#!/usr/bin/env python

import argparse
import os, sys
import datetime
import joblib
import gzip
import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import tensorflow as tf
import tensorflow_hub as hub

VERSION = 0.1
script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
TRANSFORMER_PATH = f'{script_directory}/model/enformer.finetuned.SAD.robustscaler-PCA500-robustscaler.transform.pkl'
MODEL_PATH = 'https://tfhub.dev/deepmind/enformer/1'
TARGETS_TXT = 'https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt'
SEQUENCE_LENGTH = 393216

class Enformer:
  
  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model
  
  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}
  
  @tf.function
  def contribution_input_grad(self, input_sequence,
                              target_mask, output_head='human'):
    input_sequence = input_sequence[tf.newaxis]
    
    target_mask_mass = tf.reduce_sum(target_mask)
    with tf.GradientTape() as tape:
      tape.watch(input_sequence)
      prediction = tf.reduce_sum(
          target_mask[tf.newaxis] *
          self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass
    
    input_grad = tape.gradient(prediction, input_sequence) * input_sequence
    input_grad = tf.squeeze(input_grad, axis=0)
    return tf.reduce_sum(input_grad, axis=-1)

class EnformerScoreVariantsRaw:
  
  def __init__(self, tfhub_url, organism='human'):
    self._model = Enformer(tfhub_url)
    self._organism = organism
  
  def predict_on_batch(self, inputs):
    ref_prediction = self._model.predict_on_batch(inputs['ref'])[self._organism]
    alt_prediction = self._model.predict_on_batch(inputs['alt'])[self._organism]
    
    return {
      'ref': ref_prediction,
      'alt': alt_prediction,
      'diff_raw': alt_prediction[0] - ref_prediction[0],
      'diff_mean': alt_prediction.mean(axis=1) - ref_prediction.mean(axis=1)}


class EnformerScoreVariantsNormalized:
  
  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human'):
    assert organism == 'human', 'Transforms only compatible with organism=human'
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      transform_pipeline = joblib.load(f)
    self._transform = transform_pipeline.steps[0][1]  # StandardScaler.
  
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return {
      'ref': scores['ref'],
      'alt': scores['alt'],
      'diff_raw': scores['diff_raw'],
      'diff_norm': self._transform.transform(scores['diff_mean'])}


class EnformerScoreVariantsPCANormalized:
  
  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human', num_top_features=500):
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      self._transform = joblib.load(f)
    self._num_top_features = num_top_features
  
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self._transform.transform(scores['diff_raw'])[:, :self._num_top_features]
  
class FastaStringExtractor:
  
  def __init__(self, fasta_file):
      self.fasta = pyfaidx.Fasta(fasta_file)
      self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}
  
  def extract(self, interval: Interval, **kwargs) -> str:
      # Truncate interval if it extends beyond the chromosome lengths.
      chromosome_length = self._chromosome_sizes[interval.chrom]
      trimmed_interval = Interval(interval.chrom,
                                  max(interval.start, 0),
                                  min(interval.end, chromosome_length),
                                  )
      # pyfaidx wants a 1-based interval
      sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                        trimmed_interval.start + 1,
                                        trimmed_interval.stop).seq).upper()
      # Fill truncated values with N's.
      pad_upstream = 'N' * max(-interval.start, 0)
      pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
      return pad_upstream + sequence + pad_downstream
  
  def close(self):
      return self.fasta.close()

def log(message, level="INFO"):
  #Print a log message as [timestamp] - level - message
  print(f'[{datetime.datetime.now()}] - {level} - {message}')

def variant_generator(variant_input, gzipped=False):
  """Yields a kipoiseq.dataclasses.Variant for each variant"""
  def _open(file):
    return gzip.open(file, 'rt') if gzipped else open(file)
  
  if variant_input is None:
    return []
  
  if isinstance(variant_input, str):
    # When variant_input is a string we assume a file path, read it line by line.
    with _open(variant_input) as f:
      for line in f:
        if line.startswith('#'):
          continue
        chrom, pos, id, ref, alt_list = line.split('\t')[:5]
        # Split ALT alleles and return individual variants as output.
        for alt in alt_list.split(','):
          yield kipoiseq.dataclasses.Variant(chrom=chrom, pos=pos,
                                            ref=ref, alt=alt, id=id)
  else:
    # When variant_input is a list, yield each variant.
    for v in variant_input:
      yield kipoiseq.dataclasses.Variant(chrom=v[0], pos=v[1], ref=v[2], alt=v[3], 
                                         id=v[4] if len(v) == 5 else '_'.join(v))

def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

def variant_centered_sequences(variant_input, sequence_length, fasta_file, gzipped=False,
                               chr_prefix=''):
  result = []
  seq_extractor = kipoiseq.extractors.VariantSeqExtractor(
    reference_sequence=FastaStringExtractor(fasta_file))
  
  for variant in variant_generator(variant_input, gzipped=gzipped):
    try:
      interval = Interval(chr_prefix + variant.chrom,
                          variant.pos, variant.pos)
      interval = interval.resize(sequence_length)
      center = interval.center() - interval.start
      
      reference = seq_extractor.extract(interval, [], anchor=center)
      alternate = seq_extractor.extract(interval, [variant], anchor=center)
      
      result.append( {'interval': interval,
        'inputs': {'ref': one_hot_encode(reference),
                        'alt': one_hot_encode(alternate)},
            'metadata': {'chrom': chr_prefix + variant.chrom,
                          'pos': variant.pos,
                          'id': variant.id,
                          'ref': variant.ref,
                          'alt': variant.alt}})
    except Exception as e:
      log(f'Error processing variant {variant.chrom}:{variant.pos}:{variant.ref}:{variant.alt}', 'WARNING')
      print("The following exception was generated")
      print(e)
      continue
  return result
    
def plot_tracks(tracks, interval, height=1.5):
  fig, axes = plt.subplots(len(tracks), 1, figsize=(20, height * len(tracks)), sharex=True)
  for ax, (title, y) in zip(axes, tracks.items()):
    ax.fill_between(np.linspace(interval.start, interval.end, num=len(y)), y)
    ax.set_title(title)
    sns.despine(top=True, right=True, bottom=True)
  ax.set_xlabel(str(interval))
  plt.tight_layout()
  return(fig)

def read_input(input, input_type):
  # If input is a file, read it line by line
  if os.path.isfile(input):
    with open(input) as f:
      input = f.readlines()
    # Remove newlines
    input = [x.strip() for x in input]
    if input_type == 'intervals':
      # Split intervals
      input = [x.split('\t') for x in input]
    else:
      # Split variants
      input = [x.split('_') for x in input]
  # If input is a comma-separated list, split it
  if type(input) == str:
    input = input.split(',')
    input = [x.split('_') for x in input]
  
  # Log the number of inputs
  if len(input) == 0:
    raise ValueError(f'No {input_type} found')
  log(f'Number of {input_type} from list: {len(input)}')
  return input

def main():
  log(f'== ENFORMER PREDICTION SCRIPT {VERSION} ==')
  #Parse arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('--vcf', type=str, help='VCF.gz file with variants to predict')
  parser.add_argument('--variants', type=str, help='List of variants id to predict in the foramt CHR_POS_REF_ALT(_ID). Can be a comma-separated list or a file with one variant per line.')

  parser.add_argument('--output', type=str, required=False, help='Output folder', default='enformer_output')
  parser.add_argument('--ref_genome', type=str, required=True, help='Reference genome fasta file. A fai index must be present for this file.')
  parser.add_argument('--tfhub_url', type=str, default=MODEL_PATH)
  parser.add_argument('--transform_pkl_path', type=str, default=TRANSFORMER_PATH)
  parser.add_argument('--plot_features', type=str, default="none", help='Comma separated list of feature indexes to use in output. Use all to plot for all features.')
  parser.add_argument('--add_chr_prefix', action='store_true', help='add chr prefix to chromosome names when reading input variants')
  args = parser.parse_args()

  log(f'== TEST GPU IS AVAILABLE ==')
  assert tf.config.list_physical_devices('GPU'), 'No GPU found. Please enable GPUs to run this script.'

  # Load targets meta-data
  df_targets = pd.read_csv(TARGETS_TXT, sep='\t')

  # Print targets meta-data if requested
  # if args.list_features:
  #   print('== AVAIABLE FEATURES ==')
  #   # Print the whole df_targets dataframe, limited to columns index and description
  #   print(df_targets[['index', 'description']].to_markdown(index=False))
  #   sys.exit()

  # Set up variables
  model_path = args.tfhub_url
  transform_path = args.transform_pkl_path
  vcf_input = args.vcf
  fasta_file = args.ref_genome
  output_folder = args.output
  os.makedirs(output_folder, exist_ok=True)
  chr_prefix = 'chr' if args.add_chr_prefix else ''

  # Init procedures
  log(f'Model path: {model_path}')
  log(f'Transform path: {transform_path}')
  
  if not os.path.isfile(fasta_file):
    raise ValueError(f'Fasta file not found: {fasta_file}')
  if not os.path.isfile(fasta_file + '.fai'):
    raise ValueError(f'Fasta index not found for {fasta_file}. Please run samtools faidx {fasta_file}.')
  
  log(f'Fasta file: {fasta_file}')
  
  gzipped = False
  if vcf_input:
    if vcf_input.endswith('.gz'):
      gzipped = True
    if not os.path.isfile(vcf_input):
      raise ValueError(f'VCF file not found: {vcf_input}')
  variants_input = read_input(args.variants, 'variants') if args.variants else []
  if len(variants_input) == 0 and not vcf_input:
    raise ValueError('No input variants provided. Please provide a VCF file or a list of variants.')

  plot_features_idx = []
  if args.plot_features != 'none':
    if args.plot_features == 'all':
      plot_features_idx = df_targets.index
      log(f'Plotting all features ({len(plot_features_idx)})')
    else:
      plot_features_idx = [int(x) for x in args.plot_features.split(',')]
      log(f'Plotting {len(plot_features_idx)} features: {plot_features_idx}')

  log("== START PREDICTION ==")
  # Predict for variants
  vars_from_vcf = variant_centered_sequences(vcf_input, sequence_length=SEQUENCE_LENGTH,
                                gzipped=gzipped, chr_prefix=chr_prefix, fasta_file=fasta_file)
  vars_from_list = variant_centered_sequences(variants_input, sequence_length=SEQUENCE_LENGTH,
                                chr_prefix=chr_prefix, fasta_file=fasta_file) 
  
  it = vars_from_vcf + vars_from_list

  log(f'Total number of variants to predict: {len(it)}')

  # Extract scores using the first 20 PCs
  log("== EXTRACT SCORES USING THE FIRST 20 PCs ==")
  example_list = []
  enformer_score_variants = EnformerScoreVariantsPCANormalized(model_path, transform_path, num_top_features=20)
  for i, example in enumerate(it):
    variant_scores = enformer_score_variants.predict_on_batch(
        {k: v[tf.newaxis] for k,v in example['inputs'].items()})[0]
    variant_scores = {f'PC{i}': score for i, score in enumerate(variant_scores)}
    example_list.append({**example['metadata'],
                        **variant_scores})
    if (i+1) % 10 == 0:
      print(f'{i+1} variants processed')
  df = pd.DataFrame(example_list)
  df.to_csv(f'{output_folder}/variant_scores.PC20.tsv', sep='\t', index=False)
  log(f'Saved scores to {output_folder}/variant_scores.PC20.tsv')

  # Extract all scores for each variant
  log("== EXTRACT ALL SCORES FOR EACH VARIANT ==")
  if len(plot_features_idx) > 0:
    os.makedirs(f'{output_folder}/plots', exist_ok=True)
    log(f'Plots will be saved for each var to {output_folder}/plots')

  example_list = []
  #plotting_list = []
  enformer_score_variants_all = EnformerScoreVariantsNormalized(model_path, transform_path)

  for i, example in enumerate(it):
    variant_scores = enformer_score_variants_all.predict_on_batch(
        {k: v[tf.newaxis] for k,v in example['inputs'].items()})
    variant_scores_df = {f'{i}_{name[:20]}': score for i, (name, score) in enumerate(zip(df_targets.description, variant_scores['diff_norm'][0]))}
    if len(plot_features_idx) > 0:
      plotting_scores_ref = {f'REF_{name[:20]}': score for (name, score) in zip(df_targets.description[plot_features_idx], variant_scores['ref'][0,:,plot_features_idx])}
      plotting_scores_diff = {f'DIFF_{name[:20]}': score for (name, score) in zip(df_targets.description[plot_features_idx], variant_scores['diff_raw'][np.newaxis][0,:,plot_features_idx])}
      variant_track = np.zeros_like(variant_scores['ref'][0,:, 0], dtype=bool)
      variant_track[variant_track.shape[0] // 2] = True
      plot_figure = plot_tracks({'variant':variant_track, **plotting_scores_ref, **plotting_scores_diff}, example['interval'].resize(variant_scores['ref'].shape[0] * 128), height=1)
      plot_figure.savefig(f'{output_folder}/plots/{example["metadata"]["id"]}.png')
    #plotting_list.append({**example['metadata'],**plotting_scores})
    example_list.append({**example['metadata'],
                        **variant_scores_df})
    if (i+1) % 10 == 0:
      print(f'{i+1} variants processed')
  df = pd.DataFrame(example_list)
  df.to_csv(f'{output_folder}/variant_scores.all.tsv', sep='\t', index=False)
  log(f'Saved scores to {output_folder}/variant_scores.all.tsv')

  log("== ALL DONE ==")

if __name__ == '__main__':
  main()
