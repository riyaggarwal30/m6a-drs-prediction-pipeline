#!/bin/bash
# Extract per-read m6A probabilities and base qualities at the target sites.
# Input:  merged BAM from 01_, target site BED from 03_/04_
# Output: per-read table used to compute the probability and base quality features
# Tool:   modkit v0.4.2
#
# Run on the server; commands were originally run interactively.

modkit extract full \
  hsap_combined.bam \
  per_read_mod_calls_4k.tsv \
  --include-bed target_sites_4k.bed \
  --mapped-only