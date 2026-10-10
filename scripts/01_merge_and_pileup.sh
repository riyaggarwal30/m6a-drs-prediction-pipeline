# 1. merging bam files.
samtools merge -@ 8 data/processed/hsap_combined.bam \
  data/processed/hsap_run1_merged.bam data/processed/hsap_run2_merged.bam

#indexing the file
samtools index -@ 8 data/processed/hsap_combined.bam

#running modkit pileup 
modkit pileup data/processed/hsap_combined.bam \
  data/processed/hsap_modkit_pileup.bed \
  --ref /projects/renlab/people/kgj130/data/genome/combined.fa \
  --motif DRACH 2 \
  --threads 8 \
  --log-filepath data/processed/modkit_combined_run.log
