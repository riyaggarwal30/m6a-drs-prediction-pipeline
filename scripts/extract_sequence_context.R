library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(tidyverse)

# 1. Read exported master dataframe from Python
df <- read_tsv("master_df_temp.tsv")

# 2. Construct 1-based GRanges window [start - 2, start + 2]
gr <- GRanges(
  seqnames = df$chrom,
  ranges = IRanges(start = df$start - 2, end = df$start + 2),
  strand = df$strand
)

# 3. Extract 5-mer sequences (handles reverse complement for '-' strand automatically)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
df$kmer_5mer <- as.character(getSeq(hg38, gr))

# 4. Extract position-specific columns (-2, -1, +1, +2)
df$pos_minus2 <- substr(df$kmer_5mer, 1, 1)
df$pos_minus1 <- substr(df$kmer_5mer, 2, 2)
# Index 3 is central Adenosine ('A')
df$pos_plus1  <- substr(df$kmer_5mer, 4, 4)
df$pos_plus2  <- substr(df$kmer_5mer, 5, 5)

# 5. Sanity check: Ensure central base (Index 3) is always 'A'
central_bases <- unique(substr(df$kmer_5mer, 3, 3))
cat("Central target bases found:", central_bases, "\n")
stopifnot(all(central_bases == "A"))

# 6. Save complete feature dataset
write_tsv(df, "final_m6a_feature_matrix.tsv")

cat("\n--- Sequence Context Extraction Complete! ---\n")
cat("Dimensions:", dim(df), "\n")
print(head(df[, c("chrom", "start", "strand", "kmer_5mer", "pos_minus2", "pos_minus1", "pos_plus1", "pos_plus2")]))