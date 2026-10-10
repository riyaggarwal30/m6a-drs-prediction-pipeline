awk -F'\t' 'BEGIN{OFS="\t"} {
  if ($1 ~ /^HSAP_/) {
    sub(/^HSAP_/, "chr", $1)
    if ($1 == "chrMT") $1 = "chrM"
  }
  print
}' data/processed/hsap_modkit_pileup.bed > data/processed/hsap_modkit_pileup_chr.bed