# =====================================================================
# Chapter 14 - DNA and Protein Sequences (R / Bioconductor)
# The same pipeline as seq_analysis.py, using Biostrings (+ Peptides).
#
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("Biostrings")
#   install.packages("Peptides")          # protein properties
#   source("seq_analysis.R")
#
# Data: genome.fasta, genes.fasta, proteins.fasta, plasmid.fasta.
# Expected (verified): genome 1629 bp / 49.6% GC; gene1 61% GC / 150 aa;
# gene2 35% GC; 6 ORFs >= 80 aa (longest 154 aa).
# =====================================================================
suppressMessages(library(Biostrings))

cat("\n== 1. Read the genome and basic statistics ==\n")
genome <- readDNAStringSet("genome.fasta")
chrom  <- genome[[1]]                       # a DNAString
gc <- function(x) sum(letterFrequency(x, "GC")) / length(x) * 100
cat(sprintf("id=%s  length=%d bp  GC=%.1f%%\n", names(genome)[1], length(chrom), gc(chrom)))

cat("\n== 2. The central dogma on a short sequence ==\n")
dna <- DNAString("ATGGCATTAGACTAA")
cat("DNA           :", as.character(dna), "\n")
cat("complement    :", as.character(complement(dna)), "\n")
cat("rev-complement:", as.character(reverseComplement(dna)), "\n")     # TTAGTCTAATGCCAT
cat("transcribe    :", as.character(RNAString(dna)), "\n")             # AUGGCAUUAGACUAA
cat("translate     :", as.character(translate(dna)), "\n")            # MALD*

cat("\n== 3. GC content and translation of each gene ==\n")
genes <- readDNAStringSet("genes.fasta")
for (i in seq_along(genes)) {
  g <- genes[[i]]
  prot <- translate(g)                       # AAString, '*' = stop
  prot <- subseq(prot, 1, nchar(gsub("\\*.*$", "", as.character(prot))))  # to first stop
  cat(sprintf("%s: %d bp  GC=%.1f%%  protein=%d aa\n",
              names(genes)[i], length(g), gc(g), nchar(as.character(prot))))
}

cat("\n== 4. Find open reading frames (>= 80 aa) ==\n")
find_orfs <- function(dna, min_aa = 50) {
  out <- list()
  for (strand in c(1, -1)) {
    s <- if (strand == 1) dna else reverseComplement(dna)
    for (frame in 0:2) {
      w <- ((length(s) - frame) %/% 3) * 3
      prot <- as.character(translate(subseq(s, frame + 1, width = w)))
      for (piece in strsplit(prot, "\\*", fixed = FALSE)[[1]]) {
        m <- regexpr("M", piece)
        if (m > 0) {
          orf <- substr(piece, m, nchar(piece))
          if (nchar(orf) >= min_aa)
            out[[length(out) + 1]] <- list(strand = strand, frame = frame,
                                           len = nchar(orf), prot = orf)
        }
      }
    }
  }
  out[order(-sapply(out, function(o) o$len))]
}
orfs <- find_orfs(chrom, 80)
cat(sprintf("%d ORFs found; longest %d aa\n", length(orfs), orfs[[1]]$len))
for (i in seq_along(orfs))
  cat(sprintf("  ORF%d: %d aa, strand %+d\n", i, orfs[[i]]$len, orfs[[i]]$strand))

cat("\n== 5. Motifs and restriction sites ==\n")
cat("Pribnow box TATAAT at:",
    start(matchPattern("TATAAT", chrom)), "\n")
for (site in c(EcoRI = "GAATTC", BamHI = "GGATCC", HindIII = "AAGCTT"))
  cat(sprintf("  %s at: %s\n", site, paste(start(matchPattern(site, chrom)), collapse = ", ")))

cat("\n== 6. k-mers and base composition ==\n")
print(head(sort(oligonucleotideFrequency(chrom, width = 3), decreasing = TRUE), 3))
print(round(letterFrequency(chrom, c("A", "C", "G", "T"), as.prob = TRUE) * 100, 1))

cat("\n== 7. Protein properties (gene1, needs the Peptides package) ==\n")
prot1 <- as.character(translate(genes[[1]]))
prot1 <- gsub("\\*.*$", "", prot1)            # trim at first stop
if (requireNamespace("Peptides", quietly = TRUE)) {
  cat(sprintf("length=%d aa  MW=%.1f kDa  pI=%.1f  GRAVY=%.2f\n",
              nchar(prot1), Peptides::mw(prot1) / 1000,
              Peptides::pI(prot1), Peptides::hydrophobicity(prot1)))
} else {
  cat("install.packages('Peptides') for MW, pI, hydrophobicity\n")
}

cat("\nDone.\n")
