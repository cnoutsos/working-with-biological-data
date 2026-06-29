# =====================================================================
# Chapter 15 - Genes and Genomes (R / Bioconductor)
# Read GFF/BED, build gene models, compute overlaps, query NCBI/Ensembl.
#
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("rtracklayer","GenomicFeatures","GenomicRanges",
#                          "txdbmaker"))   # txdbmaker: makeTxDbFromGFF (Bioc >= 3.19)
#   install.packages(c("rentrez","jsonlite"))   # NCBI;  biomaRt via Bioconductor
#   source("Ch15_genome_annotation.R")
#
# Data: annotation.gff3, genes.bed, peaks.bed, ensembl_lookup.json.
# =====================================================================

cat("\n== 1. Read the GFF3 and count features ==\n")
cols <- c('seqid','source','type','start','end','score','strand','phase','attr')
gff <- read.delim('annotation.gff3', comment.char='#', header=FALSE, col.names=cols)
print(table(gff$type))                       # exon 14, CDS 14, gene 5, ...
genes_df <- subset(gff, type == 'gene')
genes_df$span <- genes_df$end - genes_df$start + 1
cat("strands:", paste(names(table(genes_df$strand)), table(genes_df$strand)),
    " longest span:", max(genes_df$span), "\n")

cat("\n== 2. Coordinate systems (GFF 1-based vs BED 0-based) ==\n")
cat("GFF 100-150  ->  BED start 99, end 150;  length 51 either way\n")

cat("\n== 3. Gene models & ranges (rtracklayer / GenomicFeatures) ==\n")
suppressMessages({
  ok <- requireNamespace("rtracklayer",     quietly = TRUE) &&
        requireNamespace("GenomicFeatures", quietly = TRUE) &&
        requireNamespace("GenomicRanges",   quietly = TRUE)
})
if (ok) {
  suppressMessages({
    library(rtracklayer); library(GenomicFeatures); library(GenomicRanges)
  })

  # makeTxDbFromGFF() moved from GenomicFeatures to the txdbmaker package in
  # Bioconductor 3.19. Use whichever provides it so this runs on old and new installs.
  build_txdb <- function(f) {
    if (requireNamespace("txdbmaker", quietly = TRUE))
      txdbmaker::makeTxDbFromGFF(f)
    else
      GenomicFeatures::makeTxDbFromGFF(f)
  }

  # The whole ranges section degrades gracefully if a package/step is missing.
  res <- tryCatch({
    gr   <- import('annotation.gff3')        # a GRanges of all features
    # suppressWarnings() hides the harmless GFF-parsing notes (phase column,
    # orphan transcripts) that makeTxDbFromGFF prints on small toy files.
    txdb <- suppressWarnings(build_txdb('annotation.gff3'))
    ex_by_gene <- exonsBy(txdb, by = 'gene')
    print(lengths(ex_by_gene))               # exons per gene

    cat("\n== 4. Gene-peak overlaps ==\n")
    genes <- import('genes.bed'); peaks <- import('peaks.bed')
    # put both objects on the same seqlevels so overlap ops don't warn
    common <- union(seqlevels(genes), seqlevels(peaks))
    seqlevels(genes) <- common; seqlevels(peaks) <- common

    cat("peaks inside a gene:", sum(countOverlaps(peaks, genes) > 0), "/", length(peaks), "\n")
    hits <- findOverlaps(peaks, genes)
    cat("genes overlapping a peak:",
        paste(unique(genes$name[subjectHits(hits)]), collapse = ", "), "\n")

    # which gene is a variant in?  (build the point on the genes' own seqlevels)
    chr <- as.character(seqnames(genes))[1]
    v <- GRanges(chr, IRanges(5301, 5301))
    seqlevels(v) <- seqlevels(genes)
    hit_gene <- genes$name[countOverlaps(genes, v) > 0]
    cat("variant at ", chr, ":5301 is in: ",
        if (length(hit_gene)) paste(hit_gene, collapse = ", ") else "(intergenic)",
        "\n", sep = "")
    TRUE
  }, error = function(e) {
    cat("(ranges step skipped:", conditionMessage(e), ")\n"); FALSE
  })
} else {
  cat("(install rtracklayer + GenomicFeatures + txdbmaker for ranges/overlaps)\n")
}

cat("\n== 5. Query Ensembl ==\n")
# live REST (needs internet):
#   library(httr); r <- GET("https://rest.ensembl.org/lookup/id/ENSG00000139618",
#                           accept_json()); content(r)
# or biomaRt:
#   library(biomaRt); mart <- useMart('ensembl', 'hsapiens_gene_ensembl')
#   getBM(c('hgnc_symbol','chromosome_name','start_position','end_position'),
#         'hgnc_symbol', 'BRCA2', mart)
# offline cache:
if (requireNamespace("jsonlite", quietly = TRUE)) {
  info <- jsonlite::fromJSON('ensembl_lookup.json')
  cat(sprintf("%s  %s  chr%s:%d-%d\n", info$display_name, info$biotype,
              info$seq_region_name, info$start, info$end))
} else {
  cat("(install jsonlite to read the cached Ensembl lookup)\n")
}

cat("\n== 6. Query NCBI Entrez (needs internet) ==\n")
cat("library(rentrez)
hits <- entrez_search(db='nucleotide', term='BRCA2[Gene] AND human[Organism]')
gb <- entrez_fetch(db='nucleotide', id=hits$ids[1], rettype='gb', retmode='text')
cat(substr(gb, 1, 400))\n")

cat("\nDone.\n")
