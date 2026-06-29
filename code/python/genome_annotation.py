"""
Chapter 15 - Genes and Genomes (Python)
Read GFF/BED, build gene models, compute overlaps, query NCBI/Ensembl.

    pip install pandas gffutils pyranges biopython requests
    python genome_annotation.py

Data: annotation.gff3, genes.bed, peaks.bed, ensembl_lookup.json.
The pandas/overlap parts run with no internet; the NCBI/Ensembl calls need it
(a cached Ensembl response is parsed as a fallback).
"""
import json
import pandas as pd

def hdr(t): print("\n" + "="*56 + "\n" + t + "\n" + "-"*56)

# ---------------------------------------------------- read GFF as a dataframe
hdr("1. Read the GFF3 and count features")
cols = ['seqid','source','type','start','end','score','strand','phase','attr']
gff = pd.read_csv('annotation.gff3', sep='\t', comment='#', names=cols)
print(gff['type'].value_counts().to_dict())
genes = gff[gff.type == 'gene'].copy()
genes['span'] = genes.end - genes.start + 1
print("strands:", genes.strand.value_counts().to_dict(),
      " longest gene span:", int(genes.span.max()))

# ---------------------------------------------------- coordinate conversion
hdr("2. Coordinate systems (GFF 1-based vs BED 0-based)")
g_start, g_end = 100, 150                      # a GFF feature
print(f"GFF {g_start}-{g_end}  ->  BED start {g_start-1}, end {g_end}")
print(f"length: BED {g_end-(g_start-1)} == GFF {g_end-g_start+1}")

# ---------------------------------------------------- gene models (gffutils)
hdr("3. Gene model via gffutils")
try:
    import gffutils
    db = gffutils.create_db('annotation.gff3', ':memory:', force=True,
                            merge_strategy='create_unique', keep_order=True)
    for g in db.features_of_type('gene'):
        n_exon = sum(1 for _ in db.children(g, featuretype='exon'))
        print(f"{g.id}: {g.end-g.start+1} bp, {n_exon} exons, strand {g.strand}")
except ImportError:
    print("(install gffutils for hierarchical gene models)")

# ---------------------------------------------------- overlaps (pyranges)
hdr("4. Gene-peak overlaps")
try:
    import pyranges as pr
    genes_r = pr.read_bed('genes.bed'); peaks_r = pr.read_bed('peaks.bed')
    n_in = len(peaks_r.overlap(genes_r))
    print(f"peaks inside a gene: {n_in}/{len(peaks_r)}")
    print("genes overlapping a peak:", sorted(set(peaks_r.join(genes_r).df['Name'])))
except ImportError:
    # pure-python fallback so the script always runs
    def load_bed(p): return [(c.split('\t')[0], int(c.split('\t')[1]), int(c.split('\t')[2]),
                              c.split('\t')[3]) for c in open(p)]
    G, P = load_bed('genes.bed'), load_bed('peaks.bed')
    ov = lambda a, b: a[1] < b[2] and b[1] < a[2]
    n_in = sum(any(ov(p, g) for g in G) for p in P)
    print(f"peaks inside a gene: {n_in}/{len(P)} (pure-python fallback)")

# ---------------------------------------------------- query Ensembl
hdr("5. Query Ensembl (REST), with cached fallback")
try:
    import requests
    url = 'https://rest.ensembl.org/lookup/id/ENSG00000139618'
    info = requests.get(url, headers={'Content-Type':'application/json'}, timeout=8).json()
except Exception:
    info = json.load(open('ensembl_lookup.json'))      # offline cache
print(info['display_name'], info['biotype'],
      f"chr{info['seq_region_name']}:{info['start']}-{info['end']}")

# ---------------------------------------------------- query NCBI Entrez
hdr("6. Query NCBI Entrez (needs internet + email)")
print("""from Bio import Entrez, SeqIO
Entrez.email = 'you@example.org'
ids = Entrez.read(Entrez.esearch(db='nucleotide',
        term='BRCA2[Gene] AND human[Organism]'))['IdList']
rec = SeqIO.read(Entrez.efetch(db='nucleotide', id=ids[0],
        rettype='gb', retmode='text'), 'genbank')
print(rec.id, len(rec.seq), rec.description)""")

print("\nDone.")
