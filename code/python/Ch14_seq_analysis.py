"""
Chapter 14 - DNA and Protein Sequences (Biopython)
A complete sequence-analysis pipeline on the simulated genome.

    pip install biopython
    python Ch14_seq_analysis.py

Data: genome.fasta, genes.fasta, proteins.fasta, plasmid.fasta.
Expected output (verified): genome 1629 bp / 49.6% GC; gene1 61% GC / 150 aa;
gene2 35% GC; 6 ORFs >= 80 aa (longest 154 aa).
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
import re

def hdr(t): print("\n" + "="*56 + "\n" + t + "\n" + "-"*56)

# ---------------------------------------------------- read & basic stats
hdr("1. Read the genome and basic statistics")
chrom = SeqIO.read("genome.fasta", "fasta")
print(f"id={chrom.id}  length={len(chrom.seq)} bp  GC={gc_fraction(chrom.seq)*100:.1f}%")

hdr("2. The central dogma on a short sequence")
dna = Seq("ATGGCATTAGACTAA")
print("DNA           :", dna)
print("complement    :", dna.complement())
print("rev-complement:", dna.reverse_complement())
print("transcribe    :", dna.transcribe())
print("translate     :", dna.translate())          # MALD*

# ---------------------------------------------------- per-gene GC & translation
hdr("3. GC content and translation of each gene")
for rec in SeqIO.parse("genes.fasta", "fasta"):
    prot = rec.seq.translate(to_stop=True)
    print(f"{rec.id}: {len(rec.seq)} bp  GC={gc_fraction(rec.seq)*100:.1f}%  "
          f"protein={len(prot)} aa  MW={ProteinAnalysis(str(prot)).molecular_weight()/1000:.1f} kDa")

# ---------------------------------------------------- six-frame ORF finder
hdr("4. Find open reading frames (>= 80 aa)")
def find_orfs(seq, min_aa=50):
    orfs = []
    for strand, s in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            window = s[frame:]
            window = window[:len(window) // 3 * 3]      # trim to whole codons (no warning)
            prot = str(window.translate())
            for piece in prot.split("*"):
                if "M" in piece:
                    orf = piece[piece.index("M"):]
                    if len(orf) >= min_aa:
                        orfs.append((strand, frame, len(orf), orf))
    return sorted(orfs, key=lambda o: -o[2])

orfs = find_orfs(chrom.seq, 80)
print(f"{len(orfs)} ORFs found; longest {orfs[0][2]} aa")
for i, (strand, frame, n, prot) in enumerate(orfs, 1):
    mw = ProteinAnalysis(prot).molecular_weight()
    print(f"  ORF{i}: {n} aa, strand {strand:+d}, {mw/1000:.1f} kDa")

# ---------------------------------------------------- motifs & restriction sites
hdr("5. Motifs and restriction sites")
seq = str(chrom.seq)
print("Pribnow box TATAAT at:", [m.start() for m in re.finditer("TATAAT", seq)])
try:
    from Bio.Restriction import EcoRI, BamHI, HindIII
    print("EcoRI sites :", EcoRI.search(chrom.seq))
    print("BamHI sites :", BamHI.search(chrom.seq))
    print("HindIII     :", HindIII.search(chrom.seq))
except Exception as e:
    print("(Bio.Restriction:", e, ")")

# ---------------------------------------------------- k-mers & composition
hdr("6. k-mers and base composition")
k = 3
kmers = Counter(seq[i:i+k] for i in range(len(seq)-k+1))
print("most common 3-mers:", kmers.most_common(3))
print("base composition  :", {b: round(seq.count(b)/len(seq)*100, 1) for b in "ACGT"})

# ---------------------------------------------------- protein properties (ALL genes)
hdr("7. Protein properties (every gene)")
for rec in SeqIO.parse("genes.fasta", "fasta"):          # loop over all records
    prot = rec.seq.translate(to_stop=True)
    pa = ProteinAnalysis(str(prot))
    print(f"{rec.id}: length={len(prot)} aa  MW={pa.molecular_weight()/1000:.1f} kDa  "
          f"pI={pa.isoelectric_point():.1f}  GRAVY={pa.gravy():.2f}  "
          f"instability={pa.instability_index():.1f}")

print("\nDone.")
