# Practice Datasets
### *Working with Biological Data in Python and R*

A **fresh practice dataset for every chapter that uses data** (Chapters 3–18 and 20).
Each one is a *different biological scenario* from the worked example in the book but
exercises the **same skills**, so you can practise on data you haven't seen analysed.

All data are **simulated, biological** (no clinical or public-health data) and fully
reproducible from `make_practice_data.py` and `make_practice_bio.py`. Because they are
simulated, the "right answer" is known — noted below where relevant.

> Chapters 1, 2, 19, and 21 are conceptual (setup, reproducibility, next steps) and have
> no dataset.

---

## Chapter 3 — Reading Data
**Files:** `ch03_pollinators.csv`, `ch03_pollinators.tsv`, `ch03_pollinators_eu.csv`
Pollinator visit counts at four habitats. The **same table in three formats**: a comma CSV,
a tab-separated TSV, and a European-style file (`;` separator, `,` decimal) — the classic
import trap.
Columns: `site, species, visits, temp_C, sugar_pct`.
**Try this:** read all three correctly into one identical data frame.

## Chapter 4 — Inspecting Data Frames
**File:** `ch04_bird_survey.csv` — a bird mist-netting survey (140 birds, mixed types, a few
missing values to discover).
Columns: `ring_id, species, sex, age, wing_mm, mass_g, site`.
**Try this:** `head`/`str`/`summary` (R) or `head`/`info`/`describe` (Python); count missing
values per column.

## Chapter 5 — Cleaning Data
**File:** `ch05_amphibian_survey_messy.csv` — a **deliberately messy** amphibian capture log.
Issues to fix: inconsistent species spellings & capitalisation, mixed date formats, sex coded
many ways (`M/male/F/female/?`), lengths with stray units (`75mm`) and decimal commas, blank
/`NA`/`-` missing values, trailing spaces, and duplicate rows.
Columns: `capture_id, date, species, sex, length_mm, mass_g, site, notes`.
**Try this:** produce a clean, tidy table with consistent types and 3 canonical species.

## Chapter 6 — Correlation & Allometry
**File:** `ch06_snail_allometry.csv` — land-snail shell measurements (160 snails).
Strong **log–log (allometric)** relationships between size traits.
Columns: `snail_id, site, shell_length_mm, shell_width_mm, aperture_mm, body_mass_g, whorls`.
**Try this:** correlate traits; fit `log(mass) ~ log(length)` and interpret the slope.

## Chapter 7 — Exploratory Visualisation
**File:** `ch07_coral_transects.csv` — long-format coral genus counts across reef zones & depth.
Columns: `transect, zone, depth_m, genus, colony_count, bleached`.
**Try this:** bar charts by genus, boxplots by zone, count vs depth; spot the depth trend.

## Chapter 8 — ggplot2 (a designed experiment)
**File:** `ch08_maize_nitrogen.csv` — maize under three nitrogen levels (control/low/high, n=50).
Real dose response in height, leaf count, and an N-responsive gene.
Columns: `treatment, plant_height_cm, leaf_count, N_gene_expr, housekeeping_expr, tasseled`.
**Try this:** boxplots/violins by treatment; colour by `tasseled`; facet the gene vs housekeeping.

## Chapter 9 — matplotlib / seaborn (morphometrics)
**File:** `ch09_bumblebee_morphometrics.csv` — three bumblebee species, five traits each (165 bees).
*B. pascuorum* is long-tongued — the classes separate.
Columns: `species, body_length_mm, wing_length_mm, tongue_length_mm, thorax_width_mm, corbicula_mm`.
**Try this:** pairplot / scatter by species; histogram of tongue length.

## Chapter 10 — Describing Data (distributions)
**File:** `ch10_aphid_population.csv` — 450 aphids with traits from four distribution shapes:
normal (`body_length_mm`, `walk_speed_mm_s`), **log-normal/right-skewed** (`body_mass_mg`),
**Poisson count** (`offspring`), and **Bernoulli** (`winged`, ~37%).
**Try this:** compute mean/median/SD/IQR; show why mean ≠ median for body mass.

## Chapter 11 — Group Comparison (t-test / ANOVA)
**File:** `ch11_antibiotic_zones.csv` — inhibition-zone diameters for control + two antibiotics
(40 plates each), with a `strain` factor.
Columns: `treatment, strain, zone_diameter_mm`.
True means ≈ control 8, penicillin 17.5, streptomycin 21.
**Try this:** one-way ANOVA across treatments; Welch t-test for two of them; check assumptions.

## Chapter 12 — Relationships & Regression
**Files:**
- `ch12_yeast_growth.csv` — `time_h, density_OD` → **logistic growth curve** (plateau ≈ 1.3).
- `ch12_enzyme_kinetics.csv` — `substrate_uM, rate_umol_min` → **Michaelis–Menten** (Vmax≈2.4, Km≈6).
- `ch12_beetle_survival.csv` — `body_size_mm, overwinter_survived` → **logistic regression**
  (survival rises from ~0.1 at small size to ~1.0 at large size).
**Try this:** fit a linear model, a nonlinear curve, and a logistic GLM respectively.

## Chapter 13 — Statistical Pitfalls
**Files:**
- `ch13_metabolite_screen.csv` — `metabolite, p_value, log2FC, is_true_hit` for 300 metabolites.
  There are **18 true hits**, but ~31 metabolites have raw *p* < 0.05 — practise **multiple-testing
  correction** (BH/FDR) and compare recovered hits against `is_true_hit`.
- `ch13_mesocosm_growth.csv` — `tank, treatment, fish_id, growth_mm`. Treatment is assigned at the
  **tank** level (8 tanks) with 12 fish each — a **pseudoreplication** trap; the experimental unit
  is the tank, not the fish.
**Try this:** apply `p.adjust`/`multipletests`; then analyse the mesocosm at tank vs fish level and
see how the p-value changes.

## Chapter 14 — Sequences
**Files:** `ch14_minigenome.fasta` (one ~2.7 kb contig with three genes and promoters),
`ch14_genes.fasta` (the three CDS: `dnaK, rpoB, gyrA` — each starts ATG, no internal stops,
ends in a stop codon), `ch14_plasmid.fasta` (a reporter plasmid with EcoRI/BamHI/HindIII/PstI sites).
**Try this:** GC content, reverse complement, translate the CDS, find restriction sites.

## Chapter 15 — Genome Annotation
**Files:** `ch15_annotation.gff3` (4 genes → mRNA → exon/CDS on `chr2`),
`ch15_genes.bed` (gene spans, 0-based), `ch15_peaks.bed` (10 ChIP-style peaks, some over promoters).
**Try this:** parse the GFF3, convert coordinates, and find which peaks overlap a promoter.

## Chapter 16 — RNA-seq / Differential Expression
**Files:** `ch16_counts.csv` (600 genes × 8 samples), `ch16_metadata.csv` (mock vs infected, 4+4),
`ch16_truth.csv` (`gene, is_DE, true_log2FC` — **60 truly DE genes**).
**Try this:** normalise, run a DE test (DESeq2 / PyDESeq2 or a moderated t-test), and check your
hits and FDR against `is_DE`. Observed |log2FC| ≈ 1.6 for DE genes vs ≈ 0.4 for the rest.

## Chapter 17 — Ecology & Biodiversity
**Files:** `ch17_soil_arthropods.csv` (18 sites × 30 taxa count matrix),
`ch17_sites_meta.csv` (`site, land_use, pH` — forest / grassland / cropland, 6 each).
Cropland is the least rich; the three land uses favour different taxa.
**Try this:** Shannon/Simpson diversity per site, a rarefaction or Bray–Curtis PCoA, and PERMANOVA
by land use.

## Chapter 18 — Phylogenetic Trees
**Files:** `ch18_bird_tree.nwk` (13 songbird taxa, branch lengths + support, three clades),
`ch18_tips_meta.csv` (`tip, clade, migratory`).
**Try this:** read and plot the tree (ape / ete3), colour tips by clade, map the migratory trait,
compute a patristic distance.

## Chapter 20 — Machine Learning (classification)
**File:** `ch20_wheat_expression.csv` — 150 samples × 18 genes + `cultivar`
(three wheat cultivars: Maris / Cadenza / Paragon, 50 each).
**Six marker genes** (`gene01`–`gene06`) carry the signal; the other 12 are noise.
**Try this:** train/test split, scale, fit a classifier, evaluate with a confusion matrix and
cross-validation; rank features and confirm gene01–gene06 come out on top.

---

### Regenerate
```bash
python make_practice_data.py     # tabular datasets (Ch 3–13, 16, 17, 20)
python make_practice_bio.py      # sequences, annotation, tree (Ch 14, 15, 18)
```
Seeds are fixed, so the files are identical every time.
