# Working with Biological Data in Python and R — Code & Data

Companion code and datasets for the book ***Working with Biological Data in Python and R: An Undergraduate Introduction*** by Christos Noutsos.

Everything here is **free and open**: clone or download it, run the scripts, and follow along chapter by chapter. All datasets are **simulated** (no clinical or human data) so they are safe to share and reuse.

## What's here

```
code/python/    Python scripts (pandas, matplotlib/seaborn, Biopython, scikit-learn, ...)
code/R/         R scripts (tidyverse, ggplot2, Bioconductor, vegan, ape, ...)
data/           Simulated datasets, one set per chapter (see data/README.md)
```

`data/README.md` is a full data dictionary: it lists every file, the columns it contains, the chapter it belongs to, and a short "try this" exercise.

## Quick start

**Python** (3.10+):
```bash
pip install pandas numpy matplotlib seaborn scikit-learn biopython
python code/python/Ch7_plot_examples.py
```

**R** (4.x):
```r
install.packages(c("tidyverse", "ggplot2"))
# Bioconductor / specialised packages are installed by the scripts that need them
source("code/R/Ch17_biodiversity.R")
```

Run each script from the repository root (or from the folder that holds the data
file it reads) so the relative paths resolve. Several scripts are written to
auto-create their data or fall back gracefully if an optional package is missing.

## Chapters → files (selected)

| Topic | Python | R |
|---|---|---|
| Figures / EDA | `Ch7_plot_examples.py` | `explore_examples.R` |
| Cleaning data | `clean_examples.py` | `clean_examples.R` |
| Relationships / regression | `relationships.py` | `relationships.R` |
| Sequences (Biopython / Bioconductor) | `Ch14_seq_analysis.py` | `seq_analysis.R` |
| Genome annotation | `genome_annotation.py` | `Ch15_genome_annotation.R` |
| RNA-seq differential expression | — | `Ch16_rnaseq_de.R` |
| Ecology & biodiversity | — | `Ch17_biodiversity.R` |
| Machine learning (classification) | — | `Ch20_classify.R` |

## License

- **Code** is released under the **MIT License** (see `LICENSE`) — use it freely, including in your own work.
- **Datasets** in `data/` are released under **CC-BY-4.0**: free to use and adapt with attribution.

## Citation

If these materials help your teaching or research, please cite the book:

> Noutsos, C. *Working with Biological Data in Python and R: An Undergraduate Introduction.*

Found a bug or have a suggestion? Open an issue or a pull request — contributions are welcome.
