"""
Chapter 5 - Cleaning Messy Data (Python)
Runnable examples on a messy reef-fish field dataset that ships with the book.

    python clean_examples.py

Needs: pandas, numpy   (conda install pandas numpy)
"""
import pandas as pd, numpy as np

# =====================================================================
# 1. READ CORRECTLY  (turn placeholders into real missing values)
# =====================================================================
raw = pd.read_csv("specimens_messy.csv")
print("raw shape:", raw.shape)                       # (123, 5)
print("messy 'habitat' spellings:", raw["habitat"].nunique())   # 18

fish = pd.read_csv("specimens_messy.csv", na_values=["", "NA", "-99"])
print("dtypes after na_values:\n", fish.dtypes, sep="")

# =====================================================================
# 2. MISSING VALUES
# =====================================================================
print("\nmissing per column:\n", fish.isna().sum(), sep="")
print("rows with any missing:", fish.isna().any(axis=1).sum())        # 5
print("pct missing:\n", (fish.isna().mean()*100).round(1), sep="")

# =====================================================================
# 3. DUPLICATES
# =====================================================================
print("\nduplicate rows:", fish.duplicated().sum())                   # 3
dup_ids = sorted(fish.loc[fish["specimen_id"].duplicated(keep=False), "specimen_id"].unique())
print("duplicate IDs:", dup_ids)                                      # F010 F025 F060
fish = fish.drop_duplicates()
print("after drop_duplicates:", fish.shape)                           # (120, 5)

# =====================================================================
# 4. STANDARDISE TEXT COLUMNS
# =====================================================================
fish["habitat"] = fish["habitat"].str.strip().str.title()
fish["sex"] = fish["sex"].str.upper().str[0]
print("\nhabitats now:\n", fish["habitat"].value_counts(), sep="")
print("sex now:\n", fish["sex"].value_counts(), sep="")

# =====================================================================
# 5. OUTLIERS  (IQR vs z-score vs MAD)
# =====================================================================
L = fish["length_mm"]
q1, q3 = L.quantile(.25), L.quantile(.75); iqr = q3 - q1
low, high = q1 - 1.5*iqr, q3 + 1.5*iqr
iqr_out = fish[(L < low) | (L > high)]
print("\nIQR outliers:", iqr_out["specimen_id"].tolist())            # F026 F050 F099
print("z-score (>3) outliers:", int(((L - L.mean())/L.std()).abs().gt(3).sum()))  # 1
med = L.median(); mad = (L - med).abs().median()
print("MAD (>3.5) outliers:", int((0.6745*(L-med)/mad).abs().gt(3.5).sum()))      # 2
print("mean length WITH outliers:", round(L.mean(), 1))             # 127.0

# null out impossible values, then impute with the median
for col, lo, hi in [("length_mm", 50, 300), ("mass_g", 5, 150)]:
    fish.loc[(fish[col] < lo) | (fish[col] > hi), col] = np.nan
    fish[col] = fish[col].fillna(fish[col].median())
print("mean length AFTER cleaning:", round(fish["length_mm"].mean(), 1))  # ~119

# =====================================================================
# 6. DERIVED COLUMNS + VALIDATION  (Fulton's condition factor K)
# =====================================================================
fish = fish.rename(columns={"length_mm": "length", "mass_g": "mass"})
fish["K"] = (100000 * fish["mass"] / fish["length"]**3).round(2)
fish["K_class"] = pd.cut(fish["K"], bins=[0, 1.0, 1.5, 2.0, 100],
                         labels=["poor", "fair", "good", "high"])
assert fish["length"].between(50, 300).all()
assert fish["sex"].isin(["M", "F"]).all()
assert fish["specimen_id"].is_unique
print("\nK class counts:\n", fish["K_class"].value_counts(), sep="")
print("K max:", round(fish["K"].max(), 2))                          # ~4.2

# =====================================================================
# 7. RESHAPE  (wide tissue expression -> long)
# =====================================================================
wide = pd.read_csv("expression_wide.csv")
long = wide.melt(id_vars="specimen_id", var_name="tissue", value_name="expression")
print("\nwide", wide.shape, "-> long", long.shape)                  # (120,4) -> (360,3)

# =====================================================================
# 8. JOINS  (specimens + sequencing metadata)
# =====================================================================
gen = pd.read_csv("genetics.csv")
print("inner:", fish.merge(gen, on="specimen_id", how="inner").shape[0])   # 100
left = fish.merge(gen, on="specimen_id", how="left")
print("left :", left.shape[0], "| specimens with no sequencing record:",
      int(left["read_count"].isna().sum()))                         # 120 | 20
print("outer:", fish.merge(gen, on="specimen_id", how="outer").shape[0])   # 140

# =====================================================================
# 9. A COMBINED ANSWER  (reshape + join + group)
# =====================================================================
expr = long.merge(fish[["specimen_id", "habitat"]], on="specimen_id")
print("\nmean expression per habitat:\n",
      expr.groupby("habitat")["expression"].mean().round(1), sep="")

# =====================================================================
# 10. SAVE THE CLEAN RESULT (never overwrite the raw file)
# =====================================================================
clean = fish.merge(gen, on="specimen_id", how="left")
clean.to_csv("specimens_clean.csv", index=False)
print("\nwrote specimens_clean.csv", clean.shape)
