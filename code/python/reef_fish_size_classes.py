"""
Reef fish: assign length-based size classes and get cumulative counts.
Cut-offs are derived from the ACTUAL distribution (quartiles), so the four
groups are balanced instead of nearly empty at the ends.

Run:  python reef_fish_size_classes.py
"""
import pandas as pd

# --- 1. load (edit the path to wherever your file is) ---
df = pd.read_csv("reef_fish.csv")
L = df["length_mm"]

names = ["small", "medium", "large", "xlarge"]

# --- 2. data-driven cut-offs = the length quartiles ---
# qcut splits into 4 groups of ~equal size using the 25/50/75% points.
df["size_class"] = pd.qcut(L, q=4, labels=names)

# show the actual length boundaries it used (rounded for reporting)
edges = L.quantile([0, .25, .50, .75, 1.0]).round(1).tolist()
print("length quartile cut-offs (mm):", edges)
#   e.g. [89.3, 111.9, 121.5, 128.5, 156.6]
#   small  = 89-112 | medium = 112-121 | large = 121-129 | xlarge = 129-157

# --- 3. counts and cumulative counts (now never NaN) ---
vc = df["size_class"].value_counts().reindex(names)
print("\ncounts per class:")
print(vc.to_string())
print("\ncumulative counts:")


print(vc.cumsum().to_string())

# --- 4. (optional) a quick check that everything was classified ---
print("\ntotal classified:", int(vc.sum()), "of", len(df))

# ----------------------------------------------------------------------
# ALTERNATIVE: fixed, biologically interpretable cut-offs (whole mm).
# Unbalanced for this dataset, but easier to explain in a methods section.
# import numpy as np
# df["size_class"] = pd.cut(
#     L, bins=[-np.inf, 100, 120, 140, np.inf], labels=names)
# df["size_class"].value_counts().reindex(names).cumsum()
