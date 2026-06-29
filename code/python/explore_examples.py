"""
Chapter 6 - Exploring Your Data (Python)
Summary statistics, frequency tables, and quick checks on a reef-fish dataset.

    python explore_examples.py

Needs: pandas, numpy   (scipy optional, for trim_mean / gmean / shapiro)
"""
import pandas as pd, numpy as np

fish = pd.read_csv("reef_fish.csv")

# ===================================================================
# A FIRST LOOK
# ===================================================================
print("shape:", fish.shape)                          # (300, 7)
print(fish.describe().round(1))

# ===================================================================
# 1. MEASURES OF LOCATION
# ===================================================================
M, P = fish["mass_g"], fish["parasite_count"]
def gmean(x): return float(np.exp(np.mean(np.log(x))))
def tmean(x, p=0.1):
    x = np.sort(np.asarray(x, float)); k = int(len(x)*p)
    return float(np.mean(x[k:len(x)-k]))

print("\n--- LOCATION ---")
print("length mean/median:", round(fish["length_mm"].mean(),1), fish["length_mm"].median())
print("mass mean/median/gmean/trim:", round(M.mean(),1), M.median(),
      round(gmean(M),1), round(tmean(M),1))          # 32.4 31.4 30.7 31.8
print("parasite mean/median/gmean:", round(P.mean(),1), P.median(), round(gmean(P),1))  # 6.7 5.0 4.5
print("mode parasite / habitat:", P.mode()[0], fish["habitat"].mode()[0])               # 1, Reef
print("weighted mean length (by mass):", round(np.average(fish["length_mm"], weights=M),1))  # 123.9
print("harmonic mean parasite:", round(len(P)/(1/P).sum(),2))                            # 2.96
print("mean length by habitat:\n", fish.groupby("habitat")["length_mm"].mean().round(1))

# ===================================================================
# 2. MEASURES OF SPREAD
# ===================================================================
L = fish["length_mm"]
print("\n--- SPREAD ---")
print("range:", round(L.min(),1), "-", round(L.max(),1), "=", round(L.max()-L.min(),1))
print("var/sd:", round(L.var(),1), round(L.std(),1))
print("quartiles:", L.quantile([.25,.5,.75]).round(1).tolist(), "IQR:", round(L.quantile(.75)-L.quantile(.25),1))
print("5th/95th pct:", round(L.quantile(.05),1), round(L.quantile(.95),1))
print("MAD (scaled):", round((L-L.median()).abs().median()*1.4826,1))
print("SEM:", round(L.std()/np.sqrt(len(L)),2))
print("CV length/mass %:", round(100*L.std()/L.mean(),1), round(100*M.std()/M.mean(),1))
print("skew length/mass/parasite:", round(L.skew(),2), round(M.skew(),2), round(P.skew(),2))
print("kurt length/parasite:", round(L.kurt(),2), round(P.kurt(),2))
print("parasite skew raw vs log:", round(P.skew(),2), round(np.log(P).skew(),2))         # 2.88 -> -0.15

# ===================================================================
# 3. COUNTS AND FREQUENCIES
# ===================================================================
print("\n--- COUNTS ---")
print(fish["habitat"].value_counts())
print((fish["habitat"].value_counts(normalize=True)*100).round(1))
print("crosstab sex x habitat:\n", pd.crosstab(fish["sex"], fish["habitat"]))
print("female % per habitat:\n",
      (pd.crosstab(fish["habitat"], fish["sex"], normalize="index")["F"]*100).round(1))
classes = pd.cut(fish["length_mm"], [0,100,120,140,1000],
                 labels=["small","medium","large","xlarge"])
print("size classes:\n", classes.value_counts().reindex(["small","medium","large","xlarge"]))
# Shannon diversity of habitats
p = fish["habitat"].value_counts(normalize=True); H = -(p*np.log(p)).sum()
print("Shannon H / evenness:", round(H,3), round(H/np.log(len(p)),3))                    # 1.684 0.94
# group profile table
print("habitat profile:\n", fish.groupby("habitat").agg(
    n=("specimen_id","size"), mean_length=("length_mm","mean"),
    mean_mass=("mass_g","mean"), median_parasites=("parasite_count","median")).round(1))

# ===================================================================
# 4. QUICK CHECKS BEFORE ANALYSIS
# ===================================================================
print("\n--- QUICK CHECKS ---")
print("missing:", int(fish.isna().sum().sum()), "duplicates:", int(fish.duplicated().sum()))
print("group sizes (min,max):", fish["habitat"].value_counts().min(), fish["habitat"].value_counts().max())
print("correlations:\n", fish.select_dtypes("number").corr().round(2))
# z-score outliers and a first effect size
z = (L - L.mean())/L.std()
print("|z|>2 lengths:", int((z.abs()>2).sum()))
m = fish[fish.sex=="M"]["mass_g"]; f = fish[fish.sex=="F"]["mass_g"]
sp = np.sqrt(((len(m)-1)*m.var() + (len(f)-1)*f.var())/(len(m)+len(f)-2))
print("Cohen's d mass M vs F:", round((m.mean()-f.mean())/sp, 2))                        # 0.08
