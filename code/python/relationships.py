"""
Chapter 12 - Relationships Between Variables (Python)
Correlation, linear/multiple regression, and non-linear regression, end to end.

    pip install pandas numpy scipy statsmodels
    python relationships.py

Data: reef_fish.csv, growth.csv, kinetics.csv, survival.csv.
#> lines show the verified fitted values.
"""
import pandas as pd, numpy as np
import scipy.stats as st
from scipy.optimize import curve_fit
import statsmodels.formula.api as smf

fish   = pd.read_csv("reef_fish.csv")
growth = pd.read_csv("growth.csv")
kin    = pd.read_csv("kinetics.csv")
surv   = pd.read_csv("survival.csv")

def hdr(t): print("\n" + "="*58 + "\n" + t + "\n" + "-"*58)

# ---------------------------------------------------- correlation
hdr("Correlation: length vs mass")
print("Pearson :", round(st.pearsonr(fish.length_mm, fish.mass_g)[0], 3))   #> 0.908
print("Spearman:", round(st.spearmanr(fish.length_mm, fish.mass_g)[0], 3))  #> 0.917
print("Kendall :", round(st.kendalltau(fish.length_mm, fish.mass_g)[0], 3)) #> 0.761
print("\ncorrelation matrix:\n",
      fish[['length_mm','mass_g','parasite_count','expression']].corr().round(2))

# ---------------------------------------------------- simple regression
hdr("Simple linear regression: mass ~ length")
m = smf.ols("mass_g ~ length_mm", data=fish).fit()
print(f"intercept={m.params['Intercept']:.1f}  slope={m.params['length_mm']:.3f}  R2={m.rsquared:.3f}")
#> -63.7, 0.799, 0.824

# ---------------------------------------------------- multiple regression
hdr("Multiple regression: mass ~ length + parasites + sex")
m2 = smf.ols("mass_g ~ length_mm + parasite_count + C(sex)", data=fish).fit()
print(m2.params.round(3).to_string())
print(f"R2={m2.rsquared:.3f}  adjR2={m2.rsquared_adj:.3f}  AIC={m2.aic:.0f}")

# ---------------------------------------------------- polynomial
hdr("Polynomial regression: density ~ time")
print("quadratic coeffs:", np.round(np.polyfit(growth.time_h, growth.density_OD, 2), 4))

# ---------------------------------------------------- exponential / power
hdr("Exponential & power (via log transforms)")
b_pow, _ = np.polyfit(np.log(fish.length_mm), np.log(fish.mass_g), 1)
print("allometric exponent (mass ~ length):", round(b_pow, 2))   #> 3.06

# ---------------------------------------------------- non-linear: logistic growth
hdr("Logistic growth (curve_fit)")
logi = lambda t, K, r, t0: K / (1 + np.exp(-r*(t - t0)))
p, _ = curve_fit(logi, growth.time_h, growth.density_OD, p0=[1, 0.5, 12])
print(f"K={p[0]:.3f}  r={p[1]:.3f}  t0={p[2]:.2f}")              #> 1.06, 0.46, 11.1

# ---------------------------------------------------- non-linear: Michaelis-Menten
hdr("Michaelis-Menten (curve_fit)")
mm = lambda S, Vmax, Km: Vmax*S/(Km + S)
p, _ = curve_fit(mm, kin.substrate_mM, kin.velocity, p0=[2, 1])
print(f"Vmax={p[0]:.3f}  Km={p[1]:.3f}")                        #> 2.18, 1.57

# ---------------------------------------------------- logistic regression
hdr("Logistic regression: survived ~ body_size")
lr = smf.logit("survived ~ body_size_mm", data=surv).fit(disp=False)
print(f"intercept={lr.params['Intercept']:.2f}  slope={lr.params['body_size_mm']:.3f}")
print(f"odds ratio per mm = {np.exp(lr.params['body_size_mm']):.2f}")

print("\nDone.")
