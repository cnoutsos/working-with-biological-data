"""
Chapter 7 - Telling a Story with Figures (Python)
Ready to run on NumPy 2.x. One figure per chart type, colourblind-safe and labelled.
Saves PNGs into a 'figures/' folder next to this script.

    python Ch7_plot_examples.py

If your environment has the classic NumPy 1.x/2.x clash (pandas/matplotlib refuse
to import because pyarrow/numexpr/bottleneck were built against NumPy 1.x), this
script DETECTS it on startup, upgrades only those packages to their NumPy-2 builds
(NumPy itself stays 2.x), and relaunches itself once. No manual fixing needed.
"""

# =====================================================================
# 0. Self-heal the NumPy 2.0 ABI clash, then relaunch once if needed.
# =====================================================================
import os, sys, subprocess

# Every compiled package that pandas/matplotlib pull in and that may have been
# built against NumPy 1.x in an older Anaconda install:
_NUMPY2_REBUILD = ["pandas", "matplotlib", "pyarrow", "numexpr", "bottleneck"]

def _heavy_imports_ok():
    try:
        import pandas              # noqa: F401
        import matplotlib.pyplot   # noqa: F401
        return True, None
    except Exception as err:
        return False, err

def _bootstrap_numpy2():
    ok, err = _heavy_imports_ok()
    if ok:
        return                       # everything imports -> nothing to do

    if os.environ.get("CH7_REPAIRED") == "1":
        # We already upgraded once and it STILL fails -> stop the whack-a-mole
        # and tell the user the one fix that always works in a conda base env.
        print("\n" + "="*70)
        print("Still failing after upgrading packages.")
        print("Your Anaconda base env has more modules built for NumPy 1.x than")
        print("can be patched piecemeal. The GUARANTEED one-line fix is to pin")
        print("NumPy below 2 so it matches the rest of your conda stack:\n")
        print('    conda install "numpy<2"')
        print("\nThen rerun:  python Ch7_plot_examples.py")
        print("\n(Or make a clean env that uses NumPy 2 throughout:")
        print('    conda create -n biobook -c conda-forge python=3.12 "numpy>=2" \\')
        print("            pandas matplotlib")
        print("    conda activate biobook")
        print("    python Ch7_plot_examples.py )")
        print("="*70)
        sys.exit(1)

    text = repr(err)
    if not any(s in text for s in
               ("NumPy", "numpy", "_ARRAY_API", "multiarray", "initialization failed")):
        raise err                    # an unrelated import error -> show it

    print("NumPy 1.x/2.x ABI clash detected — upgrading", ", ".join(_NUMPY2_REBUILD),
          "to their NumPy-2 builds (NumPy itself stays 2.x)...", flush=True)
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-U", *_NUMPY2_REBUILD])
    env = dict(os.environ, CH7_REPAIRED="1")
    sys.exit(subprocess.call([sys.executable, os.path.abspath(__file__)], env=env))

_bootstrap_numpy2()

# =====================================================================
# 1. Imports + robust paths (find data & save figures next to this file)
# =====================================================================
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                     # safe headless backend; works everywhere
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
def here(name): return os.path.join(HERE, name)
os.makedirs(here("figures"), exist_ok=True)

# =====================================================================
# 2. Make the data files if they aren't already present (so it just runs)
# =====================================================================
def ensure_data():
    if not os.path.exists(here("reef_fish.csv")):
        rng = np.random.default_rng(6); N = 300
        habs = ["Reef","Seagrass","Mangrove","Estuary","Open Water","Rocky Shore"]
        hw = [.30,.20,.20,.15,.10,.05]
        length = rng.normal(120, 12, N).round(1)
        mass = (0.000018*length**3*np.exp(rng.normal(0,0.12,N))).round(1)
        parasites = np.clip(np.round(np.exp(rng.normal(1.5,1.0,N))),1,None).astype(int)
        expression = rng.normal(120, 18, N).round(1)
        pd.DataFrame({"specimen_id":[f"F{i:03d}" for i in range(1,N+1)],
                      "length_mm":length, "mass_g":mass, "parasite_count":parasites,
                      "expression":expression, "sex":rng.choice(["M","F"],N),
                      "habitat":rng.choice(habs,N,p=hw)}).to_csv(here("reef_fish.csv"), index=False)
    if not os.path.exists(here("abundance.csv")):
        rng = np.random.default_rng(11)
        months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
        reef = (40+18*np.sin((np.arange(12)-3)/12*2*np.pi)+rng.normal(0,3,12)).round().astype(int)
        seagrass = (30+10*np.sin((np.arange(12)-5)/12*2*np.pi)+rng.normal(0,3,12)).round().astype(int)
        pd.DataFrame({"month":months,"Reef":reef,"Seagrass":seagrass}).to_csv(here("abundance.csv"), index=False)

ensure_data()
fish = pd.read_csv(here("reef_fish.csv"))

# Okabe-Ito colourblind-safe palette
BLUE, GREEN, ORANGE, SKY, VERM = "#0072B2", "#009E73", "#E69F00", "#56B4E9", "#D55E00"
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                     "font.size": 12, "axes.titleweight": "bold"})

# =====================================================================
# 3. DISTRIBUTION -> histogram
# =====================================================================
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(fish["length_mm"], bins=20, color=BLUE, edgecolor="white")
ax.axvline(fish["length_mm"].median(), ls="--", color=VERM, lw=2,
           label=f"median {fish['length_mm'].median():.0f} mm")
ax.set_xlabel("Standard length (mm)"); ax.set_ylabel("Number of fish")
ax.set_title("Distribution of body length"); ax.legend(frameon=False)
fig.savefig(here("figures/histogram.png"), dpi=300, bbox_inches="tight"); plt.close(fig)

# =====================================================================
# 4. COMPARE COUNTS -> bar chart (bars start at zero!)
# =====================================================================
fig, ax = plt.subplots(figsize=(6, 4))
vc = fish["habitat"].value_counts()
ax.bar(vc.index, vc.values, color=GREEN)
ax.set_ylim(0, vc.max() * 1.12)
ax.set_xlabel("Habitat"); ax.set_ylabel("Number of fish")
ax.set_title("Sample size per habitat")
plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
fig.savefig(here("figures/bar.png"), dpi=300, bbox_inches="tight"); plt.close(fig)

# =====================================================================
# 5. COMPARE DISTRIBUTIONS -> boxplot   (tick_labels: matplotlib >= 3.9)
# =====================================================================
fig, ax = plt.subplots(figsize=(6.5, 4))
order = fish.groupby("habitat")["mass_g"].median().sort_values().index
ax.boxplot([fish[fish.habitat == h]["mass_g"] for h in order], tick_labels=list(order),
           patch_artist=True, boxprops=dict(facecolor=SKY),
           medianprops=dict(color="black"))
ax.set_xlabel("Habitat"); ax.set_ylabel("Body mass (g)")
ax.set_title("Body mass by habitat")
plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
fig.savefig(here("figures/boxplot.png"), dpi=300, bbox_inches="tight"); plt.close(fig)

# =====================================================================
# 6. RELATIONSHIP -> scatterplot (colour AND shape = accessible)
# =====================================================================
fig, ax = plt.subplots(figsize=(6, 4.2))
for sex, col, mk in [("F", ORANGE, "o"), ("M", BLUE, "^")]:
    sub = fish[fish.sex == sex]
    ax.scatter(sub["length_mm"], sub["mass_g"], c=col, marker=mk, alpha=0.75, label=sex)
ax.set_xlabel("Standard length (mm)"); ax.set_ylabel("Body mass (g)")
ax.set_title("Mass increases with length"); ax.legend(title="Sex", frameon=False)
fig.savefig(here("figures/scatter.png"), dpi=300, bbox_inches="tight"); plt.close(fig)

# =====================================================================
# 7. TREND OVER TIME -> line chart
# =====================================================================
ab = pd.read_csv(here("abundance.csv"))
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(ab["month"], ab["Reef"], marker="o", color=VERM, lw=2, label="Reef")
ax.plot(ab["month"], ab["Seagrass"], marker="s", color=BLUE, lw=2, label="Seagrass")
ax.set_ylim(0, None)
ax.set_xlabel("Month"); ax.set_ylabel("Mean catch per trip")
ax.set_title("Seasonal abundance by habitat"); ax.legend(frameon=False)
fig.savefig(here("figures/line.png"), dpi=300, bbox_inches="tight"); plt.close(fig)

# =====================================================================
# 8. MANY GROUPS -> small multiples (faceting)
# =====================================================================
habs = list(fish["habitat"].unique())
ncol = 3; nrow = int(np.ceil(len(habs) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(8, 2.6*nrow), sharex=True, sharey=True)
axes = np.atleast_1d(axes).ravel()
for ax, h in zip(axes, habs):
    sub = fish[fish.habitat == h]
    ax.scatter(sub["length_mm"], sub["mass_g"], s=12, color=BLUE, alpha=0.7)
    ax.set_title(h, fontsize=11)
for ax in axes[len(habs):]:            # hide any unused panels
    ax.set_visible(False)
fig.supxlabel("Length (mm)"); fig.supylabel("Mass (g)")
fig.tight_layout()
fig.savefig(here("figures/facets.png"), dpi=300, bbox_inches="tight"); plt.close(fig)

print(f"Saved 6 figures to {here('figures')}  — all colourblind-safe, labelled, vector-ready.")
print("NumPy", np.__version__, "| pandas", pd.__version__, "| matplotlib", matplotlib.__version__)
