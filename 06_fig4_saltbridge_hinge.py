"""Figure 4A & 4B: Catalytic Salt-Bridge Stability and Hinge Flexibility Decoupling."""

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
import seaborn as sns

Path("figures").mkdir(exist_ok=True)
plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 11,
})

df = pd.read_csv("data/processed_metrics.csv")

# 統計検定 (WT vs D754A 触媒塩橋距離)
wt_sb = df[df["Variant"] == "WT"]["Catalytic_Salt_Bridge_Dist"].dropna()
mut_sb = df[df["Variant"] == "D754A"]["Catalytic_Salt_Bridge_Dist"].dropna()
stat, p_val = mannwhitneyu(wt_sb, mut_sb)
print(
    f"Mann-Whitney U Test on Salt Bridge: stat = {stat}, p-value = {p_val:.4f}"
)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

# ----------------- Panel A: Scatter -----------------
sns.scatterplot(
    data=df,
    x="Hinge_E656_CA_Z",
    y="Catalytic_Salt_Bridge_Dist",
    hue="Variant",
    palette={"WT": "#5c768d", "D754A": "#d95f02"},
    alpha=0.8,
    s=40,
    ax=axes[0],
)
axes[0].set_title(
    "A. Salt Bridge Stability vs. Hinge Flexibility", fontweight="bold"
)
axes[0].set_xlabel("Hinge E656 Cα Z-coordinate (Å)", fontweight="bold")
axes[0].set_ylabel("K609(NZ) - E626(OE) Min Distance (Å)", fontweight="bold")
axes[0].grid(True, linestyle="--", alpha=0.5)

# ----------------- Panel B: Violin Plot -----------------
sns.violinplot(
    data=df,
    x="Variant",
    y="Catalytic_Salt_Bridge_Dist",
    palette={"WT": "#5c768d", "D754A": "#d95f02"},
    inner="box",
    cut=0,
    ax=axes[1],
)
sns.stripplot(
    data=df,
    x="Variant",
    y="Catalytic_Salt_Bridge_Dist",
    color="black",
    alpha=0.35,
    jitter=0.2,
    size=3,
    ax=axes[1],
)
axes[1].set_title(
    "B. Distribution of Catalytic Salt-Bridge Distance", fontweight="bold"
)
axes[1].set_xlabel("Variant", fontweight="bold")
axes[1].set_ylabel("K609(NZ) - E626(OE) Min Distance (Å)", fontweight="bold")
axes[1].set_ylim(2.35, 3.65)
axes[1].grid(axis="y", linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig(
    "figures/Figure4_SaltBridge_Hinge.png", dpi=300, bbox_inches="tight"
)
plt.close()
print("Figure 4A & 4B を figures/ に出力しました。")