"""Figure 1C: Hierarchical Agglomerative Clustering and Population Shifts."""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import chi2_contingency, fisher_exact

Path("figures").mkdir(exist_ok=True)
plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 11,
})

df = pd.read_csv("data/processed_metrics.csv")

# A-loop RMSDに基づくWard法クラスタリング (閾値 1.5 Å)
for var in ["WT", "D754A"]:
  rmsds = df[df["Variant"] == var]["Aloop_RMSD"].values.reshape(-1, 1)
  Z = linkage(pdist(rmsds), method="ward")
  clusters = fcluster(Z, t=1.5, criterion="distance")
  # RMSD平均値が小さいクラスターをClosed (1)、大きい方をOpen (2) と定義
  c_means = {c: np.mean(rmsds[clusters == c]) for c in np.unique(clusters)}
  closed_c = min(c_means, key=c_means.get)
  df.loc[df["Variant"] == var, "State"] = np.where(
      clusters == closed_c, "Closed", "Open"
  )

# 状態別集計
counts = pd.crosstab(df["Variant"], df["State"])
table = [
    [counts.loc["WT", "Closed"], counts.loc["WT", "Open"]],
    [counts.loc["D754A", "Closed"], counts.loc["D754A", "Open"]],
]

# 統計検定
odds, p_fisher = fisher_exact(table)
chi2, p_chi2, _, _ = chi2_contingency(table)
print(f"Fisher's Exact Test p-value: {p_fisher:.4e}, Odds ratio: {odds:.2f}")
print(f"Chi-square Test p-value: {p_chi2:.4e}, Chi2: {chi2:.2f}")

# Stacked Bar Plot 描画
wt_closed = (df[df["Variant"] == "WT"]["State"] == "Closed").mean() * 100
wt_open = 100.0 - wt_closed
mut_closed = (df[df["Variant"] == "D754A"]["State"] == "Closed").mean() * 100
mut_open = 100.0 - mut_closed

plt.figure(figsize=(6, 5))
x = ["Wild-Type (WT)", "D754A Mutant"]
p_closed = [wt_closed, mut_closed]
p_open = [wt_open, mut_open]

plt.bar(
    x,
    p_closed,
    label="Closed State (Autoinhibited)",
    color="#2b5c8f",
    edgecolor="black",
    width=0.55,
)
plt.bar(
    x,
    p_open,
    bottom=p_closed,
    label="Open State (Active-like)",
    color="#d95f02",
    edgecolor="black",
    width=0.55,
)

# パーセンテージのアノテーション
plt.text(
    0,
    wt_closed / 2,
    f"{wt_closed:.1f}%",
    ha="center",
    va="center",
    color="white",
    fontweight="bold",
)
plt.text(
    0,
    wt_closed + wt_open / 2,
    f"{wt_open:.1f}%",
    ha="center",
    va="center",
    color="white",
    fontweight="bold",
)
plt.text(
    1,
    mut_closed / 2,
    f"{mut_closed:.1f}%",
    ha="center",
    va="center",
    color="white",
    fontweight="bold",
)
plt.text(
    1,
    mut_closed + mut_open / 2,
    f"{mut_open:.1f}%",
    ha="center",
    va="center",
    color="white",
    fontweight="bold",
)

plt.title(
    "Conformational State Populations via Hierarchical Clustering\n(Ward's"
    " Linkage, Threshold = 1.5 Å)",
    fontsize=11,
    pad=12,
)
plt.ylabel("Population Fraction (%)", fontweight="bold")
plt.ylim(0, 105)
plt.legend(frameon=True, loc="upper right")
plt.grid(axis="y", linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig(
    "figures/Figure1C_Population_Shift.png", dpi=300, bbox_inches="tight"
)
plt.close()
print("Figure 1C を figures/ に出力しました。")