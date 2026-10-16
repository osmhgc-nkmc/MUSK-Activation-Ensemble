"""Figure 1B: Ensemble-Size Sensitivity and Convergence Analysis."""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

Path("figures").mkdir(exist_ok=True)
plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 11,
    "axes.linewidth": 1.2,
})

df = pd.read_csv("data/processed_metrics.csv")
wt_rmsd = df[df["Variant"] == "WT"]["Aloop_RMSD"].values
d754a_rmsd = df[df["Variant"] == "D754A"]["Aloop_RMSD"].values

np.random.seed(42)
n_iterations = 1000
subset_sizes = [25, 50, 75, 100]
records = []
summary_stats = []

for var_name, data in [("WT", wt_rmsd), ("D754A", d754a_rmsd)]:
  for N in subset_sizes:
    sub_ratios = []
    for _ in range(n_iterations):
      # 復元抽出モンテカルロブートストラップ
      sample = np.random.choice(data, size=N, replace=True)
      ratio = np.mean(sample >= 5.0) * 100.0
      sub_ratios.append(ratio)
      records.append({
          "Variant": var_name,
          "Subset Sample Size": f"N = {N}",
          "Estimated Open State Population (%)": ratio,
      })
    summary_stats.append({
        "Variant": var_name,
        "Subset_N": N,
        "Mean": np.mean(sub_ratios),
        "SD": np.std(sub_ratios),
        "CI_2.5": np.percentile(sub_ratios, 2.5),
        "CI_97.5": np.percentile(sub_ratios, 97.5),
    })

pd.DataFrame(summary_stats).to_csv(
    "figures/Table_S2_Convergence_Stats.csv", index=False
)

# 描画
plt.figure(figsize=(8, 4.8))
sns.boxplot(
    data=pd.DataFrame(records),
    x="Subset Sample Size",
    y="Estimated Open State Population (%)",
    hue="Variant",
    palette={"WT": "#5c768d", "D754A": "#b22222"},
    width=0.55,
    fliersize=2.5,
    linewidth=1.2,
)

plt.title(
    "Ensemble-Size Sensitivity & Convergence Analysis (MuSK K-ESPRIT)",
    fontweight="bold",
    pad=12,
)
plt.xlabel("Subset Sample Size", fontweight="bold")
plt.ylabel("Estimated Open State Population (%)", fontweight="bold")
plt.ylim(-5, 105)
plt.legend(title="Variant", frameon=True, loc="upper right")
sns.despine(top=True, right=True)

plt.tight_layout()
plt.savefig(
    "figures/Figure1B_Convergence_Boxplot.png", dpi=300, bbox_inches="tight"
)
plt.close()
print("Figure 1B を figures/ に出力しました。")