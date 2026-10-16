"""Figure 2A & 2B: Bivariate Correlation and Gaussian KDE Conformational Landscape."""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

Path("figures").mkdir(exist_ok=True)
plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 11,
})

df = pd.read_csv("data/processed_metrics.csv")
df_wt = df[df["Variant"] == "WT"].copy()

# ----------------- Figure 2A: Scatter Plot -----------------
plt.figure(figsize=(7.5, 4.8))
scatter = plt.scatter(
    df_wt["Aloop_RMSD"],
    df_wt["RS3_RS2_Centroid_Dist"],
    c=df_wt["RS3_RS2_Centroid_Dist"],
    cmap="viridis_r",
    edgecolors="black",
    linewidths=0.6,
    s=45,
    zorder=3,
)

# 代表重心モデルの強調表示
rep_models = {
    30: ("#30 (Closed)", 0.00, 9.98),
    45: ("#45 (Intermediate)", 3.14, 6.77),
    37: ("#37 (Open)", 12.54, 6.71),
}
for mid, (label, x_val, y_val) in rep_models.items():
  plt.scatter(
      x_val,
      y_val,
      s=120,
      facecolors="none",
      edgecolors="red",
      linewidths=2.0,
      zorder=4,
  )

plt.axhline(
    6.0,
    color="red",
    linestyle="--",
    linewidth=1.5,
    label="R-spine Packing Threshold (~6.0 Å)",
)
plt.title(
    "MuSK WT: A-loop RMSD vs. RS3-RS2 Distance", fontweight="bold", pad=12
)
plt.xlabel("A-loop Cα-RMSD vs. 1LUF (Å)", fontweight="bold")
plt.ylabel("RS3(M630) - RS2(F744) Centroid Distance (Å)", fontweight="bold")
cbar = plt.colorbar(scatter)
cbar.set_label("RS3-RS2 Dist. (Å)", fontweight="bold")
plt.grid(True, linestyle="-", alpha=0.4)

plt.tight_layout()
plt.savefig(
    "figures/Figure2A_Bivariate_Scatter.png", dpi=300, bbox_inches="tight"
)
plt.close()

# ----------------- Figure 2B: 2D Gaussian KDE -----------------
plt.figure(figsize=(7.5, 4.8))
sns.kdeplot(
    data=df_wt,
    x="Aloop_RMSD",
    y="RS3_RS2_Centroid_Dist",
    cmap="Greens",
    fill=True,
    levels=15,
    thresh=0.01,
    bw_method="scott",
)
plt.scatter(
    df_wt["Aloop_RMSD"],
    df_wt["RS3_RS2_Centroid_Dist"],
    color="black",
    alpha=0.55,
    s=20,
    edgecolors="white",
    linewidths=0.3,
)

plt.axhline(
    6.0,
    color="red",
    linestyle="--",
    linewidth=1.5,
    label="R-spine Packing Threshold (~6.0 Å)",
)
plt.title("Conformational Density Landscape (KDE)", fontweight="bold", pad=12)
plt.xlabel("A-loop Cα-RMSD vs. 1LUF (Å)", fontweight="bold")
plt.ylabel("RS3(M630) - RS2(F744) Centroid Distance (Å)", fontweight="bold")
plt.legend(loc="upper right", frameon=True)
plt.grid(True, linestyle="-", alpha=0.4)

plt.tight_layout()
plt.savefig("figures/Figure2B_KDE_Landscape.png", dpi=300, bbox_inches="tight")
plt.close()
print("Figure 2A & 2B を figures/ に出力しました。")