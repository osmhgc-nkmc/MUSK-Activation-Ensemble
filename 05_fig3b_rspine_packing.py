"""Figure 3B: R-Spine Packing Distances Across Transition Representative Models."""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

Path("figures").mkdir(exist_ok=True)
plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 11,
})

# 代表3モデルの実測最短重原子間距離 (Å)
pairs = [
    "RS1–RS2 (723–744)",
    "RS2–RS3 (744–630)",
    "RS3–RS4 (630–641)",
]
d_closed = [7.04, 6.41, 4.13]  # Model #30
d_inter = [5.20, 3.40, 4.05]  # Model #45
d_open = [4.02, 3.08, 4.00]  # Model #37

x = np.arange(len(pairs))
width = 0.25

plt.figure(figsize=(8, 4.5))
rects1 = plt.bar(
    x - width,
    d_closed,
    width,
    label="Model #30 (Closed)",
    color="#3b528b",
    edgecolor="black",
)
rects2 = plt.bar(
    x,
    d_inter,
    width,
    label="Model #45 (Intermediate)",
    color="#21918c",
    edgecolor="black",
)
rects3 = plt.bar(
    x + width,
    d_open,
    width,
    label="Model #37 (Open)",
    color="#d95f02",
    edgecolor="black",
)

plt.axhline(
    4.5,
    color="gray",
    linestyle="--",
    linewidth=1.2,
    label="Packing Threshold (4.5 Å)",
)


def autolabel(rects):
  for rect in rects:
    height = rect.get_height()
    plt.annotate(
        f"{height:.2f}",
        xy=(rect.get_x() + rect.get_width() / 2, height),
        xytext=(0, 3),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=9,
    )


autolabel(rects1)
autolabel(rects2)
autolabel(rects3)

plt.title(
    "R-Spine Packing Across Transition States (#30, #45, #37)",
    fontsize=12,
    pad=12,
)
plt.xlabel("R-Spine Adjacent Pair", fontweight="bold")
plt.ylabel("Minimum Heavy-Atom Distance (Å)", fontweight="bold")
plt.xticks(x, pairs)
plt.ylim(0, 8.5)
plt.legend(title="Conformational State", frameon=True, loc="upper right")

plt.tight_layout()
plt.savefig("figures/Figure3B_RSpine_Bar.png", dpi=300, bbox_inches="tight")
plt.close()
print("Figure 3B を figures/ に出力しました。")