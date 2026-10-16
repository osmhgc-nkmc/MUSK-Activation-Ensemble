"""MuSK Kinase Structural Metrics Extraction Pipeline.

Extracts conformational coordinates from AlphaFold3 ensemble CIF files.
"""

from pathlib import Path
from Bio.PDB import MMCIFParser, Superimposer
import numpy as np
import pandas as pd

# 定義パラメータ
A_LOOP_RES = list(range(743, 766))  # 活性化ループ (743-765)
CORE_RES = list(range(570, 605)) + list(range(640, 720))  # 不変コア残基
RS1_NUM, RS2_NUM, RS3_NUM, RS4_NUM = 723, 744, 630, 641  # R-spine
K_RES, E_RES = 609, 626  # 触媒塩橋
HINGE_RES = 656  # ヒンジ残基


def get_ca_atoms(structure, res_nums):
  atoms = []
  for model in structure:
    for chain in model:
      for res in chain:
        if res.id[1] in res_nums and "CA" in res:
          atoms.append(res["CA"])
  return atoms


def get_residue(structure, res_num):
  for model in structure:
    for chain in model:
      for res in chain:
        if res.id[1] == res_num:
          return res
  return None


def min_heavy_distance(res1, res2):
  if res1 is None or res2 is None:
    return np.nan
  atoms1 = [
      a
      for a in res1.get_atoms()
      if not a.element.startswith("H") and a.name != "H"
  ]
  atoms2 = [
      a
      for a in res2.get_atoms()
      if not a.element.startswith("H") and a.name != "H"
  ]
  min_d = np.inf
  for a1 in atoms1:
    for a2 in atoms2:
      d = np.linalg.norm(a1.coord - a2.coord)
      if d < min_d:
        min_d = d
  return min_d


def calc_centroid(res):
  atoms = [
      a
      for a in res.get_atoms()
      if not a.element.startswith("H") and a.name not in ["N", "CA", "C", "O"]
  ]
  if not atoms:
    atoms = [a for a in res.get_atoms() if not a.element.startswith("H")]
  coords = np.array([a.coord for a in atoms])
  return np.mean(coords, axis=0)


def process_cif_files(cif_dir, ref_struct, variant_label):
  parser = MMCIFParser(QUIET=True)
  ref_core = get_ca_atoms(ref_struct, CORE_RES)
  ref_aloop = get_ca_atoms(ref_struct, A_LOOP_RES)

  records = []
  cif_paths = sorted(list(Path(cif_dir).glob("*.cif")))

  for idx, p in enumerate(cif_paths):
    struct = parser.get_structure(f"{variant_label}_{idx+1}", str(p))
    mobile_core = get_ca_atoms(struct, CORE_RES)

    # 1. 不変コア骨格への剛体重ね合わせ
    sup = Superimposer()
    sup.set_atoms(ref_core, mobile_core)
    sup.apply(struct.get_atoms())

    # 2. A-loop Cα-RMSD
    mobile_aloop = get_ca_atoms(struct, A_LOOP_RES)
    aloop_rmsd = np.sqrt(
        np.mean([
            np.sum((r.coord - m.coord) ** 2)
            for r, m in zip(ref_aloop, mobile_aloop)
        ])
    )

    # 3. RS3-RS2 重心間距離
    rs3 = get_residue(struct, RS3_NUM)
    rs2 = get_residue(struct, RS2_NUM)
    rs3_rs2_dist = np.linalg.norm(calc_centroid(rs3) - calc_centroid(rs2))

    # 4. R-spine 最短原子間距離
    rs1 = get_residue(struct, RS1_NUM)
    rs4 = get_residue(struct, RS4_NUM)
    d_rs1_rs2 = min_heavy_distance(rs1, rs2)
    d_rs2_rs3 = min_heavy_distance(rs2, rs3)
    d_rs3_rs4 = min_heavy_distance(rs3, rs4)

    # 5. 触媒塩橋 K609(NZ) - E626(OE)
    k_res = get_residue(struct, K_RES)
    e_res = get_residue(struct, E_RES)
    d_sb = np.nan
    if k_res and e_res and "NZ" in k_res:
      nz_coord = k_res["NZ"].coord
      oe_atoms = [
          e_res[a].coord for a in ["OE1", "OE2"] if a in e_res.child_dict
      ]
      if oe_atoms:
        d_sb = min([np.linalg.norm(nz_coord - c) for c in oe_atoms])

    # 6. ヒンジ E656 Cα Z座標
    hinge_res = get_residue(struct, HINGE_RES)
    hinge_z = hinge_res["CA"].coord[2] if (hinge_res and "CA" in hinge_res) else np.nan

    records.append({
        "Model_ID": idx + 1,
        "Variant": variant_label,
        "Aloop_RMSD": aloop_rmsd,
        "RS3_RS2_Centroid_Dist": rs3_rs2_dist,
        "Dist_RS1_RS2": d_rs1_rs2,
        "Dist_RS2_RS3": d_rs2_rs3,
        "Dist_RS3_RS4": d_rs3_rs4,
        "Catalytic_Salt_Bridge_Dist": d_sb,
        "Hinge_E656_CA_Z": hinge_z,
    })
  return records


if __name__ == "__main__":
  parser = MMCIFParser(QUIET=True)
  # 基準モデル: WT Model #30
  ref = parser.get_structure("WT_30", "data/raw_cif/WT_model_30.cif")

  data_wt = process_cif_files("data/raw_cif/WT/", ref, "WT")
  data_d754a = process_cif_files("data/raw_cif/D754A/", ref, "D754A")

  df_all = pd.DataFrame(data_wt + data_d754a)
  df_all.to_csv("data/processed_metrics.csv", index=False)
  print(f"特徴量抽出完了: {len(df_all)} モデルを data/processed_metrics.csv に保存しました。")