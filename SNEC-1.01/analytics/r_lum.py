from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# ===============================
# xg ファイルを time(文字列) ごとに読む
# ===============================

def load_xg_by_time_prefix(filename, value_name):
    data = {}
    current_time_prefix = None
    rows = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # ★ ここを修正 ★
            if "Time" in line:
                if current_time_prefix is not None and rows:
                    data[current_time_prefix] = pd.DataFrame(
                        rows, columns=["mass", value_name]
                    )

                t = float(line.split("=")[1])
                current_time_prefix = int(t)   # ← ★ここが変更点

                rows = []
                continue

            parts = line.split()
            if len(parts) == 2:
                rows.append([float(parts[0]), float(parts[1])])

    if current_time_prefix is not None and rows:
        data[current_time_prefix] = pd.DataFrame(
            rows, columns=["mass", value_name]
        )

    return data

# ===============================
# ファイル読み込み
# ===============================
# base_dir = Path(".")  # 必要に応じて変更

radius_data = load_xg_by_time_prefix(
    r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data\radius.xg", "radius"
)
lum_data = load_xg_by_time_prefix(
    r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data\lum.xg", "lum"
)

# print("radius times:", sorted(radius_data.keys())[:10])
# print("lum times   :", sorted(lum_data.keys())[:10])

# ===============================
# 使いたい time（上5桁）を指定
# ===============================
time_prefix = 2159000   # 例： "0", "17000", "34000", など

if time_prefix not in radius_data or time_prefix not in lum_data:
    raise ValueError(f"time prefix {time_prefix} not found")

df_r = radius_data[time_prefix]
df_lum = lum_data[time_prefix]

# ===============================
# mass でマージ（安全のため）
# ===============================
df = pd.merge(df_r, df_lum, on="mass", how="inner")

# ===============================
# プロット
# ===============================
plt.figure(figsize=(7, 5))
plt.plot(df["radius"], df["lum"], lw=2)

plt.yscale("log")
plt.xscale("log")
plt.xlabel("Radius [cm]")
plt.ylabel("Luminosity [erg/s]")
plt.title(f"Radius–Luminosity relation (Time = {time_prefix})")
plt.grid(True)

plt.tight_layout()
plt.show()
