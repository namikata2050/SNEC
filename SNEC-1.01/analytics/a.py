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

radius_data1 = load_xg_by_time_prefix(
    r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11_no\radius.xg", "radius"
)
rho_data1 = load_xg_by_time_prefix(
    r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11_no\rho.xg", "density"
)

radius_data2 = load_xg_by_time_prefix(
    r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11\radius.xg", "radius"
)
rho_data2 = load_xg_by_time_prefix(
    r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11\rho.xg", "density"
)

# print("radius times:", sorted(radius_data.keys())[:10])
# print("rho times   :", sorted(rho_data.keys())[:10])

# ===============================
# 使いたい time（上5桁）を指定
# ===============================
time_prefix = 0   # 例： "0", "17000", "34000", など

if time_prefix not in radius_data1 or time_prefix not in rho_data1:
    raise ValueError(f"time prefix {time_prefix} not found")

df_r1 = radius_data1[time_prefix]
df_rho1 = rho_data1[time_prefix]

if time_prefix not in radius_data2 or time_prefix not in rho_data2:
    raise ValueError(f"time prefix {time_prefix} not found")

df_r2 = radius_data2[time_prefix]
df_rho2 = rho_data2[time_prefix]

# ===============================
# mass でマージ（安全のため）
# ===============================
df1 = pd.merge(df_r1, df_rho1, on="mass", how="inner")
df2 = pd.merge(df_r2, df_rho2, on="mass", how="inner")

# ===============================
# プロット
# ===============================
plt.figure(figsize=(7, 5))
plt.plot(df1["radius"], df1["density"], lw=2)
plt.plot(df2["radius"], df2["density"], lw=2)


plt.yscale("log")
# plt.xscale("log")
plt.xlabel("Radius [cm]")
plt.ylabel("Density [g/cm^3]")
plt.title(f"Radius–Density relation (Time = {time_prefix})")
plt.grid(True)

plt.tight_layout()
plt.show()
