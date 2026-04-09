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
    r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11\radius.xg", "radius"
)
vel_data = load_xg_by_time_prefix(
    r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11\vel.xg", "velocity"
)

# print("radius times:", sorted(radius_data.keys())[:10])
# print("vel times   :", sorted(vel_data.keys())[:10])

# ===============================
# 使いたい time（上5桁）を指定
# ===============================
time_prefix = 0   # 例： 0, 17000, 34000, など

if time_prefix not in radius_data or time_prefix not in vel_data:
    raise ValueError(f"time prefix {time_prefix} not found")

df_r = radius_data[time_prefix]
df_v = vel_data[time_prefix]

# ===============================
# mass でマージ（安全のため）
# ===============================
df = pd.merge(df_r, df_v, on="mass", how="inner")

# ===============================
# プロット
# ===============================
fig, ax = plt.subplots(figsize=(7, 5))

# 黒色の線でプロット
ax.plot(df["radius"], df["velocity"], lw=1.5, color='black')

# r = 2*10^13 の位置に赤い点線を追加
ax.axvline(x=2e13, color='red', linestyle='--', linewidth=1.5)

ax.set_yscale("log")
ax.set_xscale("log")

# 軸ラベルを大きく
ax.set_xlabel("Radius (cm)", fontsize=18)
ax.set_ylabel(r"Velocity (cm s$^{-1}$)", fontsize=18)
ax.set_title(f"Radius–Velocity relation (Time = {time_prefix})", fontsize=20)

# 目盛りを内向きに、短くする
ax.tick_params(axis='both', which='both', direction='in', 
               length=4, width=1, labelsize=14,
               top=True, right=True)

# y軸の目盛りを設定
from matplotlib.ticker import FixedLocator, LogFormatterMathtext, NullFormatter, NullLocator

# 主目盛り（ラベル付き）：各桁で表示
ax.yaxis.set_major_locator(FixedLocator([1e6, 1e7, 1e8, 1e9, 1e10]))
ax.yaxis.set_major_formatter(LogFormatterMathtext())

# 副目盛り（ラベルなし）：各桁で表示
minor_ticks_y = [10**(i) for i in range(5, 12)]
ax.yaxis.set_minor_locator(FixedLocator(minor_ticks_y))
ax.yaxis.set_minor_formatter(NullFormatter())

# 副目盛りのスタイル（短め）
ax.tick_params(axis='y', which='minor', direction='in', length=2, width=0.5)

# x軸の副目盛りを非表示
ax.xaxis.set_minor_locator(NullLocator())

# グリッドを非表示
ax.grid(False)

# 枠線の太さを設定
for spine in ax.spines.values():
    spine.set_linewidth(1)

plt.tight_layout()
plt.show()

