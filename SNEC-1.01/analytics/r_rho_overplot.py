from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# ===============================
# xg ファイルを time(文字列->数値) ごとに読む関数
# ===============================
def load_xg_by_time_prefix(filename, value_name):
    data = {}
    current_time_prefix = None
    rows = []

    try:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if "Time" in line:
                    # 前のブロックがあれば保存
                    if current_time_prefix is not None and rows:
                        data[current_time_prefix] = pd.DataFrame(
                            rows, columns=["mass", value_name]
                        )
                    
                    # Time = xxxx の行から時間を取得して整数化（キーにする）
                    try:
                        t_str = line.split("=")[1]
                        t_val = float(t_str)
                        current_time_prefix = int(t_val) 
                    except (IndexError, ValueError):
                        current_time_prefix = None

                    rows = []
                    continue

                # データ行の読み込み
                parts = line.split()
                if len(parts) == 2:
                    rows.append([float(parts[0]), float(parts[1])])

        # 最後のブロックを保存
        if current_time_prefix is not None and rows:
            data[current_time_prefix] = pd.DataFrame(
                rows, columns=["mass", value_name]
            )
            
    except FileNotFoundError:
        print(f"エラー: ファイルが見つかりません -> {filename}")
        return {}

    return data

# ===============================
# ファイルパスの設定
# ===============================
# ※ご自身の環境に合わせてパスを確認してください
path_radius = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11\radius.xg"
path_rho    = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11\rho.xg"
# path_lum    = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_cvisc10\lum.xg"

# ===============================
# データ読み込み
# ===============================
radius_data = load_xg_by_time_prefix(path_radius, "radius")
rho_data    = load_xg_by_time_prefix(path_rho, "density")
# lum_data    = load_xg_by_time_prefix(path_lum, "luminosity")

# データが空の場合は終了
if not radius_data or not rho_data:
    print("データの読み込みに失敗しました。パスを確認してください。")
    exit()
# if not radius_data or not lum_data:
#     print("データの読み込みに失敗しました。パスを確認してください。")
#     exit()

# ===============================
# ★ここに重ね書きしたい時刻(整数)をリストで指定★
# ===============================
# 例: 初期(0), 衝突前(680000), 衝突中(765000), 衝突後(969000) など
target_times = [0, 17006, 85008, 867000, 3366000]
# target_times = [0, 17005, 85012, 204000, 425000, 867000, 2414000, 4794000]
# target_times = [17000, 969000]

# ===============================
# プロット処理
# ===============================
fig, ax = plt.subplots(figsize=(10, 7))

# 指定した時刻をループしてプロット
for t in target_times:
    # データが存在するかチェック
    if t not in radius_data:
        print(f"Warning: Time {t} not found in radius data. Skipping.")
        continue
    if t not in rho_data:
        print(f"Warning: Time {t} not found in rho data. Skipping.")
        continue

    # データを取得
    df_r = radius_data[t]
    df_rho = rho_data[t]

    # mass をキーにしてマージ
    df = pd.merge(df_r, df_rho, on="mass", how="inner")

    # プロット (label引数で凡例用の名前を設定)
    ax.plot(df["radius"], df["density"], lw=1.5, label=f"Time = {t} s")

# ===============================
# グラフの装飾
# ===============================
# ax.set_ylim(1e2, 1e10)
ax.set_yscale("log")
ax.set_xscale("log")

# 軸ラベルを大きく
ax.set_xlabel("Radius (cm)", fontsize=18)
ax.set_ylabel(r"Density (g cm$^{-3}$)", fontsize=18)
# ax.set_ylabel(r"Velocity (cm s$^{-1}$)", fontsize=18)
ax.set_title("Radius–Density Evolution", fontsize=20)
# ax.set_title("Radius–Velocity Evolution", fontsize=20)

# 目盛りを内向きに、短くする
ax.tick_params(axis='both', which='both', direction='in', 
               length=4, width=1, labelsize=14,
               top=True, right=True)

# 副目盛りを非表示、y軸のラベルは偶数次数のみ
from matplotlib.ticker import NullLocator, FixedLocator, LogFormatterMathtext
ax.xaxis.set_minor_locator(NullLocator())
ax.yaxis.set_minor_locator(NullLocator())

# y軸の主目盛りを偶数次数のみに設定（10^2, 10^4, 10^6, 10^8, 10^10）
# ax.yaxis.set_major_locator(FixedLocator([1e2, 1e4, 1e6, 1e8, 1e10]))
# ax.yaxis.set_major_formatter(LogFormatterMathtext())

# y軸の目盛りを設定
from matplotlib.ticker import FixedLocator, LogFormatterMathtext, NullFormatter, NullLocator

# 主目盛り（ラベル付き）：4つ飛ばしで表示（10^-5, 10^-9, 10^-13, 10^-17）
ax.yaxis.set_major_locator(FixedLocator([1e-5, 1e-9, 1e-13, 1e-17]))
ax.yaxis.set_major_formatter(LogFormatterMathtext())

# 副目盛り（ラベルなし）：各桁で表示（10^-5 ～ 10^-22）
minor_ticks = [10**(-i) for i in range(5, 23)]
ax.yaxis.set_minor_locator(FixedLocator(minor_ticks))
ax.yaxis.set_minor_formatter(NullFormatter())

# 副目盛りのスタイル（短め）
ax.tick_params(axis='y', which='minor', direction='in', length=2, width=0.5)

# グリッドを非表示
ax.grid(False)

# 枠線の太さを設定
for spine in ax.spines.values():
    spine.set_linewidth(1)

# 凡例を左下に配置
ax.legend(fontsize=14, loc='lower left')

plt.tight_layout()
plt.show()