import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# パスの設定 (ご自身の環境に合わせて変更してください)
# 通常、SNECは "shock_radius.dat" と "photosphere.dat" (あるいは magnitues.datなど) を出力します
# ここでは output ディレクトリにある一般的なファイル名を想定します
path_shock = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_ej_sin\velshock_index.dat" 
path_photo = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_ej_sin\rad_photo.dat" # もしこのファイルがあれば

# ファイルが存在しない場合、magnitudes.dat などのカラムに含まれている場合があります
# とりあえず、shock_radius だけでも確認しましょう

try:
    # shock_radius.dat の読み込み (time, radius, velocity などの形式が多い)
    df_shock = pd.read_csv(path_shock, delim_whitespace=True, names=["grid","time", "radius", "velocity", "unused", "optical_depth"], comment="#")
    
    df_photo = pd.read_csv(path_photo, delim_whitespace=True, names=["time", "radius"], comment="#")

    # 時間を "日" に変換
    t_days_shock = df_shock["time"] / 86400.0
    r_shock = df_shock["radius"]

    t_days_photo = df_photo["time"] / 86400.0
    r_photo = df_photo["radius"]

    plt.figure(figsize=(8, 6))
    plt.plot(t_days_shock, r_shock, label=f"Shock Radius", color="red", linewidth=2)
    plt.plot(t_days_photo, r_photo, label=f"Photosphere Radius", color="blue", linewidth=2)

    # 光球半径のデータがあれば重ねる（もし radius.xg から計算する必要がある場合は複雑ですが、
    # 多くのバージョンのSNECは photosphere の情報を info.dat や専用ファイルに出します）
    # ここでは概念的な確認のため、衝撃波半径だけプロットし、
    # 「8日目付近で半径の挙動に異常がないか（＝衝撃波自体はスムーズか）」を確認します。
    
    plt.xlabel("Time [days]")
    plt.ylabel("Radius [cm]")
    plt.title("Shock / Photosphere Radius Evolution")
    plt.grid(True)
    plt.legend()
    plt.yscale("log") # ログスケールの方が見やすい場合も
    plt.show()
    
except Exception as e:
    print(f"Error: {e}")