import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sun_mag = 4.75
sun_lum = 3.828e33

path = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11\lum_observed.dat"

try:
    df = pd.read_csv(path, delim_whitespace=True, names=["time", "L_obs"], comment="#")

    # 時間を "日" に変換
    t_days = df["time"] / 86400.0

    df['M_bol'] = sun_mag - 2.5 * np.log10(df['L_obs'] / sun_lum)

    plt.figure(figsize=(8, 6))
    plt.plot(t_days, df['M_bol'], label=f"Observed Luminosity", linewidth=2)

    plt.xlabel("Time [days]")
    plt.ylabel("M_bol")
    # plt.xlim(0, 16)
    plt.ylim(-23, -19)
    plt.gca().invert_yaxis()
    plt.xlabel('Time (days)')
    plt.ylabel('Bolometric Magnitude (M_bol)')
    plt.title('SNEC Supernova Light Curve')
    plt.grid(True)
    plt.legend()
    plt.show()
    
except Exception as e:
    print(f"Error: {e}")