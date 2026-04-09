import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

path_shock = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data\velshock_index.dat" 
path_lum_photo = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data\lum_photo.dat"
path_lum_observed = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data\lum_observed.dat"

try:
    df_shock = pd.read_csv(path_shock, delim_whitespace=True, names=["grid","time", "radius", "velocity", "unused", "optical_depth"], comment="#")

    df_lum_photo = pd.read_csv(path_lum_photo, delim_whitespace=True, names=["time", "lum_photo"], comment="#")
    df_lum_observed = pd.read_csv(path_lum_observed, delim_whitespace=True, names=["time", "lum_observed"], comment="#")

    t_days_shock = df_shock["time"] / 86400.0
    v_shock = df_shock["velocity"]

    t_days_lum_photo = df_lum_photo["time"] / 86400.0
    lum_photo = df_lum_photo["lum_photo"]

    t_days_lum_observed = df_lum_observed["time"] / 86400.0
    lum_observed = df_lum_observed["lum_observed"]

    plt.figure(figsize=(8, 6))
    # plt.plot(t_days_shock, v_shock, label=f"Shock Velosity", color="red", linewidth=2)
    plt.plot(t_days_lum_photo, lum_photo, label=f"Photo Luminosity", color="blue", linewidth=2)
    plt.plot(t_days_lum_observed, lum_observed, label=f"Observed Luminosity", color="green", linewidth=2)
    plt.xlabel("Time [days]")
    # plt.ylabel("Velocity [cm/s]")
    plt.ylabel("Luminosity [erg/s]")
    # plt.title("Shock Velosity Evolution")
    plt.title("Luminosity Evolution")
    plt.grid(True)
    plt.legend()
    # plt.yscale("log") # ログスケールの方が見やすい場合も
    plt.show()

except Exception as e:
    print(f"Error: {e}")