import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sun_mag = 4.75
sun_lum = 3.828e33

data1 = pd.read_csv(
    r'C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11\lum_observed.dat', 
    delim_whitespace=True, 
    names=['time', 'L_obs']
)

data2 = pd.read_csv(
    r'C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_7d-11_no\lum_observed.dat', 
    delim_whitespace=True, 
    names=['time', 'L_obs']
)

data1['time_days'] = data1['time'] / (60 * 60 * 24)
data2['time_days'] = data2['time'] / (60 * 60 * 24)

data1['M_bol'] = sun_mag - 2.5 * np.log10(data1['L_obs'] / sun_lum)
data2['M_bol'] = sun_mag - 2.5 * np.log10(data2['L_obs'] / sun_lum)

fig, ax = plt.subplots(figsize=(10, 6))

# 黒色の線でプロット
ax.plot(data1['time_days'], data1['M_bol'], lw=1.5, color='black', label=r"$\delta\rho = 0.4$")
ax.plot(data2['time_days'], data2['M_bol'], lw=1.5, color='red', label=r"$\delta\rho = 0$")

ax.set_xlim(0, 39)
ax.set_ylim(-23, -18)
ax.invert_yaxis()

# 軸ラベルを大きく
ax.set_xlabel('Time (day)', fontsize=18)
ax.set_ylabel(r'Bolometric Magnitude ($M_{\rm bol}$)', fontsize=18)
ax.set_title('Supernova Light Curve', fontsize=20)

# 目盛りを内向きに、短くする
ax.tick_params(axis='both', which='both', direction='in', 
               length=4, width=1, labelsize=14,
               top=True, right=True)

# x軸とy軸の副目盛りを非表示
from matplotlib.ticker import NullLocator
ax.xaxis.set_minor_locator(NullLocator())
ax.yaxis.set_minor_locator(NullLocator())

# グリッドを非表示
ax.grid(False)

# 枠線の太さを設定
for spine in ax.spines.values():
    spine.set_linewidth(1)

ax.legend(fontsize=18, loc='lower right')
plt.tight_layout()
plt.show()