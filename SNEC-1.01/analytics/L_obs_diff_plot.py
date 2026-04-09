import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

sun_mag = 4.75
sun_lum = 3.828e33

data1 = pd.read_csv(
    r'C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_1d-9\lum_observed.dat', 
    delim_whitespace=True, 
    names=['time', 'L_obs']
)

data2 = pd.read_csv(
    r'C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data_1d-9_no\lum_observed.dat', 
    delim_whitespace=True, 
    names=['time', 'L_obs']
)

data1['time_days'] = data1['time'] / 86400.0
data2['time_days'] = data2['time'] / 86400.0

data1['M_bol'] = sun_mag - 2.5 * np.log10(data1['L_obs'] / sun_lum)
data2['M_bol'] = sun_mag - 2.5 * np.log10(data2['L_obs'] / sun_lum)

max_days = 33.0
min_days = 0.0

data1 = data1[(data1['time_days'] <= max_days) & (data1['time_days'] >= min_days)]
data2 = data2[(data2['time_days'] <= max_days) & (data2['time_days'] >= min_days)]

M_bol_data2_interp = np.interp(data1['time_days'], data2['time_days'], data2['M_bol'])

diff_M_bol = data1['M_bol'] - M_bol_data2_interp

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(data1['time_days'], diff_M_bol, label=r'$\Delta M_{\rm bol}$', color='black', linewidth=1.5)
    
ax.axhline(0, color='gray', linestyle='--', linewidth=1)

ax.set_xlim(min_days, max_days)
ax.set_ylim(-0.5, 0.5)

# 軸ラベルを大きく
ax.set_xlabel("Time (day)", fontsize=18)
ax.set_ylabel(r"Bolometric Magnitude (mag)", fontsize=18)
ax.set_title("Difference in Bolometric Magnitude", fontsize=20)

# 目盛りを内向きに、短くする
ax.tick_params(axis='both', which='both', direction='in', 
               length=4, width=1, labelsize=14,
               top=True, right=True)

# 副目盛りを非表示
from matplotlib.ticker import NullLocator
ax.xaxis.set_minor_locator(NullLocator())
ax.yaxis.set_minor_locator(NullLocator())

# グリッドを非表示
ax.grid(False)

# 枠線の太さを設定
for spine in ax.spines.values():
    spine.set_linewidth(1)

ax.legend(fontsize=14, loc='lower right')

plt.tight_layout()
plt.show()