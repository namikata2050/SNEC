import matplotlib.pyplot as plt
import pandas as pd
import io
import numpy as np

filename = r"C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01_sub1\profiles\simple1d_plus_csm.short"

df = pd.read_csv(
    filename,
    delim_whitespace=True,
    skiprows=1,
    names=['i', 'mass(i)', 'rad(i)', 'temp(i)', 'rho(i)', 'vel(i)', 'ye(i)', 'omega(i)']
)

# ==========================================
# プロット
# ==========================================
plt.figure(figsize=(10, 6))

# 密度プロファイルを描画

plt.plot(df['rad(i)'], df['rho(i)'], color='blue', linewidth=1.5)
# plt.plot(df['rad(i)'], df['vel(i)'], color='blue', linewidth=1.5)


# 軸の設定（密度は対数軸）
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'log Radius (cm)', fontsize=14)
plt.ylabel(r'log Density (g/cm$^3$)', fontsize=14)
#plt.ylabel(r'log velocity (cm/s)', fontsize=14)
#plt.title('Initial Density Profile (Ejecta + CSM)', fontsize=16)

# グリッドと凡例
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()

# 表示
plt.tight_layout()
plt.show()