import pandas as pd
import matplotlib.pyplot as plt

# スペース区切りのデータを読み込む
# (カラム名はSNECの出力に合わせて確認してください)
data = pd.read_csv(
    r'C:\Users\hiroto yamada\OneDrive - Kyoto University\ドキュメント\Programming\Fortran\SNEC-1.01\Data\magnitudes.dat', 
    delim_whitespace=True, 
    # ↓↓↓ この行を修正しました (13個の列名を指定) ↓↓↓
    names=['time', 'T_eff', 'M_bol', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'Band_12', 'Band_13']
)

# 時間を秒から日（days）に変換
data['time_days'] = data['time'] / (60 * 60 * 24)

# グラフを作成 (例：ボロメトリック等級 vs 時間)
plt.figure(figsize=(10, 6))
plt.plot(data['time_days'], data['M_bol'])

# plt.xlim(0,40)
# plt.ylim(-21, -15)

plt.ylim(-23, -16)

# グラフの軸を反転させる (等級は数字が小さいほど明るいため)
plt.gca().invert_yaxis()

plt.xlabel('Time (days)')
plt.ylabel('Bolometric Magnitude (M_bol)')
plt.title('SNEC Supernova Light Curve')
plt.grid(True)
plt.show()