import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import io
import os

# --- Helper function to find lines AFTER the last data block marker ---
def get_last_data_block_lines(filename):
    """
    Reads the file and returns only the lines after the last '"' marker.
    """
    lines_after_marker = []
    last_marker_line_index = -1
    all_lines = []
    try:
        with open(filename, 'r') as f:
            all_lines = f.readlines()
        # Find the index of the last line starting with '"'
        for i, line in reversed(list(enumerate(all_lines))):
            if line.strip().startswith('"'):
                last_marker_line_index = i
                break
        # Get lines after the marker
        if last_marker_line_index != -1:
            lines_after_marker = all_lines[last_marker_line_index + 1:]
        else: # If no marker found, assume all lines are data (might need adjustment)
             print(f"警告: '{filename}' 内に '\"' で始まる行が見つかりませんでした。ファイル全体をデータとして扱います。")
             lines_after_marker = all_lines
        # Filter out potentially empty lines at the very end
        lines_after_marker = [line for line in lines_after_marker if line.strip()]
        return "".join(lines_after_marker) # Return as a single string for pandas
    except FileNotFoundError:
        return None # Indicate file not found

# --- Define file paths ---
snec_output_dir = 'Data_sedov'
snec_rho_file = os.path.join(snec_output_dir, 'rho.xg')
snec_radius_file = os.path.join(snec_output_dir, 'radius.xg')
analytical_file = 'sedov_rho_rad.xg'

# --- 1. SNECシミュレーション結果の読み込み ---
try:
    snec_rho_lines = get_last_data_block_lines(snec_rho_file)
    snec_radius_lines = get_last_data_block_lines(snec_radius_file)

    if snec_rho_lines is None or snec_radius_lines is None:
        raise FileNotFoundError

    # Use io.StringIO to read the extracted lines string as if it were a file
    snec_density_data = pd.read_csv(
        io.StringIO(snec_rho_lines), sep='\s+', header=None,
        names=['_ignore', 'density'], usecols=['density']
    )
    snec_radius_data = pd.read_csv(
        io.StringIO(snec_radius_lines), sep='\s+', header=None,
        names=['_ignore', 'radius'], usecols=['radius']
    )

    snec_radius = snec_radius_data['radius']
    snec_density = snec_density_data['density']

    if len(snec_density) != len(snec_radius):
        print(f"警告: 読み込んだ '{snec_rho_file}' と '{snec_radius_file}' の行数が一致しません。")
        min_len = min(len(snec_density), len(snec_radius))
        snec_radius = snec_radius[:min_len]
        snec_density = snec_density[:min_len]

except FileNotFoundError:
    print(f"エラー: SNECの結果ファイル '{snec_rho_file}' または '{snec_radius_file}' が見つかりません。")
    exit()
except Exception as e:
    print(f"SNEC結果ファイルの読み込み中にエラーが発生しました: {e}")
    exit()

# --- 2. sedov.py が出力した解析解ファイルの読み込み ---
try:
    analytical_lines = get_last_data_block_lines(analytical_file)
    if analytical_lines is None:
         # Optional: try running sedov.py if file not found
        print(f"'{analytical_file}' が見つかりません。sedov.py を実行します...")
        try:
            subprocess.run(['python', 'sedov.py'], check=True)
            analytical_lines = get_last_data_block_lines(analytical_file)
            if analytical_lines is None: raise FileNotFoundError
        except Exception as run_e:
             print(f"sedov.py の実行またはファイル検索中にエラー: {run_e}")
             exit()

    # Read the extracted lines string
    sedov_data = pd.read_csv(
        io.StringIO(analytical_lines),
        sep='\s+',
        header=None,
        names=['radius_sedov', 'density_sedov'] # Expecting 2 columns
    )
    analytical_radius = sedov_data['radius_sedov']
    analytical_density = sedov_data['density_sedov']

    # --- DEBUG Print ---
    print("--- Analytical Solution Data Read by Pandas (first 5 rows): ---")
    print(sedov_data.head())
    print("---------------------------------------------------------------")
    # --- END DEBUG Print ---

except FileNotFoundError:
    print(f"エラー: 解析解ファイル '{analytical_file}' が見つかりません。先に python sedov.py を実行してください。")
    exit()
except Exception as e:
    print(f"解析解ファイル '{analytical_file}' の読み込み中にエラーが発生しました: {e}")
    exit()

# --- 3. グラフのプロット ---
plt.figure(figsize=(10, 7))

plt.plot(snec_radius, snec_density, label=f'SNEC Simulation (t={4.0} s)', marker='o', linestyle='-', markersize=4, color='blue')
# Make analytical solution visible
plt.plot(analytical_radius, analytical_density, label='Analytical Solution (Sedov)', linestyle='--', color='green', linewidth=4)

plt.xlabel('Radius (cm)')
plt.ylabel('Density (g/cm^3)')
plt.title('Sedov Blast Wave Test: SNEC vs Analytical Solution')
plt.legend()
plt.grid(True)

plt.yscale('log')
plt.xscale('log')

plt.show()