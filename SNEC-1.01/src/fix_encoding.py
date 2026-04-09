import os
import glob

# 現在のディレクトリ内の .F90 ファイルをすべて取得
files = glob.glob("*.F90") + glob.glob("*.f90")

print(f"Found files: {files}")

for filename in files:
    try:
        # バイナリモードで読み込む
        with open(filename, 'rb') as f:
            content = f.read()
        
        # 0xC2 0xA0 (No-break space) を 0x20 (通常のスペース) に置換
        # また、UTF-8のBOMなどが混入している場合も考慮
        new_content = content.replace(b'\xc2\xa0', b'\x20')
        
        # 全角スペースなどが混ざっている場合も考慮して置換 (E3 80 80 -> 20 20)
        new_content = new_content.replace(b'\xe3\x80\x80', b'\x20\x20')

        if content != new_content:
            print(f"Fixing: {filename}")
            with open(filename, 'wb') as f:
                f.write(new_content)
        else:
            print(f"Clean: {filename}")

    except Exception as e:
        print(f"Error processing {filename}: {e}")

print("Done.")