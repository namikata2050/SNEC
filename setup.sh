#!/bin/bash

echo "Installing required packages for SNEC..."

# Fortranコンパイラ、Make、LAPACKライブラリのインストール
sudo apt-get update
sudo apt-get install -y gfortran make liblapack-dev libblas-dev

# Pythonとそのパッケージマネージャのインストール
sudo apt-get install -y python3 python3-pip

# 解析スクリプト用のPythonライブラリをインストール
# (仮想環境を使用する場合は適宜変更してください)
pip3 install pandas matplotlib numpy

echo "Setup completed successfully."
echo "To compile SNEC, run: make"