# SNEC-1.01 — SuperNova Explosion Code

**SNEC (SuperNova Explosion Code)** は、超新星爆発の流体力学・輻射輸送を一次元球対称ラグランジュ法で解く数値シミュレーションコードです。  
本リポジトリ (`SNEC-1.01`) は、特に **超新星噴出物(Ejecta)と星周物質(CSM: Circumstellar Medium)の相互作用** 、または **マグネター(Magnetar)駆動** をシミュレーションするために拡張・利用されたバージョンです。

---

## 📖 目次

- [概要](#概要)
- [物理モデル](#物理モデル)
- [CSM相互作用シミュレーション](#csm相互作用シミュレーション)
- [SNEC-1.01\_las（マグネター駆動版）について](#snec-101_lasマグネター駆動版について)
- [ディレクトリ構成](#ディレクトリ構成)
- [動作要件](#動作要件)
- [ビルド方法](#ビルド方法)
- [使い方](#使い方)
- [パラメータ設定](#パラメータ設定)
- [出力データ](#出力データ)
- [参考文献](#参考文献)

---

## 概要

SNEC は、Morozova et al. (2015) によって開発された、超新星爆発のライトカーブ（光度曲線）を計算するためのオープンソースコードです。Version 1.01 では、Version 1.00 に比べて以下の性能改善が施されています：

- **OPALオパシティテーブルの事前構築**: 各タイムステップでOPALルーチンを呼び出す代わりに、シミュレーション開始時に一度だけテーブルを構築
- **行列要素計算の最適化** (`arrays.F90`)
- **`nickel.F90` および `simple_saha.F90` ルーチンの最適化**
- **⁵⁶Co の平均寿命と収率の修正** (Nadyozhin, ApJS 92, 527, 1994 に基づく)

これらの改善により、SNEC-1.00 と比較して約 **2倍** の高速化を実現しています。

---

## 物理モデル

SNEC-1.01 が解く物理は以下の通りです：

### 流体力学（ラグランジュ法）
- **運動方程式**: 重力、圧力勾配、人工粘性を考慮した速度更新
- **質量保存**: ラグランジュ座標上での密度更新（球殻体積の変化から計算）
- **エネルギー方程式**: Newton-Raphson 反復法による陰解法で、温度・圧力・内部エネルギーを同時に求解

### 輻射輸送
- **フラックス制限拡散近似 (Flux-Limited Diffusion)**: 光子の輸送を拡散方程式で近似し、フラックスリミッター `λ` で光学的に薄い領域への遷移を適切に処理
- **陰的三重対角行列ソルバー**: LAPACK の `dgbsv` を使用して温度補正ベクトル δT を解く
- **光学的深さ (τ)**: 星表面から内部に向かって積分し、光球 (τ = 2/3) の位置を特定

### 状態方程式 (EOS)
- **理想気体 EOS** (`eoskey = 1`): Sedov テスト問題用
- **Paczynski EOS** (`eoskey = 2`): 輻射圧、電子の縮退圧、Saha 電離平衡を考慮した現実的な EOS

### エネルギー源
- **熱爆弾 (Thermal Bomb)**: 指定した質量範囲にエネルギーを注入して爆発を駆動
- **ピストン爆発 (Piston Explosion)**: 内側境界を指定速度で押し出す
- **⁵⁶Ni → ⁵⁶Co → ⁵⁶Fe 崩壊加熱**: Swartz et al. (1995) に基づくガンマ線輸送を含む放射性崩壊エネルギーの沈着計算

### オパシティ
- **OPALオパシティテーブル**: 高温領域（log T > 3.75）の詳細なオパシティを化学組成ごとにテーブル化
- **オパシティ下限値 (Floor)**: コアとエンベロープで異なる下限値を設定し、数値的安定性を確保

---

## CSM相互作用シミュレーション

本バージョンの主要な特徴は、**超新星噴出物とCSMの相互作用** のシミュレーションです。

### CSMとは
星周物質 (CSM) は、超新星爆発の前に親星（前駆体天体）が大質量放出や恒星風によって周囲に形成した物質のことです。超新星の噴出物がこのCSMに衝突すると、強い衝撃波が発生し、運動エネルギーが熱と光に変換されます。これにより、通常の超新星に比べて明るく長期間にわたる光度曲線（IIn型超新星など）が観測されます。

### 本コードでの実装
SNEC-1.01 では、CSM相互作用を以下の方法で記述します：

1. **CSM付き初期プロファイル**: 親星の外側にCSMの密度分布を付加した初期モデル（`simple1d_plus_csm.short` など）を使用。プロファイルには密度、温度、速度、電子分率(Ye) などが含まれます。

2. **流体力学的な相互作用**: 爆発によって噴出物がCSMに衝突すると、コード内の運動方程式と人工粘性によって衝撃波が自動的に捕捉・解像されます。衝撃波面での運動エネルギーから熱エネルギーへの変換が物理的に計算されます。

3. **衝撃波追跡 (`shock_capture.F90`)**: 衝撃波の現在位置を自動追跡し、ショック・ブレイクアウト（衝撃波が星表面に到達する瞬間）を検知します。

4. **CSM電離効果**: `csm_ionization` フラグにより、CSM中の電離状態に応じたオパシティを考慮することが可能です。これにより、光球の位置と光度曲線がCSMの電離状態の影響を正確に反映します。

5. **観測量の計算 (`analysis.F90`)**: 
   - 光球半径・光度・有効温度の時間発展
   - ボロメトリック補正による多バンド等級計算
   - 拡散時間・膨張時間の診断

### 典型的なCSMシミュレーション設定例

```
profile_name     = "profiles/simple1d_plus_csm.short"
comp_profile_name = "profiles/simple1d_plus_csm.iso.dat"
initial_data     = ""              # 追加のエネルギー注入なし（CSMの運動エネルギー変換で駆動）
imax             = 400
gridding         = "same_with_input"
mass_excision    = 0
radiation        = 1
eoskey           = 2
Ni_switch        = 1
Ni_mass          = 0.05
tend             = 8000000.0d0     # ~93日
```

---

## SNEC-1.01\_las（マグネター駆動版）について

同リポジトリ内の `SNEC-1.01_las` ディレクトリには、**マグネター（強磁場中性子星）からのスピンダウン放射** によるエネルギー注入を組み込んだ別バージョンが含まれています。

### 主な追加機能
- **`magnetar_heat.F90`**: マグネターからのスピンダウン光度を計算し、指定したゾーン範囲にエネルギーを注入するルーチン
- スピンダウン光度は以下のモデルに従います：  
  $L_{\mathrm{sd}}(t) = L_{\mathrm{sd,ini}} \left(1 + \frac{t}{t_{\mathrm{sd}}}\right)^{-2}$
- 周期的な変動を加えたモデル（正弦波の振幅・周期を設定可能）もコメントアウトとして含まれています

### マグネター関連パラメータ

| パラメータ | 説明 | 例 |
|-----------|------|----|
| `Lsd_ini` | 初期スピンダウン光度 [erg/s] | `1.0d44` |
| `tsd` | スピンダウンタイムスケール [日] | `1.3d1` |
| `kout` | エネルギー注入対象のセル数 | `400` |

---

## ディレクトリ構成

```
SNEC-1.01/
├── src/                    # Fortran ソースコード
│   ├── snec.F90            # メインプログラム（時間積分ループ）
│   ├── blmod.F90           # グローバル変数モジュール定義
│   ├── problem.F90         # 初期条件の設定
│   ├── hydro_rad.F90       # 輻射流体力学ソルバー（陰解法）
│   ├── hydro.F90           # 純流体力学ソルバー（陽解法）
│   ├── blstep.F90          # 1タイムステップの実行（リトライ機能付き）
│   ├── analysis.F90        # 光球・観測量の計算
│   ├── read_profile.F90    # 初期プロファイルの読み込み・グリッドへのマッピング
│   ├── nickel.F90          # ⁵⁶Ni 崩壊加熱計算
│   ├── opacity.F90         # オパシティ計算
│   ├── opal_opacity.F90    # OPAL オパシティテーブルの構築
│   ├── luminosity.F90      # 光度・光学的深さの計算
│   ├── input_parser.F90    # パラメータファイルの解析
│   ├── arrays.F90          # ヤコビアン行列要素の計算
│   ├── eos_content.F90     # 状態方程式（理想気体 / Paczynski）
│   ├── shock_capture.F90   # 衝撃波位置の追跡
│   ├── output.F90          # データ出力ルーチン
│   ├── Makefile            # ビルドルール
│   └── ...
├── profiles/               # 初期モデル（親星プロファイル）
│   ├── simple1d_plus_csm.short    # CSM付きモデル
│   ├── simple1d_plus_csm.iso.dat  # CSM付きモデル（化学組成）
│   ├── 15Msol_RSG.short           # 15太陽質量 赤色超巨星モデル
│   ├── stripped_star.short        # エンベロープ除去モデル
│   └── sedov.short                # Sedov テスト用モデル
├── tables/                 # 物理テーブル
│   ├── BolCorr.dat         # ボロメトリック補正テーブル
│   ├── GridPattern.dat     # グリッドパターン定義
│   └── codata[a-e]         # OPAL オパシティデータ
├── parameters              # シミュレーション設定ファイル
├── make.inc                # コンパイラ・ライブラリ設定
├── Makefile                # トップレベル Makefile
├── Data/                   # 出力データディレクトリ（要事前作成）
└── setup_*.f               # 初期モデル生成用 Fortran プログラム
```

---

## 動作要件

| 項目 | 要件 |
|------|------|
| **Fortranコンパイラ** | gfortran (GCC) 推奨 |
| **数値計算ライブラリ** | LAPACK, BLAS |
| **OS** | Linux, macOS, Windows (WSL or MinGW) |
| **Python** (オプション) | プロット・解析スクリプト用 |

---

## ビルド方法

### 1. `make.inc` の設定

環境に合わせてコンパイラとライブラリのパスを設定します。

**Linux / WSL の場合:**
```makefile
export F90=gfortran
export F90FLAGS=-g -O3 -Warray-bounds -fbounds-check
export LAPACKLIBS=-llapack -lblas
```

**macOS の場合:**
```makefile
export F90=gfortran
export F90FLAGS=-g -O3 -Warray-bounds -fbounds-check
export LAPACKLIBS=-framework Accelerate
```

### 2. コンパイル

```bash
cd src
make clean
make
```

ビルドが成功すると、親ディレクトリに `snec` (Linux/macOS) または `snec.exe` (Windows) が生成されます。

---

## 使い方

### 1. 出力ディレクトリの作成

```bash
mkdir Data
```

### 2. パラメータファイルの編集

`parameters` ファイルを編集して、使用するプロファイルや爆発条件を設定します。

### 3. シミュレーションの実行

```bash
./snec
```

### 4. 結果の確認

出力ディレクトリ（例: `Data/`）に時系列データが書き出されます。

---

## パラメータ設定

`parameters` ファイルの主要な設定項目は以下の通りです：

### 初期モデル
| パラメータ | 説明 |
|-----------|------|
| `profile_name` | 初期プロファイルファイル（密度・温度・速度など） |
| `comp_profile_name` | 化学組成プロファイルファイル |

### 爆発設定
| パラメータ | 説明 |
|-----------|------|
| `initial_data` | 爆発の種類（`"Thermal_Bomb"`, `"Piston_Explosion"`, `""` (なし)） |
| `final_energy` | 熱爆弾の注入エネルギー [erg] |
| `bomb_tstart` / `bomb_tend` | エネルギー注入の開始・終了時刻 [s] |

### グリッド
| パラメータ | 説明 |
|-----------|------|
| `imax` | グリッドゾーン数 |
| `gridding` | グリッド生成法（`"uniform_in_mass"`, `"from_file_by_mass"`, `"same_with_input"`） |
| `mass_excision` | 中心質量切り捨てフラグ (0 or 1) |
| `mass_excised` | 切り捨てる質量 [太陽質量] |

### 物理設定
| パラメータ | 説明 |
|-----------|------|
| `radiation` | 輻射輸送の有効化 (0 or 1) |
| `eoskey` | 状態方程式の選択（1: 理想気体, 2: Paczynski） |
| `Ni_switch` | ⁵⁶Ni 崩壊加熱の有効化 (0 or 1) |
| `Ni_mass` | ⁵⁶Ni の総質量 [太陽質量] |
| `opacity_floor_envelope` | エンベロープのオパシティ下限値 |
| `opacity_floor_core` | コアのオパシティ下限値 |

### 時間制御
| パラメータ | 説明 |
|-----------|------|
| `tend` | シミュレーション終了時刻 [s] |
| `dtout` | プロファイルデータの出力時間間隔 [s] |
| `dtout_scalar` | スカラ量（光度曲線など）の出力時間間隔 [s] |
| `dtmin` / `dtmax` | 時間刻み幅の許容範囲 [s] |

---

## 出力データ

シミュレーション結果は指定した出力ディレクトリに以下の形式で書き出されます：

### プロファイルデータ (`.xg` ファイル)
時系列の空間プロファイルを質量座標に対して出力します。

| ファイル名 | 内容 |
|-----------|------|
| `vel.xg` | 速度 [cm/s] |
| `rho.xg` | 密度 [g/cm³] |
| `temp.xg` | 温度 [K] |
| `press.xg` | 圧力 [dyn/cm²] |
| `lum.xg` | 光度 [erg/s] |
| `tau.xg` | 光学的深さ |
| `radius.xg` | 半径 [cm] |
| `kappa.xg` | オパシティ [cm²/g] |

### スカラ量データ (`.dat` ファイル)
時系列のスカラ量（光度曲線など）を出力します。

| ファイル名 | 内容 |
|-----------|------|
| `lum_observed.dat` | 観測ボロメトリック光度 [erg/s] |
| `lum_photo.dat` | 光球光度 [erg/s] |
| `T_eff.dat` | 有効温度 [K] |
| `rad_photo.dat` | 光球半径 [cm] |
| `vel_photo.dat` | 光球速度 [cm/s] |
| `magnitudes.dat` | 多バンド等級 (U, B, V, R, I 等) |
| `Ni_total_luminosity.dat` | ⁵⁶Ni 崩壊の総光度 [erg/s] |

---

## 参考文献

- **SNEC 原論文**: Morozova, V., Piro, A. L., Renzo, M., Ott, C. D., Clausen, D., Couch, S. M., Ellis, J., & Roberts, L. F. (2015). *Light Curves of Core-Collapse Supernovae with Substantial Mass Loss using the New Open-Source SuperNova Explosion Code (SNEC)*. ApJ, 814, 63.
- **Paczynski EOS**: Paczynski, B. (1983). ApJ, 267, 315.
- **Ni崩壊加熱**: Swartz, D. A., et al. (1995). ApJ, 446, 766.
- **⁵⁶Co 修正**: Nadyozhin, D. K. (1994). ApJS, 92, 527.
- **OPALオパシティ**: Iglesias, C. A. & Rogers, F. J. (1996). ApJ, 464, 943.
- **CSM相互作用**: Moriya, T. J. & Morozova, V. (2018) に関連する議論を参照。

---

## ライセンス

SNEC はオープンソースコードです。使用・改変・再配布の条件については、原著論文の著者グループのガイドラインに従ってください。
