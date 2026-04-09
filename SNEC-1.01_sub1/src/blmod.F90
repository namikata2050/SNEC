! blmod.F90: グローバル変数定義モジュール
!
! このファイルは、SNECシミュレーションの「状態」を保持する中心的な役割を果たします。
! プログラム全体から参照・更新される重要なデータ（物理量、パラメータ、定数）が
! 4つの異なるモジュールに分けて定義されています。

!############################# BLMOD MODULE ###################################
!
! module blmod: シミュレーションの「状態変数」モジュール
!
! シミュレーションの「今」の状態を保持する、最も重要なモジュールです。
! 物理量（密度、温度、速度など）の現在値と1ステップ前の値、
! シミュレーションの時間、グリッド情報などが含まれます。
!
module blmod

 implicit none

! ===== インフラ・制御フラグ =====
! some infrastructure stuff
 logical :: wipe_outdir = .true.    ! used in input_parser.F90 (出力ディレクトリを計算開始時に空にするか)
 logical :: switch1,switch2,switch3 ! for output control (汎用のデバッグ用スイッチ)
 
! physics controls
 logical :: scratch_step 		 ! タイムステップ計算が失敗した時に .true. になり、やり直しをトリガーするフラグ
 integer :: gravity_switch 		 ! 重力計算のON/OFF (0=OFF, 1=ON)。テスト問題(Sedovなど)で使用。

! ===== 初期プロファイル（星）のパラメータ =====
!some profile parameters
 integer :: ncomps  	 	 	 ! number of isotopes (計算に含む元素・同位体の数)
 real*8  :: Rstar  	 	 	 	 ! radius of the model (初期モデルの星の半径)
 real*8  :: total_initial_energy  	 ! initial energy of the model (初期モデルの総エネルギー)

! ===== 爆発のパラメータ =====
!explosion parameters
 logical :: do_piston 			 ! ピストンによる爆発を実行するか
 logical :: do_bomb 			 ! 熱爆弾による爆発を実行するか
 real*8  :: bomb_total_energy 	 ! 爆弾で注入する総エネルギー [erg]
 integer :: bomb_spread 		 ! 爆弾のエネルギーを何ゾーンに広げるか
 real*8, allocatable :: bomb_heating(:) ! 各ゾーンの爆弾加熱率 [erg/s/g] (bomb_pattern.F90で計算)
 
 ! Time centering term for radiative transfer (輻射輸送の陰解法における時間中心化パラメータ)
 ! theta = 0.5d0 は、時間的に2次の精度を持つ Crank-Nicolson 法を意味します。
 real*8, parameter :: theta = 0.5d0

! envelope metallicity is the nominal metallicity of the OPAL Type II tables
! changing the tables, don't forget to change this number
! OPALオパシティテーブルで使用する基準金属量
 real*8, parameter :: envelope_metallicity = 0.02d0
 
 ! ===== 時間・ステップ制御変数 =====
 ! evolution vars
 real*8 :: time 	 ! 現在のシミュレーション時刻 [s]
 real*8 :: dtime 	 ! 現在の時間刻み幅 (Δt) [s]
 real*8 :: dtime_p 	 ! 1ステップ前の時間刻み幅 (previous) [s]

 real*8 :: tstart 	 ! シミュレーション開始時刻 [s]
 real*8 :: tdump_check 	 ! 次にチェックポイントデータを出力する時刻
 real*8 :: tdump 	 ! 次に詳細データを出力する時刻
 real*8 :: tdump_scalar  ! 次にスカラ量(光度曲線など)を出力する時刻

 real*8 :: dtfac = 0.4d0         ! CFL safety factor (defaults to 0.4)
 real*8 :: cvisc = 10.0d0         ! Artificial viscosity coefficient (defaults to 2.0)

 integer :: nt 		 ! timestep (現在のタイムステップ数)
 integer :: ntstart 	 ! 開始時のタイムステップ数

! _p variables contain data from previous timestep
! ( _p が付く変数は、1ステップ前の(previous)値を保持します )

! ===== グリッド変数の宣言 (スタガード・メッシュ) =====
! allocatable :: ...(:) は、サイズが後で決まる「動的配列」を宣言します。
! 
! ------------------------------------------------------------------
! variables defined at cell inner boundaries (セルの「境界」で定義される量)
! ( i 番目の境界は、i-1 番目のゾーンと i 番目のゾーンの間にあります)
! ------------------------------------------------------------------
real*8, allocatable :: mass(:)  	 	 	 	 	 ! mass coordinate (ラグランジュ質量座標 m(i)) [g]
real*8, allocatable :: delta_cmass(:) 	 		 ! 境界間の質量 (m(i+1) - m(i)) ? (詳細要確認)
real*8, allocatable :: r(:), r_p(:)  	 	 		 ! radius (境界の半径) [cm]
real*8, allocatable :: lambda(:)  	 	 		 ! flux-limiter (輻射のフラックス・リミッター)
real*8, allocatable :: lum(:)  	 	 	 		 ! luminosity (境界を通過する光度) [erg/s]
real*8, allocatable :: vel(:), vel_p(:)  	 	 ! velocity (境界の速度) [cm/s]
real*8, allocatable :: inv_kappa(:)  	 		 ! inverse opacity (オパシティの逆数)

! ------------------------------------------------------------------
! variables defined at cell centers (セルの「中心」で定義される量)
! ( i 番目のゾーンは、i 番目の境界と i+1 番目の境界の間にあります)
! ------------------------------------------------------------------
real*8, allocatable :: cmass(:) 	 			 ! ゾーン中心の質量座標
real*8, allocatable :: delta_mass(:) 			 ! ゾーンの質量 (Δm_i) [g]
real*8, allocatable :: cr(:),  	cr_p(:) 		 ! ゾーン中心の半径 [cm]
real*8, allocatable :: p(:),  	p_p(:)  		 ! pressure (圧力) [dyn/cm^2]
real*8, allocatable :: eps(:),   eps_p(:)  		 ! specific internal energy (比内部エネルギー) [erg/g]
real*8, allocatable :: temp(:), temp_p(:)  		 ! temperature (温度) [K]
real*8, allocatable :: rho(:),   rho_p(:)  		 ! density (密度) [g/cm^3]
real*8, allocatable :: kappa_table(:)  	 		 ! tabular opacity (オパシティテーブルから引いた生のオパシティ)
real*8, allocatable :: kappa(:), kappa_p(:)  	 ! opacity (floor applied) (下限値(floor)を適用したオパシティ)
real*8, allocatable :: dkappadt(:) 			 	 ! オパシティの温度微分 (dκ/dT)
real*8, allocatable :: ye(:)  	 	 	 		 ! electron fraction (電子分率 Y_e)
real*8, allocatable :: abar(:) 				 	 ! 平均原子量
real*8, allocatable :: cs2(:)  	 	 	 		 ! speed of sound squared (音速の2乗) [cm^2/s^2]
real*8, allocatable :: comp(:,:)  	 	 		 ! mass fractions of elem-s (化学組成の質量分率 comp(zone, species))
real*8, allocatable :: comp_details(:,:)  		 ! mass and atomic numbers (化学組成の詳細: [species, (質量数, 原子番号)])
real*8, allocatable :: metallicity(:)  	 		 ! metallicity (金属量)

real*8, allocatable :: Q(:)  	 	 	 		 ! artificial viscosity (人工粘性 Q)
! term in the energy equation with artificial viscosity
real*8, allocatable :: Qterm(:) 				 ! 人工粘性による加熱項 [erg/s/g]

! derivatives of energy/pressure with respect to temperature
! (陰解法ソルバー(Newton-Raphson)のヤコビアン行列に必要な微分係数)
real*8, allocatable :: dedt(:) 				 	 ! d(eps)/dT (比内部エネルギーの温度微分)
real*8, allocatable :: dpdt(:) 				 	 ! d(p)/dT (圧力の温度微分)

real*8, allocatable :: entropy(:)  	 			 ! entropy (エントロピー)
real*8, allocatable :: tau(:)  	 	 			 ! optical depth (光学的深さ)
real*8, allocatable :: p_rad(:)  	 			 ! radiation pressure (輻射圧)

 real*8, allocatable :: delta_time(:) 			 ! (用途要確認) おそらく各ゾーンのタイムステップ制限

! ===== 電離状態 (Saha EOS) 関連 =====
 real*8, allocatable :: zav(:,:) 		! average charge of heavy particles (重粒子の平均電荷)
 real*8, allocatable :: ion_fractions(:,:,:) ! ionization fractions (電離度 (species:state:gridpoint))
 real*8, allocatable :: free_electron_frac(:) ! 自由電子の割合

! ===== [追加] CSMの電離スイッチ =====
  logical :: csm_ionization = .false.  ! ここで .true. か .false. を直接指定する
  ! ==================================

 ! ===== 光球 (Photosphere) 関連 (analysis.F90で計算) =====
 ! quantities at the photosphere
 integer :: index_photo 					! 光球が存在するグリッドのインデックス
 real*8 :: lum_photo 						! 光球での光度
 real*8 :: mass_photo 					! 光球での質量座標
 real*8 :: vel_photo 						! 光球での速度
 real*8 :: rad_photo 						! 光球での半径
 integer :: photosphere_fell_on_the_center ! 光球が中心(i=1)に達したかどうかのフラグ
 integer :: Ni_contributes_five_percents 	! Ni崩壊の寄与が5%を超えたかどうかのフラグ
 real*8, allocatable :: photosphere_tracer(:) ! 光球の位置を可視化するためのトレーサー (index_photo のみ 1)

 integer :: index_lumshell 		! index of the luminosity shell (輝度シェルのインデックス)
 real*8 :: mass_lumshell  	 	! mass coordinate of the luminosity shell (輝度シェルの質量座標)

 real*8 :: lum_observed  	 	! observed luminosity (観測される光度 = 光球からの光度 + 外部のNi崩壊)

! ===== 化学組成のインデックス =====
! numbers of some important elements in the composition profile
! (comp(:,:) 配列の何番目がどの元素に対応するかを保持する)
 integer :: H_number, He_number, C_number, O_number, Ni_number

! ===== 解析用 (analysis.F90) 変数 =====
! variables, used in analysis
 real*8, allocatable :: E_shell(:)  	! internal energy of the shells (i番目のシェルより外側の総内部エネルギー)
 real*8, allocatable :: time_diff(:)  	! characteristic diffusion time (輻射の拡散時間)
 real*8, allocatable :: time_exp(:)  	! characteristic expansion time (膨張時間)

! ===== 衝撃波 (Shock) 関連 (shock_capture.F90で計算) =====
! variables, needed for the shock capturing (shock_capture.F90)
 integer :: shockpos 			! 衝撃波の現在の位置（グリッドインデックス）
 integer :: shockpos_prev 		! 衝撃波の1ステップ前の位置
 integer :: shockpos_stop 		! 衝撃波の追跡を停止するかのフラグ
 integer :: breakoutflag 		! ショック・ブレイクアウトが発生したかのフラグ
 real*8 :: radshock 			! 衝撃波の現在の半径
 real*8 :: radshock_prev 		! 衝撃波の1ステップ前の半径
 real*8 :: time_prev 			! 衝撃波の位置を記録した1ステップ前の時刻

! ===== 等級計算 (analysis.F90) 関連 =====
! quantities needed for the color magnitude calculations (analysis.F90)
 real*8, allocatable :: bol_corr(:,:) 	! ボロメトリック補正 (BC) のテーブル
 real*8, allocatable :: temp_bol_corr(:)! BCテーブルの温度軸
 integer :: nlines_bol_corr 			! BCテーブルの行数
 real*8 :: magnitudes(11) 				! 計算された等級 (U, B, V, R, I など11バンド分)
 real*8 :: T_eff 						! 有効温度 (T_effective)

 real*8 :: T_col 						! カラー温度 (T_color)を追加

! ===== Ni-56 崩壊加熱 (nickel_heating.F90) 関連 =====
! quantities related to the heating by radioactive Ni
 real*8, allocatable :: Ni_deposit_function(:) ! ガンマ線エネルギー沈着関数 (実効的なNi質量分率)
 real*8, allocatable :: Ni_heating(:) 			! 各ゾーンのNiによる加熱率 [erg/s/g]
 real*8 :: Ni_total_luminosity 				! 沈着したNi崩壊の総光度 [erg/s]
 real*8 :: time_Ni 							! 次にNi加熱を再計算する時刻
 real*8 :: Ni_energy_rate 					! 放射性物質1gあたりのエネルギー放出率 [erg/s/g_Ni]

! ===== オパシティ (opacity.F90) 関連 =====
! quantities used for the calculations of the opacity (see opacity.F90)
 real*8, allocatable :: xxc(:), xxo(:) 		! OPALテーブル用の炭素(C)・酸素(O)のエンハンスメント量
 real*8, allocatable :: logR_op(:) 		! オパシティテーブルの座標 log(R_op) (R_op = rho/T_6^3)
 real*8, allocatable :: logT(:) 			! オパシティテーブルの座標 log(T)
 real*8, allocatable :: opacity_floor(:) 	! オパシティの下限値 (数値的安定性のため)
 integer :: opacity_corrupted 				! オパシティ計算がテーブル範囲外になったかのフラグ

 integer :: rpoints, tpoints 				! オパシティテーブルの解像度 (logR, logT 方向)
 ! オパシティのルックアップテーブル (3D配列: (化学組成k, 温度i, 密度j))
 real*8, allocatable :: opacity_tables(:,:,:) ! log kappa (gridpoint:logT:LogR)
 real*8, allocatable :: logT_array(:) 		! values of log T in the tables (テーブルの logT 軸)
 real*8, allocatable :: logR_array(:) 		! values of log R in the tables (テーブルの logR 軸)
 ! OPALテーブルの有効データ範囲（長方形ではないため）を定義する
 integer, allocatable :: op_rows(:,:) 
 integer, allocatable :: op_cols(:,:)

! other
 real*8 :: eos_gamma1 						! ideal_eos で使う断熱指数 γ

end module blmod

!############################# PARAMETERS MODULE ##############################
!
! module parameters: シミュレーションの「設定パラメータ」モジュール
!
! ユーザーが input_parser.F90 を介して設定する、シミュレーションの振る舞いを
! 決定するパラメータを格納します。
!
module parameters

   implicit none

!-------------------- launch (実行設定) --------------------------------
  character(len=256) :: outdir 				! 出力ディレクトリ名

!-------------------- profile (初期モデル設定) -------------------------------
  character(len=256) :: profile_name 			! 星の初期モデル（プロファイル）のファイル名
  character(len=256) :: composition_profile_name ! 化学組成プロファイルのファイル名

!------------------- explosion (爆発設定) ------------------------------
  character(len=256) :: initial_data 		! (用途要確認)

  !piston stuff (ピストン爆発のパラメータ)
  real*8  :: piston_vel 					! ピストンの速度
  real*8  :: piston_tstart 				! ピストン開始時刻
  real*8  :: piston_tend 					! ピストン終了時刻

  !thermal bomb stuff (熱爆弾のパラメータ)
  real*8  :: final_energy 					! 注入する総エネルギー
  real*8  :: bomb_tstart 					! エネルギー注入開始時刻
  real*8  :: bomb_tend 					! エネルギー注入終了時刻
  real*8  :: bomb_mass_spread 				! エネルギーを注入する質量範囲 [M_sun] (使われていない可能性あり)
  integer :: bomb_start_point 				! エネルギー注入を開始するグリッド番号
  integer :: bomb_mode 					! 爆弾エネルギーの定義モード (1 or 2)

!-------------------- grid (グリッド設定) ------------------------------------
  integer :: imax 							! number of zones (グリッドのゾーン(セル)の総数)
  character(len=256) :: gridding 			! グリッドの生成方法

  logical :: mass_excision 				! 中心の質量を切り捨てるか
  real*8  :: mass_excised 					! 切り捨てる質量 [M_sun]

!-------------------- evolution (物理計算の設定) -------------------------------
  logical :: radiation 					! 輻射輸送(radiation transport)を計算に含めるか
  integer :: eoskey 						! eos to use (使用する状態方程式(EOS)のキー: 1=ideal, 2=Paczynski)

  integer :: Ni_switch 					! Ni-56崩壊による加熱計算を行うか (0=OFF, 1=ON)
  integer :: Ni_by_hand 					! Ni-56の分布を手動で設定するか
  real*8 :: Ni_mass 						! 手動で設定するNi-56の総質量 [M_sun]
  real*8 :: Ni_boundary_mass 				! 手動で設定するNi-56の分布境界の質量座標 [M_sun]
  real*8 :: Ni_period 						! Ni加熱を再計算する時間間隔 [s]

  integer :: saha_ncomps 					! Sahaの電離公式で考慮する元素の数

  logical :: boxcar_smoothing 			! (用途要確認) スムージングのフラグ

  real*8 :: of_env 						! (用途要確認)
  real*8 :: of_core 					! (用途要確認)

!--------------------- timing (時間・出力設定) -----------------------------------
  integer :: ntmax 						! 最大計算ステップ数

  real*8 :: tend 							! シミュレーション終了時刻 [s]

  real*8 :: dtout, dtout_scalar, dtout_check ! 各データを出力する「時間」間隔 [s]
  integer :: ntout, ntout_scalar, ntout_check ! 各データを出力する「ステップ数」間隔

  integer :: ntinfo 						! 画面に進捗情報を出力するステップ数間隔

  real*8 :: dtmin 							! 許容される最小の時間刻み幅 [s]
  real*8 :: dtmax 							! 許容される最大の時間刻み幅 [s]

!---------------------- test (テスト問題設定) -------------------------------------
  logical :: sedov = .false. 				! Sedov-Taylor 爆風問題のテストモードで実行するか

end module parameters

!############################# EOS MODULE #####################################
!
! module eosmodule: 状態方程式(EOS)関連データ モジュール
!
! paczynski_eos などで使用する、サハの電離公式に必要な物理データ（電離ポテンシャル）
! を格納し、初期化します。
!
module eosmodule

 implicit none

 ! array for ionization potential
 ! xxip(j, k): 原子番号 j (行) の元素の k 番目 (列) の電離ポテンシャル
 real*8 :: xxip(30,30)
 ! サハの式で使う統計的重みの比
 real*8,allocatable :: stat_weight_p1_ratio(:,:)

 contains
  	! init_ionpot: 電離ポテンシャル配列(xxip)を初期化するサブルーチン
  	subroutine init_ionpot
  !#############################################################################
  !
  ! These data are taken from the Timmes EOS with Saha ionization:
  ! (これらのデータは Timmes EOS (論文参照) から引用されたものである)
  !
  ! http://cococubed.asu.edu/code_pages/eos_ionize.shtml
  !
  ! For more information on the Timmes EOS see:
  !  	 timmes & arnett, apj supp. 125, 277, 1999
  !  	 timmes & swesty, apj supp. 126, 501, 2000
  !
  !#############################################################################

  	 	use blmod, only: comp_details,ncomps
  	 	use parameters, only: saha_ncomps
  	 	implicit none
  	 	integer :: k,j,i
  	 	real*8, parameter :: ev2erg  = 1.602d-12 	! eV を erg に変換する係数
  	 	real*8, parameter :: off_table = 1.0d30 	! テーブルの無効な値を示すマーカー
  	 	real*8  :: zion
  	 	real*8  :: stat_weight ! function!
  	 	integer :: izion

  	 	xxip(:,:) = 0.0d0 ! 配列を 0 で初期化

! プリプロセッサディレクティブ (通常 #if 1 で有効)
#if 1
  	 	! set up ratio of statistical weights needed
  	 	! in simple_saha.F90
  	 	! (simple_saha.F90 で必要となる統計的重みの比を計算・確保する)
  	 	allocate(stat_weight_p1_ratio(ncomps,30))
  	 	stat_weight_p1_ratio(:,:) = 0.0d0
  	 	do j=1,saha_ncomps
  	 	 	zion = comp_details(j,2)
  	 	 	izion = int(zion)
  	 	 	do i=1,izion
  	 	 	 	stat_weight_p1_ratio(j,i) = stat_weight(zion,i+1)/stat_weight(zion,i)
  	 	 	enddo
  	 	enddo
#endif

		! ===== 電離ポテンシャル [eV] のハードコード =====
		! (Z=1 の水素から Z=30 の亜鉛まで)
  	 	!hydrogen
  	 	xxip(1,1)  	 = 13.59844d0
  	 	xxip(1,2:30) = off_table

  	 	!helium
  	 	xxip(2,1)  	 = 24.58741d0
  	 	xxip(2,2)  	 = 54.41778d0
  	 	xxip(2,3:30) = off_table

  	 	!lithium
  	 	xxip(3,1)  	 = 5.39172d0
  	 	xxip(3,2)  	 = 75.64018d0
  	 	xxip(3,3)  	 = 122.45429d0
  	 	xxip(3,4:30) = off_table


      !berrylium
      xxip(4,1)    = 9.3227d0
      xxip(4,2)    = 18.21116d0
      xxip(4,3)    = 153.89661d0
      xxip(4,4)    = 217.71865d0
      xxip(4,5:30) = off_table

      !boron
      xxip(5,1)    = 8.29803d0
      xxip(5,2)    = 25.15484d0
      xxip(5,3)    = 37.93064d0
      xxip(5,4)    = 259.37521d0
      xxip(5,5)    = 340.22580d0
      xxip(5,6:30) = off_table

      !carbon
      xxip(6,1)    = 11.26030d0
      xxip(6,2)    = 24.38332d0
      xxip(6,3)    = 47.8878d0
      xxip(6,4)    = 64.4939d0
      xxip(6,5)    = 392.087d0
      xxip(6,6)    = 489.99334d0
      xxip(6,7:30) = off_table

      !nitrogen
      xxip(7,1)    = 14.53414d0
      xxip(7,2)    = 29.6013d0
      xxip(7,3)    = 47.44924d0
      xxip(7,4)    = 77.4735d0
      xxip(7,5)    = 97.8902d0
      xxip(7,6)    = 552.0718d0
      xxip(7,7)    = 667.046d0
      xxip(7,8:30) = off_table

      !oxygen
      xxip(8,1)    = 13.61806d0
      xxip(8,2)    = 35.11730d0
      xxip(8,3)    = 54.9355d0
      xxip(8,4)    = 77.41353d0
      xxip(8,5)    = 113.8990d0
      xxip(8,6)    = 138.1197d0
      xxip(8,7)    = 739.29d0
      xxip(8,8)    = 871.4101d0
      xxip(8,9:30) = off_table

      !fluorine
      xxip(9,1)      = 17.42282d0
      xxip(9,2)      = 34.97082d0
      xxip(9,3)      = 62.7084d0
      xxip(9,4)      = 87.1398d0
      xxip(9,5)      = 114.2428d0
      xxip(9,6)      = 157.1651d0
      xxip(9,7)      = 185.18d0
      xxip(9,8)      = 953.9112d0
      xxip(9,9)      = 1103.1176d0
      xxip(9,10:30)  = off_table

      !neon
      xxip(10,1)      = 21.5646d0
      xxip(10,2)      = 40.96328d0
      xxip(10,3)      = 63.45d0
      xxip(10,4)      = 97.12d0
      xxip(10,5)      = 126.21d0
      xxip(10,6)      = 157.93d0
      xxip(10,7)      = 207.2759d0
      xxip(10,8)      = 239.0989d0
      xxip(10,9)      = 1195.8286d0
      xxip(10,10)     = 1362.1995d0
      xxip(10,11:30)  = off_table

      !sodium
      xxip(11,1)      = 5.13908d0
      xxip(11,2)      = 47.2864d0
      xxip(11,3)      = 71.6200d0
      xxip(11,4)      = 98.91d0
      xxip(11,5)      = 138.40d0
      xxip(11,6)      = 172.18d0
      xxip(11,7)      = 208.50d0
      xxip(11,8)      = 264.25d0
      xxip(11,9)      = 299.864d0
      xxip(11,10)     = 1465.121d0
      xxip(11,11)     = 1648.702d0
      xxip(11,12:30)  = off_table

      !magnesium
      xxip(12,1)      = 7.64624d0
      xxip(12,2)      = 15.03528d0
      xxip(12,3)      = 80.1437d0
      xxip(12,4)      = 109.2655d0
      xxip(12,5)      = 141.27d0
      xxip(12,6)      = 186.76d0
      xxip(12,7)      = 225.02d0
      xxip(12,8)      = 265.96d0
      xxip(12,9)      = 328.06d0
      xxip(12,10)     = 367.50d0
      xxip(12,11)     = 1761.805d0
      xxip(12,12)     = 1962.6650d0
      xxip(12,13:30)  = off_table

      !aluminum
      xxip(13,1)      = 5.98577d0
      xxip(13,2)      = 18.82856d0
      xxip(13,3)      = 28.44765d0
      xxip(13,4)      = 119.992d0
      xxip(13,5)      = 153.825d0
      xxip(13,6)      = 190.49d0
      xxip(13,7)      = 241.76d0
      xxip(13,8)      = 284.66d0
      xxip(13,9)      = 330.13d0
      xxip(13,10)     = 398.75d0
      xxip(13,11)     = 442.00d0
      xxip(13,12)     = 2085.98d0
      xxip(13,13)     = 2304.1410d0
      xxip(13,14:30)  = off_table

      !silicon
      xxip(14,1)      = 8.15169d0
      xxip(14,2)      = 16.34585d0
      xxip(14,3)      = 33.49302d0
      xxip(14,4)      = 45.14181d0
      xxip(14,5)      = 166.767d0
      xxip(14,6)      = 205.27d0
      xxip(14,7)      = 246.5d0
      xxip(14,8)      = 303.54d0
      xxip(14,9)      = 351.12d0
      xxip(14,10)     = 401.37d0
      xxip(14,11)     = 476.36d0
      xxip(14,12)     = 523.42d0
      xxip(14,13)     = 2437.63d0
      xxip(14,14)     = 2673.182d0
      xxip(14,15:30)  = off_table

      !phosphorous
      xxip(15,1)      = 10.48669d0
      xxip(15,2)      = 19.7694d0
      xxip(15,3)      = 30.2027d0
      xxip(15,4)      = 51.4439d0
      xxip(15,5)      = 65.0251d0
      xxip(15,6)      = 220.421d0
      xxip(15,7)      = 263.57d0
      xxip(15,8)      = 309.60d0
      xxip(15,9)      = 372.13d0
      xxip(15,10)     = 424.4d0
      xxip(15,11)     = 479.46d0
      xxip(15,12)     = 560.8d0
      xxip(15,13)     = 611.74d0
      xxip(15,14)     = 2816.91d0
      xxip(15,15)     = 3069.842d0
      xxip(15,16:30)  = off_table

      !sulfur
      xxip(16,1)      = 10.36001d0
      xxip(16,2)      = 23.3379d0
      xxip(16,3)      = 34.79d0
      xxip(16,4)      = 47.222d0
      xxip(16,5)      = 72.5945d0
      xxip(16,6)      = 88.0530d0
      xxip(16,7)      = 280.948d0
      xxip(16,8)      = 328.7d0
      xxip(16,9)      = 379.55d0
      xxip(16,10)     = 447.5d0
      xxip(16,11)     = 504.8d0
      xxip(16,12)     = 564.44d0
      xxip(16,13)     = 652.2d0
      xxip(16,14)     = 707.01d0
      xxip(16,15)     = 3223.78d0
      xxip(16,16)     = 3494.1892d0
      xxip(16,17:30)  = off_table

      !chlorine
      xxip(17,1)      = 12.96764d0
      xxip(17,2)      = 23.814d0
      xxip(17,3)      = 39.61d0
      xxip(17,4)      = 53.4652d0
      xxip(17,5)      = 67.8d0
      xxip(17,6)      = 97.03d0
      xxip(17,7)      = 114.1958d0
      xxip(17,8)      = 348.28d0
      xxip(17,9)      = 400.06d0
      xxip(17,10)     = 455.63d0
      xxip(17,11)     = 529.28d0
      xxip(17,12)     = 591.99d0
      xxip(17,13)     = 656.71d0
      xxip(17,14)     = 749.76d0
      xxip(17,15)     = 809.40d0
      xxip(17,16)     = 3658.521d0
      xxip(17,17)     = 3946.2960d0
      xxip(17,18:30)  = off_table

      !argon
      xxip(18,1)      = 15.75962d0
      xxip(18,2)      = 27.62967d0
      xxip(18,3)      = 40.74d0
      xxip(18,4)      = 59.81d0
      xxip(18,5)      = 75.02d0
      xxip(18,6)      = 91.009d0
      xxip(18,7)      = 124.323d0
      xxip(18,8)      = 143.460d0
      xxip(18,9)      = 422.45d0
      xxip(18,10)     = 478.69d0
      xxip(18,11)     = 538.96d0
      xxip(18,12)     = 618.26d0
      xxip(18,13)     = 686.10d0
      xxip(18,14)     = 755.74d0
      xxip(18,15)     = 854.77d0
      xxip(18,16)     = 918.03d0
      xxip(18,17)     = 4120.8857d0
      xxip(18,18)     = 4426.2296d0
      xxip(18,19:30)  = off_table

      !pottasium
      xxip(19,1)      = 4.34066d0
      xxip(19,2)      = 31.63d0
      xxip(19,3)      = 45.806d0
      xxip(19,4)      = 60.91d0
      xxip(19,5)      = 82.66d0
      xxip(19,6)      = 99.4d0
      xxip(19,7)      = 117.56d0
      xxip(19,8)      = 154.88d0
      xxip(19,9)      = 175.8174d0
      xxip(19,10)     = 503.8d0
      xxip(19,11)     = 564.7d0
      xxip(19,12)     = 629.4d0
      xxip(19,13)     = 714.6d0
      xxip(19,14)     = 786.6d0
      xxip(19,15)     = 861.1d0
      xxip(19,16)     = 968.0d0
      xxip(19,17)     = 1033.4d0
      xxip(19,18)     = 4610.8d0
      xxip(19,19)     = 4934.046d0
      xxip(19,20:30)  = off_table

      !calcium
      xxip(20,1)      = 6.11316d0
      xxip(20,2)      = 11.87172d0
      xxip(20,3)      = 50.9131d0
      xxip(20,4)      = 67.27d0
      xxip(20,5)      = 84.50d0
      xxip(20,6)      = 108.78d0
      xxip(20,7)      = 127.2d0
      xxip(20,8)      = 147.24d0
      xxip(20,9)      = 188.54d0
      xxip(20,10)     = 211.275d0
      xxip(20,11)     = 591.9d0
      xxip(20,12)     = 657.2d0
      xxip(20,13)     = 726.6d0
      xxip(20,14)     = 817.6d0
      xxip(20,15)     = 894.5d0
      xxip(20,16)     = 974.0d0
      xxip(20,17)     = 1087.0d0
      xxip(20,18)     = 1157.8d0
      xxip(20,19)     = 5128.8d0
      xxip(20,20)     = 5469.864d0
      xxip(20,21:30)  = off_table

      !scandium
      xxip(21,1)      = 6.5615d0
      xxip(21,2)      = 12.79967d0
      xxip(21,3)      = 24.75666d0
      xxip(21,4)      = 73.4894d0
      xxip(21,5)      = 91.65d0
      xxip(21,6)      = 110.68d0
      xxip(21,7)      = 138.0d0
      xxip(21,8)      = 158.1d0
      xxip(21,9)      = 180.03d0
      xxip(21,10)     = 225.18d0
      xxip(21,11)     = 249.798d0
      xxip(21,12)     = 687.36d0
      xxip(21,13)     = 756.7d0
      xxip(21,14)     = 830.8d0
      xxip(21,15)     = 927.5d0
      xxip(21,16)     = 1009.0d0
      xxip(21,17)     = 1094.0d0
      xxip(21,18)     = 1213.0d0
      xxip(21,19)     = 1287.97d0
      xxip(21,20)     = 5674.8d0
      xxip(21,21)     = 6033.712d0
      xxip(21,22:30)  = off_table

      !titanium
      xxip(22,1)      = 6.8281d0
      xxip(22,2)      = 13.5755d0
      xxip(22,3)      = 27.4917d0
      xxip(22,4)      = 43.2672d0
      xxip(22,5)      = 99.30d0
      xxip(22,6)      = 119.53d0
      xxip(22,7)      = 140.8d0
      xxip(22,8)      = 170.4d0
      xxip(22,9)      = 192.1d0
      xxip(22,10)     = 215.92d0
      xxip(22,11)     = 265.07d0
      xxip(22,12)     = 291.500d0
      xxip(22,13)     = 787.84d0
      xxip(22,14)     = 863.1d0
      xxip(22,15)     = 941.9d0
      xxip(22,16)     = 1044.0d0
      xxip(22,17)     = 1131.0d0
      xxip(22,18)     = 1221.0d0
      xxip(22,19)     = 1346.0d0
      xxip(22,20)     = 1425.4d0
      xxip(22,21)     = 6249.0d0
      xxip(22,22)     = 6625.82d0
      xxip(22,24:30)  = off_table

      !vanadium
      xxip(23,1)      = 6.7463d0
      xxip(23,2)      = 14.66d0
      xxip(23,3)      = 29.311d0
      xxip(23,4)      = 46.709d0
      xxip(23,5)      = 65.2817d0
      xxip(23,6)      = 128.13d0
      xxip(23,7)      = 150.6d0
      xxip(23,8)      = 173.4d0
      xxip(23,9)      = 205.8d0
      xxip(23,10)     = 230.5d0
      xxip(23,11)     = 255.7d0
      xxip(23,12)     = 308.1d0
      xxip(23,13)     = 336.277d0
      xxip(23,14)     = 896.0d0
      xxip(23,15)     = 976.0d0
      xxip(23,16)     = 1060.0d0
      xxip(23,17)     = 1168.0d0
      xxip(23,18)     = 1260.0d0
      xxip(23,19)     = 1355.0d0
      xxip(23,20)     = 1486.0d0
      xxip(23,21)     = 1569.6d0
      xxip(23,22)     = 6851.3d0
      xxip(23,23)     = 7246.12d0
      xxip(23,24:30)  = off_table

      !chromium
      xxip(24,1)      = 6.7665d0
      xxip(24,2)      = 16.4857d0
      xxip(24,3)      = 30.96d0
      xxip(24,4)      = 49.16d0
      xxip(24,5)      = 69.46d0
      xxip(24,6)      = 90.6349d0
      xxip(24,7)      = 160.18d0
      xxip(24,8)      = 184.7d0
      xxip(24,9)      = 209.3d0
      xxip(24,10)     = 244.4d0
      xxip(24,11)     = 270.8d0
      xxip(24,12)     = 298.0d0
      xxip(24,13)     = 354.8d0
      xxip(24,14)     = 384.168d0
      xxip(24,15)     = 1010.6d0
      xxip(24,16)     = 1097.0d0
      xxip(24,17)     = 1185.0d0
      xxip(24,18)     = 1299.0d0
      xxip(24,19)     = 1396.0d0
      xxip(24,20)     = 1496.0d0
      xxip(24,21)     = 1634.0d0
      xxip(24,22)     = 1721.4d0
      xxip(24,23)     = 7481.7d0
      xxip(24,24)     = 7894.81d0
      xxip(24,25:30)  = off_table

      !manganese
      xxip(25,1)      = 7.43402d0
      xxip(25,2)      = 15.63999d0
      xxip(25,3)      = 33.668d0
      xxip(25,4)      = 51.2d0
      xxip(25,5)      = 72.4d0
      xxip(25,6)      = 95.6d0
      xxip(25,7)      = 119.203d0
      xxip(25,8)      = 194.5d0
      xxip(25,9)      = 221.8d0
      xxip(25,10)     = 248.3d0
      xxip(25,11)     = 286.0d0
      xxip(25,12)     = 314.4d0
      xxip(25,13)     = 343.6d0
      xxip(25,14)     = 403.0d0
      xxip(25,15)     = 435.163d0
      xxip(25,16)     = 1134.7d0
      xxip(25,17)     = 1224.0d0
      xxip(25,18)     = 1317.0d0
      xxip(25,19)     = 1437.0d0
      xxip(25,20)     = 1539.0d0
      xxip(25,21)     = 1644.0d0
      xxip(25,22)     = 1788.0d0
      xxip(25,23)     = 1879.9d0
      xxip(25,24)     = 8140.6d0
      xxip(25,25)     = 8571.94d0
      xxip(25,26:30)  = off_table

      !iron
      xxip(26,1)      = 7.9024d0
      xxip(26,2)      = 16.1878d0
      xxip(26,3)      = 30.652d0
      xxip(26,4)      = 54.8d0
      xxip(26,5)      = 75.0d0
      xxip(26,6)      = 99.1d0
      xxip(26,7)      = 124.98d0
      xxip(26,8)      = 151.06d0
      xxip(26,9)      = 233.6d0
      xxip(26,10)     = 262.1d0
      xxip(26,11)     = 290.2d0
      xxip(26,12)     = 330.8d0
      xxip(26,13)     = 361.0d0
      xxip(26,14)     = 392.2d0
      xxip(26,15)     = 457.0d0
      xxip(26,16)     = 489.256d0
      xxip(26,17)     = 1266.0d0
      xxip(26,18)     = 1358.0d0
      xxip(26,19)     = 1456.0d0
      xxip(26,20)     = 1582.0d0
      xxip(26,21)     = 1689.0d0
      xxip(26,22)     = 1799.0d0
      xxip(26,23)     = 1950.0d0
      xxip(26,24)     = 2023.0d0
      xxip(26,25)     = 8828.0d0
      xxip(26,26)     = 9277.69d0
      xxip(26,27:30)  = off_table

      !cobalt
      xxip(27,1)      = 7.8810d0
      xxip(27,2)      = 17.083d0
      xxip(27,3)      = 33.50d0
      xxip(27,4)      = 51.3d0
      xxip(27,5)      = 79.5d0
      xxip(27,6)      = 102.0d0
      xxip(27,7)      = 128.9d0
      xxip(27,8)      = 157.8d0
      xxip(27,9)      = 186.13d0
      xxip(27,10)     = 275.4d0
      xxip(27,11)     = 305.0d0
      xxip(27,12)     = 336.0d0
      xxip(27,13)     = 379.0d0
      xxip(27,14)     = 411.0d0
      xxip(27,15)     = 444.0d0
      xxip(27,16)     = 511.96d0
      xxip(27,17)     = 546.58d0
      xxip(27,18)     = 1397.2d0
      xxip(27,19)     = 1504.6d0
      xxip(27,20)     = 1603.0d0
      xxip(27,21)     = 1735.0d0
      xxip(27,22)     = 1846.0d0
      xxip(27,23)     = 1962.0d0
      xxip(27,24)     = 2119.0d0
      xxip(27,25)     = 2219.0d0
      xxip(27,26)     = 9544.1d0
      xxip(27,27)     = 10012.12d0
      xxip(27,28:30)  = off_table

      !nickel
      xxip(28,1)      = 7.6398d0
      xxip(28,2)      = 18.16884d0
      xxip(28,3)      = 35.19d0
      xxip(28,4)      = 54.9d0
      xxip(28,5)      = 76.06d0
      xxip(28,6)      = 108.0d0
      xxip(28,7)      = 133.0d0
      xxip(28,8)      = 162.0d0
      xxip(28,9)      = 193.0d0
      xxip(28,10)     = 224.6d0
      xxip(28,11)     = 321.0d0
      xxip(28,12)     = 352.0d0
      xxip(28,13)     = 384.0d0
      xxip(28,14)     = 430.0d0
      xxip(28,15)     = 464.0d0
      xxip(28,16)     = 499.0d0
      xxip(28,17)     = 571.08d0
      xxip(28,18)     = 607.06d0
      xxip(28,19)     = 1541.0d0
      xxip(28,20)     = 1648.0d0
      xxip(28,21)     = 1756.0d0
      xxip(28,22)     = 1894.0d0
      xxip(28,23)     = 2011.0d0
      xxip(28,24)     = 2131.0d0
      xxip(28,25)     = 2295.0d0
      xxip(28,26)     = 2399.2d0
      xxip(28,27)     = 10288.8d0
      xxip(28,28)     = 10775.40d0
      xxip(28,29:30)  = off_table

      !copper
      xxip(29,1)     = 7.72638d0
      xxip(29,2)     = 20.29240d0
      xxip(29,3)     = 36.841d0
      xxip(29,4)     = 57.38d0
      xxip(29,5)     = 79.8d0
      xxip(29,6)     = 103.0d0
      xxip(29,7)     = 139.0d0
      xxip(29,8)     = 166.0d0
      xxip(29,9)     = 199.0d0
      xxip(29,10)    = 232.0d0
      xxip(29,11)    = 265.3d0
      xxip(29,12)    = 369.0d0
      xxip(29,13)    = 401.0d0
      xxip(29,14)    = 435.0d0
      xxip(29,15)    = 484.0d0
      xxip(29,16)    = 520.0d0
      xxip(29,17)    = 557.0d0
      xxip(29,18)    = 633.0d0
      xxip(29,19)    = 670.588d0
      xxip(29,20)    = 1697.0d0
      xxip(29,21)    = 1804.0d0
      xxip(29,22)    = 1916.0d0
      xxip(29,23)    = 2060.0d0
      xxip(29,24)    = 2182.0d0
      xxip(29,25)    = 2308.0d0
      xxip(29,26)    = 2478.0d0
      xxip(29,27)    = 2587.5d0
      xxip(29,28)    = 11062.38d0
      xxip(29,29)    = 11567.617d0
      xxip(29,30)    = off_table

  	 	!zinc
  	 	!the last nine are not in the cited reference
  	 	xxip(30,1)  	= 9.3942d0
  	 	xxip(30,2)  	= 17.96440d0
  	 	xxip(30,3)  	= 39.723d0
  	 	xxip(30,4)  	= 59.4d0
  	 	xxip(30,5)  	= 82.6d0
  	 	xxip(30,6)  	= 108.0d0
  	 	xxip(30,7)  	= 134.0d0
  	 	xxip(30,8)  	= 174.0d0
  	 	xxip(30,9)  	= 203.0d0
  	 	xxip(30,10) 	= 238.0d0
  	 	xxip(30,11) 	= 274.0d0
  	 	xxip(30,12) 	= 310.8d0
  	 	xxip(30,13) 	= 419.7d0
  	 	xxip(30,14) 	= 454.0d0
  	 	xxip(30,15) 	= 490.0d0
  	 	xxip(30,16) 	= 542.0d0
  	 	xxip(30,17) 	= 579.0d0
  	 	xxip(30,18) 	= 619.0d0
  	 	xxip(30,19) 	= 698.0d0
  	 	xxip(30,20) 	= 738.0d0
  	 	xxip(30,21) 	= 1856.0d0
  	 	xxip(30,22) 	= 1920.0d0
  	 	xxip(30,23) 	= 2075.0d0
  	 	xxip(30,24) 	= 2190.0d0
  	 	xxip(30,25) 	= 2350.0d0
  	 	xxip(30,26) 	= 2500.0d0
  	 	xxip(30,27) 	= 2600.0d0
  	 	xxip(30,28) 	= 11080.0d0
  	 	xxip(30,29) 	= 11600.0d0
  	 	xxip(30,30) 	= 12000.0d0

  	 	! ===== 単位の変換 [eV] -> [erg] =====
  	 	! 読み込んだすべての電離ポテンシャルを [eV] から CGS 単位系の [erg] に変換する
  	 	do k=1,30
  	 	 	 do j=1,30
  	 	 	 	if (xxip(j,k) .ne. off_table) xxip(j,k) = xxip(j,k) * ev2erg
  	 	 	 enddo
  	 	enddo


  	end subroutine init_ionpot

end module eosmodule


!######################## PHYSICAL CONSTANTS MODULE ###########################
!
! module physical_constants: 物理定数モジュール
!
! シミュレーション全体で使用する、変更不可能な物理定数を CGS 単位系で定義します。
! `parameter` 属性により、これらは定数として扱われます。
!
module physical_constants

  	implicit none
  	
  	! ===== 基本的な物理定数 (CGS) =====
  	real*8, parameter :: msun = 1.98892d33  	 	 ! solar mass (太陽質量) [g]
  	real*8, parameter :: rsun = 6.96d10  	 	 	 ! solar radius (太陽半径) [cm]
  	real*8, parameter :: clite = 2.99792458d10  	 ! speed of light (光速) [cm/s]
  	real*8, parameter :: ggrav = 6.6742d-8  	 	 ! gravitational constant (万有引力定数) [cgs]
  	real*8, parameter :: kboltz = 1.380662d-16 	 ! ボルツマン定数 [erg/K]
  	real*8, parameter :: mev_to_erg = 1.6022d-6 	 ! MeV から erg への変換
  	real*8, parameter :: pi = 3.14159265358979d0 	 ! 円周率
  	real*8, parameter :: emev = 1.60219d-6 		 ! (mev_to_erg とほぼ同じ)
  	real*8, parameter :: avo_real = 6.0221415d23 	 ! アボガドロ数
  	real*8, parameter :: h_cgs = 6.626058d-27 	 ! プランク定数 [erg s]
  	real*8, parameter :: a_rad = 7.5657d-15 		 ! 輻射定数 (a = 4σ/c) [erg cm^-3 K^-4]
  	real*8, parameter :: mproton = 1.67262178d-24 	 ! 陽子の質量 [g]
  	real*8, parameter :: melectron = 9.1093897d-28 ! 電子の質量 [g]
  	real*8, parameter :: kdnr = 9.91d12 		 	 ! 縮退電子圧(非相対論的)の係数
  	real*8, parameter :: kdr = 1.231d15 		 	 ! 縮退電子圧(相対論的)の係数
  	real*8, parameter :: sigma_SB = 5.6704d-5 	 ! シュテファン・ボルツマン定数 [erg cm^-2 s^-1 K^-4]
  	
  	! ===== 放射性崩壊 (Ni, Co) の時定数 =====
  	real*8, parameter :: tau_Ni = 760320.0d0 		 ! 8.8 days in seconds (Ni-56 の時定数 τ = T_1/2 / ln(2))
  	real*8, parameter :: tau_Co = 9616320.0d0  	 ! 111.3 days in seconds (Co-56 の時定数)
  	real*8, parameter :: overtau_Ni = 1.0d0/760320.0d0 ! 1 / tau_Ni
  	real*8, parameter :: overtau_Co = 1.0d0/9616320.0d0 ! 1 / tau_Co
  	
  	! ===== 太陽・組成の基準値 =====
  	! relative solar mass fraction of C (G&N'93) (太陽のCの質量分率)
  	real*8, parameter :: C_frac_sol = 0.173285d0  
  	! relative solar mass fraction of O (G&N'93) (太陽のOの質量分率)
  	real*8, parameter :: O_frac_sol = 0.482273d0  
  	! absolute bolometric magnitude of sun (太陽の絶対等級)
  	real*8, parameter :: sun_mag = 4.75d0  	 	  
  	! bolometric luminosity of sun in erg/s (太陽の光度) [erg/s]
  	real*8, parameter :: sun_lum = 3.846d33

  	! ===== EOS・ヘルパー定数 =====
  	! constants used in Saha solver
  	! サハの電離公式の係数 2 * (2π m_e k_B / h^2)^(3/2)
  	real*8, parameter :: saha_coeff = 2.0d0*(2.0d0*pi*melectron*kboltz/(h_cgs**2))**(1.5d0)

  	! helpers (よく使う分数)
  	real*8, parameter :: fivethirds = 5.0d0/3.0d0 	! 5/3 (非相対論的 断熱指数 γ)
  	real*8, parameter :: fourthirds = 4.0d0/3.0d0 	! 4/3 (相対論的 断熱指数 γ)
  	real*8, parameter :: overthree  = 1.0d0/3.0d0 	! 1/3

end module physical_constants