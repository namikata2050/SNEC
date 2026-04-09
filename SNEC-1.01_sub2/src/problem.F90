! subroutine problem: シミュレーションの「初期問題」を設定するルーチン
!
! program snec (メインプログラム) からシミュレーション開始時に一度だけ呼び出される。
! シミュレーションを開始する(t=0)ために必要な、すべての初期状態を設定する。
! (例: グリッドの構築、星の初期モデルの読み込み、爆発エネルギーの準備、
!      物理テーブルの構築など)
subroutine problem

! ===== モジュールの読み込み =====
! blmod: グローバル変数を読み込む。これから「初期値」を設定する配列(mass, r, rho, temp など)が主。
  use blmod, only: ntstart, tstart, ncomps, mass, r, Rstar, &
  	 	 	 	 opacity_floor, metallicity, envelope_metallicity, do_piston, &
  	 	 	 	 do_bomb, total_initial_energy, bomb_total_energy, bomb_spread, &
  	 	 	 	 rho, kappa, kappa_table, dkappadt, tau, temp, p, delta_mass, &
  	 	 	 	 lambda, inv_kappa, lum, eos_gamma1
! parameters: ユーザーが設定したパラメータ(imax, eoskey, profile_nameなど)を読み込む。
  use parameters
! eosmodule: 状態方程式(EOS)関連のモジュール(init_ionpotなど)。
  use eosmodule
! physical_constants: 物理定数(msun, rsunなど)を読み込む。
  use physical_constants
  implicit none

! ===== ローカル変数の宣言 =====
  real*8 :: buffer(imax) ! プロファイル読み込み時などに使う一時的な作業用配列
  integer :: i
  character(len=256) :: filename ! 出力ファイル名用の文字列変数

! for OPAL interpolation routine
! (OPALオパシティルーチン(古いF77コード)とデータをやり取りするための COMMON ブロック)
  real*4 :: opact,dopact,dopacr,dopactd
  common/e/ opact,dopact,dopacr,dopactd

!------------------------------------------------------------------------------
!
! ===== ステップ 1: 時間とステップ数の初期化 =====
!
  ntstart  	= 0 	  ! 開始ステップ番号を 0 に設定
  tstart  	= 0.0d0 ! 開始時刻を 0.0 秒に設定

!****************************** EOS setup *************************************
!
! ===== ステップ 2: 状態方程式(EOS)の選択と設定 =====
! input_parserで読み込んだ `eoskey` に基づいて、使用するEOSを最終決定する。
!
  if(eoskey.eq.1) then
  	  ! (1 = 理想気体EOS)
  	  write(6,*) "Using the ideal EOS."
  	  write(6,*) "WARNING: this EOS does not return the radiation pressure."
  	  ! Sedovテスト問題用の断熱指数 γ = 1.4 を設定
  	  eos_gamma1 = 1.4d0 !for Sedov
  else if(eoskey.eq.2) then
  	  ! (2 = Paczynski EOS。輻射圧、電離、縮退を考慮する、より現実的なEOS)
  	  write(6,*) "Using the Paczynski EOS!"
  else
  	  ! (上記以外 = 未実装)
  	  stop "Choice of EOS not available, check the parameter eoskey."
  endif

!********************* Allocate and initialize variables **********************
!
! ===== ステップ 3: 変数配列のメモリ確保と初期化 =====
!
! 化学組成プロファイルファイルを一度スキャンし、元素(同位体)の数 `ncomps` を取得する。
  call get_ncomps_from_profile(composition_profile_name,ncomps)

! `imax` (グリッド数) と `ncomps` (元素数) に基づいて、
! blmod で宣言された全ての動的配列 (rho(:), temp(:), comp(:,:) など) のメモリを確保(allocate)する。
  call allocate_vars

! 確保した全ての配列の値を 0.0d0 などで初期化する。
  call initialize_vars

!***************************** Grid setup *************************************
!
! ===== ステップ 4: グリッドの構築 =====
! シミュレーションで使用するラグランジュ座標（質量メッシュ）を構築する。
!
! プロファイルファイルを読み、モデルの最内殻と最外殻の「質量座標」を取得する。
  call get_inner_outer_mass_from_profile(profile_name,mass(1),mass(imax))

! もし質量切り捨て(mass_excision)が .true. (パラメータで指定) なら、
! 最内殻の質量座標 mass(1) を、パラメータで指定された値(mass_excised)に上書きする。
  if(mass_excision) then
  	  mass(1) = mass_excised*msun
  endif

! プロファイルファイルを読み、最内殻と最外殻の「半径」を取得する。
! (注: この半径はまだ仮のもので、後の `integrate_radius_initial` で再計算される)
  call get_inner_outer_radius_from_profile(profile_name,mass(1),r(1),r(imax))

! 最外殻の半径を Rstar (星の半径)として保存する。
  Rstar=r(imax)

! info.dat ファイルに、モデルの総質量と初期半径を記録する。
  open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
  	   status="unknown",form='formatted',position="append")
  write(666,*) 'Mass of the model = ', mass(imax)/msun, 'solar masses'
  write(666,*) 'Initial radius = ', Rstar/rsun, 'solar radii'
  close(666)

  ! set up the grid
! 'grid' サブルーチンを呼び出し、mass(1) と mass(imax) の間に
! imax 個のラグランジュ・グリッド（質量殻）を配置する。
! (ここで mass(:) 配列の 1 から imax までの全ての値が設定される)
  call grid

!*************************** Read the profile *********************************
!
! ===== ステップ 5: 初期モデル(プロファイル)の読み込み =====
!
  write(*,*) "Profile file: ",trim(profile_name)
! 'read_profile' サブルーチンを呼び出す。
! ステップ4で作成した質量グリッド(mass(:))上に、プロファイルファイルから
! 物理量（密度、温度、化学組成など）を補間し、
! rho(:), temp(:), comp(:,:) などの配列にマッピング(格納)する。
  call read_profile(profile_name)

  ! set up radius coordinates based on mass and density
! 'integrate_radius_initial' サブルーチンを呼び出す。
! 読み込んだ密度 rho(:) と質量グリッド mass(:) から、
! 質量保存の式 (dm = 4π r^2 ρ dr) を積分し、各グリッドの「正しい」初期半径 r(:) を計算し、確定させる。
  call integrate_radius_initial

  ! ! ★★★ 修正箇所 ★★★
  ! ! もし「入力ファイルと同じグリッド(same_with_input)」を使うなら、
  ! ! 半径はファイルから読み込んだ正しい値を使いたいので、再計算(上書き)してはいけない。
  ! if (gridding .ne. 'same_with_input') then
  !     call integrate_radius_initial
  ! else
  !     write(*,*) "Skipping integrate_radius_initial because gridding is same_with_input."
  ! end if
  ! ! ★★★ 修正ここまで ★★★

!************************ Set up the opacity floor ****************************
!
! ===== ステップ 6: オパシティ下限値(Floor)の設定 =====
! 数値計算の安定性のため、オパシティが極端に低くならないように下限値を設定する。
!
  do i=1, imax
  	! 各グリッドの金属量(metallicity)に基づいて、
  	! パラメータで指定されたコア用(of_core)とエンベロープ用(of_env)の下限値を線形補間し、
  	! opacity_floor(i) 配列に格納する。
  	opacity_floor(i) = (envelope_metallicity*of_core - of_env  	 &
  	 	 	 	 	 + metallicity(i)*(of_env - of_core))/(envelope_metallicity - 1)
  end do

  ! 設定したオパシティ下限値のプロファイルをデバッグ用にファイル出力する。
  filename = trim(adjustl(outdir))//"/opacity_floor.dat"
  call output_screenshot(opacity_floor,filename,imax)


!******************** Set up the energy of the thermal bomb *******************
!
! ===== ステップ 7: 爆発エネルギーの準備 =====
!
! 読み込んだ初期モデルの「総エネルギー」(内部E + 重力E + 運動E)を計算する。
  call conservation_compute_energies

! 'initial_data' の設定に応じて、爆発フラグを立てる。
  if(initial_data.eq."Piston_Explosion") then
  	  ! ピストン爆発フラグを立てる
  	  do_piston = .true.

  else if(initial_data.eq."Thermal_Bomb") then
  	  ! 熱爆弾フラグを立てる
  	  do_bomb = .true.

  	  ! --- 熱爆弾の総エネルギー `bomb_total_energy` を決定 ---
  	  if(bomb_mode.eq.1) then
  	 	  ! (モード1): `final_energy` を「シミュレーションの総エネルギーの目標値」とする。
  	 	  !            注入するエネルギー = (目標値) - (初期モデルの総エネルギー)
  	 	  bomb_total_energy = final_energy - total_initial_energy
  	  else if (bomb_mode.eq.2) then
  	 	  ! (モード2): `final_energy` の値を、そのまま「注入するエネルギー量」とする。
  	 	  bomb_total_energy = final_energy
  	  else
  	 	  write(6,"(A10,I3,A12)") "bomb_mode ",bomb_mode," not defined!"
  	 	  stop
  	  endif
  	  write(6,"(A60)") "***************************************************************************"
  	  write(6,"(A31,I3)") "Operating in Thermal Bomb Mode", bomb_mode

  	  ! --- 爆弾エネルギーを広げるグリッド数 `bomb_spread` を決定 ---
  	  do i=1, imax
  	 	  ! パラメータで指定された開始点(bomb_start_point)の質量 + 広げる質量(bomb_mass_spread)
  	 	  ! を超える最初のグリッド i を見つける。
  	 	  if(mass(i).ge.(mass(bomb_start_point)+bomb_mass_spread*msun)) then
  	 	 	  ! 注入開始点から i までのグリッド数を `bomb_spread` に設定
  	 	 	  bomb_spread = i - bomb_start_point
  	 	 	  ! ログファイルに記録
  	 	 	  open(unit=666, &
  	 	 	 	 	 file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
  	 	 	 	 	 status="unknown",form='formatted',position="append")
  	 	 	  write(666,*) 'Bomb energy is spread over = ', bomb_spread, 'points'
  	 	 	  write(6,*) 'Bomb energy is spread over = ', bomb_spread, 'points'
  	 	 	  close(666)
  	 	 	  exit
  	 	  end if
  	  end do

  	  ! 計算された初期エネルギーと爆弾エネルギーをログファイルに記録する。
  	  open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
  	 	   status="unknown",form='formatted',position="append")
  	  write(666,*) 'Total energy of the model = ', total_initial_energy, ' ergs'
  	  write(666,*) 'Total energy of the bomb = ', bomb_total_energy, ' ergs'
  	  write(6,*) 'Total energy of the model = ', total_initial_energy, ' ergs'
  	  write(6,*) 'Total energy of the bomb = ', bomb_total_energy, ' ergs'
  	  close(666)

  else
  	  ! "Piston_Explosion" でも "Thermal_Bomb" でもない場合
  	  ! stop "Wrong type of explosion, check the parameter 'initial_data'"
      write(6,*) "No additoinal energy injection"
  endif
  write(6,"(A60)") "***************************************************************************"


!****************** initialize some vairables *********************************
!
! ===== ステップ 8: 輻射・観測関連の初期化 =====
!
  ! call compose_opacity_tables_OPAL (古いコード、コメントアウトされている)

  ! もし輻射輸送をオン(radiation = .true.)にするなら、
  if(radiation) then
  	  ! 'compose_opacity_tables_OPAL' を呼び出し、
  	  ! OPALライブラリを使ってオパシティの「早見表」をメモリ上に構築する。
  	  ! (注: これはシミュレーションのセットアップで最も時間がかかる処理の一つ)
  	  call compose_opacity_tables_OPAL
  endif

  ! 'opacity' を呼び出し、t=0 の初期オパシティ kappa(:) をテーブルから計算する。
  call opacity(rho(:),temp(:),kappa(:),kappa_table(:),dkappadt(:))

  ! 'optical_depth' を呼び出し、t=0 の初期光学的深さ tau(:) を計算する。
  call optical_depth(rho(:), r(:), kappa_table(:), tau(:))

  ! 'luminosity' を呼び出し、t=0 の初期光度 lum(:) を計算する。
  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))

  ! 'read_BolCorr' を呼び出し、等級計算に使う「ボロメトリック補正」の
  ! テーブルファイル (BolCorr.datなど) を読み込む。
  call read_BolCorr

!*********** output the initial values of some variables for analysis *********
!
! ===== ステップ 9: 初期状態のファイル出力 (デバッグ・解析用) =====
! 'output_screenshot' は、配列データを単純なテキストファイルに出力するヘルパー関数。
!
  filename = trim(adjustl(outdir))//"/rho_initial.dat"
  call output_screenshot(rho,filename,imax) ! 初期密度

  filename = trim(adjustl(outdir))//"/rad_initial.dat"
  call output_screenshot(r,filename,imax) ! 初期半径

  filename = trim(adjustl(outdir))//"/mass_initial.dat"
  call output_screenshot(mass,filename,imax) ! 質量グリッド

  filename = trim(adjustl(outdir))//"/delta_mass_initial.dat"
  call output_screenshot(delta_mass,filename,imax) ! 各ゾーンの質量

  filename = trim(adjustl(outdir))//"/press_initial.dat"
  call output_screenshot(p,filename,imax) ! 初期圧力

  ! density as a function of distance from the surface inwards
  ! (表面からの距離 vs 密度 のプロファイルを出力)
  filename = trim(adjustl(outdir))//"/density_profile.dat"
  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
  	 	form='formatted',position="append")
  do i=1, imax-1
  	  ! (表面からの距離), (密度), (表面からの質量)
  	  write(666,"(1P20E19.10E3)") r(imax)-r(imax-i), &
  	 	  rho(imax-i),mass(imax)-mass(imax-i)
  enddo
  close(666)

end subroutine problem