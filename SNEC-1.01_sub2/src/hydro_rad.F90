! subroutine hydro_rad: 輻射流体力学(Radiation-Hydrodynamics)ソルバー
!
! blstepから呼び出される、シミュレーションの「心臓部」。
! 物理方程式を解き、状態を 1 タイムステップ (t -> t + dtime) 進める。
!
! 手法:
! 1. 陽解法(Explicit): 運動方程式(速度)、半径、密度を更新する。
! 2. 陰解法(Implicit): エネルギー方程式と輻射輸送を Newton-Raphson 法で連立して解き、
!    温度、内部エネルギー、圧力を更新する。
subroutine hydro_rad

! ===== モジュールの読み込み =====
! blmod: ほぼ全てのグローバル変数(物理状態)を読み込む
  use blmod
! parameters: imax(グリッド数)などの設定を読み込む
  use parameters
! physical_constants: 物理定数を読み込む
  use physical_constants
  implicit none

! ===== ローカル変数の宣言 =====
  integer :: i,k 		! ループカウンター
  integer :: keytemp,keyerr 	! EOSサブルーチンに渡す制御キー
  real*8 :: dtv 		! 時間中心化された時間刻み幅 (dtime + dtime_p)/2

  ! --- 陰解法(Newton-Raphson)ソルバー用の変数 ---
  ! Variables for inversion of Jacobian (ヤコビアン行列の逆行列計算用)
  external :: dgbsv 	! LAPACK(線形代数ライブラリ)の「帯行列ソルバー」を使う宣言
  integer :: info 		! LAPACKソルバーが返すエラーコード
  ! kl, ku: 行列のバンド幅。kl=1, ku=1 は「三重対角行列(Tri-diagonal)」を意味する。
  integer, parameter :: kl = 1
  integer, parameter :: ku = 1
  integer, parameter :: ldab=2*kl+ku+1 ! LAPACKに渡す行列配列のサイズ
  ! ab: ヤコビアン行列 J をLAPACKの特殊な形式で格納する配列
  ! b:  右辺ベクトル D (残差) を格納。ソルバー実行後は解 (δT) が格納される。
  real*8 :: ab(ldab,imax-1), b(imax-1)
  integer :: ipiv(imax-1) 	! LAPACKがピボット操作に使う作業用配列
  
  real*8 :: delta_max 		! 収束判定用: 最大の相対補正量 |δT / T|
  integer :: location_max 	! 最大補正が起きたグリッドのインデックス
  
  ! A, B, C: 三重対角行列 J の (A=上, B=対角, C=下) 成分を一時的に格納する配列
  ! D:       右辺ベクトル D (残差) を一時的に格納する配列
  real*8 :: Aarray(imax-1), Barray(imax-1), Carray(imax-1), Darray(imax-1)
  ! ..._temp: Newton-Raphson の反復計算中に使う「現在の推測値」を格納する作業用配列
  real*8 :: p_temp(imax), eps_temp(imax), lum_temp(imax), temp_temp(imax)
  real*8 :: lambda_temp(imax)

  ! --- 陰解法ソルバーの収束判定パラメータ ---
  real*8, parameter :: EPSTOL = 1.0d-7 	! 許容する相対誤差 (Epsilon Tolerance)
  integer, parameter :: ITMAX = 300 	! 最大反復回数 (Iteration Maximum)


!------------------------------------------------------------------------------
!
! ===== ステップ 0: 状態のバックアップと時間刻み幅の定義 =====
!
  ! follow dt convention of MB93: (Mihalas & Burrows 1993 などの慣習に従う)
  ! 時間中心化された時間刻み幅 (dtv) を計算する。 (dtime + dtime_p)/2
  dtv = 0.5d0 * (dtime + dtime_p)

  ! copy over data into _p arrays
  ! 現在の状態(t=n)を「1ステップ前(_p)」の配列にコピー(バックアップ)する。
  ! これにより、rho, vel, r などの配列を「新しい状態(t=n+1)」として
  ! 安心して上書き更新できるようになる。
  rho_p(:) = rho(:)
  vel_p(:) = vel(:)
  r_p(:)   = r(:)
  cr_p(:)  = cr(:)
  eps_p(:) = eps(:)
  p_p(:) = p(:)
  temp_p(:) = temp(:)
  kappa_p(:) = kappa(:)


!------------------------------------------------------------------------------
!
! ===== ステップ 1 (陽解法 / Explicit): 流体力学の更新 =====
!
!---------------------------- update velocities (速度の更新) -------------------------------
! 運動方程式: v_new = v_old + a * dt
! (注意: ここで使う圧力 p や半径 r は 1ステップ前(t=n) の値)

  ! 内側境界の速度をピストンとして強制的に設定する
  if(do_piston .and. time.ge.piston_tstart .and. time.le.piston_tend) then
  	  vel(1) = piston_vel
  	  vel(2) = piston_vel
  endif

  ! 2番目からimax番目のグリッドまでループ
  do i=2,imax
  	! 新しい速度 vel(i) = 1ステップ前の速度 vel_p(i) + (各種の力) * dtv
  	vel(i) = vel_p(i) &

  	  ! gravity (重力による加速度項)
  	  - dtv * ggrav*mass(i) / r(i)**2 *gravity_switch   &

  	  ! pressure (圧力勾配による加速度項 (ラグランジュ座標系))
  	  - dtv * 4.0d0*pi*r(i)**2 * (p(i) - p(i-1)) / delta_cmass(i-1)   &

  	  ! artificial viscosity (人工粘性による加速度項)
  	  - dtv * 4.0d0*pi * (cr(i)**2 * Q(i) - cr(i-1)**2 * Q(i-1))/delta_cmass(i-1)
  enddo

  ! 内側境界(i=1)の速度を 0 に設定する (ピストン終了後、または爆弾の場合)
  if(do_piston.and.time.ge.piston_tend) then
  	vel(1) = 0.0d0
  else if(do_bomb) then
  	vel(1) = 0.0d0
  endif

!----------------------- update the radial coordinates (半径の更新) -------------------------
! r_new = r_old + v_new * dt
  do i=1,imax
  	! 新しい半径 r(i) = 1ステップ前の半径 r_p(i) + (さっき計算した新しい速度 vel(i)) * dtime
  	r(i) = r_p(i) + dtime * vel(i)
  	
  	! ★重要: グリッド交差(Grid Crossing)のチェック
  	! ラグランジュ法では、質量殻(グリッド)が追い越してはならない。
  	! もし外側の殻 i が内側の殻 i-1 より内側に来たら、計算は破綻している。
  	if(i.gt.1) then
  	 	 if (r(i).lt.r(i-1)) then
  	 	 	  write(*,*) 'radius of a gridpoint', i, 'is less than preceding'
  	 	 	  stop
  	 	 end if
  	end if
  enddo

!------------------------- update the zone densities (密度の更新) --------------------------
! rho_new = M / V_new
  do i=1,imax-1
  	  ! 新しい密度 rho(i) = ゾーンの質量 delta_mass(i) / (球殻の体積 V)
  	  ! 体積 V = (4/3) * pi * (r_outer^3 - r_inner^3)
  	  ! (さっき計算した新しい半径 r(i+1) と r(i) を使う)
  	  rho(i) = delta_mass(i) / (4.0d0*pi * (r(i+1)**3 - r(i)**3)/3.0d0)
  enddo
  rho(imax) = 0.0d0 !passive boundary condition (最外殻の密度は0)


!------------------------- update zone center radius (ゾーン中心半径の更新) --------------------------
  do i=1,imax-1
  	! ゾーン中心の半径 cr(i) を、体積平均で計算する
  	cr(i) = ( ( r(i)**3 + r(i+1)**3 ) / 2.0d0 )**(1.0d0/3.0d0)
  enddo
  cr(imax) = r(imax) + (r(imax) - cr(imax-1))
  !passive boundary condition, ... (最外殻のゾーン中心半径。人工粘性Qが0なので使われない)


  ! update the artificial viscosity
  ! 新しい rho, vel に基づいて、人工粘性 Q を再計算する
  call artificial_viscosity


!------------------------------------------------------------------------------
!
! ===== ステップ 2 (陰解法 / Implicit): 熱・輻射の更新 =====
!
!----------- update the temperature, pressure and internal energy -------------

  ! calculate heating term due to Ni
  ! (Ni崩壊による加熱率を計算/更新する)
  if(time.ge.time_Ni) then
  	  time_Ni = time_Ni + Ni_period
  	  call nickel_heating
  endif

  ! calculate heating term due to bomb
  ! (熱爆弾による加熱率を計算する)
  if(do_bomb .and. time.ge.bomb_tstart .and. time.le.bomb_tend) then
  	  call bomb_pattern
  else
  	  bomb_heating(:) = 0.0d0
  endif

  ! --- Newton-Raphson 法の準備 ---
  ! Initial guess for the quantities at the next time step
  ! (次のタイムステップ(t=n+1)の値の「初期推測値」として、
  !  現在(t=n)の値を作業用配列(..._temp)にコピーする)
  p_temp(1:imax) = p(1:imax)
  eps_temp(1:imax) = eps(1:imax)
  temp_temp(1:imax) = temp(1:imax)

  ! --- Newton-Raphson 反復ループ (最大 ITMAX 回) ---
  do k=1, ITMAX

  	! --- ステップ 2a: 状態方程式(EOS)の計算 ---
  	! 現在の「推測値」 (temp_temp, rho) を使って、
  	! 「推測値」に対応する圧力(p_temp)、エネルギー(eps_temp)、
  	! および陰解法ソルバーに必要な微分係数(dpdt, dedt)を計算する。
  	keytemp = 1
  	call eos(rho(1:imax-1),temp_temp(1:imax-1),ye(1:imax-1), &
  	 	 	abar(1:imax-1),p_temp(1:imax-1),eps_temp(1:imax-1), &
  	 	 	cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), entropy(1:imax-1), &
  	 	 	p_rad(1:imax-1),keyerr,keytemp,eoskey)

  	! --- ステップ 2b: 光度(Luminosity)の計算 ---
  	! 現在の「推測値」 (temp_temp) を使って、
  	! 「推測値」に対応する光度(lum_temp)を計算する。
  	call luminosity(r(:), temp_temp(:), kappa_p(:), &
  	 	 	 	 	lambda_temp(:), inv_kappa(:), lum_temp(:))

  	! --- ステップ 2c: ヤコビアン行列と残差ベクトルの作成 ---
  	! calculate the coefficients of the equation:
  	! A(i)*\delta T(i+1) + B(i)*\delta T(i) + C(i)*\delta T(i-1) = D(i)
  	! (エネルギー方程式を T について線形化し、
  	!  ヤコビアン行列 J (A, B, C 配列) と
  	!  残差ベクトル D (Darray 配列) を作成する)
  	call matrix_arrays(temp_temp(:), lambda_temp(:), inv_kappa(:), &
  	 	 	eps_temp(:), p_temp(:), lum_temp(:), &
  	 	 	Aarray(:), Barray(:), Carray(:), Darray(:))

  	 	 	 	
  	! --- ステップ 2d: 行列の組み立て ---
  	! assemble the matrix in the form used by lapack
  	! (A, B, C, D 配列を、LAPACKソルバー(dgbsv)が要求する
  	!  特殊な格納形式 (ab, b) に詰め替える)
  	ab(2,2:imax-1) = Aarray(1:imax-2) ! 上対角(Upper)
  	ab(3,1:imax-1) = Barray(1:imax-1) ! 主対角(Main)
  	ab(4,1:imax-2) = Carray(2:imax-1) ! 下対角(Lower)
 
  	b(1:imax-1) = Darray(1:imax-1) 	  ! 右辺ベクトル D (残差)

  	! --- ステップ 2e: 線形方程式 J * (δT) = D を解く ---
  	! invert the matrix (行列を解く)
  	! if the inversion fails, ... (もし失敗したら 'failed_matrix.dat' に書き出す)
  	info = 0
  	! LAPACKソルバーを呼び出す。
  	! 成功すると、b 配列の中身が D から 解(δT) (補正量) に置き換わる。
  	call dgbsv(imax-1,kl,ku,1,ab,ldab,ipiv,b,imax-1,info)

  	! LAPACK エラーチェック
  	if(info.ne.0) then
  	 	open(unit=666, &
  	 	 	 file=trim(adjustl(trim(adjustl(outdir))//"/failed_matrix.dat")), &
  	 	 	 status="unknown",form='formatted',position="append")
  	 	do i=1,imax-1
  	 	 	write(666,*) ab(2,i), ab(3,i), ab(4,i), Darray(i)
  	 	enddo
  	 	close(666)
  	 	stop "problem in the matrix inversion (see Data/failed_matrix.dat)"
  	endif

  	! --- ステップ 2f: 収束判定 ---
  	! check if the iteration procedure converged
  	! (反復計算が収束したかチェックする)
  	delta_max = 0.0d0
  	! 最大の「相対」補正量 |δT / T| を探す
  	do i=1,imax-1
  	 	 if(abs(b(i)/temp_temp(i)).gt.delta_max) then
  	 	 	  delta_max = abs(b(i)/temp_temp(i))
  	 	 	  location_max = i
  	 	 endif
  	enddo
  	! もし最大の相対補正量が許容誤差(EPSTOL)以下なら、収束したとみなし、
  	! 成功ラベル(101)へジャンプする。
  	if(delta_max.le.EPSTOL) goto 101

  	! --- ステップ 2g: 推測値の更新 (未収束の場合) ---
  	! add the increment to the temperature
  	! (温度に補正量(δT)を加える)
  	do i=1, imax-1
  	 	 ! T_new_guess = T_old_guess + δT
  	 	 ! (b(i) には δT が格納されている)
  	 	 temp_temp(i) = temp_temp(i) + b(i)
  	 	 ! 安全装置: もし温度が負になったら、計算が発散したとみなし、
  	 	 ! 失敗ラベル(100)へジャンプする。
  	 	 if(temp_temp(i).lt.0.0d0) then
  	 	 	  goto 100
  	 	 end if
  	end do

  enddo ! Newton-Raphson ループ (k=1, ITMAX) の終端

  ! ----- ループ終了後の処理 (失敗) -----
  ! (ITMAX 回反復しても収束しなかった場合、ここに到達する)
  100 continue

  ! 失敗を画面に表示
  write(6,*) "EOS problem", delta_max, location_max
  ! ★重要: 失敗フラグ `scratch_step` を .true. にセットする。
  ! これにより、blstep はこのステップを「失敗」と認識し、
  ! dtime を半分にしてやり直す。
  scratch_step = .true.

  ! ----- ループ終了後の処理 (成功) -----
  ! (goto 101 でここにジャンプしてきた場合)
  101 continue

  	  
  ! ===== ステップ 3: 最終状態の確定 (成功した場合) =====
  !
  ! 収束した「推測値」(..._temp)を、
  ! 正式な「次のステップ(t=n+1)の状態」としてグローバル変数(eps, p, temp)にコピーする。
  eps(1:imax-1) = eps_temp(1:imax-1)
  p(1:imax-1)   = p_temp(1:imax-1)
  temp(1:imax-1) = temp_temp(1:imax-1)


  ! passive boundary conditions, do not participate in the evolution
  ! (最外殻の境界条件を設定)
  temp(imax) = 0.0d0
  eps(imax) = 0.0d0

  ! active boundary condition, used in the velocity update
  ! (最外殻の圧力は0。これは次のステップの vel(imax) の計算で使われる)
  p(imax) = 0.0d0

  ! ===== ステップ 4: 最終的な物理量の更新 =====
  ! 
  ! 確定した「新しい温度」temp_temp(:) を使って、
  ! このステップの最終的な「オパシティ」kappa(:) を計算・更新する。
  call opacity(rho(:),temp_temp(:),kappa(:),kappa_table(:),dkappadt(:))
  	
  ! 確定した「新しい温度」temp(:) と「新しいオパシティ」kappa(:) を使って、
  ! このステップの最終的な「光度」lum(:) を計算・更新する。
  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))


end subroutine hydro_rad