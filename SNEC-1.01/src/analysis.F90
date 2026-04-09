! subroutine analysis: 派生的な物理量や観測量を計算する。
!
! blstep の中で物理計算(hydro_rad)が終わった後に呼び出される。
! 更新されたシミュレーションの「生」の状態(rho, temp, r など)を元に、
! 観測される可能性のある「派生的な量」を計算・分析(Analysis)する。
! (例: 光球はどこか？ 光度や等級はいくつか？)
subroutine analysis

! ===== モジュールの読み込み =====
! blmod: ほぼ全てのグローバル変数(物理状態)を読み込む
  use blmod
! parameters: outdir(出力ディレクトリ名)などの設定を読み込む
  use parameters
! physical_constants: 物理定数(光速 clite, シュテファン・ボルツマン定数 sigma_SB など)を読み込む
  use physical_constants
  implicit none

! ===== ローカル変数の宣言 =====
  character(len=1024) :: filename 	! ファイル出力用の文字列変数
  integer :: i 						! ループカウンター

  real*8 :: T_eff_for_BC 				! ボロメトリック補正(BC)に使う有効温度 (下限値適用後)
  real*8 :: bol_corr_used(11) 			! 補間されたボロメトリック補正の値 (11バンド分)
  ! 有効温度 T_eff の下限値。BCテーブルがこれより低い温度を持っていないため。
  real*8, parameter :: T_eff_min = 5000.0d0

!KM
  real*8 :: mag_col(11)
  real*8 :: T_col_for_BC
  integer :: n_lum

!------------------------------------------------------------------------------
!
! ===== セクション 1: 光学的深さの計算と光球(Photosphere)の特定 =====
!
! 星の「表面」である光球(Photosphere)の位置を特定する。
! 光球は「光学的深さ τ = 2/3」になる場所として定義される。
!
!---------- Calculating optical depth and tracing the photosphere -------------

  ! kappa_table (without the opacity floor) is used to trace the photosphere
  ! (オパシティの下限値(floor)の影響を受けない kappa_table を使って光球を探す)
  ! 'optical_depth' サブルーチンを呼び出し、
  ! 星の表面(imax)から中心に向かって光学的深さ τ = ∫κρdr を計算し、
  ! 結果をグローバル変数 tau(:) 配列に格納する。
  ! call optical_depth(rho(:), r(:), kappa_table(:), tau(:))

  !KM
  if(csm_ionization) then
    call optical_depth(rho(:), r(:), kappa(:), tau(:))
  else
  !kappa_table (without the opacity floor) is used to trace the photosphere
    call optical_depth(rho(:), r(:), kappa_table(:), tau(:))
  end if
!

  ! find the grid point, where the photosphere is located
  ! (光球が存在するグリッド地点を探す)
  ! 星の外側(imax-1)から内側(1)に向かってループ
  do i=imax-1, 1, -1
  	  ! tau(i) が初めて 2/3 (0.66d0) を超えた場所を見つける
  	  if(tau(i).gt.0.66d0) then
  	 	  ! その1つ外側のグリッド(i+1)を光球の位置(index_photo)として記録する
  	 	  index_photo = i + 1
  	 	  ! 見つかったらループを抜ける
  	 	  exit
  	  endif
  enddo

  ! fix the moment, when the photosphere reaches the inner boundary, if it does
  ! (もし光球が内部境界(中心)に達した場合の処理)
  ! 中心(i=1)まで積分しても τ < 2/3 なら、星全体がほぼ透明になったことを意味する。
  if(tau(1).lt.0.66d0) then
  	  ! 光球の位置を中心(1)とする
  	  index_photo = 1
  	  ! photosphere_fell_on_the_center フラグがまだ 0 なら (＝最初に到達した瞬間なら)
  	  if(photosphere_fell_on_the_center.eq.0) then
  	 	  ! フラグを 1 に立てる
  	 	  photosphere_fell_on_the_center = 1
  	 	  ! info.dat ファイルに「光球が中心に到達した」時刻を記録する
  	 	  open(unit=666, &
  	 	 	 	 file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
  	 	 	 	 status="unknown",form='formatted',position="append")
  	 	  write(666,*) 'Photosphere reached the center at ', time, 'seconds'
  	 	  close(666)
  	  end if
  end if

! ===== セクション 2: 光球での物理量の計算 (補間) =====
! index_photo はおおよその位置。τ = 2/3 の「正確な」位置での物理量を補間で求める。

  ! find the values of some variables at the photosphere
  ! (光球でのいくつかの変数の値を求める)
  ! 光球がまだ中心に落ちていない場合
  if(photosphere_fell_on_the_center.eq.0) then
  	  ! 'map_map' (カスタム補間ルーチン) を呼び出す。
  	  ! 「tau 配列をX軸、lum 配列をY軸として、X=0.66 の時のYの値を求めよ」
  	  ! (配列を逆順(imax:1:-1)で渡しているのは、tau が外側で 0 に近いため)
  	  call map_map(lum_photo,  0.66d0, lum(imax:1:-1),  tau(imax:1:-1),imax)
  	  call map_map(mass_photo, 0.66d0, mass(imax:1:-1), tau(imax:1:-1),imax)
  	  call map_map(vel_photo,  0.66d0, vel(imax:1:-1),  tau(imax:1:-1),imax)
  	  call map_map(rad_photo,  0.66d0, r(imax:1:-1),  	 tau(imax:1:-1),imax)
  ! 星全体が透明になった場合
  else
  	  ! 中心(i=1)の値をそのまま光球の値とする
  	  lum_photo   = lum(1)
  	  mass_photo  = mass(1)
  	  vel_photo   = vel(1)
  	  rad_photo   = r(1)
  end if

  ! photosphere tracer is equal to 1 at the photosphere, and 0 everywhere else
  ! used for visualization
  ! (光球トレーサー配列を作成。光球の位置(index_photo)だけ 1 にする。可視化用。)
  photosphere_tracer(:) = 0.0d0
  photosphere_tracer(index_photo) = 1.0d0

  ! --- オパシティに関する警告フラグ ---
  ! check if the photosphere moves through the regions with wrong opacity
  ! log10(T) = 3.75 is the lower boundary of the OPAL opacity tables
  ! (光球が、OPALテーブルの範囲外(低温)かつ金属量が異なる領域にあるかチェック)
  ! (kappa_table > kappa は、opacity_floor が適用されたことを意味する)
  if(metallicity(index_photo).gt.envelope_metallicity .and. &
  	 	log10(temp(index_photo)).lt.3.75d0 .and. &
  	 	kappa_table(index_photo).gt.kappa(index_photo) &
  	 	.and. index_photo.gt.1 ) then
  	  opacity_corrupted = 1 ! 警告フラグを立てる
  else
  	  opacity_corrupted = 0
  end if

!------------------------------------------------------------------------------
!
! ===== セクション 3: DIAGNOSTICS - その他の分析量の計算 =====
!
!------------------------ Tracing the luminosity shell (輝度シェルの追跡) --------
! 輝度シェル(Luminosity Shell)またはトラッピング半径(Trapping Radius)を探す。
! これは、光子がガスの膨張から「逃げ出せる」ようになる境界 (τ ≈ c/v) を示す。
  index_lumshell = imax
  ! 外側から内側へループ
  do i=1, imax - 1
  	  ! tau が初めて c/v を超えた場所を見つける
  	  if(tau(imax-i).gt.(clite/vel(imax-i))) then
  	 	  index_lumshell = imax - i + 1
  	 	  exit
  	  end if
  end do
  ! もし中心近くまで τ < c/v なら、輝度シェルは中心にあるとする
  if (tau(2).le.(clite/vel(2))) then
  	  index_lumshell = 1
  end if
  ! 輝度シェルの質量座標を記録
  mass_lumshell = mass(index_lumshell)
  
  ! characteristic diffusion and expansion times for different shells
  ! (各シェルにおける、特徴的なタイムスケールを計算する)
  do i=1, imax-1
  	  ! 拡散時間 t_diff ≈ κρR^2 / c (光子が表面まで拡散する時間)
  	  time_diff(i) = kappa(i)*rho(i)*(r(imax)-r(i))**2.0/clite
  	  ! 膨張時間 t_exp ≈ R / v (ガスが表面まで膨張する時間)
  	  time_exp(i) = (r(imax)-r(i))/(vel(imax)-vel(i))
  end do

  ! internal energies of shells from the given radius out to the surface
  ! (各シェル i より「外側」のガスの総内部エネルギー E_shell(i) を計算する)
  do i=1, imax
  	  ! E_shell(i) = Σ[k=i to imax] (eps(k) * delta_mass(k))
  	  E_shell(i) = sum(eps(i:imax)*delta_mass(i:imax))
  end do

!------------------------------------------------------------------------------
!
! ===== セクション 4: 観測光度(Observed Luminosity)の計算 =====
!
! 観測される光度は、光球から漏れ出す熱(lum_photo)だけでなく、
! 光球より「外側」で起こるNi崩壊のガンマ線も直接寄与する。
!
!------------------- Calculate the observed luminosity ------------------------

  ! ! observed luminosity is the sum of lum_photosphere and Ni contribution
  ! ! (観測光度 = 光球からの光度 + Ni崩壊からの寄与)
  ! ! 光球がまだ中心に落ちていない場合
  ! if(photosphere_fell_on_the_center.eq.0) then
  ! 	  ! lum_observed = 光球光度 + Σ[k=光球から外側] (Ni加熱率 * 質量)
  ! 	  lum_observed = lum_photo + sum(Ni_energy_rate* &
  ! 	 	 	 Ni_deposit_function(index_photo:imax)*delta_mass(index_photo:imax))
  ! ! 星全体が透明になった場合
  ! else
  ! 	  ! lum_observed = 中心光度 + Σ[全領域] (Ni加熱率 * 質量) (＝Ni総光度)
  ! 	  lum_observed = lum(1) + &
  ! 	 	 	 sum(Ni_energy_rate*Ni_deposit_function(1:imax)*delta_mass(1:imax))
  ! end if

  if(rad_photo.gt.r(shockpos+10)) then
     lum_observed=lum_photo + sum(Ni_energy_rate* &
          Ni_deposit_function(index_photo:imax)*delta_mass(index_photo:imax))
  else
     n_lum=shockpos+10
     n_lum=min(n_lum,imax)
     lum_observed=lum(n_lum) + sum(Ni_energy_rate* &
          Ni_deposit_function(n_lum:imax)*delta_mass(n_lum:imax))
  end if
! KM 2018.10.18

  ! write down the time when the contribution of the Ni above the
  ! photosphere to the luminosity is greater than 5%
  ! (光球より外側のNi崩壊の寄与が、光球光度の5%を超えた瞬間を記録する)
  ! (shockpos_stop=1 は、衝撃波が星を突き抜けた後であることを確認)
  if(shockpos_stop.eq.1 .and. &
  	 	abs((lum_observed - lum_photo)/lum_photo).gt.0.05 .and. &
  	 	Ni_contributes_five_percents.eq.0) then
  	 
  	  ! フラグを立てる (一度だけ記録するため)
  	  Ni_contributes_five_percents = 1
  	  ! info.dat ファイルに時刻を記録する
  	  open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
  	 	 	 status="unknown",form='formatted',position="append")
  	  write(666,*) 'Ni contribution to the luminosity is 5% at ', time, 'seconds'
  	  close(666)

  end if

!------------------------------------------------------------------------------
!
! ===== セクション 5: 天文観測量(等級 Magnitude)への変換 =====
!
!-------- Calculation of color magnitudes using bolometric corrections --------

  ! 有効温度 T_eff を計算する
  ! シュテファン=ボルツマンの法則 L = 4π R^2 σ T_eff^4 を T_eff について解く
  T_eff = (lum_photo/(4.0d0*pi*sigma_SB*rad_photo**2))**0.25d0
  
  T_col = (lum_observed/(4.0d0*pi*sigma_SB*radshock**2))**0.25d0

  ! see Eq.(3) of Swartz et al., ApJ 374:266 (1991) and explanation there
  ! (BCテーブルの下限に合わせて、T_eff に下限値 T_eff_min (例: 5000K) を適用する)
  T_eff_for_BC = MAX(T_eff,T_eff_min)

  T_col_for_BC = MAX(T_col,T_eff_min)

  ! here in cases, when the effective temperature goes beyond the boundaries
  ! of the table BolCorr.dat, the linear extrapolation is used
  ! (BCテーブルの範囲外になった場合、map_map は線形外挿を行う)
  ! 11種類のフィルターバンドについてループ
  do i=1, 11
  	  ! 'map_map' を使い、T_eff_for_BC に対応するボロメトリック補正(BC)の値を
  	  ! テーブル(bol_corr, temp_bol_corr)から補間して求める
  	  call map_map(bol_corr_used(i),T_eff_for_BC,bol_corr(1:nlines_bol_corr,i), &
  	 	 	 temp_bol_corr,nlines_bol_corr)
  	  ! 絶対等級 M を計算する
  	  ! M = M_sun - BC - 2.5 * log10(L / L_sun)
  	  ! magnitudes(i) = sun_mag - bol_corr_used(i) - 2.5d0*log10(lum_photo/sun_lum) !default
      magnitudes(i) = sun_mag - 2.5d0*log10(lum_observed/sun_lum)
  end do

  ! write down the magnitudes in a file
  ! (等級をファイルに書き出す)
  ! スカラ量出力のタイミングかどうかをチェック
  if(time.eq.0.0d0.or.time.gt.tdump_scalar) then
  	  ! magnitudes.dat ファイルを追記モード('append')で開く
  	  open(666,file=trim(adjustl(outdir))//"/magnitudes.dat",&
  	 	 	 status='unknown',position='append')
  	  ! 時刻、有効温度、11バンドの等級を1行に書き込む
  	  write(666,"(13E18.9)") time, T_eff_for_BC, magnitudes(1:11)
  	  close(666)
  endif

end subroutine analysis