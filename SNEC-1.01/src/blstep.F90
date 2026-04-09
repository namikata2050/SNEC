! subroutine blstep: シミュレーションの「1ステップ」を実行する心臓部。
! program snec のメインループから呼び出され、
! 物理状態を t から t + dtime へと進める役割を担います。
!
! 計算が失敗した場合に備えて、ステップをやり直す「リトライ機能」も備えています。
subroutine blstep

! ===== モジュールの読み込み =====
! blmod: 主要な物理量配列(rho, temp, vel...)を読み込む。
!        特に重要なのが `scratch_step`。これは計算が失敗したことを示すフラグ。
  use blmod, only: rho, temp, ye, abar, eps, p, cs2, vel, r, cr, scratch_step, &
  	 	 	 	 time, dtime, shockpos_stop, nt
! parameters: imax (グリッド数) や radiation (輻射計算フラグ) などの設定を読み込む。
  use parameters
! physical_constants: 物理定数を読み込む。
  use physical_constants
  implicit none

! ===== ローカル変数の宣言 =====
! これらは、計算に失敗した場合に備えて、1ステップ前(t時点)の状態を
! 一時的に保存(バックアップ)するための「_save」配列です。
! imax (グリッド総数) と同じサイズで宣言されます。
  real*8 :: rho_save(imax),temp_save(imax),vel_save(imax)
  real*8 :: p_save(imax),ye_save(imax),cs2_save(imax), abar_save(imax)
  real*8 :: dedt_save(imax), r_save(imax), cr_save(imax)
  real*8 :: eps_save(imax)

! リトライ(やり直し)の回数を数えるカウンター
  integer :: iterations

!------------------------------------------------------------------------------

! ===== ステップ 1: 計算前の「状態保存」 (セーブポイントの作成) =====
  ! save old values in case we have to redo the step
! 現在(t時点)のシミュレーション状態(rho, tempなど)を、
! 上で宣言したバックアップ用のローカル配列(rho_saveなど)に丸ごとコピーします。
! これにより、この後の計算が失敗しても、この時点の状態に安全に戻ることができます。
  rho_save   = rho
  temp_save  = temp
  ye_save  	 = ye
  abar_save  = abar
  eps_save  	 = eps
  p_save  	 	= p
  cs2_save  	= cs2
  vel_save  = vel
  r_save  	 	= r
  cr_save  	 = cr

! リトライ回数を0に初期化
  iterations = 0

! ===== ステップ 2: メイン計算ループ (リトライ機能付き) =====
!
! 'scratch_step' が .true. (＝計算が失敗した)
! または 'iterations.eq.0' (＝最初の1回目の挑戦)
! である限り、この do-while ループを繰り返します。
  do while(scratch_step.or.(iterations.eq.0))

! ループの繰り返し回数を1増やす
  	  iterations = iterations + 1

! ----- リカバリー処理 (もし `scratch_step` が .true. なら) -----
! (ループの2回目以降、つまり計算が失敗してやり直す場合にのみ実行される)
  	  if(scratch_step) then
  	 	 ! redo step with half the timestep; this is 
  	 	 ! sometimes necessary if the solver does not converge
! 	 	 (ソルバーが収束しない場合、タイムステップを半分にしてやり直す)

! 	 	 安全装置: 10回やり直しても失敗する場合は、計算が破綻したとみなしプログラムを停止する。
  	 	 if (iterations.gt. 10) then
  	 	 	  stop "Stopping evolution. Time step repeated 10 times without luck"
  	 	 endif
  	 	 
! 	 	 ログ出力: 「ステップを破棄(scratch)してやり直す」ことを画面に表示
  	 	 write(6,*) "Scratching entire step...", nt
  	 	 write(6,"(A25,1P10E15.6)") "time,dt,dt_new: ",time,dtime,dtime/2.0d0
  	 	 
! 	 	 ★重要: タイムステップ(dtime)を半分にする。
! 	 	 これにより、より安定な(小さい)時間刻み幅で再挑戦する。
  	 	 dtime = dtime / 2.0d0
  	 	 
! 	 	 ★重要: 状態の復元(リストア)。
! 	 	 計算に失敗してゴミデータが入った可能性のあるグローバル変数(rhoなど)を、
! 	 	 バックアップ配列(rho_save)から書き戻し、t時点のクリーンな状態に戻す。
  	 	 rho = rho_save
  	 	 temp = temp_save
  	 	 ye = ye_save
  	 	 abar = abar_save
  	 	 eps = eps_save
  	 	 p = p_save
  	 	 cs2 = cs2_save
  	 	 vel = vel_save
  	 	 r = r_save
  	 	 cr = cr_save
  	 	 
! 	 	 やり直し処理が完了したので、フラグを .false. に戻す。
! 	 	 (もし次の `hydro_rad` の呼び出しで *また* 失敗すれば、`hydro_rad` 内部で再び .true. にセットされる)
  	 	 scratch_step = .false.
  	  endif

! ----- ステップ 3: 物理計算の実行 -----
  	  ! select between pure hydro or 
  	  ! radiation+hydro solver
! 	  (純粋な流体力学ソルバーか、輻射流体力学ソルバーかを選択)

! 	  'radiation' フラグ (parametersで設定) をチェック
  	  if(radiation) then
! 	 	 輻射 + 流体 を解くメインソルバーを呼び出す
  	 	 call hydro_rad
  	  else
! 	 	 輻射を無視し、流体力学だけを解くソルバーを呼び出す
  	 	 call hydro
  	  end if
! 	  (もし `hydro_rad` や `hydro` が収束に失敗したら、
! 	  それらの内部で `scratch_step = .true.` がセットされる)

! ----- ステップ 4: 計算後の解析処理 -----
! (物理計算が *成功* した後、新しい状態(t+dtime)に基づいて派生的な量を計算する)

! 爆発エネルギー注入が終わり(time >= bomb_tend)、
! かつ衝撃波追跡がまだ停止していなければ(shockpos_stop == 0)、
! 'shock_capture' を呼び出し、衝撃波の現在位置を特定する。
  	  if (time.ge.bomb_tend .and. shockpos_stop.eq.0) then
  	 	 call shock_capture
  	  endif

! Sedovテストモード(.true.) で *ない* 場合は、
! 'analysis' を呼び出し、光球の位置や光度、等級などを計算する。
  	  if(.not.sedov) then
  	 	 call analysis
  	  endif
  	 	 	 	 
! エネルギー保存則が守られているかをチェックするルーチンを呼び出す。
  	  if(time.gt.0.0d0) call conservation_compute_energies

! do-while ループの終端。
! もし `scratch_step` が .false. (＝計算が成功した) なら、
! ループの条件が .false. となり、ループを抜ける。
  enddo

! サブルーチンを終了し、program snec のメインループに戻る。
end subroutine blstep