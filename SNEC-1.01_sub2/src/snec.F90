! program snec: SNEC（超新星爆発コード）のメインプログラム
! シミュレーション全体の流れ（初期化、メインループ、終了処理）を制御する司令塔です。
program snec

! ===== モジュールの読み込み =====
! 必要なモジュール（グローバル変数やサブルーチンが定義されたファイル）を読み込みます。
! 
! blmod: シミュレーションの「状態」を保持するモジュール。
!        時間(time)、時間刻み幅(dtime)、ステップ数(nt)、密度(rho)など、
!        計算の核となるグローバル変数を読み込みます。
  use blmod, only: dtime, dtime_p, time, nt, ntstart, tstart,   &
  	 	tdump, tdump_scalar, rho, tdump_check
! parameters: ユーザーが設定するシミュレーションのパラメータ
!             (例: tend, ntmax, dtoutなど) を保持するモジュール。
  use parameters
! outinfomod: 画面出力(stdout)に関する変数を読み込みます。
  use outinfomod, only: outinfo_count
! implicit none: 宣言されていない変数の使用を禁止する、Fortranの重要なおまじない。
!                タイプミスなどによるバグを防ぎます。
  implicit none

! ===== ローカル変数の宣言 =====
! このプログラム(snec)の中だけで使うローカル変数を宣言します。
! logicalは .true. (真) または .false. (偽) の値を持つ型です。
! これら3つは、データを出力するタイミングかどうかを判定する「旗（フラグ）」として使われます。
  logical :: OutputFlag = .false.
  logical :: OutputFlagScalar = .false.
  logical :: OutputFlagCheck = .false.

!------------------------------------------------------------------------------

! ===== タイトルバナーの表示 =====
! プログラム実行時に、画面にタイトルを表示します。
  write(*,*)
  
  write(*,*) "***********************************"
  write(*,*) "* Supernova Explosion Code (SNEC) *"
  write(*,*) "***********************************"

  write(*,*)
! *****************************************************
! INITIALIZATION (初期化)
! *****************************************************        

! 'input_parser' サブルーチンを呼び出します。
! ユーザーが作成したパラメータファイル（param.inなど）を読み込み、
! tend (終了時刻) や dtout (出力間隔) などの変数を設定します。
  call input_parser

! 'problem' サブルーチンを呼び出します。
! シミュレーションの「初期条件」を設定します。
! (例: 爆発前の星の密度分布、温度分布などをグリッド上に設定する)
  call problem

! 'artificial_viscosity' (人工粘性) サブルーチンを呼び出します。
! 衝撃波を数値的に安定して扱うための初期設定を行います。
  call artificial_viscosity

! output before first timestep
! ===== 初期状態(t=0)の出力 =====
! 計算を1ステップも進める前の「初期条件」をデータファイルに出力します。
! (0, 1, 2 は異なる種類の出力ファイルに対応)
  call output_all(0)
  call output_all(1)
  call output_all(2)
 
! 'timestep' サブルーチンを呼び出します。
! 物理状態（CFL条件など）に基づいて、「最初(1回目)の時間刻み幅 dtime (Δt)」を計算します。
  call timestep

! ===== 時間・ステップ制御変数の初期化 =====
! tstart (開始時刻) や dtout (出力間隔) を使って、
! 「次に各データを出力すべき時刻」を初期設定します。
  tdump_check = tstart+dtout_check
  tdump_scalar = tstart+dtout_scalar
  tdump = tstart+dtout
! 現在時刻(time)と現在ステップ数(nt)を、開始時の値(tstart, ntstart)に設定します。
  time = tstart
  nt = ntstart
  

! *****************************************************
! MAIN LOOP (メインループ)
! *****************************************************
! シミュレーションの心臓部。
! 'exit' が呼ばれるまで、ここの処理が繰り返し実行されます。

  IntegrationLoop: do
  	  
! ----- 時間刻み幅(dt)の計算 -----
! 
! 現在の時間刻み幅 dtime を dtime_p (previous = 前の) に保存します。
  	  dtime_p = dtime
  	  ! determine dt
! 'timestep' サブルーチンを呼び出し、現在の物理状態に基づいて
! 「次」の計算に使う時間刻み幅 dtime (Δt) を決定します。
  	  call timestep
  	  

! ----- 進捗状況の画面出力 -----
! もし ntinfo (例: 100) ステップごとに画面出力する設定なら、
! 'outinfo' サブルーチンを呼び出し、現在の時刻やステップ数を画面に表示します。
! mod(nt, ntinfo) は nt を ntinfo で割った余りを計算します。
  	  if(ntinfo.gt.0) then
  	 	 if(mod(nt,ntinfo).eq.0) then
  	 	 	  ! print useful info to stdout
  	 	 	  call outinfo
  	 	 endif
  	  endif
  	  
! ----- 最終ステップの調整 -----
! もし次のステップ (time + dtime) が終了時刻 (tend) を超えてしまう場合、
! dtime を (tend - time) に調整します。
! これにより、シミュレーションは正確に tend の時刻で停止します。
  	  if((time+dtime).gt.tend) dtime = tend-time

! ----- 物理計算の実行 (1ステップ進める) -----
  	  ! actual integration step
! 'blstep' サブルーチンを呼び出します。
! これがシミュレーションの「心臓部」であり、物理方程式（流体力学など）を解いて、
! シミュレーションの状態を t から t + dtime へと1ステップ進めます。
  	  call blstep
  	  
! ----- ステップ数のインクリメント -----
  	  ! increment timestep
! タイムステップ数を1増やします。
  	  nt = nt + 1
 
! ----- 出力タイミングの判定 (フラグ立て) -----
  	  ! various output related things
! 【ステップ数に基づく判定】
! もし ntout (例: 1000) ステップごとに出力する設定で、
! 現在の nt が ntout の倍数になったら、OutputFlag を .true. (出力する) にします。
  	  if (ntout.gt.0) then
  	 	 if ( mod(nt,ntout) .eq. 0) OutputFlag = .true.
  	  endif
  	
  	  if (ntout_scalar.gt.0) then
  	 	 if ( mod(nt,ntout_scalar) .eq. 0 ) OutputFlagScalar = .true.
  	  endif

  	  if (ntout_check.gt.0) then
  	 	 if ( mod(nt,ntout_check) .eq. 0 ) OutputFlagCheck = .true.
  	  endif
  	  
! 【シミュレーション時間に基づく判定】
! もし現在の time が、次に出力すべき時刻 tdump を超えていたら、
  	  if ( time.ge.tdump) then
! 	 1. 「次」の出力時刻を tdump + dtout (例: 10秒後 + 10秒 = 20秒後) に更新します。
  	 	 tdump=tdump+dtout
! 	 2. OutputFlag を .true. (出力する) にします。
  	 	 OutputFlag = .true.
  	  endif

  	  if ( time.ge.tdump_scalar) then
  	 	 tdump_scalar=tdump_scalar+dtout_scalar
  	 	 OutputFlagScalar = .true.
  	  endif

  	  if ( time.ge.tdump_check) then
  	 	 tdump_check=tdump_check+dtout_check
  	 	 OutputFlagCheck = .true.
  	  endif
  	  
! ----- 時間のインクリメント -----
  	  ! increment time
! 'blstep' での計算が終わったので、シミュレーションの現在時刻を dtime だけ進めます。
  	  time = time+dtime
  	  
! ----- 実際の出力処理 (フラグ実行) -----
! ループの最後で、立てられたフラグをチェックします。
! もし OutputFlag が .true. なら、'output_all' を呼び出してデータを書き出します。
! 書き出しが終わったら、フラグを .false. に戻し、次のタイミングに備えます。
  	  if (OutputFlag) then
  	 	 call output_all(0)
  	 	 call output_all(1)
  	 	 OutputFlag = .false.
  	  endif
  	  
  	  if (OutputFlagScalar) then
  	 	 call output_all(2)
  	 	 OutputFlagScalar = .false.
  	  endif
  	  
  	  if (OutputFlagCheck) then
  	 	 OutputFlagCheck = .false.
  	  endif

! ----- 終了条件の判定 -----
! 
! 【1. 終了時刻に到達】
! もしシミュレーションが指定された終了時刻 tend にピッタリ到達したら、
  	  if (time.eq.tend) then
  	 	 write(*,*) "Done! :-) tend reached"
! 	 	 最後の状態をファイルに出力し、
  	 	 call output_all(0)
  	 	 call output_all(1)
  	 	 call output_all(2)
! 	 	 'exit' コマンドで 'IntegrationLoop' から脱出します。
  	 	 exit
  	  endif

! 【2. 最大ステップ数に到達】
! もし終了時刻に達する前に、指定された最大ステップ数 ntmax に達してしまったら、
  	  if (nt.ge.ntmax) then
  	 	 write(*,*) "Done! :-) ntmax reached"
! 	 	 同様に、最後の状態を出力し、
  	 	 call output_all(0)
  	 	 call output_all(1)
  	 	 call output_all(2)
! 	 	 'exit' コマンドで 'IntegrationLoop' から脱出します。(安全停止)
  	 	 exit
  	  endif
  	
! メインループの終端。'exit' されなければ、'IntegrationLoop: do' に戻ります。
  enddo IntegrationLoop

! メインプログラムの終了。
end program snec