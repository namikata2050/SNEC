! subroutine input_parser: パラメータファイルの解析ルーチン
!
! このルーチンは、シミュレーションの挙動を制御する設定ファイル "parameters" を読み込み、
! blmod と parameters モジュール内の様々な変数（パラメータやフラグ）を設定します。
subroutine input_parser
!
! This routine parses the "parameters" file and sets various flags
! that are kept in the modules blmod and parameters.
!
! ===== モジュールの読み込み =====
  use blmod, only: gravity_switch, wipe_outdir, dtfac, cvisc ! blmod から変数を読み込む
  use parameters 							 ! parameters モジュール全体を読み込む (パラメータ変数を設定するため)
  implicit none

! ===== ローカル変数の宣言 =====
  character(len=128) :: cpstring 	! ファイルコピー用のコマンド文字列
  character(len=500) :: rmstring 	! ファイル削除用のコマンド文字列
  logical :: opt 					! パラメータがオプション(省略可能)かどうかを示すフラグ (.false.=必須, .true.=オプション)
  logical :: outdirthere 			! 出力ディレクトリが存在するかどうかを示すフラグ

!------------------------------------------------------------------------------

! opt フラグを .false. (必須パラメータ) に初期化
  opt = .false.

!****************************** LAUNCH (実行設定) ****************************************
! 出力ディレクトリ名 'outdir' をパラメータファイルから読み込む (必須)
  call get_string_parameter('outdir',outdir,opt)

!************************* STELLAR PROFILE (初期モデル設定) ************************************
! 星の初期モデルファイル名 'profile_name' を読み込む (必須)
  call get_string_parameter('profile_name',profile_name,opt)
! 化学組成プロファイル名 'comp_profile_name' を読み込む (必須)
  call get_string_parameter('comp_profile_name',composition_profile_name,opt)


!***************************** EXPLOSION (爆発設定) **************************************
! 爆発の種類 'initial_data' ("Piston_Explosion" or "Thermal_Bomb") を読み込む (必須)
  call get_string_parameter('initial_data',initial_data,opt)

! もしピストン爆発なら、関連パラメータを読み込む (必須)
  if(initial_data.eq."Piston_Explosion") then
  	  call get_double_parameter('piston_vel',piston_vel,opt) 		! ピストン速度
  	  call get_double_parameter('piston_tstart',piston_tstart,opt) ! ピストン開始時刻
  	  call get_double_parameter('piston_tend',piston_tend,opt) 	! ピストン終了時刻
  endif

! もし熱爆弾なら、関連パラメータを読み込む
  if(initial_data.eq."Thermal_Bomb") then
  	  call get_double_parameter('final_energy',final_energy,opt) 		! 総注入エネルギー (必須)
  	  call get_double_parameter('bomb_tstart',bomb_tstart,opt) 		! 注入開始時刻 (必須)
  	  call get_double_parameter('bomb_tend',bomb_tend,opt) 			! 注入終了時刻 (必須)
  	  call get_double_parameter('bomb_mass_spread',bomb_mass_spread,opt)! エネルギー注入質量範囲 (必須だが使われていない可能性あり)
  	  call get_integer_parameter('bomb_start_point',bomb_start_point,opt)! 注入開始グリッド (必須)
  	  ! bomb_mode はオプション(.true.)。見つからなければ -666 が返る。
  	  call get_integer_parameter('bomb_mode',bomb_mode,.true.)
  	  ! もし bomb_mode が見つからなかったら、デフォルト値 1 を設定する。
  	  if(bomb_mode.eq.-666) then
  	 	 bomb_mode = 1 ! default
  	  endif
  endif

!******************************** GRID (グリッド設定) ****************************************
! グリッドゾーン数 'imax' を読み込む (必須)
  call get_integer_parameter('imax',imax,opt)
! グリッド生成法 'gridding' を読み込む (必須)
  call get_string_parameter('gridding',gridding,opt)
! 質量切り捨てフラグ 'mass_excision' を読み込む (必須)
  call get_logical_parameter('mass_excision',mass_excision,opt)
! もし質量切り捨てを行うなら、切り捨てる質量 'mass_excised' を読み込む (必須)
  if(mass_excision) then
  	  call get_double_parameter('mass_excised',mass_excised,opt)
  end if

!****************************** EVOLUTION (物理計算設定) *************************************
! 輻射輸送フラグ 'radiation' を読み込む (必須)
  call get_logical_parameter('radiation',radiation,opt)
! 状態方程式キー 'eoskey' を読み込む (必須)
  call get_integer_parameter('eoskey',eoskey,opt)
! Ni崩壊加熱スイッチ 'Ni_switch' を読み込む (必須)
  call get_integer_parameter('Ni_switch',Ni_switch,opt)
! Ni質量 'Ni_mass' を読み込む (必須だが Ni_by_hand=1 の場合のみ使われる)
  call get_double_parameter('Ni_mass',Ni_mass,opt)
! Ni分布境界質量 'Ni_boundary_mass' を読み込む (必須だが Ni_by_hand=1 の場合のみ使われる)
  call get_double_parameter('Ni_boundary_mass',Ni_boundary_mass,opt)
! Ni加熱再計算周期 'Ni_period' を読み込む (必須)
  call get_double_parameter('Ni_period',Ni_period,opt)
! Ni分布を手動設定するかのフラグ 'Ni_by_hand' (オプション, デフォルト=1)
  call get_integer_parameter('Ni_by_hand',Ni_by_hand,.true.)
  if (Ni_by_hand.eq.-666) Ni_by_hand = 1 ! set to default value

! Saha電離で考慮する元素数 'saha_ncomps' を読み込む (必須)
  call get_integer_parameter('saha_ncomps',saha_ncomps,opt)
  
! Boxcarスムージングフラグ 'boxcar_smoothing' を読み込む (必須)
  call get_logical_parameter('boxcar_smoothing',boxcar_smoothing,opt)

! オパシティ下限値(エンベロープ) 'opacity_floor_envelope' を読み込む (必須)
  call get_double_parameter('opacity_floor_envelope',of_env,opt)
! オパシティ下限値(コア) 'opacity_floor_core' を読み込む (必須)
  call get_double_parameter('opacity_floor_core',of_core,opt)

!********************** WHEN TO DO THINGS (時間・出力設定) *************************************
! 最大ステップ数 'ntmax' を読み込む (必須)
  call get_integer_parameter('ntmax',ntmax,opt)
! 終了時刻 'tend' を読み込む (必須)
  call get_double_parameter('tend',tend,opt)
! 詳細データ出力時間間隔 'dtout' を読み込む (必須)
  call get_double_parameter('dtout',dtout,opt)
! スカラ量出力時間間隔 'dtout_scalar' を読み込む (必須)
  call get_double_parameter('dtout_scalar',dtout_scalar,opt)
! チェックポイント出力時間間隔 'dtout_check' を読み込む (必須)
  call get_double_parameter('dtout_check',dtout_check,opt)
! 詳細データ出力ステップ間隔 'ntout' を読み込む (必須)
  call get_integer_parameter('ntout',ntout,opt)
! スカラ量出力ステップ間隔 'ntout_scalar' を読み込む (必須)
  call get_integer_parameter('ntout_scalar',ntout_scalar,opt)
! チェックポイント出力ステップ間隔 'ntout_check' を読み込む (必須)
  call get_integer_parameter('ntout_check',ntout_check,opt)
! 画面情報出力ステップ間隔 'ntinfo' を読み込む (必須)
  call get_integer_parameter('ntinfo',ntinfo,opt)
! 最小時間刻み幅 'dtmin' を読み込む (必須)
  call get_double_parameter('dtmin',dtmin,opt)
! 最大時間刻み幅 'dtmax' を読み込む (必須)
  call get_double_parameter('dtmax',dtmax,opt)

  ! CFL safety factor (optional)
  call get_double_parameter('dtfac',dtfac,.true.)

  ! Artificial viscosity coefficient (optional)
  call get_double_parameter('cvisc',cvisc,.true.)
  
!********************************** TEST (テスト設定) **************************************
! Sedovテストモードフラグ 'sedov' を読み込む (必須)
  call get_logical_parameter('sedov',sedov,opt)
! もし Sedov テストモードなら、重力スイッチを 0 (OFF) にする
  if(sedov) then
  	  gravity_switch = 0
  else
  	  gravity_switch = 1
  end if

!******************************************************************************
! ===== 出力ディレクトリのチェックと準備 =====

! check if output directory exists (出力ディレクトリが存在するか確認)
! (コンパイラ依存の処理: Intel Compiler とそれ以外で inquire の書式が違う)
#if __INTEL_COMPILER
  inquire(directory=trim(adjustl(outdir)),exist=outdirthere)
#else
  inquire(file=trim(adjustl(outdir)),exist=outdirthere)
#endif
! もし出力ディレクトリが存在しなければ、エラーメッセージを表示して停止する
  if(.not.outdirthere) then
  	  write(6,*) "*** Output directory does not exist."
  	  write(6,*) "Please create the output directory: ", trim(adjustl(outdir))
  	  stop
  endif

! wipe output dir if requested: (もし要求されていれば出力ディレクトリを空にする)
! wipe_outdir は blmod で .true. に初期化されている
  if(wipe_outdir) then
  	  write(*,*) "Removing output directory contents: ", trim(outdir)
  	  ! OS依存のファイル削除コマンドを組み立てる (Windows の del コマンド)
  	  write(rmstring,*) "del /s /q ", trim(adjustl(outdir)), '\*'
  	  ! system コマンドで外部コマンドを実行する
  	  call system(rmstring)
  endif
  
! copy parameter file (パラメータファイルを結果確認用にコピーする)
  ! cpstring="cp parameters "//trim(adjustl(outdir)) ! Linux/macOS 用
  ! OS依存のファイルコピーコマンドを組み立てる (Windows の copy コマンド)
  cpstring="copy parameters "//trim(adjustl(outdir))
  call system(cpstring)


!######## subroutines to get parse the parameters file ########################
!
! ===== ヘルパーサブルーチン (パラメータファイル解析用) =====
!
contains
! --------------------------------------------------------------------------
! get_string_parameter: 文字列パラメータをファイルから読み込む基本ルーチン
! --------------------------------------------------------------------------
  	subroutine get_string_parameter(parname,par,opt)

  	 	implicit none
  	 	logical opt 			! 入力: .true. ならオプション(見つからなくてもエラーにしない)
  	 	character*(*) parname 	! 入力: 探すパラメータ名 (例: 'outdir')
  	 	character*(*) par 		! 出力: 見つかったパラメータの値 (文字列)
  	 	character*(200) line_string ! ファイルから読み込んだ1行
  	 	integer i,j,l,ll 		! 作業用変数
  	 	character*(200) temp_string ! 作業用変数

  	 	! パラメータファイル 'parameters' をユニット番号 27 で開く
  	 	open(unit=27,file='parameters',status='unknown')

  	 	! ----- ファイルを行ごとに読み込み、目的のパラメータを探すループ -----
  	 	10 continue
  	 	! 1行読み込む (ファイル終端なら 19 へ、エラーなら 10 へ)
  	 	read(27,'(a)',end=19,err=10) line_string
  	 	! separator is an equal sign '=', # is comment (区切り文字は'=', '#'以降はコメント)
  	 	i = index(line_string,'=') ! '=' の位置を探す
  	 	j = index(line_string,'#') ! '#' の位置を探す

  	 	! '=' がない行、または行頭が '#' の行はスキップ
  	 	if (i.eq.0.or.j.eq.1) goto 10
  	 	!   if the whole line is a comment or there is no
  	 	!   equal sign, then go on to the next line  	

  	 	! '#' が '=' より前にある場合 (例: # param = value) もスキップ
  	 	if(j.gt.0.and.j.lt.i) goto 10
  	 	!   if there is an equal sign, but it is in a comment
  	 	!   then go on to the next line

  	 	! is this the right parameter? If not, cycle
  	 	! (これが探しているパラメータか？)
  	 	! '=' より前の部分を抽出し、空白を除去
  	 	temp_string=trim(adjustl(line_string(1:i-1)))
  	 	l=len(parname)
  	 	! パラメータ名が一致しなければ次の行へ
  	 	if(parname.ne.temp_string(1:l)) goto 10

  	 	! ----- パラメータが見つかった場合の処理 -----
  	 	!  If there is a comment in the line, exclude it!
  	 	! (もし行の途中に '#' があれば、それ以降を除外する)
  	 	l = len(line_string)
  	 	if (j.gt.0) l = j - 1

  	 	! '=' より後ろの部分 (パラメータ値) を抽出
  	 	par = line_string(i+1:l)
  	 	! now remove potential crap! (不要な文字を除去)
  	 	do ll=1,len(par)
  	 	 	 if(par(ll:ll).eq.'\t') par(ll:ll) = ' ' ! タブをスペースに
  	 	 	 if(par(ll:ll).eq.'"') par(ll:ll) = ' ' ! ダブルクォートを除去
  	 	 	 if(par(ll:ll).eq."'") par(ll:ll) = ' ' ! シングルクォートを除去
  	 	enddo
  	 	! adjust left... (先頭の空白を除去)
  	 	par = adjustl(par)
  	 	! get rid of trailing blanks (末尾の空白を除去)
  	 	j = index(par," ")
  	 	if (j.gt.0) par = par(1:j-1) ! 最初の空白以降を切り捨てる


  	 	! now look for " or ' and remove them (再度クォートをチェックし、あればエラー)
  	 	j=index(par,'"')
  	 	if(j.ne.0) stop "No quotes in my strings, please!"

  	 	j=index(par,"'")
  	 	if(j.ne.0) stop "No quotes in my strings, please!"

  	 	! ファイルを閉じて正常終了
  	 	close(27)
  	 	return

  	 	! ----- パラメータが見つからなかった場合の処理 (ファイルの終端に到達) -----
  	 	19 continue
  	 	! もしパラメータが必須 (.not.opt) なら、エラーメッセージを表示して停止
  	 	if(.not.opt) then
  	 	 	 write(6,*) "Fatal problem in input parser:"
  	 	 	 write(6,*) "Parameter ",parname
  	 	 	 write(6,*) "could not be read!"
  	 	 	 write(6,*) 
  	 	 	 call flush(6) ! メッセージを確実に出力
  	 	 	 stop
  	 	! もしパラメータがオプション (.true.) なら、
  	 	else
  	 	 	 ! 特殊な文字列 "NOTTHERE" を返して終了
  	 	 	 par = "NOTTHERE"
  	 	 	 close(27)
  	 	endif

  	end subroutine get_string_parameter

! --------------------------------------------------------------------------
! get_double_parameter: 倍精度実数(real*8)パラメータを読み込む
! --------------------------------------------------------------------------
  	subroutine get_double_parameter(parname,par,opt)

  	 	implicit none
  	 	logical opt
  	 	character(*) parname
  	 	character*256 line_string ! 文字列として読み込むための中間変数
  	 	real*8 par 				! 出力: 読み込んだ倍精度実数値

  	 	! まず、文字列としてパラメータ値を読み込む
  	 	call get_string_parameter(parname,line_string,opt)

        if((opt).and.trim(adjustl(line_string)) &
                 .eq."NOTTHERE") then
            return
        endif

  	 	! 簡単なフォーマットチェック: '.' がなければエラー (整数と区別するため？)
  	 	if(index(line_string,'.').eq.0) then
  	 	 	 write(6,*) "Uh. Bad double parameter ",trim(parname)
  	 	 	 write(6,*) "Please check input file!"
  	 	 	 call flush(6)
  	 	 	 stop
  	 	endif

  	 	! 文字列を倍精度実数(指数形式 e20.15)として読み取り、par に格納
  	 	read(line_string,"(e20.15)") par

  	end subroutine get_double_parameter

! --------------------------------------------------------------------------
! get_integer_parameter: 整数(integer)パラメータを読み込む
! --------------------------------------------------------------------------
  	subroutine get_integer_parameter(parname,par,opt)

  	 	implicit none
  	 	logical opt
  	 	character(*) parname
  	 	character*256 line_string ! 文字列として読み込むための中間変数
  	 	integer par 				! 出力: 読み込んだ整数値

  	 	! まず、文字列としてパラメータ値を読み込む
  	 	call get_string_parameter(parname,line_string,opt)

  	 	! オプションパラメータで、かつ見つからなかった場合 ("NOTTHERE" が返ってきた場合)
  	 	if((opt).and.trim(adjustl(line_string)) &
  	 	 	 	 .eq."NOTTHERE") then
             par = -666 ! Default marker
             return     ! Return early, do not attempt to read
        endif
  	 	 	 
  	 	read(line_string,"(i10)") par

  	end subroutine get_integer_parameter

! --------------------------------------------------------------------------
! get_logical_parameter: 論理値(logical)パラメータを読み込む
! --------------------------------------------------------------------------
  	subroutine get_logical_parameter(parname,par,opt)

  	 	implicit none
  	 	logical opt
  	 	character*(*) parname
  	 	character*(50) value_string ! 文字列として読み込むための中間変数
  	 	integer temp_par 			! 一時的に整数として読み込む
  	 	logical par 				! 出力: 読み込んだ論理値 (.true. or .false.)


  	 	! まず、文字列としてパラメータ値を読み込む
  	 	call get_string_parameter(parname,value_string,opt)

  	 	! オプションパラメータの場合の処理
  	 	if(opt) then
  	 	 	 ! don't try to set the parameter if it is not
  	 	 	 ! in the input file
  	 	 	 ! (パラメータがファイルになくても設定しようとしない)
  	 	 	 ! もし "NOTTHERE" が返ってきたら
  	 	 	 if(value_string .eq. "NOTTHERE") then
  	 	 	 	 ! デフォルト値を使う旨のメッセージを表示して、値を変更せずに終了
  	 	 	 	 write(6,*) "*** Parameter ",trim(adjustl(parname)), &
  	 	 	 	 	 	 "not found in input file. Using default value."
  	 	 	 	 return
  	 	 	 endif
  	 	endif
  	 	! 文字列を一時的に整数として読み込む (0 or 1 を想定)
  	 	read(value_string,"(i10)") temp_par

  	 	! 整数値に基づいて論理値を設定
  	 	if(temp_par.ne.0) then ! 0 以外なら .true.
  	 	 	 par = .true.
  	 	else 				   ! 0 なら .false.
  	 	 	 par = .false.
  	 	endif

  	end subroutine get_logical_parameter

! input_parser サブルーチン本体の終了
end subroutine input_parser