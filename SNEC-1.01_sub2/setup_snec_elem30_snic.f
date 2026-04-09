      program setup_snec_elem30
      implicit none

      real*8 msun,rsun,day,yr,pi,arad,kb, mwd
      parameter(msun=1.989d33)
      parameter(rsun=6.96d10)
      parameter(day=24.*60.*60.) !day in sec
      parameter(yr=365.*day) !year in sec
      parameter(pi=3.141592)
      parameter(arad=7.56d-15)
      parameter(kb=1.38d-16)
      parameter(mwd=1.3776068d0) !in Msun

      real*8 rhocsm1,rcsm0,rcsmin,powcsm1
      real*8 rhocsm2,rhocsm3,rhocsm4,powcsm2,powcsm3,powcsm4

      real*8 rt12,rt23,rt34

c 追加した
      real*8 :: r_gap_end
      integer :: n_gap, n_wave
      real*8 :: dr_gap, dr_wave


c final structure
c      parameter(rhocsm1=1.5d-15*30.0*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
c      parameter(powcsm1=1.0)
c      parameter(rhocsm2=5.0d21*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
c      parameter(powcsm2=3.3)
c      parameter(rhocsm3=0.85*25.*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
c      parameter(powcsm3=2.0)
c      parameter(rhocsm4=2.7d-10*25.0*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
c      parameter(powcsm4=1.4)
c final structure - 2020.12.28
c      parameter(rhocsm1=1.4*1.5d-15*30.0*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
c      parameter(powcsm1=1.0)
c      parameter(rhocsm2=1.4*5.0d21*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
c      parameter(powcsm2=3.3)
c      parameter(rhocsm3=1.6*0.85*25.*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
c      parameter(powcsm3=2.0)
c      parameter(rhocsm4=1.4*2.7d-10*25.0*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
c      parameter(powcsm4=1.4)
c tmp
      parameter(rhocsm1=1.4*1.5d-15*30.0*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
      parameter(powcsm1=1.0)
      parameter(rhocsm2=1.4*5.5d24*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
      parameter(powcsm2=3.5)
      parameter(rhocsm3=1.4*0.85*25.*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
      parameter(powcsm3=2.0)
      parameter(rhocsm4=1.4*2.7d-10*25.0*5.0d11) !in Msun ! rho_csm=rhocsm1*r**-powcsm
      parameter(powcsm4=1.4)

!1.2 for 1.0Msun
!0.006 for 0.005Msun
!0.012 for 0.01Msun
!0.06 for 0.05Msun
!0.12 for 0.1Msun
!0.35 for 0.3Msun
!0.68 for 0.6Msun 
      parameter(rcsm0=3.0d16) !in cm
      real*8 rhocsm0,rfac

      integer nn0,nn,nn_in,nncsm,nn0_nuc
      parameter(nn=200) !200
      parameter(nncsm=200)      !0 for no CSM, 100 for def, all meash = nn+nncsm, original : 200
! depending on nncsm, chenage the rfac in the main program

! ===== ↓↓ ここから追加 ↓↓ =====
! CSMの基本設定
      real*8 :: csm_rho_base, csm_power
      parameter(csm_rho_base=1.0d-9) ! CSM内縁(イジェクタ表面)での密度 [g/cm^3] original:2.0d-16は低すぎ
      ! r0=1.0d14の場合で連続にするなら2.0d-13
      ! 3.0d-12に上げる時、r0=5.0d13の場合で連続
      ! Mdot=1(6.3d25 g/s)の場合、v_wind=1.0d7 cm/sでr0=2.5d13 cm
      parameter(csm_power=2.0)        ! 密度のべき指数 (ρ ∝ r^-2 の風モデル)

c この場合光学的厚さは15で厚め。

! sin波のパラメータ
      real*8 :: csm_A, csm_B
      parameter(csm_A=0.0)            ! ゆらぎの振幅 (0.2 = ±20%)
      parameter(csm_B=1.0d-14)        ! ゆらぎの波数 (2π/B が周期 [cm] になる) default:1.0d-15
! ===== ↑↑ ここまで追加 ↑↑ =====

      integer nelm
      parameter(nelm=30) !30 for snec 

      real*8 zsol
      real*8 metallicity !0.02 for solar
      parameter(metallicity=0.02)

c parames:
c      real*8 vsn
c      parameter(vsn=2.5d9)
      integer nhyd_sn, nhyd_csm
      parameter(nhyd_sn=nn)
      parameter(nhyd_csm=nncsm)
      real*8 tempwd0,temp0
      parameter(tempwd0=1.0d2,temp0=1.0d2) ! original : temp01.0d2

c      real*8 rmin,rmax
c      parameter(rmax=8.4e15)


      integer i,j,k,m,l

      real*8 radwd(nn),encmwd(nn),rhowd(nn),tempwd(nn),velwd(nn)
      real*8 mfwd(nn,nelm),encmwd_nuc(nn)
      real*8 mf56niwd(nn)

      real*8 zz(nelm),aa(nelm),mf0(nelm)
      real*8 zz56ni,aa56ni

      real*8 mass(nn+nncsm),rad(nn+nncsm)
     &     ,temp(nn+nncsm)
     &     ,rho(nn+nncsm),vel(nn+nncsm)
     &     ,ye(nn+nncsm),omega(nn+nncsm)
     &     ,mf(nn+nncsm,nelm),mf56ni(nn+nncsm)
      real*8 check(nn+nncsm)

      real*8 masselm(nelm),masselm_wd(nelm)
      real*8 mass56ni,mass56niwd
      real*8 ene,mass_wd

      real*8 msum,dum,esum,eksum,etsum



c composition:

      do k=1,nelm
         zz(k)=real(k)
      end do   
      do k=2,nelm
         aa(k)=2.0*real(k)
      end do
      aa(1)=1.0d0 !H
      aa(28)=58.0 !stable ni
 
      zz56ni=28.0
      aa56ni=56.0

 

!solar
c solar abundance:
      mf0(1)=0.0 !H
      mf0(2)=0.0 !He
      mf0(3)=1.05d-8 !Li
      mf0(4)=1.74d-10 !Be
      mf0(5)=6.05d-9 !B
      mf0(6)=3.23d-3 !C
      mf0(7)=1.17d-3 !N
      mf0(8)=1.01d-2 !O
      mf0(9)=4.27d-7 !F
      mf0(10)=1.85d-3 !Ne
      mf0(11)=3.52d-5 !Na
      mf0(12)=6.95d-4 !Mg
      mf0(13)=6.10d-5 !Al
      mf0(14)=7.49d-4 !Si
      mf0(15)=8.59d-6 !P
      mf0(16)=4.40d-4 !S
      mf0(17)=3.57d-6 !Cl
      mf0(18)=9.78d-5 !Ar
      mf0(19)=3.93d-6 !K
      mf0(20)=6.53d-5 !Ca
      mf0(21)=4.10d-8 !Sc
      mf0(22)=3.07d-6 !Ti
      mf0(23)=3.98d-7 !V
      mf0(24)=1.87d-5 !Cr
      mf0(25)=1.40d-5 !Mn
      mf0(26)=1.34d-3 !Fe  
      mf0(27)=3.54d-6 !Co
      mf0(28)=7.73d-5 !Ni
      mf0(29)=8.85d-7 !Cu
      mf0(30)=2.20d-6 !Zn

      dum=0.0d0
      do j=1,nelm
         dum=dum+mf0(j)
      end do
      mf0(1)=1.0-dum !scaled by Hydrogen
      zsol=0.0
      do j=11,nelm !beyond Nr
         zsol=zsol+mf0(j)
      end do
      open(unit=21,status="old",file=
     &     "inmodel_simple1d.dat")

c begin:

      do i=1,nn
         read(21,*) encmwd(i),radwd(i),rhowd(i),dum
     &        ,dum,velwd(i),mf56niwd(i)
     &        ,(mfwd(i,k),k=1,nelm)
         tempwd(i)=tempwd0
      end do
      close(21)


c setup:
      do i=1,nn
         rad(i)=radwd(i)
         vel(i)=velwd(i)
         rho(i)=rhowd(i)
         temp(i)=tempwd(i)
         mf56ni(i)=mf56niwd(i)
         do k=1,nelm
            mf(i,k)=mfwd(i,k)
         end do
      end do
      mass(1)=4.0*pi/3.0*rad(1)**3.*rho(1)
      ye(1)=0.5
      omega(1)=0.0d0
      do i=2,nn
         mass(i)=mass(i-1)
     &        +4.0*pi/3.0*(rad(i)**3.-rad(i-1)**3.)*rho(i)
         ye(i)=0.5
         omega(i)=0.0d0
      end do
c normalized
      do i=1,nn
         dum=mf56ni(i)
         do j=1,nelm
            dum=dum+mf(i,j)
         end do
         mf56ni(i)=mf56ni(i)/dum
         do k=1,nelm
            mf(i,k)=mf(i,k)/dum
         end do
      end do




      rcsmin=rad(nn)
      if(rcsm0.le.rcsmin) then
         write(*,*) "csm wrong"
         stop
      end if



c      rhocsm0=masscsm0*msun/4.0/pi/(rcsmin**3.)/log(rcsm0/rcsmin)
c      write(*,*) "CSM demnsity scale", rhocsm0
      rfac=(rcsm0/rcsmin)**(1.0/real(nncsm))
c      for no CSM
c      rfac=1.0
      dum=log10(rhocsm2/rhocsm1)*1.0/(powcsm2-powcsm1)
      rt12=10.0**dum
      dum=log10(rhocsm3/rhocsm2)*1.0/(powcsm3-powcsm2)
      rt23=10.0**dum
       dum=log10(rhocsm4/rhocsm3)*1.0/(powcsm4-powcsm3)
      rt34=10.0**dum
      

      do i=nn+1,nn+nncsm
         rad(i)=rad(i-1)*rfac
         vel(i)=0.0d0
c         rho(i)=rhocsm0*(rad(i)/rcsmin)**(-3.)

! ========== 修正箇所ここから ==========
         
!         if(rad(i).le.rt12) then
!            rho(i)=rhocsm1*rad(i)**(-1.0*powcsm1)
!         elseif(rad(i).le.rt23) then
!            rho(i)=rhocsm2*rad(i)**(-1.0*powcsm2)
!         elseif(rad(i).le.rt34) then
!            rho(i)=rhocsm3*rad(i)**(-1.0*powcsm3)
!         else
!            rho(i)=rhocsm4*rad(i)**(-1.0*powcsm4)
!         end if

         ! 新しいロジック (単調減少 + sin波)
         ! dum 変数を一時的に使って基本密度を計算
         ! (イジェクタ表面 rad(nn) で密度が csm_rho_base になるようにつなぐ)
         dum = csm_rho_base * (rad(i)/rad(nn))**(-1.0*csm_power)
         
         ! sin波のゆらぎを乗算で加える
         rho(i) = dum * (1.0d0 + csm_A * sin(csm_B * rad(i)))

         ! ========== 修正箇所ここまで ==========
        
         temp(i)=temp0

         mf56ni(i)=0.0d0  ! ★これを追加：CSM領域のNickelを0で初期化する

         ! 0密度の部分を除く
         do j=1,nelm
            mf(i,j)=0.0d0
         end do
c         if(rad(i).gt.2.0*rcsmin) then
         mf(i,6)=0.5-metallicity*0.5
         mf(i,8)=0.5-metallicity*0.5
         mf(i,10)=metallicity-zsol
         do j=11,nelm
            mf(i,j)=mf(i,j)+mf0(j)
         end do

         dum=0.0d0
         do j=1,nelm
            dum=dum+mf(i,j)
         end do
         do k=1,nelm
            mf(i,k)=mf(i,k)/dum
         end do
      end do
      do i=nn+1,nn+nncsm
         mass(i)=mass(i-1)
     &        +4.0*pi/3.0*(rad(i)**3.-rad(i-1)**3.)*rho(i)
         ye(i)=0.5
         omega(i)=0.0d0
      end do

      open(unit=21,status="unknown",
     &file="profiles/simple1d_plus_csm.short")
      write(21,*) nn+nncsm
      do i=1,nn+nncsm
         write(21,100) i,mass(i),rad(i),temp(i),rho(i),vel(i)
     &        ,ye(i),omega(i)
      end do
 100  format(i6,7(4x,1pe16.9))

      open(unit=22,status="unknown",
     & file="profiles/simple1d_plus_csm.iso.dat")
      write(22,*) nn+nncsm,nelm+1
      write(22,150) (aa(j),j=1,nelm), aa56ni
      write(22,150) (zz(j),j=1,nelm), zz56ni
      do i=1,nn+nncsm
         write(22,200) mass(i),rad(i),(mf(i,j),j=1,nelm),mf56ni(i)
      end do
 150  format(31(2x,f5.1))
 200  format(33(3x,1pe12.6))

      open(unit=23,status="unknown",
     & file="GridPattern.dat")
c      write(23,*) nn
      write(23,*) "0"
c      do i=3,ncut+nncsm-1
c         write(23,250) mass(i)/mass(ncut+nncsm)
c      end do
      do i=1,nhyd_sn
         write(23,250) real(i)/real(nhyd_sn)
     &        *mass(nn)/mass(nn+nncsm)
      end do
      dum=mass(nn)/mass(nn+nncsm)
      do i=1,nhyd_csm
         write(23,250) real(i)/real(nhyd_sn+1)*(1.0-dum)+dum
      end do
      write(23,*) "1"
 250  format(1pe18.10)
      close(23)

      open(unit=23,status="old",file="GridPattern.dat")
      open(unit=24,status="unknown",file="checkgrid.dat")
c      do i=2,ncut+nncsm

      mass56ni=mass(1)*mf56ni(1)
      do k=1,nelm
         masselm(k)=mass(1)*mf(1,k)
      end do
      do i=1,nhyd_sn+nhyd_csm
         read(23,*) check(i)
         if(i.gt.2) then
            write(24,*) check(i)-check(i-1)
            if(check(i)-check(i-1).le.0.0) then
               write(*,*) "below 0 at ", i
            end if
         end if
         if(i.gt.1) then
            mass56ni=mass56ni+(mass(i)-mass(i-1))*mf56ni(i)
            do k=1,nelm
               masselm(k)=masselm(k)+(mass(i)-mass(i-1))*mf(i,k)
            end do
         end if
      end do


      write(*,*) "Ejecta mass", mass(nn)/msun
      write(*,*) "CSM mass1: ", (mass(nn+nncsm)-mass(nn))/msun

      mass56ni=mass56ni/msun
      do k=1,nelm
         masselm(k)=masselm(k)/msun
      end do
      write(*,*) "M(56Ni): ", mass56ni
      do k=1,nelm
         write(*,*) "mass elm", k, masselm(k)
      end do

      msum=0.0d0
      eksum=0.0d0
      etsum=0.0d0
      msum=4.0*pi/3.0*(rad(1)**3.)*rho(1)
      eksum=0.5*4.0*pi/3.0*(rad(1)**3.)*rho(1)*vel(1)**2.
      etsum=4.0*pi/3.0*(rad(1)**3.)*7.56d-15*(temp(1)**4.)
      do i=2,nhyd_sn     
         msum=msum
     &        +4.0*pi/3.0*(rad(i)**3.-rad(i-1)**3.)*rho(i)
         eksum=eksum+0.5*4.0*pi/3.0*(rad(i)**3.-rad(i-1)**3.)
     &        *rho(i)*vel(i)**2.
         etsum=etsum+4.0*pi/3.0*(rad(i)**3.-rad(i-1)**3.)
     &        *7.56d-15*(temp(i)**4.)
      end do
      write(*,*) "Ejecta mass", msum/msun
      write(*,*) "Ek",eksum/1.0d51
      write(*,*) "Eth", etsum/1.0d51
      write(*,*) "Ek + ETh", (eksum+etsum)/1.0d51

      stop
      end