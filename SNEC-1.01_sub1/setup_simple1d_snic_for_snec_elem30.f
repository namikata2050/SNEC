      program setup_simple1d_snic_for_snec
      implicit none
      
      integer nn
      parameter(nn=200) 
      real*8 r0  
      parameter(r0=2.2d13) ! original : 1.0d14 5.0d13
      real*8 vmax
      parameter(vmax=3.0d9)
      integer t0 !in sec
      parameter(t0=r0/vmax)
      real*8 zmetal
      parameter(zmetal=1.0)

      real*8 pow
      parameter(pow=10.0) !outer ejecta , inner=0.0
      real*8 ek
      parameter(ek=10.0*0.9) !original trial 1.6, default 1.0*0.9
      real*8 mece,m56ni,mime,mo,mej
      parameter(mece=0.0)
      parameter(m56ni=2.0*0.9) !default 0.1*0.9
      parameter(mime=4.0*0.9) !default 0.2*0.9
      parameter(mo=14.0*0.9) !default 0.7*0.9
      
      real*8 mix_ece,mix_56ni,mix_ime !mixed region in msun
      parameter(mix_ece=0.0) !0.1 for w7, 0 otherwise
      parameter(mix_56ni=0.05) !0.1 as a standard
      parameter(mix_ime=0.05) !0.05 as a standard, 0 if O layer thin

c fixed parames:
      integer nelm
      parameter(nelm=30)
c change the snec to allow 30 stable isotopes

c model constants:
      real*8 encmece0,encmece,encmece1
      real*8 encm56ni0,encm56ni,encm56ni1
      real*8 encmime0,encmime,encmime1
      real*8 encmo

      real*8 den(nn),vel(nn),rad(nn),encm(nn),dv
      real*8 mf(nn,nelm),mfni(nn)

      real*8 mf0_ece(nelm),ni0_ece,enuc0_ece
      real*8 mf0_56ni(nelm),ni0_56ni,enuc0_56ni
      real*8 mf0_ime(nelm),ni0_ime,enuc0_ime
      real*8 mf0_o(nelm),ni0_o,enuc0_o

c internal:
      integer i,j,k,l,m,n
      real*8 dum,dum2
      real*8 ve,enuc,rho0
      real*8 mdum,edum

c constants:
      real*8 msun,rsun,gg,pi,foe
      parameter(msun=1.989d33)
      parameter(rsun=6.96d10)
      parameter(gg=6.67d-8)
      parameter(pi=3.141592)
      parameter(foe=1.0d51)

      open(unit=31,status="unknown",file="inmodel_simple1d.dat")
      open(unit=32,status="unknown",file="log_inmodel_simple1d.dat")
c composition set:
c ece:
      do m=1,nelm
         mf0_ece(m)=0.0d0
      end do
      ni0_ece=0.1
      mf0_ece(26)=0.45
      mf0_ece(28)=0.45
      enuc0_ece=1.76
c 56ni:
      do m=1,nelm
         mf0_56ni(m)=0.0d0
      end do
      mf0_56ni(26)=0.05*zmetal
      mf0_56ni(28)=0.05*zmetal
      ni0_56ni=1.0-mf0_56ni(26)-mf0_56ni(28)
      enuc0_56ni=1.57
c ime:
      do m=1,nelm
         mf0_ime(m)=0.0d0
      end do
      ni0_ime=0.0d0
      mf0_ime(16)=0.29
      mf0_ime(18)=0.04
      mf0_ime(20)=0.03
      mf0_ime(22)=2.0d-6
      mf0_ime(24)=1.0d-5
      mf0_ime(26)=0.06
      mf0_ime(28)=5.0d-5
      dum=ni0_ime
      do m=1,nelm
         dum=dum+mf0_ime(m)
      end do
      mf0_ime(14)=1.0-dum
      write(*,*) "X(Si) in IME: ", mf0_ime(14)
      enuc0_ime=1.25
c o:     
      do m=1,nelm
         mf0_o(m)=0.0d0
      end do
      ni0_o=0.0d0
      mf0_o(6)=0.01 
      mf0_o(12)=0.03
      mf0_o(14)=0.03
      mf0_o(16)=0.01
      mf0_o(18)=0.01
      dum=ni0_o
      do m=1,nelm
         dum=dum+mf0_o(m)
      end do
      mf0_o(8)=1.0-dum
      write(*,*) "X(O) in O: ", mf0_o(8)
      enuc0_o=0.56


c density distribution: input=ek,mej,vmax, output=
      enuc=enuc0_ece*mece + enuc0_56ni*m56ni + enuc0_ime*mime
     &     + enuc0_o*mo
      mej=mece+m56ni+mime+mo
      ve=sqrt(2.0*ek*foe/mej/msun) !initiaql guess
      rho0=mej*msun/(4.0*pi/3.0*ve**3.*t0**3.) !initial guess
      dv=vmax/real(nn)
      mdum=0.0d0
      edum=0.0d0
      do i=1,nn
         vel(i)=dv*real(i)
         rad(i)=vel(i)*t0
         if(vel(i).le.ve) then
            den(i)=rho0
         else
            den(i)=rho0*(vel(i)/ve)**(-1.0*pow)
         end if
         if(i.eq.1) then
c     dum2=dum2+4.0*pi/3.0*den(i)*rad(i)**3.0*vel(i)**2.*0.5
            mdum=4.0*pi/3.0*den(i)*rad(i)**3.0/msun
            edum=1.0/2.0*den(i)*vel(i)**2.
     &           *4.0*pi/3.0*rad(i)**3.0/foe
         else
c     dum2=dum2+4.0*pi/3.0*den(i)
c     &           *(rad(i)**3.0-rad(i-1)**3.0)*vel(i)**2.*0.5
            mdum=mdum+den(i)
     &           *4.0*pi/3.0*(rad(i)**3.0-rad(i-1)**3.0)/msun
            edum=edum+1.0/2.0*den(i)*(0.5*(vel(i-1)+vel(i)))**2.0
     &           *4.0*pi/3.0*(rad(i)**3.0-rad(i-1)**3.0)/foe
         end if
      end do
      write(*,*) Mej, mdum
      write(*,*) ek, edum

c rescaling: 
      ve=ve*(ek/edum)**0.5*(mej/mdum)**(-0.5) !updated guess
      rho0=rho0*(ek/edum)**(-1.5)*(mej/mdum)**2. !initial guess
      dv=vmax/real(nn)
      mdum=0.0d0
      edum=0.0d0
      do i=1,nn
         vel(i)=dv*real(i)
         rad(i)=vel(i)*t0
         if(vel(i).le.ve) then
            den(i)=rho0
         else
            den(i)=rho0*(vel(i)/ve)**(-1.0*pow)
         end if
         if(i.eq.1) then
c     dum2=dum2+4.0*pi/3.0*den(i)*rad(i)**3.0*vel(i)**2.*0.5
            mdum=4.0*pi/3.0*den(i)*rad(i)**3.0/msun
            encm(i)=mdum
            edum=1.0/2.0*den(i)*vel(i)**2.
     &           *4.0*pi/3.0*rad(i)**3.0/foe
         else
c     dum2=dum2+4.0*pi/3.0*den(i)
c     &           *(rad(i)**3.0-rad(i-1)**3.0)*vel(i)**2.*0.5
            mdum=mdum+den(i)
     &           *4.0*pi/3.0*(rad(i)**3.0-rad(i-1)**3.0)/msun
            encm(i)=mdum
            edum=edum+1.0/2.0*den(i)*(0.5*(vel(i-1)+vel(i)))**2.0
     &           *4.0*pi/3.0*(rad(i)**3.0-rad(i-1)**3.0)/foe
         end if
      end do

      write(*,*) "Model Summary: "
      write(*,*) "Mej (aimed, result): ", mej, mdum
      write(*,*) "Ek (aimed, result): ", Ek, edum
      write(*,*) "M (ECE): ", mece
      write(*,*) "M (56Ni): ", m56ni
      write(*,*) "M (IME): ", mime
      write(*,*) "M (O): ", mo
      write(*,*) "Enuc: ", enuc
      write(*,*) "rho0", rho0
      write(*,*) "Ve (km/s): ", ve/1.0d5

      write(32,*) "Model Summary: "
      write(32,*) "Mej (aimed, result): ", mej, mdum
      write(32,*) "Ek (aimed, result): ", Ek, edum
      write(32,*) "M (ECE): ", mece
      write(32,*) "M (56Ni): ", m56ni
      write(32,*) "M (IME): ", mime
      write(32,*) "M (O): ", mo
      write(32,*) "Enuc: ", enuc
      write(32,*) "rho0", rho0
      write(32,*) "Ve (km/s): ", ve/1.0d5


c mixing
      encmece=mece
      encmece0=encmece-mix_ece
      encmece1=encmece+mix_ece
      encm56ni=encmece+m56ni
      encm56ni0=encm56ni-mix_56ni
      encm56ni1=encm56ni+mix_56ni
      encmime=encm56ni+mime
      encmime0=encmime-mix_ime
      encmime1=encmime+mix_ime
      encmo=encmime+mo
      write(*,*) "mece0: ", encmece0
      write(*,*) "mece: ", encmece
      write(*,*) "mece1: ", encmece1
      write(*,*) "m56ni0: ", encm56ni0
      write(*,*) "m56ni: ", encm56ni
      write(*,*) "m56ni1: ", encm56ni1
      Write(*,*) "mime0: ", encmime0
      write(*,*) "mime: ", encmime
      write(*,*) "mime1: ", encmime1
      write(*,*) "mo:",encmo

      write(32,*) "mece0: ", encmece0
      write(32,*) "mece: ", encmece
      write(32,*) "mece1: ", encmece1
      write(32,*) "m56ni0: ", encm56ni0
      write(32,*) "m56ni: ", encm56ni
      write(32,*) "m56ni1: ", encm56ni1
      Write(32,*) "mime0: ", encmime0
      write(32,*) "mime: ", encmime
      write(32,*) "mime1: ", encmime1
      write(32,*) "mo:",encmo

      write(*,*) "OK to proceed???"
      read(*,*)
      do i=1,nn
         if(encm(i).lt.encmece0) then
            mfni(i)=ni0_ece
            do k=1,nelm
               mf(i,k)=mf0_ece(k)
            end do
         elseif(encm(i).lt.encmece1) then
            mfni(i)=(ni0_56ni-ni0_ece)/(encmece1-encmece0)
     &           *(encm(i)-encmece0)+ni0_ece
            do k=1,nelm
               mf(i,k)=(mf0_56ni(k)-mf0_ece(k))
     &              /(encmece1-encmece0)
     &              *(encm(i)-encmece0)+mf0_ece(k)
            end do
         elseif(encm(i).lt.encm56ni0) then
            mfni(i)=ni0_56ni
            do k=1,nelm
               mf(i,k)=mf0_56ni(k)
            end do
         elseif(encm(i).lt.encm56ni1) then
            mfni(i)=(ni0_ime-ni0_56ni)/(encm56ni1-encm56ni0)
     &           *(encm(i)-encm56ni0)+ni0_56ni
            do k=1,nelm
               mf(i,k)=(mf0_ime(k)-mf0_56ni(k))
     &              /(encm56ni1-encm56ni0)
     &              *(encm(i)-encm56ni0)+mf0_56ni(k)
            end do
         elseif(encm(i).lt.encmime0) then
            mfni(i)=ni0_ime
            do k=1,nelm
               mf(i,k)=mf0_ime(k)
            end do
         elseif(encm(i).lt.encmime1) then
            mfni(i)=(ni0_o-ni0_ime)/(encmime1-encmime0)
     &           *(encm(i)-encmime0)+ni0_ime
            do k=1,nelm
               mf(i,k)=(mf0_o(k)-mf0_ime(k))
     &              /(encmime1-encmime0)
     &              *(encm(i)-encmime0)+mf0_ime(k)
            end do
         elseif(mo.gt.0.0) then
            mfni(i)=ni0_o
            do k=1,nelm
               mf(i,k)=mf0_o(k)
            end do
         else
            mfni(i)=ni0_ime
            do k=1,nelm
               mf(i,k)=mf0_ime(k)
            end do
         end if
      end do
      
      dum=0.0d0
      do i=1,nn
         write(31,100) encm(i),rad(i),den(i),dum,dum,vel(i),
     &        mfni(i),(mf(i,k),k=1,nelm)
      end do
      close(31)
 100  format(40(1x,1pe12.5))

      stop
      end