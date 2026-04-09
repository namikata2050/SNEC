subroutine magnetar_heat


  use blmod, only: magnetar_heating, time, delta_mass, mass, nt, magnetar_total_luminosity  
  use parameters

  use physical_constants
 
  implicit none

!local:
  real*8 :: Lsd, tsd_sec
  real*8 :: coef_a, coef_b
  real*8, parameter :: ratio_mass = 10000.0d0 ! 100.0d0
  integer :: k
  real*8 :: exponent_array(kout)
  
!----------------------------
  
  !local rate of gamma-ray energy deposition at a given grid point
  
  tsd_sec = tsd*86400d0 ! Spin-down timescale (sec)
  
 Lsd = Lsd_ini*(1d0 + time/tsd_sec)**(-2)
  
  ! Injecting the heat into the cells from 1 to kout.

  magnetar_heating(:) = 0
  
  !  Uniformly inject heat.
!!  if (time .ge. 10.) then  ! turn on 10 sec after explosion.
     do k=1, imax
        if (k .le. kout) then
           magnetar_heating(k) = Lsd/sum(delta_mass(1:kout)) ! (mass(kout) - mass_excised*msun)
        end if
     end do
!!  end if
  

  ! Inject heat exponentially in mass.

  coef_a = log(ratio_mass)/(mass(kout) - mass(1))
  
  do k =1, kout 
     exponent_array(k) = delta_mass(k) * exp( - coef_a * mass(k) )
  end do
  
  coef_b = Lsd/sum(exponent_array(1:kout))
  
 ! do k= 1, imax
 !    magnetar_heating(k) = coef_b*exp(-coef_a*mass(k))
 ! end do

   magnetar_total_luminosity = dot_product(magnetar_heating, delta_mass)
  
  ! Write Lsd to the terminal
  if (mod(nt, 5000) .eq. 0) then
     write(*,*) "Lsd (erg/s) = ", Lsd
     write(*,*) "L_inject (erg/s) = ",  magnetar_total_luminosity 
  end if
  
  end subroutine magnetar_heat
