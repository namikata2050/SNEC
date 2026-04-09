subroutine artificial_viscosity

  use blmod, only: vel, rho, Q, cr, cr_p, delta_mass, dtime, Qterm, cvisc
  use parameters
  use physical_constants
  implicit none

  integer :: i
  ! real*8, parameter :: cvisc = 4.0d0 ! original : 2.0d0 10.0d0? timestep.f90と合わせること

!------------------------------------------------------------------------------
! Richtmeyer-Von Neuman artificial viscosity

  do i=1,imax-1          
      if(vel(i+1).lt.vel(i)) then
          Q(i) = cvisc*rho(i) &
           * (vel(i+1) - vel(i))**2
      else
          Q(i) = 0
      endif
  enddo

  !outer boundary
  Q(imax) = 0.0d0


  !the term with viscosity at the right hand side of the energy equation
  do i=1,imax-1
      Qterm(i) = - 4.0d0*pi*dtime* (0.5d0 * (cr(i) + cr_p(i)))**2 &
          * Q(i) *(vel(i+1) - vel(i))/delta_mass(i)
  enddo
  Qterm(imax) = 0.0d0


end subroutine artificial_viscosity
