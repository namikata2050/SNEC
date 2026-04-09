subroutine timestep
  ! this routine sets the timestep

  use blmod, only: cs2, r, vel, dtime, dtime_p, nt, delta_time
  use parameters
  use physical_constants
  implicit none

  integer :: i

  real*8 :: sound
  real*8 :: dttrans
!  real*8, parameter :: dtfac = 0.95d0
  real*8, parameter :: dtfac = 0.4d0
! KM ... check artifical viscocity as well, if you want to change here.
  real*8, parameter :: cvisc = 2.0d0

  real*8 :: delv
  real*8 :: dtc,dtq
!------------------------------------------------------------------------------

  dttrans = dtmax

  if(nt.le.1) then
    dttrans = 1.0d-8
    dtime = 1.0d-8
  endif

  !copy old timestep
  dtime_p = dtime


  dtime = dtmax
  do i=1,imax-1
!KM
    sound = sqrt(cs2(i))
    delv=abs(vel(i+1)-vel(i))
    delv=max(delv,1.0) !avoid zero 
    dtc=(r(i+1)-r(i))/sound
    dtq=(r(i+1)-r(i))/4.0/cvisc/delv
    delta_time(i)=min(dtc,dtq)  
    dtime=min(dtime,dtc,dtq)
!    delta_time(i) = (r(i+1) - r(i)) / &
!         max(abs(vel(i)+sound),abs(vel(i)-sound))
!    dtime = min(dtime, (r(i+1) - r(i)) / &
!         max(abs(vel(i)+sound),abs(vel(i)-sound)))
!
!KM 2018.03.02
!    if(delta_time(i).le.1.0d-2) then
!      write(*,*) i,delta_time(i),vel(i),sound,r(i+1)-r(i)
!    endif
 end do

  dtime = dtfac * dtime

  dtime = min(1.025d0*dtime_p,dtime)

  if(dttrans.lt.0.75d0*dtime) then
    dtime = dttrans
  endif

  dtime = min(max(dtime,dtmin),dtmax)


end subroutine timestep
