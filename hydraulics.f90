! Fortran 95: one-node soil-tree RC circuit, MPa, L, hours.
module hydraulics
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)
  real(dp) :: c=0.5_dp, r_other=0.2_dp, r_xylem=0.8_dp, height=1.0_dp
  real(dp) :: leaf_area=2.0_dp, pressure=101.325_dp
  real(dp) :: gmin=0.002_dp, gmax=0.20_dp, topt=25.0_dp, twidth=15.0_dp
  real(dp) :: d0=1.5_dp, gs50=-1.5_dp, sg=0.25_dp, p50=-3.6_dp, sx=0.35_dp
  ! van Genuchten: alpha in m^-1, n dimensionless; theta in m3/m3.
  real(dp) :: theta_r=0.05_dp, theta_s=0.45_dp, alpha=3.6_dp, vn=1.56_dp
  real(dp) :: max_step=0.01_dp, initial_psi=-0.2_dp, initial_loss=0.0_dp
  real(dp) :: loss_floor=1.e-6_dp
  namelist /model/ c,r_other,r_xylem,height,leaf_area,pressure,gmin,gmax, &
    topt,twidth,d0,gs50,sg,p50,sx,theta_r,theta_s,alpha,vn, &
    max_step,initial_psi,initial_loss,loss_floor
contains
  subroutine check(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not.ok) then
      write(*,'(A)') 'ERROR: '//message
      stop 1
    end if
  end subroutine check

  real(dp) function logistic(x)
    real(dp), intent(in) :: x
    if (x>=0.0_dp) then
      logistic=exp(-x)/(1.0_dp+exp(-x))
    else
      logistic=1.0_dp/(1.0_dp+exp(x))
    end if
  end function logistic

  real(dp) function soil_potential(theta)
    real(dp), intent(in) :: theta
    real(dp) :: se,m
    call check(theta>theta_r .and. theta<=theta_s, &
      'soil moisture must be > theta_r and <= theta_s')
    se=(theta-theta_r)/(theta_s-theta_r)
    m=1.0_dp-1.0_dp/vn
    soil_potential=-0.00981_dp*((se**(-1.0_dp/m)-1.0_dp)**(1.0_dp/vn))/alpha
  end function soil_potential

  real(dp) function vulnerability(psi)
    real(dp), intent(in) :: psi
    vulnerability=logistic((psi-p50)/sx)
  end function vulnerability

  subroutine fluxes(psi,loss,t,vpd,day,soil,g,e,q,r)
    real(dp), intent(in) :: psi,loss,t,vpd,day,soil
    real(dp), intent(out) :: g,e,q,r
    g=gmin+(gmax-gmin)*day*exp(-((t-topt)/twidth)**2)* &
      logistic((gs50-psi)/sg)/(1.0_dp+vpd/d0)
    e=64.8_dp*leaf_area*g*vpd/pressure
    r=r_other+r_xylem/max(1.0_dp-loss,loss_floor)
    q=(soil-0.00981_dp*height-psi)/r
  end subroutine fluxes

  subroutine step(psi,loss,dt,t,vpd,day,soil,e,q)
    ! Backward Euler water balance. Damage uses endpoint potential and
    ! running maximum; only one scalar nonlinear solve is needed.
    real(dp), intent(inout) :: psi,loss
    real(dp), intent(in) :: dt,t,vpd,day,soil
    real(dp), intent(out) :: e,q
    real(dp) :: old,lo,hi,mid,f,g,r,lnew,emax,source
    integer :: k
    old=psi
    source=soil-0.00981_dp*height
    emax=64.8_dp*leaf_area*gmax*vpd/pressure
    lo=min(old,source)-dt*emax/c-1.0_dp
    hi=max(old,source)+1.0_dp
    ! At the lower bracket storage loss exceeds the maximum possible
    ! transpiration; at the upper bracket inflow is negative.
    ! The residual need not be globally monotone with embolism, but the
    ! continuous sign bracket locates a root. Check timestep convergence.
    do k=1,100
      mid=0.5_dp*(lo+hi)
      lnew=max(loss,vulnerability(mid))
      call fluxes(mid,lnew,t,vpd,day,soil,g,e,q,r)
      f=c*(mid-old)/dt-q+e
      if (f>0.0_dp) then
        hi=mid
      else
        lo=mid
      end if
      if (hi-lo<1.e-11_dp*max(1.0_dp,abs(mid))) exit
    end do
    psi=0.5_dp*(lo+hi)
    loss=max(loss,vulnerability(psi))
    call fluxes(psi,loss,t,vpd,day,soil,g,e,q,r)
  end subroutine step
end module hydraulics

