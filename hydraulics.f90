! Version 2.0 -- F95. Hydraulic units MPa, L, hours; biological rates per DAY.
module hydraulics
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)
  real(dp) :: c=0.5_dp, r_other=0.2_dp, r_xylem=0.8_dp, height=1.0_dp
  real(dp) :: leaf_area=2.0_dp, pressure=101.325_dp
  real(dp) :: gmin=0.002_dp, gmax=0.20_dp, topt=25.0_dp, twidth=15.0_dp
  real(dp) :: d0=1.5_dp, gs50=-1.5_dp, sg=0.25_dp, p50=-3.6_dp, sx=0.35_dp
  real(dp) :: theta_r=0.05_dp, theta_s=0.45_dp, alpha=3.6_dp, vn=1.56_dp
  real(dp) :: max_step=0.01_dp, initial_psi=-0.2_dp, initial_loss=0.0_dp
  real(dp) :: loss_floor=1.e-6_dp
  ! Defaults zero preserve version 1 physics. Explicit file values may override.
  real(dp) :: growth_rate_per_day=0._dp, turnover_rate_per_day=0._dp
  real(dp) :: growth_soil_min=-0.5_dp, growth_psi_min=-1.0_dp
  real(dp) :: growth_temp_min=5._dp, growth_temp_max=40._dp
  real(dp) :: repair_rate_per_day=0._dp
  real(dp) :: repair_soil_min=-0.2_dp, repair_psi_min=-0.5_dp
  real(dp) :: repair_temp_min=5._dp, repair_temp_max=35._dp
  namelist /model/ c,r_other,r_xylem,height,leaf_area,pressure,gmin,gmax, &
    topt,twidth,d0,gs50,sg,p50,sx,theta_r,theta_s,alpha,vn, &
    max_step,initial_psi,initial_loss,loss_floor,growth_rate_per_day, &
    turnover_rate_per_day,growth_soil_min,growth_psi_min,growth_temp_min, &
    growth_temp_max,repair_rate_per_day,repair_soil_min,repair_psi_min, &
    repair_temp_min,repair_temp_max
contains
  subroutine check(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not.ok) then
      write(*,'(A)') 'ERROR: '//message
      stop 1
    end if
  end subroutine check

  subroutine validate_biology()
    call check(growth_rate_per_day>=0 .and. growth_rate_per_day<=1,'invalid growth rate')
    call check(turnover_rate_per_day>=0 .and. turnover_rate_per_day<=1,'invalid turnover rate')
    call check(repair_rate_per_day>=0 .and. repair_rate_per_day<=10,'invalid repair rate')
    call check(growth_temp_min<growth_temp_max,'invalid growth temperature range')
    call check(repair_temp_min<repair_temp_max,'invalid repair temperature range')
    call check(growth_soil_min<=0 .and. growth_psi_min<=0,'invalid growth water thresholds')
    call check(repair_soil_min<=0 .and. repair_psi_min<=0,'invalid repair water thresholds')
  end subroutine validate_biology

  real(dp) function logistic(x)
    real(dp), intent(in) :: x
    if (x>=0._dp) then
      logistic=exp(-x)/(1._dp+exp(-x))
    else
      logistic=1._dp/(1._dp+exp(x))
    end if
  end function logistic

  real(dp) function soil_potential(theta)
    real(dp), intent(in) :: theta
    real(dp) :: se,m
    call check(theta>theta_r .and. theta<=theta_s,'soil moisture outside retention bounds')
    se=(theta-theta_r)/(theta_s-theta_r)
    m=1._dp-1._dp/vn
    soil_potential=-0.00981_dp*((se**(-1._dp/m)-1._dp)**(1._dp/vn))/alpha
  end function soil_potential

  real(dp) function vulnerability(psi)
    real(dp), intent(in) :: psi
    vulnerability=logistic((psi-p50)/sx)
  end function vulnerability

  subroutine read_forcing(line,time,t,theta,vpd,day,growth,leaf)
    ! Strict numeric CSV: exactly 5 legacy fields or 7 extended fields.
    ! No embedded commas, empty fields, comments, slash terminators or NaNs.
    character(*), intent(in) :: line
    real(dp), intent(out) :: time,t,theta,vpd,day,growth,leaf
    real(dp) :: fields(7)
    integer :: j,start,n,ios,last
    character(2048) :: token
    n=0; start=1; last=len_trim(line)
    call check(last>0,'blank forcing row')
    do j=1,last+1
      if(j<=last) then
        if(line(j:j)/=',') cycle
      end if
      n=n+1
      call check(n<=7,'too many forcing columns')
      call check(j>start,'empty forcing field')
      token=adjustl(line(start:j-1))
      call check(len_trim(token)>0,'empty forcing field')
      call check(verify(trim(token),'0123456789+-.eEdD')==0,'invalid numeric forcing field')
      read(token,*,iostat=ios) fields(n)
      call check(ios==0,'cannot parse forcing field')
      start=j+1
    end do
    call check(n==5 .or. n==7,'forcing must have exactly 5 or 7 columns')
    time=fields(1); t=fields(2); theta=fields(3); vpd=fields(4); day=fields(5)
    growth=0._dp; leaf=1._dp
    if(n==7) then
      growth=fields(6); leaf=fields(7)
    end if
    call check(growth>=0 .and. growth<=1,'growth_factor must be in [0,1]')
    call check(leaf>=0 .and. leaf<=100,'leaf_area_factor must be in [0,100]')
  end subroutine read_forcing

  subroutine fluxes(psi,kpot,kfun,t,vpd,day,soil,leaf,g,e,q,r)
    real(dp), intent(in) :: psi,kpot,kfun,t,vpd,day,soil,leaf
    real(dp), intent(out) :: g,e,q,r
    g=gmin+(gmax-gmin)*day*exp(-((t-topt)/twidth)**2)* &
      logistic((gs50-psi)/sg)/(1._dp+vpd/d0)
    e=64.8_dp*leaf_area*leaf*g*vpd/pressure
    r=r_other+1._dp/max(kfun,loss_floor*kpot)
    q=(soil-0.00981_dp*height-psi)/r
  end subroutine fluxes

  subroutine capacity_update(psi,dt,t,soil,growth,kpot,kfun,knext,fstar,added,repaired,retired)
    ! Operator split: proportional retirement, new functional capacity, repair.
    ! Environmental gates use START-of-substep stem Psi and current forcing.
    real(dp), intent(in) :: psi,dt,t,soil,growth,kpot,kfun
    real(dp), intent(out) :: knext,fstar,added,repaired,retired
    real(dp) :: decay,dday,fp
    dday=dt/24._dp
    decay=exp(-turnover_rate_per_day*dday)
    retired=kpot*(1._dp-decay)
    knext=kpot*decay
    fp=kfun*decay
    added=0._dp
    if(soil>=growth_soil_min .and. psi>=growth_psi_min .and. &
       t>=growth_temp_min .and. t<=growth_temp_max) then
      ! Absolute addition relative to INITIAL potential conductance, not
      ! exponential compounding and not a biomass or carbon-budget model.
      added=(1._dp/r_xylem)*growth_rate_per_day*growth*dday
    end if
    knext=knext+added; fp=fp+added
    call check(knext>tiny(1._dp)*100._dp,'potential capacity depleted; revise horizon/turnover')
    repaired=0._dp
    if(soil>=repair_soil_min .and. psi>=repair_psi_min .and. &
       t>=repair_temp_min .and. t<=repair_temp_max) then
      repaired=max(0._dp,knext-fp)*(1._dp-exp(-repair_rate_per_day*dday))
    end if
    fstar=min(knext,max(0._dp,fp+repaired))
  end subroutine capacity_update

  subroutine step(psi,loss,kpot,kfun,dt,t,vpd,day,soil,growth,leaf,e,q,added,repaired,retired,damaged)
    real(dp), intent(inout) :: psi,loss,kpot,kfun
    real(dp), intent(in) :: dt,t,vpd,day,soil,growth,leaf
    real(dp), intent(out) :: e,q,added,repaired,retired,damaged
    real(dp) :: old,lo,hi,mid,f,g,r,knext,fstar,fnew,emax,source
    integer :: it
    call capacity_update(psi,dt,t,soil,growth,kpot,kfun,knext,fstar,added,repaired,retired)
    old=psi; source=soil-0.00981_dp*height
    emax=64.8_dp*leaf_area*leaf*gmax*vpd/pressure
    lo=min(old,source)-dt*emax/c-1._dp; hi=max(old,source)+1._dp
    do it=1,100
      mid=0.5_dp*(lo+hi)
      ! KEY P50 LINK: no more function than vulnerability allows at new Psi.
      fnew=min(fstar,knext*(1._dp-vulnerability(mid)))
      call fluxes(mid,knext,fnew,t,vpd,day,soil,leaf,g,e,q,r)
      f=c*(mid-old)/dt-q+e
      if(f>0._dp) then
        hi=mid
      else
        lo=mid
      end if
      if(hi-lo<1.e-11_dp*max(1._dp,abs(mid))) exit
    end do
    call check(it<=100,'hydraulic bisection did not converge')
    psi=0.5_dp*(lo+hi)
    kpot=knext
    kfun=min(fstar,kpot*(1._dp-vulnerability(psi)))
    damaged=max(0._dp,fstar-kfun)
    loss=min(1._dp,max(0._dp,1._dp-kfun/kpot))
    call fluxes(psi,kpot,kfun,t,vpd,day,soil,leaf,g,e,q,r)
  end subroutine step
end module hydraulics
