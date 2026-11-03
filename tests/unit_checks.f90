program unit_checks
  use hydraulics
  implicit none
  real(dp) :: kp,kf,kn,fs,ad,rp,rt,k0,psi,loss,e,q,damage,expected
  k0=1._dp/r_xylem
  growth_soil_min=-10;growth_psi_min=-10;growth_temp_min=0;growth_temp_max=40
  repair_soil_min=-10;repair_psi_min=-10;repair_temp_min=0;repair_temp_max=40
  kp=k0;kf=.25_dp*k0
  growth_rate_per_day=.1_dp;repair_rate_per_day=0;turnover_rate_per_day=0
  call capacity_update(-.2_dp,24._dp,20._dp,-.1_dp,1._dp,kp,kf,kn,fs,ad,rp,rt)
  call check(abs(kn-1.1_dp*k0)<1.e-12_dp,'growth potential capacity')
  call check(abs(fs-.35_dp*k0)<1.e-12_dp,'growth functional capacity')
  call check(1._dp-fs/kn<.75_dp,'growth must lower fractional loss when wet')
  growth_rate_per_day=0;repair_rate_per_day=.2_dp
  call capacity_update(-.2_dp,24._dp,20._dp,-.1_dp,0._dp,kp,kf,kn,fs,ad,rp,rt)
  expected=kf+(kp-kf)*(1._dp-exp(-.2_dp))
  call check(abs(fs-expected)<1.e-12_dp,'analytic repair')
  repair_psi_min=-.5_dp
  call capacity_update(-1._dp,24._dp,20._dp,-.1_dp,0._dp,kp,kf,kn,fs,ad,rp,rt)
  call check(abs(fs-kf)<1.e-12_dp,'dry stem blocks repair')
  call capacity_update(-.2_dp,24._dp,-5._dp,-.1_dp,0._dp,kp,kf,kn,fs,ad,rp,rt)
  call check(abs(fs-kf)<1.e-12_dp,'cold blocks repair')
  repair_rate_per_day=0;turnover_rate_per_day=.1_dp
  call capacity_update(-.2_dp,24._dp,20._dp,-.1_dp,0._dp,kp,kf,kn,fs,ad,rp,rt)
  call check(abs(kn-kp*exp(-.1_dp))<1.e-12_dp,'analytic turnover')
  call check(abs(fs/kn-kf/kp)<1.e-12_dp,'uniform turnover preserves loss fraction')
  ! At P50 repair cannot restore more than 50% functional capacity.
  turnover_rate_per_day=0;repair_rate_per_day=10;repair_psi_min=-10
  psi=p50;kp=k0;kf=.1_dp*k0;loss=.9_dp
  call step(psi,loss,kp,kf,24._dp,20._dp,0._dp,0._dp, &
    p50+0.00981_dp*height,0._dp,1._dp,e,q,ad,rp,rt,damage)
  call check(abs(psi-p50)<1.e-9_dp,'constant potential hydraulic balance')
  call check(loss>=.5_dp-1.e-10_dp,'P50 must bound repair')
  call check(abs(loss-(1._dp-kf/kp))<1.e-12_dp,'loss capacity identity')
  ! Worse historical damage remains at P50 in the absence of recovery.
  repair_rate_per_day=0;psi=p50;kp=k0;kf=.25_dp*k0;loss=.75_dp
  call step(psi,loss,kp,kf,1._dp,20._dp,0._dp,0._dp, &
    p50+0.00981_dp*height,0._dp,1._dp,e,q,ad,rp,rt,damage)
  call check(abs(loss-.75_dp)<1.e-10_dp,'historical damage must persist without recovery')
  write(*,'(A)') 'PASS module checks: growth repair gates turnover P50 and damage memory'
end program unit_checks
