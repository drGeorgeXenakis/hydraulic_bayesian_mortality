! F95 Bayesian importance sampling, conditional on fixed hydraulic traits
! except for a categorical P50 prior. All time integrals for mortality use DAYS.
module mortality_tools
  use hydraulics
  implicit none
  integer, parameter :: max_grid=8
  integer :: exposure_mode=1,n_particles=8000,seed_value=271828,nq=3,np50=3
  logical :: use_observations=.false.
  real(dp) :: psi_scale=1.0_dp,loss_threshold=0.5_dp
  real(dp) :: lambda0_median=1.e-4_dp,lambda0_logsd=1.5_dp
  real(dp) :: k_median=0.05_dp,k_logsd=1.0_dp
  real(dp) :: q_values(max_grid)=(/1._dp,2._dp,3._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
  real(dp) :: q_prior(max_grid)=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp/)
  real(dp) :: p50_offsets(max_grid)=(/-.2_dp,0._dp,.2_dp,0._dp,0._dp,0._dp,0._dp,0._dp/)
  real(dp) :: p50_prior(max_grid)=(/1._dp,2._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp/)
  namelist /bayes/ exposure_mode,psi_scale,loss_threshold,n_particles,seed_value, &
    use_observations,lambda0_median,lambda0_logsd,k_median,k_logsd, &
    nq,q_values,q_prior,np50,p50_offsets,p50_prior
contains
  subroutine commas(line)
    character(*), intent(inout) :: line
    integer :: j
    do j=1,len_trim(line)
      if(line(j:j)==',') line(j:j)=' '
    end do
  end subroutine commas

  subroutine init_rng()
    integer :: n,j
    integer, allocatable :: seeds(:)
    call random_seed(size=n)
    allocate(seeds(n))
    do j=1,n
      seeds(j)=mod(seed_value,1000000)+37*j
    end do
    call random_seed(put=seeds)
    deallocate(seeds)
  end subroutine init_rng

  real(dp) function normal_draw()
    real(dp) :: u,v
    call random_number(u)
    call random_number(v)
    normal_draw=sqrt(-2._dp*log(max(u,tiny(1._dp))))*cos(6.283185307179586_dp*v)
  end function normal_draw

  integer function category(prob,n)
    integer, intent(in) :: n
    real(dp), intent(in) :: prob(n)
    real(dp) :: u,s
    integer :: j
    call random_number(u)
    s=0
    category=n
    do j=1,n
      s=s+prob(j)
      if(u<s) then
        category=j
        return
      end if
    end do
  end function category

  real(dp) function death_probability(haz)
    real(dp), intent(in) :: haz
    if(haz<1.e-5_dp) then
      death_probability=haz*(1._dp-haz/2._dp+haz*haz/6._dp)
    else
      death_probability=1._dp-exp(-haz)
    end if
  end function death_probability

  real(dp) function log_death(haz)
    real(dp), intent(in) :: haz
    if(haz<=0._dp) then
      log_death=-huge(1._dp)/100._dp
    else
      log_death=log(max(death_probability(haz),tiny(1._dp)))
    end if
  end function log_death

  ! Sort values and weights together, then report discrete weighted quantiles.
  recursive subroutine sort_pairs(a,w,left,right)
    real(dp), intent(inout) :: a(:),w(:)
    integer, intent(in) :: left,right
    integer :: i,j
    real(dp) :: pivot,tmp
    if(left>=right) return
    pivot=a((left+right)/2)
    i=left; j=right
    do
      do while(a(i)<pivot)
        i=i+1
      end do
      do while(a(j)>pivot)
        j=j-1
      end do
      if(i<=j) then
        tmp=a(i); a(i)=a(j); a(j)=tmp
        tmp=w(i); w(i)=w(j); w(j)=tmp
        i=i+1; j=j-1
      end if
      if(i>j) exit
    end do
    if(left<j) call sort_pairs(a,w,left,j)
    if(i<right) call sort_pairs(a,w,i,right)
  end subroutine sort_pairs

  subroutine summary(values,weights,mean,lo,median,hi)
    real(dp), intent(in) :: values(:),weights(:)
    real(dp), intent(out) :: mean,lo,median,hi
    real(dp) :: a(size(values)),w(size(values)),cw
    integer :: j,n
    n=size(values)
    mean=sum(values*weights)
    a=values; w=weights
    call sort_pairs(a,w,1,n)
    lo=a(n);median=a(n);hi=a(n);cw=0
    do j=1,n
      cw=cw+w(j)
      if(cw>=.025_dp) then
        lo=a(j);exit
      end if
    end do
    cw=0
    do j=1,n
      cw=cw+w(j)
      if(cw>=.5_dp) then
        median=a(j);exit
      end if
    end do
    cw=0
    do j=1,n
      cw=cw+w(j)
      if(cw>=.975_dp) then
        hi=a(j);exit
      end if
    end do
  end subroutine summary
end module mortality_tools

program hydraulic_mortality
  use mortality_tools
  implicit none
  integer :: ios,j,k,z,l,b,ncase,maxrow,ii,nn,nbin,nobs
  integer, allocatable :: counts(:),ntrees(:),ndead(:),iq(:),ip(:)
  character(64), allocatable :: ids(:)
  character(256), allocatable :: paths(:)
  character(2048) :: line
  character(32) :: inference
  ! hyd columns: time,T,theta,VPD,day,soil,psi,loss,g,E,Q,R,water_out,water_in
  real(dp), allocatable :: hyd(:,:,:,:),dose(:,:,:,:),lam(:),kap(:),weights(:)
  real(dp), allocatable :: logw(:),vals(:),uniform(:),effective_p50(:)
  real(dp) :: base_p50,time,temp,theta,vpd,day,soil,psi,loss,previous,dt,h
  real(dp) :: e,q,g,r,water_out,water_in,exposure,cumdose(max_grid),duration
  real(dp) :: haz,ess,maxlog,mu,lo,med,hi,logz,balance,max_balance,elapsed

  open(10,file='parameters.nml',status='old',action='read',iostat=ios)
  call check(ios==0,'cannot open parameters.nml')
  read(10,nml=model,iostat=ios)
  call check(ios==0,'invalid hydraulic namelist')
  close(10)
  open(10,file='bayes_settings.nml',status='old',action='read',iostat=ios)
  call check(ios==0,'cannot open bayes_settings.nml')
  read(10,nml=bayes,iostat=ios)
  call check(ios==0,'invalid Bayesian namelist')
  close(10)
  call check(c>0 .and. r_other>=0 .and. r_xylem>0,'invalid capacitance/resistance')
  call check(height>=0 .and. leaf_area>0 .and. pressure>0,'invalid height/leaf area/pressure')
  call check(gmin>=0 .and. gmax>=gmin,'invalid conductance bounds')
  call check(twidth>0 .and. d0>0 .and. sg>0 .and. sx>0,'invalid response width')
  call check(theta_r>=0 .and. theta_s<=1 .and. theta_s>theta_r,'invalid moisture bounds')
  call check(alpha>0 .and. vn>1 .and. max_step>0,'invalid soil parameters or step')
  call check(initial_loss>=0 .and. initial_loss<=1,'invalid initial loss')
  call check(initial_psi<=0 .and. initial_psi>-1000,'invalid initial potential')
  call check(loss_floor>0 .and. loss_floor<1,'invalid conductance floor')
  call check(n_particles>=100 .and. n_particles<=1000000,'invalid particle count')
  call check(nq>=1 .and. nq<=max_grid .and. np50>=1 .and. np50<=max_grid,'invalid grid sizes')
  call check(all(q_values(1:nq)>0) .and. all(q_values(1:nq)<=10),'require 0 < q <= 10')
  call check(all(q_prior(1:nq)>=0) .and. sum(q_prior(1:nq))>0,'invalid q prior')
  call check(all(p50_prior(1:np50)>=0) .and. sum(p50_prior(1:np50))>0,'invalid P50 prior')
  call check(all(p50+p50_offsets(1:np50)<0),'P50 members must be negative')
  call check(exposure_mode==1 .or. exposure_mode==2,'exposure_mode must be 1 or 2')
  call check(psi_scale>0 .and. psi_scale<=100,'invalid MPa exposure scale')
  call check(loss_threshold>=0 .and. loss_threshold<1,'invalid PLC exposure threshold')
  call check(lambda0_median>0 .and. k_median>0,'prior medians must be positive')
  call check(lambda0_logsd>=0 .and. lambda0_logsd<=5,'invalid baseline log SD')
  call check(k_logsd>=0 .and. k_logsd<=5,'invalid stress log SD')
  call check(seed_value>=0,'seed must be nonnegative')
  q_prior(1:nq)=q_prior(1:nq)/sum(q_prior(1:nq))
  p50_prior(1:np50)=p50_prior(1:np50)/sum(p50_prior(1:np50))
  base_p50=p50

  open(11,file='cases.csv',status='old',action='read',iostat=ios)
  call check(ios==0,'cannot open cases.csv')
  read(11,'(A)',iostat=ios) line
  call check(ios==0,'missing cases header')
  ncase=0
  do
    read(11,'(A)',iostat=ios) line
    if(ios<0) exit
    call check(ios==0 .and. len_trim(line)>0,'invalid cases row')
    ncase=ncase+1
  end do
  call check(ncase>0,'no cases')
  rewind(11)
  read(11,'(A)') line
  allocate(counts(ncase),ntrees(ncase),ndead(ncase),ids(ncase),paths(ncase))
  do j=1,ncase
    read(11,'(A)') line
    call commas(line)
    read(line,*,iostat=ios) ids(j),paths(j),ntrees(j),ndead(j)
    call check(ios==0,'invalid cases columns')
    call check(len_trim(ids(j))>0,'empty case ID')
    if(j>1) call check(.not.any(ids(1:j-1)==ids(j)),'duplicate case ID')
    call check(ntrees(j)>=0 .and. ndead(j)>=0 .and. ndead(j)<=ntrees(j),'invalid tree counts')
    open(12,file=trim(paths(j)),status='old',action='read',iostat=ios)
    call check(ios==0,'cannot open forcing: '//trim(paths(j)))
    read(12,'(A)',iostat=ios) line
    call check(ios==0,'missing forcing header')
    counts(j)=0
    do
      read(12,'(A)',iostat=ios) line
      if(ios<0) exit
      call check(ios==0 .and. len_trim(line)>0,'invalid forcing row')
      counts(j)=counts(j)+1
    end do
    close(12)
    call check(counts(j)>=2,'need at least two forcing rows')
  end do
  close(11)
  nobs=sum(ntrees)
  if(use_observations) call check(nobs>0,'posterior mode needs survival counts')
  inference='prior_predictive'
  if(use_observations) inference='posterior'
  maxrow=maxval(counts)
  allocate(hyd(14,maxrow,ncase,np50),dose(nq,maxrow,ncase,np50))
  allocate(effective_p50(np50))
  hyd=0; dose=0; max_balance=0
  do b=1,np50
    p50=base_p50+p50_offsets(b)
    effective_p50(b)=p50
    do j=1,ncase
      open(12,file=trim(paths(j)),status='old',action='read')
      read(12,'(A)') line
      psi=initial_psi; loss=max(initial_loss,vulnerability(psi))
      water_out=0; water_in=0; cumdose=0; previous=0
      do k=1,counts(j)
        read(12,'(A)') line
        call commas(line)
        read(line,*,iostat=ios) time,temp,theta,vpd,day
        call check(ios==0,'invalid forcing columns: '//trim(paths(j)))
        call check(time>=0 .and. time<1.e7_dp,'invalid forcing time')
        call check(temp>-100 .and. temp<100 .and. vpd>=0 .and. vpd<100,'invalid T/VPD')
        call check(day>=0 .and. day<=1,'invalid day factor')
        soil=soil_potential(theta)
        if(k>1) then
          call check(time>previous,'forcing times must strictly increase')
          dt=time-previous
          call check(dt/max_step<1.e7_dp,'interval has too many substeps')
          nn=ceiling(dt/max_step); h=dt/real(nn,dp)
          do z=1,nn
            call step(psi,loss,h,temp,vpd,day,soil,e,q)
            water_out=water_out+h*e; water_in=water_in+h*q
            if(exposure_mode==1) then
              exposure=max(0._dp,(p50-psi)/psi_scale)
            else
              exposure=max(0._dp,(loss-loss_threshold)/(1._dp-loss_threshold))
            end if
            do l=1,nq
              cumdose(l)=cumdose(l)+(h/24._dp)*exposure**q_values(l)
            end do
          end do
        end if
        call fluxes(psi,loss,temp,vpd,day,soil,g,e,q,r)
        hyd(:,k,j,b)=(/time,temp,theta,vpd,day,soil,psi,loss,g,e,q,r,water_out,water_in/)
        dose(:,k,j,b)=cumdose(1:nq)
        balance=c*(psi-initial_psi)-(water_in-water_out)
        max_balance=max(max_balance,abs(balance))
        previous=time
      end do
      close(12)
    end do
  end do
  p50=base_p50

  allocate(lam(n_particles),kap(n_particles),weights(n_particles),logw(n_particles))
  allocate(vals(n_particles),uniform(n_particles),iq(n_particles),ip(n_particles))
  uniform=1._dp/real(n_particles,dp)
  call init_rng()
  do ii=1,n_particles
    lam(ii)=exp(log(lambda0_median)+lambda0_logsd*normal_draw())
    kap(ii)=exp(log(k_median)+k_logsd*normal_draw())
    iq(ii)=category(q_prior(1:nq),nq)
    ip(ii)=category(p50_prior(1:np50),np50)
    logw(ii)=0
    if(.not.use_observations) cycle
    do j=1,ncase
      if(ntrees(j)==0) cycle
      k=counts(j)
      duration=(hyd(1,k,j,ip(ii))-hyd(1,1,j,ip(ii)))/24._dp
      haz=lam(ii)*duration+kap(ii)*dose(iq(ii),k,j,ip(ii))
      ! Binomial combinatorial factor is constant in parameters, omitted.
      if(ndead(j)>0) logw(ii)=logw(ii)+real(ndead(j),dp)*log_death(haz)
      logw(ii)=logw(ii)-real(ntrees(j)-ndead(j),dp)*haz
    end do
  end do
  maxlog=maxval(logw)
  weights=exp(logw-maxlog)
  call check(sum(weights)>0,'all importance weights vanished')
  logz=maxlog+log(sum(weights)/real(n_particles,dp))
  weights=weights/sum(weights)
  ess=1._dp/sum(weights*weights)
  open(13,file='diagnostics.txt',status='replace',action='write')
  write(13,'(A)') 'Inference: '//trim(inference)
  write(13,'(A,I10)') 'Particles: ',n_particles
  write(13,'(A,F14.3)') 'Importance ESS: ',ess
  write(13,'(A,F14.6)') 'ESS / particles: ',ess/real(n_particles,dp)
  write(13,'(A,ES16.7)') 'Maximum weight: ',maxval(weights)
  write(13,'(A,ES16.7)') 'Log mean likelihood (binomial constants omitted): ',logz
  write(13,'(A,ES16.7)') 'Maximum absolute hydraulic water balance residual L: ',max_balance
  write(13,'(A)') 'Bands condition on fixed hydraulic traits except the specified P50 grid.'
  write(13,'(A)') 'Importance sampling is an approximation; repeat with more particles and other seeds.'
  if(ess<max(200._dp,.05_dp*real(n_particles,dp))) then
    write(13,'(A)') 'WARNING: low ESS. Posterior summaries may be unreliable.'
    write(*,'(A)') 'WARNING: low importance ESS; see diagnostics.txt.'
  end if
  close(13)

  open(14,file='parameter_particles.csv',status='replace',action='write')
  write(14,'(A)') 'particle,lambda0_per_day,k_per_day,q,P50_MPa,prior_weight,weight'
  do ii=1,n_particles
    write(14,'(I10,6(",",ES23.15))') ii,lam(ii),kap(ii),q_values(iq(ii)), &
      effective_p50(ip(ii)),uniform(ii),weights(ii)
  end do
  close(14)

  open(14,file='parameter_summary.csv',status='replace',action='write')
  write(14,'(A)') 'parameter,mean,p025,median,p975'
  call summary(lam,weights,mu,lo,med,hi)
  write(14,'(A,4(",",ES23.15))') 'lambda0_per_day',mu,lo,med,hi
  call summary(kap,weights,mu,lo,med,hi)
  write(14,'(A,4(",",ES23.15))') 'k_per_day',mu,lo,med,hi
  do ii=1,n_particles
    vals(ii)=q_values(iq(ii))
  end do
  call summary(vals,weights,mu,lo,med,hi)
  write(14,'(A,4(",",ES23.15))') 'q',mu,lo,med,hi
  do ii=1,n_particles
    vals(ii)=effective_p50(ip(ii))
  end do
  call summary(vals,weights,mu,lo,med,hi)
  write(14,'(A,4(",",ES23.15))') 'P50_MPa',mu,lo,med,hi
  close(14)

  open(14,file='hydraulic_trajectories.csv',status='replace',action='write')
  write(14,'(A)') 'case_id,p50_member,P50_MPa,time_h,T_C,theta,VPD_kPa,day_factor,'// &
    'soil_MPa,stem_MPa,PLC_percent,gsw_mol_m2_s,E_L_h,Qin_L_h,R_MPa_h_L,transpired_L,uptake_L'
  open(15,file='exposure_doses.csv',status='replace',action='write')
  write(15,'(A)') 'case_id,p50_member,q,time_h,dose_days'
  do j=1,ncase
    do b=1,np50
      do k=1,counts(j)
        write(14,'(A,",",I4,15(",",ES23.15))') trim(ids(j)),b,effective_p50(b), &
          hyd(1:7,k,j,b),100._dp*hyd(8,k,j,b),hyd(9:14,k,j,b)
        do l=1,nq
          write(15,'(A,",",I4,3(",",ES23.15))') trim(ids(j)),b,q_values(l), &
            hyd(1,k,j,b),dose(l,k,j,b)
        end do
      end do
    end do
  end do
  close(14);close(15)

  open(14,file='risk_timeseries.csv',status='replace',action='write')
  write(14,'(A)') 'case_id,inference,time_h,elapsed_days,mean_Pdeath,p025,median,p975,mean_survival'
  open(15,file='endpoint_predictions.csv',status='replace',action='write')
  write(15,'(A)') 'case_id,inference,n_trees,n_dead,mean_Pdeath,p025,median,p975,'// &
    'replicate_deaths_p025,replicate_deaths_median,replicate_deaths_p975'
  do j=1,ncase
    do k=1,counts(j)
      elapsed=(hyd(1,k,j,1)-hyd(1,1,j,1))/24._dp
      do ii=1,n_particles
        haz=lam(ii)*elapsed+kap(ii)*dose(iq(ii),k,j,ip(ii))
        vals(ii)=death_probability(haz)
      end do
      call summary(vals,weights,mu,lo,med,hi)
      write(14,'(A,",",A,7(",",ES23.15))') trim(ids(j)),trim(inference), &
        hyd(1,k,j,1),elapsed,mu,lo,med,hi,1._dp-mu
    end do
    write(15,'(A,",",A,2(",",I10),4(",",ES23.15))',advance='no') &
      trim(ids(j)),trim(inference),ntrees(j),ndead(j),mu,lo,med,hi
    ! Replicate endpoint counts: mixture of binomials, unconditioned on realised
    ! counts in each replicate. Sample from weighted parameter particles.
    if(ntrees(j)>0) then
      do ii=1,n_particles
        b=category(weights,n_particles)
        haz=lam(b)*elapsed+kap(b)*dose(iq(b),counts(j),j,ip(b))
        exposure=death_probability(haz)
        nbin=0
        do z=1,ntrees(j)
          call random_number(h)
          if(h<exposure) nbin=nbin+1
        end do
        vals(ii)=real(nbin,dp)
      end do
      call summary(vals,uniform,mu,lo,med,hi)
      write(15,'(3(",",ES23.15))') lo,med,hi
    else
      write(15,'(A)') ',NA,NA,NA'
    end if
  end do
  close(14);close(15)
  write(*,'(A)') 'Finished: '//trim(inference)//'. See diagnostics.txt and risk_timeseries.csv.'
end program hydraulic_mortality
