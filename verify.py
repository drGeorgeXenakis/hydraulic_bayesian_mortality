#!/usr/bin/env python3
"""Compile and check the actual Fortran executable; standard library only."""
import csv
import math
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def rows(folder, name):
    with (folder/name).open() as f:
        return list(csv.DictReader(f))

def run(folder):
    subprocess.run([str(folder/'hydraulic_mortality')], cwd=folder, check=True)

def check_outputs(folder):
    trajectories = defaultdict(list)
    for row in rows(folder, 'hydraulic_trajectories.csv'):
        trajectories[row['case_id'], row['p50_member']].append(row)
        psi = float(row['stem_MPa'])
        balance = .5*(psi+.2)-(float(row['uptake_L'])-float(row['transpired_L']))
        assert abs(balance) < 2e-6, balance
    for data in trajectories.values():
        plc = [float(r['PLC_percent']) for r in data]
        assert all(-1e-10 <= x <= 100+1e-10 for x in plc)
        # PLC can decline with growth or repair in version 2. Every row
        # must still respect both the capacity identity and P50 curve.
        for r in data:
            kp=float(r['Kpot_L_h_MPa']);kf=float(r['Kfun_L_h_MPa'])
            loss=float(r['PLC_percent'])/100
            assert kp>0 and -1e-12<=kf<=kp+1e-12
            assert abs(loss-(1-kf/kp))<1e-11
            x=(float(r['stem_MPa'])-float(r['P50_MPa']))/.35
            curve=math.exp(-x)/(1+math.exp(-x)) if x>=0 else 1/(1+math.exp(x))
            assert loss>=curve-1e-9
            # Potential-capacity accounting separates addition and retirement.
            expected=1/.8+float(r['growth_added_K'])-float(r['retired_potential_K'])
            assert abs(kp-expected)<1e-8
    risks = defaultdict(list)
    for row in rows(folder, 'risk_timeseries.csv'):
        risks[row['case_id']].append(row)
        mean, lo, med, hi = [float(row[x]) for x in ('mean_Pdeath','p025','median','p975')]
        assert 0 <= mean <= 1 and 0 <= lo <= med <= hi <= 1
        assert abs(mean+float(row['mean_survival'])-1) < 1e-12
    for data in risks.values():
        p = [float(r['mean_Pdeath']) for r in data]
        assert p[0] == 0
        assert all(y >= x-1e-10 for x, y in zip(p, p[1:]))
    return risks

def main():
    compiler = shutil.which('gfortran')
    if not compiler:
        raise SystemExit('gfortran is required; no executable checks were run.')
    with tempfile.TemporaryDirectory(prefix='hydraulic_mortality_check_') as td:
        work=Path(td)
        for p in ROOT.iterdir():
            if p.suffix in ('.f90','.nml','.csv'):
                shutil.copy(p,work/p.name)
        subprocess.run([compiler,'-std=f95','-Wall','-Wextra','-O2','-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow','hydraulics.f90','mortality.f90',
                        '-o','hydraulic_mortality'], cwd=work,check=True)
        settings=(work/'bayes_settings.nml').read_text().replace('n_particles=8000','n_particles=2000')
        (work/'bayes_settings.nml').write_text(settings)
        # Direct analytic module tests exercise the new equations.
        subprocess.run([compiler,'-std=f95','-O2','-fcheck=all',
                        'hydraulics.f90',str(ROOT/'tests/unit_checks.f90'),
                        '-o','unit_checks'],cwd=work,check=True)
        subprocess.run([str(work/'unit_checks')],cwd=work,check=True)
        run(work)
        prior=check_outputs(work)
        particles=rows(work,'parameter_particles.csv')
        assert abs(sum(float(r['weight']) for r in particles)-1) < 1e-10
        assert all(abs(float(r['weight'])-1/2000) < 1e-12 for r in particles)
        # Actual endpoint posterior weights versus independent binomial likelihood.
        shutil.copy(work/'cases_SYNTHETIC_calibration.csv', work/'cases.csv')
        (work/'bayes_settings.nml').write_text(settings.replace('use_observations=.false.',
                                                               'use_observations=.true.'))
        run(work)
        check_outputs(work)
        posts=rows(work,'parameter_particles.csv')
        assert any(abs(float(r['weight'])-1/2000)>1e-5 for r in posts)
        doses={}
        durations={}
        for r in rows(work,'exposure_doses.csv'):
            doses[r['case_id'],int(r['p50_member']),float(r['q'])]=float(r['dose_days'])
        for r in rows(work,'risk_timeseries.csv'):
            durations[r['case_id']]=float(r['elapsed_days'])
        p50s=sorted({float(r['P50_MPa']) for r in posts})
        counts={'control':(20,0),'moderate':(20,1),'severe':(20,6)}
        logs=[]
        for r in posts:
            ll=0
            b=min(range(len(p50s)), key=lambda i:abs(p50s[i]-float(r['P50_MPa'])))+1
            for name,(n,d) in counts.items():
                H=float(r['lambda0_per_day'])*durations[name]+float(r['k_per_day'])*doses[name,b,float(r['q'])]
                if d: ll+=d*math.log(-math.expm1(-H))
                ll-=(n-d)*H
            logs.append(ll)
        a=max(logs);ws=[math.exp(x-a) for x in logs];total=sum(ws)
        assert max(abs(w/total-float(r['weight'])) for w,r in zip(ws,posts)) < 1e-10
        # Halve the hydraulic timestep using identical prior draws.
        (work/'bayes_settings.nml').write_text(settings)
        (work/'cases.csv').write_text((ROOT/'cases.csv').read_text())
        hp=(work/'parameters.nml').read_text()
        (work/'parameters.nml').write_text(hp.replace('max_step=0.01','max_step=0.005'))
        run(work)
        fine=check_outputs(work)
        diff=max(abs(float(a['mean_Pdeath'])-float(b['mean_Pdeath']))
                 for name in prior for a,b in zip(prior[name],fine[name]))
        assert diff < .01, diff
        # Analytic no-excess-stress survival: a fixed lambda0=0.01 /day.
        (work/'bayes_settings.nml').write_text('''&bayes
 exposure_mode=1,n_particles=200,seed_value=1,use_observations=.false.,
 lambda0_median=0.01,lambda0_logsd=0,k_median=0.2,k_logsd=0,
 nq=1,q_values=2,q_prior=1,np50=1,p50_offsets=0,p50_prior=1
/\n''')
        (work/'parameters.nml').write_text(hp)
        (work/'cases.csv').write_text('case_id,forcing_file,n_trees,n_dead\ncontrol,forcing_control.csv,0,0\n')
        run(work)
        for r in rows(work,'risk_timeseries.csv'):
            expected=-math.expm1(-.01*float(r['elapsed_days']))
            assert abs(float(r['mean_Pdeath'])-expected) < 1e-11
        # Wet, damaged cohort: growth, repair, and seasonal gates.
        (work/'cases.csv').write_text('case_id,forcing_file,n_trees,n_dead\nwet,wet.csv,0,0\n')
        with (work/'wet.csv').open('w') as f:
            f.write('time_h,T_C,theta,VPD_kPa,day_factor,growth_factor,leaf_area_factor\n')
            for hour in range(0,24*30+1,6):
                f.write(f'{hour},20,0.35,0,0,1,1\n')
        wet_hp=hp.replace('initial_loss=0.0','initial_loss=0.7')
        (work/'parameters.nml').write_text(wet_hp)
        (work/'bayes_settings.nml').write_text(settings)
        run(work);check_outputs(work)
        wet=rows(work,'hydraulic_trajectories.csv')
        for b in ('1','2','3'):
            part=[r for r in wet if int(r['p50_member'])==int(b)]
            assert float(part[-1]['PLC_percent'])<float(part[0]['PLC_percent'])
            assert float(part[-1]['growth_added_K'])>0
            assert float(part[-1]['repair_added_K'])==0
        # Turn growth off and enable explicit repair; potential capacity is fixed.
        rp_hp=wet_hp.replace('growth_rate_per_day=0.001','growth_rate_per_day=0.0')
        rp_hp=rp_hp.replace('repair_rate_per_day=0.0','repair_rate_per_day=0.1')
        (work/'parameters.nml').write_text(rp_hp)
        run(work);check_outputs(work)
        repaired=rows(work,'hydraulic_trajectories.csv')
        assert max(float(r['repair_added_K']) for r in repaired)>0
        assert max(abs(float(r['Kpot_L_h_MPa'])-1/.8) for r in repaired)<1e-10
        # Both recovery mechanisms off: persistent PLC under wet forcing.
        (work/'parameters.nml').write_text(rp_hp.replace('repair_rate_per_day=0.1','repair_rate_per_day=0.0'))
        run(work)
        no_repair=rows(work,'hydraulic_trajectories.csv')
        assert max(abs(float(r['PLC_percent'])-70) for r in no_repair)<1e-8
        # Strict parsing must reject a six-column record instead of silently truncating it.
        (work/'wet.csv').write_text('time_h,T_C,theta,VPD_kPa,day_factor,growth_factor\n'
                                   '0,20,.35,0,0,1\n24,20,.35,0,0,1\n')
        failed=subprocess.run([str(work/'hydraulic_mortality')],cwd=work,capture_output=True)
        assert failed.returncode!=0
        (work/'parameters.nml').write_text(hp)
        # Alternative PLC exposure: finite monotone risk and memory across rewetting.
        (work/'bayes_settings.nml').write_text(settings.replace('exposure_mode=1','exposure_mode=2'))
        (work/'cases.csv').write_text((ROOT/'cases.csv').read_text())
        run(work)
        check_outputs(work)
        print('PASS: compiled F95, water balance, persistent damage, monotone probabilities,')
        print('Bayesian likelihood weights, analytic baseline survival, PLC mode, timestep convergence.')
        print('PASS additional v2 growth repair capacity-accounting and malformed-input checks.')
        print('Maximum mean-risk difference after halving timestep:',diff)

if __name__=='__main__':
    main()
