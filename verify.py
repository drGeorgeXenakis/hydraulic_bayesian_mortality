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
        assert all(y >= x-1e-9 for x, y in zip(plc, plc[1:]))
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
        # Alternative PLC exposure: finite monotone risk and memory across rewetting.
        (work/'bayes_settings.nml').write_text(settings.replace('exposure_mode=1','exposure_mode=2'))
        (work/'cases.csv').write_text((ROOT/'cases.csv').read_text())
        run(work)
        check_outputs(work)
        print('PASS: compiled F95, water balance, persistent damage, monotone probabilities,')
        print('Bayesian likelihood weights, analytic baseline survival, PLC mode, timestep convergence.')
        print('Maximum mean-risk difference after halving timestep:',diff)

if __name__=='__main__':
    main()
