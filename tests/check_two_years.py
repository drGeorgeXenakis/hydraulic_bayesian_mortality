#!/usr/bin/env python3
"""Executable two-year growth and timestep checks; Python standard library only."""
import csv,math,shutil,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def rows(d,name):
    with (d/name).open() as f:return list(csv.DictReader(f))
def run(d):subprocess.run([str(d/'hydraulic_mortality')],cwd=d,check=True)
def main():
    compiler=shutil.which('gfortran')
    if not compiler:raise SystemExit('gfortran is required')
    with tempfile.TemporaryDirectory(prefix='two_year_hydraulics_') as tmp:
        d=Path(tmp)
        for f in ['hydraulics.f90','mortality.f90']:shutil.copy(ROOT/f,d/f)
        for f in (ROOT/'examples/two_years').iterdir():shutil.copy(f,d/f.name)
        settings=d/'bayes_settings.nml'
        settings.write_text(settings.read_text().replace('n_particles=2000','n_particles=500'))
        subprocess.run([compiler,'-std=f95','-O2','-fcheck=all','hydraulics.f90',
                        'mortality.f90','-o','hydraulic_mortality'],cwd=d,check=True)
        run(d)
        hyd=rows(d,'hydraulic_trajectories.csv');risk=rows(d,'risk_timeseries.csv')
        balance=0
        for r in hyd:
            kp=float(r['Kpot_L_h_MPa']);kf=float(r['Kfun_L_h_MPa'])
            balance=max(balance,abs(.5*(float(r['stem_MPa'])+.2)-
                                   (float(r['uptake_L'])-float(r['transpired_L']))))
            assert 0<=kf<=kp
            assert abs(float(r['PLC_percent'])/100-(1-kf/kp))<1e-10
            assert abs(kp-(1/.8+float(r['growth_added_K'])-float(r['retired_potential_K'])))<1e-8
        assert balance<1e-5
        assert any(float(r['growth_added_K'])>0 for r in hyd)
        assert all(float(r['repair_added_K'])==0 for r in hyd)
        for member in (1,2,3):
            data=[r for r in hyd if int(r['p50_member'])==member]
            assert any(float(y['PLC_percent'])<float(x['PLC_percent'])-1e-8
                       for x,y in zip(data,data[1:]))
        hp=(d/'parameters.nml').read_text()
        (d/'parameters.nml').write_text(hp.replace('max_step=0.10','max_step=0.05'))
        run(d);fine=rows(d,'risk_timeseries.csv')
        difference=max(abs(float(a['mean_Pdeath'])-float(b['mean_Pdeath']))
                       for a,b in zip(risk,fine))
        assert difference<.01
        (d/'parameters.nml').write_text(hp.replace('growth_rate_per_day=0.001','growth_rate_per_day=0.0'))
        run(d);off=rows(d,'hydraulic_trajectories.csv')
        assert max(abs(float(r['Kpot_L_h_MPa'])-1/.8) for r in off)<1e-10
        assert any(float(a['Kfun_L_h_MPa'])>float(b['Kfun_L_h_MPa'])+1e-6
                   for a,b in zip(hyd,off))
        p=[float(r['mean_Pdeath']) for r in risk]
        assert all(b>=a-1e-10 for a,b in zip(p,p[1:]))
        print('PASS two-year growth, capacity accounting, water balance and continuous mortality')
        print('Maximum balance residual L:',balance)
        print('Maximum mean risk difference after halving timestep:',difference)
if __name__=='__main__':main()
