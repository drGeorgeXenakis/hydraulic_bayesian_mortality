# Hydraulic Bayesian mortality model version 2.0

Updated 3 November 2026. This Fortran 95 model adds xylem growth, optional
embolism repair, optional proportional xylem retirement and prescribed seasonal
leaf area. Actual conductance loss remains bounded by the P50 vulnerability
curve. The Bayesian hazard/survival equations are retained.

All supplied rates, thresholds, priors, environmental scenarios and synthetic
counts are illustrative. This is not a calibrated species mortality model.

## Compile and run

```sh
gfortran -std=f95 -Wall -Wextra -O2 -fcheck=all \
  hydraulics.f90 mortality.f90 -o hydraulic_mortality
./hydraulic_mortality
python3 verify.py
```

Or use `make run` and `make check`. Python is used only for verification, not to
run the Fortran model. Inputs are read from the working directory and output
files are replaced there. There are no command-line input arguments.

## Inputs

- `parameters.nml`: hydraulic traits, initial states, growth/repair/retirement.
- `bayes_settings.nml`: exposure formulation, priors and inference settings.
- `cases.csv`: unique IDs, forcing paths, cohort size and endpoint death count.
- `forcing_*.csv`: one header followed by numeric rows.

Five-column legacy forcing:

`time_h,T_C,theta,VPD_kPa,day_factor`

Seven-column seasonal forcing:

`time_h,T_C,theta,VPD_kPa,day_factor,growth_factor,leaf_area_factor`

Five columns imply growth_factor=0 and leaf_area_factor=1. Six columns, blank
fields and nonnumeric fields are rejected. Growth requires both a seasonal
factor in [0,1] and warm/wet eligibility. Leaf area equals the base parameter
times leaf_area_factor (0 to 100). Drivers are constant over the interval ending
at each row; the first row only supplies initial output. No calendar or
phenology model is embedded.

## New equations and parameters

Initial potential conductance K0=1/r_xylem. Actual loss L=1-Kfun/Kpot.
For an internal step of d days:

1. Retire both conductance pools proportionally: multiply by exp(-mu*d).
2. Add new functional conductance G=K0*gamma*growth_factor*growth_gate*d to both.
3. Restore J=(Kpot-Kfun)*(1-exp(-r*repair_gate*d)) to functional capacity only.
4. Solve the hydraulic balance and enforce:
   Kfun_new=min(Kfun_provisional,Kpot_new*(1-Lcurve(Psi_new))).
5. Update L=1-Kfun_new/Kpot_new, fluxes, cumulative exposure and hazard.

Lcurve(Psi)=1/(1+exp((Psi-P50)/sx)). Thus L_new is never below the loss implied
by current Psi and P50. Growth/repair can lower actual L under wetter conditions,
but do not reset mortality hazard accumulated in earlier years.

New entries in `parameters.nml`:

| Name | Supplied value | Meaning |
|---|---|---|
| growth_rate_per_day | 0.001 | Fraction of initial capacity added per eligible day |
| turnover_rate_per_day | 0.0 | Proportional retirement per day |
| growth_soil_min | -0.5 MPa | Soil potential threshold |
| growth_psi_min | -1.0 MPa | Old stem potential threshold |
| growth_temp_min / growth_temp_max | 5 / 40 C | Inclusive growth temperature range |
| repair_rate_per_day | 0.0 | Repair per day of missing functional capacity |
| repair_soil_min | -0.2 MPa | Soil potential threshold for repair |
| repair_psi_min | -0.5 MPa | Old stem potential threshold for repair |
| repair_temp_min / repair_temp_max | 5 / 35 C | Inclusive repair temperature range |

Module defaults for all new rates are zero, preserving old namelist behaviour.
The shipped file explicitly selects growth_rate_per_day=0.001; legacy drivers
still disable growth. Repair is OFF unless deliberately enabled. Rates and
thresholds are fixed hydraulic parameters, not inferred by the Bayesian routine.

Capacitance, height and non-xylem resistance remain fixed. Growth is conductance
addition, not biomass production; water/carbon costs, xylem age cohorts and
changes in vulnerability are not simulated. Retirement removes damaged and
functional xylem proportionally. Seasonal leaf area is externally prescribed.

## Examples

The root example retains the original 14-day, five-column drivers. It defaults
to prior predictive mode, with no real mortality observations. The separately
labelled synthetic calibration file contains invented counts (0/20,1/20,6/20).

After compiling in the root directory, run the two-year scenario using:

```sh
cd examples/two_years
../../hydraulic_mortality
```

That folder contains a complete set of settings and a 730-day seasonal driver.
It has two repeated synthetic droughts, seasonal growth and leaf-area changes;
repair and retirement are disabled. It uses coarse six-hourly drivers for a
small demonstration, not a validated representation of any real experiment.
The internal hydraulic timestep is 0.10 h. Test forcing and timestep resolution
before research use. One continuous case retains both conductance pools and
cumulative hazard across years; separate cases or runs reset initial states.

## Outputs and compatibility

All seven output filenames from version 1 are retained. The first 17 columns of
`hydraulic_trajectories.csv` retain their order. Eight columns are appended:

- Kpot_L_h_MPa, Kfun_L_h_MPa: potential and functional conductance.
- growth_factor, leaf_area_factor: supplied seasonal factors.
- growth_added_K, repair_added_K: cumulative gross capacity additions.
- retired_potential_K: cumulative potential capacity removed by retirement.
- new_damage_K: cumulative functional capacity removed by new embolism.

Capacity totals use conductance units L h^-1 MPa^-1, not litres. Gross repair
can be followed by immediate embolism, so it is not necessarily a net gain.
Actual PLC may now decrease, but cumulative death probability cannot decrease.
Check Kpot=K0+growth_added_K-retired_potential_K and L=1-Kfun/Kpot.

The Bayesian calculation is prior importance sampling over lambda0, k, a
categorical q and a categorical P50. The same P50 controls embolism and mortality.
`use_observations=.false.` gives prior predictive risk. To calibrate, enter actual
endpoint counts in cases.csv and set it to true. Prediction-only cases use N=0.
Do not enter repeated annual endpoints from the same cohort as independent cases.
Recalibrate mortality after changing recovery assumptions.

Manual copies are included in `docs/`.

See the updated Word/PDF manual for full theory, 28 numbered equations, all
original and new parameters, references, outputs, numerical assumptions and
calibration instructions. `VERIFICATION.md` records checks actually executed.
