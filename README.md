# Hydraulic–Bayesian mortality model — Fortran 95

This package links a soil–tree resistance–capacitance (RC) model to the Bayesian
hydraulic survival model proposed in the earlier discussion. It includes
stomatal closure, residual transpiration, persistent embolism, hydraulic exposure
integration, prior predictive simulation and Bayesian updating from observed
endpoint survival counts.

**All example inputs and prior parameters are illustrative. The outputs are not
calibrated species mortality probabilities.** Default mode does not read any
survival observations into the likelihood. Synthetic calibration counts are
provided separately and must never be reported as experimental data.

## 1. Compile and run

Requires a Fortran 95 compiler. No external Fortran libraries are needed.
Run the following from this extracted directory:

```sh
gfortran -std=f95 -Wall -Wextra -O2 -fcheck=all \
  hydraulics.f90 mortality.f90 -o hydraulic_mortality
./hydraulic_mortality
```

Alternatively, `make run`. Output files in the working directory are replaced.
Use a separate directory for each experiment and archive inputs with outputs.
The code uses fixed input filenames to remain within Fortran 95; no
GET_COMMAND_ARGUMENT extension is required. Python is needed only for the
optional verification script, not to run the model.

## 2. How the models are linked

The hydraulic module predicts stem potential Psi(t) and fractional loss of
xylem conductance L(t) from temperature, moisture, VPD and a daylight factor.
The mortality module then integrates a hazard along the living-tree trajectory.
It does not change the hydraulic trajectory when a tree dies or randomly kill
individual trees during the simulation.

### Original P50 depth–duration formulation (default exposure_mode=1)

x(t) = max[0, (P50 - Psi_stem(t)) / psi_scale]

D_q(t) = integral x(u)^q du, with u in DAYS

lambda(t) = lambda0 + k*x(t)^q

H(t) = lambda0*t + k*D_q(t)

S(t) = exp[-H(t)]

P(death by t) = 1 - S(t)

psi_scale defaults to 1 MPa. With negative water potentials, Psi=-4 MPa and
P50=-3 MPa gives x=1 when psi_scale=1 MPa. More negative than P50 means greater
stress. P50 is a predictor threshold here, not a universal biological death
threshold. lambda0 and k have units day^-1 because x is dimensionless.
q controls the severity response, while integration supplies duration.
No extra multiplication of an exceedance probability by duration is needed.

The baseline hazard is constant; the total hazard is time varying because x
changes. The survival/cumulative-hazard identity is standard survival analysis:
https://mc-stan.org/docs/stan-users-guide/survival.html
The particular hydraulic exposure function is a modelling hypothesis to test.

### Optional persistent PLC exposure (exposure_mode=2)

x(t) = max[0, (L(t) - loss_threshold) / (1-loss_threshold)]

Use the same D_q, hazard and survival equations. The default threshold is
L=0.5 (50% loss), giving x between zero and one. Because L remembers the worst
stress, this model can retain excess hazard after rewetting. The original
potential-based hazard returns to baseline once Psi rises above P50, although
cumulative mortality risk never decreases. Select the exposure formulation
before calibration and use priors appropriate to its scale. The two exposure
models are alternatives; their terms are not added together automatically.

No recovery/repair of embolism, explicit delayed-death compartment, carbon
starvation, growth, pests, or competition is represented.

## 3. Bayesian inference implemented here

The program draws independent parameter particles from these priors:

- lambda0 ~ LogNormal(log(lambda0_median), lambda0_logsd)
- k ~ LogNormal(log(k_median), k_logsd)
- q ~ categorical(q_values, q_prior)
- P50 ~ categorical(base_p50 + p50_offsets, p50_prior)

The baseline and stress coefficients have continuous priors. The q and P50
priors are **finite discrete distributions**, not approximated continuous
posteriors. At most eight members are supported per grid. The same sampled
P50 controls both hydraulic embolism and the mortality exposure threshold.
Each P50 member receives a separate hydraulic simulation. Other hydraulic
traits, soil retention parameters and forcing inputs are fixed.

Priors are independent before observing survival outcomes. The P50 prior
represents uncertainty about a shared species parameter, not individual-tree
variation. All cases in one run share the mortality parameters and the P50
member. Run separate directories for different species; this is not a
hierarchical multispecies model. The chosen hydraulic traits should describe
the species and size class under study.

For observed case j, all N_j trees are alive at the first forcing time and
follow the same representative environment. Given final cumulative hazard H_j,

p_j = 1-exp(-H_j)

n_dead_j ~ Binomial(N_j, p_j)

log likelihood = sum_j [n_dead_j*log(p_j) - (N_j-n_dead_j)*H_j]

The binomial coefficient, constant in the parameters, is omitted. Trees are
assumed conditionally independent; block effects and individual frailty are
not represented. This models
deaths at unknown times between the first and last row, and survival until the
last row for the remaining trees. Different cases may have different end times.
Exact death times and repeated survival inspections are not used. Do not enter
repeated observations of the same cohort as independent cases. For individual
trajectories, use one case per tree with N=1 and n_dead=0 or 1.

Posterior particle weight is proportional to this likelihood, since particles
were sampled from the prior. Weights are normalised with a log-sum-exp
calculation. This is Bayesian **prior importance sampling**, not MCMC. No
resampling or hidden fitting package is used. Summary quantiles use weighted
empirical inverse CDFs. The output uncertainty bands describe uncertainty in
risk conditional on the specified model and priors, not empirical confidence
limits or a guarantee about future outcomes.

Default `use_observations=.false.` gives equal weights and explicitly labels
outputs `prior_predictive`. Setting it to `.true.` requires positive observed
tree counts, incorporates those observations and labels outputs `posterior`.

Check diagnostics.txt. ESS=1/sum(weight^2) measures weight concentration.
ESS below max(200, 0.05*n_particles) triggers a warning. This threshold is a
heuristic, not a guarantee of accuracy. Repeat with more particles and several
seeds and examine stability; poor prior coverage cannot be repaired simply by
trusting one large weight. Strongly concentrated posteriors may require a more
efficient sampler. Very sparse data may not distinguish lambda0, k, q and P50;
use independently measured traits and well-justified priors. This implementation
does not fit C, R, gsw response parameters, or forcing uncertainty.

## 4. Files you need to supply

| File | Purpose |
|---|---|
| parameters.nml | Shared hydraulic, plant and soil traits |
| bayes_settings.nml | Exposure choice, uncertainty priors and Monte Carlo settings |
| cases.csv | Case IDs, forcing paths and optional observed endpoint counts |
| forcing_*.csv | Environmental time series for each case |

### Forcing drivers

Exactly one header followed by five numeric columns in this order:

| Column | Units / meaning |
|---|---|
| time_h | Strictly increasing elapsed hours; any nonnegative starting value |
| T_C | Air temperature, degrees Celsius |
| theta | Root-zone volumetric water content, m3/m3, NOT percentage |
| VPD_kPa | Air vapour pressure deficit, kPa |
| day_factor | 0 at night, 1 in daylight, or supplied light response in [0,1] |

Every row after the first supplies constant forcing for the preceding interval
ending at that row's time. The first row only establishes the initial output.
No interpolation is performed. Mortality elapsed time is measured from the
first forcing time; cumulative hazard is zero there regardless of that timestamp.
Use enough temporal resolution to capture daily extremes and drying events.
Derive day_factor externally from sunrise/sunset and location, or calibrate a
light response from radiation. Soil moisture must satisfy
`theta_r < theta <= theta_s`; invalid moisture is rejected rather than clipped.
Missing values, blank rows, quoted numeric fields and embedded commas are
unsupported. Dates must be converted to elapsed hours before running.

The supplied control, moderate and severe drivers cover 14 days in half-hourly
steps, with eight days of drying, two days held at minimum moisture, then rewetting. These are synthetic
scenarios, not reproductions of your experimental treatment schedules. Their
minimum moisture values are 0.35, 0.057 and 0.055 m3/m3 respectively. A prescribed
moisture trajectory does not require a rainfall or drainage model.

### Case table

`case_id,forcing_file,n_trees,n_dead`

Case IDs must be unique, without spaces or commas. Use simple relative filenames
without spaces; quote paths containing `/` in the CSV because Fortran
list-directed input treats an unquoted slash as a terminator. For example:

```csv
case_id,forcing_file,n_trees,n_dead
calibration_A,"drivers/forcing_A.csv",30,4
prediction_B,forcing_B.csv,0,0
```

N=0 marks a prediction-only case and requires n_dead=0. N>0 supplies observed
counts only when use_observations is true. All cases are still predicted in
either mode. Prediction-only cases have no replicate-death-count interval,
because no prediction cohort size is supplied. For observed cohorts, predictive
counts describe a NEW replicate cohort of the same size, not unobserved extra
deaths in the already observed cohort.

Follow-up should include an appropriate post-rewatering assessment period.
If survival is assessed later than the forcing ends, extend the forcing through
that assessment date. Baseline risk continues over the whole follow-up period.

## 5. Hydraulic parameters and equations

The original one-compartment RC circuit is retained:

C*dPsi/dt = (Psi_soil - 0.00981*height - Psi)/R - E

R = r_other + r_xylem / max(1-L, loss_floor)

L_new = max[L_old, 1/(1+exp((Psi-P50)/sx))]

gsw = gmin + (gmax-gmin)*day_factor*exp(-((T-topt)/twidth)^2)
      * 1/(1+VPD/d0) * 1/(1+exp((gs50-Psi)/sg))

E = 64.8*leaf_area*gsw*VPD/pressure

| Parameter | Units / meaning |
|---|---|
| c | Constant water capacitance, L/MPa |
| r_other | Constant non-xylem series resistance, MPa h/L |
| r_xylem | Fully functional xylem resistance, MPa h/L |
| height | Effective stem-node elevation above soil, metres |
| leaf_area | Total tree leaf area, m2 |
| pressure | Atmospheric pressure, kPa |
| gmin, gmax | Residual and maximum effective conductance, mol/m2/s |
| topt, twidth | Conductance temperature optimum and Gaussian width, degrees C |
| d0 | VPD response scale, kPa |
| gs50, sg | Stomatal half-response potential and positive width, MPa |
| p50, sx | Base embolism half-loss potential and positive width, MPa |
| theta_r, theta_s | Residual and saturated volumetric moisture, m3/m3 |
| alpha, vn | van Genuchten parameters, m^-1 and dimensionless |
| max_step | Maximum hydraulic integration substep, HOURS |
| initial_psi | Initial stem water potential, MPa |
| initial_loss | Initial fractional conductance loss [0,1] |
| loss_floor | Numerical lower bound on remaining conductance |

For mortality PLC exposure, `loss_threshold` in bayes_settings.nml is used.
The initial damage is max(initial_loss, vulnerability(initial_psi)). Each
P50 member can therefore have slightly different initial damage. Initial
potential, storage and tree size are otherwise identical across all cases.

The inverse van Genuchten retention relation is:

Se=(theta-theta_r)/(theta_s-theta_r); m=1-1/vn

Psi_soil= -0.00981/alpha * [Se^(-1/m)-1]^(1/vn)

Fit the soil retention curve using your measurements; replace the
soil_potential function with interpolation if preferred. Matric potential only
is represented. Soil osmotic potential is omitted. Soil/root resistance is
constant in this version, despite the effect of drying on real soil conductance.

Leaf temperature is assumed to equal air temperature; leaf boundary-layer
resistance is neglected. VPD already contains temperature's effect on vapour
pressure; the additional temperature multiplier is an empirical direct
conductance response. Keep gs50 separate from P50. Both sx and sg must be
constrained: a P50 value alone does not determine vulnerability-curve width.
No additional physiological stomatal lag is included.

Hydraulic steps use hours, while exposure/hazard use days. The conversion occurs
at each substep (h/24). Backward Euler solves the water balance by bisection with
new-state stomatal closure and damage. Re-run with max_step halved and compare
mortality probabilities as well as hydraulic trajectories, especially near
rapid embolism. The residual is continuous but need not be globally monotone
under strong embolism feedback. Capacitance has no finite-storage bound;
near-total conductance loss depends on the numerical floor. Such extreme
trajectories require physiological validation before mortality interpretation.
Signed uptake permits reverse exchange with the single soil compartment;
this is not a multilayer hydraulic redistribution model.

## 6. Bayesian settings

| Parameter | Meaning |
|---|---|
| exposure_mode | 1: original P50 deficit; 2: persistent PLC exposure |
| psi_scale | MPa scale making potential deficit dimensionless (mode 1) |
| loss_threshold | Fractional PLC threshold, [0,1), used in mode 2 |
| n_particles | Number of prior importance particles, 100 to 1,000,000 |
| seed_value | Nonnegative RNG seed; repeatable with the same compiler/runtime |
| use_observations | False: prior predictive; true: endpoint Bayesian updating |
| lambda0_median | Prior median baseline hazard, day^-1 |
| lambda0_logsd | Natural-log SD, dimensionless; 0 fixes lambda0 |
| k_median | Prior median stress-hazard coefficient, day^-1 |
| k_logsd | Natural-log SD, dimensionless; 0 fixes k |
| nq | Number of categorical q values, 1 to 8 |
| q_values | Positive exposure exponents, at most 10 |
| q_prior | Nonnegative relative probabilities, automatically normalised |
| np50 | Number of P50 prior members, 1 to 8 |
| p50_offsets | Offsets in MPa added to p50 in parameters.nml |
| p50_prior | Nonnegative relative probabilities, automatically normalised |

Supply only the first nq/np50 entries or complete arrays as needed. All
retained entries must be specified intentionally. The code's default prior
medians and widths are arbitrary demonstration values. Elicit baseline risk
from controls and stress risk from independent mortality evidence; do prior
predictive checks before fitting. Changing exposure_mode or psi_scale changes
the meaning of k and requires reconsidering its prior.

To fix independently measured P50, set `np50=1, p50_offsets=0.0, p50_prior=1.0`.
To fix q=2, set `nq=1, q_values=2.0, q_prior=1.0`.

## 7. Outputs

| File | Contents |
|---|---|
| hydraulic_trajectories.csv | Each case and P50 member: potentials, PLC, gsw, fluxes, resistance, water volumes |
| exposure_doses.csv | Cumulative dimensionless-stress exposure in days for every q and P50 member |
| parameter_particles.csv | Prior-drawn lambda0, k, q, P50, prior weights and updated weights |
| parameter_summary.csv | Weighted parameter means, medians and 95% equal-tail intervals |
| risk_timeseries.csv | Mean P(death), median, 95% risk interval and mean survival at each forcing time |
| endpoint_predictions.csv | Endpoint risks; replicate death-count intervals where N>0 |
| diagnostics.txt | Inference mode, importance ESS, weight concentration and water balance error |

Risk intervals are distributions of parameter-conditional probabilities, NOT
intervals for a binary individual outcome. Replicate-death-count intervals add
binomial outcome variability to parameter uncertainty. In posterior mode,
predictions for observed cases reuse their calibration data; they are
in-sample checks, not held-out validation. The prior_weight column remains
uniform even when posterior weights concentrate.

The diagnostics log mean likelihood omits binomial constants; do not treat it
as a fully normalised marginal likelihood for comparing different data sets.
Risk bands are pointwise, not simultaneous bands across time.

## 8. Suggested workflow with your experiment

1. Constrain soil retention, leaf area, resistance, capacitance and stomatal
   responses from hydraulic observations, keeping held-out episodes for checks.
2. Choose one species per run. Insert a measured P50 and uncertainty appropriate
   to its estimate; constrain vulnerability-curve width independently.
3. Supply the actual forcing for each treatment or individual. Extend through
   the final survival assessment after rewatering.
4. Examine prior predictions with use_observations=false. The default example
   has no real mortality observations and cannot establish absolute risk.
5. Add actual endpoint counts in cases.csv and set use_observations=true.
   Use N=0 for held-out trajectories. Run and inspect ESS and posterior stability.
6. Validate against independent cohorts or episodes. Compare potential-based
   and PLC-based exposure as separate hypotheses using appropriate priors.

For the **synthetic demonstration only**, copy cases_SYNTHETIC_calibration.csv
to cases.csv, set use_observations=true and rerun. The hand-written counts
(0/20, 1/20, 6/20) are invented to exercise Bayesian updating; they were not
sampled from a fitted model or observed in your greenhouse experiment.

Identifiability is limited with only a few deaths or one treatment. If no deaths
were observed, controls and treated cohorts can constrain upper risk bounds,
but do not establish a well-determined mortality threshold. More hydraulic
measurements do not substitute for survival observations when calibrating death.

## 9. Verification

Run `python3 verify.py` (or `make check`) where gfortran is available. It compiles
with F95 compliance and runtime checking, runs prior and synthetic-posterior
modes in temporary folders, and checks water balance, probability bounds,
monotonic cumulative risk, weighted endpoint likelihoods, an analytical no-stress
survival case, and timestep sensitivity. Only Python's standard library is used.
See VERIFICATION.md for what was actually executed in the creation environment.
