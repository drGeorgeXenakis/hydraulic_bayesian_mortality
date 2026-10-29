# Verification record

## Source checks actually completed

- Both Fortran files passed fparser's strict free-format parser and static
  analyser (and the fparser2 grammar parser).
- Maximum line length is below the standard 132-character free-format limit.
- The included verification script passed Python bytecode syntax compilation.
- No Fortran compiler was available in the creation environment. Therefore
  the Fortran executable has NOT been compiled or run here. Static parsing is
  not a substitute for compiler or executable testing.

## Independent equation-level checks actually completed

An independent Python implementation reproduced the hydraulic equations,
backward-Euler stepping and Bayesian likelihood on the supplied synthetic
forcing. Its random generator differs from Fortran's, so numerical Monte Carlo
summaries will not match exactly between implementations.

- Maximum absolute hydraulic water balance residual: 6.13e-9 L.
- Fractional embolism remained finite, within [0,1], and nondecreasing.
- Synthetic posterior weights normalised correctly; importance ESS was about
  1,702 of 8,000 prior particles using the supplied lognormal medians and grids.
- Halving the internal step from 0.01 h to 0.005 h in the severe scenario
  (P50=-3.6 MPa) changed stem potential by at most 0.00226 MPa.
- For the same severe scenario, fixed q=2, lambda0=0.0001/day and k=0.2/day,
  endpoint death probability changed by 3.67e-5 after halving the step. This
  particular sensitivity check uses k=0.2, not the prior median k=0.05.
- PLC-based exposure continued to accumulate after rewetting.
- Analytical no-stress check: lambda0=0.01/day over 14 days gives mortality
  probability 1-exp(-0.14)=0.1306417646.

These are internal numerical checks, not empirical validation. No experimental
survival or mortality data were used. The example's synthetic death counts
were manually specified for demonstration.

## Executable checks supplied for your environment

Run `python3 verify.py` from the package directory with gfortran installed.
It compiles the actual source with strict F95 and runtime checks, then tests:

1. Hydraulic mass balance and embolism memory.
2. Probability bounds and nondecreasing cumulative mortality.
3. Equal prior weights and correct posterior endpoint-binomial weights.
4. Analytic baseline-only survival with fixed hazard.
5. Timestep sensitivity of prior mean mortality.
6. The alternative PLC-based exposure mode.

The script uses temporary directories and does not overwrite your working
inputs or outputs. Its tests assume the supplied example parameter values;
keep an unmodified package copy for regression testing.
