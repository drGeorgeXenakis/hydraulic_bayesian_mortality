# Verification of version 2.0

Date: 3 November 2026. Compiler: GNU Fortran 13.3.0 (Ubuntu build).

## Executable tests completed

The actual Fortran source was compiled with `-std=f95 -Wall -Wextra -O2
-fcheck=all`; the primary regression build also used
`-ffpe-trap=invalid,zero,overflow`. No compiler warnings were emitted.
`python3 verify.py` completed successfully with 2000 importance particles in its
short-example runs, plus a fixed-hazard analytical test with 200 particles.

Passed checks:

- Analytic additions to potential and functional conductance by growth.
- Exponential repair of missing capacity and blocked repair in dry/cold states.
- Proportional retirement and its preservation of fractional loss by itself.
- P50-constrained functional conductance after provisional repair.
- Historical damage persistence when growth and repair are off.
- Water balance, positive potential capacity, bounded functional capacity,
  L=1-Kfun/Kpot, and L at least the current vulnerability-curve value.
- Potential-capacity accounting: Kpot=K0+growth-added-retired-potential.
- Prior weight normalisation, independent endpoint-binomial posterior weight
  calculation, analytic baseline-only survival and the alternative PLC mode.
- Bounded and nondecreasing cumulative mortality, including during recovery.
- Malformed six-column forcing rejection.
- Halving the internal step in the original short example: maximum difference
  in prior mean mortality probability 8.17e-5 (0.01 h versus 0.005 h).

The source also passed strict free-format static parser/analysis checks.
Source lines fit the 132-character free-format limit.

## Backward compatibility test completed

Both version 1 and version 2 were compiled and run on the original five-column
example with 500 identical-seed prior particles. The new growth factor defaults
to zero; repair and retirement were off. Maximum differences:

- Stem water potential: 0 MPa.
- PLC: 2.14e-14 percentage points.
- Hydraulic resistance: 1.14e-13 MPa h/L.
- Prior mean mortality probability: 0.

## Two-year checks completed

A continuous 730-day example was run with all three P50 members and 500 prior
particles; repair and turnover were zero. A second run halved the internal
step from 0.10 h to 0.05 h, and a third disabled growth.

- Maximum absolute water balance error in the 0.10 h run: 1.16e-8 L.
- Potential conductance increased from 1.25 to about 1.372615 L h^-1 MPa^-1.
- Actual PLC declined during eligible recovery periods despite persistent
  damage when recovery was disabled.
- Maximum change in prior mean mortality probability after halving the step:
  9.99e-4 (less than 0.001 in absolute probability).
- Maximum stem-potential change in that timestep comparison: 0.00804 MPa.
- The growth-off comparison confirmed fixed potential capacity and less
  functional capacity than the growth-enabled scenario at recovery times.

The root `verify.py` and `tests/check_two_years.py` reproduce the principal
executable checks in temporary folders. Keep an unmodified example copy since
regression assertions use its known parameters. The legacy comparison used the
archived version 1 source and is reported here rather than bundled as a second
model implementation.

## Limits of these results

These are numerical and implementation checks, not physiological validation.
The forcing scenarios, rates, thresholds, priors and mortality counts are
illustrative. More particles/different seeds and finer driver/timestep
resolution are required to establish numerical precision for a particular
application. The two-year timestep result is an example, not a universal error
bound. No real survival data were used.

Testing exposed and corrected an inherited Fortran CSV issue: exponents with
three digits could omit the letter E. All scientific-format CSV outputs now
explicitly reserve three exponent digits, including tiny posterior weights.
