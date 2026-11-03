# Version 2.0 — 3 November 2026

- Replace the maximum-loss-only state with potential and functional conductance.
- Add seasonal new functional xylem, optional gated conduit repair and optional
  proportional retirement, while retaining the common P50 vulnerability cap.
- Add optional growth and leaf-area drivers; retain five-column input support.
- Append conductance/capacity accounting fields to hydraulic output.
- Keep continuous cumulative mortality across years and allow instantaneous
  hydraulic hazard to decline as functional capacity recovers.
- Add a two-year example and analytic/regression checks for the new mechanisms.
- Retain fixed capacitance, height and non-xylem resistance; no carbon budget,
  xylem age cohorts, physiological phenology model or restart file is introduced.
- Existing mortality posteriors must be reassessed if recovery assumptions change.

- Use explicit three-digit scientific exponents in CSV outputs so very small
  weights remain valid numeric tokens in R and Python.
