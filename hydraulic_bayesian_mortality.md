Yes. **The code combines established hydraulic and statistical relationships with several modelling choices introduced in our discussion. The complete coupled model is not taken from one published paper.** For your manuscript, distinguish the sources supporting each component from the assumptions you will need to describe and test.

The references below include sources cited earlier and additional references checked now to give you clearer attribution.

| Model component | References to cite | Attribution |
|---|---|---|
| Resistance–capacitance water balance and stem potential | Ruffault et al. (2022); Zhang et al. (2025) | Established framework; our single-compartment implementation is a simplification |
| Soil moisture–water potential conversion | van Genuchten (1980) | Direct mathematical foundation |
| Transpiration from conductance and VPD | Campbell & Norman (1998) | Standard vapour-transfer relationship |
| Multiplicative environmental controls on stomata | Jarvis (1976) | Framework; our particular response functions are empirical choices |
| Water-potential-dependent stomatal closure | Tuzet et al. (2003) | Related physiological formulation; our logistic function is an adaptation |
| Sigmoidal vulnerability curve and \\(P_{50}\\) | Pammenter & Van der Willigen (1998); Duursma & Choat (2017) | Mathematical foundation |
| Conductance loss, residual transpiration and hydraulic failure | Ruffault et al. (2022) | Physiological and modelling foundation |
| Hazard–survival relationship | Klein & Moeschberger (2003); Stan documentation | Standard survival-analysis mathematics |
| Bayesian importance weighting and effective sample size | Owen (2013) | Statistical computational method |

**References and their relevance**

1. **Ruffault, J., Pimont, F., Cochard, H., Dupuy, J.-L., & Martin-StPaul, N. (2022).** SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations of plant water status and drought-induced mortality at the ecosystem level. *Geoscientific Model Development, 15*, 5593–5626. [https://doi.org/10.5194/gmd-15-5593-2022](https://doi.org/10.5194/gmd-15-5593-2022).  
   The strongest single reference for the hydraulic framework: compartment water balances, capacitance, vulnerability-dependent conductance, stomatal closure and residual transpiration. Its equations 15–16 support the logistic PLC curve and conductance proportional to remaining functional xylem. Our model is substantially simpler than SurEau-Ecos. [gmd.copernicus.org](https://gmd.copernicus.org/articles/15/5593/2022/?utm_source=chatgpt.com)

2. **Zhang, Z., Guan, H., Veneklaas, E., Singha, K., & Batelaan, O. (2025).** Revealing seasonal plasticity of whole-plant hydraulic properties using sap-flow and stem water-potential monitoring. *Hydrology and Earth System Sciences, 29*, 3975–3992. [https://doi.org/10.5194/hess-29-3975-2025](https://doi.org/10.5194/hess-29-3975-2025).  
   Supports representing whole-plant water uptake and storage with an RC model and relating effective plant potential to measured stem potential. [hess.copernicus.org](https://hess.copernicus.org/articles/29/3975/2025/?utm_source=chatgpt.com)

3. **van Genuchten, M. Th. (1980).** A closed-form equation for predicting the hydraulic conductivity of unsaturated soils. *Soil Science Society of America Journal, 44*, 892–898. [https://doi.org/10.2136/sssaj1980.03615995004400050002x](https://doi.org/10.2136/sssaj1980.03615995004400050002x).  
   Cite for the soil water-retention relation. The code algebraically inverts this relation to obtain matric potential from volumetric moisture; it does not implement the paper’s complete unsaturated-conductivity model. [Wiley Online Library](https://doi.org/10.2136/sssaj1980.03615995004400050002x?trk=article-ssr-frontend-pulse_little-text-block\\&utm_source=chatgpt.com)

4. **Campbell, G. S., & Norman, J. M. (1998).** *An introduction to environmental biophysics* (2nd ed.). Springer. [https://doi.org/10.1007/978-1-4612-1626-1](https://doi.org/10.1007/978-1-4612-1626-1).  
   A standard reference for water-vapour transport, conductance and the physical basis of \\(E_{\\mathrm{leaf}}\\approx g_wD/P\\). The factor \\(64.8\\) in our code is a unit conversion from mol s-1 to approximately L h-1, rather than an empirical coefficient. [Springer Nature Link](https://link.springer.com/book/10.1007/978-1-4612-1626-1?utm_source=chatgpt.com)

5. **Jarvis, P. G. (1976).** The interpretation of the variations in leaf water potential and stomatal conductance found in canopies in the field. *Philosophical Transactions of the Royal Society of London. B, Biological Sciences, 273*, 593–610. [https://doi.org/10.1098/rstb.1976.0035](https://doi.org/10.1098/rstb.1976.0035).  
   Cite for the empirical framework describing stomatal conductance through environmental response functions. Our Gaussian temperature function and particular VPD function should be identified as chosen parameterisations, not exact reproductions of Jarvis’s equations. [electronicsandbooks.com](https://electronicsandbooks.com/edt/manual/Magazine/P/Philosophical%20Transactions%20of%20the%20Royal%20Society%20B%20UK/1976%20Volume%20273/i927_593-610.pdf?utm_source=chatgpt.com)

6. **Tuzet, A., Perrier, A., & Leuning, R. (2003).** A coupled model of stomatal conductance, photosynthesis and transpiration. *Plant, Cell & Environment, 26*, 1097–1116. [https://doi.org/10.1046/j.1365-3040.2003.01035.x](https://doi.org/10.1046/j.1365-3040.2003.01035.x).  
   Supports coupling stomatal regulation to plant water potential and hydraulic supply. Our code uses stem potential and a simplified logistic closure response; it does not implement the full Tuzet model. [Wiley Online Library](https://onlinelibrary.wiley.com/doi/full/10.1046/j.1365-3040.2003.01035.x?utm_source=chatgpt.com)

7. **Pammenter, N. W., & Van der Willigen, C. (1998).** A mathematical and statistical analysis of the curves illustrating vulnerability of xylem to cavitation. *Tree Physiology, 18*, 589–593. [https://doi.org/10.1093/treephys/18.8-9.589](https://doi.org/10.1093/treephys/18.8-9.589).  
   A primary reference for the exponential sigmoidal vulnerability curve, including separate parameters for its position and slope. This supports our \\(P_{50}\\)–width parameterisation. [Oxford Academic](https://academic.oup.com/treephys/article-abstract/18/8-9/589/1632591?utm_source=chatgpt.com)

8. **Duursma, R. A., & Choat, B. (2017).** fitplc—an R package to fit hydraulic vulnerability curves. *Journal of Plant Hydraulics, 4*, e002. [https://doi.org/10.20870/jph.2017.e002](https://doi.org/10.20870/jph.2017.e002).  
   Supports fitting sigmoidal vulnerability curves and estimating \\(P_{50}\\), curve slope and uncertainty. Useful when replacing the illustrative vulnerability parameters with measured estimates. [Journal of Plant Hydraulics](https://jplanthydro.org/article/view/1541/0?utm_source=chatgpt.com)

9. **Klein, J. P., & Moeschberger, M. L. (2003).** *Survival analysis: Techniques for censored and truncated data* (2nd ed.). Springer. [https://doi.org/10.1007/b97377](https://doi.org/10.1007/b97377).  
   A conventional paper reference for hazard, cumulative hazard, survival and censoring. It supports the statistical framework, rather than our particular hydraulic exposure function. [Springer Nature Link](https://link.springer.com/book/10.1007/b97377?utm_source=chatgpt.com)

10. **Owen, A. B. (2013).** *Monte Carlo theory, methods and examples*. Chapter 9: Importance sampling. [Author’s chapter](https://artowen.su.domains/mc/Ch-var-is.pdf).  
    Supports self-normalised importance sampling and weight-based effective sample size, as used in the Fortran Bayesian calculation. [artowen.su.domains](https://artowen.su.domains/mc/Ch-var-is.pdf?utm_source=chatgpt.com)

11. **Stan Development Team.** *Stan User’s Guide: Survival models*. [Online documentation](https://mc-stan.org/docs/stan-users-guide/survival.html).  
    This was the statistical source explicitly consulted when constructing the linked package. It explains \\(S(t)=\\exp[-H(t)]\\), survival likelihoods and examples of Bayesian priors. The Fortran code does not use Stan. [mc-stan.org](https://mc-stan.org/docs/stan-users-guide/survival.html?utm_source=chatgpt.com)

Two additional papers appeared in the earlier replies, but I would treat them as **supporting context rather than essential derivation references**:

12. **Gao, Q., Zhao, P., Zeng, X., Cai, X., & Shen, W. (2002).** A model of stomatal conductance to quantify the relationship between leaf transpiration, microclimate and soil water stress. *Plant, Cell & Environment, 25*, 1373–1381. [https://doi.org/10.1046/j.1365-3040.2002.00926.x](https://doi.org/10.1046/j.1365-3040.2002.00926.x). Our conductance equation is not a reproduction of this model. [Wiley Online Library](https://onlinelibrary.wiley.com/doi/full/10.1046/j.1365-3040.2002.00926.x?utm_source=chatgpt.com)

13. **Emilio, T., et al. (2019).** Embolism resistance in petioles and leaflets of palms. *Annals of Botany, 124*, 1173–1183. [https://doi.org/10.1093/aob/mcz104](https://doi.org/10.1093/aob/mcz104). This was an example of vulnerability research; Pammenter & Van der Willigen and Ruffault et al. provide more direct references for our equations. [Oxford Academic](https://academic.oup.com/aob/issue/124/7?utm_source=chatgpt.com)

**The parts you should present as assumptions of your model**

The exact hydraulic mortality linkage we introduced,

\\[
x(t)=\\max\\!\\left[0,\\frac{P_{50}-\\Psi_{\\mathrm{stem}}(t)}{\\Psi_{\\mathrm{scale}}}\\right],
\\qquad
\\lambda(t)=\\lambda_0+kx(t)^q,
\\]

is a **proposed exposure–hazard relationship**, not an equation established by the references above. The same applies to the alternative PLC exposure threshold, the no-repair rule \\(L_t=\\max[L_{t-1},L_{\\mathrm{curve}}(\\Psi_t)]\\), and the illustrative priors.

In the paper, describe these explicitly as hypotheses or simplifying assumptions. Cite survival-analysis references for converting hazard to mortality probability, and hydraulic references for the physiological predictors; neither establishes that crossing \\(P_{50}\\) universally initiates mortality.
