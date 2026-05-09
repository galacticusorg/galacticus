# Inline Citation Report

## Summary
- Total inline citations found: 328
- Citations with confirmed BibTeX key: 282
- Citations with NO BibTeX match (need to be added): 46

## Citations to fix (grouped by file)

### source/XRay_Absorption_ISM_Wilms2000.F90
- Line 74: "Wilms et al. (2000)" -> `\cite{wilms_absorption_2000}`
  - Context [inline-comment]: `! Specify values of the parameters to be passed to the dotbvabs function. The default values from Wilms et al. (2000) are used.`

### source/accretion.halo.Naoz_Barkana_2007.F90
- Line 114: "Gnedin (2000)" -> `\cite{gnedin_effect_2000}`
  - Context [inline-comment]: `! Virial density contrast definition used by Gnedin (2000) to define halos, and therefore used in the filtering mass fitting functions.`
- Line 284: "Naoz & Barkana (2007" -> `\cite{naoz_formation_2007}`
  - Context [inline-comment]: `! Evaluate the filtering mass suppression fitting formula as defined by Naoz & Barkana (2007;`
- Line 320: "Naoz & Barkana (2007" -> `\cite{naoz_formation_2007}`
  - Context [inline-comment]: `! Evaluate the rate of change of the filtering mass suppression fitting formula as defined by Naoz & Barkana (2007;`

### source/accretion.halo.cold_mode.F90
- Line 538: "Birnboim & Dekel (2003)" -> `\cite{birnboim_virial_2003}`
  - Context [inline-comment]: `! Compute the shock stability parameter from Birnboim & Dekel (2003).`
- Line 549: "Benson & Bower (2011)" -> `\cite{benson_galaxy_2010-1}`
  - Context [inline-comment]: `! Compute the cold fraction using the model from eqn. (2) of Benson & Bower (2011). The original form does not allow the`

### source/accretion_disks.Shakura_Sunyaev.F90
- Line 124: "Meier (2001)" -> `\cite{meier_association_2001}`
  - Context [inline-comment]: `! Get the black hole spin and dimensionless accretion rate and mass as defined by Meier (2001).`

### source/accretion_disks.spectra.Hopkins2007.F90
- Line 245: "Hopkins et al. (2007)" -> `\cite{hopkins_observational_2007}`
  - Context [string-doc]: `call file%writeAttribute("Hopkins et al. (2007)"                                                                              ,"reference"   )`

### source/black_holes.binaries.separation_growth_rate.standard.F90
- Line 278: "Volonteri et al. (2003)" -> `\cite{volonteri_assembly_2003}`
  - Context [inline-comment]: `! If it does change, we first compute the fraction of that change according to Volonteri et al. (2003) else we set it as the`

### source/black_holes.seeds.Vergara2023.F90
- Line 236: "Binney & Tremaine (2008" -> `\cite{binney_galactic_2008}`
  - Context [inline-comment]: `! Safronov number defined by Binney & Tremaine (2008, https://ui.adsabs.harvard.edu/abs/2008gady.book.....B/abstract)`

### source/chemical.reaction_rates.hydrogen.F90
- Line 255: "Tegmark et al. (1997)" -> `\cite{tegmark_small_1997}`
  - Context [inline-comment]: `! Compute rates. References after each call refer to the rate coefficient in Tegmark et al. (1997) and the equation number in`
- Line 256: "Abel et al. (1997)" -> `\cite{abel_modeling_1997}`
  - Context [inline-comment]: `! Abel et al. (1997) respectively.`
- Line 1587: "Safranek-Shrader et al. (2012" -> `\cite{safranek-shrader_star_2012}`
  - Context [inline-comment]: `! Apply the self-shielding model from Safranek-Shrader et al. (2012;`

### source/cooling.cooling_function.molecular_hydrogen_Galli_Palla.F90
- Line 546: "Galli & Palla (1998)" -> `\cite{galli_chemistry_1998}`
  - Context [inline-comment]: `! The expression from Galli & Palla (1998), assumes an equilibrium (1:3) ratio of para:ortho.`
- Line 559: "Hollenbach & McKee (1979" -> `\cite{hollenbach_molecule_1979}`
  - Context [inline-comment]: `! Rotational and vibrational cooling functions from Hollenbach & McKee (1979; their equations 6.37 and 6.38).`

### source/cooling.time_available.Benson_Bower_2010.F90
- Line 245: "Benson & Bower 2010" -> `\cite{benson_galaxy_2010-1}`
  - Context [inline-comment]: `! Compute the time available for cooling (Benson & Bower 2010; eqn. 16).`

### source/dark_matter_halos.mass_accretion_history.Wechsler2002.F90
- Line 290: "Wechsler et al. (2002" -> `\cite{wechsler_concentrations_2002}`
  - Context [inline-comment]: `double precision                                                , parameter     :: haloMassFraction   =0.015d0 ! Wechsler et al. (2002;  Astrophysical Journal, 568:52-70).`

### source/dark_matter_halos.mass_accretion_history.Zhao2009.F90
- Line 165: "Zhao et al. (2009)" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `sObserved=sigmaObserved*10.0d0**dSigmadMassLogarithmicObserved ! Equation 8 from Zhao et al. (2009).`
- Line 169: "Zhao et al. (2009)" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `wObserved=deltaCriticalObserved/sObserved ! Equation 7 from Zhao et al. (2009).`
- Line 171: "Zhao et al. (2009)" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `pObserved=0.5d0*wObserved/(1.0d0+(0.25d0*wObserved)**6) ! Equation 11 from Zhao et al. (2009).`
- Line 207: "Zhao et al. (2009)" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `&              /                                                                  self%linearGrowth_      %value                               ( time=nowTime(1))   ! Equation 8 from Zhao et al. (2009)`
- Line 217: "Zhao et al. (2009)" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `wNow=deltaCriticalNow/sNow      ! Equation 7 from Zhao et al. (2009).`
- Line 219: "Zhao et al. (2009)" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `pNow=pObserved*max(0.0d0,1.0d0-log10(deltaCriticalNow/deltaCriticalObserved)*wObserved/0.272d0) ! Equation 10 from Zhao et al. (2009).`
- Line 221: "Zhao et al. (2009)" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `dSigmadDeltaCriticalLogarithmic=(wNow-pNow)/5.85d0 ! Equation 12 from Zhao et al. (2009).`

### source/dark_matter_halos.spins.F90
- Line 92: "Bullock et al. (2001" -> `\cite{bullock_profiles_2001}`
  - Context [inline-comment]: `! Use the halo angular momentum scale used in the Bullock et al. (2001; http://adsabs.harvard.edu/abs/2001ApJ...555..240B)`

### source/dark_matter_profiles.structure.concentration.Brown2021.F90
- Line 196: "Brown et al. (2021)" -> `\cite{brown_towards_2022}`
  - Context [inline-comment]: `! Evaluate the concentration using equation (20) of Brown et al. (2021).`

### source/dark_matter_profiles.structure.concentration.Correa2015.F90
- Line 239: "Correa et al. (2015)" -> `\cite{correa_accretion_2015}`
  - Context [inline-comment]: `! Evaluate left-hand side of eqn. 18 of Correa et al. (2015).`
- Line 249: "Correa et al. (2015)" -> `\cite{correa_accretion_2015}`
  - Context [inline-comment]: `! Evaluate eqn. 17 of Correa et al. (2015), rearranged to give formation redshift.`
- Line 265: "Correa et al. (2015)" -> `\cite{correa_accretion_2015}`
  - Context [inline-comment]: `! Evaluate right-hand side of eqn. 19 of Correa et al. (2015).`

### source/dark_matter_profiles.structure.concentration.NFW.F90
- Line 261: "Navarro et al. (1996)" -> `\cite{navarro_structure_1996}`
  - Context [inline-comment]: `! definition used by Navarro et al. (1996).`

### source/dark_matter_profiles.structure.concentration.Schneider2015.F90
- Line 244: "Schneider et al. (2015)" -> `\cite{schneider_structure_2015}`
  - Context [inline-comment]: `! definition used by Schneider et al. (2015).)`
- Line 316: "Schneider et al. (2015)" -> `\cite{schneider_structure_2015}`
  - Context [inline-comment]: `! definition used by Schneider et al. (2015).)`

### source/dark_matter_profiles.structure.concentration.WDM.Bose2016.F90
- Line 139: "Bose et al. (2016)" -> `\cite{bose_copernicus_2016}`
  - Context [inline-comment]: `! Parameters of Bose et al. (2016)'s fitting formula.`

### source/dark_matter_profiles.structure.concentration.WDM.F90
- Line 132: "Schneider et al. (2012)" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! Parameters of Schneider et al. (2012)'s fitting formula.`

### source/dark_matter_profiles.structure.concentration.Zhao2009.F90
- Line 201: "Zhao et al. (2009)" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `! Compute the concentration from the formation time using the Zhao et al. (2009) fitting formula.`

### source/dark_matter_profiles.structure.scale.Ludlow2016.F90
- Line 173: "Ludlow et al. (2016)" -> `\cite{ludlow_mass-concentration-redshift_2016}`
  - Context [inline-comment]: `! Find the density contrast as used to define masses by Ludlow et al. (2016).`

### source/dark_matter_profiles.structure.scale.Ludlow2016.analytic.F90
- Line 142: "Ludlow et al. (2016)" -> `\cite{ludlow_mass-concentration-redshift_2016}`
  - Context [inline-comment]: `! Find the density contrast as used to define masses by Ludlow et al. (2016).`

### source/dark_matter_profiles_DMO.SIDM.coreNFW.F90
- Line 81: "Jiang et al. (2022)" -> `\cite{jiang_semi-analytic_2023}`
  - Context [xml:defaultSource]: `<defaultSource>Jiang et al. (2022)</defaultSource>`

### source/dark_matter_profiles_DMO.accretion_flow.DiemerKravtsov2014.F90
- Line 225: "Diemer & Kravtsov (2014)" -> `\cite{diemer_dependence_2014}`
  - Context [inline-comment]: `! forms which fit the plots in figure 18 of Diemer & Kravtsov (2014). There is no guarantee that these fits will perform`
- Line 253: "Diemer & Kravtsov (2014" -> `\cite{diemer_dependence_2014}`
  - Context [inline-comment]: `! Compute the transition radius following Diemer & Kravtsov (2014; equation 6).`

### source/dark_matter_profiles_DMO.accretion_flow.Shi2016.F90
- Line 196: "Diemer & Kravtsov (2014" -> `\cite{diemer_dependence_2014}`
  - Context [inline-comment]: `! Compute the transition radius following Diemer & Kravtsov (2014; equation 6).`

### source/dark_matter_profiles_DMO.accretion_flow.correlation_function.F90
- Line 205: "Diemer & Kravtsov (2014" -> `\cite{diemer_dependence_2014}`
  - Context [inline-comment]: `! Compute the transition radius following Diemer & Kravtsov (2014; equation 6).`

### source/dark_matter_profiles_DMO.soliton_NFW.F90
- Line 450: "Chan et al. (2022" -> `\cite{chan_diversity_2022}`
  - Context [inline-comment]: `radiusCore         =+5.5d6                           & ! Equation (14) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).`
- Line 455: "Schive et al. (2014" -> `\cite{schive_understanding_2014}`
  - Context [inline-comment]: `densityCore       =+massCore                         & ! Equation (3) of Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).`

### source/dark_matter_profiles_DMO.soliton_NFW_heated.F90
- Line 635: "Chan et al. (2022" -> `\cite{chan_diversity_2022}`
  - Context [inline-comment]: `radiusCore         =+5.5d6                           & ! Equation (14) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).`
- Line 640: "Schive et al. (2014" -> `\cite{schive_understanding_2014}`
  - Context [inline-comment]: `densityCore       =+massCore                         & ! Equation (3) of Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).`

### source/geometry.surveys.Martin-2010-ALFALFA.F90
- Line 154: "Martin et al. (2010)" -> `\cite{martin_arecibo_2010}`
  - Context [inline-comment]: `! The signal-to-noise limit used by Martin et al. (2010).`
- Line 180: "Martin et al. (2010)" -> `\cite{martin_arecibo_2010}`
  - Context [inline-comment]: `! Compute the limiting integrated flux using equation (A1) of Martin et al. (2010).`
- Line 203: "Martin et al. (2010" -> `\cite{martin_arecibo_2010}`
  - Context [inline-comment]: `double precision                                 , parameter               :: solidAngleSurvey=0.79415674617213461d0 ! Computed from survey bounds in Martin et al. (2010; ApJ; 723; 1359)`
- Line 222: "Haynes et al. (2011" -> `\cite{haynes_arecibo_2011}`
  - Context [inline-comment]: `! Survey geometry from Haynes et al. (2011; http://adsabs.harvard.edu/abs/2011AJ....142..170H).`

### source/geometry.surveys.Moustakas-2013-PRIMUS.F90
- Line 238: "Moustakas et al. (2013" -> `\cite{moustakas_primus:_2013}`
  - Context [inline-comment]: `! Find the limiting redshift for this mass completeness limits from Moustakas et al. (2013; Table 2). (See`

### source/halo_model.power_spectrum_modifier.triaxiality.F90
- Line 58: "Smith et al. (2005)" -> `\cite{smith_triaxial_2005}`
  - Context [inline-comment]: `! Tabulated results read from figures in Smith et al. (2005).`

### source/intergalactic_medium.filtering_mass.Gnedin2000.F90
- Line 476: "Naoz & Barkana (2007" -> `\cite{naoz_formation_2007}`
  - Context [doc-comment]: `!! These represent the 2nd order ODE given in equation 11 of Naoz & Barkana (2007, MNRAS, 377, 667;`

### source/kinematic_distributions.Lam2013.F90
- Line 170: "Lam et al. (2013)" -> `\cite{lam_modeling_2013}`
  - Context [inline-comment]: `! Evaluate the mass in the shell outside the halo virial radius using equation (B4) of Lam et al. (2013).`
- Line 183: "Lam et al. (2013)" -> `\cite{lam_modeling_2013}`
  - Context [inline-comment]: `! Compute the nonlinear density contrast using equation (B1) of Lam et al. (2013).`
- Line 195: "Lam et al. (2013)" -> `\cite{lam_modeling_2013}`
  - Context [inline-comment]: `! Evaluate the inflow velocity in the spherical collapse model using equation (B2) of Lam et al. (2013).`

### source/mass_distributions.spherical.Hernquist.F90
- Line 319: "Hernquist (1990)" -> `\cite{hernquist_analytical_1990}`
  - Context [inline-comment]: `! Use the analytic result - equation (37) of Hernquist (1990).`

### source/mass_distributions.spherical.PatejLoeb2015.F90
- Line 244: "Patej & Loeb (2015)" -> `\cite{patej_simple_2015}`
  - Context [inline-comment]: `! Evaluate density using equation 12 of Patej & Loeb (2015).`
- Line 289: "Patej & Loeb 2015)" -> `\cite{patej_simple_2015}`
  - Context [inline-comment]: `! Compute the enclosed mass (eqn. 4 of Patej & Loeb 2015).`

### source/mass_distributions.spherical.SIDM.coreNFW.F90
- Line 28: "Jiang et al. (2022)" -> `\cite{jiang_semi-analytic_2023}`
  - Context [xml:description]: `on the model of Jiang et al. (2022). The profile is defined by the enclosed mass, with \citep{jiang_semi-analytic_2023}:`

### source/mass_distributions.spherical.accretion_flow.Shi2016.F90
- Line 33: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `!        * Scaled    - these correspond to the scaled, self-similar variables used in Appendix A of Shi (2016) - i.e. column 3 of Table A1, "y" for radius, etc.`
- Line 34: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `!        * Original  - these correspond to the original variables used in Appendix A of Shi (2016) - i.e. column 2 of Table A1, "R" for radius, etc.`
- Line 332: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `self%timeNowScaled             =+sqrt(1.0d0-self %cosmologyFunctions_%OmegaMatterEpochal(self%time              )) & ! Equation A5 from Shi (2016).`
- Line 340: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `! Choose a value of the scaled overdensity, "β" in the notation of Shi (2016), identifying a mass shell. We avoid β=1 because such a shell never collapses.`
- Line 356: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `self%radiusTurnaroundScaled(i)=+2.0d0**(2.0d0/3.0d0)                                            & ! Equation A9 from Shi (2016).`
- Line 359: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `bigA                          =+1.0d0                                                           & ! Text after equation of A12 from Shi (2016).`
- Line 361: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `bigB                          =+2.0d0                                                           & ! Text after equation of A12 from Shi (2016).`
- Line 363: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `bigC                          =+                sqrt(1.0d0+4.0d0*bigA)                          & ! Text after equation of A12 from Shi (2016).`
- Line 365: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `self%timeTurnaroundScaled  (i)=+    (    +1.0d0+sqrt(1.0d0+4.0d0*bigA)      )                   & ! Equation A12 of Shi (2016).`
- Line 390: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `! with what Shi (2016) assumes, and has the nice feature that it ensures the mass-epoch relation is the simple scale-free`
- Line 413: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `! not clear if this is precisely what Shi (2016) chose, but it should not matter.`
- Line 439: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `radiusGrowthRateInitialScaled=+sqrt(                                                                                      & ! Equation A8 from Shi (2016).`
- Line 445: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [inline-comment]: `timeInitialScaled            =+2.0d0                                                                                      & ! Equation A6 from Shi (2016).`
- Line 739: "Shi (2016)" -> `\cite{shi_outer_2016}`
  - Context [doc-comment]: `!! Velocity rate of change is given by equation (A7) of Shi (2016).`

### source/mass_distributions.spherical.soliton.Schive2014.F90
- Line 32: "Schive et al. (2014" -> `\cite{schive_understanding_2014}`
  - Context [inline-comment]: `double precision, parameter, public :: coefficientCore=0.091d0 ! Schive et al. (2014; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S; equation 3).`

### source/merger_trees.branching_probability.Fakhouri2010.F90
- Line 220: "Boylan-Kolchin (2010)" -> `\cite{boylan-kolchin_theres_2010}`
  - Context [inline-comment]: `! Evaluate the fitting function of Fakhouri, Ma, & Boylan-Kolchin (2010): d²N/dξ/dz.`
- Line 239: "Boylan-Kolchin (2010)" -> `\cite{boylan-kolchin_theres_2010}`
  - Context [inline-comment]: `! Convert the Fakhouri, Ma, & Boylan-Kolchin (2010) merger rate, d²N/dξ/dz, to the form we want, d²N/dB/dm, where B=δ/σ is the barrier.`

### source/merger_trees.construct.builder.Cole2000.F90
- Line 767: "Cole et al. (2000)" -> `\cite{cole_hierarchical_2000}`
  - Context [inline-comment]: `! In this case we're using the original Cole et al. (2000) algorithm.`

### source/nBody.operator.identify_flybys.MansfieldKravtsov2020.F90
- Line 226: "Mansfield & Kravtsov (2020" -> `\cite{mansfield_three_2020}`
  - Context [inline-comment]: `! Check conditions from Mansfield & Kravtsov (2020; §2.3.1).`

### source/nodes.operators.physics.black_holes.triple_interaction.F90
- Line 129: "Hoffman and Loeb (2007" -> `\cite{hoffman_dynamics_2007}`
  - Context [inline-comment]: `! three body interaction occurs using the radial condition derived in Hoffman and Loeb (2007;`

### source/nodes.operators.physics.dark_matter_profile.prompt_cusp.F90
- Line 453: "Delos (2025)" -> `\cite{delos_cusp-halo_2025}`
  - Context [inline-comment]: `! Compute a new scale radius to obtain the target r₋₂, using equations (23) from Delos (2025).`

### source/nodes.operators.physics.dark_matter_profile.soliton.F90
- Line 70: "Chan et al. (2022" -> `\cite{chan_diversity_2022}`
  - Context [inline-comment]: `! Best-fitting parameters from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).`
- Line 137: "Schive et al. (2014" -> `\cite{schive_understanding_2014}`
  - Context [inline-comment]: `! Using Mc,min ~ 0.25 × Mmin,halo (soliton-dominated limit), from Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).`
- Line 210: "Chan et al. (2022" -> `\cite{chan_diversity_2022}`
  - Context [inline-comment]: `! first term of Equation 15 from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).`
- Line 215: "Chan et al. (2022" -> `\cite{chan_diversity_2022}`
  - Context [inline-comment]: `! second term of Equation 15 from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).`
- Line 269: "Chan et al. (2022" -> `\cite{chan_diversity_2022}`
  - Context [inline-comment]: `! Evaluate the core mass using equation (15) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).`

### source/nodes.operators.physics.dark_matter_profile_axis_ratios.MenkerBenson2022.F90
- Line 225: "Menker & Benson 2022" -> `\cite{menker_random_2022}`
  - Context [inline-comment]: `! Sum energy tensors (Menker & Benson 2022, equation 19).`
- Line 246: "Menker & Benson (2022" -> `\cite{menker_random_2022}`
  - Context [doc-comment]: `!! Menker & Benson (2022; equation 20).`
- Line 255: "Menker & Benson (2022" -> `\cite{menker_random_2022}`
  - Context [doc-comment]: `!! Menker & Benson (2022; equation 23).`
- Line 393: "Menker & Benson (2022)" -> `\cite{menker_random_2022}`
  - Context [doc-comment]: `!! NOTE: This is specific to an NFW profile as in Menker & Benson (2022). Ideally this should be generalized to an arbitrary`

### source/nodes.operators.physics.position.interpolated.F90
- Line 964: "Merson et al. (2013)" -> `\cite{merson_lightcone_2013}`
  - Context [doc-comment]: `!! it and the relative position vector at the final time. This is how we enforce the choice of Merson et al. (2013) to`

### source/nodes.operators.physics.tidal_mass_loss.soliton.F90
- Line 105: "Du et al. (2018" -> `\cite{du_tidal_2018}`
  - Context [inline-comment]: `! Parameters appearing in the fitting function for soliton energy from equation (7) of Du et al. (2018; PRD; 97; 3507;`
- Line 164: "Du et al. (2018" -> `\cite{du_tidal_2018}`
  - Context [inline-comment]: `! radius. Note that equation (7) of Du et al. (2018; PRD; 97; 3507; https://ui.adsabs.harvard.edu/abs/2018PhRvD..97f3507D) is`
- Line 173: "Du et al. (2018" -> `\cite{du_tidal_2018}`
  - Context [inline-comment]: `! Compute the imaginary part of the energy eigenvalue E, using the fitting formula, equation (7) of Du et al. (2018; PRD; 97;`
- Line 180: "Du et al. (2018)" -> `\cite{du_tidal_2018}`
  - Context [inline-comment]: `! Set the soliton core mass evolution rate following Du et al. (2018). From Eq. (17), the core mass obeys (1/M_c) dM_c/dt =`

### source/nodes.property_extractor.Vergara_2023.F90
- Line 189: "Vergara et al. (2024)" -> `\cite{vergara_global_2023}`
  - Context [string-doc]: `descriptions(1)=var_str('Redshift at the formation of the black hole seed in the Vergara et al. (2024) model.'                         )`
- Line 190: "Vergara et al. (2024)" -> `\cite{vergara_global_2023}`
  - Context [string-doc]: `descriptions(2)=var_str('Black hole seed mass in the Vergara et al. (2024) model [M⊙].'                                                )`

### source/nodes.property_extractor.halo_collapse_epoch.F90
- Line 160: "Schneider et al. (2015)" -> `\cite{schneider_structure_2015}`
  - Context [inline-comment]: `! definition used by Schneider et al. (2015).)`

### source/nodes.property_extractor.luminosity_stellar_from_SED.F90
- Line 206: "Hogg et al. (2002)" -> `\cite{hogg_k_2002}`
  - Context [inline-comment]: `! Hogg et al. (2002) and assume that the filter response gives the fraction of incident photons received by the detector at`

### source/nodes.property_extractor.mass_accretion_history.Hearin2021.F90
- Line 150: "Hearin et al. (2021)" -> `\cite{hearin_differentiable_2021}`
  - Context [string-doc]: `descriptions(1)='Early-time power-law index in the Hearin et al. (2021) mass accretion history for this tree.'`
- Line 151: "Hearin et al. (2021)" -> `\cite{hearin_differentiable_2021}`
  - Context [string-doc]: `descriptions(2)='Late-time power-law index in the Hearin et al. (2021) mass accretion history for this tree.'`
- Line 152: "Hearin et al. (2021)" -> `\cite{hearin_differentiable_2021}`
  - Context [string-doc]: `descriptions(3)='Parameter log₁₀t₀ in the Hearin et al. (2021) mass accretion history for this tree.'`

### source/nodes.property_extractor.spin_Bullock.F90
- Line 193: "Bullock et al. (2001)" -> `\cite{bullock_profiles_2001}`
  - Context [string-doc]: `descriptions    (1)=var_str('Spin parameter of the halo under the Bullock et al. (2001) definition [].'                   )`
- Line 195: "Bullock et al. (2001)" -> `\cite{bullock_profiles_2001}`
  - Context [string-doc]: `descriptions(2)=var_str('x-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].')`
- Line 196: "Bullock et al. (2001)" -> `\cite{bullock_profiles_2001}`
  - Context [string-doc]: `descriptions(3)=var_str('y-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].')`
- Line 197: "Bullock et al. (2001)" -> `\cite{bullock_profiles_2001}`
  - Context [string-doc]: `descriptions(4)=var_str('z-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].')`

### source/nodes.property_extractor.velocity_dispersion.F90
- Line 385: "Cappellari et al. (2007" -> `\cite{cappellari_sauron_2007}`
  - Context [inline-comment]: `! The "lambdaR" parameter of Cappellari et al. (2007; MNRAS; 379; 418)`

### source/objects.nodes.components.disk.standard.F90
- Line 1123: "Cole et al. (2000)" -> `\cite{cole_hierarchical_2000}`
  - Context [inline-comment]: `! If using the Cole et al. (2000) method for disk radii, adjust the specific angular momentum to account for the`

### source/output.analyses.HI_vs_halo_mass_relation.ALFALFA_Padmanabhan_2017.F90
- Line 340: "Bryan & Norman (1998)" -> `\cite{bryan_statistical_1998}`
  - Context [inline-comment]: `! Create a halo scale object from which to compute virial velocities. Padmanabhan & Refrigier use the Bryan & Norman (1998)`

### source/output.analyses.Sunyaev-Zeldovich_Planck2013.F90
- Line 221: "Li & White (2009)" -> `\cite{li_distribution_2009}`
  - Context [inline-comment]: `! Construct survey geometry. We use the Li & White (2009) SDSS geometry here, since Planck Intermediate Results XI uses the`

### source/output.analyses.black_hole_vs_halo_mass_relation.F90
- Line 68: "Leauthaud et al. (2012)" -> `\cite{leauthaud_new_2011}`
  - Context [xml:description]: `\item \mono{reference}: a reference for the dataset suitable for inclusion in figures, e.g. \mono{Leauthaud et al. (2012)}.`

### source/output.analyses.concentration_distribution.CDM.COCO.F90
- Line 130: "Ludlow et al. (2016)" -> `\cite{ludlow_mass-concentration-redshift_2016}`
  - Context [inline-comment]: `! Create a virial density contrast object matched to the definition used by Ludlow et al. (2016).`

### source/output.analyses.concentration_vs_mass_relation.CDM.Ludlow_2016.F90
- Line 153: "Ludlow et al. (2016)" -> `\cite{ludlow_mass-concentration-redshift_2016}`
  - Context [inline-comment]: `! Construct mass bins matched to those used by Ludlow et al. (2016).`
- Line 192: "Ludlow et al. (2016)" -> `\cite{ludlow_mass-concentration-redshift_2016}`
  - Context [inline-comment]: `! Create a virial density contrast object matched to the definition used by Ludlow et al. (2016).`

### source/output.analyses.correlation_function.Hearin2014_SDSS.F90
- Line 201: "Hearin et al. (2013)" -> `\cite{hearin_dark_2013}`
  - Context [inline-comment]: `! Build the SDSS survey geometry of Hearin et al. (2013) with their imposed redshift limits.`
- Line 281: "Hearin et al. (2013)" -> `\cite{hearin_dark_2013}`
  - Context [string-doc]: `&                                   var_str('Correlation function for the Hearin et al. (2013) SDSS analysis')                                                , &`

### source/output.analyses.galaxy_sizes_SDSS.F90
- Line 175: "Shen et al. (2003)" -> `\cite{shen_size_2003}`
  - Context [inline-comment]: `!  Random error (in dex) on galaxy stellar masses in the Shen et al. (2003) sample. Shen et al. (2003) quote 95% confidence`
- Line 175: "Shen et al. (2003)" -> `\cite{shen_size_2003}`
  - Context [inline-comment]: `!  Random error (in dex) on galaxy stellar masses in the Shen et al. (2003) sample. Shen et al. (2003) quote 95% confidence`
- Line 194: "Shen et al. (2003)" -> `\cite{shen_size_2003}`
  - Context [inline-comment]: `! Construct sizes matched to those used by  Shen et al. (2003). Also read stellar mass and Sersic index ranges.`
- Line 426: "Shen et al. (2003)" -> `\cite{shen_size_2003}`
  - Context [string-doc]: `&                                                     targetLabel     =var_str('Shen et al. (2003)')                       , &`

### source/output.analyses.luminosity_function.Halpha.Gunawardhana2013.SDSS.F90
- Line 239: "Gunawardhana et al. (2013)" -> `\cite{gunawardhana_galaxy_2013}`
  - Context [inline-comment]: `! Build the SDSS survey geometry of Gunawardhana et al. (2013).`
- Line 296: "Gunawardhana et al. (2013)" -> `\cite{gunawardhana_galaxy_2013}`
  - Context [string-doc]: `&                                  var_str('H$\alpha$ luminosity function for the Gunawardhana et al. (2013) SDSS analysis')                               , &`

### source/output.analyses.luminosity_function.Halpha.Sobral2013.HiZELS.F90
- Line 329: "Sobral et al. (2013)" -> `\cite{sobral_large_2013}`
  - Context [string-doc]: `&                                        var_str('H$\alpha$ luminosity function for the Sobral et al. (2013) HiZELS analysis; redshift interval ')//redshiftInterval , &`

### source/output.analyses.luminosity_function.Montero-Dorta2009.SDSS.F90
- Line 244: "Montero-Dorta & Prada (2009)" -> `\cite{montero-dorta_sdss_2009}`
  - Context [inline-comment]: `! Build the SDSS survey geometry of Montero-Dorta & Prada (2009).`
- Line 301: "Montero-Dorta & Prada (2009)" -> `\cite{montero-dorta_sdss_2009}`
  - Context [string-doc]: `&                                               var_str(band//'-band luminosity function for the Montero-Dorta & Prada (2009) SDSS analysis')                                 , &`

### source/output.analyses.luminosity_function.Stefanon-Marchesini2013.F90
- Line 346: "Stefanon & Marchesini (2013)" -> `\cite{stefanon_evolution_2013}`
  - Context [string-doc]: `&                                               var_str(band//'-band luminosity function for the Stefanon & Marchesini (2013) analysis'), &`

### source/output.analyses.mass_function_HI.ALFALFA_Martin2010.F90
- Line 210: "Martin et al. (2010)" -> `\cite{martin_arecibo_2010}`
  - Context [inline-comment]: `! Build the ALFALFA survey geometry of Martin et al. (2010).`
- Line 253: "Martin et al. (2010)" -> `\cite{martin_arecibo_2010}`
  - Context [string-doc]: `&                              var_str('HI mass function for the Martin et al. (2010) ALFALFA analysis')                              , &`

### source/output.analyses.mass_function_stellar.Bernardi_SDSS.F90
- Line 217: "Bernardi et al. (2013)" -> `\cite{bernardi_massive_2013}`
  - Context [inline-comment]: `! Build the survey geometry of Bernardi et al. (2013) with their imposed redshift limits.`
- Line 274: "Bernardi et al. (2013)" -> `\cite{bernardi_massive_2013}`
  - Context [string-doc]: `&                                   var_str('Stellar mass function for the Bernardi et al. (2013) SDSS analysis')                                     , &`

### source/output.analyses.mass_function_stellar.GAMA.F90
- Line 245: "Baldry et al. (2012)" -> `\cite{baldry_galaxy_2012}`
  - Context [inline-comment]: `! Build the survey geometry of Baldry et al. (2012) with their imposed redshift limits.`
- Line 315: "Baldry et al. (2012)" -> `\cite{baldry_galaxy_2012}`
  - Context [string-doc]: `&                                   var_str('Stellar mass function for the Baldry et al. (2012) GAMA analysis')                                   , &`

### source/output.analyses.mass_function_stellar.PRIMUS.F90
- Line 287: "Moustakas et al. (2013)" -> `\cite{moustakas_primus:_2013}`
  - Context [inline-comment]: `! Build the PRIMUS survey geometry of Moustakas et al. (2013) with their imposed redshift limits.`
- Line 344: "Moustakas et al. 2013)" -> `\cite{moustakas_primus:_2013}`
  - Context [string-doc]: `&                                   var_str('$')//redshiftLabelLow//' < z < '//redshiftLabelHigh//'$ stellar mass function from PRIMUS (Moustakas et al. 2013)', &`

### source/output.analyses.mass_function_stellar.SDSS.F90
- Line 243: "Li & White (2009)" -> `\cite{li_distribution_2009}`
  - Context [inline-comment]: `! Build the SDSS survey geometry of Li & White (2009) with their imposed redshift limits.`
- Line 300: "Li \& White (2009)" -> `\cite{li_distribution_2009}`
  - Context [string-doc]: `&                                   var_str('Stellar mass function for the Li \& White (2009) SDSS analysis')                                         , &`

### source/output.analyses.mass_function_stellar.UKIDSS_UDS.F90
- Line 265: "Caputi et al. (2011)" -> `\cite{caputi_stellar_2011}`
  - Context [inline-comment]: `! Build the UKIDSSUDS survey geometry of Caputi et al. (2011) with their imposed redshift limits.`
- Line 322: "Caputi et al. (2011)" -> `\cite{caputi_stellar_2011}`
  - Context [string-doc]: `&                                   var_str('Stellar mass function for the Caputi et al. (2011) UKIDSS UDS analysis')   , &`

### source/output.analyses.mass_function_stellar.ULTRAVISTA.F90
- Line 247: "Muzzin et al. (2013)" -> `\cite{muzzin_evolution_2013}`
  - Context [string-doc]: `description=var_str('Stellar mass function for the Muzzin et al. (2013) ULTRAVISTA analysis, ')`
- Line 285: "Muzzin et al. (2013)" -> `\cite{muzzin_evolution_2013}`
  - Context [inline-comment]: `! Build the ULTRAVISTA survey geometry of Muzzin et al. (2013) with their imposed redshift limits.`

### source/output.analyses.mass_function_stellar.VIPERS.F90
- Line 265: "Davidzon et al. (2013)" -> `\cite{davidzon_vimos_2013}`
  - Context [inline-comment]: `! Build the VIPERS survey geometry of Davidzon et al. (2013) with their imposed redshift limits.`
- Line 322: "Davidzon et al. (2013)" -> `\cite{davidzon_vimos_2013}`
  - Context [string-doc]: `&                                   var_str('Stellar mass function for the Davidzon et al. (2013) VIPERS analysis')     , &`

### source/output.analyses.mass_function_stellar.ZFOURGE.F90
- Line 298: "Tomczak et al. (2014)" -> `\cite{tomczak_galaxy_2014}`
  - Context [inline-comment]: `! Build the ZFOURGE survey geometry of Tomczak et al. (2014) with their imposed redshift limits.`
- Line 355: "Tomczak et al. 2014)" -> `\cite{tomczak_galaxy_2014}`
  - Context [string-doc]: `&                                   var_str('$')//redshiftLabelLow//' < z < '//redshiftLabelHigh//'$ stellar mass function from ZFOURGE (Tomczak et al. 2014)', &`

### source/output.analyses.mass_metallicity_relation.Andrews2013.F90
- Line 406: "Andrews \& Martini (2013)" -> `\cite{andrews_mass-metallicity_2013}`
  - Context [string-doc]: `&                                                     targetLabel     =var_str('Andrews \& Martini (2013)'                   ), &`

### source/output.analyses.mass_size_relation.Shen2003.F90
- Line 180: "Shen et al. (2003)" -> `\cite{shen_size_2003}`
  - Context [inline-comment]: `! Evaluate the parameters of the Shen et al. (2003) distribution function. (Convert from kpc to Mpc.)`
- Line 182: "Shen et al. 2003" -> `\cite{shen_size_2003}`
  - Context [inline-comment]: `! Late-type galaxy (Shen et al. 2003; eq. 18).`
- Line 185: "Shen et al. 2003" -> `\cite{shen_size_2003}`
  - Context [inline-comment]: `! Early-type galaxy (Shen et al. 2003; eq. 17).`
- Line 188: "Shen et al. (2003" -> `\cite{shen_size_2003}`
  - Context [inline-comment]: `scatterLogarithmic   =sigma2+(sigma1-sigma2)/(1.0d0+(mass/massZeroPoint)**2) ! Shen et al. (2003; eq. 19).`
- Line 274: "Shen et al. (2003)" -> `\cite{shen_size_2003}`
  - Context [string-doc]: `call analysisGroup%writeAttribute('The mass-size relation of Shen et al. (2003).','description'  )`

### source/output.analyses.molecular_ratio.Obreschkow2009.F90
- Line 212: "Obreschkow et al. (2009)" -> `\cite{obreschkow_simulation_2009}`
  - Context [inline-comment]: `! Find the H₂/HI ratio using eqn. (15) of Obreschkow et al. (2009).`

### source/output.analyses.morphological_fraction.GAMA_Moffett2016.F90
- Line 366: "Moffett et al. (2016)" -> `\cite{moffett_galaxy_2016}`
  - Context [string-doc]: `&                                                     targetLabel     =var_str('Moffett et al. (2016)'     ),  &`

### source/output.analyses.quiescent_fraction.Wagner2016.F90
- Line 236: "Wagner et al. (2016)" -> `\cite{wagner_evolution_2016}`
  - Context [string-doc]: `description='Quiescent fraction from Wagner et al. (2016) for galaxies with '`
- Line 237: "Wagner et al. (2016)" -> `\cite{wagner_evolution_2016}`
  - Context [doc-comment]: `!! Determine redshift range properties. Host halo mass threshold is judged approximately from Figure 1 of Wagner et al. (2016).`

### source/output.analyses.spin_distribution.Bett2007.F90
- Line 134: "Bett et al. (2007)" -> `\cite{bett_spin_2007}`
  - Context [inline-comment]: `! Create a virial density contrast class to match the friends-of-friends halo definition used by Bett et al. (2007).`

### source/output.analyses.star_formation_rate_function.Robotham2011.F90
- Line 288: "Robotham et al. 2011" -> `\cite{robotham_galaxy_2011}`
  - Context [string-doc]: `&                                         var_str('Star formation rate function for the Robotham et al. 2011'                                                          ), &`

### source/output.analyses.star_forming_main_sequence.Schreiber2015.F90
- Line 199: "Schreiber et al. (2015)" -> `\cite{schreiber_herschel_2015}`
  - Context [string-doc]: `description='Mean sSFR sequence from Schreiber et al. (2015) for galaxies with '`
- Line 201: "Schreiber et al. (2015)" -> `\cite{schreiber_herschel_2015}`
  - Context [doc-comment]: `!! Schreiber et al. (2015).`
- Line 290: "Schreiber et al. (2015)" -> `\cite{schreiber_herschel_2015}`
  - Context [inline-comment]: `! Build the survey geometry. A more elaborate model could be used here accounting for the different fields and depths used by Schreiber et al. (2015).`

### source/output.analyses.star_forming_main_sequence.Wagner2016.F90
- Line 243: "Wagner et al. (2016)" -> `\cite{wagner_evolution_2016}`
  - Context [string-doc]: `description='Mean sSFR sequence from Wagner et al. (2016) for '`
- Line 260: "Wagner et al. (2016)" -> `\cite{wagner_evolution_2016}`
  - Context [doc-comment]: `!! Determine redshift range properties. Host halo mass threshold is judged approximately from Figure 1 of Wagner et al. (2016).`

### source/output.analyses.stellar_vs_halo_mass_relation.COSMOS_Leauthaud2012.F90
- Line 526: "Leauthaud et al. (2012)" -> `\cite{leauthaud_new_2011}`
  - Context [string-doc]: `&                                                       targetLabel     =var_str('Leauthaud et al. (2012)'              ), &`

### source/output.analyses.stellar_vs_halo_mass_relation.F90
- Line 68: "Leauthaud et al. (2012)" -> `\cite{leauthaud_new_2011}`
  - Context [xml:description]: `\item \mono{reference}: a reference for the dataset suitable for inclusion in figures, e.g. \mono{Leauthaud et al. (2012)}.`

### source/output.analyses.tidal_tracks.velocity_maximum.F90
- Line 163: "Penarrubia et al. (2010)" -> `\cite{penarrubia_impact_2010}`
  - Context [inline-comment]: `! Evaluate the target value. Uses the Penarrubia et al. (2010) fitting function.`

### source/satellites.deceleration_SIDM.Kummer2018.F90
- Line 192: "Kummer et al. (2018)" -> `\cite{kummer_effective_2018}`
  - Context [inline-comment]: `! Note that Kummer et al. (2018) suggest computing an effective escape velocity which is the escape velocity averaged over the`
- Line 234: "Kummer et al. (2018)" -> `\cite{kummer_effective_2018}`
  - Context [inline-comment]: `! Kummer et al. (2018).`
- Line 301: "Kummer et al. (2018)" -> `\cite{kummer_effective_2018}`
  - Context [inline-comment]: `thetaCritical        =+acos(                  & ! Equation (2) from Kummer et al. (2018).`

### source/satellites.dynamical_friction.acceleration.Kaur2018.F90
- Line 117: "Kaur & Sridhar (2018)" -> `\cite{kaur_stalling_2018}`
  - Context [inline-comment]: `! of Figure 8 of Kaur & Sridhar (2018). Specifically we tabulate log₁₀(S) vs. r/r*.`

### source/satellites.evaporation_SIDM.Kummer2018.F90
- Line 190: "Kummer et al. (2018)" -> `\cite{kummer_effective_2018}`
  - Context [inline-comment]: `! Note that Kummer et al. (2018) suggest computing an effective escape velocity which is the escape velocity averaged over the`
- Line 295: "Kummer et al. (2018)" -> `\cite{kummer_effective_2018}`
  - Context [inline-comment]: `thetaCritical        =+acos(                  & ! Equation (2) from Kummer et al. (2018).`

### source/satellites.merging.remnant_sizes.Cole2000.F90
- Line 366: "Cole et al. (2000)" -> `\cite{cole_hierarchical_2000}`
  - Context [inline-comment]: `! Apply the Cole et al. (2000) algorithm to compute the size of the new remnant.`

### source/satellites.merging.remnant_sizes.Covington2008.F90
- Line 338: "Covington et al. (2008)" -> `\cite{covington_predicting_2008}`
  - Context [inline-comment]: `! Apply the Covington et al. (2008) algorithm to compute the size of the new remnant.`

### source/satellites.merging.timescale.Lacey-Cole_Tormen.F90
- Line 100: "Cole et al. (2000)" -> `\cite{cole_hierarchical_2000}`
  - Context [inline-comment]: `double precision                                               , parameter     :: orbitalFactorDistributionSigma=+0.26d0                   !   Cole et al. (2000).`
- Line 101: "Cole et al. (2000)" -> `\cite{cole_hierarchical_2000}`
  - Context [inline-comment]: `double precision                                               , parameter     :: orbitalFactorDistributionMean =-0.14d0                   !   Cole et al. (2000).`

### source/satellites.merging.timescale.Poulton2021.F90
- Line 207: "Poulton et al. 2021" -> `\cite{poulton_extracting_2021}`
  - Context [inline-comment]: `! reduced mass correction (Poulton et al. 2021, eq. 13 for R_peri):`
- Line 224: "Poulton et al. (2021)" -> `\cite{poulton_extracting_2021}`
  - Context [inline-comment]: `! Evaluate the Poulton et al. (2021) fitting formula.`

### source/satellites.merging.timescale.Wetzel-White.F90
- Line 125: "Wetzel & White (2010)" -> `\cite{wetzel_what_2010}`
  - Context [inline-comment]: `double precision                                           , parameter     :: timeScaleNormalization=0.2d0        !   C_dyn from Wetzel & White (2010).`
- Line 149: "Wetzel & White (2010)" -> `\cite{wetzel_what_2010}`
  - Context [inline-comment]: `! Compute dynamical friction timescale using eqn. (2) from Wetzel & White (2010).`

### source/satellites.merging.virial_orbits.Benson2005.F90
- Line 185: "Benson (2005)" -> `\cite{benson_orbital_2005}`
  - Context [inline-comment]: `! Find mass, radius, and velocity in the host corresponding to the Benson (2005) virial density contrast definition.`

### source/satellites.merging.virial_orbits.Jiang2014.F90
- Line 489: "Jiang et al. (2014)" -> `\cite{jiang_generating_2014}`
  - Context [inline-comment]: `! Find virial density contrast under Jiang et al. (2014) definition.`
- Line 493: "Jiang et al. (2014)" -> `\cite{jiang_generating_2014}`
  - Context [inline-comment]: `! Find mass, radius, and velocity in the host and satellite corresponding to the Jiang et al. (2014) virial density contrast`
- Line 495: "Jiang et al. (2014)" -> `\cite{jiang_generating_2014}`
  - Context [inline-comment]: `! mass under the definition of Jiang et al. (2014) can sometimes lead to large satellite masses if the satellite is moving to`

### source/satellites.merging.virial_orbits.Li2020.F90
- Line 297: "Li et al. (2020)" -> `\cite{li_orbital_2020}`
  - Context [inline-comment]: `! Find virial density contrast under Li et al. (2020) definition.`
- Line 301: "Li et al. (2020)" -> `\cite{li_orbital_2020}`
  - Context [inline-comment]: `! Find mass, radius, and velocity in the host and satellite corresponding to the Li et al. (2020) virial density contrast`
- Line 426: "Li et al. (2020)" -> `\cite{li_orbital_2020}`
  - Context [inline-comment]: `! Find virial density contrast under Li et al. (2020) definition.`
- Line 430: "Li et al. (2020)" -> `\cite{li_orbital_2020}`
  - Context [inline-comment]: `! Find mass, radius, and velocity in the host and satellite corresponding to the Li et al. (2020) virial density contrast`

### source/satellites.merging.virial_orbits.Wetzel2010.F90
- Line 259: "Wetzel (2010)" -> `\cite{wetzel_what_2010}`
  - Context [inline-comment]: `! Find virial density contrast under Wetzel (2010) definition.`
- Line 263: "Wetzel (2010)" -> `\cite{wetzel_what_2010}`
  - Context [inline-comment]: `! Find mass, radius, and velocity in the host corresponding to the Wetzel (2010) virial density contrast definition.`

### source/satellites.merging.virial_orbits.loss_cone.F90
- Line 1277: "Sheth & Diaferio (2001)" -> `\cite{sheth_linear_2001}`
  - Context [inline-comment]: `double precision, parameter     :: radiusReferenceExternal=10.0d0 ! Reference radius (in units of Mpc/h) in the Sheth & Diaferio (2001) model for the environmental dependence of velocity dispersion.`

### source/star_formation.rate.nuclear_star_clusters.Krumholz2009.F90
- Line 167: "Krumholz et al. (2009" -> `\cite{krumholz_star_2009}`
  - Context [inline-comment]: `! Computation of the timescale (in units of Gyr⁻¹) given by Krumholz et al. (2009;`

### source/star_formation.rate_surface_density.disks.Blitz2006.F90
- Line 342: "Blitz & Rosolowsky (2006)" -> `\cite{blitz_role_2006}`
  - Context [inline-comment]: `! Compute the pressure ratio that Blitz & Rosolowsky (2006) use to compute the molecular fraction.`
- Line 1198: "Blitz & Rosolowsky (2006)" -> `\cite{blitz_role_2006}`
  - Context [inline-comment]: `! Compute the pressure ratio that Blitz & Rosolowsky (2006) use to compute the molecular fraction.`

### source/stellar_astrophysics.supernovae_type_Ia.Nagashima2005.F90
- Line 61: "Nagashima et al. (2005" -> `\cite{nagashima_metal_2005}`
  - Context [inline-comment]: `! Parameters of the distribution of binaries from Nagashima et al. (2005; MNRAS; 358; 1427; eqn. 17).`
- Line 131: "Nagashima et al. 2005)" -> `\cite{nagashima_metal_2005}`
  - Context [inline-comment]: `! maximum mass is the largest single-star mass for which the endpoint is a C-O white dwarf (Nagashima et al. 2005).`
- Line 193: "Nagashima et al. (2005)" -> `\cite{nagashima_metal_2005}`
  - Context [inline-comment]: `! Evaluate the integrand. Nagashima et al. (2005) give this integrand in the 2D space of (m_b,μ). Here, the quantity we`
- Line 196: "Nagashima et al. (2005)" -> `\cite{nagashima_metal_2005}`
  - Context [inline-comment]: `!   1. divide the integrand given by Nagashima et al. (2005) by φ(m₂) (as it will later be multiplied by it), and;`

### source/stellar_populations.broad_band_luminosities.standard.F90
- Line 659: "Hogg et al. (2002)" -> `\cite{hogg_k_2002}`
  - Context [inline-comment]: `! -c/λ² dλ. Note that we follow the convention of Hogg et al. (2002) and assume that the filter response gives the fraction`

### source/stellar_populations.spectra.postprocess.birth_clouds.Lacey2016.F90
- Line 101: "Lacey et al. (2016)" -> `\cite{lacey_unified_2016}`
  - Context [inline-comment]: `! Implement equation (A5) of Lacey et al. (2016).`

### source/structure_formation.gravitational_lensing.Takahashi_2011.F90
- Line 454: "Takahashi et al. (2011)" -> `\cite{takahashi_probability_2011}`
  - Context [string-doc]: `call parametersFile%writeDataset(tableNKappa             ,"NKappa"             ,"Parameter N_kappa from Takahashi et al. (2011)"    )`
- Line 455: "Takahashi et al. (2011)" -> `\cite{takahashi_probability_2011}`
  - Context [string-doc]: `call parametersFile%writeDataset(tableAKappa             ,"AKappa"             ,"Parameter A_kappa from Takahashi et al. (2011)"    )`
- Line 456: "Takahashi et al. (2011)" -> `\cite{takahashi_probability_2011}`
  - Context [string-doc]: `call parametersFile%writeDataset(tableOmegaKappa         ,"omegaKappa"         ,"Parameter omega_kappa from Takahashi et al. (2011)")`

### source/structure_formation.halo_bias.Tinker2005.F90
- Line 123: "Tinker et al. (2005" -> `\cite{tinker_mass--light_2005}`
  - Context [inline-comment]: `! Fitting function from eqn. (B7) of Tinker et al. (2005; ApJ; 631; 41).`

### source/structure_formation.halo_bias.Tinker2010.F90
- Line 150: "Tinker et al. 2010)" -> `\cite{tinker_large_2010}`
  - Context [inline-comment]: `! Compute parameters as a function of halo overdensity (from Table 2 of Tinker et al. 2010)`
- Line 155: "Tinker et al. 2010)" -> `\cite{tinker_large_2010}`
  - Context [inline-comment]: `! Compute halo bias (equation 6 of Tinker et al. 2010).`

### source/structure_formation.halo_environment.normal.F90
- Line 206: "Mo & White (1996" -> `\cite{mo_analytic_1996}`
  - Context [inline-comment]: `! region has not collapsed on any larger scale. The resulting distribution is given by eqn. (9) of Mo & White (1996; MNRAS;`

### source/structure_formation.halo_mass_function.Reed2007.F90
- Line 121: "Reed et al. (2007)" -> `\cite{reed_halo_2007}`
  - Context [inline-comment]: `! Parameter values from Reed et al. (2007), text after equations (11) and (12).`
- Line 139: "Reed et al. (2007" -> `\cite{reed_halo_2007}`
  - Context [doc-comment]: `!! Scaled peak height defined by Reed et al. (2007; test after equation 9).`
- Line 145: "Reed et al. (2007" -> `\cite{reed_halo_2007}`
  - Context [doc-comment]: `!! Reed et al. (2007; equation 12).`
- Line 148: "Reed et al. (2007" -> `\cite{reed_halo_2007}`
  - Context [doc-comment]: `!! Reed et al. (2007; equation 13).`
- Line 155: "Reed et al. (2007" -> `\cite{reed_halo_2007}`
  - Context [doc-comment]: `!! Reed et al. (2007; equation 12). Note that the published version has some typos. Specifically, in the exponential term`

### source/structure_formation.halo_mass_function.Tinker2008.F90
- Line 219: "Tinker et al. (2008" -> `\cite{tinker_towardhalo_2008}`
  - Context [inline-comment]: `! Extrapolate to higher redshift using redshift scalings given by Tinker et al. (2008; eqns. 5-8).`

### source/structure_formation.halo_mass_function.friends_of_friends_bias.F90
- Line 177: "More et al. (2011" -> `\cite{more_overdensity_2011}`
  - Context [inline-comment]: `double precision                         , parameter               :: correctionCoefficient                                 =  0.220000d+0 ! More et al. (2011; eqn. B11).`
- Line 178: "More et al. (2011" -> `\cite{more_overdensity_2011}`
  - Context [inline-comment]: `double precision                         , parameter               :: fofFractionalAccuracyExponent                         =  1.330000d+0 ! More et al. (2011; nu).`
- Line 179: "More et al. (2011" -> `\cite{more_overdensity_2011}`
  - Context [inline-comment]: `double precision                         , parameter               :: numberDensityPercolation                              =  0.652960d+0 ! More et al. (2011; n_c).`
- Line 233: "More et al. 2011)" -> `\cite{more_overdensity_2011}`
  - Context [inline-comment]: `! Compute friends-of-friends fractional accuracy parameter (Lsize from More et al. 2011).`
- Line 287: "More et al. (2011" -> `\cite{more_overdensity_2011}`
  - Context [inline-comment]: `! utilizes the percolation theory-motivated results of More et al. (2011; their equation B11). The minus sign below arises`

### source/structure_formation.power_spectrum.nonlinear.Smith2003.F90
- Line 185: "Smith et al. (2003" -> `\cite{smith_stable_2003}`
  - Context [inline-comment]: `powerSpectrumQuasiLinear=+self%powerSpectrum_%power(wavenumber,time)       & ! Smith et al. (2003; eqn. 48).`
- Line 194: "Smith et al. (2003" -> `\cite{smith_stable_2003}`
  - Context [inline-comment]: `deltaHaloPrimedSquared  =+   self%a        *y **(3.0d0*self%f1)    & ! Smith et al. (2003; eqn. 50).`
- Line 200: "Smith et al. (2003" -> `\cite{smith_stable_2003}`
  - Context [inline-comment]: `deltaHaloSquared        =+deltaHaloPrimedSquared & ! Smith et al. (2003; eqn. 51).`
- Line 260: "Smith et al. 2003" -> `\cite{smith_stable_2003}`
  - Context [inline-comment]: `! Compute effective spectral index (Smith et al. 2003; eqn. 59) and spectral curvature (Smith et al. 2003; eqn. 60).`
- Line 260: "Smith et al. 2003" -> `\cite{smith_stable_2003}`
  - Context [inline-comment]: `! Compute effective spectral index (Smith et al. 2003; eqn. 59) and spectral curvature (Smith et al. 2003; eqn. 60).`
- Line 294: "Smith et al. (2003" -> `\cite{smith_stable_2003}`
  - Context [inline-comment]: `! Evaluate fitting function coefficients - using the fits from Smith et al. (2003; eqns C9 through C16).`

### source/structure_formation.transfer_function.BBKS.WDM.F90
- Line 216: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/structure_formation.transfer_function.Bode2001.F90
- Line 191: "Barkana et al. (2001" -> `\cite{barkana_constraints_2001}`
  - Context [inline-comment]: `! This uses equation (4) from Barkana et al. (2001; http://adsabs.harvard.edu/abs/2001ApJ...558..482B), with the`
- Line 207: "Bode et al. (2001" -> `\cite{bode_halo_2001}`
  - Context [inline-comment]: `! This uses equation (A9) from Bode et al. (2001; https://ui.adsabs.harvard.edu/abs/2001ApJ...556...93B).`
- Line 224: "Viel et al. (2005" -> `\cite{viel_constraining_2005}`
  - Context [inline-comment]: `! This uses equation (7) from Viel et al. (2005; https://ui.adsabs.harvard.edu/abs/2005PhRvD..71f3534V).`
- Line 391: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`
- Line 419: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`
- Line 448: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/structure_formation.transfer_function.ETHOS.F90
- Line 562: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/structure_formation.transfer_function.Hu2000.FDM.F90
- Line 231: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`
- Line 262: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/structure_formation.transfer_function.Murgia2017.F90
- Line 233: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`
- Line 274: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/structure_formation.transfer_function.envelope.F90
- Line 316: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/structure_formation.transfer_function.file.F90
- Line 760: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/structure_formation.transfer_function.file.fuzzy_dark_matter.F90
- Line 329: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`
- Line 376: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/structure_formation.transfer_function.fuzzyDM.Passaglia2022.F90
- Line 136: "Passaglia & Hu (2022" -> `\cite{passaglia_accurate_2022}`
  - Context [inline-comment]: `self%A                 =+2.22d0*     self%m22**(1.0d0/25.0d0-1.0d0/1000.0d0*log(self%m22))  ! Equation 46 of Passaglia & Hu (2022, PRD, 105, 123529; https://ui.adsabs.harvard.edu/abs/2022PhRvD.105l352`
- Line 137: "Passaglia & Hu (2022" -> `\cite{passaglia_accurate_2022}`
  - Context [inline-comment]: `self%B                 =+0.16d0/     self%m22**(1.0d0/20.0d0                             )  ! Equation 46 of Passaglia & Hu (2022, PRD, 105, 123529; https://ui.adsabs.harvard.edu/abs/2022PhRvD.105l352`
- Line 138: "Passaglia & Hu (2022" -> `\cite{passaglia_accurate_2022}`
  - Context [inline-comment]: `self%wavenumberJeans   =+9.00d0*sqrt(self%m22                                             ) ! Equation 45 of Passaglia & Hu (2022, PRD, 105, 123529; https://ui.adsabs.harvard.edu/abs/2022PhRvD.105l352`
- Line 270: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/structure_formation.transfer_function.half_mode_slope.F90
- Line 134: "Murgia et al. (2017)" -> `\cite{murgia_non-cold_2017}`
  - Context [inline-comment]: `! Compute the parameters for the underlying Murgia et al. (2017) transfer function.`

### source/tasks.power_spectrum.F90
- Line 350: "Schneider et al. (2012" -> `\cite{schneider_non-linear_2012}`
  - Context [inline-comment]: `! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].`

### source/tests.Zhao2009_algorithms.EdS.F90
- Line 87: "Zhao et al. 2009" -> `\cite{zhao_accurate_2009}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Zhao et al. 2009 algorithms: Einstein-de Sitter cosmology")`
- Line 88: "Zhao et al. 2009" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `! Test Zhao et al. 2009 algorithms in an Einstein-de Sitter universe.`

### source/tests.Zhao2009_algorithms.dark_energy.F90
- Line 87: "Zhao et al. 2009" -> `\cite{zhao_accurate_2009}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Zhao et al. 2009 algorithms: dark energy cosmology")`
- Line 88: "Zhao et al. 2009" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `! Test Zhao et al. 2009 algorithms in a dark energy universe.`

### source/tests.Zhao2009_algorithms.open.F90
- Line 87: "Zhao et al. 2009" -> `\cite{zhao_accurate_2009}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Zhao et al. 2009 algorithms: open cosmology")`
- Line 88: "Zhao et al. 2009" -> `\cite{zhao_accurate_2009}`
  - Context [inline-comment]: `! Test Zhao et al. 2009 algorithms in an open universe.`

### source/tests.biases.F90
- Line 94: "Sheth et al. (2001)" -> `\cite{sheth_ellipsoidal_2001}`
  - Context [string-doc]: `modelName           (2)='Sheth et al. (2001)'`
- Line 98: "Tinker et al. (2010" -> `\cite{tinker_large_2010}`
  - Context [string-doc]: `modelName           (3)='Tinker et al. (2010; 200c)'`

### source/tests.concentration.Correa2015.F90
- Line 66: "Correa et al. 2015" -> `\cite{correa_accretion_2015}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Correa et al. 2015 concentration-mass relation")`
- Line 67: "Correa et al. 2015" -> `\cite{correa_accretion_2015}`
  - Context [inline-comment]: `! Test Correa et al. 2015 algorithm.`

### source/tests.concentrations.F90
- Line 105: "Diemer & Kravtsov (2015)" -> `\cite{diemer_universal_2014}`
  - Context [string-doc]: `modelName           (1)='Diemer & Kravtsov (2015)'`
- Line 111: "Dutton et al. (2014" -> `\cite{dutton_cold_2014}`
  - Context [string-doc]: `modelName           (2)='Dutton et al. (2014; 200c)'`
- Line 115: "Dutton et al. (2014" -> `\cite{dutton_cold_2014}`
  - Context [string-doc]: `modelName           (3)='Dutton et al. (2014; vir)'`
- Line 119: "Ludlow et al. (2016)" -> `\cite{ludlow_mass-concentration-redshift_2016}`
  - Context [string-doc]: `modelName           (4)='Ludlow et al. (2016)'`
- Line 123: "Prada et al. (2012)" -> `\cite{prada_halo_2011}`
  - Context [string-doc]: `modelName           (5)='Prada et al. (2012)'`
- Line 127: "Diemer & Joyce (2019)" -> `\cite{diemer_accurate_2019}`
  - Context [string-doc]: `modelName           (6)='Diemer & Joyce (2019)'`
- Line 131: "Diemer & Joyce (2019" -> `\cite{diemer_accurate_2019}`
  - Context [string-doc]: `modelName           (7)='Diemer & Joyce (2019; vir [converted])'`
- Line 554: "Ludlow et al. (2016)" -> `\cite{ludlow_mass-concentration-redshift_2016}`
  - Context [inline-comment]: `! For the Ludlow et al. (2016) model, also check the direct (non-fitting function) calculation.`

### source/tests.dark_matter_profiles.F90
- Line 63: "Jiang et al. 2022)" -> `\cite{jiang_semi-analytic_2023}`
  - Context [inline-comment]: `! Mass and concentration of Pippin halos (Jiang et al. 2022).`
- Line 1563: "Jiang et al. (2022)" -> `\cite{jiang_semi-analytic_2023}`
  - Context [inline-comment]: `! Interaction radius estimated from Figure 2 of Jiang et al. (2022).`
- Line 1696: "Jiang et al. (2022)" -> `\cite{jiang_semi-analytic_2023}`
  - Context [doc-comment]: `!! Target values were measured from Figure A1 of Jiang et al. (2022).`

### source/tests.dark_matter_profiles.adiabaticGnedin2004.F90
- Line 167: "Gnedin et al. (2004)" -> `\cite{gnedin_response_2004}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group('Gnedin et al. (2004) dark matter profile')`

### source/tests.gaunt_factors.F90
- Line 48: "Hoof et al. (2014)" -> `\cite{van_hoof_accurate_2014}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("van Hoof et al. (2014) fitting function:")`

### source/tests.halo_mass_function.Reed2007.F90
- Line 76: "Reed et al. (2007)" -> `\cite{reed_halo_2007}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Halo mass function: Reed et al. (2007)")`

### source/tests.halo_mass_function.Tinker.F90
- Line 70: "Tinker et al. (2008)" -> `\cite{tinker_towardhalo_2008}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Halo mass function: Tinker et al. (2008)")`

### source/tests.initial_mass_functions.F90
- Line 167: "Nagashima et al. (2005)" -> `\cite{nagashima_metal_2005}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group('Nagashima et al. (2005)')`
- Line 182: "Nagashima et al. (2005)" -> `\cite{nagashima_metal_2005}`
  - Context [inline-comment]: `! The tolerance for this test is (very) low. Nagashima et al. (2005) do not specify what they use for M(t) - the mass of a`

### source/tests.mass_accretion_history.Correa2015.F90
- Line 67: "Correa et al. 2015" -> `\cite{correa_accretion_2015}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Correa et al. 2015 mass accretion history algorithms")`
- Line 68: "Correa et al. 2015" -> `\cite{correa_accretion_2015}`
  - Context [inline-comment]: `! Test Correa et al. 2015 algorithm.`

### source/tests.mass_accretion_history.Hearin2021.F90
- Line 58: "Correa et al. 2015" -> `\cite{correa_accretion_2015}`
  - Context [inline-comment]: `! Test Correa et al. 2015 algorithm.`

### source/tests.mass_accretion_history.Hearin2021_stochastic.F90
- Line 85: "Correa et al. 2015" -> `\cite{correa_accretion_2015}`
  - Context [inline-comment]: `! Test Correa et al. 2015 algorithm.`

### source/tests.mass_distributions.F90
- Line 53: "Mazure & Capelato (2001)" -> `\cite{mazure_exact_2002}`
  - Context [inline-comment]: `! Mass targets for Sersic profile from Mazure & Capelato (2001).`
- Line 55: "Mazure & Capelato (2001)" -> `\cite{mazure_exact_2002}`
  - Context [inline-comment]: `! Density targets for Sersic profile from Mazure & Capelato (2001).`
- Line 563: "Patej & Loeb (2015)" -> `\cite{patej_simple_2015}`
  - Context [inline-comment]: `! Patej & Loeb (2015) profile.`

### source/tests.spectra.postprocess.Inoue2014.F90
- Line 43: "Inoue et al. (2014)" -> `\cite{inoue_updated_2014}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Inoue et al. (2014) IGM attenuation model")`

### source/tests.stellar_populations.F90
- Line 122: "Chabrier (2001)" -> `\cite{chabrier_galactic_2001}`
  - Context [inline-comment]: `! the interval. This is compared to a previously computed value for the Chabrier (2001) IMF.`

### source/tests.stellar_populations.luminosities.F90
- Line 231: "Bruzual & Charlot (2003)" -> `\cite{bruzual_stellar_2003}`
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Comparison with Bruzual & Charlot (2003) code")`

### source/tests.transfer_functions.F90
- Line 151: "Eisenstein & Hu (1999)" -> `\cite{eisenstein_power_1999}`
  - Context [inline-comment]: `! We expect agreement between Eisenstein & Hu (1999) and CAMB over only a limited range of wavenumbers. Outside of that range force them to be equal to avoid failed assertions.`

## Citations with NO BibTeX match (need new entries)

- `source/XRay_Absorption_ISM_Wilms2000.F90:113`: "Allen & McCray (2000" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `call outputFile%writeAttribute("Wilms, Allen & McCray (2000; ApJ, 542, 914; http://adsabs.harvard.edu/abs/2000ApJ...542..914W)","source")`
- `source/atomic.rates.recombination.radiative.Verner.F90:588`: "Fergusen & Ferland (1997)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Coefficients in the fitting function of Fergusen & Ferland (1997).`
- `source/atomic.rates.recombination.radiative.Verner.F90:676`: "Fergusen & Ferland (1997)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Fergusen & Ferland (1997).`
- `source/atomic.rates.recombination.radiative.Verner.F90:731`: "Storey & Hummer (1995" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Fit from Aparna Venkatesan, via Mike Shull. This is a fit to the data in Table 1 of Storey & Hummer (1995; MNRAS; 272; 41) for Z=2 and Ne=100.`
- `source/black_holes.seeds.Vergara2023.F90:245`: "Landau & Lifshitz (1980" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Probabilistic mean free path defined as in Landau & Lifshitz (1980,`
- `source/galactic_structure.radius_definition.F90:67`: "Ensellem et al. (2007" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `<entry label="lambdaR"                    description="The λᵣ parameter of Ensellem et al. (2007; https://ui.adsabs.harvard.edu/abs/2007MNRAS.379..401E) is computed"/>`
- `source/geometry.surveys.Li-White-2009-SDSS.F90:261`: "Percival et al. (2010" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `double precision                               , parameter               :: solidAngleSurvey=2.1901993d0 ! From Percival et al. (2010; MNRAS; 401; 2148)`
- `source/geometry.surveys.Local_Group_DES.F90:123`: "Drlica-Wagner et al. (2015" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! units in the V band) to the detection efficiencies reported by Drlica-Wagner et al. (2015, ApJ, 813, 109;`
- `source/geometry.surveys.Local_Group_SDSS.F90:101`: "Peter & Hargis (2018" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Find the limiting distance for this mass completeness limits. We adopt the model of Kim, Peter & Hargis (2018,`
- `source/models.likelihoods.halo_mass_function.F90:806`: "Cole et al. (2008" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `fractionMassProgenitors=+0.4d0                                     & ! "Global fit" from Cole et al. (2008; https://ui.adsabs.harvard.edu/abs/2008MNRAS.383..546C).`
- `source/nodes.operators.physics.dark_matter_profile.soliton.F90:138`: "Hui et al. (2017" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! The minimum halo-mass scaling is mentioned in the abstract of Hui et al. (2017; PRD; 95; 3541; https://ui.adsabs.harvard.edu/abs/2017PhRvD..95d3541H).`
- `source/nodes.property_extractor.luminosity_emission_line.Panuzzo2003.F90:320`: "Savage & Mathis 1979)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `double precision                                                 , parameter               :: AVToEBV                       =+3.10d+00                   ! (A_V/E(B-V); Savage & Mathis 1979)`
- `source/nodes.property_extractor.luminosity_emission_line.Panuzzo2003.F90:321`: "Savage & Mathis 1979)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `double precision                                                 , parameter               :: NHToEBV                       =+5.80d+21                   ! (N_H/E(B-V); atoms/cm²/mag; Savage & Mathis 1`
- `source/nodes.property_extractor.luminosity_emission_line.Panuzzo2003.F90:411`: "Tumlinson (2009)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Compute the (logarithm of) hydrogen density, based on the model of Krumholz, McKee, & Tumlinson (2009) for molecular cloud`
- `source/nodes.property_extractor.luminosity_emission_line_AGN.F90:439`: "Gutkin & Charlot (2016" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Calculate the Strömgren radius (equation 3 of Feltre, Gutkin & Charlot (2016; MNRAS; 456; 3354;`
- `source/nodes.property_extractor.luminosity_stellar.dust.CharlotFall2000.F90:211`: "Savage & Mathis 1979)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `double precision                                        , parameter               :: AVToEBV                  =+3.10d+00            ! (A_V/E(B-V); Savage & Mathis 1979)`
- `source/nodes.property_extractor.luminosity_stellar.dust.CharlotFall2000.F90:212`: "Savage & Mathis 1979)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `double precision                                        , parameter               :: NHToEBV                  =+5.80d+21            ! (N_H/E(B-V); atoms/cm²/mag; Savage & Mathis 1979)`
- `source/output.analyses.HI_vs_halo_mass_relation.ALFALFA_Padmanabhan_2017.F90:330`: "Mo & Tormen (2001)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! friends-of-friends algorithm with linking length parameter b=0.2 since that is what was used by Sheth, Mo & Tormen (2001) in`
- `source/output.analyses.HI_vs_halo_mass_relation.ALFALFA_Padmanabhan_2017.F90:331`: "Padmanabhan & Refregier 2017)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! their original calibration of their halo mass function (as used by Padmanabhan & Refregier 2017).`
- `source/output.analyses.HI_vs_halo_mass_relation.ALFALFA_Padmanabhan_2017.F90:367`: "Padmanabhan & Refrigier (2017)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! contrast used by Padmanabhan & Refrigier (2017).`
- `source/output.analyses.HI_vs_halo_mass_relation.ALFALFA_Padmanabhan_2017.F90:413`: "Refrigier (2017)" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `&                                                     targetLabel     =var_str('Padmanabhan \\& Refrigier (2017)'           ), &`
- `source/output.analyses.Local_Group.occupation_fraction.F90:373`: "Nadler et al. (2020)" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `&                                                      targetLabel     =var_str('Nadler et al. (2020)'              ), &`
- `source/output.analyses.Local_Group.occupation_fraction.F90:549`: "Nadler et al. (2020)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Nadler et al. (2020) report that the 50% mass is log₁₀(M₅₀/M☉)=7.51⁺⁰˙²¹₋₀.₀₀ - that is it is unconstrained below`
- `source/output.analyses.Local_Group.stellar_mass_halo_mass_relation.F90:371`: "Nadler et al. (2020)" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `&                                                      targetLabel     =var_str('Nadler et al. (2020)'                                ), &`
- `source/output.analyses.Sunyaev-Zeldovich_Planck2013.F90:326`: "XI (2013)" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `&                                                     targetLabel     =var_str('Planck Intermediate Results XI (2013)'), &`
- `source/output.analyses.black_hole_vs_halo_mass_relation.F90:60`: "Norman (1998)" -- no matching entry in Galacticus.bib
  - Context [xml:description]: `\item \mono{`Bryan \&amp; Norman (1998)'}: halos are defined as have mean density contrasts given by the`
- `source/output.analyses.color_distribution.SDSS.F90:147`: "Baldry et al. (2004)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Construct colors matched to those used by Baldry et al. (2004). Also read magnitude range.`
- `source/output.analyses.color_distribution.SDSS.F90:189`: "Baldry et al. (2004)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Build the SDSS survey geometry of Baldry et al. (2004) with their imposed redshift limits.`
- `source/output.analyses.color_distribution.SDSS.F90:282`: "Baldry et al. (2004)" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `&                                                     targetLabel     =var_str('Baldry et al. (2004)')         , &`
- `source/output.analyses.galaxy_sizes_SDSS.F90:240`: "Shen et al. (2013)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Build the SDSS survey geometry of Shen et al. (2013) with their imposed redshift limits.`
- `source/output.analyses.mass_metallicity_relation.Blanc2019.F90:410`: "Blanc et al. (2017)" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `&                                                     targetLabel     =var_str('Blanc et al. (2017)'                         ), &`
- `source/output.analyses.size_vs_stellar_mass_relation.F90:81`: "Wel et al. (2014)" -- no matching entry in Galacticus.bib
  - Context [xml:description]: `\item \mono{reference}: a reference for the dataset suitable for inclusion in figures, e.g. \mono{van der Wel et al. (2014)}.`
- `source/output.analyses.stellar_vs_halo_mass_relation.COSMOS_Leauthaud2012.F90:456`: "More et al. (2009" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! on the constraint on σ_{log₁₀L}=0.16±0.04 from More et al. (2009; MNRAS; 392; 801) for SDSS galaxies.`
- `source/output.analyses.stellar_vs_halo_mass_relation.COSMOS_Leauthaud2012.F90:479`: "More et al. (2009)" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `&                                                       targetLabel     =var_str('More et al. (2009)'                            ), &`
- `source/output.analyses.stellar_vs_halo_mass_relation.F90:60`: "Norman (1998)" -- no matching entry in Galacticus.bib
  - Context [xml:description]: `\item \mono{`Bryan \&amp; Norman (1998)'}: halos are defined as have mean density contrasts given by the`
- `source/star_formation.rate.nuclear_star_clusters.Krumholz2009.F90:168`: "Sesana et al. (2014" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! https://ui.adsabs.harvard.edu/abs/2009ApJ...699..850K/abstract) and Sesana et al. (2014,`
- `source/structure_formation.power_spectrum.nonlinear.PeacockDodds1996.F90:172`: "Press & Turner 1992)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Compute growth factor using same fitting function as Peacock & Dodds (from Carroll, Press & Turner 1992).`
- `source/structure_formation.spherical_collapse.solver.collisionlessMatter_cosmologicalConstant.F90:449`: "Lahav et al. 1991)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Find the η-factor (see Lahav et al. 1991) which measures the dark energy contribution to the energy of the`
- `source/structure_formation.transfer_function.Bode2001.F90:241`: "Vogel & Azabajian (2023" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! This uses equation (9) from Vogel & Azabajian (2023; https://ui.adsabs.harvard.edu/abs/2023PhRvD.108d3520V) with parameters for spin-1/2 particles.`
- `source/structure_formation.transfer_function.Bode2001.F90:258`: "Vogel & Azabajian (2023" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! This uses equation (9) from Vogel & Azabajian (2023; https://ui.adsabs.harvard.edu/abs/2023PhRvD.108d3520V) with parameters for spin-3/2 particles.`
- `source/structure_formation.transfer_function.Eisenstein_Hu1999.F90:220`: "Komatsu et al. (2007" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Specify properties of neutrinos. Mass fraction formula is from Komatsu et al. (2007; http://adsabs.harvard.edu/abs/2010arXiv1001.4538K).`
- `source/tests.dark_matter_profiles.F90:1548`: "Elbert et al. (2015" -- no matching entry in Galacticus.bib
  - Context [doc-comment]: `!! Set properties to match the Pippin halos of Elbert et al. (2015; https://ui.adsabs.harvard.edu/abs/2015MNRAS.453...29E).`
- `source/tests.dark_matter_profiles.fuzzy_dark_matter.F90:182`: "Chowdhury et al. (2021" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Construct a node with properties matched to that in Chowdhury et al. (2021; https://ui.adsabs.harvard.edu/abs/2021ApJ...916...27D).`
- `source/tests.mass_distributions.F90:57`: "Young (1976)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Potential targets for Sersic profile from Young (1976).`
- `source/tests.mass_distributions.F90:151`: "Young (1976)" -- no matching entry in Galacticus.bib
  - Context [inline-comment]: `! Convert to reduced potential used by Young (1976).`
- `source/tests.mass_distributions.F90:564`: "Patej-Loeb (2015)" -- no matching entry in Galacticus.bib
  - Context [string-doc]: `call Unit_Tests_Begin_Group("Patej-Loeb (2015) profile")`
