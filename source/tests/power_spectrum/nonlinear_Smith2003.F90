!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Contains a program to test the :cite:t:`smith_stable_2003` nonlinear power spectrum.
  !!}

program Test_Power_Spectrum_Nonlinear_Smith2003
  !!{RST
  Test the :galacticus-class:`powerSpectrumNonlinearSmith2003` nonlinear matter power spectrum - the ``halofit`` algorithm -
  against the independent implementation in CAMB, run with ``halofit_version = 1``, which selects the original form of
  :cite:t:`smith_stable_2003` rather than any of its later revisions.

  The comparison is made on the ratio :math:`P_\mathrm{nonlinear}(k)/P_\mathrm{linear}(k)` rather than on the power spectrum
  itself. That ratio is what ``halofit`` actually supplies, and taking it removes any difference in the overall normalization of
  the linear spectrum between the two codes. The linear spectrum itself is taken from CAMB here, through
  :galacticus-class:`transferFunctionCAMB`, because ``halofit`` is not a function of :math:`k` alone: it depends on the nonlinear
  scale, and on the effective index and curvature of the linear spectrum there, so a comparison made with a different linear
  spectrum would not be testing the same thing.

  The algorithm writes the nonlinear power as the sum of a quasi-linear term, which carries the linear spectrum damped on small
  scales, and a halo term, which supplies the power from collapsed objects. The class can be built with either term alone, and
  the two are checked to sum to the whole, and separately to dominate where each should: the quasi-linear term on large scales,
  where the ratio must approach unity, and the halo term on small scales.

  The cosmology, and the spectral index and :math:`\sigma_8` below, are those of the CAMB run which produced the reference
  values, with :math:`\sigma_8` taken from what CAMB itself reported for that run rather than imposed.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard   , displayMessage
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Power_Spectra                       , only : powerSpectrumStandard
  use :: Power_Spectra_Nonlinear             , only : powerSpectrumNonlinearSmith2003
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionCAMB                    , transferFunctionTypeTotal
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group   , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple               ), pointer      :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda          ), pointer      :: cosmologyFunctions_
  type            (darkMatterParticleCDM                   ), pointer      :: darkMatterParticle_
  type            (linearGrowthCollisionlessMatter         ), pointer      :: linearGrowth_
  type            (powerSpectrumPrimordialPowerLaw         ), pointer      :: powerSpectrumPrimordial_
  type            (transferFunctionCAMB                    ), pointer      :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple), pointer      :: powerSpectrumPrimordialTransferred_
  type            (powerSpectrumWindowFunctionTopHat       ), pointer      :: powerSpectrumWindowFunction_
  type            (cosmologicalMassVarianceFilteredPower   ), pointer      :: cosmologicalMassVariance_
  type            (powerSpectrumStandard                   ), pointer      :: powerSpectrum_
  type            (powerSpectrumNonlinearSmith2003         ), pointer      :: nonlinear_                        , nonlinearQuasiLinear_, &
       &                                                                      nonlinearHalo_
  ! The cosmology of the CAMB run which produced the reference values: Ω_b h² = 0.0226, Ω_c h² = 0.112, h = 0.67, and
  ! a flat universe.
  double precision                                          , parameter    :: HubbleConstant       =67.0000d0   , OmegaBaryon    =0.0503450d0, &
       &                                                                      OmegaMatter          = 0.2998440d0, OmegaDarkEnergy=0.7001560d0, &
       &                                                                      temperatureCMB       = 2.7255d0   , indexPrimordial=0.9600000d0, &
       &                                                                      sigma8               = 0.7722d0
  ! Wavenumbers [Mpc⁻¹] and the ratio of the nonlinear to the linear power spectrum given by CAMB at each.
  integer                                                   , parameter    :: countWavenumbers     =9
  double precision                                          , dimension(9) :: wavenumbers         =[0.00670000d0,0.03350000d0,0.06700000d0,0.20100000d0,0.33500000d0,0.67000000d0,1.34000000d0,3.35000000d0,6.70000000d0]
  double precision                                          , dimension(9) :: ratioReference      =[0.9939648453d0,0.9768489238d0,0.9951035101d0,1.4202527540d0,2.2625702250d0,5.4373947440d0,12.6195313966d0,25.4813838837d0,33.3795600800d0]
  ! The two terms of the algorithm separately, from an independent implementation of Appendix C of Smith et al. (2003) driven
  ! by the same CAMB linear spectrum (`halofitDecompose.py` in galacticusDevTools). That implementation reproduces the total
  ! CAMB reports above to better than 0.4%, which is what licenses its use for the decomposition. For this linear spectrum it
  ! finds a nonlinear scale k_sigma = 0.39839153 h/Mpc, an effective index n_eff = -1.79706366, and a curvature C = 0.32249769.
  !
  ! The quasi-linear term is compared only where it is not utterly negligible: it carries exp(-y²/8), so beyond k ~ 1/Mpc it
  ! falls below 10⁻¹⁰ of the linear spectrum and a relative comparison is meaningless. The halo term is compared from the point
  ! at which it rises above a per cent of the linear spectrum.
  integer                                                   , parameter    :: indexQuasiLinearLast =6         , indexHaloFirst=3
  double precision                                          , dimension(9) :: ratioQuasiLinearReference=[9.9355258005d-01,9.6003004697d-01,9.0696838061d-01,6.8051476724d-01,4.8764757689d-01,1.7226408970d-01,7.4327137875d-03,6.0378772087d-11,5.0110540068d-38]
  double precision                                          , dimension(9) :: ratioHaloReference       =[4.1929170682d-04,1.7321783167d-02,8.9490157537d-02,7.4211429283d-01,1.7752140594d+00,5.2494678839d+00,1.2564312461d+01,2.5403827896d+01,3.3296547504d+01]
  double precision                                          , dimension(9) :: ratio                              , ratioQuasiLinear     , &
       &                                                                      ratioHalo                          , powerLinear
  ! Tolerances.
  !! The two codes evaluate the same algorithm, but on linear spectra which agree only to the accuracy of the tabulation each
  !! builds from CAMB, and `halofit` amplifies differences in the linear spectrum through the nonlinear scale and the effective
  !! index and curvature there. The agreement achieved is better than 0.6% at every wavenumber - see the values reported by the
  !! test itself - so this leaves a factor of a few in hand for a different CAMB build.
  double precision                                          , parameter   :: tolerance             =2.0d-2
  !! The decomposition into quasi-linear and halo terms is exact.
  double precision                                          , parameter   :: toleranceDecomposition=1.0d-6
  !! The independent decomposition reproduces the CAMB total to better than 0.4%, which sets the accuracy of each term.
  double precision                                          , parameter   :: toleranceTerms        =2.0d-2
  double precision                                                        :: time
  character       (len=256                                 )              :: message
  integer                                                                 :: i

  call displayVerbositySet(verbosityLevelStandard)
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  call Unit_Tests_Begin_Group("Nonlinear power spectrum: Smith et al. (2003)")
  allocate(cosmologyParameters_               )
  allocate(cosmologyFunctions_                )
  allocate(darkMatterParticle_                )
  allocate(linearGrowth_                      )
  allocate(powerSpectrumPrimordial_           )
  allocate(transferFunction_                  )
  allocate(powerSpectrumPrimordialTransferred_)
  allocate(powerSpectrumWindowFunction_       )
  allocate(cosmologicalMassVariance_          )
  allocate(powerSpectrum_                     )
  allocate(nonlinear_                         )
  allocate(nonlinearQuasiLinear_              )
  allocate(nonlinearHalo_                     )
  cosmologyParameters_               =cosmologyParametersSimple               (OmegaMatter=OmegaMatter,OmegaBaryon=OmegaBaryon,OmegaDarkEnergy=OmegaDarkEnergy,temperatureCMB=temperatureCMB,HubbleConstant=HubbleConstant)
  cosmologyFunctions_                =cosmologyFunctionsMatterLambda          (cosmologyParameters_=cosmologyParameters_)
  darkMatterParticle_                =darkMatterParticleCDM                   ()
  linearGrowth_                      =linearGrowthCollisionlessMatter         (cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  powerSpectrumPrimordial_           =powerSpectrumPrimordialPowerLaw         (index_=indexPrimordial,running=0.0d0,runningRunning=0.0d0,wavenumberReference=1.0d0,runningSmallScalesOnly=.false.)
  transferFunction_                  =transferFunctionCAMB                    (darkMatterParticle_=darkMatterParticle_,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,transferFunctionType=transferFunctionTypeTotal,redshift=0.0d0,cambCountPerDecade=0)
  powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferredSimple(powerSpectrumPrimordial_=powerSpectrumPrimordial_,transferFunction_=transferFunction_,linearGrowth_=linearGrowth_)
  powerSpectrumWindowFunction_       =powerSpectrumWindowFunctionTopHat       (cosmologyParameters_=cosmologyParameters_)
  cosmologicalMassVariance_          =cosmologicalMassVarianceFilteredPower   (sigma8=sigma8,tolerance=1.0d-6,toleranceTopHat=1.0d-6,rootVarianceLogarithmicGradientTolerance=1.0d-9,integrationFailureIsFatal=.true.,storeTabulations=.true.,nonMonotonicIsFatal=.true.,monotonicInterpolation=.false.,truncateAtParticleHorizon=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,linearGrowth_=linearGrowth_,powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_,powerSpectrumWindowFunction_=powerSpectrumWindowFunction_)
  powerSpectrum_                     =powerSpectrumStandard                   (cosmologicalMassVariance_=cosmologicalMassVariance_,powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_)
  ! The full algorithm, and each of its two terms in isolation.
  nonlinear_                         =powerSpectrumNonlinearSmith2003         (includePeacockCorrection=.false.,includeQuasiLinearPower=.true. ,includeHaloPower=.true. ,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,powerSpectrum_=powerSpectrum_)
  nonlinearQuasiLinear_              =powerSpectrumNonlinearSmith2003         (includePeacockCorrection=.false.,includeQuasiLinearPower=.true. ,includeHaloPower=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,powerSpectrum_=powerSpectrum_)
  nonlinearHalo_                     =powerSpectrumNonlinearSmith2003         (includePeacockCorrection=.false.,includeQuasiLinearPower=.false.,includeHaloPower=.true. ,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,powerSpectrum_=powerSpectrum_)
  time                               =cosmologyFunctions_%cosmicTime(1.0d0)
  do i=1,countWavenumbers
     powerLinear     (i)=powerSpectrum_       %power(wavenumbers(i),time)
     ratio           (i)=nonlinear_           %value(wavenumbers(i),time)/powerLinear(i)
     ratioQuasiLinear(i)=nonlinearQuasiLinear_%value(wavenumbers(i),time)/powerLinear(i)
     ratioHalo       (i)=nonlinearHalo_       %value(wavenumbers(i),time)/powerLinear(i)
  end do
  ! Report the comparison, so that the margin by which it passes is visible and the tolerance above can be seen to be
  ! justified rather than simply generous.
  do i=1,countWavenumbers
     write (message,'(a,f9.5,a,f12.6,a,f12.6,a,f7.3,a)') '  k = ',wavenumbers(i),' /Mpc: P_nl/P_lin = ',ratio(i),', CAMB ',ratioReference(i),' (',100.0d0*(ratio(i)/ratioReference(i)-1.0d0),'%)'
     call displayMessage(trim(message))
  end do
  call Assert('nonlinear to linear ratio against CAMB',ratio,ratioReference,relTol=tolerance)
  ! The two terms must sum to the whole.
  call Assert('quasi-linear and halo terms sum to the total',ratioQuasiLinear+ratioHalo,ratio,relTol=toleranceDecomposition)
  ! Each term against the independent decomposition. Checking the total alone is not enough: the two terms overlap, so an error
  ! in one is diluted in the sum - a ten per cent error in the coefficient of the quasi-linear damping moves the total by only
  ! one per cent, but moves the quasi-linear term itself by up to thirteen.
  call Assert('quasi-linear term',ratioQuasiLinear(             1:indexQuasiLinearLast),ratioQuasiLinearReference(             1:indexQuasiLinearLast),relTol=toleranceTerms)
  call Assert('halo term'        ,ratioHalo       (indexHaloFirst:countWavenumbers    ),ratioHaloReference       (indexHaloFirst:countWavenumbers    ),relTol=toleranceTerms)
  ! On the largest scale the halo term is negligible and the quasi-linear term carries the linear spectrum, so the ratio
  ! approaches unity; on the smallest it is the halo term which dominates.
  call Assert('the halo term is negligible on large scales' ,ratioHalo       (1)                < 1.0d-2*ratioQuasiLinear(1               ),.true.              )
  call Assert('the quasi-linear term tends to linear'       ,ratioQuasiLinear(1)                                                           ,1.0d0 ,relTol=1.0d-2)
  call Assert('the halo term dominates on small scales'     ,ratioHalo       (countWavenumbers) > 1.0d+1*ratioQuasiLinear(countWavenumbers),.true.              )
  ! Clean up.
  deallocate(nonlinear_           )
  deallocate(nonlinearQuasiLinear_)
  deallocate(nonlinearHalo_       )
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Power_Spectrum_Nonlinear_Smith2003
