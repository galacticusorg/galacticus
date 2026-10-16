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
  Contains a program to test the :cite:t:`bhattacharya_mass_2011` halo mass function.
  !!}

program Test_Halo_Mass_Function_Bhattacharya2011
  !!{RST
  Test the :galacticus-class:`haloMassFunctionBhattacharya2011` halo mass function against the independent implementation in the
  ``colossus`` package :cite:p:`diemer_colossus_2018`, over a range of redshift.

  The comparison is made in terms of the multiplicity function, :math:`f(\sigma)`, defined by

  .. math::

     \frac{\mathrm{d}n}{\mathrm{d}M} = \frac{\bar{\rho}_\mathrm{m}}{M^2} f(\sigma) \left| \frac{\mathrm{d}\ln\sigma}{\mathrm{d}\ln
     M} \right|,

  which is what :cite:t:`bhattacharya_mass_2011` actually fit. Recovering :math:`f(\sigma)` from the differential mass function
  divides out the matter density and the gradient of the variance, leaving only the fitting function itself.

  The masses at which it is evaluated are those at which :math:`\sigma` takes a set of chosen values, obtained by inverting
  :math:`\sigma(M)`. The comparison is therefore made at exactly the same :math:`\sigma` as the reference, and so does not depend
  on the power spectrum, the transfer function, or the normalization of the variance - any of which could change without this
  test either noticing or caring. The critical overdensity is held fixed at 1.68647, the value ``colossus`` uses, since
  :math:`f` depends on it through :math:`\nu = (\delta_\mathrm{c}/\sigma)^2`.

  The point of testing over redshift is that two of the fit parameters evolve: :math:`\bar{A} \propto (1+z)^{-0.11}` and
  :math:`\bar{a} \propto (1+z)^{-0.01}`. That evolution was absent from this class until it was added alongside this test;
  without it the mass function is too high by 8% at :math:`z=2` and 17% at :math:`z=3`. Setting the two exponents to zero must
  recover the un-evolving fit, which is checked here too, since that is how the reference parameter files are configured - the
  parameters there having been calibrated without any redshift dependence.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower   , criticalOverdensityFixed
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Halo_Mass_Functions                 , only : haloMassFunctionBhattacharya2011
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group  , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple               ), pointer        :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda          ), pointer        :: cosmologyFunctions_
  type            (darkMatterParticleCDM                   ), pointer        :: darkMatterParticle_
  type            (linearGrowthCollisionlessMatter         ), pointer        :: linearGrowth_
  type            (powerSpectrumPrimordialPowerLaw         ), pointer        :: powerSpectrumPrimordial_
  type            (transferFunctionEisensteinHu1999        ), pointer        :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple), pointer        :: powerSpectrumPrimordialTransferred_
  type            (powerSpectrumWindowFunctionTopHat       ), pointer        :: powerSpectrumWindowFunction_
  type            (cosmologicalMassVarianceFilteredPower   ), pointer        :: cosmologicalMassVariance_
  type            (criticalOverdensityFixed                ), pointer        :: criticalOverdensity_
  type            (haloMassFunctionBhattacharya2011        ), pointer        :: haloMassFunction_                     , haloMassFunctionNoEvolution_
  ! The critical overdensity used by colossus, held fixed here so that both sides form the same ν.
  double precision                                          , parameter      :: criticalOverdensityValue  =1.68647d0
  ! Values of sigma at which the multiplicity function is compared, and the redshifts.
  integer                                                   , parameter      :: countSigma                =5        , countRedshift            =4
  double precision                                          , dimension(5  ) :: sigmas                =[0.4d0,0.7d0,1.0d0,1.5d0,2.5d0]
  double precision                                          , dimension(  4) :: redshifts             =[0.0d0,1.0d0,2.0d0,3.0d0]
  ! The multiplicity function given by colossus for the bhattacharya11 model at those σ and redshifts.
  double precision                                          , dimension(5,4) :: multiplicityReference=reshape(                                                             &
       &                                                                       [2.8859933239541d-03,1.3659755242365d-01,2.7196260748573d-01,3.2231174064359d-01,2.9086422460145d-01, &
       &                                                                        2.7909338461091d-03,1.2795004628102d-01,2.5286569866108d-01,2.9865751196321d-01,2.6921857315884d-01, &
       &                                                                        2.7363652337684d-03,1.2314171137574d-01,2.4231733110677d-01,2.8563169735633d-01,2.5731108138665d-01, &
       &                                                                        2.6981117135873d-03,1.1983756526209d-01,2.3509884268217d-01,2.7673483865755d-01,2.4918327131625d-01],[5,4])
  double precision                                          , dimension(5,4) :: multiplicity                       , multiplicityNoEvolution
  double precision                                          , dimension(5,4) :: rootVarianceRecovered              , sigmasTarget
  ! Tolerance. Both sides evaluate the same closed-form expression, so only the precision to which the reference values are
  ! written, and the accuracy of the inversion of σ(M), enter.
  double precision                                          , parameter      :: tolerance                 =1.0d-6
  !! The accuracy to which σ(M) is inverted, which the comparison above amplifies.
  double precision                                          , parameter      :: toleranceInversion        =1.0d-10
  double precision                                                           :: time                               , mass                , &
       &                                                                        rootVariance                       , rootVarianceGradient, &
       &                                                                        densityMatter
  integer                                                                    :: i                                  , j                   , &
       &                                                                        k
  !! The number of Newton refinements applied to the inversion of σ(M).
  integer                                                   , parameter      :: countRefinements          =4

  call displayVerbositySet(verbosityLevelStandard)
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  call Unit_Tests_Begin_Group("Halo mass function: Bhattacharya et al. (2011)")
  allocate(cosmologyParameters_               )
  allocate(cosmologyFunctions_                )
  allocate(darkMatterParticle_                )
  allocate(linearGrowth_                      )
  allocate(powerSpectrumPrimordial_           )
  allocate(transferFunction_                  )
  allocate(powerSpectrumPrimordialTransferred_)
  allocate(powerSpectrumWindowFunction_       )
  allocate(cosmologicalMassVariance_          )
  allocate(criticalOverdensity_               )
  allocate(haloMassFunction_                  )
  allocate(haloMassFunctionNoEvolution_       )
  cosmologyParameters_               =cosmologyParametersSimple               (OmegaMatter=0.3153d0,OmegaBaryon=0.0493d0,OmegaDarkEnergy=0.6847d0,temperatureCMB=2.72548d0,HubbleConstant=67.36d0)
  cosmologyFunctions_                =cosmologyFunctionsMatterLambda          (cosmologyParameters_=cosmologyParameters_)
  darkMatterParticle_                =darkMatterParticleCDM                   ()
  linearGrowth_                      =linearGrowthCollisionlessMatter         (cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  powerSpectrumPrimordial_           =powerSpectrumPrimordialPowerLaw         (index_=0.9649d0,running=0.0d0,runningRunning=0.0d0,wavenumberReference=1.0d0,runningSmallScalesOnly=.false.)
  transferFunction_                  =transferFunctionEisensteinHu1999        (neutrinoNumberEffective=3.046d0,neutrinoMassSummed=0.0d0,darkMatterParticle_=darkMatterParticle_,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferredSimple(powerSpectrumPrimordial_=powerSpectrumPrimordial_,transferFunction_=transferFunction_,linearGrowth_=linearGrowth_)
  powerSpectrumWindowFunction_       =powerSpectrumWindowFunctionTopHat       (cosmologyParameters_=cosmologyParameters_)
  cosmologicalMassVariance_          =cosmologicalMassVarianceFilteredPower   (sigma8=0.8111d0,tolerance=1.0d-9,toleranceTopHat=1.0d-9,rootVarianceLogarithmicGradientTolerance=1.0d-9,integrationFailureIsFatal=.true.,storeTabulations=.true.,nonMonotonicIsFatal=.true.,monotonicInterpolation=.false.,truncateAtParticleHorizon=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,linearGrowth_=linearGrowth_,powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_,powerSpectrumWindowFunction_=powerSpectrumWindowFunction_)
  criticalOverdensity_               =criticalOverdensityFixed                (criticalOverdensity_=criticalOverdensityValue,linearGrowth_=linearGrowth_,cosmologyFunctions_=cosmologyFunctions_,cosmologicalMassVariance_=cosmologicalMassVariance_)
  ! The fit as published, with the redshift dependence of its normalization and `a` parameter.
  haloMassFunction_                  =haloMassFunctionBhattacharya2011        (cosmologyParameters_=cosmologyParameters_,cosmologicalMassVariance_=cosmologicalMassVariance_,criticalOverdensity_=criticalOverdensity_,cosmologyFunctions_=cosmologyFunctions_,a=0.788d0,b=1.0d0,c=1.0d0,p=0.807d0,q=1.795d0,normalization=0.333d0,exponentRedshiftA=-0.01d0,exponentRedshiftNormalization=-0.11d0)
  ! The same fit with that evolution switched off, as the reference parameter files configure it.
  haloMassFunctionNoEvolution_       =haloMassFunctionBhattacharya2011        (cosmologyParameters_=cosmologyParameters_,cosmologicalMassVariance_=cosmologicalMassVariance_,criticalOverdensity_=criticalOverdensity_,cosmologyFunctions_=cosmologyFunctions_,a=0.788d0,b=1.0d0,c=1.0d0,p=0.807d0,q=1.795d0,normalization=0.333d0,exponentRedshiftA=+0.00d0,exponentRedshiftNormalization=+0.00d0)
  densityMatter                      =+cosmologyParameters_%OmegaMatter    () &
       &                              *cosmologyParameters_%densityCritical()
  do j=1,countRedshift
     time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshifts(j)))
     do i=1,countSigma
        ! Find the mass at which the variance takes the required value, so that the comparison is made at exactly the σ of
        ! the reference rather than at whatever σ this power spectrum happens to give.
        ! `mass()` inverts a tabulation of σ(M) and so lands within a few parts in 10,000 of the target. That is not close
        ! enough here: the multiplicity varies as steeply as σ¹² at the smallest σ used, so such an error would show up
        ! at the 0.1% level in the comparison below. Refine by Newton iteration on ln(σ) against ln(M), for which the
        ! logarithmic gradient is already to hand; this converges to round-off in a step or two.
        mass=cosmologicalMassVariance_%mass(sigmas(i),time)
        do k=1,countRefinements
           call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(mass,time,rootVariance,rootVarianceGradient)
           mass=+mass                              &
                &*exp(                             &
                &     +log(sigmas(i)/rootVariance) &
                &     /rootVarianceGradient        &
                &    )
        end do
        call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(mass,time,rootVariance,rootVarianceGradient)
        rootVarianceRecovered(i,j)=rootVariance
        sigmasTarget         (i,j)=sigmas(i)
        ! Recover the multiplicity function from the differential mass function.
        multiplicity           (i,j)=+haloMassFunction_           %differential(time,mass) &
             &                       *mass                                             **2 &
             &                       /densityMatter                                        &
             &                       /abs(rootVarianceGradient)
        multiplicityNoEvolution(i,j)=+haloMassFunctionNoEvolution_%differential(time,mass) &
             &                       *mass                                             **2 &
             &                       /densityMatter                                        &
             &                       /abs(rootVarianceGradient)
     end do
  end do
  ! The comparison below is made at the sigma of the reference, so first confirm that inverting sigma(M) really landed there.
  ! The multiplicity varies as steeply as σ¹² at the smallest sigma used here, so an error in the inversion is amplified
  ! by that factor in the comparison which follows.
  call Assert('inversion of sigma(M) recovers the target',rootVarianceRecovered,sigmasTarget         ,relTol=toleranceInversion)
  call Assert('multiplicity function against colossus'   ,multiplicity         ,multiplicityReference,relTol=tolerance         )
  ! With the evolution switched off the multiplicity must be independent of redshift, and equal to the published fit at z=0.
  do j=2,countRedshift
     call Assert('no redshift evolution when the exponents vanish',multiplicityNoEvolution(:,j),multiplicityNoEvolution(:,1),relTol=tolerance)
  end do
  call Assert('vanishing exponents recover the z=0 fit',multiplicityNoEvolution(:,1),multiplicityReference(:,1),relTol=tolerance)
  ! The published evolution suppresses the mass function with increasing redshift.
  call Assert('the published fit is suppressed at high redshift',all(multiplicity(:,countRedshift) < multiplicityNoEvolution(:,countRedshift)),.true.)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Halo_Mass_Function_Bhattacharya2011
