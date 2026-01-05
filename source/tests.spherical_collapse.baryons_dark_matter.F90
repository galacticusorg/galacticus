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

program Tests_Spherical_Collapse_Baryons_Dark_Matter
  !!{
  Tests linear growth calculations.
  !!}
  use :: Cosmological_Density_Field           , only : cosmologicalMassVarianceFilteredPower                     , criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy      , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Cosmology_Functions                  , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                 , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles                , only : darkMatterParticleCDM
  use :: Display                              , only : displayVerbositySet                                       , verbosityLevelStandard
  use :: Events_Hooks                         , only : eventsHooksInitialize
  use :: Intergalactic_Medium_Filtering_Masses, only : intergalacticMediumFilteringMassGnedin2000
  use :: Intergalactic_Medium_State           , only : intergalacticMediumStateSimple
  use :: Linear_Growth                        , only : componentDarkMatter                                       , linearGrowthBaryonsDarkMatter                                 , linearGrowthCollisionlessMatter
  use :: Power_Spectra_Primordial             , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred , only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions      , only : powerSpectrumWindowFunctionSharpKSpace
  use :: Spherical_Collapse_Solvers           , only : cllsnlssMttrDarkEnergyFixedAtTurnaround
  use :: Transfer_Functions                   , only : transferFunctionIdentity
  use :: Unit_Tests                           , only : Assert                                                    , Unit_Tests_Begin_Group                                        , Unit_Tests_End_Group                                        , Unit_Tests_Finish
  use :: Virial_Density_Contrast              , only : virialDensityContrastSphericalCollapseBrynsDrkMttrDrkEnrgy, virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  implicit none
  double precision                                                                , dimension(5), parameter :: redshift                                                   =[0.0d0,1.0d0,3.0d0,9.0d0,30.0d0]
  integer                                                                                       , parameter :: countFractionBaryon                                        =25
  double precision                                                                              , parameter :: fractionBaryonsMaximum                                     =0.5d0
  type            (cosmologyParametersSimple                                     )                          :: cosmologyParametersDMO_                                                                     , cosmologyParametersBaryons_
  type            (cosmologyFunctionsMatterLambda                                )                          :: cosmologyFunctionsMatterLambda_
  type            (cosmologicalMassVarianceFilteredPower                         )                          :: cosmologicalMassVarianceFilteredPower_
  type            (linearGrowthBaryonsDarkMatter                                 )                          :: linearGrowthBaryonsDarkMatter_
  type            (linearGrowthCollisionlessMatter                               )                          :: linearGrowthCollisionlessMatter_
  type            (intergalacticMediumStateSimple                                )                          :: intergalacticMediumState_
  type            (darkMatterParticleCDM                                         )                          :: darkMatterParticleCDM_
  type            (powerSpectrumWindowFunctionSharpKSpace                        )                          :: powerSpectrumWindowFunctionSharpKSpace_
  type            (powerSpectrumPrimordialPowerLaw                               )                          :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionIdentity                                      )                          :: transferFunctionIdentity_
  type            (powerSpectrumPrimordialTransferredSimple                      )                          :: powerSpectrumPrimordialTransferredSimple_
  type            (intergalacticMediumFilteringMassGnedin2000                    )                          :: intergalacticMediumFilteringMassGnedin2000_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt  )                          :: criticalOverdensitySphrclCllpsCllsnlssMttrCsmlgclCnstnt_
  type            (criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy      )                          :: criticalOverdensitySphrclCllpsBrynsDrkMttrDrkEnrgy_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)                          :: virialDensityContrastSphrclCllpsCllsnlssMttrCsmlgclCnstnt_
  type            (virialDensityContrastSphericalCollapseBrynsDrkMttrDrkEnrgy    )                          :: virialDensityContrastSphrclCllpsBrynsDrkMttrDrkEnrgy_
  character       (len=1024                                                      )                          :: message                                                                                     , outputFractions
  integer                                                                                                   :: i                                                                                           , j                         , &
       &                                                                                                       outputFile
  double precision                                                                                          :: expansionFactor                                                                             , time                      , &
       &                                                                                                       criticalOverdensityBaryons                                                                  , criticalOverdensityDMO    , &
       &                                                                                                       virialDensityContrastBaryons                                                                , virialDensityContrastDMO  , &
       &                                                                                                       radiusTurnaroundBaryons                                                                     , radiusTurnaroundDMO       , &
       &                                                                                                       fractionBaryons

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Get argument.
  if (Command_Argument_Count() > 0) then
     call Get_Command_Argument(1,outputFractions)
  else
     outputFractions="no"
  end if
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Open output file for z=0 results.
  if (trim(outputFractions) == "yes") then
     open(newUnit=outputFile,file='testSuite/outputs/sphericalCollapseBaryonsDarkMatter.txt',status='unknown',form='formatted')
     write (outputFile,'(a)') '# Spherical collapse solutions at z=0 as a function of baryon fraction'
     write (outputFile,'(a)') '#'
     write (outputFile,'(a)') '# fb        δc        Δᵥᵢᵣ      rₜₐ/rᵥᵢᵣ'
  end if
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: dark matter-only universe")
  ! Iterate over baryon fraction.
  do j=1,countFractionBaryon
     fractionBaryons=fractionBaryonsMaximum*dble(j-1)/dble(countFractionBaryon-1)
     if (trim(outputFractions) /= "yes" .and. fractionBaryons > 0.0d0) exit
     ! Construct model.
     darkMatterParticleCDM_                                     =darkMatterParticleCDM                                          (                                                                                   &
          &                                                                                                                     )
     cosmologyParametersDMO_                                    =cosmologyParametersSimple                                      (                                                                                   &
          &                                                                                                                      OmegaMatter                          = 0.2750d0                                  , &
          &                                                                                                                      OmegaBaryon                          = 0.2750d0*fractionBaryons                  , &
          &                                                                                                                      OmegaDarkEnergy                      = 0.7250d0                                  , &
          &                                                                                                                      temperatureCMB                       = 2.7800d0                                  , &
          &                                                                                                                      HubbleConstant                       =70.2000d0                                    &
          &                                                                                                                     )
     cosmologyParametersBaryons_                                =cosmologyParametersSimple                                      (                                                                                   &
          &                                                                                                                      OmegaMatter                          = 0.2750d0                                  , &
          &                                                                                                                      OmegaBaryon                          = 0.0458d0                                  , &
          &                                                                                                                      OmegaDarkEnergy                      = 0.7250d0                                  , &
          &                                                                                                                      temperatureCMB                       = 2.7800d0                                  , &
          &                                                                                                                      HubbleConstant                       =70.2000d0                                    &
          &                                                                                                                     )
     cosmologyFunctionsMatterLambda_                            =cosmologyFunctionsMatterLambda                                 (                                                                                   &
          &                                                                                                                      cosmologyParametersBaryons_                                                        &
          &                                                                                                                     )
     intergalacticMediumState_                                  =intergalacticMediumStateSimple                                 (                                                                                   &
          &                                                                                                                      reionizationRedshift                 = 8.00d0                                    , &
          &                                                                                                                      reionizationTemperature              = 1.00d4                                    , &
          &                                                                                                                      preReionizationTemperature           = 1.00d4                                    , &
          &                                                                                                                      cosmologyFunctions_                  =cosmologyFunctionsMatterLambda_            , &
          &                                                                                                                      cosmologyParameters_                 =cosmologyParametersBaryons_                  &
          &                                                                                                                     )
     linearGrowthBaryonsDarkMatter_                             =linearGrowthBaryonsDarkMatter                                  (                                                                                   &
          &                                                                                                                      redshiftInitial                      =160.0d0                                    , &
          &                                                                                                                      redshiftInitialDelta                 =  1.6d0                                    , &
          &                                                                                                                      cambCountPerDecade                   =  0                                        , &
          &                                                                                                                      darkMatterOnlyInitialConditions      =.false.                                    , &
          &                                                                                                                      cosmologyParameters_                 =cosmologyParametersBaryons_                , &
          &                                                                                                                      cosmologyParametersInitialConditions_=cosmologyParametersBaryons_                , &
          &                                                                                                                      cosmologyFunctions_                  =cosmologyFunctionsMatterLambda_            , &
          &                                                                                                                      intergalacticMediumState_            =intergalacticMediumState_                    &
          &                                                                                                                     )
     linearGrowthCollisionlessMatter_                           =linearGrowthCollisionlessMatter                                (                                                                                   &
          &                                                                                                                      cosmologyParameters_                 =cosmologyParametersBaryons_                , &
          &                                                                                                                      cosmologyFunctions_                  =cosmologyFunctionsMatterLambda_              &
          &                                                                                                                     )
     powerSpectrumPrimordialPowerLaw_                           =powerSpectrumPrimordialPowerLaw                                (                                                                                   &
          &                                                                                                                      index_                               =-1.0d0                                     , &
          &                                                                                                                      running                              =+0.0d0                                     , &
          &                                                                                                                      runningRunning                       =+0.0d0                                     , &
          &                                                                                                                      wavenumberReference                  =+1.0d0                                     , &
          &                                                                                                                      runningSmallScalesOnly               =.false.                                      &
          &                                                                                                                     )
     transferFunctionIdentity_                                  =transferFunctionIdentity                                       (                                                                                   &
          &                                                                                                                      cosmologyParameters_                 =cosmologyParametersBaryons_                , &
          &                                                                                                                      time                                 =13.8d0                                       &
          &                                                                                                                     )
     powerSpectrumPrimordialTransferredSimple_                  =powerSpectrumPrimordialTransferredSimple                       (                                                                                   &
          &                                                                                                                      powerSpectrumPrimordial_             =powerSpectrumPrimordialPowerLaw_           , &
          &                                                                                                                      transferFunction_                    =transferFunctionIdentity_                  , &
          &                                                                                                                      linearGrowth_                        =linearGrowthBaryonsDarkMatter_               &
          &                                                                                                                     )
     powerSpectrumWindowFunctionSharpKSpace_                    =powerSpectrumWindowFunctionSharpKSpace                         (                                                                                   &
          &                                                                                                                      cosmologyParameters_                 =cosmologyParametersBaryons_                , &
          &                                                                                                                      normalization                        =0.0d0                                        &
          &                                                                                                                     )

     cosmologicalMassVarianceFilteredPower_                     =cosmologicalMassVarianceFilteredPower                          (                                                                                   &
          &                                                                                                                      sigma8                               =1.0d+0                                     , &
          &                                                                                                                      tolerance                            =1.0d-4                                     , &
          &                                                                                                                      toleranceTopHat                      =1.0d-4                                     , &
          &                                                                                                                      nonMonotonicIsFatal                  =.true.                                     , &
          &                                                                                                                      monotonicInterpolation               =.false.                                    , &
          &                                                                                                                      truncateAtParticleHorizon            =.false.                                    , &
          &                                                                                                                      cosmologyParameters_                 =cosmologyParametersBaryons_                , &
          &                                                                                                                      cosmologyFunctions_                  =cosmologyFunctionsMatterLambda_            , &
          &                                                                                                                      linearGrowth_                        =linearGrowthBaryonsDarkMatter_             , &
          &                                                                                                                      powerSpectrumPrimordialTransferred_  =powerSpectrumPrimordialTransferredSimple_  , &
          &                                                                                                                      powerSpectrumWindowFunction_         =powerSpectrumWindowFunctionSharpKSpace_      &
          &                                                                                                                     )
     intergalacticMediumFilteringMassGnedin2000_                =intergalacticMediumFilteringMassGnedin2000                     (                                                                                   &
          &                                                                                                                      timeTooEarlyIsFatal                  =.true.                                     , &
          &                                                                                                                      cosmologyParameters_                 =cosmologyParametersBaryons_                , &
          &                                                                                                                      cosmologyFunctions_                  =cosmologyFunctionsMatterLambda_            , &
          &                                                                                                                      linearGrowth_                        =linearGrowthBaryonsDarkMatter_             , &
          &                                                                                                                      intergalacticMediumState_            =intergalacticMediumState_                    &
          &                                                                                                                     )
     criticalOverdensitySphrclCllpsCllsnlssMttrCsmlgclCnstnt_   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt   (                                                                                   &
          &                                                                                                                      linearGrowth_                        =linearGrowthCollisionlessMatter_           , &
          &                                                                                                                      cosmologyFunctions_                  =cosmologyFunctionsMatterLambda_            , &
          &                                                                                                                      cosmologicalMassVariance_            =cosmologicalMassVarianceFilteredPower_     , &
          &                                                                                                                      darkMatterParticle_                  =darkMatterParticleCDM_                     , &
          &                                                                                                                      tableStore                           =.true.                                       &
          &                                                                                                                     )
     criticalOverdensitySphrclCllpsBrynsDrkMttrDrkEnrgy_        =criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy       (                                                                                   &
          &                                                                                                                      cosmologyParameters_                 =cosmologyParametersDMO_                    , &
          &                                                                                                                      cosmologyFunctions_                  =cosmologyFunctionsMatterLambda_            , &
          &                                                                                                                      cosmologicalMassVariance_            =cosmologicalMassVarianceFilteredPower_     , &
          &                                                                                                                      darkMatterParticle_                  =darkMatterParticleCDM_                     , &
          &                                                                                                                      intergalacticMediumFilteringMass_    =intergalacticMediumFilteringMassGnedin2000_, &
          &                                                                                                                      tablePointsPerOctave                 =300                                        , &
          &                                                                                                                      tableStore                           =.false.                                    , &
          &                                                                                                                      energyFixedAt                        =cllsnlssMttrDarkEnergyFixedAtTurnaround    , &
          &                                                                                                                      normalization                        =1.0d0                                        &
          &                                                                                                                     )
     virialDensityContrastSphrclCllpsCllsnlssMttrCsmlgclCnstnt_ =virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                                    &
          &                                                                                                                      tableStore                           =.false.                                    , &
          &                                                                                                                      cosmologyFunctions_                  =cosmologyFunctionsMatterLambda_              &
          &                                                                                                                     )

     virialDensityContrastSphrclCllpsBrynsDrkMttrDrkEnrgy_      =virialDensityContrastSphericalCollapseBrynsDrkMttrDrkEnrgy     (                                                                                   &
          &                                                                                                                      tableStore                           =.true.                                     , &
          &                                                                                                                      tablePointsPerOctave                 =300                                        , &
          &                                                                                                                      energyFixedAt                        =cllsnlssMttrDarkEnergyFixedAtTurnaround    , &
          &                                                                                                                      cosmologyParameters_                 =cosmologyParametersDMO_                    , &
          &                                                                                                                      cosmologyFunctions_                  =cosmologyFunctionsMatterLambda_            , &
          &                                                                                                                      intergalacticMediumFilteringMass_    =intergalacticMediumFilteringMassGnedin2000_  &
          &                                                                                                                     )
     ! Iterate over redshifts.
     do i=1,size(redshift)
        expansionFactor             =cosmologyFunctionsMatterLambda_                           %expansionFactorFromRedshift(     redshift       (i)           )
        time                        =cosmologyFunctionsMatterLambda_                           %cosmicTime                 (     expansionFactor              )
        criticalOverdensityDMO      =criticalOverdensitySphrclCllpsCllsnlssMttrCsmlgclCnstnt_  %value                      (time=time              ,mass=1.0d8)
        criticalOverdensityBaryons  =criticalOverdensitySphrclCllpsBrynsDrkMttrDrkEnrgy_       %value                      (time=time              ,mass=1.0d8)
        virialDensityContrastDMO    =virialDensityContrastSphrclCllpsCllsnlssMttrCsmlgclCnstnt_%densityContrast            (time=time              ,mass=1.0d8)
        virialDensityContrastBaryons=virialDensityContrastSphrclCllpsBrynsDrkMttrDrkEnrgy_     %densityContrast            (time=time              ,mass=1.0d8)
        radiusTurnaroundDMO         =virialDensityContrastSphrclCllpsCllsnlssMttrCsmlgclCnstnt_%turnAroundOverVirialRadii  (time=time              ,mass=1.0d8)
        radiusTurnaroundBaryons     =virialDensityContrastSphrclCllpsBrynsDrkMttrDrkEnrgy_     %turnAroundOverVirialRadii  (time=time              ,mass=1.0d8)
        if (fractionBaryons <= 0.0d0) then
           write (message,'(a,f6.1,a)') "critical overdensity for collapse [z = ",redshift(i),"]"
           call Assert(trim(message),criticalOverdensityBaryons  ,criticalOverdensityDMO  ,relTol=1.0d-3)
           write (message,'(a,f6.1,a)') "virial density contrast           [z = ",redshift(i),"]"
           call Assert(trim(message),virialDensityContrastBaryons,virialDensityContrastDMO,relTol=1.0d-3)
           write (message,'(a,f6.1,a)') "turnaround radius                 [z = ",redshift(i),"]"
           call Assert(trim(message),radiusTurnaroundBaryons     ,radiusTurnaroundDMO     ,relTol=1.0d-3)
        end if
        if (trim(outputFractions) == "yes" .and. redshift(i) == 0.0d0) then
           write (outputFile,'(2x,4(f9.5,1x))') fractionBaryons,criticalOverdensityBaryons,virialDensityContrastBaryons,radiusTurnaroundBaryons
           if (fractionBaryons > 0.0d0) exit
        end if
     end do
  end do
  if (trim(outputFractions) == "yes") close(outputFile)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Spherical_Collapse_Baryons_Dark_Matter
