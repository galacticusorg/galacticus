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

!!{RST
Tests of the large-:math:`x` (:math:`x>x_\mathrm{critical}=100`) asymptotic deceleration factor, :math:`\chi_\mathrm{d}`, of the :cite:t:`kummer_effective_2018` self-interacting dark matter satellite deceleration model, for both the velocity-dependent and constant cross-section models. In this regime the factor has a closed analytic form, dispatched on the dark matter particle class; the reference values are evaluated directly from those forms.
!!}

program Tests_SIDM_Kummer_Deceleration
  use :: Cosmology_Parameters             , only : cosmologyParametersSimple
  use :: Cosmology_Functions              , only : cosmologyFunctionsMatterLambda
  use :: Virial_Density_Contrast          , only : virialDensityContrastBryanNorman1998
  use :: Dark_Matter_Halo_Scales          , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profiles_DMO         , only : darkMatterProfileDMONFW
  use :: Dark_Matter_Particles            , only : darkMatterParticleCDM                             , darkMatterParticleSIDMVelocityDependent, darkMatterParticleSelfInteractingDarkMatterConstant
  use :: Satellite_Deceleration_SIDM      , only : satelliteDecelerationSIDMKummer2018
  use :: Display                          , only : displayVerbositySet                               , verbosityLevelStandard
  use :: Events_Hooks                     , only : eventsHooksInitialize
  use :: Functions_Global_Utilities       , only : Functions_Global_Set
  use :: Galacticus_Nodes                 , only : nodeClassHierarchyInitialize
  use :: Node_Components                  , only : Node_Components_Initialize                        , Node_Components_Thread_Initialize      , Node_Components_Thread_Uninitialize                , Node_Components_Uninitialize
  use :: Input_Parameters                 , only : inputParameters
  use :: Unit_Tests                       , only : Assert                                            , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group                               , Unit_Tests_Finish
  implicit none
  type(inputParameters                                    )          :: parameters
  type(cosmologyParametersSimple                          ), pointer :: cosmologyParameters_
  type(cosmologyFunctionsMatterLambda                     ), pointer :: cosmologyFunctions_
  type(virialDensityContrastBryanNorman1998               ), pointer :: virialDensityContrast_
  type(darkMatterHaloScaleVirialDensityContrastDefinition ), pointer :: darkMatterHaloScale_
  type(darkMatterProfileDMONFW                            ), pointer :: darkMatterProfileDMO_
  type(darkMatterParticleCDM                              ), pointer :: darkMatterParticleCDM_
  type(darkMatterParticleSIDMVelocityDependent            ), pointer :: particleVelocityDependent
  type(darkMatterParticleSelfInteractingDarkMatterConstant), pointer :: particleConstant
  type(satelliteDecelerationSIDMKummer2018                )          :: decelerationVelocityDependent, decelerationConstant

  ! Set verbosity level and initialize the minimal global state needed to construct the objects.
  call displayVerbositySet  (verbosityLevelStandard)
  call eventsHooksInitialize(                      )
  call Functions_Global_Set (                      )
  ! Initialize the node-component hierarchy (the NFW profile required to construct the deceleration object needs the dark matter
  ! profile component, with a settable scale radius, to be active).
  parameters=inputParameters('testSuite/parameters/nodes/nodes_SIDM_parametric.xml')
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("SIDM Kummer (2018) deceleration factor")
  ! Build the dependency stack (a cosmology, halo scale, and NFW profile are required to construct the deceleration object, but
  ! are not used by the large-x asymptotic branch under test).
  allocate(cosmologyParameters_     )
  allocate(cosmologyFunctions_      )
  allocate(virialDensityContrast_   )
  allocate(darkMatterHaloScale_     )
  allocate(darkMatterProfileDMO_    )
  allocate(darkMatterParticleCDM_   )
  allocate(particleVelocityDependent)
  allocate(particleConstant         )
  cosmologyParameters_     =cosmologyParametersSimple                          (OmegaMatter=0.2815d0,OmegaBaryon=0.0465d0,OmegaDarkEnergy=0.7185d0,temperatureCMB=2.78d0,HubbleConstant=70.0d0)
  cosmologyFunctions_      =cosmologyFunctionsMatterLambda                     (cosmologyParameters_=cosmologyParameters_)
  virialDensityContrast_   =virialDensityContrastBryanNorman1998               (allowUnsupportedCosmology=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  darkMatterHaloScale_     =darkMatterHaloScaleVirialDensityContrastDefinition (cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,virialDensityContrast_=virialDensityContrast_)
  darkMatterProfileDMO_    =darkMatterProfileDMONFW                            (velocityDispersionUseSeriesExpansion=.false.,darkMatterHaloScale_=darkMatterHaloScale_)
  darkMatterParticleCDM_   =darkMatterParticleCDM                              ()
  particleVelocityDependent=darkMatterParticleSIDMVelocityDependent            (velocityCharacteristic=50.0d0,sigma0=10.0d0,darkMatterParticle_=darkMatterParticleCDM_)
  particleConstant         =darkMatterParticleSelfInteractingDarkMatterConstant(crossSectionSelfInteraction=10.0d0,darkMatterParticle_=darkMatterParticleCDM_)
  ! Build one deceleration object per cross-section model.
  decelerationVelocityDependent=satelliteDecelerationSIDMKummer2018(cosmologyParameters_,particleVelocityDependent,darkMatterHaloScale_,darkMatterProfileDMO_)
  decelerationConstant         =satelliteDecelerationSIDMKummer2018(cosmologyParameters_,particleConstant         ,darkMatterHaloScale_,darkMatterProfileDMO_)
  ! Velocity-dependent large-x deceleration factor: q = speedOrbital/v_c with v_c = 50 km/s.
  call Assert("velocity-dependent deceleration factor, x=200, q=2",decelerationVelocityDependent%decelerationFactor(200.0d0,100.0d0),9.998333516667d-01,relTol=1.0d-9)
  call Assert("velocity-dependent deceleration factor, x=150, q=1",decelerationVelocityDependent%decelerationFactor(150.0d0, 50.0d0),9.998814920165d-01,relTol=1.0d-9)
  ! Constant large-x deceleration factor: independent of the orbital speed.
  call Assert("constant deceleration factor, x=200",decelerationConstant%decelerationFactor(200.0d0,100.0d0),9.999333353333d-01,relTol=1.0d-9)
  call Assert("constant deceleration factor, x=150",decelerationConstant%decelerationFactor(150.0d0, 50.0d0),9.998814878025d-01,relTol=1.0d-9)
  ! Clean up the node-component hierarchy.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_SIDM_Kummer_Deceleration
