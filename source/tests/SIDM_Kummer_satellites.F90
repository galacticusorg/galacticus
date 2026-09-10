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
Tests of the :cite:t:`kummer_effective_2018` self-interacting dark matter satellite evaporation and deceleration models, as
applied to an explicitly constructed satellite/host pair. (The deceleration factor, :math:`\chi_\mathrm{d}`, is tested in
isolation by ``tests.SIDM_Kummer_deceleration.exe``.)

First, the tabulated evaporation factor, :math:`\chi_\mathrm{e}`, is compared against the closed form
:math:`\chi_\mathrm{e}=(1-x^2)/(1+x^2)` which it takes for an isotropic, constant cross section. The mass loss rate and
acceleration computed for the satellite are then used to pin the construction of :math:`x` itself. Since
:math:`x=v_\mathrm{esc}/v_\mathrm{rel}`, in the limit of arbitrarily large orbital speed :math:`x \rightarrow 0`, so that the
evaporation factor implied by the mass loss rate must approach unity while that implied by the acceleration must approach
:math:`\chi_\mathrm{d}(0)=-1`; and for a satellite bound deeply enough that :math:`v_\mathrm{esc} > v_\mathrm{rel}`
evaporation must switch off entirely while the drag must remain a genuine deceleration.
!!}

program Tests_SIDM_Kummer_Satellites
  use :: Coordinates                     , only : coordinateCartesian                               , assignment(=)
  use :: Cosmology_Parameters            , only : cosmologyParametersSimple
  use :: Cosmology_Functions             , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Particles           , only : darkMatterParticleCDM                             , darkMatterParticleSelfInteractingDarkMatterConstant
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMONFW
  use :: Display                         , only : displayVerbositySet                               , verbosityLevelStandard
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Functions_Global_Utilities      , only : Functions_Global_Set
  use :: Galacticus_Nodes                , only : mergerTree                                        , nodeClassHierarchyInitialize                       , nodeComponentBasic    , nodeComponentDarkMatterProfile, &
       &                                          nodeComponentSatellite                            , treeNode
  use :: Input_Parameters                , only : inputParameters
  use :: Mass_Distributions              , only : massDistributionClass
  use :: Node_Components                 , only : Node_Components_Initialize                        , Node_Components_Thread_Initialize                  , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Numerical_Constants_Astronomical, only : gigaYear                                          , massSolar                                          , megaParsec
  use :: Numerical_Constants_Prefixes    , only : centi                                             , kilo                                               , milli
  use :: Satellite_Deceleration_SIDM     , only : satelliteDecelerationSIDMKummer2018
  use :: Satellite_Evaporation_SIDM      , only : satelliteEvaporationSIDMKummer2018
  use :: Unit_Tests                      , only : Assert                                            , Unit_Tests_Begin_Group                             , Unit_Tests_End_Group  , Unit_Tests_Finish
  use :: Virial_Density_Contrast         , only : virialDensityContrastBryanNorman1998
  implicit none
  type            (inputParameters                                    )               :: parameters
  type            (cosmologyParametersSimple                          ), pointer      :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                     ), pointer      :: cosmologyFunctions_
  type            (virialDensityContrastBryanNorman1998               ), pointer      :: virialDensityContrast_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition ), pointer      :: darkMatterHaloScale_
  type            (darkMatterProfileDMONFW                            ), pointer      :: darkMatterProfileDMO_
  type            (darkMatterParticleCDM                              ), pointer      :: darkMatterParticleCDM_
  type            (darkMatterParticleSelfInteractingDarkMatterConstant), pointer      :: particleConstant
  type            (satelliteEvaporationSIDMKummer2018                 )               :: evaporation
  type            (satelliteDecelerationSIDMKummer2018                )               :: deceleration
  type            (mergerTree                                         )               :: tree
  type            (treeNode                                           ), pointer      :: node                        , nodeHost                    , &
       &                                                                                 nodeDense
  class           (nodeComponentBasic                                 ), pointer      :: basic                       , basicHost                   , &
       &                                                                                 basicDense
  class           (nodeComponentDarkMatterProfile                     ), pointer      :: darkMatterProfile           , darkMatterProfileHost       , &
       &                                                                                 darkMatterProfileDense
  class           (nodeComponentSatellite                             ), pointer      :: satellite                   , satelliteDense
  class           (massDistributionClass                              ), pointer      :: massDistributionHost_
  double precision                                                     , parameter    :: crossSectionConstant=10.0d0                                   ! [cm² g⁻¹]
  double precision                                                     , parameter    :: massHost            = 1.0d12, radiusScaleHost      =2.0d-2, & ! [M☉], [Mpc]
       &                                                                                 massSatellite       = 1.0d09, radiusScaleSatellite =2.0d-3, & ! [M☉], [Mpc]
       &                                                                                 massSatelliteDense  = 1.0d11, radiusScaleDense     =1.0d-3    ! [M☉], [Mpc]
  double precision                                                     , parameter    :: radiusOrbital       = 5.0d-2                                  ! [Mpc]
  double precision                                                     , parameter    :: speedOrbitalDense   = 1.0d1                                   ! [km/s]
  double precision                                                     , dimension(7) :: xTabulated          =[0.10d0,0.30d0,0.50d0,0.70d0,0.90d0,1.05d0,2.00d0]
  double precision                                                     , dimension(7) :: factorTabulated             , factorAnalytic
  double precision                                                     , dimension(3) :: speedOrbital        =[5.0d1,2.0d2,5.0d3]                      ! [km/s]
  double precision                                                     , dimension(3) :: factorEvaporation           , factorDeceleration
  double precision                                                     , dimension(3) :: acceleration
  type            (coordinateCartesian                                )               :: coordinates
  double precision                                                                    :: densityHost                 , rateNormalization         , &
       &                                                                                 massLossRate                , massBoundary              , &
       &                                                                                 fractionDarkMatter          , massLossRateDense         , &
       &                                                                                 factorDecelerationDense
  integer                                                                             :: i

  ! Set verbosity level and initialize the minimal global state needed to construct the objects.
  call displayVerbositySet  (verbosityLevelStandard)
  call eventsHooksInitialize(                      )
  call Functions_Global_Set (                      )
  ! Initialize the node-component hierarchy. The basic and dark matter profile components supply the mass, time, and scale radius
  ! of each halo, while the orbiting satellite component carries the position and velocity of the satellite relative to its host.
  parameters=inputParameters('testSuite/parameters/nodes/nodes_SIDM_Kummer.xml')
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("SIDM Kummer (2018) satellite evaporation and deceleration")
  ! Build the dependency stack. These must match the objects built from the parameter file by the dark matter profile component,
  ! so that the mass distributions seen by the evaporation object and by the nodes are the same.
  allocate(cosmologyParameters_  )
  allocate(cosmologyFunctions_   )
  allocate(virialDensityContrast_)
  allocate(darkMatterHaloScale_  )
  allocate(darkMatterProfileDMO_ )
  allocate(darkMatterParticleCDM_)
  allocate(particleConstant      )
  cosmologyParameters_  =cosmologyParametersSimple                          (OmegaMatter=0.2815d0,OmegaBaryon=0.0465d0,OmegaDarkEnergy=0.7185d0,temperatureCMB=2.78d0,HubbleConstant=70.0d0)
  cosmologyFunctions_   =cosmologyFunctionsMatterLambda                     (cosmologyParameters_=cosmologyParameters_)
  virialDensityContrast_=virialDensityContrastBryanNorman1998               (allowUnsupportedCosmology=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  darkMatterHaloScale_  =darkMatterHaloScaleVirialDensityContrastDefinition (cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,virialDensityContrast_=virialDensityContrast_)
  darkMatterProfileDMO_ =darkMatterProfileDMONFW                            (velocityDispersionUseSeriesExpansion=.false.,darkMatterHaloScale_=darkMatterHaloScale_)
  darkMatterParticleCDM_=darkMatterParticleCDM                              ()
  particleConstant      =darkMatterParticleSelfInteractingDarkMatterConstant(crossSectionSelfInteraction=crossSectionConstant,darkMatterParticle_=darkMatterParticleCDM_)
  evaporation           =satelliteEvaporationSIDMKummer2018                 (cosmologyParameters_,particleConstant,darkMatterHaloScale_,darkMatterProfileDMO_)
  deceleration          =satelliteDecelerationSIDMKummer2018                (cosmologyParameters_,particleConstant,darkMatterHaloScale_,darkMatterProfileDMO_)
  ! For an isotropic, constant cross section the tabulated integral over the range of scattering angles which result in net
  ! evaporation, ∫dσ/dθ dθ from π-θ_c to θ_c with θ_c=acos([x²-1]/[x²+1]), has the closed form (1-x²)/(1+x²) given by Kummer et
  ! al. (2018). Evaporation switches off for x≥1, where that range of angles is empty. The x values tested fall on the pinned
  ! lattice on which the factor is tabulated, so that this tests the tabulated integral itself and not the linear interpolation
  ! between lattice points.
  do i=1,size(xTabulated)
     factorTabulated(i)=evaporation%evaporationFactor(xTabulated(i),1.0d2)
     if (xTabulated(i) < 1.0d0) then
        factorAnalytic(i)=+(1.0d0-xTabulated(i)**2) &
             &            /(1.0d0+xTabulated(i)**2)
     else
        factorAnalytic(i)=+0.0d0
     end if
  end do
  call Assert("evaporation factor matches the analytic (1-x²)/(1+x²)",factorTabulated,factorAnalytic,absTol=1.0d-6)
  ! Construct a satellite orbiting within a host halo, and evaluate the mass loss rate. Since the only quantity in that rate
  ! which is not directly known here is the evaporation factor itself, it can be recovered from the rate, and used to pin the
  ! construction of x.
  nodeHost              => treeNode                                    (hostTree=tree    )
  node                  => treeNode                                    (hostTree=tree    )
  tree%nodeBase         => nodeHost
  nodeHost%firstChild   => node
  node    %parent       => nodeHost
  basicHost             => nodeHost%basic            (autoCreate=.true.)
  basic                 => node    %basic            (autoCreate=.true.)
  darkMatterProfileHost => nodeHost%darkMatterProfile(autoCreate=.true.)
  darkMatterProfile     => node    %darkMatterProfile(autoCreate=.true.)
  satellite             => node    %satellite        (autoCreate=.true.)
  call basicHost            %massSet (massHost                       )
  call basicHost            %timeSet (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic                %massSet (massSatellite                  )
  call basic                %timeSet (cosmologyFunctions_%cosmicTime(1.0d0))
  call darkMatterProfileHost%scaleSet(radiusScaleHost                )
  call darkMatterProfile    %scaleSet(radiusScaleSatellite           )
  call satellite            %boundMassSet(massSatellite              )
  call satellite            %positionSet ([radiusOrbital,0.0d0,0.0d0])
  ! Evaluate the factors, other than the evaporation factor itself, which enter the mass loss rate.
  fractionDarkMatter=+(                                    &
       &               +cosmologyParameters_%OmegaMatter() &
       &               -cosmologyParameters_%OmegaBaryon() &
       &              )                                    &
       &             /  cosmologyParameters_%OmegaMatter()
  massBoundary      =+massSatellite      &
       &             *fractionDarkMatter
  coordinates           =  [radiusOrbital,0.0d0,0.0d0]
  massDistributionHost_ => nodeHost             %massDistribution(           )
  densityHost           =  massDistributionHost_%density         (coordinates)
  !![
  <objectDestructor name="massDistributionHost_"/>
  !!]
  rateNormalization=+crossSectionConstant       & ! [cm² g⁻¹]
       &            *centi    **2/milli         & ! Convert cross-section from cm² g⁻¹ to m² kg⁻¹.
       &            *kilo                       & ! Convert velocity from km s⁻¹ to m s⁻¹.
       &            *massSolar   /megaParsec**3 & ! Convert density from M☉ Mpc⁻³ to kg m⁻³.
       &            *gigaYear                     ! Convert rate from s⁻¹ to Gyr⁻¹.
  ! Recover the evaporation factor from the mass loss rate at each of a range of orbital speeds.
  do i=1,size(speedOrbital)
     call satellite%velocitySet([0.0d0,speedOrbital(i),0.0d0])
     massLossRate         =evaporation %massLossRate(node)
     acceleration         =deceleration%acceleration(node)
     factorEvaporation (i)=-massLossRate         &
          &                /massBoundary         &
          &                /speedOrbital     (i) &
          &                /densityHost          &
          &                /rateNormalization
     ! Recover the product of the deceleration factor and the velocity dispersion correction factor. The acceleration returned
     ! is -v χ_d f_σ (ρ σ v), so the orbital speed appears in it twice, and must be divided out twice: once as the velocity
     ! vector which sets the direction of the drag, and once as the flux factor within the scattering rate. The velocity is
     ! directed along the y-axis, so only that component of the acceleration is non-zero.
     factorDeceleration(i)=-acceleration     (2) &
          &                /speedOrbital     (i) & ! The velocity which the drag acts along.
          &                /speedOrbital     (i) & ! The flux factor within the scattering rate.
          &                /densityHost          &
          &                /rateNormalization
  end do
  ! Evaporation must occur at every one of these speeds: the satellite here is far more weakly bound than the speed at which it
  ! moves through its host. Prior to the fix of the inverted definition of x, the constructed x was bounded below by unity and so
  ! the mass loss rate was identically zero.
  call Assert("mass loss rate is non-zero for a weakly bound satellite"           ,all(factorEvaporation      > 0.0d0                 ),.true.)
  ! The evaporation factor must increase with orbital speed: x=v_esc/v_rel decreases as the satellite moves more quickly, and
  ! χ_e=(1-x²)/(1+x²) decreases with x. Were x inverted, the ordering would be reversed.
  call Assert("evaporation factor increases with orbital speed"                   ,all(factorEvaporation(2:3) > factorEvaporation(1:2)),.true.)
  ! In the limit of an arbitrarily large orbital speed x→0 and hence χ_e→1: every scattering evaporates a particle. This fixes
  ! the normalization and the sense of x without reference to the details of either halo.
  call Assert("evaporation factor tends to unity at large orbital speed"          ,factorEvaporation(3) , 1.0d0,absTol=5.0d-3)
  ! The same limit fixes the deceleration factor. As x→0 the velocity dispersion correction factor of Appendix A of Kummer et
  ! al. (2018) also tends to unity, since it depends on the dispersion only through its ratio to the orbital speed, so the
  ! product recovered from the acceleration must tend to χ_d(0)=-1. For a weakly bound subhalo the scattered particles carry
  ! away more momentum than they bring in, and the sign of the drag reverses.
  call Assert("deceleration factor tends to -1 at large orbital speed"            ,factorDeceleration(3),-1.0d0,absTol=5.0d-3)
  ! Finally, a satellite bound deeply enough that its escape velocity exceeds the speed at which it moves relative to its host
  ! must suffer no evaporation at all: every scattering which ejects one of its own particles also captures a host particle. A
  ! second node is built for this, rather than the properties of the first being changed: mass distributions are memoized
  ! against the unique ID of the node, so a node modified in place can be seen with a stale mass distribution.
  nodeDense                  => treeNode                   (hostTree=tree    )
  nodeHost  %firstChild      => nodeDense
  nodeDense %sibling         => node
  nodeDense %parent          => nodeHost
  basicDense                 => nodeDense%basic            (autoCreate=.true.)
  darkMatterProfileDense     => nodeDense%darkMatterProfile(autoCreate=.true.)
  satelliteDense             => nodeDense%satellite        (autoCreate=.true.)
  call basicDense            %massSet     (massSatelliteDense                   )
  call basicDense            %timeSet     (cosmologyFunctions_%cosmicTime(1.0d0))
  call darkMatterProfileDense%scaleSet    (radiusScaleDense                     )
  call satelliteDense        %boundMassSet(massSatelliteDense                   )
  call satelliteDense        %positionSet ([radiusOrbital,0.0d0,0.0d0]          )
  call satelliteDense        %velocitySet ([0.0d0,speedOrbitalDense,0.0d0]      )
  massLossRateDense=evaporation%massLossRate(nodeDense)
  call Assert("mass loss rate is zero for a deeply bound, slowly moving satellite",massLossRateDense,0.0d0)
  ! In that same regime the drag must be a genuine deceleration: nothing escapes the subhalo, so momentum transfer from the host
  ! is complete and χ_d is positive.
  acceleration           =deceleration%acceleration(nodeDense)
  factorDecelerationDense=-acceleration     (2) &
       &                  /speedOrbitalDense    & ! The velocity which the drag acts along.
       &                  /speedOrbitalDense    & ! The flux factor within the scattering rate.
       &                  /densityHost          &
       &                  /rateNormalization
  call Assert("deceleration is a drag for a deeply bound, slowly moving satellite",factorDecelerationDense > 0.0d0,.true.)
  ! Clean up.
  call node     %destroy()
  call nodeDense%destroy()
  call nodeHost %destroy()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_SIDM_Kummer_Satellites
