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
Contains a program which tests orbit calculations.
!!}

program Test_Orbits
  !!{RST
  Tests dark matter profiles.
  !!}
  use :: Calculations_Resets       , only : Calculations_Reset
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOIsothermal                                , darkMatterProfileDMONFW
  use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Display                   , only : displayVerbositySet                                           , verbosityLevelStandard
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize     , nodeComponentBasic                 , treeNode                                       , &
       &                                    nodeComponentDarkMatterProfile
  use :: Input_Parameters          , only : inputParameters
  use :: Node_Components           , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Kepler_Orbits             , only : keplerOrbit                                                   , keplerOrbitMasses                , keplerOrbitRadius
  use :: Unit_Tests                , only : Assert                                                        , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  use :: Satellite_Orbits          , only : Satellite_Orbit_Extremum_Phase_Space_Coordinates              , extremumPericenter               , extremumApocenter                  , Satellite_Orbit_Equivalent_Circular_Orbit_Radius
  use :: ISO_Varying_String        , only : var_str
  implicit none
  type            (treeNode                                                      ), pointer      :: node
  class           (nodeComponentBasic                                            ), pointer      :: basic
  class           (nodeComponentDarkMatterProfile                                ), pointer      :: darkMatterProfile
  double precision                                                                , parameter    :: massHost                =1.0d12      , massSatellite               =0.00d0, &
       &                                                                                            velocityFractionalRadial=0.9d00      , velocityFractionalTangential=0.75d0, &
       &                                                                                            fractionMassHostOriginal=2.0d00/3.0d0, concentration               =1.00d1
  type            (darkMatterProfileDMOIsothermal                                ), pointer      :: darkMatterProfileDMO_
  type            (darkMatterProfileDMONFW                                       ), pointer      :: darkMatterProfileDMONFW_
  type            (cosmologyParametersSimple                                     ), pointer      :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                ), pointer      :: cosmologyFunctions_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            ), pointer      :: darkMatterHaloScale_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer      :: virialDensityContrast_
  type            (keplerOrbit                                                   )               :: orbit
  type            (inputParameters                                               )               :: parameters
  double precision                                                                               :: radiusPericenter                     , velocityPericenter                 , &
       &                                                                                            radiusApocenter                      , velocityApocenter                  , &
       &                                                                                            radiusCircular                       , fractionRadiusVirialOriginal
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks and global functions.
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Initialize node components.
  parameters=inputParameters(var_str('testSuite/parameters/orbits.xml'))
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build required objects.
  allocate(cosmologyParameters_  )
  allocate(cosmologyFunctions_   )
  allocate(virialDensityContrast_)
  allocate(darkMatterHaloScale_  )
  allocate(darkMatterProfileDMO_ )
  allocate(darkMatterProfileDMONFW_)
  !![
  <referenceConstruct object="cosmologyParameters_"  >
   <constructor>
    cosmologyParametersSimple                                     (                                               &amp;
     &amp;                                                         OmegaMatter           = 0.2815d0             , &amp;
     &amp;                                                         OmegaBaryon           = 0.0000d0             , &amp;
     &amp;                                                         OmegaDarkEnergy       = 0.7185d0             , &amp;
     &amp;                                                         temperatureCMB        = 2.7800d0             , &amp;
     &amp;                                                         HubbleConstant        =69.3000d0               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"   >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                               &amp;
     &amp;                                                         cosmologyParameters_  =cosmologyParameters_    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_">
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                               &amp;
     &amp;                                                         tableStore            =.true.                , &amp;
     &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_     &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"  >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                               &amp;
     &amp;                                                         cosmologyParameters_  =cosmologyParameters_  , &amp;
     &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_   , &amp;
     &amp;                                                         virialDensityContrast_=virialDensityContrast_  &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMO_" >
   <constructor>
    darkMatterProfileDMOIsothermal                                (                                               &amp;
     &amp;                                                         darkMatterHaloScale_  =darkMatterHaloScale_    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMONFW_">
   <constructor>
    darkMatterProfileDMONFW                                       (                                                 &amp;
     &amp;                                                         velocityDispersionUseSeriesExpansion=.false.   , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_ &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  ! Create a node.
  node              => treeNode                  (                 )
  ! Create components.
  basic             => node    %basic            (autoCreate=.true.)
  darkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
  ! Set node properties.
  call basic%timeSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic%massSet(massHost                           )
  ! Create an orbit.
  orbit=keplerOrbit()
  call orbit%massesSet(                                          &
       &                                    massSatellite      , &
       &                                    massHost             &
       &              )
  call orbit%radiusSet(                                          &
       &               darkMatterHaloScale_%radiusVirial (node)  &
       &              )
  ! Begin unit tests.
  call Unit_Tests_Begin_Group('Orbits')
  ! Non-circular orbit. Solutions for peri- and apo-centric radii found using Mathematica.
  call orbit%reset                (keep=[keplerOrbitMasses,keplerOrbitRadius]                            )
  call orbit%velocityRadialSet    (darkMatterHaloScale_%velocityVirial(node)*velocityFractionalRadial    )
  call orbit%velocityTangentialSet(darkMatterHaloScale_%velocityVirial(node)*velocityFractionalTangential)
  call Satellite_Orbit_Extremum_Phase_Space_Coordinates(node,orbit,extremumPericenter,radiusPericenter,velocityPericenter,darkMatterHaloScale_)
  call Satellite_Orbit_Extremum_Phase_Space_Coordinates(node,orbit,extremumApocenter ,radiusApocenter ,velocityApocenter ,darkMatterHaloScale_)
  call Assert('non-circular orbit, pericenter radius',radiusPericenter,0.428095d0*darkMatterHaloScale_%radiusVirial(node),relTol=1.0d-5)
  call Assert('non-circular orbit,  apocenter radius',radiusApocenter ,1.825500d0*darkMatterHaloScale_%radiusVirial(node),relTol=1.0d-5)
  ! Circular orbit.
  call orbit%reset                (keep=[keplerOrbitMasses,keplerOrbitRadius]                            )
  call orbit%velocityRadialSet    (                     0.0d0                                            )
  call orbit%velocityTangentialSet(darkMatterHaloScale_%velocityVirial(node)                             )
  call Satellite_Orbit_Extremum_Phase_Space_Coordinates(node,orbit,extremumPericenter,radiusPericenter,velocityPericenter,darkMatterHaloScale_)
  call Satellite_Orbit_Extremum_Phase_Space_Coordinates(node,orbit,extremumApocenter ,radiusApocenter ,velocityApocenter ,darkMatterHaloScale_)
  call Assert('    circular orbit, pericenter radius',radiusPericenter,1.000000d0*darkMatterHaloScale_%radiusVirial(node),relTol=1.0d-5)
  call Assert('    circular orbit,  apocenter radius',radiusApocenter ,1.000000d0*darkMatterHaloScale_%radiusVirial(node),relTol=1.0d-5)
  ! Equivalent circular orbit of the non-circular orbit. For the isothermal halo, Φ(r)=V²[ln(r/r_vir)-1], so this is at
  ! r=r_vir exp(e/V²+1/2), where e=(v_r²+v_t²)/2-V² is the orbital energy.
  call orbit%reset                (keep=[keplerOrbitMasses,keplerOrbitRadius]                            )
  call orbit%velocityRadialSet    (darkMatterHaloScale_%velocityVirial(node)*velocityFractionalRadial    )
  call orbit%velocityTangentialSet(darkMatterHaloScale_%velocityVirial(node)*velocityFractionalTangential)
  radiusCircular=Satellite_Orbit_Equivalent_Circular_Orbit_Radius(node,orbit,darkMatterHaloScale_,darkMatterProfileDMO_)
  call Assert('non-circular orbit, equivalent circular orbit radius',radiusCircular,1.204723d0*darkMatterHaloScale_%radiusVirial(node),relTol=1.0d-5)
  call Unit_Tests_End_Group()
  ! Orbit defined in a host of lower mass at the same epoch (e.g. a host which has since grown, or a satellite which has moved to
  ! a more massive host). The orbit is defined at the virial radius of that original host, which is smaller by a factor
  ! f^{1/3}, where f is the ratio of original to current host mass, with velocities in units of its virial velocity (also
  ! smaller by f^{1/3}). Converting to the potential of the current host gives the energy e=(v_r²+v_t²)/2+Φ(r_vir,0), with
  ! Φ(r_vir)=-V_vir² for the current host. Solutions found by root finding (scipy), in units of the current virial radius.
  call Unit_Tests_Begin_Group('Orbits in a more massive host')
  fractionRadiusVirialOriginal=fractionMassHostOriginal**(1.0d0/3.0d0)
  call orbit%reset                (                                                                                                    )
  call orbit%massesSet            (                                          massSatellite                                            ,  &
       &                                                                     massHost                    *fractionMassHostOriginal       &
       &                          )
  call orbit%radiusSet            (darkMatterHaloScale_%radiusVirial  (node)                             *fractionRadiusVirialOriginal )
  call orbit%velocityRadialSet    (darkMatterHaloScale_%velocityVirial(node)*velocityFractionalRadial    *fractionRadiusVirialOriginal )
  call orbit%velocityTangentialSet(darkMatterHaloScale_%velocityVirial(node)*velocityFractionalTangential*fractionRadiusVirialOriginal )
  call Satellite_Orbit_Extremum_Phase_Space_Coordinates(node,orbit,extremumPericenter,radiusPericenter,velocityPericenter,darkMatterHaloScale_)
  call Satellite_Orbit_Extremum_Phase_Space_Coordinates(node,orbit,extremumApocenter ,radiusApocenter ,velocityApocenter ,darkMatterHaloScale_)
  radiusCircular=Satellite_Orbit_Equivalent_Circular_Orbit_Radius(node,orbit,darkMatterHaloScale_,darkMatterProfileDMO_)
  call Assert('non-circular orbit, pericenter radius'                 ,radiusPericenter,0.331139d0*darkMatterHaloScale_%radiusVirial(node),relTol=1.0d-5)
  call Assert('non-circular orbit,  apocenter radius'                 ,radiusApocenter ,1.347652d0*darkMatterHaloScale_%radiusVirial(node),relTol=1.0d-5)
  call Assert('non-circular orbit, equivalent circular orbit radius'  ,radiusCircular  ,0.894538d0*darkMatterHaloScale_%radiusVirial(node),relTol=1.0d-5)
  call Unit_Tests_End_Group()
  ! Equivalent circular orbit found in a dark matter only profile (here NFW, with concentration c) which differs from the
  ! (isothermal) mass distribution of the host. The radius solves Φ(r)+V_c²(r)/2=e, with e=(v_r²+v_t²)/2-V² as above, and
  ! Φ(r)=-V²{[ln(1+cr/r_vir)r_vir/r-ln(1+c)]/μ(c)+1}, where μ(c)=ln(1+c)-c/(1+c). Solution found by root finding (mpmath), in
  ! units of the virial radius. It differs from that in the isothermal mass distribution of the host (1.204723), so detects use
  ! of the wrong profile.
  call Unit_Tests_Begin_Group('Equivalent circular orbit in the dark matter only profile')
  call darkMatterProfile%scaleSet(darkMatterHaloScale_%radiusVirial(node)/concentration)
  call orbit%reset                (                                                                      )
  call orbit%massesSet            (                                          massSatellite             , &
       &                                                                     massHost                    &
       &                          )
  call orbit%radiusSet            (darkMatterHaloScale_%radiusVirial  (node)                             )
  call orbit%velocityRadialSet    (darkMatterHaloScale_%velocityVirial(node)*velocityFractionalRadial    )
  call orbit%velocityTangentialSet(darkMatterHaloScale_%velocityVirial(node)*velocityFractionalTangential)
  radiusCircular=Satellite_Orbit_Equivalent_Circular_Orbit_Radius(node,orbit,darkMatterHaloScale_,darkMatterProfileDMONFW_)
  call Assert('non-circular orbit, equivalent circular orbit radius',radiusCircular,1.294712d0*darkMatterHaloScale_%radiusVirial(node),relTol=1.0d-5)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Uninitialize node components.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  ! Clean up objects.
  !![
  <objectDestructor name="cosmologyParameters_"  />
  <objectDestructor name="cosmologyFunctions_"   />
  <objectDestructor name="virialDensityContrast_"/>
  <objectDestructor name="darkMatterHaloScale_"  />
  <objectDestructor name="darkMatterProfileDMO_" />
  <objectDestructor name="darkMatterProfileDMONFW_"/>
  !!]
end program Test_Orbits
