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

!!{
Contains a program to test calculations for tidal track dark matter profiles.
!!}

program Test_Dark_Matter_Profiles_Tidal_Tracks
  !!{
  Test calculations for tidal track dark matter profiles.
  !!}
  use :: Calculations_Resets       , only : Calculations_Reset
  use :: Coordinates               , only : coordinateSpherical                                           , assignment(=)
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOPenarrubia2010                            , darkMatterProfileDMONFW
  use :: Display                   , only : displayMessage                                                , displayVerbositySet              , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize     , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
          &                                 treeNode                                                      , nodeComponentSatellite
  use :: Input_Parameters          , only : inputParameters
  use :: Node_Components           , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Mass_Distributions        , only : massDistributionClass
  use :: Unit_Tests                , only : Assert                                                        , Skip                             , Unit_Tests_Begin_Group             , Unit_Tests_End_Group          , &
          &                                 Unit_Tests_Finish
  implicit none
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            ), pointer   :: darkMatterHaloScale_
  type            (cosmologyParametersSimple                                     ), pointer   :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                ), pointer   :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer   :: virialDensityContrast_
  type            (darkMatterProfileDMONFW                                       ), pointer   :: darkMatterProfileNFW_
  type            (darkMatterProfileDMOPenarrubia2010                            ), pointer   :: darkMatterProfilePenarrubia2010_
  type            (treeNode                                                      ), pointer   :: node_
  class           (nodeComponentBasic                                            ), pointer   :: basic_
  class           (nodeComponentDarkMatterProfile                                ), pointer   :: darkMatterProfile_
  class           (nodeComponentSatellite                                        ), pointer   :: satellite_
  class           (massDistributionClass                                         ), pointer   :: massDistributionNFW_                   , massDistributionPenarrubia2010_
  double precision                                                                , parameter :: concentration                   =+8.0d0, massVirial                     =+1.0d12, &
       &                                                                                         muRadius                        =-0.3d0, etaRadius                      =+0.4d-0, &
       &                                                                                         muVelocity                      =+0.4d0, etaVelocity                    =+0.3d-0
  type            (inputParameters                                               )            :: parameters
  double precision                                                                            :: radiusScale                            , radiusVirial                           , &
       &                                                                                         velocityMaximumInitial                 , radiusMaximumInitial                   , &
       &                                                                                         velocityMaximumTidalTrack              , radiusMaximumTidalTrack                , &
       &                                                                                         fractionMassBound
  type            (coordinateSpherical                                           )            :: coordinatesScale                       , coordinatesVirial

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Tidal track dark matter profiles")
  parameters=inputParameters('testSuite/parameters/darkMatterProfilesTidalTracks.xml')
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  node_                            =>  treeNode                              (                 )
  basic_                           =>  node_               %basic            (autoCreate=.true.)
  darkMatterProfile_               =>  node_               %darkMatterProfile(autoCreate=.true.)
  satellite_                       =>  node_               %satellite        (autoCreate=.true.)
  allocate(darkMatterHaloScale_            )
  allocate(cosmologyParameters_            )
  allocate(cosmologyFunctions_             )
  allocate(virialDensityContrast_          )
  allocate(darkMatterProfileNFW_           )
  allocate(darkMatterProfilePenarrubia2010_)
  !![
  <referenceConstruct object="cosmologyParameters_"            >
   <constructor>
    cosmologyParametersSimple                                     (                                                             &amp;
     &amp;                                                         OmegaMatter                         = 0.30d0               , &amp;
     &amp;                                                         OmegaBaryon                         = 0.00d0               , &amp;
     &amp;                                                         OmegaDarkEnergy                     = 0.70d0               , &amp;
     &amp;                                                         temperatureCMB                      = 2.78d0               , &amp;
     &amp;                                                         HubbleConstant                      =70.00d0                 &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"             >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                                             &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_"          >
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                             &amp;
     &amp;                                                         tableStore                          =.true.,                 &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_     &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"            >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                                             &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_  , &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_   , &amp;
     &amp;                                                         virialDensityContrast_              =virialDensityContrast_  &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileNFW_"           >
   <constructor>
    darkMatterProfileDMONFW                                       (                                                             &amp;
     &amp;                                                         velocityDispersionUseSeriesExpansion=.false.               , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfilePenarrubia2010_">
   <constructor>
    darkMatterProfileDMOPenarrubia2010                            (                                                             &amp;
     &amp;                                                         alpha                               =+1.0d0                , &amp;
     &amp;                                                         beta                                =+3.0d0                , &amp;
     &amp;                                                         betaStripped                        =+5.0d0                , &amp;
     &amp;                                                         gamma                               =+1.0d0                , &amp;
     &amp;                                                         muRadius                            =muRadius              , &amp;
     &amp;                                                         etaRadius                           =etaRadius             , &amp;
     &amp;                                                         muVelocity                          =muVelocity            , &amp;
     &amp;                                                         etaVelocity                         =etaVelocity           , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  call basic_    %timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %massSet            (massVirial                           )
  call satellite_%boundMassSet       (massVirial                           )
  radiusVirial =+darkMatterHaloScale_%radiusVirial(node_)
  radiusScale  =+radiusVirial &
       &        /concentration
  coordinatesScale =[radiusScale ,0.0d0,0.0d0]
  coordinatesVirial=[radiusVirial,0.0d0,0.0d0]
  call darkMatterProfile_%scaleSet(radiusScale)
  call Calculations_Reset(node_)
  ! Store the initial values of rmax and Vmax.
  massDistributionNFW_            => darkMatterProfileNFW_           %get                         (node_)
  massDistributionPenarrubia2010_ => darkMatterProfilePenarrubia2010_%get                         (node_)
  radiusMaximumInitial            =  massDistributionNFW_            %  radiusRotationCurveMaximum(     )
  velocityMaximumInitial          =  massDistributionNFW_            %velocityRotationCurveMaximum(     )
  ! Begin tests.
  call Unit_Tests_Begin_Group("Unstripped profile matches NFW"    )
  call Assert("Density at scale radius" ,massDistributionNFW_%density(coordinates=coordinatesScale ),massDistributionPenarrubia2010_%density(coordinates=coordinatesScale ),relTol=1.0d-6)
  call Assert("Density at virial radius",massDistributionNFW_%density(coordinates=coordinatesVirial),massDistributionPenarrubia2010_%density(coordinates=coordinatesVirial),relTol=1.0d-6)
  call Unit_Tests_End_Group  (                                    )
  !![
  <objectDestructor name="massDistributionNFW_"           />
  <objectDestructor name="massDistributionPenarrubia2010_"/>
  !!]
  call Unit_Tests_Begin_Group("Stripped profile (95%) tidal track")
  call satellite_%boundMassSet(0.95d0*massVirial)
  call Calculations_Reset(node_)
  massDistributionNFW_            => darkMatterProfileNFW_           %get                         (node_)
  massDistributionPenarrubia2010_ => darkMatterProfilePenarrubia2010_%get                         (node_)
  fractionMassBound        =+satellite_%boundMass() &
       &                    /basic_    %     mass()
  velocityMaximumTidalTrack=velocityMaximumInitial*2.0d0**muVelocity*fractionMassBound**etaVelocity/(1.0d0+fractionMassBound)**muVelocity
  radiusMaximumTidalTrack  =  radiusMaximumInitial*2.0d0**muRadius  *fractionMassBound**etaRadius  /(1.0d0+fractionMassBound)**muRadius
  call Assert("Tidal track rmax",radiusMaximumTidalTrack  ,massDistributionPenarrubia2010_%  radiusRotationCurveMaximum(),relTol=1.0d-6)
  call Assert("Tidal track Vmax",velocityMaximumTidalTrack,massDistributionPenarrubia2010_%velocityRotationCurveMaximum(),relTol=1.0d-6)
  call Unit_Tests_End_Group               ()
  !![
  <objectDestructor name="massDistributionNFW_"           />
  <objectDestructor name="massDistributionPenarrubia2010_"/>
  !!]
  call Unit_Tests_Begin_Group("Stripped profile (30%) tidal track")
  call satellite_%boundMassSet(0.30d0*massVirial)
  call Calculations_Reset(node_)
  massDistributionNFW_            => darkMatterProfileNFW_           %get                         (node_)
  massDistributionPenarrubia2010_ => darkMatterProfilePenarrubia2010_%get                         (node_)
  fractionMassBound        =+satellite_%boundMass() &
       &                    /basic_    %     mass()
  velocityMaximumTidalTrack=velocityMaximumInitial*2.0d0**muVelocity*fractionMassBound**etaVelocity/(1.0d0+fractionMassBound)**muVelocity
  radiusMaximumTidalTrack  =  radiusMaximumInitial*2.0d0**muRadius  *fractionMassBound**etaRadius  /(1.0d0+fractionMassBound)**muRadius
  call Assert("Tidal track rmax",radiusMaximumTidalTrack  ,massDistributionPenarrubia2010_%  radiusRotationCurveMaximum(),relTol=1.0d-6)
  call Assert("Tidal track Vmax",velocityMaximumTidalTrack,massDistributionPenarrubia2010_%velocityRotationCurveMaximum(),relTol=1.0d-6)
  call Unit_Tests_End_Group               ()
  !![
  <objectDestructor name="massDistributionNFW_"           />
  <objectDestructor name="massDistributionPenarrubia2010_"/>
  !!]
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Clean up.
  call node_%destroy()
  deallocate(node_)
  !![
  <objectDestructor name="darkMatterHaloScale_"            />
  <objectDestructor name="cosmologyParameters_"            />
  <objectDestructor name="cosmologyFunctions_"             />
  <objectDestructor name="virialDensityContrast_"          />
  <objectDestructor name="darkMatterProfileNFW_"           />
  <objectDestructor name="darkMatterProfilePenarrubia2010_"/>
  !!]
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Tidal_Tracks
