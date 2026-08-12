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
  Contains a program to test the limiting of satellite tidal radii by the radius of a solitonic core.
  !!}

program Test_Satellite_Tidal_Stripping_Radius_Soliton_Limited
  !!{RST
  Test that the :galacticus-class:`satelliteTidalStrippingRadiusSolitonLimited` class limits the tidal radius returned by the model
  which it wraps to be no smaller than the radius of the solitonic core.

  The wrapped model is :galacticus-class:`satelliteTidalStrippingRadiusKing1962`, applied to an isolated node which is not a
  satellite. For such a node that class returns the virial radius immediately, giving a value which this test can compute
  independently, so that the assertions below test only the limiting behavior and not the tidal radius calculation itself.
  !!}
  use :: Cosmology_Functions            , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters           , only : cosmologyParametersSimple                                     , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales        , only : darkMatterHaloScaleVirialDensityContrastDefinition            , darkMatterHaloScaleClass
  use :: Display                        , only : displayVerbositySet                                           , verbosityLevelStandard
  use :: Events_Hooks                   , only : eventsHooksInitialize
  use :: Functions_Global_Utilities     , only : Functions_Global_Set
  use :: Galacticus_Nodes               , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize        , &
       &                                         nodeComponentBasic                                            , nodeComponentDarkMatterProfile      , &
       &                                         treeNode                                                      , defaultDarkMatterProfileComponent
  use :: Input_Parameters               , only : inputParameters
  use :: ISO_Varying_String             , only : var_str
  use :: Node_Components                , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize   , &
       &                                         Node_Components_Thread_Uninitialize                           , Node_Components_Uninitialize
  use :: Satellite_Tidal_Stripping_Radii, only : satelliteTidalStrippingRadiusKing1962                         , satelliteTidalStrippingRadiusSolitonLimited
  use :: Satellites_Tidal_Fields        , only : satelliteTidalFieldNull                                       , satelliteTidalFieldClass
  use :: Unit_Tests                     , only : Assert                                                        , Unit_Tests_Begin_Group              , &
       &                                         Unit_Tests_End_Group                                          , Unit_Tests_Finish
  use :: Virial_Density_Contrast        , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  implicit none
  type            (treeNode                                                      ), pointer   :: node
  class           (nodeComponentBasic                                            ), pointer   :: basic
  class           (nodeComponentDarkMatterProfile                                ), pointer   :: darkMatterProfile
  type            (cosmologyParametersSimple                                     ), pointer   :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                ), pointer   :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer   :: virialDensityContrast_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            ), pointer   :: darkMatterHaloScale_
  type            (satelliteTidalFieldNull                                       ), pointer   :: satelliteTidalField_
  type            (satelliteTidalStrippingRadiusKing1962                         ), pointer   :: radiusKing1962_
  type            (satelliteTidalStrippingRadiusSolitonLimited                          ), pointer   :: radiusLimited_
  type            (inputParameters                                               )            :: parameters
  double precision                                                                            :: radiusVirial                 , radiusSoliton
  integer                                                                                     :: radiusSolitonID
  double precision                                                                , parameter :: massVirial            =1.0d12, timeNode     =13.8d0

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Satellite tidal stripping radius: soliton-limited")
  parameters=inputParameters('testSuite/parameters/satellites/tidalStrippingRadiusSolitonLimited.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Register the soliton radius meta-property as a creator. The class under test registers it as a non-creator, so without a
  ! creator the meta-property has no storage and can not be written. This must be done *before* the dark matter profile
  ! component is created below, as the meta-property storage is allocated at component creation.
  radiusSolitonID=defaultDarkMatterProfileComponent%addFloatRank0MetaProperty(                                                   &
       &                                                                      var_str('solitonRadiusSoliton'                  ), &
       &                                                                              'darkMatterProfile:solitonRadiusSoliton' , &
       &                                                                      isCreator  =.true.                               , &
       &                                                                      isEvolvable=.false.                                &
       &                                                                     )
  ! Build the objects required by the tidal radius models.
  allocate(cosmologyParameters_  )
  allocate(cosmologyFunctions_   )
  allocate(virialDensityContrast_)
  allocate(darkMatterHaloScale_  )
  allocate(satelliteTidalField_  )
  allocate(radiusKing1962_       )
  allocate(radiusLimited_        )
  cosmologyParameters_  = cosmologyParametersSimple                                     (OmegaMatter=0.3153d0,OmegaBaryon=0.0493d0,OmegaDarkEnergy=0.6847d0,temperatureCMB=2.72548d0,HubbleConstant=67.36d0)
  cosmologyFunctions_   = cosmologyFunctionsMatterLambda                                (cosmologyParameters_=cosmologyParameters_)
  virialDensityContrast_= virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(tableStore=.true.,cosmologyFunctions_=cosmologyFunctions_)
  darkMatterHaloScale_  = darkMatterHaloScaleVirialDensityContrastDefinition            (cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,virialDensityContrast_=virialDensityContrast_)
  satelliteTidalField_  = satelliteTidalFieldNull                                       ()
  radiusKing1962_       = satelliteTidalStrippingRadiusKing1962                         (efficiencyCentrifugal=1.0d0,applyPreInfall=.false.,cosmologyParameters_=cosmologyParameters_,darkMatterHaloScale_=darkMatterHaloScale_,satelliteTidalField_=satelliteTidalField_)
  radiusLimited_        = satelliteTidalStrippingRadiusSolitonLimited                          (satelliteTidalStrippingRadius_=radiusKing1962_)
  ! Build a node. It is not a satellite, so the wrapped King (1962) model returns the virial radius.
  node              => treeNode                  (                 )
  basic             => node    %basic            (autoCreate=.true.)
  darkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
  call basic%massSet(massVirial)
  call basic%timeSet(timeNode  )
  radiusVirial=darkMatterHaloScale_%radiusVirial(node)
  ! With no solitonic core the tidal radius must be returned unchanged. The soliton profiles store a negative radius in this
  ! case.
  call darkMatterProfile%floatRank0MetaPropertySet(radiusSolitonID,-1.0d0)
  call Assert('no soliton: tidal radius returned unchanged' ,radiusLimited_%radius(node),radiusVirial ,relTol=1.0d-6)
  ! With a solitonic core smaller than the tidal radius the tidal radius must again be returned unchanged.
  radiusSoliton=0.5d0*radiusVirial
  call darkMatterProfile%floatRank0MetaPropertySet(radiusSolitonID,radiusSoliton)
  call Assert('soliton smaller than tidal radius: unchanged',radiusLimited_%radius(node),radiusVirial ,relTol=1.0d-6)
  ! With a solitonic core larger than the tidal radius the tidal radius must be limited to the soliton radius, so that stripping
  ! of the outer halo can not remove material from within the core.
  radiusSoliton=2.0d0*radiusVirial
  call darkMatterProfile%floatRank0MetaPropertySet(radiusSolitonID,radiusSoliton)
  call Assert('soliton larger than tidal radius: limited'   ,radiusLimited_%radius(node),radiusSoliton,relTol=1.0d-6)
  ! Clean up.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
end program Test_Satellite_Tidal_Stripping_Radius_Soliton_Limited
