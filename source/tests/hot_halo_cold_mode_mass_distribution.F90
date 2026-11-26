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
Contains a program which tests the mass distribution of the ``coldMode`` hot halo component.
!!}

!+    Contributions to this file made by: Claude.

program Test_Hot_Halo_Cold_Mode_Mass_Distribution
  !!{RST
  Tests the mass distribution returned by the ``coldMode`` hot halo component. When both the hot-mode and the cold-mode
  reservoirs hold gas, a request for the gaseous mass distribution of all components must return the composite of the two, so
  the mass enclosed within the outer radius (at which both distributions are truncated) must be the sum of the hot-mode and
  cold-mode masses. The mass enclosed at a smaller radius (where the two distributions differ, as their core radii differ) must
  be the sum of the masses enclosed by the hot-mode and cold-mode distributions requested separately. Requests for the
  hot-mode, or the cold-mode, component alone must enclose that component's mass only.
  !!}
  use :: Abundances_Structure      , only : zeroAbundances
  use :: Cosmology_Functions       , only : cosmologyFunctions                 , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale                , darkMatterHaloScaleClass
  use :: Display                   , only : displayVerbositySet                , verbosityLevelStandard
  use :: Error                     , only : Error_Handler_Register
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galactic_Structure_Options, only : componentTypeColdHalo              , componentTypeHotHalo             , &
       &                                    massTypeGaseous
  use :: Galacticus_Nodes          , only : nodeClassHierarchyInitialize       , nodeComponentBasic               , &
       &                                    nodeComponentHotHalo               , treeNode
  use :: IO_HDF5                   , only : ioHDF5AccessInitialize
  use :: ISO_Varying_String        , only : assignment(=)                      , varying_string
  use :: Input_Parameters          , only : inputParameters
  use :: Mass_Distributions        , only : massDistributionClass
  use :: Node_Components           , only : Node_Components_Initialize         , Node_Components_Thread_Initialize, &
       &                                    Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                , only : Assert                             , Unit_Tests_Begin_Group           , &
       &                                    Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  class           (cosmologyFunctionsClass ), pointer   :: cosmologyFunctions_
  class           (darkMatterHaloScaleClass), pointer   :: darkMatterHaloScale_
  type            (treeNode                ), pointer   :: node
  class           (nodeComponentBasic      ), pointer   :: basic
  class           (nodeComponentHotHalo    ), pointer   :: hotHalo
  class           (massDistributionClass   ), pointer   :: massDistributionAll             , massDistributionHotMode       , &
       &                                                   massDistributionColdMode
  type            (varying_string          )            :: parameterFile
  type            (inputParameters         )            :: parameters
  ! Properties of the test halo.
  double precision                          , parameter :: massHalo                =1.0d12, massHotMode             =1.0d11, &
       &                                                   massColdMode            =5.0d10
  double precision                                      :: time                           , radiusOuter                    , &
       &                                                   massAllOuter                   , massHotModeOuter               , &
       &                                                   massColdModeOuter              , massAllInner                   , &
       &                                                   massHotModeInner               , massColdModeInner

  ! Establish error handlers, so that errors reported by the GSL are trapped rather than aborting.
  call Error_Handler_Register()
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize the lock used to serialize access to the HDF5 library.
  call ioHDF5AccessInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cold mode hot halo mass distribution")
  ! Read in controlling parameters and initialize the node component hierarchy.
  parameterFile='testSuite/parameters/hotHaloColdModeMassDistribution.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
  <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
  !!]
  time=cosmologyFunctions_%cosmicTime(1.0d0)
  ! Build an isolated halo with both hot-mode and cold-mode gas in its hot halo, extending to half of the virial radius.
  node  => treeNode()
  basic => node%basic(autoCreate=.true.)
  call basic%massSet            (massHalo)
  call basic%timeSet            (time    )
  call basic%timeLastIsolatedSet(time    )
  hotHalo     => node%hotHalo(autoCreate=.true.)
  radiusOuter =  0.5d0*darkMatterHaloScale_%radiusVirial(node)
  call hotHalo%massSet          (massHotMode   )
  call hotHalo%massColdSet      (massColdMode  )
  call hotHalo%abundancesSet    (zeroAbundances)
  call hotHalo%abundancesColdSet(zeroAbundances)
  call hotHalo%outerRadiusSet   (radiusOuter   )
  ! Get the mass distributions directly from the component (bypassing the node's memoized distributions).
  massDistributionAll      => hotHalo%massDistribution(                                    massType=massTypeGaseous)
  massDistributionHotMode  => hotHalo%massDistribution(componentType=componentTypeHotHalo ,massType=massTypeGaseous)
  massDistributionColdMode => hotHalo%massDistribution(componentType=componentTypeColdHalo,massType=massTypeGaseous)
  massAllOuter     =massDistributionAll     %massEnclosedBySphere(      radiusOuter)
  massHotModeOuter =massDistributionHotMode %massEnclosedBySphere(      radiusOuter)
  massColdModeOuter=massDistributionColdMode%massEnclosedBySphere(      radiusOuter)
  massAllInner     =massDistributionAll     %massEnclosedBySphere(0.2d0*radiusOuter)
  massHotModeInner =massDistributionHotMode %massEnclosedBySphere(0.2d0*radiusOuter)
  massColdModeInner=massDistributionColdMode%massEnclosedBySphere(0.2d0*radiusOuter)
  call Assert("hot-mode mass within the outer radius"            ,massHotModeOuter ,massHotMode                       ,relTol=1.0d-6)
  call Assert("cold-mode mass within the outer radius"           ,massColdModeOuter,massColdMode                      ,relTol=1.0d-6)
  call Assert("total gaseous mass within the outer radius"       ,massAllOuter     ,massHotMode     +massColdMode     ,relTol=1.0d-6)
  call Assert("total gaseous mass within 0.2 of the outer radius",massAllInner     ,massHotModeInner+massColdModeInner,relTol=1.0d-6)
  !![
  <objectDestructor name="massDistributionAll"     />
  <objectDestructor name="massDistributionHotMode" />
  <objectDestructor name="massDistributionColdMode"/>
  !!]
  call node%destroy()
  deallocate(node)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  !![
  <objectDestructor name="cosmologyFunctions_" />
  <objectDestructor name="darkMatterHaloScale_"/>
  !!]
end program Test_Hot_Halo_Cold_Mode_Mass_Distribution
