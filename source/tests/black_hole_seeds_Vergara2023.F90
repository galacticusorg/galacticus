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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a program which tests the :galacticus-class:`blackHoleSeedsVergara2023` black hole seeds class.
!!}

program Tests_Black_Hole_Seeds_Vergara2023
  !!{RST
  Tests the :galacticus-class:`blackHoleSeedsVergara2023` black hole seeds class. Its ``mass`` method is called during
  differential evolution, so must not change the nuclear star cluster: any change made there is overwritten by the evolution
  (or, for the collapse flag, persists even if the trial state at which it was made is never reached). The changes that follow
  from the collapse (marking the cluster as collapsed, recording its properties, and removing the mass of the seed from its
  stars, with their abundances and luminosities) must instead be made by ``seedFormed``, which is called once the seed has been placed.
  !!}
  use :: Abundances_Structure      , only : abundances                  , Abundances_Property_Count   , operator(*)
  use :: Black_Hole_Seeds          , only : blackHoleSeedsClass
  use :: Display                   , only : displayVerbositySet         , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes          , only : nodeClassHierarchyInitialize, nodeComponentBasic          , nodeComponentNSC                 , treeNode                           , &
       &                                    nodeClassHierarchyFinalize  , defaultNSCComponent
  use :: Input_Paths               , only : inputPath                   , pathTypeExec
  use :: Input_Parameters          , only : inputParameters
  use :: ISO_Varying_String        , only : assignment(=)               , operator(//)                , varying_string                   , var_str
  use :: Node_Components           , only : Node_Components_Initialize  , Node_Components_Uninitialize, Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
  use :: Stellar_Luminosities_Structure, only : stellarLuminosities
  use :: Unit_Tests                , only : Assert                      , Unit_Tests_Begin_Group      , Unit_Tests_End_Group             , Unit_Tests_Finish
  implicit none
  type            (treeNode           ), pointer                   :: node
  class           (nodeComponentBasic ), pointer                   :: basic
  class           (nodeComponentNSC   ), pointer                   :: nuclearStarCluster
  class           (blackHoleSeedsClass), pointer                   :: blackHoleSeeds_
  type            (inputParameters    )                            :: parameters
  type            (varying_string     )                            :: parameterFile
  type            (abundances         )                            :: abundancesStellar
  type            (stellarLuminosities)                            :: luminositiesStellar
  ! Mass [M☉] of the nuclear star cluster, the time [Gyr] at which it is tested, and the mean time [Gyr] at which its stars formed
  ! (giving an age of 1 Gyr). With the radius efficiency of the parameter file the critical mass for collapse is about 10⁶ M☉, so
  ! the cluster collapses.
  double precision                     , parameter                 :: massStellar        =1.0d7, time                   =10.0d0, &
       &                                                              timeFormation      =9.0d0
  ! The seed mass efficiency, as set in the parameter file.
  double precision                     , parameter                 :: massEfficiency     =0.1d0
  double precision                     , allocatable, dimension(:) :: arrayAbundances           , arrayLuminosities
  double precision                                                 :: massSeed                 , massSeedRepeat
  integer                                                          :: stellarMassFormedID      , timeStellarMassFormedID       , &
       &                                                              blackHoleSeedMassID

  ! Initialize.
  parameterFile=inputPath(pathTypeExec)//'testSuite/parameters/blackHoleSeedsVergara2023.xml'
  parameters   =inputParameters(parameterFile)
  call displayVerbositySet              (verbosityLevelStandard)
  call eventsHooksInitialize            (                      )
  call Functions_Global_Set             (                      )
  call nodeClassHierarchyInitialize     (parameters            )
  call Node_Components_Initialize       (parameters            )
  call Node_Components_Thread_Initialize(parameters            )
  ! Create the meta-properties holding the star formation history of the nuclear star cluster (in a model these are created by
  ! the "ages" node operator), and the one holding the seed mass formed (created by the black hole seeds class).
  stellarMassFormedID    =defaultNSCComponent%addFloatRank0MetaProperty(var_str('agesStellarMassFormed'    ),'NSC:agesStellarMassFormed'    ,isCreator=.true.,isEvolvable=.true.)
  timeStellarMassFormedID=defaultNSCComponent%addFloatRank0MetaProperty(var_str('agesTimeStellarMassFormed'),'NSC:agesTimeStellarMassFormed',isCreator=.true.,isEvolvable=.true.)
  !![
  <objectBuilder class="blackHoleSeeds" name="blackHoleSeeds_" source="parameters"/>
  !!]
  blackHoleSeedMassID    =defaultNSCComponent%addFloatRank0MetaProperty(var_str('blackHoleSeedMassFormed'  ),'NSC:blackHoleSeedMassFormed'  ,isCreator=.false.                 )
  ! Build a node with a nuclear star cluster.
  node               => treeNode      (                 )
  basic              => node    %basic(autoCreate=.true.)
  nuclearStarCluster => node    %NSC  (autoCreate=.true.)
  call basic%timeSet            (time)
  call basic%timeLastIsolatedSet(time)
  allocate(arrayAbundances(Abundances_Property_Count()))
  arrayAbundances=massStellar
  call abundancesStellar%deserialize(arrayAbundances)
  allocate(arrayLuminosities(luminositiesStellar%serializeCount()))
  arrayLuminosities=massStellar
  call luminositiesStellar%deserialize(arrayLuminosities)
  call nuclearStarCluster%massStellarSet           (massStellar                                    )
  call nuclearStarCluster%abundancesStellarSet     (abundancesStellar                              )
  call nuclearStarCluster%luminositiesStellarSet   (luminositiesStellar                            )
  call nuclearStarCluster%floatRank0MetaPropertySet(stellarMassFormedID    ,massStellar              )
  call nuclearStarCluster%floatRank0MetaPropertySet(timeStellarMassFormedID,massStellar*timeFormation)
  call Unit_Tests_Begin_Group("Vergara et al. (2023) black hole seeds")
  ! Compute the seed mass. This must not change the nuclear star cluster.
  massSeed      =blackHoleSeeds_%mass(node)
  massSeedRepeat=blackHoleSeeds_%mass(node)
  call Assert('nuclear star cluster collapses'                                ,massSeed > 0.0d0                ,.true.                                   )
  call Assert('seed mass is the mass efficiency times the stellar mass'       ,massSeed                        ,massEfficiency*massStellar,relTol=1.0d-12)
  call Assert('computing the seed mass does not mark the cluster as collapsed',nuclearStarCluster%isCollapsed(),.false.                                  )
  call Assert('computing the seed mass does not change the stellar mass'      ,nuclearStarCluster%massStellar(),massStellar               ,relTol=1.0d-12)
  call Assert('computing the seed mass again gives the same mass'             ,massSeedRepeat                  ,massSeed                  ,relTol=1.0d-12)
  ! Place the seed. The changes that follow from the collapse must now be made.
  call blackHoleSeeds_%seedFormed(node)
  abundancesStellar=nuclearStarCluster%abundancesStellar()
  call abundancesStellar%serialize(arrayAbundances)
  luminositiesStellar=nuclearStarCluster%luminositiesStellar()
  call luminositiesStellar%serialize(arrayLuminosities)
  call Assert('placing the seed marks the cluster as collapsed'              ,nuclearStarCluster%isCollapsed()                                 ,.true.                                                                           )
  call Assert('placing the seed removes its mass from the stars'             ,nuclearStarCluster%massStellar()                                 ,       (1.0d0-massEfficiency)*massStellar                         ,relTol=1.0d-12)
  call Assert('placing the seed removes its share of the stellar abundances' ,arrayAbundances                                                  ,spread((1.0d0-massEfficiency)*massStellar,1,size(arrayAbundances)),relTol=1.0d-12)
  call Assert('one stellar luminosity is tracked'                           ,size(arrayLuminosities)                                         ,1                                                                                )
  call Assert('placing the seed removes its share of the stellar luminosities',arrayLuminosities                                              ,spread((1.0d0-massEfficiency)*massStellar,1,size(arrayLuminosities)),relTol=1.0d-12)
  call Assert('placing the seed records the seed mass'                       ,nuclearStarCluster%floatRank0MetaPropertyGet(blackHoleSeedMassID),massSeed                                                          ,relTol=1.0d-12)
  call Assert('no further seed forms once the cluster has collapsed'         ,blackHoleSeeds_%mass(node)                                       ,0.0d0                                                                            )
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Clean up.
  call node%destroy()
  deallocate(node)
  !![
  <objectDestructor name="blackHoleSeeds_"/>
  !!]
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Tests_Black_Hole_Seeds_Vergara2023
