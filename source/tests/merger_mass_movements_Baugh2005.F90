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
Contains a program which tests the burst criterion of the :galacticus-class:`mergerMassMovementsBaugh2005` merger mass movements
class.
!!}

program Test_Merger_Mass_Movements_Baugh2005
  !!{RST
  Tests the burst criterion of the :galacticus-class:`mergerMassMovementsBaugh2005` class in minor mergers. Following
  :cite:t:`baugh_can_2005`, a minor merger triggers a burst if the satellite to central galaxy mass ratio exceeds
  :math:`f_\mathrm{burst}` and the gas fraction of the central galaxy exceeds :math:`f_\mathrm{gas,crit}`. The test uses two
  mergers for which this criterion and the bulge-to-total criterion (selectable through ``[burstCriterion]``) disagree:

  * a satellite with 10% of the central's mass, merging with a gas-rich central with a bulge-to-total ratio of 0.1: a burst under
    the mass ratio criterion, but not under the bulge-to-total criterion;
  * a satellite with 1% of the central's mass, merging with a gas-rich central with a bulge-to-total ratio of 0.002: no burst under
    the mass ratio criterion, but a burst under the bulge-to-total criterion.

  In a burst the gas of the central moves to its spheroid; otherwise it does not move. A fresh pair of nodes is used for each
  configuration, since the class memoizes its results by node unique ID.
  !!}
  use :: Display                         , only : displayVerbositySet           , verbosityLevelStandard
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Functions_Global_Utilities      , only : Functions_Global_Set
  use :: Galacticus_Nodes                , only : nodeClassHierarchyInitialize  , nodeComponentBasic               , nodeComponentDisk                  , nodeComponentSpheroid       , &
       &                                          treeNode
  use :: ISO_Varying_String              , only : assignment(=)                 , varying_string
  use :: Input_Parameters                , only : inputParameters
  use :: Node_Components                 , only : Node_Components_Initialize    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Satellite_Merging_Mass_Movements, only : destinationMergerSpheroid     , destinationMergerUnmoved         , enumerationDestinationMergerType   , mergerMassMovementsBaugh2005, &
       &                                          baugh2005BurstCriterionMassRatio, baugh2005BurstCriterionBulgeToTotal, enumerationBaugh2005BurstCriterionType
  use :: Unit_Tests                      , only : Assert                        , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type   (varying_string                        )               :: parameterFile
  type   (inputParameters                       )               :: parameters
  type   (enumerationBaugh2005BurstCriterionType), dimension(2) :: criteria             =[baugh2005BurstCriterionMassRatio,baugh2005BurstCriterionBulgeToTotal]
  character(len=14                              ), dimension(2) :: criteriaLabel        =['mass ratio    ','bulge-to-total']
  ! Expected outcomes (burst or not) for each criterion (first index) and configuration (second index).
  logical                                        , dimension(2,2) :: burstExpected      =reshape([.true.,.false.,.false.,.true.],[2,2])
  integer                                                       :: i                    , j
  logical                                                       :: burst

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Read in controlling parameters and initialize the node component hierarchy.
  parameterFile='testSuite/parameters/diskSpheroidComponents.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Baugh et al. (2005) merger burst criterion")
  do i=1,2
     do j=1,2
        burst=triggersBurst(criteria(i),j)
        call Assert(                                                                           &
             &      trim(criteriaLabel(i))//' criterion: configuration '//char(ichar('0')+j), &
             &      burst                                                                    , &
             &      burstExpected(i,j)                                                         &
             &     )
     end do
  end do
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()

contains

  logical function triggersBurst(burstCriterion,configuration)
    !!{RST
    Build a satellite and central galaxy in the given ``configuration``, and return whether their merger triggers a burst under
    the given ``burstCriterion``.
    !!}
    implicit none
    type            (enumerationBaugh2005BurstCriterionType), intent(in   ) :: burstCriterion
    integer                                                 , intent(in   ) :: configuration
    type            (mergerMassMovementsBaugh2005          )                :: mergerMassMovements_
    type            (treeNode                              ), pointer       :: nodeHost               , nodeSatellite
    class           (nodeComponentBasic                    ), pointer       :: basic
    class           (nodeComponentDisk                     ), pointer       :: diskHost
    class           (nodeComponentSpheroid                 ), pointer       :: spheroidHost           , spheroidSatellite
    type            (enumerationDestinationMergerType      )                :: destinationGasSatellite, destinationStarsSatellite, &
         &                                                                     destinationGasHost     , destinationStarsHost
    logical                                                                 :: mergerIsMajor
    ! Masses [M☉] of the central's disk gas, disk stars, and spheroid stars, and of the satellite's spheroid stars, for each
    ! configuration. Both centrals have a gas fraction (0.8 and 0.9) above f_gas,crit=0.75.
    double precision                                        , dimension(2)  :: massGasDiskHost        =[8.0d9,9.0d9], massStellarDiskHost=[1.0d9,9.8d8], &
         &                                                                     massStellarSpheroidHost=[1.0d9,2.0d7], massSatellite      =[1.0d9,1.0d8]

    ! Use the defaults of Baugh et al. (2005) for f_burst and f_gas,crit, and the default major merger mass ratio, so that both
    ! mergers are minor.
    mergerMassMovements_=mergerMassMovementsBaugh2005(massRatioMajorMerger=0.25d0,destinationGasMinorMerger=destinationMergerSpheroid,ratioMassBurst=0.05d0,fractionGasCriticalBurst=0.75d0,burstCriterion=burstCriterion)
    nodeHost                  => treeNode              (                 )
    nodeSatellite             => treeNode              (                 )
    basic                     => nodeHost     %basic   (autoCreate=.true.)
    basic                     => nodeSatellite%basic   (autoCreate=.true.)
    diskHost                  => nodeHost     %disk    (autoCreate=.true.)
    spheroidHost              => nodeHost     %spheroid(autoCreate=.true.)
    spheroidSatellite         => nodeSatellite%spheroid(autoCreate=.true.)
    nodeSatellite%mergeTarget => nodeHost
    call diskHost         %massGasSet    (massGasDiskHost        (configuration))
    call diskHost         %massStellarSet(massStellarDiskHost    (configuration))
    call diskHost         %radiusSet     (3.0d-3                                )
    call spheroidHost     %massStellarSet(massStellarSpheroidHost(configuration))
    call spheroidHost     %radiusSet     (1.0d-3                                )
    call spheroidSatellite%massStellarSet(massSatellite          (configuration))
    call spheroidSatellite%radiusSet     (5.0d-4                                )
    call mergerMassMovements_%get(nodeSatellite,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    call Assert('merger is minor',mergerIsMajor,.false.)
    ! In a burst the gas of the central moves to its spheroid, while its stars stay in place. Otherwise neither moves.
    triggersBurst=destinationGasHost == destinationMergerSpheroid
    call Assert('central stars do not move',destinationStarsHost == destinationMergerUnmoved,.true.)
    call nodeSatellite%destroy()
    call nodeHost     %destroy()
    deallocate(nodeSatellite)
    deallocate(nodeHost     )
    return
  end function triggersBurst

end program Test_Merger_Mass_Movements_Baugh2005
