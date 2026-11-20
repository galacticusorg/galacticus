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
Contains a program which tests the dissipative correction of the :galacticus-class:`mergerRemnantSizeCovington2008` merger remnant
size class.
!!}

program Test_Merger_Remnant_Size_Covington2008
  !!{RST
  Tests the dissipative correction, :math:`(1+f_\mathrm{gas}C_\mathrm{rad})`, of the :galacticus-class:`mergerRemnantSizeCovington2008`
  class. The form of this correction used, from :cite:t:`covington_effects_2011`, was calibrated for major mergers only, and so is
  applied only to those. The test compares remnant radii computed with :math:`C_\mathrm{rad}=2.75` to those computed with
  :math:`C_\mathrm{rad}=0` for gas-rich progenitors:

  * in a major merger (mass ratio 0.5), the remnant radius must be smaller by exactly the factor
    :math:`1/(1+f_\mathrm{gas}C_\mathrm{rad})`, where :math:`f_\mathrm{gas}` is the gas fraction of the remnant spheroid;
  * in a minor merger (mass ratio :math:`10^{-3}`), the remnant radius must not depend on :math:`C_\mathrm{rad}`, and must be close
    to the half mass radius of the host spheroid.

  No dark matter component is active, so the remnant size depends only on the progenitor masses and radii. A fresh pair of nodes is
  used for each calculation, since the classes memoize their results by node unique ID.
  !!}
  use :: Dark_Matter_Halo_Scales                , only : darkMatterHaloScaleClass
  use :: Display                                , only : displayVerbositySet               , verbosityLevelStandard
  use :: Events_Hooks                           , only : eventsHooksInitialize
  use :: Functions_Global_Utilities             , only : Functions_Global_Set
  use :: Galacticus_Nodes                       , only : nodeClassHierarchyInitialize      , nodeComponentBasic               , nodeComponentSpheroid              , treeNode
  use :: ISO_Varying_String                     , only : assignment(=)                     , varying_string
  use :: Input_Parameters                       , only : inputParameters
  use :: Node_Components                        , only : Node_Components_Initialize        , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Satellite_Merging_Mass_Movements       , only : destinationMergerSpheroid         , mergerMassMovementsSimple
  use :: Satellite_Merging_Progenitor_Properties, only : mergerProgenitorPropertiesStandard
  use :: Satellite_Merging_Remnant_Sizes        , only : mergerRemnantSizeCovington2008
  use :: Unit_Tests                             , only : Assert                            , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (varying_string                    )          :: parameterFile
  type            (inputParameters                   )          :: parameters
  class           (darkMatterHaloScaleClass          ), pointer :: darkMatterHaloScale_
  type            (mergerMassMovementsSimple         )          :: mergerMassMovements_
  type            (mergerProgenitorPropertiesStandard)          :: mergerProgenitorProperties_
  ! Radiative efficiency, and progenitor properties: spheroid stellar and gas masses [M☉] of the host, the stellar mass of the
  ! satellites in the major and minor mergers, the gas fraction of all progenitor spheroids, and spheroid scale radii [Mpc].
  double precision                                    , parameter :: efficiencyRadiative=2.75d0   , massStellarHost        =1.0d10, &
       &                                                             massSatelliteMajor =0.5d10   , massSatelliteMinor     =1.0d07, &
       &                                                             fractionGas        =0.2d0    , radiusHost             =1.0d-3, &
       &                                                             radiusSatellite    =5.0d-4
  ! The half mass radius of a Hernquist profile in units of its scale radius.
  double precision                                    , parameter :: radiusHalfMassHernquist=2.4142135624d0
  double precision                                                :: radiusDissipative           , radiusNonDissipative

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
  !![
  <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
  <referenceConstruct object="mergerMassMovements_"       >
   <constructor>
    mergerMassMovementsSimple(massRatioMajorMerger=0.25d0,destinationGasMinorMerger=destinationMergerSpheroid,destinationStarsMinorMerger=destinationMergerSpheroid)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="mergerProgenitorProperties_">
   <constructor>
    mergerProgenitorPropertiesStandard(mergerMassMovements_=mergerMassMovements_)
   </constructor>
  </referenceConstruct>
  !!]
  call Unit_Tests_Begin_Group("Covington et al. (2008) merger remnant sizes")
  ! Major merger: all mass moves to the remnant spheroid, so its gas fraction is that of the progenitors.
  radiusDissipative   =remnantRadius(efficiencyRadiative,massSatelliteMajor)
  radiusNonDissipative=remnantRadius(0.0d0              ,massSatelliteMajor)
  call Assert('major merger: radius reduced by 1/(1+f_gas C_rad)',radiusDissipative,radiusNonDissipative/(1.0d0+fractionGas*efficiencyRadiative),relTol=1.0d-10)
  ! Minor merger: the dissipative correction must not be applied.
  radiusDissipative   =remnantRadius(efficiencyRadiative,massSatelliteMinor)
  radiusNonDissipative=remnantRadius(0.0d0              ,massSatelliteMinor)
  call Assert('minor merger: radius independent of C_rad'        ,radiusDissipative,radiusNonDissipative                                        ,relTol=1.0d-10)
  call Assert('minor merger: radius close to that of the host'   ,radiusDissipative,radiusHalfMassHernquist*radiusHost                          ,relTol=1.0d-02)
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  !![
  <objectDestructor name="darkMatterHaloScale_"/>
  !!]
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()

contains

  double precision function remnantRadius(efficiencyRadiative,massSatellite)
    !!{RST
    Compute the remnant radius for a merger with a satellite of the given stellar mass, using the given radiative efficiency.
    !!}
    implicit none
    double precision                                , intent(in   ) :: efficiencyRadiative    , massSatellite
    type            (mergerRemnantSizeCovington2008)                :: mergerRemnantSize_
    type            (treeNode                      ), pointer       :: nodeHost               , nodeSatellite
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentSpheroid         ), pointer       :: spheroidHost           , spheroidSatellite
    double precision                                                :: velocityCircular       , angularMomentumSpecific

    mergerRemnantSize_=mergerRemnantSizeCovington2008(energyOrbital=1.0d0,efficiencyRadiative=efficiencyRadiative,darkMatterHaloScale_=darkMatterHaloScale_,mergerProgenitorProperties_=mergerProgenitorProperties_,mergerMassMovements_=mergerMassMovements_)
    nodeHost                  => treeNode              (                 )
    nodeSatellite             => treeNode              (                 )
    basic                     => nodeHost     %basic   (autoCreate=.true.)
    basic                     => nodeSatellite%basic   (autoCreate=.true.)
    spheroidHost              => nodeHost     %spheroid(autoCreate=.true.)
    spheroidSatellite         => nodeSatellite%spheroid(autoCreate=.true.)
    nodeSatellite%mergeTarget => nodeHost
    call spheroidHost     %massStellarSet    (massStellarHost                           )
    call spheroidHost     %massGasSet        (massStellarHost*fractionGas/(1.0d0-fractionGas))
    call spheroidHost     %radiusSet         (radiusHost                                )
    call spheroidHost     %angularMomentumSet(0.0d0                                     )
    call spheroidSatellite%massStellarSet    (massSatellite                             )
    call spheroidSatellite%massGasSet        (massSatellite  *fractionGas/(1.0d0-fractionGas))
    call spheroidSatellite%radiusSet         (radiusSatellite                           )
    call spheroidSatellite%angularMomentumSet(0.0d0                                     )
    call mergerRemnantSize_%get(nodeSatellite,remnantRadius,velocityCircular,angularMomentumSpecific)
    call nodeSatellite%destroy()
    call nodeHost     %destroy()
    deallocate(nodeSatellite)
    deallocate(nodeHost     )
    return
  end function remnantRadius

end program Test_Merger_Remnant_Size_Covington2008
