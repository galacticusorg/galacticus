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
Contains a program which tests the chemicals accreted in each mode by the ``coldMode`` halo accretion class.
!!}

program Tests_Accretion_Halo_Cold_Mode_Chemicals
  !!{RST
  Tests the masses (and rates) of chemicals accreted in each mode by the :galacticus-class:`accretionHaloColdMode` class. The
  chemicals accreted in the hot and cold modes must sum to those accreted in total, and the chemicals in each mode must scale
  with the mass accreted in that mode. The halo mass is chosen so that both modes contribute substantially.
  !!}
  use :: Accretion_Halos              , only : accretionHaloClass             , accretionModeTotal               , accretionModeHot                   , accretionModeCold
  use :: Chemical_Abundances_Structure, only : chemicalAbundances             , Chemicals_Property_Count
  use :: Cosmology_Functions          , only : cosmologyFunctionsClass
  use :: Display                      , only : displayVerbositySet            , verbosityLevelStandard
  use :: Events_Hooks                 , only : eventsHooksInitialize
  use :: Functions_Global_Utilities   , only : Functions_Global_Set
  use :: Galacticus_Nodes             , only : nodeClassHierarchyInitialize   , nodeComponentBasic               , treeNode                           , nodeClassHierarchyFinalize
  use :: Input_Paths                  , only : inputPath                      , pathTypeExec
  use :: Input_Parameters             , only : inputParameters
  use :: ISO_Varying_String           , only : assignment(=)                  , operator(//)                     , varying_string
  use :: Node_Components              , only : Node_Components_Initialize     , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                   , only : Assert                         , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode               ), pointer                   :: node
  class           (nodeComponentBasic     ), pointer                   :: basic
  class           (accretionHaloClass     ), pointer                   :: accretionHalo_
  class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_
  type            (chemicalAbundances     )                            :: chemicalsTotal       , chemicalsHot        , &
       &                                                                  chemicalsCold
  double precision                         , allocatable, dimension(:) :: arrayTotal           , arrayHot            , &
       &                                                                  arrayCold
  ! Halo mass [M☉] and accretion rate [M☉/Gyr], chosen to lie in the transition between cold and hot mode accretion.
  double precision                         , parameter                 :: massHalo      =3.0d11, rateAccretionHalo=3.0d10
  double precision                                                     :: massTotal            , massHot             , &
       &                                                                  massCold             , rateTotal           , &
       &                                                                  rateHot              , rateCold
  type            (inputParameters        )                            :: parameters
  type            (varying_string         )                            :: parameterFile

  ! Initialize.
  parameterFile=inputPath(pathTypeExec)//'testSuite/parameters/accretionHaloColdModeChemicals.xml'
  parameters   =inputParameters(parameterFile)
  call displayVerbositySet              (verbosityLevelStandard)
  call eventsHooksInitialize            (                      )
  call Functions_Global_Set             (                      )
  call nodeClassHierarchyInitialize     (parameters            )
  call Node_Components_Initialize       (parameters            )
  call Node_Components_Thread_Initialize(parameters            )
  !![
  <objectBuilder class="accretionHalo"      name="accretionHalo_"      source="parameters"/>
  <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
  !!]
  ! Build an isolated node at z=1.
  node  => treeNode      (                 )
  basic => node    %basic(autoCreate=.true.)
  call basic%timeSet            (cosmologyFunctions_%cosmicTime(0.5d0))
  call basic%timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(0.5d0))
  call basic%massSet            (massHalo                             )
  call basic%accretionRateSet   (rateAccretionHalo                    )
  allocate(arrayTotal(Chemicals_Property_Count()))
  allocate(arrayHot  (Chemicals_Property_Count()))
  allocate(arrayCold (Chemicals_Property_Count()))
  call Unit_Tests_Begin_Group("Cold mode accretion: chemicals")
  ! Masses.
  massTotal     =accretionHalo_%accretedMass         (node,accretionModeTotal)
  massHot       =accretionHalo_%accretedMass         (node,accretionModeHot  )
  massCold      =accretionHalo_%accretedMass         (node,accretionModeCold )
  chemicalsTotal=accretionHalo_%accretedMassChemicals(node,accretionModeTotal)
  chemicalsHot  =accretionHalo_%accretedMassChemicals(node,accretionModeHot  )
  chemicalsCold =accretionHalo_%accretedMassChemicals(node,accretionModeCold )
  call chemicalsTotal%serialize(arrayTotal)
  call chemicalsHot  %serialize(arrayHot  )
  call chemicalsCold %serialize(arrayCold )
  call Assert('both modes contribute to the accreted mass'          ,min(massHot,massCold) > 0.05d0*massTotal,.true.                          )
  call Assert('accreted masses in the two modes sum to the total'    ,massHot+massCold                        ,massTotal          ,relTol=1.0d-9)
  call Assert('accreted chemicals in the two modes sum to the total' ,arrayHot+arrayCold                      ,arrayTotal         ,relTol=1.0d-9)
  ! The ionization state of each mode depends on its temperature, but the mass of hydrogen (neutral plus ionized) scales with
  ! the mass accreted in that mode. (Only to the extent that the mass of an electron is negligible, as the mass of ionized
  ! hydrogen excludes it: the two modes differ in ionization, so these agree to about 10⁻⁶.)
  call Assert('hot mode hydrogen mass scales with the hot mode mass'  ,arrayHot (1)+arrayHot (2),(arrayTotal(1)+arrayTotal(2))*massHot /massTotal,relTol=1.0d-5)
  call Assert('cold mode hydrogen mass scales with the cold mode mass',arrayCold(1)+arrayCold(2),(arrayTotal(1)+arrayTotal(2))*massCold/massTotal,relTol=1.0d-5)
  ! Rates.
  rateTotal     =accretionHalo_%accretionRate         (node,accretionModeTotal)
  rateHot       =accretionHalo_%accretionRate         (node,accretionModeHot  )
  rateCold      =accretionHalo_%accretionRate         (node,accretionModeCold )
  chemicalsTotal=accretionHalo_%accretionRateChemicals(node,accretionModeTotal)
  chemicalsHot  =accretionHalo_%accretionRateChemicals(node,accretionModeHot  )
  chemicalsCold =accretionHalo_%accretionRateChemicals(node,accretionModeCold )
  call chemicalsTotal%serialize(arrayTotal)
  call chemicalsHot  %serialize(arrayHot  )
  call chemicalsCold %serialize(arrayCold )
  call Assert('accretion rates in the two modes sum to the total'          ,rateHot+rateCold  ,rateTotal ,relTol=1.0d-9)
  call Assert('chemical accretion rates in the two modes sum to the total' ,arrayHot+arrayCold,arrayTotal,relTol=1.0d-9)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Clean up.
  call node%destroy()
  deallocate(node)
  !![
  <objectDestructor name="accretionHalo_"     />
  <objectDestructor name="cosmologyFunctions_"/>
  !!]
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Tests_Accretion_Halo_Cold_Mode_Chemicals
