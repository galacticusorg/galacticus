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
Contains a program which tests the limit on the energy radiated by the hot halo in the :cite:t:`benson_galaxy_2010-1` cooling model.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Cooling_Energy_Radiated_Cap
  !!{RST
  Tests the limit applied by the :galacticus-class:`nodeOperatorCoolingEnergyRadiated` node operator to the energy radiated by the
  hot halo, :math:`\mathcal{E}`, in the :cite:t:`benson_galaxy_2010-1` cooling model. The time available for cooling (from
  :galacticus-class:`coolingTimeAvailableBensonBower2010`) is :math:`\mathcal{E}/\dot{\mathcal{E}}`, and the limit is such that this
  can not exceed the cooling time of the hot halo gas at the outer radius of the hot halo.

  For a set of halos, each given a hot atmosphere containing one tenth of its mass with an outer radius equal to the virial
  radius, the radiated energy is first set to a value far in excess of the limit, and the operator's post-step limiter
  applied. The test checks that the limiter reports that it changed the node state (requiring the ODE solver to recompute
  derivatives), and that the resulting time available for cooling equals the cooling time at the outer radius, evaluated directly
  from the :galacticus-class:`coolingTime` class using the density and temperature of the hot halo at that radius. The radiated
  energy is then set to half of its limited value, and the test checks that the limiter leaves it, and the node state, unchanged.
  !!}
  use :: Abundances_Structure         , only : zeroAbundances                         , abundances
  use :: Chemical_Abundances_Structure, only : chemicalAbundances
  use :: Cooling_Times                , only : coolingTime                            , coolingTimeClass
  use :: Cooling_Times_Available      , only : coolingTimeAvailable                   , coolingTimeAvailableClass
  use :: Coordinates                  , only : coordinateSpherical                    , assignment(=)
  use :: Cosmology_Functions          , only : cosmologyFunctions                     , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales      , only : darkMatterHaloScale                    , darkMatterHaloScaleClass
  use :: Display                      , only : displayVerbositySet                    , verbosityLevelStandard
  use :: Error                        , only : Error_Handler_Register
  use :: Events_Hooks                 , only : eventsHooksInitialize
  use :: Functions_Global_Utilities   , only : Functions_Global_Set
  use :: Galactic_Structure_Options   , only : componentTypeHotHalo                   , massTypeGaseous
  use :: Galacticus_Nodes             , only : nodeClassHierarchyInitialize           , nodeComponentBasic               , &
       &                                       nodeComponentHotHalo                   , treeNode
  use :: IO_HDF5                      , only : ioHDF5AccessInitialize
  use :: ISO_Varying_String           , only : assignment(=)                          , varying_string
  use :: Input_Parameters             , only : inputParameters
  use :: Interface_GSL                , only : GSL_Success                            , GSL_Failure
  use :: Mass_Distributions           , only : massDistributionClass                  , kinematicsDistributionClass
  use :: Node_Components              , only : Node_Components_Initialize             , Node_Components_Thread_Initialize, &
       &                                       Node_Components_Thread_Uninitialize    , Node_Components_Uninitialize
  use :: Nodes_Operators              , only : nodeOperator                           , nodeOperatorClass
  use :: Radiation_Fields             , only : radiationFieldCosmicMicrowaveBackground
  use :: Unit_Tests                   , only : Assert                                 , Unit_Tests_Begin_Group           , &
       &                                       Unit_Tests_End_Group                   , Unit_Tests_Finish
  implicit none
  class           (cosmologyFunctionsClass                ), pointer               :: cosmologyFunctions_
  class           (darkMatterHaloScaleClass               ), pointer               :: darkMatterHaloScale_
  class           (coolingTimeClass                       ), pointer               :: coolingTime_
  class           (coolingTimeAvailableClass              ), pointer               :: coolingTimeAvailable_
  class           (nodeOperatorClass                      ), pointer               :: nodeOperator_
  class           (massDistributionClass                  ), pointer               :: massDistribution_
  class           (kinematicsDistributionClass            ), pointer               :: kinematicsDistribution_
  type            (radiationFieldCosmicMicrowaveBackground), pointer               :: radiation
  type            (treeNode                               ), pointer               :: node
  class           (nodeComponentBasic                     ), pointer               :: basic
  class           (nodeComponentHotHalo                   ), pointer               :: hotHalo
  type            (coordinateSpherical                    )                        :: coordinates
  type            (varying_string                         )                        :: parameterFile
  type            (inputParameters                        )                        :: parameters
  type            (abundances                             )                        :: abundancesGas
  type            (chemicalAbundances                     )                        :: densityChemicals
  ! The fraction of the halo mass placed in the hot atmosphere.
  double precision                                         , parameter             :: massHotFraction        =1.0d-1
  ! A radiated energy (in units of ergs s⁻¹ Gyr) far in excess of the limit for any halo considered here (for which it is of order
  ! 10⁴³ in these units).
  double precision                                         , parameter             :: energyRadiatedExcessive=1.0d60
  ! Test cases: halo masses and redshifts.
  integer                                                  , parameter             :: countCases             =4
  double precision                                         , dimension(countCases) :: massHalo               =[1.0d11,1.0d12,1.0d12,1.0d13] , &
       &                                                                              redshift               =[0.0d00,0.0d00,1.0d00,0.0d00]
  double precision                                         , dimension(countCases) :: timeAvailableLimited          , timeCoolingOuter      , &
       &                                                                              energyRadiatedUnlimited       , energyRadiatedHalf
  logical                                                  , dimension(countCases) :: limitedStateChanged           , unlimitedStateChanged
  integer                                                                          :: energyRadiatedID              , status                , &
       &                                                                              i
  double precision                                                                 :: time                          , radiusOuter           , &
       &                                                                              density                       , temperature

  ! Establish error handlers, so that errors reported by the GSL are trapped rather than aborting.
  call Error_Handler_Register()
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize the lock used to serialize access to the HDF5 library (needed by the tabulated cooling function).
  call ioHDF5AccessInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cooling energy radiated limit")
  ! Read in controlling parameters and initialize the node component hierarchy.
  parameterFile='testSuite/parameters/coolingEnergyRadiatedCap.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
  <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
  <objectBuilder class="coolingTime"          name="coolingTime_"          source="parameters"/>
  <objectBuilder class="nodeOperator"         name="nodeOperator_"         source="parameters"/>
  <objectBuilder class="coolingTimeAvailable" name="coolingTimeAvailable_" source="parameters"/>
  <addMetaProperty component="hotHalo" name="energyRadiatedBensonBower2010" id="energyRadiatedID" isEvolvable="yes"/>
  !!]
  allocate(radiation)
  !![
  <referenceConstruct object="radiation" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
  !!]
  ! Build a halo with a hot atmosphere at each test mass and redshift.
  do i=1,countCases
     time    =  cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift(i)))
     node    => treeNode                      (                                                            )
     basic   => node               %basic     (autoCreate=.true.                                           )
     hotHalo => node               %hotHalo   (autoCreate=.true.                                           )
     call basic  %massSet            (                massHalo(i)   )
     call basic  %timeSet            (                time          )
     call basic  %timeLastIsolatedSet(                time          )
     radiusOuter=darkMatterHaloScale_%radiusVirial(node)
     call hotHalo%massSet            (massHotFraction*massHalo(i)   )
     call hotHalo%outerRadiusSet     (                radiusOuter   )
     call hotHalo%abundancesSet      (                zeroAbundances)
     ! Find the cooling time of the hot halo gas at the outer radius directly.
     coordinates             =  [radiusOuter,0.0d0,0.0d0]
     massDistribution_       => node                   %massDistribution      (componentTypeHotHalo,massTypeGaseous)
     kinematicsDistribution_ => massDistribution_      %kinematicsDistribution(                                    )
     density                 =  massDistribution_      %density               (coordinates                         )
     temperature             =  kinematicsDistribution_%temperature           (coordinates                         )
     !![
     <objectDestructor name="massDistribution_"      />
     <objectDestructor name="kinematicsDistribution_"/>
     !!]
     abundancesGas=hotHalo%abundances()
     call abundancesGas%massToMassFraction(hotHalo%mass())
     call radiation%timeSet(time)
     timeCoolingOuter(i)=coolingTime_%time(node,temperature,density,abundancesGas,densityChemicals,radiation)
     ! Set a radiated energy far in excess of the limit, and apply the limiter.
     call hotHalo%floatRank0MetaPropertySet(energyRadiatedID,energyRadiatedExcessive)
     status=GSL_Success
     call nodeOperator_%differentialEvolutionPostStep(node,status)
     limitedStateChanged    (i)=status == GSL_Failure
     timeAvailableLimited   (i)=coolingTimeAvailable_%timeAvailable(node)
     ! Set a radiated energy below the limit, and apply the limiter.
     energyRadiatedHalf     (i)=0.5d0*hotHalo%floatRank0MetaPropertyGet(energyRadiatedID)
     call hotHalo%floatRank0MetaPropertySet(energyRadiatedID,energyRadiatedHalf(i))
     status=GSL_Success
     call nodeOperator_%differentialEvolutionPostStep(node,status)
     unlimitedStateChanged  (i)=status /= GSL_Success
     energyRadiatedUnlimited(i)=hotHalo%floatRank0MetaPropertyGet(energyRadiatedID)
     call node%destroy()
     deallocate(node)
  end do
  call Assert("limiter changes the node state"          ,all(limitedStateChanged)  ,.true.                           )
  call Assert("limited time available = t_cool(r_outer)",timeAvailableLimited      ,timeCoolingOuter  ,relTol=1.0d-06)
  call Assert("no change to the node state below limit" ,any(unlimitedStateChanged),.false.                          )
  call Assert("radiated energy unchanged below limit"   ,energyRadiatedUnlimited   ,energyRadiatedHalf,relTol=1.0d-12)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  !![
  <objectDestructor name="cosmologyFunctions_"  />
  <objectDestructor name="darkMatterHaloScale_" />
  <objectDestructor name="coolingTime_"         />
  <objectDestructor name="coolingTimeAvailable_"/>
  <objectDestructor name="nodeOperator_"        />
  <objectDestructor name="radiation"            />
  !!]
end program Test_Cooling_Energy_Radiated_Cap
