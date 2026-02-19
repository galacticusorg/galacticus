!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Contains a class which implements a task to pre-build tabulations needed for broadband luminosity calculations.
  !!}

  use :: Stellar_Population_Broad_Band_Luminosities, only : stellarPopulationBroadBandLuminositiesClass
  use :: Stellar_Populations                       , only : stellarPopulationClass

  !![
  <task name="taskBuildBroadbandLuminosityTabulations">
   <description>A task which pre-builds tabulations needed for broadband luminosity calculations.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskBuildBroadbandLuminosityTabulations
     !!{
     Implementation of a task which pre-builds tabulations needed for broadband luminosity calculations.
     !!}
     private
     class  (stellarPopulationBroadBandLuminositiesClass), pointer :: stellarPopulationBroadBandLuminosities_ => null()
     class  (stellarPopulationClass                     ), pointer :: stellarPopulation_                      => null()
     logical                                                       :: nodeComponentsInitialized               =  .false.
   contains
     final     ::                       buildBroadbandLuminosityTabulationsDestructor
     procedure :: perform            => buildBroadbandLuminosityTabulationsPerform
     procedure :: requiresOutputFile => buildBroadbandLuminosityTabulationsRequiresOutputFile
  end type taskBuildBroadbandLuminosityTabulations

  interface taskBuildBroadbandLuminosityTabulations
     !!{
     Constructors for the \refClass{taskBuildBroadbandLuminosityTabulations} task.
     !!}
     module procedure buildBroadbandLuminosityTabulationsParameters
     module procedure buildBroadbandLuminosityTabulationsInternal
  end interface taskBuildBroadbandLuminosityTabulations

contains

  function buildBroadbandLuminosityTabulationsParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskBuildBroadbandLuminosityTabulations} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Node_Components , only : Node_Components_Initialize
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    implicit none
    type (taskBuildBroadbandLuminosityTabulations    )                         :: self
    type (inputParameters                            ), intent(inout), target  :: parameters
    class(stellarPopulationBroadBandLuminositiesClass)               , pointer :: stellarPopulationBroadBandLuminosities_
    class(stellarPopulationClass                     )               , pointer :: stellarPopulation_
    type (inputParameters                            )               , pointer :: parametersRoot

    ! Ensure the nodes objects are initialized.
    parametersRoot => parameters
    do while (associated(parametersRoot%parent))
       parametersRoot => parametersRoot%parent
    end do
    call nodeClassHierarchyInitialize(parametersRoot)
    call Node_Components_Initialize  (parametersRoot)
    self%nodeComponentsInitialized=.true.
    !![
    <objectBuilder class="stellarPopulationBroadBandLuminosities" name="stellarPopulationBroadBandLuminosities_" source="parameters"/>
    <objectBuilder class="stellarPopulation"                      name="stellarPopulation_"                      source="parameters"/>
    !!]
    self=taskBuildBroadbandLuminosityTabulations(stellarPopulationBroadBandLuminosities_,stellarPopulation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarPopulationBroadBandLuminosities_"/>
    <objectDestructor name="stellarPopulation_"                     />
    !!]
    return
  end function buildBroadbandLuminosityTabulationsParameters

  function buildBroadbandLuminosityTabulationsInternal(stellarPopulationBroadBandLuminosities_,stellarPopulation_) result(self)
    !!{
    Constructor for the \refClass{taskBuildBroadbandLuminosityTabulations} task class which takes a parameter set as input.
    !!}
    implicit none
    type (taskBuildBroadbandLuminosityTabulations    )                        :: self
    class(stellarPopulationBroadBandLuminositiesClass), intent(in   ), target :: stellarPopulationBroadBandLuminosities_
    class(stellarPopulationClass                     ), intent(in   ), target :: stellarPopulation_
    !![
    <constructorAssign variables="*stellarPopulationBroadBandLuminosities_, *stellarPopulation_"/>
    !!]

    return
  end function buildBroadbandLuminosityTabulationsInternal

  subroutine buildBroadbandLuminosityTabulationsDestructor(self)
    !!{
    Destructor for the \refClass{taskBuildBroadbandLuminosityTabulations} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskBuildBroadbandLuminosityTabulations), intent(inout) :: self
    
    !![
    <objectDestructor name="self%stellarPopulationBroadBandLuminosities_"/>
    <objectDestructor name="self%stellarPopulation_"                     />
    !!]
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine buildBroadbandLuminosityTabulationsDestructor

  subroutine buildBroadbandLuminosityTabulationsPerform(self,status)
    !!{
    Builds the tabulation.
    !!}
    use :: Display                       , only : displayIndent      , displayUnindent
    use :: Error              , only : errorStatusSuccess
    use :: Abundances_Structure          , only : abundances         , metallicityTypeLinearByMassSolar
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class  (taskBuildBroadbandLuminosityTabulations), intent(inout), target   :: self
    integer                                         , intent(  out), optional :: status
    type   (abundances                             )                          :: abundancesStellar
    type   (stellarLuminosities                    )                          :: luminosities
    
    call displayIndent ('Begin task: build broadband luminosity tabulations')
    ! Set luminosities, which will trigger tabulation if necessary.
    call abundancesStellar%metallicitySet (                                                                                      &
         &                                 metallicity                            =1.0d0                                       , &
         &                                 metallicityType                        =metallicityTypeLinearByMassSolar              &
         &                                )
    call luminosities     %setLuminosities(                                                                                      &
         &                                 mass                                   =1.0d0                                       , &
         &                                 time                                   =0.0d0                                       , &
         &                                 abundancesStellar                      =abundancesStellar                           , &
         &                                 stellarPopulation_                     =self%stellarPopulation_                     , &
         &                                 stellarPopulationBroadBandLuminosities_=self%stellarPopulationBroadBandLuminosities_  &
         &                                )
    ! Finish up.
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: build broadband luminosity tabulations')
    return
  end subroutine buildBroadbandLuminosityTabulationsPerform
  
  logical function buildBroadbandLuminosityTabulationsRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskBuildBroadbandLuminosityTabulations ), intent(inout) :: self
    !$GLC attributes unused :: self

    buildBroadbandLuminosityTabulationsRequiresOutputFile=.false.
    return
  end function buildBroadbandLuminosityTabulationsRequiresOutputFile
