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
Implements an N-body data operator which applies some other operator to a selected simulation.
!!}
  
  !![
  <nbodyOperator name="nbodyOperatorSimulationSelector">
   <description>An N-body data operator which applies some other operator to a selected simulation.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSimulationSelector
     !!{
     An N-body data operator which applies some other operator to a selected simulation.
     !!}
     private
     class  (nbodyOperatorClass), pointer :: nbodyOperator_  => null()
     integer                              :: indexSimulation
   contains
     final     ::            simulationSelectorDestructor
     procedure :: operate => simulationSelectorOperate
  end type nbodyOperatorSimulationSelector

  interface nbodyOperatorSimulationSelector
     !!{
     Constructors for the {\normalfont \ttfamily simulationSelector} N-body operator class.
     !!}
     module procedure simulationSelectorConstructorParameters
     module procedure simulationSelectorConstructorInternal
  end interface nbodyOperatorSimulationSelector

contains

  function simulationSelectorConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily simulationSelector} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nbodyOperatorSimulationSelector)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    class  (nbodyOperatorClass             ), pointer       :: nbodyOperator_
    integer                                                 :: indexSimulation

    !![
    <inputParameter>
      <name>indexSimulation</name>
      <source>parameters</source>
      <description>The index of the simulation to which to apply the operator.</description>
    </inputParameter>
    <objectBuilder class="nbodyOperator" name="nbodyOperator_" source="parameters"/>
    !!]
    self=nbodyOperatorSimulationSelector(indexSimulation,nbodyOperator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nbodyOperator_"/>
    !!]
    return
  end function simulationSelectorConstructorParameters

  function simulationSelectorConstructorInternal(indexSimulation,nbodyOperator_) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily simulationSelector} N-body operator class.
    !!}
    implicit none
    type   (nbodyOperatorSimulationSelector)                        :: self
    integer                                 , intent(in   )         :: indexSimulation
    class  (nbodyOperatorClass             ), intent(in   ), target :: nbodyOperator_
    !![
    <constructorAssign variables="indexSimulation, *nbodyOperator_"/>
    !!]

    return
  end function simulationSelectorConstructorInternal

  subroutine simulationSelectorDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily simulationSelector} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorSimulationSelector), intent(inout) :: self

    !![
    <objectDestructor name="self%nbodyOperator_"/>
    !!]
    return
  end subroutine simulationSelectorDestructor

  subroutine simulationSelectorOperate(self,simulations)
    !!{
    Operate on the selected simulation.
    !!}
    use :: Display, only : displayIndent, displayUnindent, verbosityLevelStandard
    implicit none
    class  (nbodyOperatorSimulationSelector), intent(inout)                 :: self
    type   (nBodyData                      ), intent(inout), dimension(  :) :: simulations

    call displayIndent('operate on selected simulation',verbosityLevelStandard)
    call self%nbodyOperator_%operate(simulations(self%indexSimulation:self%indexSimulation))
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine simulationSelectorOperate
