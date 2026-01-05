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
Implements a filter that passes nodes only if they exist at an output time.
!!}

  use :: Output_Times, only : outputTimesClass

  !![
  <galacticFilter name="galacticFilterOutputTimes">
   <description>
     A filter that passes nodes only if they exist at an output time.
  </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterOutputTimes
     !!{
     A filter that passes nodes only if they exist at an output time.
     !!}
     private
     class           (outputTimesClass), pointer :: outputTimes_      => null()
     double precision                            :: toleranceRelative
   contains
     final     ::           outputTimesDestructor
     procedure :: passes => outputTimesPasses
  end type galacticFilterOutputTimes

  interface galacticFilterOutputTimes
     !!{
     Constructors for the \refClass{galacticFilterOutputTimes} galactic filter class.
     !!}
     module procedure outputTimesConstructorParameters
     module procedure outputTimesConstructorInternal
  end interface galacticFilterOutputTimes

contains

  function outputTimesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterOutputTimes} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterOutputTimes)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    class           (outputTimesClass         ), pointer       :: outputTimes_
    double precision                                           :: toleranceRelative

    !![
    <inputParameter>
      <name>toleranceRelative</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The fractional tolerance to allow when comparing the time at which a node exists to output times.</description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=galacticFilterOutputTimes(toleranceRelative,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function outputTimesConstructorParameters

  function outputTimesConstructorInternal(toleranceRelative,outputTimes_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterOutputTimes} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterOutputTimes)                        :: self
    double precision                           , intent(in   )         :: toleranceRelative
    class           (outputTimesClass         ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="toleranceRelative, *outputTimes_"/>
    !!]
    
    return
  end function outputTimesConstructorInternal

  subroutine outputTimesDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterOutputTimes} galactic filter class.
    !!}
    implicit none
    type(galacticFilterOutputTimes), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine outputTimesDestructor

  logical function outputTimesPasses(self,node) result(passes)
    !!{
    Filter based on whether a subhalo can impact a stream in the timestep.
    !!}
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: Galacticus_Nodes    , only : nodeComponentBasic
    use            :: Numerical_Comparison, only : Values_Agree
    implicit none
    class  (galacticFilterOutputTimes), intent(inout)          :: self
    type   (treeNode                 ), intent(inout), target  :: node
    class  (nodeComponentBasic       )               , pointer :: basic
    integer(c_size_t                 )                         :: i

    passes =  .false.
    basic  => node%basic()
    do i=1,self%outputTimes_%count()
       if (Values_Agree(basic%time(),self%outputTimes_%time(i),relTol=self%toleranceRelative)) then
          passes=.true.
          exit
       end if
    end do
    return
  end function outputTimesPasses
