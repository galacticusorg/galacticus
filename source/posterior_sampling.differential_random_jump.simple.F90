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
  Implementation of a posterior sampling differential evolution random jump class in which the jump is drawn from a fixed
  distribution.
  !!}

  !![
  <posteriorSampleDffrntlEvltnRandomJump name="posteriorSampleDffrntlEvltnRandomJumpSimple">
   <description>
    In this class, the random jumps are drawn directly from the distributions specified in the {\normalfont \ttfamily random} object
    of each \refClass{modelParameterClass} object.
   </description>
  </posteriorSampleDffrntlEvltnRandomJump>
  !!]
  type, extends(posteriorSampleDffrntlEvltnRandomJumpClass) :: posteriorSampleDffrntlEvltnRandomJumpSimple
     !!{
     Implementation of a posterior sampling differential evolution random jump class in which the jump is drawn from a fixed
     distribution.
     !!}
     private
   contains
     procedure :: sample => simpleSample
  end type posteriorSampleDffrntlEvltnRandomJumpSimple

  interface posteriorSampleDffrntlEvltnRandomJumpSimple
     !!{
     Constructors for the \refClass{posteriorSampleDffrntlEvltnRandomJumpSimple} posterior sampling differential evolution random jump class.
     !!}
     module procedure simpleConstructorParameters
  end interface posteriorSampleDffrntlEvltnRandomJumpSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleDffrntlEvltnRandomJumpSimple} posterior sampling differential evolution random jump class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(posteriorSampleDffrntlEvltnRandomJumpSimple)                 :: self
    type(inputParameters                            ), intent(inout)  :: parameters

    self=posteriorSampleDffrntlEvltnRandomJumpSimple()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function simpleConstructorParameters

  function simpleSample(self,modelParameters_,simulationState)
    !!{
    Sample from the random jump distribution.
    !!}
    implicit none
    class           (posteriorSampleDffrntlEvltnRandomJumpSimple)                                   , intent(inout) :: self
    type            (modelParameterList                         ), dimension(:)                     , intent(in   ) :: modelParameters_
    class           (posteriorSampleStateClass                  )                                   , intent(inout) :: simulationState
    double precision                                             , dimension(size(modelParameters_))                :: simpleSample
    integer                                                                                                         :: i
    !$GLC attributes unused :: self, simulationState

    simpleSample=0.0d0
    do i=1,size(modelParameters_)
       simpleSample(i)=simpleSample(i)+modelParameters_(i)%modelParameter_%randomPerturbation()
    end do
    return
  end function simpleSample
