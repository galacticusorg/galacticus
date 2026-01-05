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
  Implementation of a posterior sampling differential evolution random jump class in which the jump is drawn from an adaptive
  distribution which scales with the range spanned by the sample states.
  !!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t
  
  !![
  <posteriorSampleDffrntlEvltnRandomJump name="posteriorSampleDffrntlEvltnRandomJumpAdaptive">
   <description>
    The random jumps are drawn from the distributions specified in the {\normalfont \ttfamily random} element of each
    \refClass{modelParameterClass} object and then multiplied by the currently occupied range of each parameter (i.e. the maximum
    value of the parameter over all current chain states minus the minimum value of each parameter over all current chain
    states).
   </description>
  </posteriorSampleDffrntlEvltnRandomJump>
  !!]
  type, extends(posteriorSampleDffrntlEvltnRandomJumpClass) :: posteriorSampleDffrntlEvltnRandomJumpAdaptive
     !!{
     Implementation of a posterior sampling differential evolution random jump class in which the jump is drawn from an
     adaptive distribution which scales with the range spanned by the sample states.
     !!}
     private
     integer(c_size_t) :: countLargeJumps , countJumps
     logical           :: warnedLargeJumps
   contains
     procedure :: sample => adaptiveSample
  end type posteriorSampleDffrntlEvltnRandomJumpAdaptive

  interface posteriorSampleDffrntlEvltnRandomJumpAdaptive
     !!{
     Constructors for the \refClass{posteriorSampleDffrntlEvltnRandomJumpAdaptive} posterior sampling differential evolution random jump class.
     !!}
     module procedure adaptiveConstructorParameters
     module procedure adaptiveConstructorInternal
  end interface posteriorSampleDffrntlEvltnRandomJumpAdaptive

contains

  function adaptiveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleDffrntlEvltnRandomJumpAdaptive} posterior sampling differential evolution random jump class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(posteriorSampleDffrntlEvltnRandomJumpAdaptive)                 :: self
    type(inputParameters                              ), intent(inout)  :: parameters

    self=posteriorSampleDffrntlEvltnRandomJumpAdaptive()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function adaptiveConstructorParameters

  function adaptiveConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily adaptive} posterior sampling differential evolution random jump class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(posteriorSampleDffrntlEvltnRandomJumpAdaptive) :: self

    self%warnedLargeJumps=.false.
    self%countJumps      =0_c_size_t
    self%countLargeJumps =0_c_size_t
    return
  end function adaptiveConstructorInternal

  function adaptiveSample(self,modelParameters_,simulationState)
    !!{
    Sample from the random jump distribution.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    use :: Display      , only : displayMessage, displayMagenta, displayReset
    implicit none
    class           (posteriorSampleDffrntlEvltnRandomJumpAdaptive)                                   , intent(inout) :: self
    type            (modelParameterList                           ), dimension(:)                     , intent(in   ) :: modelParameters_
    class           (posteriorSampleStateClass                    )                                   , intent(inout) :: simulationState
    double precision                                               , dimension(size(modelParameters_))                :: adaptiveSample
    double precision                                               , dimension(size(modelParameters_))                :: parameterRange
    integer         (c_size_t                                     ), parameter                                        :: countMinimum    =100_c_size_t
    integer                                                                                                           :: i
    double precision                                                                                                  :: sizeJump
    logical                                                                                                           :: jumpIsLarge

    ! Find the current range of each parameter.
    parameterRange=+mpiSelf%maxval(simulationState%get()) &
         &         -mpiSelf%minval(simulationState%get())
    jumpIsLarge   =.false.
    adaptiveSample=0.0d0
    do i=1,size(modelParameters_)
       sizeJump         =modelParameters_(i)%modelParameter_%randomPerturbation()
       adaptiveSample(i)=adaptiveSample(i)+sizeJump*parameterRange(i)
       if (abs(sizeJump) >= 1.0d0) jumpIsLarge=.true.
    end do
    if (.not.self%warnedLargeJumps) then
       if (jumpIsLarge)                                            &
            & self%countLargeJumps=self%countLargeJumps+1_c_size_t
       self       %countJumps     =self%countJumps     +1_c_size_t
       if (self%countJumps > countMinimum .and. self%countLargeJumps > self%countJumps/10_c_size_t) then
          call displayMessage(displayMagenta()//'WARNING:'//displayReset()//': random jumps to model parameters are O(1) more than 10% of the time - consider using smaller jumps')
          self%warnedLargeJumps=.true.
       end if
    end if
    return
  end function adaptiveSample

