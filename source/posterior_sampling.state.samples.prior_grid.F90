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
  Implementation of a posterior state samples class which draws samples from a grid in the cumulative distribution of the priors.
  !!}

  !![
  <posteriorSamples name="posteriorSamplesPriorGrid">
   <description>
    A posterior state samples class which draws samples from a grid in the cumulative distribution of the priors.
   </description>
  </posteriorSamples>
  !!]
  type, extends(posteriorSamplesClass) :: posteriorSamplesPriorGrid
     !!{
     Implementation of a posterior state samples class which draws samples from a grid in the cumulative distribution of the priors.
     !!}
     private
     integer :: countGrid
   contains
     procedure :: samples  => priorGridSamples
  end type posteriorSamplesPriorGrid

  interface posteriorSamplesPriorGrid
     !!{
     Constructors for the \refClass{posteriorSamplesPriorGrid} posterior sampling state initialization class.
     !!}
     module procedure priorGridConstructorParameters
     module procedure priorGridConstructorInternal
  end interface posteriorSamplesPriorGrid

contains

  function priorGridConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSamplesPriorGrid} posterior sampling state initialization class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (posteriorSamplesPriorGrid)                :: self
    type   (inputParameters          ), intent(inout) :: parameters
    integer                                           :: countGrid
    
    !![
    <inputParameter>
      <name>countGrid</name>
      <description>The number of grid steps in each parameter.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSamplesPriorGrid(countGrid)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function priorGridConstructorParameters

  function priorGridConstructorInternal(countGrid) result(self)
    !!{
    Constructor for the \refClass{posteriorSamplesPriorGrid} posterior sampling state initialization class.
    !!}
    implicit none
    type   (posteriorSamplesPriorGrid)                :: self
    integer                           , intent(in   ) :: countGrid
    !![
    <constructorAssign variables="countGrid"/>
    !!]
    
    return
  end function priorGridConstructorInternal

  subroutine priorGridSamples(self,simulationStates,modelParameters_)
    !!{
    Generate simulation states by sampling on a grid in the cumulative distribution of priors.
    !!}
    use, intrinsic :: ISO_C_Binding , only : c_size_t
    use            :: Multi_Counters, only : multiCounter
    implicit none
    class           (posteriorSamplesPriorGrid ), intent(inout)                                                 :: self
    type            (posteriorSampleStateSimple), intent(inout), dimension(                    : ), allocatable :: simulationStates
    type            (modelParameterList        ), intent(inout), dimension(                    : )              :: modelParameters_
    double precision                                           , dimension(size(modelParameters_))              :: stateVector
    type            (multiCounter              )                                                                :: counter
    integer         (c_size_t                  )                                                                :: i               , j

    allocate(simulationStates(self%countGrid**size(modelParameters_)))
    counter=multiCounter(spread(self%countGrid,1,size(modelParameters_)))
    i      =0
    do while (counter%increment())
       i=i+1
       ! Construct the state.
       do j=1,size(modelParameters_)
          stateVector(j)=modelParameters_(j)%modelParameter_%map        (                                     &
               &         modelParameters_(j)%modelParameter_%priorInvert (                                    &
               &                                                          +(dble(counter%state    (j))-0.5d0) &
               &                                                          / dble(self   %countGrid   )        &
               &                                                         )                                    &
               &                                                        )
       end do
       ! Store the state vector.
       simulationStates(i)=posteriorSampleStateSimple(acceptedStateCount=1)
       call simulationStates(i)%parameterCountSet(size(modelParameters_))
       call simulationStates(i)%update(stateVector,.false.,.false.)
    end do
    return
  end subroutine priorGridSamples
