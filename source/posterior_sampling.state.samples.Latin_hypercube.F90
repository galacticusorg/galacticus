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
  Implementation of a posterior state samples class which draws samples from a Latin hypercube in the cumulative distribution of the priors.
  !!}
  
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <posteriorSamples name="posteriorSamplesLatinHypercube">
   <description>
    A posterior state samples class which draws samples from a Latin hypercube in the cumulative distribution of the priors.
   </description>
  </posteriorSamples>
  !!]
  type, extends(posteriorSamplesClass) :: posteriorSamplesLatinHypercube
     !!{
     Implementation of a posterior state samples class which draws samples from a Latin hypercube in the cumulative distribution of the priors.
     !!}
     private
     integer                                      :: countSamples                    , maximinTrialCount
     class  (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::             latinHypercubeDestructor
     procedure :: samples  => latinHypercubeSamples
  end type posteriorSamplesLatinHypercube

  interface posteriorSamplesLatinHypercube
     !!{
     Constructors for the \refClass{posteriorSamplesLatinHypercube} posterior sampling state initialization class.
     !!}
     module procedure latinHypercubeConstructorParameters
     module procedure latinHypercubeConstructorInternal
  end interface posteriorSamplesLatinHypercube

contains

  function latinHypercubeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSamplesLatinHypercube} posterior sampling state initialization class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (posteriorSamplesLatinHypercube)                :: self
    type   (inputParameters               ), intent(inout) :: parameters
    class  (randomNumberGeneratorClass    ), pointer       :: randomNumberGenerator_
    integer                                                :: countSamples          , maximinTrialCount
    
    !![
    <inputParameter>
      <name>countSamples</name>
      <description>The number of samples to draw.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>maximinTrialCount</name>
      <defaultValue>1000</defaultValue>
      <description>The number of trial Latin Hypercubes to construct when seeking the maximum minimum separation sample.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=posteriorSamplesLatinHypercube(countSamples,maximinTrialCount,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function latinHypercubeConstructorParameters

  function latinHypercubeConstructorInternal(countSamples,maximinTrialCount,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{posteriorSamplesLatinHypercube} posterior sampling state initialization class.
    !!}
    implicit none
    type   (posteriorSamplesLatinHypercube)                        :: self
    integer                                , intent(in   )         :: countSamples          , maximinTrialCount
    class  (randomNumberGeneratorClass    ), intent(in   ), target :: randomNumberGenerator_

    !![
    <constructorAssign variables="countSamples, maximinTrialCount, *randomNumberGenerator_"/>
    !!]
    return
  end function latinHypercubeConstructorInternal

  subroutine latinHypercubeDestructor(self)
    !!{
    Destructor for the \refClass{posteriorSamplesLatinHypercube} posterior sampling state initialization class.
    !!}
    implicit none
    type(posteriorSamplesLatinHypercube), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine latinHypercubeDestructor

  subroutine latinHypercubeSamples(self,simulationStates,modelParameters_)
    !!{
    Generate simulation states by sampling a Latin hypercube in the cumulative distribution of priors.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: MPI_Utilities, only : mpiBarrier, mpiSelf
    use            :: Sorting      , only : sortIndex
    implicit none
    class           (posteriorSamplesLatinHypercube), intent(inout)                              :: self
    type            (posteriorSampleStateSimple    ), intent(inout), dimension(:  ), allocatable :: simulationStates
    type            (modelParameterList            ), intent(inout), dimension(:  )              :: modelParameters_
    integer         (kind=c_size_t                 )               , dimension(:  ), allocatable :: order
    double precision                                               , dimension(:  ), allocatable :: x                , y
    double precision                                               , dimension(:,:), allocatable :: stateGrid        , stateGridBest
    integer                                                                                      :: i                , j                       , &
         &                                                                                          i1               , i2                      , &
         &                                                                                          k
    double precision                                                                             :: separationMinimum, separationMinimumMaximum

    ! Generate random sequence.
    allocate(        x       (0:self%countSamples-1                       ))
    allocate(        y       (0:self%countSamples-1                       ))
    allocate(    order       (0:self%countSamples-1                       ))
    allocate(stateGrid       (0:self%countSamples-1,size(modelParameters_)))
    allocate(stateGridBest   (0:self%countSamples-1,size(modelParameters_)))
    allocate(simulationStates(1:self%countSamples                         ))
    separationMinimumMaximum=0.0d0
    do k=1,self%maximinTrialCount
       do j=1,size(modelParameters_)
          x                =0.0d0
          do i=0,self%countSamples-1
             if (mod(i,mpiSelf%count()) /= mpiSelf%rank()) cycle
             x(i)=self%randomNumberGenerator_%uniformSample()
          end do
          y    =mpiSelf%sum(x)
          call mpiBarrier()
          order=sortIndex(y)-1
          do i=0,self%countSamples-1
             stateGrid(i,j)=(dble(order(i))+0.5d0)/self%countSamples
          end do
       end do
       ! Find minimum separation.
       separationMinimum=huge(1.0d0)
       do i1=0,self%countSamples-1
          do i2=i1+1,self%countSamples-1
             separationMinimum=min(separationMinimum,sum((stateGrid(i1,:)-stateGrid(i2,:))**2))
          end do
       end do
       ! Check if this is the maximum minimum separation yet found.
       if (separationMinimum > separationMinimumMaximum) then
          separationMinimumMaximum=separationMinimum
          stateGridBest           =stateGrid
       end if
    end do
    ! Convert the cumulative density to a value of the prior.
    do i=0,self%countSamples-1
       do j=1,size(modelParameters_)
          stateGridBest(i,j)=modelParameters_(j)%modelParameter_%map        (      &
               &             modelParameters_(j)%modelParameter_%priorInvert (     &
               &             stateGridBest                                    (    &
               &                                                               i,j &
               &                                                              )    &
               &                                                             )     &
               &                                                            )
       end do
       simulationStates(i+1)=posteriorSampleStateSimple(acceptedStateCount=1)
       call simulationStates(i+1)%parameterCountSet(size(modelParameters_))       
       call simulationStates(i+1)%update(stateGridBest(i,:),.false.,.false.)
    end do
    return
  end subroutine latinHypercubeSamples
