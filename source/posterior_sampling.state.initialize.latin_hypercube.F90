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
  Implementation of Latin hypercube state initializer.
  !!}

  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <posteriorSampleStateInitialize name="posteriorSampleStateInitializeLatinHypercube">
   <description>
    This class uses a \gls{latinhypercube} design (in the cumulative prior probability distribution of each parameter) to assign
    initial state vectors. In particular, a \gls{maximin} design is used in which a number of trial \glspl{latinhypercube} are
    constructed and the hypercube with the greatest minimum distance between any pair of state vectors is selected. The {\normalfont
    \ttfamily [maximinTrialCount]} parameter is used to specify the number of trial hypercubes to construct.
   </description>
  </posteriorSampleStateInitialize>
  !!]
  type, extends(posteriorSampleStateInitializeClass) :: posteriorSampleStateInitializeLatinHypercube
     !!{
     Implementation of a posterior sampling state initialization class which samples the initial state at random from the priors using Latin Hypercube sampling.
     !!}
     private
     integer                                      :: maximinTrialCount
     class  (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::                latinHypercubeDestructor
     procedure :: initialize  => latinHypercubeInitialize
  end type posteriorSampleStateInitializeLatinHypercube

  interface posteriorSampleStateInitializeLatinHypercube
     !!{
     Constructors for the \refClass{posteriorSampleStateInitializeLatinHypercube} posterior sampling state initialization class.
     !!}
     module procedure latinHypercubeConstructorParameters
     module procedure latinHypercubeConstructorInternal
  end interface posteriorSampleStateInitializeLatinHypercube

contains

  function latinHypercubeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializeLatinHypercube} posterior sampling state initialization class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (posteriorSampleStateInitializeLatinHypercube)                :: self
    type   (inputParameters                             ), intent(inout) :: parameters
    class  (randomNumberGeneratorClass                  ), pointer       :: randomNumberGenerator_
    integer                                                              :: maximinTrialCount

    !![
    <inputParameter>
      <name>maximinTrialCount</name>
      <defaultValue>1000</defaultValue>
      <description>The number of trial Latin Hypercubes to construct when seeking the maximum minimum separation sample.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=posteriorSampleStateInitializeLatinHypercube(maximinTrialCount,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function latinHypercubeConstructorParameters

  function latinHypercubeConstructorInternal(maximinTrialCount,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializeLatinHypercube} posterior sampling state initialization class.
    !!}
    implicit none
    type   (posteriorSampleStateInitializeLatinHypercube)                        :: self
    integer                                              , intent(in   )         :: maximinTrialCount
    class  (randomNumberGeneratorClass                  ), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="maximinTrialCount, *randomNumberGenerator_"/>
    !!]

    return
  end function latinHypercubeConstructorInternal

  subroutine latinHypercubeDestructor(self)
    !!{
    Destructor for the \refClass{posteriorSampleStateInitializeLatinHypercube} posterior sampling state initialization class.
    !!}
    implicit none
    type(posteriorSampleStateInitializeLatinHypercube), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine latinHypercubeDestructor

  subroutine latinHypercubeInitialize(self,simulationState,modelParameters_,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
    !!{
    Initialize simulation state by drawing at random from the parameter priors.
    !!}
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: MPI_Utilities               , only : mpiBarrier   , mpiSelf
    use            :: Models_Likelihoods_Constants, only : logImpossible
    use            :: Sorting                     , only : sortIndex
    implicit none
    class           (posteriorSampleStateInitializeLatinHypercube), intent(inout)                 :: self
    class           (posteriorSampleStateClass                   ), intent(inout)                 :: simulationState
    class           (posteriorSampleLikelihoodClass              ), intent(inout)                 :: modelLikelihood
    type            (modelParameterList                          ), intent(inout), dimension(:  ) :: modelParameters_
    double precision                                              , intent(  out)                 :: timeEvaluatePrevious, logLikelihood , &
         &                                                                                           logPosterior
    integer         (kind=c_size_t                               ), allocatable  , dimension(:  ) :: order
    double precision                                              , allocatable  , dimension(:  ) :: x                    , y
    double precision                                              , allocatable  , dimension(:,:) :: stateGrid            , stateGridBest
    integer                                                                                       :: i                    , j                       , &
         &                                                                                           i1                   , i2                      , &
         &                                                                                           k
    double precision                                                                              :: separationMinimum    , separationMinimumMaximum
    !$GLC attributes unused :: modelLikelihood

    ! No knowledge of evaluation time.
    timeEvaluatePrevious=-1.0d0
    ! We have no information about the likelihood of this state.
    logLikelihood=logImpossible
    logPosterior =logImpossible
    ! Generate random sequence.
    allocate(        x    (0:mpiSelf%count()-1                            ))
    allocate(        y    (0:mpiSelf%count()-1                            ))
    allocate(    order    (0:mpiSelf%count()-1                            ))
    allocate(stateGrid    (0:mpiSelf%count()-1,simulationState%dimension()))
    allocate(stateGridBest(0:mpiSelf%count()-1,simulationState%dimension()))
    separationMinimumMaximum=0.0d0
    do k=0,self%maximinTrialCount
       do j=1,simulationState%dimension()
          x                =0.0d0
          x(mpiSelf%rank())=self%randomNumberGenerator_%uniformSample()
          y                =mpiSelf%sum(x)
          call mpiBarrier()
          order            =sortIndex(y)-1
          do i=0,mpiSelf%count()-1
             stateGrid(i,j)=(dble(order(i))+0.5d0)/mpiSelf%count()
          end do
       end do
       ! Find minimum separation.
       separationMinimum=huge(1.0d0)
       do i1=0,mpiSelf%count()-1
          do i2=i1+1,mpiSelf%count()-1
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
    do j=1,simulationState%dimension()
       stateGridBest(mpiSelf%rank(),j)=modelParameters_(j)%modelParameter_%map        (                  &
            &                          modelParameters_(j)%modelParameter_%priorInvert (                 &
            &                          stateGridBest                                    (                &
            &                                                                            mpiSelf%rank(), &
            &                                                                            j               &
            &                                                                           )                &
            &                                                                          )                 &
            &                                                                         )
    end do
    ! Store the state.
    call simulationState%update(stateGridBest(mpiSelf%rank(),:),.false.,.false.)
    ! Clean up workspace.
    deallocate(stateGridBest)
    deallocate(stateGrid    )
    deallocate(        x    )
    deallocate(    order    )
    return
  end subroutine latinHypercubeInitialize
