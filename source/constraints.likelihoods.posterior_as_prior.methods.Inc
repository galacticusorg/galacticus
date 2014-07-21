!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines methods for a ``posterior as prior'' likelihood function.

function likelihoodPosteriorPriorConstructor(wrappedLikelihoodDefinition,chainBaseName,neighborCount,tolerance,configFileName)
  !% Constructor for ``posterior as prior'' regression likelihood class.
  use FoX_DOM
  use String_Handling
  use File_Utilities
  use ISO_Varying_String
  use Memory_Management
  implicit none
  type            (likelihoodPosteriorPrior)                              :: likelihoodPosteriorPriorConstructor
  type            (node                    ), intent(in   ), pointer      :: wrappedLikelihoodDefinition
  character       (len=*                   ), intent(in   )               :: chainBaseName
  integer                                   , intent(in   )               :: neighborCount
  double precision                          , intent(in   )               :: tolerance
  double precision                          , allocatable  , dimension(:) :: state
  type            (varying_string          ), intent(in   ), optional     :: configFileName
  character       (len=1024                )                              :: line
  integer                                                                 :: chainCount                , iChain               , &
       &                                                                     chainUnit                 , ioStatus             , &
       &                                                                     i                         , j 
  type            (varying_string          )                              :: chainFileName
  character       (len=4                   )                              :: label
  logical                                                                 :: converged
  double precision                                                        :: time                      , likelihood           , &
       &                                                                     logPriorNormalization

  ! Store wrapped likelihood, neighbor count, tolerance, and allocate workspace.
  likelihoodPosteriorPriorConstructor%wrappedLikelihoodFunction => likelihoodNew(wrappedLikelihoodDefinition,configFileName)
  likelihoodPosteriorPriorConstructor%tolerance                 =  tolerance
  likelihoodPosteriorPriorConstructor%neighborCount             =  neighborCount
  likelihoodPosteriorPriorConstructor%initialized               =  .false.
  call Alloc_Array(likelihoodPosteriorPriorConstructor%neighborIndices  ,[neighborCount])
  call Alloc_Array(likelihoodPosteriorPriorConstructor%neighborDistances,[neighborCount])
  ! Read MCMC chains.
  iChain=-1
  do while (iChain < 10000)
     iChain=iChain+1
     write (label,'(i4.4)') iChain
     chainFileName=chainBaseName//"_"//trim(adjustl(label))//".log"
     if (.not.File_Exists(chainFileName)) exit
     if (iChain == 0) then
        open(newUnit=chainUnit,file=char(chainFileName),status='old',form='formatted',iostat=ioStatus)
        read (chainUnit,'(a)') line
        likelihoodPosteriorPriorConstructor%dimCount=String_Count_Words(line)
        close(chainUnit)
    end if
  end do  
  chainCount=iChain+1
  i         = 0
  do iChain=0,chainCount-1
     write (label,'(i4.4)') iChain
     chainFileName=chainBaseName//"_"//trim(adjustl(label))//".log"
     if (iChain == 0) then
        ! Count number of converged states in the file.
        likelihoodPosteriorPriorConstructor%convergedStateCount=0
        open(newUnit=chainUnit,file=char(chainFileName),status='old',form='formatted',iostat=ioStatus)
        do while (ioStatus == 0)
           read (chainUnit,*,iostat=ioStatus) i,j,time,converged
           if (ioStatus /= 0) cycle
           if (converged) likelihoodPosteriorPriorConstructor%convergedStateCount=likelihoodPosteriorPriorConstructor%convergedStateCount+1
        end do
        close(chainUnit)
        ! Allocate storage for all converged points.
        likelihoodPosteriorPriorConstructor%convergedStateCount=likelihoodPosteriorPriorConstructor%convergedStateCount*chainCount
        call Alloc_Array(likelihoodPosteriorPriorConstructor%states,[likelihoodPosteriorPriorConstructor%convergedStateCount,likelihoodPosteriorPriorConstructor%dimCount])
        call Alloc_Array(                                    state ,[                                                        likelihoodPosteriorPriorConstructor%dimCount])
     end if
     open(newUnit=chainUnit,file=char(chainFileName),status='old',form='formatted',iostat=ioStatus)
     do while (ioStatus == 0)
        read (chainUnit,*,iostat=ioStatus) i,j,time,converged,likelihood,state
        if (ioStatus /= 0) cycle
        if (converged) then
           i=i+1
           likelihoodPosteriorPriorConstructor%states(i,:)=state
        end if
     end do
     close(chainUnit)
  end do
  return
end function likelihoodPosteriorPriorConstructor

subroutine likelihoodPosteriorPriorInitialize(self,parameterMappings)
  !% Initialize a ``posterior as prior'' likelihood object.
  use Memory_Management
  implicit none
  class           (likelihoodPosteriorPrior), intent(inout)               :: self
  type            (mappingList             ), intent(in   ), dimension(:) :: parameterMappings
  double precision                          , allocatable  , dimension(:) :: mean, meanSquared
  integer                                                                 :: i   , j

  ! Set initialization state.
  if (self%initialized) return
  self%initialized=.true.
  ! Map states.
  do i=1,self%convergedStateCount
     do j=1,self%dimCount
        self%states(i,j)=parameterMappings(j)%thisMapping%map(self%states(i,j))
     end do
  end do
  ! Estimate variance in each dimension.
  call Alloc_Array(mean             ,[self%dimCount])
  call Alloc_Array(meanSquared      ,[self%dimCount])
  call Alloc_Array(self%rootVariance,[self%dimCount])
  mean             =sum(self%states   ,dim=1)/dble(self%convergedStateCount)
  meanSquared      =sum(self%states**2,dim=1)/dble(self%convergedStateCount)
  self%rootVariance=sqrt(meanSquared-mean**2)  
  ! Scale states by root-variance.
  forall(i=1:self%convergedStateCount)
     self%states(i,:)=self%states(i,:)/self%rootVariance
  end forall
  ! Evaluate normalization factor for prior. Includes normalization for
  ! number of states used and conversion back to unscaled parameters.
  self%logPriorNormalization=-log(dble(self%convergedStateCount))-log(product(self%rootVariance))
  ! Construct the search object.  
  allocate(self%searcher)
  select type (s => self%searcher)
  type is (nearestNeighbors)
     s=nearestNeighbors(self%states)
  end select
  return
end subroutine likelihoodPosteriorPriorInitialize

subroutine likelihoodPosteriorPriorDestructor(self)
  !% Destructor for ``posterior as prior'' likelihood class.
  implicit none
  type(likelihoodPosteriorPrior), intent(inout) :: self
    
  deallocate(self%searcher)
  return
end subroutine likelihoodPosteriorPriorDestructor

double precision function likelihoodPosteriorPriorEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  !% Return the log-likelihood for a ``posterior as prior'' likelihood function.
  use Constraints_State
  use Constraints_Convergence
  use Numerical_Constants_Math
  implicit none
  class           (likelihoodPosteriorPrior), intent(inout)               :: self
  class           (state                   ), intent(in   )               :: simulationState
  type            (mappingList             ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence             ), intent(inout)               :: simulationConvergence
  double precision                          , intent(in   )               :: temperature          , logLikelihoodCurrent, &
       &                                                                     logPriorCurrent      , logPriorProposed
  double precision                          , allocatable  , dimension(:) :: stateVector
  double precision                                                        :: kernelScale          , weight
  integer                                                                 :: i
  
  ! Initialize.
  call self%initialize(parameterMappings)
  ! Get the simulation state.
  stateVector=simulationState%get()
  ! Scale the state vector.
  stateVector=stateVector/self%rootVariance
  ! Search for nearest neighbors.
  call self%searcher%search(stateVector,self%neighborCount,self%tolerance,self%neighborIndices,self%neighborDistances)
  ! Evaluate kernel density estimate of log-prior.
  kernelScale=self%neighborDistances(self%neighborCount)
  weight=0.0d0
  do i=1,self%neighborCount
     weight=weight+exp(-0.5d0*self%neighborDistances(i)/kernelScale)
  end do
  likelihoodPosteriorPriorEvaluate=                                      &
       & +self%logPriorNormalization                                     &
       & +log(weight)                                                    &
       & -0.5d0                                                          &
       & *dble(self%dimCount)                                            &
       & *log(2.0d0*Pi*kernelScale)                                      &
       & +self%wrappedLikelihoodFunction%evaluate(                       &
       &                                          simulationState      , &
       &                                          parameterMappings    , &
       &                                          simulationConvergence, &
       &                                          temperature          , &
       &                                          logLikelihoodCurrent , &
       &                                          logPriorCurrent      , &
       &                                          logPriorProposed       &
       &                                         )
  return
end function likelihoodPosteriorPriorEvaluate

subroutine likelihoodPosteriorPriorFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodPosteriorPrior), intent(inout) :: self
  
  return
end subroutine likelihoodPosteriorPriorFunctionChanged
