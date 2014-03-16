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

!% Defines methods for a latin hypercube state initializor.

function stateInitializorLatinHypercubeConstructor(maximinTrialCount)
  !% Constructor for the latin hypercube state initializor class.
  implicit none
  type   (stateInitializorLatinHypercube)                :: stateInitializorLatinHypercubeConstructor
  integer                                , intent(in   ) :: maximinTrialCount

  stateInitializorLatinHypercubeConstructor%maximinTrialCount=maximinTrialCount
  return
end function stateInitializorLatinHypercubeConstructor

subroutine stateInitializorLatinHypercubeInitialize(self,simulationState,parameterPriors)
  !% Initialize simulation state by drawing at random from the parameter priors.
  use MPI_Utilities
  use Pseudo_Random
  use Sort
  use, intrinsic :: ISO_C_Binding
  implicit none
  class           (stateInitializorLatinHypercube), intent(inout)                 :: self
  class           (state                         ), intent(inout)                 :: simulationState
  type            (prior                         ), intent(inout), dimension(:  ) :: parameterPriors
  integer         (kind=c_size_t                 ), allocatable  , dimension(:  ) :: order
  double precision                                , allocatable  , dimension(:  ) :: x
  double precision                                , allocatable  , dimension(:,:) :: stateGrid            , stateGridBest
  type            (pseudoRandom)                                                  :: randomNumberGenerator
  integer                                                                         :: i                    , j                       , &
       &                                                                             i1                   , i2                      , &
       &                                                                             k
  double precision                                                                :: separationMinimum    , separationMinimumMaximum

  ! Generate random sequence.
  allocate(        x    (0:mpiSelf%count()-1                            ))
  allocate(    order    (0:mpiSelf%count()-1                            ))
  allocate(stateGrid    (0:mpiSelf%count()-1,simulationState%dimension()))
  allocate(stateGridBest(0:mpiSelf%count()-1,simulationState%dimension()))
  separationMinimumMaximum=0.0d0
  do k=0,self%maximinTrialCount
     do j=1,simulationState%dimension()
        x                =0.0d0
        x(mpiSelf%rank())=randomNumberGenerator%sample(ompThreadOffset=.true.,mpiRankOffset=.true.)
        x                =mpiSelf%sum(x)
        order            =Sort_Index_Do(x)-1
        do i=0,mpiSelf%count()-1
           stateGrid(i,j)=(dble(order(i))+0.5d0)/mpiSelf%count()
        end do
     end do
     ! Find minimum separation.
     separationMinimum=1.0d30
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
     stateGridBest(mpiSelf%rank(),j)=parameterPriors(j)%invert(stateGridBest(mpiSelf%rank(),j))
  end do
  ! Store the state.
  call simulationState%update(stateGridBest(mpiSelf%rank(),:),.false.)
  ! Clean up workspace.
  deallocate(stateGridBest)
  deallocate(stateGrid    )
  deallocate(        x    )
  deallocate(    order    )
  return
end subroutine stateInitializorLatinHypercubeInitialize
