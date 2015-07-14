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

!% Defines methods for the ``adaptive'' differential evolution random jump class.

function deRandomJumpAdaptiveConstructor()
  !% Constructor for the ``adaptive'' differential evolution random jump class.
  implicit none
  type(deRandomJumpAdaptive) :: deRandomJumpAdaptiveConstructor

  return
end function deRandomJumpAdaptiveConstructor

function deRandomJumpAdaptiveSample(self,randomDistributions,simulationState)
  !% Sample from the random jump distribution.
  use Galacticus_Error
  use MPI_Utilities
  implicit none
  class           (deRandomJumpAdaptive)                                      , intent(inout) :: self
  type            (distributionList     ), dimension(:)                        , intent(in   ) :: randomDistributions
  class           (state                )                                      , intent(in   ) :: simulationState
  double precision                       , dimension(size(randomDistributions))                :: deRandomJumpAdaptiveSample
  double precision                       , dimension(:)                        , allocatable   :: parameterRange
  integer                                                                                      :: i

  ! Find the current range of each parameter.
  parameterRange=mpiSelf%maxval(simulationState%get())-mpiSelf%minval(simulationState%get())
  deRandomJumpAdaptiveSample=0.0d0
  do i=1,size(randomDistributions)
     select type (d => randomDistributions(i)%thisDistribution)
     class is (distribution1D)
        deRandomJumpAdaptiveSample(i)=deRandomJumpAdaptiveSample(i)+d%sample()*parameterRange(i)
     class default
        call Galacticus_Error_Report('deRandomJumpAdaptiveSample','only 1D distributions are supported')
     end select
  end do
  return
end function deRandomJumpAdaptiveSample
