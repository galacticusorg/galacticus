!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements Poisson random deviates.

module Poisson_Random
  !% Implements Poisson random deviates.
  use FGSL
  implicit none
  private
  public :: Poisson_Random_Get, Poisson_Random_Free

  logical                 :: Seed_Is_Set        =.false.
  integer                 :: poissonRandomSeed
  integer(kind=fgsl_long) :: poissonRandomSeedC=-1
  !$omp threadprivate(poissonRandomSeedC)
contains

  integer function Poisson_Random_Get(pseudoSequenceObject,mean,reset)
    !% Returns a Poisson random deviate.
    use Input_Parameters
    !$ use OMP_Lib
    implicit none
    type            (fgsl_rng), intent(inout)           :: pseudoSequenceObject
    double precision          , intent(in   )           :: mean
    logical                   , intent(inout), optional :: reset
    logical                                             :: resetActual

    ! Determine if we need to reset.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    else
       resetActual=.false.
    end if

    ! Read in the random number seed if necessary.
    !$omp critical (Poisson_Random_Get)
    if (.not.Seed_Is_Set) then
       !@ <inputParameter>
       !@   <name>poissonRandomSeed</name>
       !@   <defaultValue>348</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     A seed for the Poisson random number generator.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('poissonRandomSeed',poissonRandomSeed,defaultValue=843)
       Seed_Is_Set=.true.
    end if
    !$omp end critical (Poisson_Random_Get)

    if (resetActual.or..not.FGSL_Well_Defined(pseudoSequenceObject)) then
       pseudoSequenceObject=FGSL_RNG_Alloc(FGSL_RNG_Default)
       if (poissonRandomSeedC < 0) then
          poissonRandomSeedC=poissonRandomSeed
          !$ poissonRandomSeedC=poissonRandomSeedC+omp_get_thread_num()
       end if
       poissonRandomSeedC=poissonRandomSeedC+1
       call FGSL_RNG_Set(pseudoSequenceObject,poissonRandomSeedC)
    end if
    Poisson_Random_Get=FGSL_Ran_Poisson(pseudoSequenceObject,mean)
    return
  end function Poisson_Random_Get

  subroutine Poisson_Random_Free(pseudoSequenceObject)
    !% Frees a pseudo-random sequence object.
    implicit none
    type(fgsl_rng), intent(inout) :: pseudoSequenceObject

    call FGSL_RNG_Free(pseudoSequenceObject)
    return
  end subroutine Poisson_Random_Free

end module Poisson_Random
