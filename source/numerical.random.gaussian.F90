!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements Gaussian random deviates.

module Gaussian_Random
  !% Implements Gaussian random deviates.
  use FGSL
  use, intrinsic :: ISO_C_Binding
  private
  public :: Gaussian_Random_Get, Gaussian_Random_Free
  
  logical           :: Seed_Is_Set=.false.
  integer(c_size_t) :: gaussianRandomSeedC

contains

  double precision function Gaussian_Random_Get(pseudoSequenceObject,width,reset)
    !% Returns a Gaussian random deviate.
    use Input_Parameters
    implicit none
    type(fgsl_rng),   intent(inout)           :: pseudoSequenceObject
    double precision, intent(in)              :: width
    logical,          intent(inout), optional :: reset
    integer                                   :: gaussianRandomSeed
    logical                                   :: resetActual

    ! Determine if we need to reset.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    else
       resetActual=.false.
    end if
    
    ! Read in the random number seed if necessary.
    !$omp critical (Gaussian_Random_Get)
    if (.not.Seed_Is_Set) then
       !@ <inputParameter>
       !@   <name>gaussianRandomSeed</name>
       !@   <defaultValue>843</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     A seed for the Gaussian random number generator.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('gaussianRandomSeed',gaussianRandomSeed,defaultValue=843)
       gaussianRandomSeedC=gaussianRandomSeed
       Seed_Is_Set=.true.
    end if
    !$omp end critical (Gaussian_Random_Get)
    
    if (resetActual.or..not.FGSL_Well_Defined(pseudoSequenceObject)) then
       pseudoSequenceObject=FGSL_RNG_Alloc(FGSL_RNG_Default)
       gaussianRandomSeed=gaussianRandomSeed+1
       call FGSL_RNG_Set(pseudoSequenceObject,gaussianRandomSeedC)
    end if
    Gaussian_Random_Get=FGSL_Ran_Gaussian(pseudoSequenceObject,width)
    return
  end function Gaussian_Random_Get
  
  subroutine Gaussian_Random_Free(pseudoSequenceObject)
    !% Frees a pseudo-random sequence object.
    implicit none
    type(fgsl_rng), intent(inout) :: pseudoSequenceObject
    
    call FGSL_RNG_Free(pseudoSequenceObject)
    return
  end subroutine Gaussian_Random_Free
  
end module Gaussian_Random
