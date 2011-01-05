!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
