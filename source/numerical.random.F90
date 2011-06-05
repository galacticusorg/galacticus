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


!% Contains a module which implements pseudo-random numbers.

module Pseudo_Random
  !% Implements pseudo-random numbers.
  use FGSL
  use, intrinsic :: ISO_C_Binding
  private
  public :: Pseudo_Random_Get, Pseudo_Random_Free, Pseudo_Random_Store, Pseudo_Random_Retrieve
  
  logical           :: Seed_Is_Set=.false.
  integer(c_size_t) :: randomSeedC

contains

  double precision function Pseudo_Random_Get(pseudoSequenceObject,reset)
    !% Returns a scalar giving a pseudo-random number.
    use Input_Parameters
    implicit none
    type(fgsl_rng), intent(inout)           :: pseudoSequenceObject
    logical,        intent(inout), optional :: reset
    integer                                 :: randomSeed
    logical                                 :: resetActual

    ! Determine if we need to reset.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    else
       resetActual=.false.
    end if

    ! Read in the random number seed if necessary.
    !$omp critical (Pseudo_Random_Get)
    if (.not.Seed_Is_Set) then
       !@ <inputParameter>
       !@   <name>randomSeed</name>
       !@   <defaultValue>219</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     A seed value for the random number generator.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('randomSeed',randomSeed,defaultValue=219)
       randomSeedC=randomSeed
       Seed_Is_Set=.true.
    end if
    !$omp end critical (Pseudo_Random_Get)
    
    if (resetActual.or..not.FGSL_Well_Defined(pseudoSequenceObject)) then
       pseudoSequenceObject=FGSL_RNG_Alloc(FGSL_RNG_Default)
       randomSeed=randomSeed+1
       call FGSL_RNG_Set(pseudoSequenceObject,randomSeedC)
    end if
    Pseudo_Random_Get=FGSL_RNG_Uniform(pseudoSequenceObject)
    return
  end function Pseudo_Random_Get
  
  subroutine Pseudo_Random_Free(pseudoSequenceObject)
    !% Frees a pseudo-random sequence object.
    implicit none
    type(fgsl_rng), intent(inout) :: pseudoSequenceObject
    
    call FGSL_RNG_Free(pseudoSequenceObject)
    return
  end subroutine Pseudo_Random_Free

  subroutine Pseudo_Random_Store(pseudoSequenceObject,fgslFile)
    !% Stores a pseudo-random sequence object to file.
    implicit none
    type(fgsl_rng),  intent(in) :: pseudoSequenceObject
    type(fgsl_file), intent(in) :: fgslFile
    integer                     :: iError

    iError=FGSL_Rng_FWrite(fgslFile,pseudoSequenceObject)
    return
  end subroutine Pseudo_Random_Store

  subroutine Pseudo_Random_Retrieve(pseudoSequenceObject,fgslFile)
    !% Stores a pseudo-random sequence object to file.
    implicit none
    type(fgsl_rng),  intent(inout) :: pseudoSequenceObject
    type(fgsl_file), intent(in)    :: fgslFile
    integer                        :: iError

    if (.not.FGSL_Well_Defined(pseudoSequenceObject)) pseudoSequenceObject=FGSL_RNG_Alloc(FGSL_RNG_Default)
    iError=FGSL_Rng_FRead(fgslFile,pseudoSequenceObject)
    return
  end subroutine Pseudo_Random_Retrieve

end module Pseudo_Random
