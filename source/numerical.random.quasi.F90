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


!% Contains a module which implements quasi-random sequences.

module Quasi_Random
  !% Implements quasi-random sequences.
  use FGSL
  private
  public :: Quasi_Random_Get, Quasi_Random_Free

  interface Quasi_Random_Get
     module procedure Quasi_Random_Get_Scalar
     module procedure Quasi_Random_Get_Array
  end interface

contains

  double precision function Quasi_Random_Get_Scalar(quasiSequenceObject,reset,quasiSequenceType)
    !% Returns a scalar giving a quasi-random point in a 1-dimensional space.
    implicit none
    type(fgsl_qrng),      intent(inout)           :: quasiSequenceObject
    logical,              intent(inout), optional :: reset
    type(fgsl_qrng_type), intent(in),    optional :: quasiSequenceType
    integer(kind=4),      parameter               :: nGet=1
    double precision                              :: pointArray(1)
    logical                                       :: resetActual
    type(fgsl_qrng_type)                          :: quasiSequenceTypeActual

    ! Determine what type of quasi-random sequence to use.
    if (present(quasiSequenceType)) then
       quasiSequenceTypeActual=quasiSequenceType
    else
       quasiSequenceTypeActual=FGSL_qRng_Sobol
    end if
    
    ! Determine if we need to reset.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    else
       resetActual=.false.
    end if
    if (resetActual.or..not.FGSL_Well_Defined(quasiSequenceObject)) quasiSequenceObject=FGSL_qRng_Alloc(quasiSequenceTypeActual&
         &,1)
    pointArray=Quasi_Random_Get_Array(quasiSequenceObject,nGet,resetActual,quasiSequenceTypeActual)
    Quasi_Random_Get_Scalar=pointArray(1)
    return
  end function Quasi_Random_Get_Scalar

  function Quasi_Random_Get_Array(quasiSequenceObject,quasiSequenceDimension,reset,quasiSequenceType) result (quasiSequencePoint)
    !% Returns an array giving a quasi-random points in a {\tt quasiSequenceDimension}-dimensional space.
    implicit none
    type(fgsl_qrng),      intent(inout)                     :: quasiSequenceObject
    integer(kind=4),      intent(in)                        :: quasiSequenceDimension
    logical,              intent(inout), optional           :: reset
    type(fgsl_qrng_type), intent(in),    optional           :: quasiSequenceType
    double precision,     dimension(quasiSequenceDimension) :: quasiSequencePoint
    logical                                                 :: resetActual
    integer                                                 :: status
    type(fgsl_qrng_type)                                    :: quasiSequenceTypeActual

    ! Determine what type of quasi-random sequence to use.
    if (present(quasiSequenceType)) then
       quasiSequenceTypeActual=quasiSequenceType
    else
       quasiSequenceTypeActual=FGSL_qRng_Sobol
    end if
    
    ! Determine if we need to reset.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    else
       resetActual=.false.
    end if
    if (resetActual.or..not.FGSL_Well_Defined(quasiSequenceObject)) quasiSequenceObject=FGSL_qRng_Alloc(quasiSequenceTypeActual&
         &,quasiSequenceDimension)
    status=FGSL_qRng_Get(quasiSequenceObject,quasiSequencePoint)
    return
  end function Quasi_Random_Get_Array

  subroutine Quasi_Random_Free(quasiSequenceObject)
    !% Frees a quasi-random sequence object.
    implicit none
    type(fgsl_qrng), intent(inout) :: quasiSequenceObject
    
    call FGSL_qRng_Free(quasiSequenceObject)
    return
  end subroutine Quasi_Random_Free
  
end module Quasi_Random
