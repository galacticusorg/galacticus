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

!% Contains a module which implements quasi-random sequences.

module Quasi_Random
  !% Implements quasi-random sequences.
  use FGSL
  implicit none
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
    type            (fgsl_qrng     ), intent(inout)           :: quasiSequenceObject
    logical                         , intent(inout), optional :: reset
    type            (fgsl_qrng_type), intent(in   ), optional :: quasiSequenceType
    integer         (kind=4        ), parameter               :: nGet                      =1
    double precision                                          :: pointArray             (1)
    logical                                                   :: resetActual
    type            (fgsl_qrng_type)                          :: quasiSequenceTypeActual

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
    type            (fgsl_qrng     ), intent(inout)                               :: quasiSequenceObject
    integer         (kind=4        ), intent(in   )                               :: quasiSequenceDimension
    logical                         , intent(inout)                    , optional :: reset
    type            (fgsl_qrng_type), intent(in   )                    , optional :: quasiSequenceType
    double precision                , dimension(quasiSequenceDimension)           :: quasiSequencePoint
    logical                                                                       :: resetActual
    integer                                                                       :: status
    type            (fgsl_qrng_type)                                              :: quasiSequenceTypeActual

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
