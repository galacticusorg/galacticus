!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements quasi-random sequences.

module Quasi_Random
  !% Implements quasi-random sequences.
  use, intrinsic :: ISO_C_Binding, only : c_ptr            , c_int          , c_double
  use            :: FGSL         , only : FGSL_Well_Defined, FGSL_qRng_Alloc, FGSL_qRng_Free, FGSL_qRng_Get, &
          &                               FGSL_qRng_Sobol  , fgsl_qrng      , fgsl_qrng_type
  implicit none
  private
  public :: Quasi_Random_Get, Quasi_Random_Free

  ! Sequence types.
  integer, public, parameter :: gsl_qrng_niederreiter_2=1
  integer, public, parameter :: gsl_qrng_sobol         =2
  integer, public, parameter :: gsl_qrng_halton        =3
  integer, public, parameter :: gsl_qrng_reversehalton =4

  interface
     function gsl_qrng_alloc(T,d) bind(c,name='gsl_qrng_alloc')
       !% Template for the GSL quasi-random number generator allocator function.
       import c_ptr, c_int
       type   (c_ptr)        :: gsl_qrng_alloc
       type   (c_ptr), value :: T
       integer(c_int), value :: d
     end function gsl_qrng_alloc
     subroutine gsl_qrng_free(q) bind(c,name='gsl_qrng_free')
       !% Template for the GSL quasi-random number generator free function.
       import c_ptr
       type(c_ptr), value :: q
     end subroutine gsl_qrng_free
     function gsl_qrng_get(q,x) bind(c,name='gsl_qrng_get')
       !% Template for the GSL quasi-random number generator get function.
       import c_ptr, c_double, c_int
       integer(c_int   )               :: gsl_qrng_get
       type   (c_ptr   ), value        :: q
       real   (c_double), dimension(*) :: x
     end function gsl_qrng_get
  end interface
  
  interface Quasi_Random_Get
     module procedure Quasi_Random_Get_Scalar
     module procedure Quasi_Random_Get_Array
  end interface Quasi_Random_Get

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
    !% Returns an array giving a quasi-random points in a {\normalfont \ttfamily quasiSequenceDimension}-dimensional space.
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
