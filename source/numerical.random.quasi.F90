!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module which implements quasi-random sequences.
!!}

! Specify an explicit dependence on the interface.GSL.C.quasi_random.o object file.
!: $(BUILDPATH)/interface.GSL.C.quasi_random.o

! Add dependency on GSL library.
!; gsl

module Numerical_Quasi_Random_Sequences
  !!{
  Implements quasi-random sequences.
  !!}
  use, intrinsic :: ISO_C_Binding   , only : c_ptr          , c_int, c_double, c_null_ptr
  use            :: Resource_Manager, only : resourceManager
  implicit none
  private
  public :: quasiRandomNumberGenerator

  ! Sequence types.
  integer, public, parameter :: gsl_qrng_niederreiter_2=1
  integer, public, parameter :: gsl_qrng_sobol         =2
  integer, public, parameter :: gsl_qrng_halton        =3
  integer, public, parameter :: gsl_qrng_reversehalton =4

  interface
     function gsl_qrng_alloc(T,d) bind(c,name='gsl_qrng_alloc')
       !!{
       Template for the GSL quasi-random number generator allocator function.
       !!}
       import c_ptr, c_int
       type   (c_ptr)        :: gsl_qrng_alloc
       type   (c_ptr), value :: T
       integer(c_int), value :: d
     end function gsl_qrng_alloc
     subroutine gsl_qrng_free(q) bind(c,name='gsl_qrng_free')
       !!{
       Template for the GSL quasi-random number generator free function.
       !!}
       import c_ptr
       type(c_ptr), value :: q
     end subroutine gsl_qrng_free
     function gsl_qrng_get(q,x) bind(c,name='gsl_qrng_get')
       !!{
       Template for the GSL quasi-random number generator get function.
       !!}
       import c_ptr, c_double, c_int
       integer(c_int   )               :: gsl_qrng_get
       type   (c_ptr   ), value        :: q
       real   (c_double), dimension(*) :: x
     end function gsl_qrng_get
     function gsl_qrng_type_get(i) bind(c,name='gsl_qrng_type_get')
       !!{
       Template for GSL interface quasi-random number generator type function.
       !!}
       import c_ptr, c_int
       type   (c_ptr)                       :: gsl_qrng_type_get
       integer(c_int), intent(in   ), value :: i
     end function gsl_qrng_type_get
  end interface

  type :: gslQRNGWrapper
     !!{
     Wrapper class for managing GSL quasi-random number generators.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: gslQRNGWrapperDestructor
  end type gslQRNGWrapper
  
  type :: quasiRandomNumberGenerator
     !!{
     Type providing quasi-random number generators.
     !!}
     private
     type   (resourceManager)              :: qrngManager
     type   (gslQRNGWrapper ), pointer     :: gsl_qrng      => null()
     type   (c_ptr          ), allocatable :: gsl_qrng_type
     integer                               :: qrngType
   contains
     !![
     <methods>
       <method description="Get numbers from the sequence." method="get" />
     </methods>
     !!]
     procedure :: get => quasiRandomNumberGeneratorGet
  end type quasiRandomNumberGenerator
  
  interface quasiRandomNumberGenerator
     !!{
     Constructor for the {\normalfont \ttfamily quasiRandomNumberGenerator} class.
     !!}
     module procedure quasiRandomNumberGeneratorConstructor
  end interface quasiRandomNumberGenerator
  
contains

  function quasiRandomNumberGeneratorConstructor(qrngType) result(self)
    !!{
    Constructor for {\normalfont \ttfamily quasiRandomNumberGenerator} objects.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (quasiRandomNumberGenerator)                          :: self
    integer                            , intent(in   ), optional :: qrngType
    class  (*                         ), pointer                 :: dummyPointer_
    !![
    <optionalArgument name="qrngType" defaultsTo="gsl_qrng_sobol"/>
    !!]
    
    ! Get the interpolator type.
    self%qrngType=qrngType_
    allocate(self%gsl_qrng_type)
    self%gsl_qrng_type=gsl_qrng_type_get(qrngType_)
    ! Allocate the sequence.
    allocate(self%gsl_qrng)
    self%gsl_qrng%gsl=gsl_qrng_alloc (self%gsl_qrng_type,1)
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_    => self%gsl_qrng
    self%qrngManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    return
  end function quasiRandomNumberGeneratorConstructor

  subroutine gslQRNGWrapperDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily gslQRNGWrapper} object.
    !!}
    implicit none
    type(gslQRNGWrapper), intent(inout) :: self

    call gsl_qrng_free(self%gsl)
    return
  end subroutine gslQRNGWrapperDestructor

  double precision function quasiRandomNumberGeneratorGet(self)
    !!{
    Return the next entry in the quasi-random sequence.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class           (quasiRandomNumberGenerator), intent(inout) :: self
    double precision                            , dimension(1)  :: sequenceNext
    integer         (c_int                     )                :: status

    status=GSL_qRng_Get(self%gsl_qrng%gsl,sequenceNext)
    if (status /= GSL_Success) call Error_Report('failed to get next entry in quasi-random sequence'//{introspection:location})
    quasiRandomNumberGeneratorGet=sequenceNext(1)
    return
  end function quasiRandomNumberGeneratorGet
  
end module Numerical_Quasi_Random_Sequences
