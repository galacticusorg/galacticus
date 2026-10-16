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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which implements quasi-random sequences.
!!}

! Specify an explicit dependence on the interface.GSL.C.quasi_random.o object file.
!: $(BUILDPATH)/interface/GSL/C/quasi_random.o

! Add dependency on GSL library.
!; gsl

module Numerical_Quasi_Random_Sequences
  !!{RST
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

  ! Maximum dimensions supported by each sequence type. GSL does not enforce these limits itself: for a Sobol sequence of
  ! more than 40 dimensions `gsl_qrng_alloc` returns a generator without error, but the sequence it produces is invalid.
  integer, parameter, dimension(4) :: dimensionMaximum=[12,40,1229,1229]

  interface
     function gsl_qrng_alloc(T,d) bind(c,name='gsl_qrng_alloc')
       !!{RST
       Template for the GSL quasi-random number generator allocator function.
       !!}
       import c_ptr, c_int
       type   (c_ptr)        :: gsl_qrng_alloc
       type   (c_ptr), value :: T
       integer(c_int), value :: d
     end function gsl_qrng_alloc
     subroutine gsl_qrng_free(q) bind(c,name='gsl_qrng_free')
       !!{RST
       Template for the GSL quasi-random number generator free function.
       !!}
       import c_ptr
       type(c_ptr), value :: q
     end subroutine gsl_qrng_free
     function gsl_qrng_get(q,x) bind(c,name='gsl_qrng_get')
       !!{RST
       Template for the GSL quasi-random number generator get function.
       !!}
       import c_ptr, c_double, c_int
       integer(c_int   )               :: gsl_qrng_get
       type   (c_ptr   ), value        :: q
       real   (c_double), dimension(*) :: x
     end function gsl_qrng_get
     function gsl_qrng_type_get(i) bind(c,name='gsl_qrng_type_get')
       !!{RST
       Template for GSL interface quasi-random number generator type function.
       !!}
       import c_ptr, c_int
       type   (c_ptr)                       :: gsl_qrng_type_get
       integer(c_int), intent(in   ), value :: i
     end function gsl_qrng_type_get
  end interface

  type :: gslQRNGWrapper
     !!{RST
     Wrapper class for managing GSL quasi-random number generators.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: gslQRNGWrapperDestructor
  end type gslQRNGWrapper
  
  type :: quasiRandomNumberGenerator
     !!{RST
     Type providing quasi-random number generators.
     !!}
     private
     type   (resourceManager)              :: qrngManager
     type   (gslQRNGWrapper ), pointer     :: gsl_qrng       => null()
     type   (c_ptr          ), allocatable :: gsl_qrng_type
     integer                               :: qrngType                , countDimensions_
   contains
     !![
     <methods docformat="rst">
       <method description="Get the next number from a one-dimensional sequence." method="get"       />
       <method description="Get the next point from a sequence of any dimension." method="getVector" />
       <method description="Return the dimension of the sequence."                method="dimensions"/>
     </methods>
     !!]
     procedure :: get        => quasiRandomNumberGeneratorGet
     procedure :: getVector  => quasiRandomNumberGeneratorGetVector
     procedure :: dimensions => quasiRandomNumberGeneratorDimensions
  end type quasiRandomNumberGenerator
  
  interface quasiRandomNumberGenerator
     !!{RST
     Constructor for the ``quasiRandomNumberGenerator`` class.
     !!}
     module procedure quasiRandomNumberGeneratorConstructor
  end interface quasiRandomNumberGenerator
  
contains

  function quasiRandomNumberGeneratorConstructor(qrngType,countDimensions) result(self)
    !!{RST
    Constructor for ``quasiRandomNumberGenerator`` objects. The sequence is one-dimensional unless ``countDimensions`` is given.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : var_str     , operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    type   (quasiRandomNumberGenerator)                          :: self
    integer                            , intent(in   ), optional :: qrngType     , countDimensions
    class  (*                         ), pointer                 :: dummyPointer_
    !![
    <optionalArgument name="qrngType"        defaultsTo="gsl_qrng_sobol"/>
    <optionalArgument name="countDimensions" defaultsTo="1"             />
    !!]
    
    ! Validate the sequence type and dimension.
    if (qrngType_ < 1 .or. qrngType_ > size(dimensionMaximum))                                   &
         & call Error_Report('unknown quasi-random sequence type'//{introspection:location})
    if (countDimensions_ < 1 .or. countDimensions_ > dimensionMaximum(qrngType_))                &
         & call Error_Report(                                                                    &
         &                   var_str('quasi-random sequence dimension must be between 1 and ')// &
         &                   dimensionMaximum(qrngType_)                                      // &
         &                   ' for this sequence type'                                        // &
         &                   {introspection:location}                                            &
         &                  )
    ! Get the interpolator type.
    self%qrngType        =qrngType_
    self%countDimensions_=countDimensions_
    allocate(self%gsl_qrng_type)
    self%gsl_qrng_type=gsl_qrng_type_get(qrngType_)
    ! Allocate the sequence.
    allocate(self%gsl_qrng)
    self%gsl_qrng%gsl=gsl_qrng_alloc (self%gsl_qrng_type,countDimensions_)
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi?id=105807" docformat="rst">
      <description>
      ICE when passing a derived type component to a class(*) function argument.
      </description>
    !!]
    dummyPointer_    => self%gsl_qrng
    self%qrngManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    return
  end function quasiRandomNumberGeneratorConstructor

  subroutine gslQRNGWrapperDestructor(self)
    !!{RST
    Destroy a ``gslQRNGWrapper`` object.
    !!}
    implicit none
    type(gslQRNGWrapper), intent(inout) :: self

    call gsl_qrng_free(self%gsl)
    return
  end subroutine gslQRNGWrapperDestructor

  double precision function quasiRandomNumberGeneratorGet(self)
    !!{RST
    Return the next entry in the quasi-random sequence.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class           (quasiRandomNumberGenerator), intent(inout) :: self
    double precision                            , dimension(1)  :: sequenceNext
    integer         (c_int                     )                :: status

    if (self%countDimensions_ /= 1) call Error_Report('`get` requires a one-dimensional sequence - use `getVector` instead'//{introspection:location})
    status=GSL_qRng_Get(self%gsl_qrng%gsl,sequenceNext)
    if (status /= GSL_Success) call Error_Report('failed to get next entry in quasi-random sequence'//{introspection:location})
    quasiRandomNumberGeneratorGet=sequenceNext(1)
    return
  end function quasiRandomNumberGeneratorGet

  subroutine quasiRandomNumberGeneratorGetVector(self,sequenceNext)
    !!{RST
    Return the next point in the quasi-random sequence. Note that, as implemented by GSL, the sequence omits its first point
    (the origin), so for a Sobol sequence the first :math:`2^m-1` points returned here, together with the origin, form a
    balanced :math:`2^m`-point set.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class           (quasiRandomNumberGenerator), intent(inout)               :: self
    double precision                            , intent(  out), dimension(:) :: sequenceNext
    integer         (c_int                     )                              :: status

    if (size(sequenceNext) /= self%countDimensions_) call Error_Report('array size does not match the sequence dimension'//{introspection:location})
    status=GSL_qRng_Get(self%gsl_qrng%gsl,sequenceNext)
    if (status /= GSL_Success) call Error_Report('failed to get next entry in quasi-random sequence'//{introspection:location})
    return
  end subroutine quasiRandomNumberGeneratorGetVector

  integer function quasiRandomNumberGeneratorDimensions(self)
    !!{RST
    Return the dimension of the quasi-random sequence.
    !!}
    implicit none
    class(quasiRandomNumberGenerator), intent(in   ) :: self

    quasiRandomNumberGeneratorDimensions=self%countDimensions_
    return
  end function quasiRandomNumberGeneratorDimensions
  
end module Numerical_Quasi_Random_Sequences
