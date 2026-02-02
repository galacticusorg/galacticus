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
Contains a module which implements multidimensional minimizers.
!!}

! Specify an explicit dependence on the interface.GSL.C.multidimensional_minimization.o object file.
!: $(BUILDPATH)/interface.GSL.C.multidimensional_minimization.o

! Add dependency on GSL library.
!; gsl

module Multidimensional_Minimizer
  !!{
  Implements multidimensional minimizers.
  !!}
  use, intrinsic :: ISO_C_Binding   , only : c_ptr          , c_funptr, c_double, c_size_t, &
       &                                     c_null_ptr     , c_int
  use            :: Resource_Manager, only : resourceManager
  implicit none
  private
  public :: multiDMinimizer

  type :: gslMinimizerWrapper
     !!{
     Wrapper class for managing GSL minimizers.
     !!}
     type   (c_ptr) :: gsl        =c_null_ptr
     logical        :: useGradient
   contains
     final :: gslMinimizerWrapperDestructor
  end type gslMinimizerWrapper
  
  type :: gslMinimizerFWrapper
     !!{
     Wrapper class for managing GSL minimizer F-functions.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: gslMinimizerFWrapperDestructor
  end type gslMinimizerFWrapper
  
  type :: gslMinimizerFDFWrapper
     !!{
     Wrapper class for managing GSL minimizer FDF-functions.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: gslMinimizerFDFWrapperDestructor
  end type gslMinimizerFDFWrapper
  
  type, public :: multiDMinimizer
     !!{
     Multidimensional minimizer class.
     !!}
     private
     type     (resourceManager               )                  :: minimizerManager                    , minimizerFManager, &
          &                                                        minimizerFDFManager
     type     (gslMinimizerWrapper           ), pointer         :: minimizer_                 => null()
     type     (gslMinimizerFWrapper          ), pointer         :: minimizerFFunction         => null()
     type     (gslMinimizerFDFWrapper        ), pointer         :: minimizerFDFFunction       => null()
     type     (c_ptr                         )                  :: gslMinimizerType
     integer                                                    :: minimizerType
     integer  (c_size_t                      )                  :: countDimensions
     procedure(gslMultiminFunctionFTemplate  ), pointer, nopass :: minimizeFunction
     procedure(gslMultiminFunctionDFTemplate ), pointer, nopass :: minimizeFunctionDerivative
     procedure(gslMultiminFunctionFDFTemplate), pointer, nopass :: minimizeFunctionBoth
   contains
     !![
     <methods>
       <method description="Set the initial state of the minimizer."                                  method="set"         />
       <method description="Iterate the minimizer."                                                   method="iterate"     />
       <method description="Test the gradient of the function at the current point in the minimizer." method="testGradient"/>
       <method description="Test the size of the interval in the minimizer."                          method="testSize"    />
       <method description="Retrieve the parameter values at the current minimum from the minimizer." method="minimum"     />
       <method description="Retrieve the function value at the current minimum from the minimizer."   method="x"           />
     </methods>
     !!]
     procedure ::                  multiDMinimizerAssign
     generic   :: assignment(=) => multiDMinimizerAssign
     procedure :: set           => multiDMinimizerSet
     procedure :: iterate       => multiDMinimizerIterate
     procedure :: testGradient  => multiDMinimizerTestGradient
     procedure :: testSize      => multiDMinimizerTestSize
     procedure :: minimum       => multiDMinimizerMinimum
     procedure :: x             => multiDMinimizerX
  end type multiDMinimizer

  interface multiDMinimizer
     !!{
     Interface to multidimensional minimizer constructors.
     !!}
     module procedure multiDMinimizerConstructor
  end interface multiDMinimizer

  abstract interface
     function gslMultiminFunctionFTemplate(x)
       !!{
       Interface for {\normalfont \ttfamily gslMultiminFunctionF} type.
       !!}
       double precision                              :: gslMultiminFunctionFTemplate
       double precision, intent(in   ), dimension(:) :: x
     end function gslMultiminFunctionFTemplate

     function gslMultiminFunctionDFTemplate(x)
       !!{
       Interface for {\normalfont \ttfamily gslMultiminFunctionDF} type.
       !!}
       double precision, intent(in   ), dimension(     : ) :: x
       double precision               , dimension(size(x)) :: gslMultiminFunctionDFTemplate
     end function gslMultiminFunctionDFTemplate

     subroutine gslMultiminFunctionFDFTemplate(x,f,g)
       !!{
       Interface for {\normalfont \ttfamily gslMultiminFunctionFDF} type.
       !!}
       double precision, intent(in   ), dimension(     : ) :: x
       double precision, intent(  out)                     :: f
       double precision, intent(  out), dimension(size(x)) :: g
     end subroutine gslMultiminFunctionFDFTemplate
  end interface

  interface
     subroutine gslMultiminFunctionFdFDestructor(f) bind(c,name="gslMultiminFunctionFdFDestructor")
       !!{
       Interface to a C function which destroys a {\normalfont \ttfamily gslFunction} type.
       !!}
       import c_funptr
       type(c_funptr), value :: f
     end subroutine gslMultiminFunctionFdFDestructor

     subroutine gslMultiminFunctionFDestructor(f) bind(c,name="gslMultiminFunctionFDestructor")
       !!{
       Interface to a C function which destroys a {\normalfont \ttfamily gslFunction} type.
       !!}
       import c_funptr
       type(c_funptr), value :: f
     end subroutine gslMultiminFunctionFDestructor

     function gslMultiminFunctionFdFConstructor(n,f,df,fdf) bind(c,name="gslMultiminFunctionFdFConstructor")
       !!{
       Interface to a C function which establishes a {\normalfont \ttfamily gslMultiminFunctionFdF} type.
       !!}
       import c_ptr, c_funptr, c_size_t
       type   (c_ptr   )        :: gslMultiminFunctionFdFConstructor
       integer(c_size_t), value :: n
       type   (c_funptr), value :: f                                , df, &
            &                      fdf
     end function gslMultiminFunctionFdFConstructor

     function gsl_multimin_fdfminimizer_type_get(i) bind(c,name='gsl_multimin_fdfminimizer_type_get')
       !!{
       Interface to C function which gets a GSL minimizer type.
       !!}
       import c_ptr, c_int
       type   (c_ptr)        :: gsl_multimin_fdfminimizer_type_get
       integer(c_int), value :: i
     end function gsl_multimin_fdfminimizer_type_get
     
     function gsl_multimin_fdfminimizer_alloc(T,n) bind(c,name='gsl_multimin_fdfminimizer_alloc')
       !!{
       Template for the GSL multimin allocate function.
       !!}
       import c_ptr, c_size_t
       type   (c_ptr   )        :: gsl_multimin_fdfminimizer_alloc
       type   (c_ptr   ), value :: T
       integer(c_size_t), value :: n
     end function gsl_multimin_fdfminimizer_alloc

     function gsl_multimin_fdfminimizer_set(s,fdf,x,step_size,tol) bind(c,name='gsl_multimin_fdfminimizer_set')
       !!{
       Template for the GSL multimin set function.
       !!}
       import c_ptr, c_int, c_double
       integer(c_int   )        :: gsl_multimin_fdfminimizer_set
       type   (c_ptr   ), value :: s                            , fdf, &
            &                      x
       real   (c_double), value :: step_size                    , tol
     end function gsl_multimin_fdfminimizer_set

     subroutine gsl_multimin_fdfminimizer_free(s) bind(c,name='gsl_multimin_fdfminimizer_free')
       !!{
       Template for the GSL multimin set function.
       !!}
       import c_ptr
       type(c_ptr), value :: s
     end subroutine gsl_multimin_fdfminimizer_free

     function gsl_multimin_fdfminimizer_iterate(s) bind(c,name='gsl_multimin_fdfminimizer_iterate')
       !!{
       Template for the GSL multimin iterator function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_multimin_fdfminimizer_iterate
       type   (c_ptr), value :: s
     end function gsl_multimin_fdfminimizer_iterate

     function gsl_multimin_fdfminimizer_x(s) bind(c,name='gsl_multimin_fdfminimizer_x')
       !!{
       Template for the GSL multimin x function.
       !!}
       import c_ptr
       type(c_ptr)        :: gsl_multimin_fdfminimizer_x
       type(c_ptr), value :: s
     end function gsl_multimin_fdfminimizer_x

     function gsl_multimin_fdfminimizer_minimum(s) bind(c,name='gsl_multimin_fdfminimizer_minimum')
       !!{
       Template for the GSL multimin minimum function.
       !!}
       import c_ptr, c_double
       real(c_double)        :: gsl_multimin_fdfminimizer_minimum
       type(c_ptr   ), value :: s
     end function gsl_multimin_fdfminimizer_minimum

     function gsl_multimin_fdfminimizer_gradient(s) bind(c,name='gsl_multimin_fdfminimizer_gradient')
       !!{
       Template for the GSL multimin gradient function.
       !!}
       import c_ptr, c_int, c_double
       type(c_ptr)        :: gsl_multimin_fdfminimizer_gradient
       type(c_ptr), value :: s
     end function gsl_multimin_fdfminimizer_gradient

     function gsl_multimin_fdfminimizer_dx(s) bind(c,name='gsl_multimin_fdfminimizer_dx')
       !!{
       Template for the GSL multimin size function.
       !!}
       import c_ptr, c_int, c_double
       type(c_ptr)        :: gsl_multimin_fdfminimizer_dx
       type(c_ptr), value :: s
     end function gsl_multimin_fdfminimizer_dx

     function gsl_multimin_test_gradient(g,epsabs) bind(c,name='gsl_multimin_test_gradient')
       !!{
       Template for the GSL multimin minimum function.
       !!}
       import c_ptr, c_int, c_double
       integer(c_int   )        :: gsl_multimin_test_gradient
       type   (c_ptr   ), value :: g
       real   (c_double), value :: epsabs
     end function gsl_multimin_test_gradient

     function gsl_multimin_test_size(size,epsabs) bind(c,name='gsl_multimin_test_size')
       !!{
       Template for the GSL multimin test size function.
       !!}
       import c_ptr, c_int, c_double
       integer(c_int   )        :: gsl_multimin_test_size
       real   (c_double), value :: size                  , epsabs
     end function gsl_multimin_test_size

     function gslMultiminFunctionFConstructor(n,f) bind(c,name="gslMultiminFunctionFConstructor")
       !!{
       Interface to a C function which establishes a {\normalfont \ttfamily gslMultiminFunctionF} type.
       !!}
       import c_ptr, c_funptr, c_size_t
       type   (c_ptr   )        :: gslMultiminFunctionFConstructor
       integer(c_size_t), value :: n
       type   (c_funptr), value :: f
     end function gslMultiminFunctionFConstructor

     function gsl_multimin_fminimizer_type_get(i) bind(c,name='gsl_multimin_fminimizer_type_get')
       !!{
       Interface to C function which gets a GSL minimizer type.
       !!}
       import c_ptr, c_int
       type   (c_ptr)        :: gsl_multimin_fminimizer_type_get
       integer(c_int), value :: i
     end function gsl_multimin_fminimizer_type_get
     
     function gsl_multimin_fminimizer_alloc(T,n) bind(c,name='gsl_multimin_fminimizer_alloc')
       !!{
       Template for the GSL multimin allocate function.
       !!}
       import c_ptr, c_size_t
       type   (c_ptr   )        :: gsl_multimin_fminimizer_alloc
       type   (c_ptr   ), value :: T
       integer(c_size_t), value :: n
     end function gsl_multimin_fminimizer_alloc

     function gsl_multimin_fminimizer_set(s,f,x,step_size) bind(c,name='gsl_multimin_fminimizer_set')
       !!{
       Template for the GSL multimin set function.
       !!}
       import c_ptr, c_int, c_double
       integer(c_int   )        :: gsl_multimin_fminimizer_set
       type   (c_ptr   ), value :: s                            , f        , &
            &                      x                            , step_size
     end function gsl_multimin_fminimizer_set

     subroutine gsl_multimin_fminimizer_free(s) bind(c,name='gsl_multimin_fminimizer_free')
       !!{
       Template for the GSL multimin set function.
       !!}
       import c_ptr
       type(c_ptr), value :: s
     end subroutine gsl_multimin_fminimizer_free

     function gsl_multimin_fminimizer_iterate(s) bind(c,name='gsl_multimin_fminimizer_iterate')
       !!{
       Template for the GSL multimin iterator function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_multimin_fminimizer_iterate
       type   (c_ptr), value :: s
     end function gsl_multimin_fminimizer_iterate

     function gsl_multimin_fminimizer_x(s) bind(c,name='gsl_multimin_fminimizer_x')
       !!{
       Template for the GSL multimin x function.
       !!}
       import c_ptr
       type(c_ptr)        :: gsl_multimin_fminimizer_x
       type(c_ptr), value :: s
     end function gsl_multimin_fminimizer_x

     function gsl_multimin_fminimizer_minimum(s) bind(c,name='gsl_multimin_fminimizer_minimum')
       !!{
       Template for the GSL multimin minimum function.
       !!}
       import c_ptr, c_double
       real(c_double)        :: gsl_multimin_fminimizer_minimum
       type(c_ptr   ), value :: s
     end function gsl_multimin_fminimizer_minimum

     function gsl_multimin_fminimizer_size(s) bind(c,name='gsl_multimin_fminimizer_size')
       !!{
       Template for the GSL multimin size function.
       !!}
       import c_ptr, c_int, c_double
       real(c_double)        :: gsl_multimin_fminimizer_size
       type(c_ptr   ), value :: s
     end function gsl_multimin_fminimizer_size

  end interface

  ! Minimizer types.
  integer, public, parameter :: gsl_multimin_fdfminimizer_conjugate_fr    =1
  integer, public, parameter :: gsl_multimin_fdfminimizer_conjugate_pr    =2
  integer, public, parameter :: gsl_multimin_fdfminimizer_vector_bfgs2    =3
  integer, public, parameter :: gsl_multimin_fdfminimizer_vector_bfgs     =4
  integer, public, parameter :: gsl_multimin_fdfminimizer_steepest_descent=5
  integer, public, parameter :: gsl_multimin_fminimizer_nmsimplex2        =6
  integer, public, parameter :: gsl_multimin_fminimizer_nmsimplex2rand    =7

  ! Sub-module scope pointer to self for use in wrapper functions.
  class(multiDMinimizer), pointer :: self_
  !$omp threadprivate(self_)
  
contains

  function multiDMinimizerConstructor(countDimensions,minimizeFunction,minimizeFunctionDerivative,minimizeFunctionBoth,minimizerType) result(self)
    !!{
    Constructor for {\normalfont \ttfamily multiDMinimizer} class.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_funloc, c_loc
    use            :: Error        , only : Error_Report
    implicit none
    type     (multiDMinimizer               )                          :: self
    integer  (c_size_t                      ), intent(in   )           :: countDimensions
    procedure(gslMultiminFunctionFTemplate  )                          :: minimizeFunction
    procedure(gslMultiminFunctionDFTemplate )               , optional :: minimizeFunctionDerivative
    procedure(gslMultiminFunctionFDFTemplate)               , optional :: minimizeFunctionBoth
    integer                                  , intent(in   ), optional :: minimizerType
    integer                                                            :: minimizerType_
    class    (*                             ), pointer                 :: dummyPointer_
    !![
    <constructorAssign variables="countDimensions, *minimizeFunction, *minimizeFunctionDerivative, *minimizeFunctionBoth"/>
    !!]

    ! Determine if a gradient function is supplied.
    allocate(self%minimizer_)
    self%minimizer_%useGradient=present(minimizeFunctionDerivative)
    if (self%minimizer_%useGradient) then
       if (.not.present(minimizeFunctionBoth)) call Error_Report('a function evaluating both the function and its derivatives must be supplied'//{introspection:location})
       ! Get the interpolator type.
       minimizerType_=gsl_multimin_fdfminimizer_conjugate_pr
       if (present(minimizerType)) minimizerType_=minimizerType
       if (minimizerType_ >  gsl_multimin_fdfminimizer_steepest_descent) call Error_Report('this minimizer type is not intended for cases where a gradient function is available'//{introspection:location})
       self%minimizerType   =                                           minimizerType_
       self%gslMinimizerType=gsl_multimin_fdfminimizer_type_get(        minimizerType_                )
       ! Allocate the minimizer object.
       self%minimizer_          %gsl=gsl_multimin_fdfminimizer_alloc   (self%gslMinimizerType ,countDimensions)
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_          => self%minimizer_
       self%minimizerManager  =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
       ! Build the function object.
       allocate(self%minimizerFDFFunction)
       self%minimizerFDFFunction%gsl=gslMultiminFunctionFdFConstructor(countDimensions,c_funloc(gslMultiminFunctionFWrapper),c_funloc(gslMultiminFunctionDFWrapper),c_funloc(gslMultiminFunctionFDFWrapper))
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_            => self%minimizerFDFFunction
       self%minimizerFDFManager =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
    else
       if (     present(minimizeFunctionBoth)) call Error_Report('no function evaluating both the function and its derivatives should be supplied'//{introspection:location})
       ! Get the interpolator type.
       minimizerType_=gsl_multimin_fminimizer_nmsimplex2
       if (present(minimizerType)) minimizerType_=minimizerType
       if (minimizerType_ <  gsl_multimin_fminimizer_nmsimplex2) call Error_Report('this minimizer type is intended for cases where a gradient function is available'//{introspection:location})
       self%minimizerType   =                                         minimizerType_
       self%gslMinimizerType=gsl_multimin_fminimizer_type_get(        minimizerType_                )
       ! Allocate the minimizer object.
       self%minimizer_        %gsl=gsl_multimin_fminimizer_alloc   (self%gslMinimizerType ,countDimensions)
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_          => self%minimizer_
       self%minimizerManager  =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
       ! Build the function object.
       allocate(self%minimizerFFunction)
       self%minimizerFFunction%gsl=gslMultiminFunctionFConstructor(countDimensions,c_funloc(gslMultiminFunctionFWrapper))
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_          => self%minimizerFFunction
       self%minimizerFManager =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
    end if
    return
  end function multiDMinimizerConstructor

  subroutine gslMinimizerWrapperDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily gslMinimizerWrapper} object.
    !!}
    implicit none
    type(gslMinimizerWrapper), intent(inout) :: self

    if (self%useGradient) then
       call gsl_multimin_fdfminimizer_free(self%gsl)
    else
       call gsl_multimin_fminimizer_free  (self%gsl)
    end if
    return
  end subroutine gslMinimizerWrapperDestructor

  subroutine gslMinimizerFWrapperDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily gslMinimizerFWrapper} object.
    !!}
    implicit none
    type(gslMinimizerFWrapper), intent(inout) :: self

    call gslMultiminFunctionFDestructor(self%gsl)
    return
  end subroutine gslMinimizerFWrapperDestructor

  subroutine gslMinimizerFDFWrapperDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily gslMinimizerFDFWrapper} object.
    !!}
    implicit none
    type(gslMinimizerFDFWrapper), intent(inout) :: self

    call gslMultiminFunctionFdFDestructor(self%gsl)
    return
  end subroutine gslMinimizerFDFWrapperDestructor

  subroutine multiDMinimizerAssign(to,from)
    !!{
    Assignment operator for \refClass{multiDMinimizer} objects.
    !!}
    implicit none
    class(multiDMinimizer), intent(  out) :: to
    class(multiDMinimizer), intent(in   ) :: from

    to%minimizer_                 => from%minimizer_
    to%minimizerFFunction         => from%minimizerFFunction
    to%minimizerFDFFunction       => from%minimizerFDFFunction
    to%minimizeFunction           => from%minimizeFunction
    to%minimizeFunctionDerivative => from%minimizeFunctionDerivative
    to%minimizeFunctionBoth       => from%minimizeFunctionBoth
    to%minimizerManager           =  from%minimizerManager
    to%minimizerFManager          =  from%minimizerFManager
    to%minimizerFDFManager        =  from%minimizerFDFManager
    to%gslMinimizerType           =  from%gslMinimizerType
    to%minimizerType              =  from%minimizerType
    to%countDimensions            =  from%countDimensions
    return
  end subroutine multiDMinimizerAssign
  
  subroutine multiDMinimizerSet(self,x,stepSize,tolerance)
    !!{
    Set the initial state for a {\normalfont \ttfamily multiDMinimizer} object.
    !!}
    use :: Linear_Algebra, only : vector
    use :: Error         , only : Error_Report
    use :: Interface_GSL , only : GSL_Success , GSL_Failure
    implicit none
    class           (multiDMinimizer), intent(inout), target        :: self
    double precision                 , intent(in   ), dimension(:)  :: x
    double precision                 , intent(in   ), dimension(..) :: stepSize
    double precision                 , intent(in   ), optional      :: tolerance
    type            (vector         )                               :: vector_  , stepSize_
    integer         (c_int          )                               :: status

    self_ => self
    vector_=vector(x)
    status=GSL_Failure
    if (self%minimizer_%useGradient) then
       if (.not.present(tolerance)) call Error_Report('a tolerance must be provided'   //{introspection:location})
       select rank (stepSize)
       rank (0)
          status =GSL_MultiMin_FDFMinimizer_Set(self%minimizer_%gsl,self%minimizerFDFFunction%gsl,vector_%gslObject(),stepSize,tolerance)
       rank default
          call Error_Report('stepSize should be a scalar'//{introspection:location})
       end select
    else
       if (     present(tolerance)) call Error_Report('no tolerance should be provided'//{introspection:location})
       select rank (stepSize)
       rank (1)
          stepSize_=vector(stepSize)
          status =GSL_MultiMin_FMinimizer_Set  (self%minimizer_%gsl,self%minimizerFFunction  %gsl,vector_%gslObject(),stepSize_%gslObject())
       rank default
          call Error_Report('stepSize should be a rank-1 array'//{introspection:location})
       end select
    end if
    if (status /= GSL_Success) call Error_Report('failed to set minimizer'//{introspection:location})
    self_ => null()
    return
  end subroutine multiDMinimizerSet
  
  subroutine multiDMinimizerIterate(self,status)
    !!{
    Iterate a multidimensional minimizer.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class  (multiDMinimizer), intent(inout), target   :: self
    integer(c_int          ), intent(  out), optional :: status
    integer(c_int          )                          :: status_

    self_ => self
    if (self%minimizer_%useGradient) then
       status_=gsl_multimin_fdfminimizer_iterate(self%minimizer_%gsl)
    else
       status_=gsl_multimin_fminimizer_iterate  (self%minimizer_%gsl)
    end if
    if (present(status)) status=status_
    if (status /= GSL_Success .and. .not.present(status)) call Error_Report('failed to iterate minimizer'//{introspection:location})
    self_ => null()
    return
  end subroutine multiDMinimizerIterate

  logical function multiDMinimizerTestGradient(self,toleranceAbsolute)
    !!{
    Test for convergence in the gradient of the function for a multidimensional minimizer.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class           (multiDMinimizer), intent(inout) :: self
    double precision                 , intent(in   ) :: toleranceAbsolute

    if (.not.self%minimizer_%useGradient) call Error_Report('gradient test can not be performed for this minimizer'//{introspection:location})
    multiDMinimizerTestGradient=GSL_Multimin_Test_Gradient(GSL_Multimin_FDFMinimizer_Gradient(self%minimizer_%gsl),toleranceAbsolute) == GSL_Success
    return
  end function multiDMinimizerTestGradient
  
  logical function multiDMinimizerTestSize(self,toleranceAbsolute)
    !!{
    Test for convergence in the interval size for a multidimensional minimizer.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class           (multiDMinimizer), intent(inout) :: self
    double precision                 , intent(in   ) :: toleranceAbsolute
    !$GLC attributes unused :: self

    if (    self%minimizer_%useGradient) call Error_Report('size test can not be performed for this minimizer'//{introspection:location})
    multiDMinimizerTestSize=GSL_Multimin_Test_Size(GSL_Multimin_FMinimizer_Size(self%minimizer_%gsl),toleranceAbsolute) == GSL_Success
    return
  end function multiDMinimizerTestSize
  
  double precision function multiDMinimizerMinimum(self)
    !!{
    Return the minimized function value for a multidimensional minimizer.
    !!}
    implicit none
    class(multiDMinimizer), intent(inout) :: self

    if (self%minimizer_%useGradient) then
       multiDMinimizerMinimum=GSL_Multimin_FDFMinimizer_Minimum(self%minimizer_%gsl)
    else
       multiDMinimizerMinimum=GSL_Multimin_FMinimizer_Minimum  (self%minimizer_%gsl)
    end if
    return
  end function multiDMinimizerMinimum

  function multiDMinimizerX(self)
    !!{
    Return the parameter values at the minimum of the minimized function value for a multidimensional minimizer.
    !!}
    use :: Linear_Algebra, only : gsl_vector_get
    implicit none
    class           (multiDMinimizer), intent(inout)                   :: self
    double precision                 , dimension(self%countDimensions) :: multiDMinimizerX
    type            (c_ptr          )                                  :: x
    integer         (c_size_t       )                                  :: i
    
    if (self%minimizer_%useGradient) then
       x=GSL_Multimin_FDFMinimizer_X(self%minimizer_%gsl)
    else
       x=GSL_Multimin_FMinimizer_X  (self%minimizer_%gsl)
    end if
    do i=1,self%countDimensions
       multiDMinimizerX(i)=GSL_Vector_Get(x,i-1)
    end do
    return
  end function multiDMinimizerX

  function gslMultiminFunctionFWrapper(x,parameters)
    !!{
    Wrapper function for minimized function value.
    !!}
    use :: Linear_Algebra, only : gsl_vector_get
    implicit none
    real   (c_double)                                   :: gslMultiminFunctionFWrapper
    type   (c_ptr   ), value                            :: x                          , parameters
    real   (c_double), dimension(self_%countDimensions) :: x_
    integer(c_size_t)                                   :: i
    !$GLC attributes unused :: parameters
    
    do i=1,self_%countDimensions
       x_(i)=GSL_Vector_Get(x,i-1)
    end do
    gslMultiminFunctionFWrapper=self_%minimizeFunction(x_)
    return
  end function gslMultiminFunctionFWrapper

  subroutine gslMultiminFunctionDFWrapper(x,parameters,g)
    !!{
    Wrapper function for minimized function gradient.
    !!}
    use :: Linear_Algebra, only : gsl_vector_get, gsl_vector_set
    implicit none
    type   (c_ptr   ), value                            :: x, parameters, &
         &                                                 g
    real   (c_double), dimension(self_%countDimensions) :: x_, g_
    integer(c_size_t)                                   :: i
    !$GLC attributes unused :: parameters

    do i=1,self_%countDimensions
       x_(i)=GSL_Vector_Get(x,i-1)
    end do
    g_=self_%minimizeFunctionDerivative(x_)
    do i=1,self_%countDimensions
       call GSL_Vector_Set(g,i-1,g_(i))
    end do
    return
  end subroutine gslMultiminFunctionDFWrapper

  subroutine gslMultiminFunctionFDFWrapper(x,parameters,f,g)
    !!{
    Wrapper function for minimized function value and gradient.
    !!}
    use :: Linear_Algebra, only : gsl_vector_get, gsl_vector_set
    implicit none
    type   (c_ptr   ), value                            :: x, parameters, &
         &                                                 g
    real   (c_double)                                   :: f
    real   (c_double), dimension(self_%countDimensions) :: x_, g_
    integer(c_size_t)                                   :: i
    !$GLC attributes unused :: parameters

    do i=1,self_%countDimensions
       x_(i)=GSL_Vector_Get(x,i-1)
    end do
    call self_%minimizeFunctionBoth(x_,f,g_)
    do i=1,self_%countDimensions
       call GSL_Vector_Set(g,i-1,g_(i))
    end do
    return
  end subroutine gslMultiminFunctionFDFWrapper

end module Multidimensional_Minimizer
