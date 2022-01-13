!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Implements linear algebra calculations.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_ptr, c_funptr, c_double, c_size_t, c_int
  implicit none
  private
  public :: multiDMinimizer

  type, public :: multiDMinimizer
     !!{
     Multidimensional minimizer class.
     !!}
     private
     type     (c_ptr                         ), allocatable         :: minimizer_                , minimizerFDFFunction
     type     (c_ptr                         )                      :: gslMinimizerType
     integer                                                        :: minimizerType
     integer  (c_size_t                      )                      :: countDimensions
     procedure(gslMultiminFunctionFTemplate  ), pointer    , nopass :: minimizeFunction
     procedure(gslMultiminFunctionDFTemplate ), pointer    , nopass :: minimizeFunctionDerivative
     procedure(gslMultiminFunctionFDFTemplate), pointer    , nopass :: minimizeFunctionBoth
   contains
     !![
     <methods>
       <method description="Set the initial state of the minimizer." method="set" />
       <method description="Iterate the minimizer." method="iterate" />
       <method description="Test the gradient of the function at the current point in the minimizer." method="testGradient" />
       <method description="Retrieve the parameter values at the current minimum from the minimizer." method="minimum" />
       <method description="Retrieve the function value at the current minimum from the minimizer." method="x" />
     </methods>
     !!]
     final     ::                 multiDMinimizerDestructor
     procedure :: set          => multiDMinimizerSet
     procedure :: iterate      => multiDMinimizerIterate
     procedure :: testGradient => multiDMinimizerTestGradient
     procedure :: minimum      => multiDMinimizerMinimum
     procedure :: x            => multiDMinimizerX
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

     subroutine gslMultiminFunctionDestructor(f) bind(c,name="gslMultiminFunctionDestructor")
       !!{
       Interface to a C function which destroys a {\normalfont \ttfamily gslFunction} type.
       !!}
       import c_funptr
       type(c_funptr), value :: f
     end subroutine gslMultiminFunctionDestructor

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
       real(c_double)        ::  gsl_multimin_fdfminimizer_minimum
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

     function gsl_multimin_test_gradient(g,epsabs) bind(c,name='gsl_multimin_test_gradient')
       !!{
       Template for the GSL multimin minimum function.
       !!}
       import c_ptr, c_int, c_double
       integer(c_int   )        :: gsl_multimin_test_gradient
       type   (c_ptr   ), value :: g
       real   (c_double), value :: epsabs
     end function gsl_multimin_test_gradient
  end interface

  ! Minimizer types.
  integer, public, parameter :: gsl_multimin_fdfminimizer_conjugate_fr    =1
  integer, public, parameter :: gsl_multimin_fdfminimizer_conjugate_pr    =2
  integer, public, parameter :: gsl_multimin_fdfminimizer_vector_bfgs2    =3
  integer, public, parameter :: gsl_multimin_fdfminimizer_vector_bfgs     =4
  integer, public, parameter :: gsl_multimin_fdfminimizer_steepest_descent=5

  ! Sub-module scope pointer to self for use in wrapper functions.
  class(multiDMinimizer), pointer :: self_
  !$omp threadprivate(self_)
  
contains

  function multiDMinimizerConstructor(countDimensions,minimizeFunction,minimizeFunctionDerivative,minimizeFunctionBoth,minimizerType) result(self)
    !!{
    Constructor for {\normalfont \ttfamily multiDMinimizer} class.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_funloc
    implicit none
    type     (multiDMinimizer               )                          :: self
    integer  (c_size_t                      ), intent(in   )           :: countDimensions
    procedure(gslMultiminFunctionFTemplate  )                          :: minimizeFunction
    procedure(gslMultiminFunctionDFTemplate )                          :: minimizeFunctionDerivative
    procedure(gslMultiminFunctionFDFTemplate)                          :: minimizeFunctionBoth
    integer                                  , intent(in   ), optional :: minimizerType
    !![
    <optionalArgument name="minimizerType" defaultsTo="gsl_multimin_fdfminimizer_conjugate_pr"/>
    <constructorAssign variables="countDimensions, *minimizeFunction, *minimizeFunctionDerivative, *minimizeFunctionBoth"/>
    !!]
    
    ! Get the interpolator type.
    self%minimizerType   =                                           minimizerType_
    self%gslMinimizerType=gsl_multimin_fdfminimizer_type_get(        minimizerType_                )
    ! Allocate the minimizer object.
    allocate(self%minimizer_)
    self%minimizer_      =gsl_multimin_fdfminimizer_alloc   (self%gslMinimizerType ,countDimensions)
    ! Build the function object.
    allocate(self%minimizerFDFFunction)
    self%minimizerFDFFunction=gslMultiminFunctionFdFConstructor(countDimensions,c_funloc(gslMultiminFunctionFWrapper),c_funloc(gslMultiminFunctionDFWrapper),c_funloc(gslMultiminFunctionFDFWrapper))
    return
  end function multiDMinimizerConstructor

  subroutine multiDMinimizerDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily multiDMinimizer} class.
    !!}
    implicit none
    type(multiDMinimizer), intent(inout) :: self

    if (allocated(self%minimizer_)) then
       call gsl_multimin_fdfminimizer_free(self%minimizer_)
       deallocate(self%minimizer_)
    end if
    if (allocated(self%minimizerFDFFunction)) then
       call gslMultiminFunctionDestructor(self%minimizerFDFFunction)
       deallocate(self%minimizerFDFFunction)
    end if
    return
  end subroutine multiDMinimizerDestructor
  
  subroutine multiDMinimizerSet(self,x,stepSize,tolerance)
    !!{
    Set the initial state for a {\normalfont \ttfamily multiDMinimizer} object.
    !!}
    use :: Linear_Algebra  , only : vector
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Interface_GSL   , only : GSL_Success
    implicit none
    class           (multiDMinimizer), intent(inout), target       :: self
    double precision                 , intent(in   ), dimension(:) :: x
    double precision                 , intent(in   )               :: stepSize, tolerance
    type            (vector         )                              :: vector_
    integer         (c_int          )                              :: status

    self_ => self
    vector_=vector(x)
    status =GSL_MultiMin_FDFMinimizer_Set(self%minimizer_,self%minimizerFDFFunction,vector_%gslObject(),stepSize,tolerance)
    if (status /= GSL_Success) call Galacticus_Error_Report('failed to set minimizer'//{introspection:location})
    self_ => null()
    return
  end subroutine multiDMinimizerSet
  
  subroutine multiDMinimizerIterate(self,status)
    !!{
    Iterate a multidimensional minimizer.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Interface_GSL   , only : GSL_Success
    implicit none
    class  (multiDMinimizer), intent(inout), target   :: self
    integer(c_int          ), intent(  out), optional :: status
    integer(c_int          )                          :: status_

    self_ => self
    status_=gsl_multimin_fdfminimizer_iterate(self%minimizer_)
    if (present(status)) status=status_
    if (status /= GSL_Success .and. .not.present(status)) call Galacticus_Error_Report('failed to iterate minimizer'//{introspection:location})
    self_ => null()
    return
  end subroutine multiDMinimizerIterate

  logical function multiDMinimizerTestGradient(self,toleranceAbsolute)
    !!{
    Test for convergence in the gradient of the function for a multidimensional minimizer.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class           (multiDMinimizer), intent(inout) :: self
    double precision                 , intent(in   ) :: toleranceAbsolute
    !$GLC attributes unused :: self

    multiDMinimizerTestGradient=GSL_Multimin_Test_Gradient(GSL_Multimin_FDFMinimizer_Gradient(self%minimizer_),toleranceAbsolute) == GSL_Success
    return
  end function multiDMinimizerTestGradient

  double precision function multiDMinimizerMinimum(self)
    !!{
    Return the minimized function value for a multidimensional minimizer.
    !!}
    implicit none
    class(multiDMinimizer), intent(inout) :: self

    multiDMinimizerMinimum=GSL_Multimin_FDFMinimizer_Minimum(self%minimizer_)
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
    
    x=GSL_Multimin_FDFMinimizer_X(self%minimizer_)
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
