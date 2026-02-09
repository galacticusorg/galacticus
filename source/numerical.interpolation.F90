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
Contains a module which acts as a simple interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library}
\href{http://www.gnu.org/software/gsl/manual/html_node/Interpolation.html}{interpolation routines}.
!!}

! Specify an explicit dependence on the interface.GSL.C.interpolation.o object file.
!: $(BUILDPATH)/interface.GSL.C.interpolation.o

! Add dependency on GSL library.
!; gsl

module Numerical_Interpolation
  !!{
  A simple interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library}
  \href{http://www.gnu.org/software/gsl/manual/html_node/Interpolation.html}{interpolation routines}.
  !!}
  use, intrinsic :: ISO_C_Binding,    only : c_ptr                           , c_size_t, c_int, c_double, &
       &                                     c_null_ptr
  use            :: Table_Labels    , only : enumerationExtrapolationTypeType
  use            :: Resource_Manager, only : resourceManager
  implicit none
  private
  public :: interpolator, interpolator2D

  ! Interpolator types.
  integer, public, parameter :: gsl_interp_linear          =1
  integer, public, parameter :: gsl_interp_polynomial      =2
  integer, public, parameter :: gsl_interp_cspline         =3
  integer, public, parameter :: gsl_interp_cspline_periodic=4
  integer, public, parameter :: gsl_interp_akima           =5
  integer, public, parameter :: gsl_interp_akima_periodic  =6

  interface
     function gsl_interp_alloc(T,size) bind(c,name='gsl_interp_alloc')
       !!{
       Template for GSL interface interpolation allocation function.
       !!}
       import c_ptr, c_size_t
       type   (c_ptr   )        :: gsl_interp_alloc
       type   (c_ptr   ), value :: T
       integer(c_size_t), value :: size
     end function gsl_interp_alloc

     function gsl_interp_init(interp,xa,ya,size) bind(c,name='gsl_interp_init')
       !!{
       Template for GSL interface interpolation initialization function.
       !!}
       import c_ptr, c_size_t, c_int, c_double
       integer(c_int   )               :: gsl_interp_init
       type   (c_ptr   ), value        :: interp
       real   (c_double), dimension(*) :: xa             , ya
       integer(c_size_t), value        :: size
     end function gsl_interp_init

     function gsl_interp_eval_e(interp,xa,ya,x,acc,y) bind(c,name='gsl_interp_eval_e')
       !!{
       Template for GSL interface interpolation function.
       !!}
       import c_ptr, c_size_t, c_int, c_double
       integer(c_int   )               :: gsl_interp_eval_e
       type   (c_ptr   ), value        :: interp
       real   (c_double), dimension(*) :: xa               , ya
       real   (c_double), value        :: x
       type   (c_ptr   ), value        :: acc
       real   (c_double)               :: y
     end function gsl_interp_eval_e
     
     function gsl_interp_eval_deriv_e(interp,xa,ya,x,acc,d) bind(c,name='gsl_interp_eval_deriv_e')
       !!{
       Template for GSL interface interpolation function.
       !!}
       import c_ptr, c_size_t, c_int, c_double
       integer(c_int   )               :: gsl_interp_eval_deriv_e
       type   (c_ptr   ), value        :: interp
       real   (c_double), dimension(*) :: xa                     , ya
       real   (c_double), value        :: x
       type   (c_ptr   ), value        :: acc
       real   (c_double)               :: d
     end function gsl_interp_eval_deriv_e
     
     subroutine gsl_interp_free(interp) bind(c,name='gsl_interp_free')
       !!{
       Template for GSL interface interpolation free function.
       !!}
       import c_ptr
       type(c_ptr), value :: interp
     end subroutine gsl_interp_free

     function gsl_interp_accel_alloc() bind(c,name='gsl_interp_accel_alloc')
       !!{
       Template for GSL interface interpolation accelerator allocation function.
       !!}
       import c_ptr
       type(c_ptr) :: gsl_interp_accel_alloc
     end function gsl_interp_accel_alloc

     function gsl_interp_accel_find(a,x_array,size,x) bind(c,name='gsl_interp_accel_find')
       !!{
       Template for GSL interface interpolation accelerator allocation function.
       !!}
       import c_ptr, c_size_t, c_double
       integer(c_size_t)               :: gsl_interp_accel_find
       type   (c_ptr   ), value        :: a
       real   (c_double), dimension(*) :: x_array
       integer(c_size_t), value        :: size
       real   (c_double), value        :: x
     end function gsl_interp_accel_find

     subroutine gsl_interp_accel_free(acc) bind(c,name='gsl_interp_accel_free')
       !!{
       Template for GSL interface interpolation accelerator free function.
       !!}
       import c_ptr
       type(c_ptr), value :: acc
     end subroutine gsl_interp_accel_free

     function gsl_interp_type_get(i) bind(c,name='gsl_interp_type_get')
       !!{
       Template for GSL interface interpolation type function.
       !!}
       import c_ptr, c_int
       type   (c_ptr)                       :: gsl_interp_type_get
       integer(c_int), intent(in   ), value :: i
     end function gsl_interp_type_get
  end interface

  !![
  <stateStorable class="interpolator">
   <interpolator>
     <restoreTo  variables="initialized" state=".false."/>
     <methodCall method="GSLReallocate"/>
   </interpolator>
  </stateStorable>
  <deepCopyActions class="interpolator">
   <interpolator>
    <methodCall method="GSLReallocate"/>
   </interpolator>
  </deepCopyActions>
  !!]

  type :: gslInterpWrapper
     !!{
     Wrapper class for managing GSL interpolators.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: gslInterpWrapperDestructor
  end type gslInterpWrapper
  
  type :: gslInterpAccelWrapper
     !!{
     Wrapper class for managing GSL interpolation accelerators.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: gslInterpAccelWrapperDestructor
  end type gslInterpAccelWrapper
  
  type :: interpolator
     !!{
     Type providing interpolation in 1-D arrays.
     !!}
     private
     type            (resourceManager                 )                            :: interpManager              , interpAccelManager
     type            (gslInterpWrapper                ), pointer                   :: interp_           => null()
     type            (gslInterpAccelWrapper           ), pointer                   :: interpAccel_      => null()
     type            (c_ptr                           ), allocatable               :: gsl_interp_type
     integer                                                                       :: interpolationType
     type            (enumerationExtrapolationTypeType)             , dimension(2) :: extrapolationType
     integer         (c_size_t                        )                            :: countArray
     logical                                                                       :: initialized                , interpolatable
     double precision                                  , allocatable, dimension(:) :: x                          , y 
   contains
     !![
     <methods>
       <method description="Interpolate in the tabulated function."                                                 method="interpolate"         />
       <method description="Interpolate the derivative in the tabulated function."                                  method="derivative"          />
       <method description="Locate the position in the array corresponding to the given {\normalfont \ttfamily x}." method="locate"              />
       <method description="Return factors required to perform a linear interpolation."                             method="linearFactors"       />
       <method description="Return weights required to perform a linear interpolation."                             method="linearWeights"       />
       <method description="Allocate GSL objects."                                                                  method="gslAllocate"         />
       <method description="Reallocate GSL objects."                                                                method="gslReallocate"       />
       <method description="Initialize GSL interpolator."                                                           method="gslInitialize"       />
       <method description="Assert that the data is interpolatable."                                                method="assertInterpolatable"/>
       <method description="Assign interpolator objects."                                                           method="assignment(=)"       />
     </methods>
     !!]
     procedure ::                         interpolatorAssign
     generic   :: assignment(=)        => interpolatorAssign
     procedure ::                         interpolatorInterpolate
     procedure ::                         interpolatorInterpolateNoYa
     generic   :: interpolate          => interpolatorInterpolate           , &
          &                               interpolatorInterpolateNoYa
     procedure ::                         interpolatorDerivative
     procedure ::                         interpolatorDerivativeNoYa
     generic   :: derivative           => interpolatorDerivative            , &
          &                               interpolatorDerivativeNoYa
     procedure :: linearFactors        => interpolatorLinearFactors
     procedure :: linearWeights        => interpolatorLinearWeights
     procedure :: locate               => interpolatorLocate
     procedure :: gslAllocate          => interpolatorGSLAllocate
     procedure :: gslReallocate        => interpolatorGSLReallocate
     procedure :: gslInitialize        => interpolatorGSLInitialize
     procedure :: assertInterpolatable => interpolatorAssertInterpolatable
  end type interpolator

  interface interpolator
     !!{
     Constructor for the {\normalfont \ttfamily interpolator} class.
     !!}
     module procedure interpolatorConstructor
  end interface interpolator

  !![
  <stateStorable class="interpolator2D">
   <interpolator2D>
    <methodCall method="GSLReallocate"/>
   </interpolator2D>
  </stateStorable>
  <deepCopyActions class="interpolator2D">
   <interpolator2D>
    <methodCall method="GSLReallocate"/>
   </interpolator2D>
  </deepCopyActions>
  !!]

  type :: interpolator2D
     !!{
     Type providing interpolation in 2-D arrays.
     !!}
     private
     type            (resourceManager      )                              :: interpXManager           , interpAccelXManager          , &
          &                                                                  interpYManager           , interpAccelYManager
     type            (gslInterpWrapper     ), pointer                     :: interpX         => null(), interpY             => null()
     type            (gslInterpAccelWrapper), pointer                     :: interpAccelX    => null(), interpAccelY        => null()
     type            (c_ptr                ), allocatable                 :: gsl_interp_type
     integer         (c_size_t             )                              :: countArrayX              , countArrayY
     double precision                       , allocatable, dimension(:  ) :: x                        , y 
     double precision                       , allocatable, dimension(:,:) :: z
   contains
     !![
     <methods>
       <method description="Interpolate in the tabulated function." method="interpolate"  />
       <method description="Allocate GSL objects."                  method="gslAllocate"  />
       <method description="Reallocate GSL objects."                method="gslReallocate"/>
       <method description="Assign 2D interpolator objects."        method="assignment(=)"/>
     </methods>
     !!]
     procedure ::                  interpolator2DAssign
     generic   :: assignment(=) => interpolator2DAssign
     procedure :: interpolate   => interpolator2DInterpolate
     procedure :: gslAllocate   => interpolator2DGSLAllocate
     procedure :: gslReallocate => interpolator2DGSLReallocate
  end type interpolator2D

  interface interpolator2D
     !!{
     Constructor for the {\normalfont \ttfamily interpolator2D} class.
     !!}
     module procedure interpolator2DConstructor
  end interface interpolator2D

contains

  function interpolatorConstructor(x,y,interpolationType,extrapolationType) result(self)
    !!{
    Constructor for {\normalfont \ttfamily interpolator} objects.
    !!}
    use :: Error       , only : Error_Report
    use :: Table_Labels, only : extrapolationTypeAbort
    implicit none
    type            (interpolator                    )                                         :: self
    double precision                                  , intent(in   ), dimension(: )           :: x
    double precision                                  , intent(in   ), dimension(: ), optional :: y
    integer                                           , intent(in   )               , optional :: interpolationType
    type            (enumerationExtrapolationTypeType), intent(in   ), dimension(..), optional :: extrapolationType
    !![
    <optionalArgument name="interpolationType" defaultsTo="gsl_interp_linear"     />
    !!]
    
    ! Validate data.
    if    (size(x) <       1 ) call Error_Report('"x" must contain at least 1 datapoint'  //{introspection:location})
    if (present(y)) then
       if (size(x) ==      1 ) call Error_Report('"x" must contain at least 2 datapoints' //{introspection:location})
       if (size(x) /= size(y)) call Error_Report('"x" and "y" arrays must be of same size'//{introspection:location})
    end if
    ! Store extrapolation.
    if (present(extrapolationType)) then
       select rank (extrapolationType)
       rank (0)
          self%extrapolationType=extrapolationType
       rank (1)
          self%extrapolationType=extrapolationType
       rank default
          call Error_Report('"extrapolationType must be a scalar or rank-1"'//{introspection:location})
       end select
    else
       self%extrapolationType=extrapolationTypeAbort
    end if
    ! Get the interpolator type.
    self%interpolationType=interpolationType_
    self%gsl_interp_type=gsl_interp_type_get(interpolationType_)
    ! Store data.
    self%countArray=size(x)
    allocate   (self%x(self%countArray))
    self%x         =     x
    if (present(y)) then
       allocate(self%y(self%countArray))
       self%y=y
    end if
    ! Determine if the data is interpolatable.
    self%interpolatable=self%countArray > 1
    ! Allocate GSL interpolation objects.
    call self%GSLAllocate()
    return
  end function interpolatorConstructor

  subroutine interpolatorGSLAllocate(self)
    !!{
    Allocate GSL objects.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class(interpolator), intent(inout) :: self
    class(*           ), pointer       :: dummyPointer_

    if (.not.associated(self%interp_     ) .and. self%interpolatable) then
       allocate(self%interp_     )
       self%interp_     %gsl=gsl_interp_alloc      (self%gsl_interp_type,self%countArray)
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_           => self%interp_
       self%interpManager      =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
    end if
    if (.not.associated(self%interpAccel_)                          ) then
       allocate(self%interpAccel_)
       self%interpAccel_%gsl=gsl_interp_accel_alloc(                                    )
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_           => self%interpAccel_
       self%interpAccelManager =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
    end if
    ! Initialize.
    self%initialized=.false.
    if (allocated(self%y)) then
       ! A y-array was provided, so we can initialize normally.
       call self%GSLInitialize(self%y)
    else if (self%interpolationType /= gsl_interp_linear) then
       ! For linear interpolation only we can handle having no y-array (since the GSL linear interpolator makes no use of this
       ! during initialization). In this case we defer initialization until needed.
       call Error_Report('a y-array must be provided (except for linear interpolators)'//{introspection:location})
    end if
    return
  end subroutine interpolatorGSLAllocate

  subroutine interpolatorGSLReallocate(self)
    !!{
    Reallocate GSL objects.
    !!}
    implicit none
    class(interpolator), intent(inout) :: self
    
    if (associated(self%interp_     )) then
       call self%interpManager     %release()
       nullify(self%interp_     )
    end if
    if (associated(self%interpAccel_)) then
       call self%interpAccelManager%release()
       nullify(self%interpAccel_)
    end if
    ! If x() is not allocated, then this interpolator was never initialized - so do not attempt to allocate any GSL objects.
    if (.not.allocated(self%x)) return
    ! Re-establish the interpolation type if necessary.
    if (.not.allocated(self%gsl_interp_type)) then
       allocate(self%gsl_interp_type)
       self%gsl_interp_type=gsl_interp_type_get(self%interpolationType)
    end if
    call self%GSLAllocate()
    return
  end subroutine interpolatorGSLReallocate

  subroutine gslInterpWrapperDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily gslInterpWrapper} object.
    !!}
    implicit none
    type(gslInterpWrapper), intent(inout) :: self

    call gsl_interp_free(self%gsl)
    return
  end subroutine gslInterpWrapperDestructor

  subroutine gslInterpAccelWrapperDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily gslInterpWrapper} object.
    !!}
    implicit none
    type(gslInterpAccelWrapper), intent(inout) :: self

    call gsl_interp_accel_free(self%gsl)
    return
  end subroutine gslInterpAccelWrapperDestructor

  subroutine interpolatorAssign(self,from)
    !!{
    Perform assignment of interpolators.
    !!}
    implicit none
    class(interpolator), intent(inout) :: self
    class(interpolator), intent(in   ) :: from
    
    self%interpManager      =  from%interpManager
    self%interpAccelManager =  from%interpAccelManager
    self%interp_            => from%interp_
    self%interpAccel_       => from%interpAccel_
    self%interpolationType  =  from%interpolationType
    self%extrapolationType  =  from%extrapolationType
    self%countArray         =  from%countArray
    self%initialized        =  from%initialized
    self%interpolatable     =  from%interpolatable
    if (allocated(self%gsl_interp_type)) deallocate(self%gsl_interp_type                            )
    if (allocated(from%gsl_interp_type))   allocate(self%gsl_interp_type,source=from%gsl_interp_type)
    if (allocated(self%x              )) deallocate(self%x                                          )
    if (allocated(from%x              ))   allocate(self%x              ,source=from%x              )
    if (allocated(self%y              )) deallocate(self%y                                          )
    if (allocated(from%y              ))   allocate(self%y              ,source=from%y              )    
    return
  end subroutine interpolatorAssign
  
  subroutine interpolatorGSLInitialize(self,ya)
    !!{
    Initialize GSL interpolator.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class           (interpolator), intent(inout)               :: self
    double precision              , intent(in   ), dimension(:) :: ya
    integer         (c_int       )                              :: status

    status=gsl_interp_init(self%interp_%gsl,self%x,ya,self%countArray)
    if (status /= GSL_Success) call Error_Report('failed to initialize interpolator'//{introspection:location})
    self%initialized=.true.
    return
  end subroutine interpolatorGSLInitialize

  subroutine interpolatorAssertInterpolatable(self)
    !!{
    Assert that the data is interpolatable.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(interpolator), intent(inout) :: self

    if (.not.self%interpolatable) call Error_Report('data is not interpolatable'//{introspection:location})
    return
  end subroutine interpolatorAssertInterpolatable
  
  subroutine interpolatorLinearFactors(self,x,i,h)
    !!{
    Return interpolating factors for linear interpolation in the array {\normalfont \ttfamily xArray()} given {\normalfont \ttfamily x}.
    !!}
    use :: Error       , only : Error_Report
    use :: Table_Labels, only : extrapolationTypeAbort, extrapolationTypeExtrapolate, extrapolationTypeFix, extrapolationTypeZero
    implicit none
    class           (interpolator)                , intent(inout) :: self
    double precision                              , intent(in   ) :: x
    integer         (c_size_t    )                , intent(  out) :: i
    double precision              , dimension(0:1), intent(  out) :: h
    double precision              , parameter                     :: rangeTolerance=1.0d-6
    double precision                                              :: x_

    call self%assertInterpolatable()
    ! Handle extrapolation types.
    if (x > self%x(self%countArray)) then
       ! Extrapolate to high values.
       select case (self%extrapolationType(2)%ID)
       case (extrapolationTypeExtrapolate%ID)
          x_=     x
       case (extrapolationTypeFix        %ID)
          x_=self%x(self%countArray)
       case (extrapolationTypeZero       %ID)
          i =self%countArray-1
          h =0.0d0
          return
       case (extrapolationTypeAbort      %ID)
          if (x > self%x(self%countArray)*(1.0d0+rangeTolerance)) then
             i =0
             h =0.0d0
             call Error_Report('extrapolation is not allowed'//{introspection:location})
          else
             x_=self%x(self%countArray)
          end if
       case default
          i =0
          h =0.0d0
          call Error_Report('unknown extrapolation type'  //{introspection:location})
       end select
    else if (x < self%x(              1)) then
       ! Extrapolate to low values.
       select case (self%extrapolationType(1)%ID)
       case (extrapolationTypeExtrapolate%ID)
          x_=     x
       case (extrapolationTypeFix        %ID)
          x_=self%x(              1)
       case (extrapolationTypeZero       %ID)
          i =1
          h =0.0d0
          return
       case (extrapolationTypeAbort      %ID)
          if (x < self%x(              1)*(1.0d0-rangeTolerance)) then
             i =0
             h =0.0d0
             call Error_Report('extrapolation is not allowed'//{introspection:location})
          else
             x_=self%x(              1)
          end if
       case default
          i =0
          h =0.0d0
          call Error_Report('unknown extrapolation type'  //{introspection:location})
       end select
    else
       ! No extrapolation needed.
       x_=x
    end if
    i=self%locate(x_)
    call self%linearWeights(x_,i,h)
    return
  end subroutine interpolatorLinearFactors

  subroutine interpolatorLinearWeights(self,x,i,h)
    !!{
    Return interpolation weights factors for linear interpolation in the array {\normalfont \ttfamily xArray()} given {\normalfont \ttfamily x} and {\normalfont \ttfamily i}.
    !!}
    implicit none
    class           (interpolator)                , intent(inout) :: self
    double precision                              , intent(in   ) :: x
    integer         (c_size_t    )                , intent(in   ) :: i
    double precision              , dimension(0:1), intent(  out) :: h

    call self%assertInterpolatable()
    h(0)   =+(self%x(i+1)-     x   ) &
         &  /(self%x(i+1)-self%x(i))
    h(1)   =1.0d0-h(0)
    return
  end subroutine interpolatorLinearWeights

  double precision function interpolatorInterpolateNoYa(self,x)
    !!{
    Interpolate a function to {\normalfont \ttfamily x}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (interpolator  ), intent(inout) :: self
    double precision                , intent(in   ) :: x

    if (allocated(self%y)) then
       interpolatorInterpolateNoYa=self%interpolate(x,self%y)
    else
       interpolatorInterpolateNoYa=0.0d0
       call Error_Report("no y() array specified"//{introspection:location})       
    end if
    return
  end function interpolatorInterpolateNoYa
  
  double precision function interpolatorInterpolate(self,x,ya)
    !!{
    Interpolate a function to {\normalfont \ttfamily x}.
    !!}
    use :: Error             , only : Error_Report
    use :: Interface_GSL     , only : GSL_Success           , GSL_EDom
    use :: ISO_Varying_String, only : assignment(=)         , operator(//)                , varying_string
    use :: Table_Labels      , only : extrapolationTypeAbort, extrapolationTypeExtrapolate, extrapolationTypeFix, extrapolationTypeZero
    implicit none
    class           (interpolator                    ), intent(inout)               :: self
    double precision                                  , intent(in   )               :: x
    double precision                                  , intent(in   ), dimension(:) :: ya
    double precision                                  , parameter                   :: rangeTolerance   =1.0d-6
    type            (enumerationExtrapolationTypeType)                              :: extrapolationType
    integer         (c_size_t                        )                              :: basePoint
    type            (varying_string                  ), save                        :: message
    !$omp threadprivate(message)
    integer         (c_int                           )                              :: statusGSL
    double precision                                                                :: gradient                , x_
    character       (len=12                          )                              :: labelX                  , labelX_      , &
         &                                                                             labelXMaximum           , labelXMinimum

    call self%assertInterpolatable()
    if (.not.self%initialized) call self%GSLInitialize(ya)
    ! If extrapolation is allowed, check if this is necessary.
    if (x > self%x(self%countArray)) then
       extrapolationType=self%extrapolationType(2)
    else
       extrapolationType=self%extrapolationType(1)
    end if
    if     (                                                   &
         &   extrapolationType == extrapolationTypeExtrapolate &
         &  .and.                                              &
         &   (                                                 &
         &     x < self%x(1              )                     &
         &    .or.                                             &
         &     x > self%x(self%countArray)                     &
         &   )                                                 &
         & ) then
       if (x < self%x(1)) then
          basePoint=1_c_size_t
       else
          basePoint=self%countArray
       end if
       statusGSL=gsl_interp_eval_deriv_e(self%interp_%gsl,self%x,ya,self%x(basePoint),self%interpAccel_%gsl,gradient)
       if (statusGSL /= 0) then
          select case (statusGSL)
          case (GSL_EDom)
             message='requested point is outside of allowed range'
          case default
             message='interpolation failed for unknown reason'
          end select
          call Error_Report(message//{introspection:location})
       end if
       interpolatorInterpolate=ya(basePoint)+gradient*(x-self%x(basePoint))
    else
       ! Allow for rounding errors.
       x_=x
       select case (extrapolationType%ID)
       case (extrapolationTypeFix%ID)
          if (x < self%x(1              )                                                                              ) x_=self%x(1              )
          if (x > self%x(self%countArray)                                                                              ) x_=self%x(self%countArray)
       case default
          if (x < self%x(1              ) .and. x > self%x(1              )-rangeTolerance*abs(self%x(1              ))) x_=self%x(1              )
          if (x > self%x(self%countArray) .and. x < self%x(self%countArray)+rangeTolerance*abs(self%x(self%countArray))) x_=self%x(self%countArray)
       end select
       ! Do the interpolation.
       statusGSL=gsl_interp_eval_e(self%interp_%gsl,self%x,ya,x_,self%interpAccel_%gsl,interpolatorInterpolate)
       if (statusGSL /= GSL_Success) then
          select case (statusGSL)
          case (GSL_EDom)
             if (extrapolationType == extrapolationTypeZero) then
                ! Return zero outside of the tabulated range.
                interpolatorInterpolate=0.0d0
                return
             end if
             write (labelX       ,'(e12.6)')      x
             write (labelX_      ,'(e12.6)')      x_
             write (labelXMinimum,'(e12.6)') self%x (1              )
             write (labelXMaximum,'(e12.6)') self%x (self%countArray)
             message='requested point ['//labelX//' ⟶ '//labelX_//'] is outside of allowed range ['//labelXMinimum//':'//labelXMaximum//']'
          case default
             message='interpolation failed for unknown reason'
          end select
          call Error_Report(message//{introspection:location})
       end if
    end if
    return
  end function interpolatorInterpolate
  
  double precision function interpolatorDerivativeNoYa(self,x)
    !!{
    Interpolate the derivative of the function to {\normalfont \ttfamily x}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (interpolator  ), intent(inout) :: self
    double precision                , intent(in   ) :: x

    if (allocated(self%y)) then
       interpolatorDerivativeNoYa=self%derivative(x,self%y)
    else
       interpolatorDerivativeNoYa=0.0d0
       call Error_Report("no y() array specified"//{introspection:location})       
    end if
    return
  end function interpolatorDerivativeNoYa
  
  double precision function interpolatorDerivative(self,x,ya)
    !!{
    Interpolate the derivative of the function to {\normalfont \ttfamily x}.
    !!}
    use            :: Error             , only : Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: Interface_GSL     , only : GSL_EDom
    use            :: ISO_Varying_String, only : assignment(=)          , operator(//)                , varying_string
    use            :: Table_Labels      , only : extrapolationTypeAbort , extrapolationTypeExtrapolate, extrapolationTypeFix, extrapolationTypeZero
    implicit none
    class           (interpolator                    ), intent(inout)               :: self
    double precision                                  , intent(in   )               :: x
    double precision                                  , intent(in   ), dimension(:) :: ya
    type            (varying_string                  ), save                        :: message
    !$omp threadprivate(message)
    type            (enumerationExtrapolationTypeType)                              :: extrapolationType
    integer         (c_int                           )                              :: statusGSL
    double precision                                                                :: x_

    call self%assertInterpolatable()
    if (.not.self%initialized) call self%GSLInitialize(ya)
    ! If extrapolation is allowed, check if this is necessary.
    if (x > self%x(self%countArray)) then
       extrapolationType=self%extrapolationType(2)
    else
       extrapolationType=self%extrapolationType(1)
    end if
    x_=x
    select case (extrapolationType%ID)
    case (extrapolationTypeExtrapolate%ID,extrapolationTypeFix%ID)
       if (x < self%x(1              )) x_=self%x(1              )
       if (x > self%x(self%countArray)) x_=self%x(self%countArray)
    end select
    ! Do the interpolation.
    statusGSL=gsl_interp_eval_deriv_e(self%interp_%gsl,self%x,ya,x_,self%interpAccel_%gsl,interpolatorDerivative)
    if (statusGSL /= 0) then
       select case (statusGSL)
       case (GSL_EDom)
          if (extrapolationType == extrapolationTypeZero) then
             ! Return zero outside of the tabulated range.
             interpolatorDerivative=0.0d0
             return
          end if
          message='requested point is outside of allowed range'
       case default
          message='interpolation failed for unknown reason'
       end select
       call Error_Report(message//{introspection:location})
    end if
    return
  end function interpolatorDerivative
  
  function interpolatorLocate(self,x,closest) result(i)
    !!{
    Locate the element in the table for interpolation of {\normalfont \ttfamily x}.
    !!}
    implicit none
    integer         (c_size_t    )                          :: i
    class           (interpolator), intent(inout)           :: self
    double precision              , intent(in   )           :: x
    logical                       , intent(in   ), optional :: closest
    !![
    <optionalArgument name="closest" defaultsTo=".false."/>
    !!]
    
    ! If array has just one point, always return it.
    if (self%countArray == 1_c_size_t) then
       i=1_c_size_t
       return
    end if
    ! Do the interpolation.
    i=gsl_interp_accel_find(self%interpAccel_%gsl,self%x,self%countArray,x)+1_c_size_t
    i=max(min(i,self%countArray-1_c_size_t),1_c_size_t)
    ! If we want the closest point, find it.
    if (closest_ .and. abs(x-self%x(i+1)) < abs(x-self%x(i))) i=i+1_c_size_t
    return
  end function interpolatorLocate

  function interpolator2DConstructor(x,y,z) result(self)
    !!{
    Constructor for {\normalfont \ttfamily interpolator2D} objects.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (interpolator2D)                                :: self
    double precision                , intent(in   ), dimension(:  ) :: x   , y
    double precision                , intent(in   ), dimension(:,:) :: z
    
    ! Validate data.
    if     (                                                                                         &
         &   size(x      ) <       2                                                                 &
         & ) call Error_Report('"x" must contain at least 2 datapoints'  //{introspection:location})
    if     (                                                                                         &
         &   size(y      ) <       2                                                                 &
         & ) call Error_Report('"y" must contain at least 2 datapoints'  //{introspection:location})
    if     (                                                                                         &
         &   size(z,dim=1) /= size(x)                                                                &
         &  .or.                                                                                     &
         &   size(z,dim=2) /= size(y)                                                                &
         & ) call Error_Report('"z" must be of same shape as "x" and "y"'//{introspection:location})
    ! Get the interpolator type.
    self%gsl_interp_type=gsl_interp_type_get(gsl_interp_linear)
    ! Store data.
    self%countArrayX=size(x)
    self%countArrayY=size(y)
    allocate(self%x(self%countArrayX                 ))
    allocate(self%y(                 self%countArrayY))
    allocate(self%z(self%countArrayX,self%countArrayY))
    self%x=x
    self%y=y
    self%z=z
    ! Allocate GSL interpolation objects.
    call self%GSLAllocate()
    return
  end function interpolator2DConstructor

  subroutine interpolator2DAssign(self,from)
    !!{
    Perform assignment of 2D interpolators.
    !!}
    implicit none
    class(interpolator2D), intent(inout) :: self
    class(interpolator2D), intent(in   ) :: from
    
    self%interpXManager      =  from%interpXManager
    self%interpYManager      =  from%interpYManager
    self%interpAccelXManager =  from%interpAccelXManager
    self%interpAccelYManager =  from%interpAccelYManager
    self%interpX             => from%interpX
    self%interpY             => from%interpY
    self%interpAccelX        => from%interpAccelX
    self%interpAccelY        => from%interpAccelY
    self%countArrayX         =  from%countArrayX
    self%countArrayY         =  from%countArrayY
    if (allocated(self%gsl_interp_type)) deallocate(self%gsl_interp_type                            )
    if (allocated(from%gsl_interp_type))   allocate(self%gsl_interp_type,source=from%gsl_interp_type)
    if (allocated(self%x              )) deallocate(self%x                                          )
    if (allocated(from%x              ))   allocate(self%x              ,source=from%x              )
    if (allocated(self%y              )) deallocate(self%y                                          )
    if (allocated(from%y              ))   allocate(self%y              ,source=from%y              )    
    if (allocated(self%z              )) deallocate(self%z                                          )
    if (allocated(from%z              ))   allocate(self%z              ,source=from%z              )    
    return
  end subroutine interpolator2DAssign
  
  subroutine interpolator2DGSLAllocate(self)
    !!{
    Allocate GSL objects.
    !!}
    implicit none
    class(interpolator2D), intent(inout) :: self
    class(*             ), pointer       :: dummyPointer_

    if (.not.associated(self%interpX     )) then
       allocate(self%interpX     )
       self%interpX     %gsl=gsl_interp_alloc           (self%gsl_interp_type,self%countArrayX)
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_            => self%interpX
       self%interpXManager      =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
    end if
    if (.not.associated(self%interpY     )) then
       allocate(self%interpY     )
       self%interpY     %gsl=gsl_interp_alloc           (self%gsl_interp_type,self%countArrayY)
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_            => self%interpY
       self%interpYManager      =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
    end if
    if (.not.associated(self%interpAccelX)) then
       allocate(self%interpAccelX)
       self%interpAccelX%gsl=gsl_interp_accel_alloc(                                     )
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_            => self%interpAccelX
       self%interpAccelXManager =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
    end if
    if (.not.associated(self%interpAccelY)) then
       allocate(self%interpAccelY)
       self%interpAccelY%gsl=gsl_interp_accel_alloc(                                     )
       !![
       <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
	 <description>ICE when passing a derived type component to a class(*) function argument.</description>
       !!]
       dummyPointer_            => self%interpAccelY
       self%interpAccelYManager =  resourceManager(dummyPointer_)
       !![
       </workaround>
       !!]
    end if   
    return
  end subroutine interpolator2DGSLAllocate

  subroutine interpolator2DGSLReallocate(self)
    !!{
    Reallocate GSL objects.
    !!}
    implicit none
    class  (interpolator2D), intent(inout) :: self
    
    if (associated(self%interpX      )) then
       call self%interpXManager     %release()
       nullify(self%interpX     )
    end if
    if (associated(self%interpY      )) then
       call self%interpYManager     %release()
       nullify(self%interpY     )
    end if
    if (associated(self%interpAccelX)) then
       call self%interpAccelXManager%release()
       nullify(self%interpAccelX)
    end if
    if (associated(self%interpAccelY)) then
       call self%interpAccelYManager%release()
       nullify(self%interpAccelY)
    end if
    ! Re-establish the interpolation type if necessary.
    if (.not.allocated(self%gsl_interp_type)) then
       allocate(self%gsl_interp_type)
       self%gsl_interp_type=gsl_interp_type_get(gsl_interp_linear)
    end if
    call self%GSLAllocate()
    return
  end subroutine interpolator2DGSLReallocate

  double precision function interpolator2DInterpolate(self,x,y)
    !!{
    Interpolate a function to {\normalfont \ttfamily (x,y)}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (interpolator2D), intent(inout)  :: self
    double precision                , intent(in   )  ::  x  ,  y
    integer         (c_size_t      )                 :: ix  , iy, &
         &                                              jx  , jy
    double precision                , dimension(0:1) :: hx  , hy
    ! Check for out of range values.
    if     (                                                             &
         &   x < self%x(               1)                                &
         &  .or.                                                         &
         &   x > self%x(self%countArrayX)                                &
         &  .or.                                                         &
         &   y < self%y(               1)                                &
         &  .or.                                                         &
         &   y > self%y(self%countArrayY)                                &
         & ) call Error_Report('out of range'//{introspection:location})
    ! Find the interpolating factors.
    ix=gsl_interp_accel_find(self%interpAccelX%gsl,self%x,self%countArrayX,x)+1_c_size_t
    iy=gsl_interp_accel_find(self%interpAccelY%gsl,self%y,self%countArrayY,y)+1_c_size_t
    hx(0)  =+(self%x(ix+1)-     x   )   &
         &   /(self%x(ix+1)-self%x(ix))
    hy(0)  =+(self%y(iy+1)-     y   )   &
         &  /(self%y(iy+1)-self%y(iy))
    hx(1)  =+1.0d0-hx(0)
    hy(1)  =+1.0d0-hy(0)
    ! Do the interpolation.
    interpolator2DInterpolate=0.0d0
    do jx=0,1
       do jy=0,1
          interpolator2DInterpolate=+interpolator2DInterpolate &
               &                    +self%z(ix+jx,iy+jy)       &
               &                    *hx    (   jx      )       &
               &                    *hy    (         jy)
       end do
    end do
    return
  end function interpolator2DInterpolate
  
end module Numerical_Interpolation
