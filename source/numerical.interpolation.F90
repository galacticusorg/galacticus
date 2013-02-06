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

!% Contains a module which acts as a simple interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library}
!% \href{http://www.gnu.org/software/gsl/manual/html_node/Interpolation.html}{interpolation routines}.

module Numerical_Interpolation
  !% A simple interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library}
  !% \href{http://www.gnu.org/software/gsl/manual/html_node/Interpolation.html}{interpolation routines}.
  use FGSL
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Interpolate, Interpolate_Derivative, Interpolate_Locate, Interpolate_Done, Interpolate_Linear_Generate_Factors,&
       & Interpolate_Linear_Do, Interpolate_Linear_Generate_Gradient_Factors
  
  ! Labels for extrapolation methods.
  integer, parameter, public :: extrapolationTypeNone  =0
  integer, parameter, public :: extrapolationTypeLinear=1
  integer, parameter, public :: extrapolationTypeFixed =2

contains
  
  double precision function Interpolate_Linear_Do(nPoints,yArray,iInterpolate,interpolationFactors)
    !% Given an array index {\tt iInterpolate} and interpolating factors {\tt interpolationFactors} for array {\tt yArray}, return
    !% a linearly interpolated value.
    implicit none
    integer,          intent(in) :: nPoints,iInterpolate
    double precision, intent(in) :: yArray(:)
    double precision, intent(in) :: interpolationFactors(2)

    Interpolate_Linear_Do=yArray(iInterpolate)*interpolationFactors(1)+yArray(iInterpolate+1)*interpolationFactors(2)
    return
  end function Interpolate_Linear_Do

  function Interpolate_Linear_Generate_Factors(nPoints,xArray,iInterpolate,x)
    !% Return interpolating factors for linear interpolation in the array {\tt xArray()} given the index in the array which
    !% brackets value {\tt x}.
    implicit none
    double precision,             dimension(0:1) :: Interpolate_Linear_Generate_Factors
    integer,          intent(in)                 :: nPoints,iinterpolate
    double precision, intent(in), dimension(:)   :: xArray
    double precision, intent(in)                 :: x
    
    Interpolate_Linear_Generate_Factors(0)=(xArray(iInterpolate+1)-x)/(xArray(iInterpolate+1)-xArray(iInterpolate))
    Interpolate_Linear_Generate_Factors(1)=1.0d0-Interpolate_Linear_Generate_Factors(0)
    return
  end function Interpolate_Linear_Generate_Factors

  function Interpolate_Linear_Generate_Gradient_Factors(nPoints,xArray,iInterpolate,x)
    !% Return interpolating factors for linear interpolation in the array {\tt xArray()} given the index in the array which
    !% brackets value {\tt x}.
    implicit none
    double precision,             dimension(0:1) :: Interpolate_Linear_Generate_Gradient_Factors
    integer,          intent(in)                 :: nPoints,iinterpolate
    double precision, intent(in), dimension(:)   :: xArray
    double precision, intent(in)                 :: x
    
    Interpolate_Linear_Generate_Gradient_Factors(0)=-1.0d0/(xArray(iInterpolate+1)-xArray(iInterpolate))
    Interpolate_Linear_Generate_Gradient_Factors(1)=-Interpolate_Linear_Generate_Gradient_Factors(0)
    return
  end function Interpolate_Linear_Generate_Gradient_Factors

  double precision function Interpolate(nPoints,xArray,yArray,interpolationObject,interpolationAccelerator,x,interpolationType&
       &,extrapolationType,reset)
    !% Perform an interpolation of {\tt x} into {\tt xArray()} and return the corresponding value in {\tt yArray()}.
    use Galacticus_Error
    use ISO_Varying_String
    implicit none
    integer,                 intent(in)                  :: nPoints
    double precision,        intent(in),    dimension(:) :: xArray,yArray
    type(fgsl_interp),       intent(inout)               :: interpolationObject
    type(fgsl_interp_accel), intent(inout)               :: interpolationAccelerator
    double precision,        intent(in)                  :: x
    type(fgsl_interp_type),  intent(in),    optional     :: interpolationType
    logical,                 intent(inout), optional     :: reset
    integer,                 intent(in),    optional     :: extrapolationType
    double precision,        parameter                   :: rangeTolerance=1.0d-6
    type(fgsl_interp_type)                               :: interpolationTypeActual
    integer                                              :: status,basePoint,extrapolationTypeActual
    integer(c_size_t)                                    :: nPointsC
    logical                                              :: resetActual
    type(varying_string)                                 :: message
    integer(fgsl_int)                                    :: errorCode
    double precision                                     :: gradient,xActual

    ! Decide whether to reset.
    resetActual=.false.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    end if
    ! Check if this interpolation needs initializing.
    if (resetActual.or..not.fgsl_well_defined(interpolationAccelerator).or..not.fgsl_well_defined(interpolationObject)) then
       ! Allocate the accelerator.
       interpolationAccelerator=fgsl_interp_accel_alloc()
       ! Determine the interpolation type.
       if (present(interpolationType)) then
          interpolationTypeActual=interpolationType
       else
          interpolationTypeActual=fgsl_interp_linear
       end if
       ! Allocate the interpolation object.
       nPointsC=nPoints
       interpolationObject=fgsl_interp_alloc(interpolationTypeActual,nPointsC)
       ! Check status.
       status=fgsl_interp_init(interpolationObject,xArray,yArray,nPointsC)
       if (status /= fgsl_success) call Galacticus_Error_Report('Interpolate','interpolation initialization failed')
    end if

    ! If extrapolation is allowed, check if this is necessary.
    if (present(extrapolationType)) then
       extrapolationTypeActual=extrapolationType
    else
       extrapolationTypeActual=extrapolationTypeNone
    end if
    if     (                                                    &
         &   extrapolationTypeActual == extrapolationTypeLinear &
         &  .and.                                               &
         &   (                                                  &
         &     x < xArray(1           )                         &
         &    .or.                                              &
         &     x > xArray(size(xArray))                         &
         &   )                                                  &
         & ) then
       if (x < xArray(1)) then
          basePoint=1
       else
          basePoint=size(xArray)
       end if
       errorCode=fgsl_interp_eval_deriv_e(interpolationObject,xArray,yArray,xArray(basePoint),interpolationAccelerator,gradient)
       if (errorCode /= 0) then
          select case (errorCode)
          case (FGSL_EDOM)
             message='requested point is outside of allowed range'
          case default
             message='interpolation failed for unknown reason'
          end select
          call Galacticus_Error_Report('Interpolate',message)
       end if
       Interpolate=yArray(basePoint)+gradient*(x-xArray(basePoint))
    else
       ! Allow for rounding errors.	
       xActual=x
       select case (extrapolationTypeActual)
       case (extrapolationTypeFixed)
          if (x < xArray(1           )                                                                         ) xActual=xArray(1           )
          if (x > xArray(size(xArray))                                                                         ) xActual=xArray(size(xArray))
       case default
          if (x < xArray(1           ) .and. x > xArray(1           )-rangeTolerance*dabs(xArray(1)           )) xActual=xArray(1           )
          if (x > xArray(size(xArray)) .and. x < xArray(size(xArray))+rangeTolerance*dabs(xArray(size(xArray)))) xActual=xArray(size(xArray))
       end select
       ! Do the interpolation.
       errorCode=fgsl_interp_eval_e(interpolationObject,xArray,yArray,xActual,interpolationAccelerator,Interpolate)
       if (errorCode /= 0) then
          select case (errorCode)
          case (FGSL_EDOM)
             message='requested point is outside of allowed range'
          case default
             message='interpolation failed for unknown reason'
          end select
          call Galacticus_Error_Report('Interpolate',message)
       end if
    end if
    return
  end function Interpolate

  double precision function Interpolate_Derivative(nPoints,xArray,yArray,interpolationObject,interpolationAccelerator,x&
       &,interpolationType ,extrapolationType,reset)
    !% Perform an interpolation of {\tt x} into {\tt xArray()} and return the corresponding first derivative of {\tt yArray()}.
    use Galacticus_Error
    use ISO_Varying_String
    implicit none
    integer,                 intent(in)                  :: nPoints
    double precision,        intent(in),    dimension(:) :: xArray,yArray
    type(fgsl_interp),       intent(inout)               :: interpolationObject
    type(fgsl_interp_accel), intent(inout)               :: interpolationAccelerator
    double precision,        intent(in)                  :: x
    type(fgsl_interp_type),  intent(in),    optional     :: interpolationType
    logical,                 intent(inout), optional     :: reset
    integer,                 intent(in),    optional     :: extrapolationType
    type(fgsl_interp_type)                               :: interpolationTypeActual
    integer                                              :: status,extrapolationTypeActual
    integer(c_size_t)                                    :: nPointsC
    logical                                              :: resetActual
    type(varying_string)                                 :: message
    integer(fgsl_int)                                    :: errorCode
    double precision                                     :: xActual

    ! Decide whether to reset.
    resetActual=.false.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    end if
    ! Check if this interpolation needs initializing.
    if (resetActual.or..not.fgsl_well_defined(interpolationAccelerator).or..not.fgsl_well_defined(interpolationObject)) then
       ! Allocate the accelerator.
       interpolationAccelerator=fgsl_interp_accel_alloc()
       ! Determine the interpolation type.
       if (present(interpolationType)) then
          interpolationTypeActual=interpolationType
       else
          interpolationTypeActual=fgsl_interp_linear
       end if
       ! Allocate the interpolation object.
       nPointsC=nPoints
       interpolationObject=fgsl_interp_alloc(interpolationTypeActual,nPointsC)
       ! Check status.
       status=fgsl_interp_init(interpolationObject,xArray,yArray,nPointsC)
       if (status /= fgsl_success) call Galacticus_Error_Report('Interpolate_Derivative','interpolation initialization failed')
    end if


    ! If extrapolation is allowed, check if this is necessary.
    xActual=x
    if (present(extrapolationType)) then
       extrapolationTypeActual=extrapolationType
    else
       extrapolationTypeActual=extrapolationTypeNone
    end if
    select case (extrapolationTypeActual)
    case (extrapolationTypeLinear,extrapolationTypeFixed)
       if (x < xArray(1           )) xActual=xArray(1           )
       if (x > xArray(size(xArray))) xActual=xArray(size(xArray))
    end select
    ! Do the interpolation.
    errorCode=fgsl_interp_eval_deriv_e(interpolationObject,xArray,yArray,xActual,interpolationAccelerator,Interpolate_Derivative)
    if (errorCode /= 0) then
       select case (errorCode)
       case (FGSL_EDOM)
          message='requested point is outside of allowed range'
       case default
          message='interpolation failed for unknown reason'
       end select
       call Galacticus_Error_Report('Interpolate_Derivative',message)
    end if
     return
   end function Interpolate_Derivative

  integer function Interpolate_Locate(nPoints,xArray,interpolationAccelerator,x,reset,closest)
    !% Perform an interpolation of {\tt x} into {\tt xArray()} and return the corresponding value in {\tt yArray()}.
    use Galacticus_Error
    implicit none
    integer,                 intent(in)                  :: nPoints
    double precision,        intent(in),    dimension(:) :: xArray
    type(fgsl_interp_accel), intent(inout)               :: interpolationAccelerator
    double precision,        intent(in)                  :: x
    logical,                 intent(inout), optional     :: reset
    logical,                 intent(in),    optional     :: closest
    integer(c_size_t)                                    :: nPointsC
    logical                                              :: resetActual,closestActual

    ! Abort on non-positive sized arrays.
    if (nPoints <= 0) call Galacticus_Error_Report('Interpolate_Locate','array has non-positive size')

    ! If array has just one point, always return it.
    if (nPoints == 1) then
       Interpolate_Locate=1
       return
    end if

    ! Decide whether to reset.
    resetActual=.false.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    end if

    ! Determine if we want the closest point.
    if (present(closest)) then
       closestActual=closest
    else
       closestActual=.false.
    end if

    ! Check if this interpolation needs initializing.
    if (resetActual.or..not.fgsl_well_defined(interpolationAccelerator)) then
       ! Allocate the accelerator.
       interpolationAccelerator=fgsl_interp_accel_alloc()
    end if

    ! Do the interpolation.
    nPointsC=nPoints
    Interpolate_Locate=fgsl_interp_accel_find(interpolationAccelerator,xArray,nPointsC,x)
    Interpolate_Locate=max(min(Interpolate_Locate,nPointsC-1),1)

    ! If we want the closest point, find it.
    if (closestActual .and. dabs(x-xArray(Interpolate_Locate+1)) < dabs(x-xArray(Interpolate_Locate))) Interpolate_Locate=Interpolate_Locate+1

    return
  end function Interpolate_Locate

  subroutine Interpolate_Done(interpolationObject,interpolationAccelerator,reset)
    !% Free interpolation objects when they are no longer required.
    implicit none
    type(fgsl_interp),       intent(inout), optional :: interpolationObject
    type(fgsl_interp_accel), intent(inout), optional :: interpolationAccelerator
    logical,                 intent(in),    optional :: reset
    logical                                          :: resetActual
    
    ! Determine reset status.
    if (present(reset)) then
       resetActual=reset
    else
       resetActual=.false.
    end if
    ! If reset is true then these objects haven't been allocated and so we don't need to do anything.
    if (.not.resetActual) then
       if (present(interpolationObject     )) call fgsl_interp_free      (interpolationObject     )
       if (present(interpolationAccelerator)) call fgsl_interp_accel_free(interpolationAccelerator)
    end if
    return
  end subroutine Interpolate_Done

end module Numerical_Interpolation
