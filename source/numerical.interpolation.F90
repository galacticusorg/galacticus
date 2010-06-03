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






!% Contains a module which acts as a simple interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library}
!% \href{http://www.gnu.org/software/gsl/manual/html_node/Interpolation.html}{interpolation routines}.

module Numerical_Interpolation
  !% A simple interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library}
  !% \href{http://www.gnu.org/software/gsl/manual/html_node/Interpolation.html}{interpolation routines}.
  use FGSL
  use, intrinsic :: ISO_C_Binding
  private
  public :: Interpolate, Interpolate_Derivative, Interpolate_Locate, Interpolate_Done, Interpolate_Linear_Generate_Factors,&
       & Interpolate_Linear_Do
  
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

  double precision function Interpolate(nPoints,xArray,yArray,interpolationObject,interpolationAccelerator,x,interpolationType&
       &,reset)
    !% Perform an interpolation of {\tt x} into {\tt xArray()} and return the corresponding value in {\tt yArray()}.
    use Galacticus_Error
    implicit none
    integer,                 intent(in)                  :: nPoints
    double precision,        intent(in),    dimension(:) :: xArray,yArray
    type(fgsl_interp),       intent(inout)               :: interpolationObject
    type(fgsl_interp_accel), intent(inout)               :: interpolationAccelerator
    double precision,        intent(in)                  :: x
    type(fgsl_interp_type),  intent(in),    optional     :: interpolationType
    logical,                 intent(inout), optional     :: reset
    type(fgsl_interp_type)                               :: interpolationTypeActual
    integer                                              :: status
    integer(c_size_t)                                    :: nPointsC
    logical                                              :: resetActual

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

    ! Do the interpolation.
    Interpolate=fgsl_interp_eval(interpolationObject,xArray,yArray,x,interpolationAccelerator)
    return
  end function Interpolate

  double precision function Interpolate_Derivative(nPoints,xArray,yArray,interpolationObject,interpolationAccelerator,x&
       &,interpolationType ,reset)
    !% Perform an interpolation of {\tt x} into {\tt xArray()} and return the corresponding first derivative of {\tt yArray()}.
    use Galacticus_Error
    implicit none
    integer,                 intent(in)                  :: nPoints
    double precision,        intent(in),    dimension(:) :: xArray,yArray
    type(fgsl_interp),       intent(inout)               :: interpolationObject
    type(fgsl_interp_accel), intent(inout)               :: interpolationAccelerator
    double precision,        intent(in)                  :: x
    type(fgsl_interp_type),  intent(in),    optional     :: interpolationType
    logical,                 intent(inout), optional     :: reset
    type(fgsl_interp_type)                               :: interpolationTypeActual
    integer                                              :: status
    integer(c_size_t)                                    :: nPointsC
    logical                                              :: resetActual

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

    ! Do the interpolation.
    Interpolate_Derivative=fgsl_interp_eval_deriv(interpolationObject,xArray,yArray,x,interpolationAccelerator)
    return
  end function Interpolate_Derivative

  integer function Interpolate_Locate(nPoints,xArray,interpolationAccelerator,x,reset)
    !% Perform an interpolation of {\tt x} into {\tt xArray()} and return the corresponding value in {\tt yArray()}.
    use Galacticus_Error
    implicit none
    integer,                 intent(in)                  :: nPoints
    double precision,        intent(in),    dimension(:) :: xArray
    type(fgsl_interp_accel), intent(inout)               :: interpolationAccelerator
    double precision,        intent(in)                  :: x
    logical,                 intent(inout), optional     :: reset
    integer(c_size_t)                                    :: nPointsC
    logical                                              :: resetActual

    ! Decide whether to reset.
    resetActual=.false.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
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
