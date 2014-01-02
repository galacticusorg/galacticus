!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which performs numerical integration.

module Numerical_Integration
  !% Implements numerical integration.
  use FGSL
  implicit none
  private
  public :: Integrate, Integrate_Done

  ! Module scope error status.
  integer :: errorStatusGlobal

contains

  recursive double precision function Integrate(lowerLimit,upperLimit,integrand,parameterPointer,integrandFunction&
       &,integrationWorkspace,maxIntervals,toleranceAbsolute,toleranceRelative,hasSingularities,integrationRule,reset,errorStatus)
    !% Integrates the supplied {\tt integrand} function.
    use Galacticus_Error
    use, intrinsic :: ISO_C_Binding
    implicit none
    double precision                            , external                           :: integrand
    double precision                                       , intent(in   )           :: lowerLimit                      , upperLimit
    type            (c_ptr                     )           , intent(in   )           :: parameterPointer
    type            (fgsl_function             )           , intent(inout)           :: integrandFunction
    type            (fgsl_integration_workspace)           , intent(inout)           :: integrationWorkspace
    type            (fgsl_error_handler_t      )                                     :: integrationErrorHandler         , standardGslErrorHandler
    integer                                                , intent(in   ), optional :: integrationRule                 , maxIntervals
    double precision                                       , intent(in   ), optional :: toleranceAbsolute               , toleranceRelative
    logical                                                , intent(in   ), optional :: hasSingularities
    logical                                                , intent(inout), optional :: reset
    integer                                                , intent(  out), optional :: errorStatus
    integer                                     , parameter                          :: maxIntervalsDefault     =1000
    double precision                            , parameter                          :: toleranceAbsoluteDefault=1.0d-10, toleranceRelativeDefault=1.0d-10
    integer                                                                          :: integrationRuleActual           , status
    integer         (kind=c_size_t             )                                     :: maxIntervalsActual
    double precision                                                                 :: integrationError                , integrationValue                , &
         &                                                                              toleranceAbsoluteActual         , toleranceRelativeActual
    logical                                                                          :: hasSingularitiesActual          , resetActual

    ! Set optional parameters to specified or default values.
    if (present(maxIntervals)) then
       maxIntervalsActual=maxIntervals
    else
       maxIntervalsActual=maxIntervalsDefault
    end if
    if (present(toleranceAbsolute)) then
       toleranceAbsoluteActual=toleranceAbsolute
    else
       toleranceAbsoluteActual=toleranceAbsoluteDefault
    end if
    if (present(toleranceRelative)) then
       toleranceRelativeActual=toleranceRelative
    else
       toleranceRelativeActual=toleranceRelativeDefault
    end if
    if (present(hasSingularities)) then
       hasSingularitiesActual=hasSingularities
    else
       hasSingularitiesActual=.false.
    end if
    if (present(integrationRule)) then
       integrationRuleActual=integrationRule
    else
       integrationRuleActual=FGSL_Integ_Gauss61
    end if
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    else
       resetActual=.true.
    end if

    ! Initialize the integration variables if necessary.
    if (resetActual) then
       integrationWorkspace=FGSL_Integration_Workspace_Alloc(maxIntervalsActual)
       integrandFunction   =FGSL_Function_Init(integrand,parameterPointer)
    end if

    ! Set error handler if necessary.
    if (present(errorStatus)) then
       integrationErrorHandler=FGSL_Error_Handler_Init(Integration_GSL_Error_Handler)
       standardGslErrorHandler=FGSL_Set_Error_Handler (integrationErrorHandler      )
       errorStatusGlobal=errorStatusSuccess
    end if

    ! Do the integration
    select case (hasSingularitiesActual)
    case (.false.)
       status=FGSL_Integration_QAG(integrandFunction,lowerLimit,upperLimit,toleranceAbsoluteActual,toleranceRelativeActual &
            &,maxIntervalsActual,integrationRuleActual,integrationWorkspace,integrationValue,integrationError)
    case (.true.)
       status=FGSL_Integration_QAGS(integrandFunction,lowerLimit,upperLimit,toleranceAbsoluteActual,toleranceRelativeActual &
            &,maxIntervalsActual,integrationWorkspace,integrationValue,integrationError)
    end select
    Integrate=integrationValue

    ! Reset error handler.
    if (present(errorStatus)) then
       errorStatus            =errorStatusGlobal
       standardGslErrorHandler=FGSL_Set_Error_Handler (standardGslErrorHandler)
    end if
    return
  end function Integrate

  subroutine Integrate_Done(integrandFunction,integrationWorkspace)
    !% Frees up integration objects that are no longer required.
    implicit none
    type(fgsl_function             ), intent(inout) :: integrandFunction
    type(fgsl_integration_workspace), intent(inout) :: integrationWorkspace

    call FGSL_Function_Free(integrandFunction)
    call FGSL_Integration_Workspace_Free(integrationWorkspace)
    return
  end subroutine Integrate_Done

  subroutine Integration_GSL_Error_Handler(reason,file,line,errorNumber) bind(c)
    !% Handle errors from the GSL library during integration.
    use, intrinsic :: ISO_C_Binding
    type   (c_ptr     ), value :: file       , reason
    integer(kind=c_int), value :: errorNumber, line

    errorStatusGlobal=errorNumber
    return
  end subroutine Integration_GSL_Error_Handler

end module Numerical_Integration
