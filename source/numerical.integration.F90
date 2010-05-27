!% Contains a module which performs numerical integration.

module Numerical_Integration
  !% Implements numerical integration.
  use FGSL
  private
  public :: Integrate, Integrate_Done

contains

  recursive double precision function Integrate(lowerLimit,upperLimit,integrand,parameterPointer,integrandFunction&
       &,integrationWorkspace,maxIntervals,toleranceAbsolute,toleranceRelative,hasSingularities,reset)
    !% Integrates the supplied {\tt integrand} function.
    use, intrinsic :: ISO_C_Binding                             
    implicit none
    double precision,                 external                :: integrand
    double precision,                 intent(in)              :: lowerLimit,upperLimit
    type(c_ptr),                      intent(in)              :: parameterPointer
    type(fgsl_function),              intent(inout)           :: integrandFunction
    type(fgsl_integration_workspace), intent(inout)           :: integrationWorkspace
    integer,                          intent(in),    optional :: maxIntervals
    double precision,                 intent(in),    optional :: toleranceAbsolute,toleranceRelative
    logical,                          intent(in),    optional :: hasSingularities
    logical,                          intent(inout), optional :: reset
    integer,                          parameter               :: maxIntervalsDefault=1000
    double precision,                 parameter               :: toleranceAbsoluteDefault=1.0d-10,toleranceRelativeDefault=1.0d-10
    integer                                                   :: status
    integer(c_size_t)                                         :: maxIntervalsActual
    double precision                                          :: toleranceAbsoluteActual,toleranceRelativeActual,integrationValue&
         &,integrationError
    logical                                                   :: hasSingularitiesActual,resetActual
    
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

    ! Do the integration
    select case (hasSingularitiesActual)
    case (.false.)
       status=FGSL_Integration_QAG(integrandFunction,lowerLimit,upperLimit,toleranceAbsoluteActual,toleranceRelativeActual &
            &,maxIntervalsActual,FGSL_INTEG_GAUSS61,integrationWorkspace,integrationValue,integrationError)
    case (.true.)
       status=FGSL_Integration_QAGS(integrandFunction,lowerLimit,upperLimit,toleranceAbsoluteActual,toleranceRelativeActual &
            &,maxIntervalsActual,integrationWorkspace,integrationValue,integrationError)
    end select
    Integrate=integrationValue
    return
  end function Integrate

  subroutine Integrate_Done(integrandFunction,integrationWorkspace)
    !% Frees up integration objects that are no longer required.
    implicit none
    type(fgsl_function),              intent(inout) :: integrandFunction
    type(fgsl_integration_workspace), intent(inout) :: integrationWorkspace

    call FGSL_Function_Free(integrandFunction)
    call FGSL_Integration_Workspace_Free(integrationWorkspace)
    return
  end subroutine Integrate_Done

end module Numerical_Integration
