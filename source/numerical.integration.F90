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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
