!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which defines the radiation structure data type, used to describe radiation fields. (Currently only includes
!% CMB radiation temperature.)

module Radiation_Structure
  !% Defines the radiation structure data type, used to describe radiation fields. (Currently only includes
  !% CMB radiation temperature.)
  implicit none
  private
  public :: radiationStructure, Radiation_Temperature, Radiation_Flux

  ! Include auto-generated labels for different radiation types.
  !# <include directive="radiationLabel" type="label" prefix="radiationType">
  include 'objects.radiation.labels.inc'
  !# </include>

  type radiationData
     !% A structure used to store data for components of radiation objects.
     double precision, allocatable, dimension(:) :: properties
  end type radiationData

  type radiationStructure
     !% The radiation structure data type, used to describe radiation fields.
     private
     double precision                               :: timeValue
     integer,             allocatable, dimension(:) :: radiationType
     type(radiationData), allocatable, dimension(:) :: components
   contains
     !@ <objectMethods>
     !@   <object>radiationStructure</object>
     !@   <objectMethod>
     !@     <method>define</method>
     !@     <description>Define the radiation components active in a given radiation object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isDefined</method>
     !@     <description>Return true if the radiation component is defined, false otherwise.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>set</method>
     !@     <description>Set the radiation components in the radiation object.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: isDefined      => Radiation_Is_Defined
     procedure :: define         => Radiation_Define
     procedure :: set            => Radiation_Set
     !@ <objectMethods>
     !@   <object>radiationStructure</object>
     !@   <objectMethod>
     !@     <method>temperature</method>
     !@     <description>Return the temperature (in units of Kelvin) of the given radiation structure.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>flux</method>
     !@     <description>Return the flux (in units of ergs cm$^2$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$) of the given radiation structure.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>integrateOverCrossSection</method>
     !@     <description>
     !@       Integrates the flux (in units of ergs cm$^2$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$) of the given radiation structure between
     !@       the wavelengths given in {\tt wavelengthRange} over a cross section specified by the function {\tt crossSectionFunction}.
     !@     </description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>time</method>
     !@     <description>The cosmic time at which this radiation object was set.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: temperature               => Radiation_Temperature
     procedure :: flux                      => Radiation_Flux
     procedure :: integrateOverCrossSection => Radiation_Integrate_Over_Cross_Section
     procedure :: time                      => Radiation_Time
  end type radiationStructure

  ! Module global variables for use in integrand routines.
  type(radiationStructure)             :: radiationGlobal
  procedure(double precision), pointer :: crossSectionFunctionGlobal
  !$omp threadprivate(radiationGlobal,crossSectionFunctionGlobal)

contains

  logical function Radiation_Is_Defined(radiation)
    !% Return true if the radiation object has been defined, false otherwise.
    implicit none
    class(radiationStructure), intent(in) :: radiation

    Radiation_Is_Defined=allocated(radiation%radiationType)
    return
  end function Radiation_Is_Defined

  subroutine Radiation_Define(radiation,radiationTypes)
    !% Define which radiation fields are active in this {\tt radiation} object.
    use Memory_Management
    implicit none
    class(radiationStructure), intent(inout)              :: radiation
    integer,                   intent(in),   dimension(:) :: radiationTypes

    ! Allocate the array of radiation types to the correct size.
    if (allocated(radiation%radiationType)) call Dealloc_Array(radiation%radiationType)
    if (allocated(radiation%components   )) then
       deallocate(radiation%components)
       call Memory_Usage_Record(sizeof(radiation%components),addRemove=-1)
    end if
    call Alloc_Array(radiation%radiationType,shape(radiationTypes))
    allocate(radiation%components(size(radiationTypes)))
    call Memory_Usage_Record(sizeof(radiation%components))

    ! Populate the array of radiation types.
    radiation%radiationType=radiationTypes
    return
  end subroutine Radiation_Define

  subroutine Radiation_Set(radiation,thisNode)
    !% Set the {\tt radiation} field as specified.
    !# <include directive="radiationSet" type="moduleUse">
    include 'objects.radiation.set.modules.inc'
    !# </include>
    use Tree_Nodes
    implicit none
    class(radiationStructure), intent(inout)          :: radiation
    type(treeNode),            intent(inout), pointer :: thisNode         
    integer                                           :: iComponent

    ! For an unallocated radiation object, return immediately.
    if (.not.allocated(radiation%radiationType)) return

    ! Set the time.
    radiation%timeValue=Tree_Node_Time(thisNode)

    ! Loop over all radiation components.
    do iComponent=1,size(radiation%radiationType)
       ! Call the appropriate routine to set the component.
       !# <include directive="radiationSet" type="code" action="subroutine">
       !#  <subroutineArgs>radiation%radiationType(iComponent)==radiationType#label,thisNode,radiation%components(iComponent)%properties</subroutineArgs>
       include 'objects.radiation.set.inc'
       !# </include>
    end do
    return
  end subroutine Radiation_Set

  double precision function Radiation_Time(radiation)
    !% Return the time of the {\tt radiation} object.
    implicit none
    class(radiationStructure), intent(in) :: radiation

    Radiation_Time=radiation%timeValue
    return
  end function Radiation_Time

  double precision function Radiation_Temperature(radiation,radiationType)
    !% Return the temperature of the {\tt radiation} object.
    !# <include directive="radiationTemperature" type="moduleUse">
    include 'objects.radiation.temperature.modules.inc'
    !# </include>
    implicit none
    class(radiationStructure), intent(in)                         :: radiation
    integer,                   intent(in), optional, dimension(:) :: radiationType
    integer                                                       :: iComponent

    ! Loop over all radiation components.
    Radiation_Temperature=0.0d0
    do iComponent=1,size(radiation%radiationType)
       ! Call the appropriate routine to get the temperature.
       !# <include directive="radiationTemperature" type="code" action="subroutine">
       !#  <subroutineArgs>radiation%radiationType(iComponent),radiationType#label,radiation%components(iComponent)%properties,Radiation_Temperature,radiationType</subroutineArgs>
       include 'objects.radiation.temperature.inc'
       !# </include>
    end do
    return
  end function Radiation_Temperature

  double precision function Radiation_Flux(radiation,wavelength,radiationType)
    !% Return the flux of the {\tt radiation} object in units of ergs cm$^2$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$ at the specified {\tt
    !% wavelength} (in \AA).
    !# <include directive="radiationFlux" type="moduleUse">
    include 'objects.radiation.flux.modules.inc'
    !# </include>
    implicit none
    class(radiationStructure), intent(in)                        :: radiation
    double precision,         intent(in)                         :: wavelength
    integer,                  intent(in), optional, dimension(:) :: radiationType
    integer                                                      :: iComponent

    ! Loop over all radiation components.
    Radiation_Flux=0.0d0
    do iComponent=1,size(radiation%radiationType)
       ! Call the appropriate routine to get the flux.
       !# <include directive="radiationFlux" type="code" action="subroutine">
       !#  <subroutineArgs>radiation%radiationType(iComponent),radiationType#label,radiation%components(iComponent)%properties,wavelength,Radiation_Flux,radiationType</subroutineArgs>
       include 'objects.radiation.flux.inc'
       !# </include>
    end do
    return
  end function Radiation_Flux

  double precision function Radiation_Integrate_Over_Cross_Section(radiation,crossSectionFunction,wavelengthRange)
    !% Integrate the photon number of the radiation field over a given cross-section function (which should return the cross
    !% section in units of cm$^2$), i.e.:
    !% \begin{equation}
    !% {4 \pi \over {\rm h}} \int_{\lambda_1}^{\lambda_2} \sigma(\lambda) j_{\nu}(\lambda) {{\rm d}\lambda \over \lambda},
    !% \end{equation}
    !% where $j_{\nu}$ is the flux of energy per unit area per unit solid angle and per unit frequency.
    use, intrinsic :: ISO_C_Binding                             
    use Numerical_Integration
    use FGSL
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    class(radiationStructure), intent(in)              :: radiation
    double precision,         intent(in), dimension(2) :: wavelengthRange
    double precision,         external                 :: crossSectionFunction
    type(c_ptr)                                        :: parameterPointer
    type(fgsl_function)                                :: integrandFunction
    type(fgsl_integration_workspace)                   :: integrationWorkspace

    ! Copy the radiation object to a module global copy for use in the integrand routine.
    select type (radiation)
    type is (radiationStructure)
       radiationGlobal=radiation
    end select

    ! Copy the procedure pointer to a module global copy for use in the integrand routine.
    crossSectionFunctionGlobal => crossSectionFunction
    
    ! Perform the integration.
    Radiation_Integrate_Over_Cross_Section=Integrate( wavelengthRange(1),wavelengthRange(2)            &
         &                                           ,Cross_Section_Integrand,parameterPointer         &
         &                                           ,integrandFunction,integrationWorkspace           &
         &                                           ,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3 &
         &                                           ,integrationRule =FGSL_Integ_Gauss15&
         &                                          )
    call Integrate_Done(integrandFunction,integrationWorkspace)

    ! Scale result by multiplicative prefactors to give answer in units of inverse seconds.
    Radiation_Integrate_Over_Cross_Section=Radiation_Integrate_Over_Cross_Section*4.0d0*Pi*ergs/plancksConstant

    return
  end function Radiation_Integrate_Over_Cross_Section

  function Cross_Section_Integrand(wavelength,parameterPointer) bind(c)
    !% Integrand function use in integrating a radiation field over a cross section function.
    use, intrinsic :: ISO_C_Binding                             
    real(c_double)        :: Cross_Section_Integrand
    real(c_double), value :: wavelength
    type(c_ptr),    value :: parameterPointer

    if (wavelength > 0.0d0) then
       Cross_Section_Integrand=crossSectionFunctionGlobal(wavelength)*Radiation_Flux(radiationGlobal,wavelength)/wavelength
    else
       Cross_Section_Integrand=0.0d0
    end if
    return
  end function Cross_Section_Integrand

end module Radiation_Structure
