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
     double precision                                           :: timeValue
     integer                        , allocatable, dimension(:) :: radiationType
     type            (radiationData), allocatable, dimension(:) :: components
   contains
     !@ <objectMethods>
     !@   <object>radiationStructure</object>
     !@   <objectMethod>
     !@     <method>define</method>
     !@     <description>Define the radiation components active in a given radiation object.</description>
     !@     <type>\void</type>
     !@     <arguments>\intone\ radiationTypes\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isDefined</method>
     !@     <description>Return true if the radiation component is defined, false otherwise.</description>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>set</method>
     !@     <description>Set the radiation components in the radiation object.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless *type(treeNode)\textgreater thisNode}</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>temperature</method>
     !@     <description>Return the temperature (in units of Kelvin) of the given radiation structure.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\intone\ [radiationType]\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>flux</method>
     !@     <description>Return the flux (in units of ergs cm$^2$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$) of the given radiation structure.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ wavelength\argin, \intone\ [radiationType]\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>integrateOverCrossSection</method>
     !@     <description>
     !@       Integrates the flux (in units of ergs cm$^2$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$) of the given radiation structure between
     !@       the wavelengths given in {\tt wavelengthRange} over a cross section specified by the function {\tt crossSectionFunction}.
     !@     </description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ wavelength\argin, \intone\ [radiationType]\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>time</method>
     !@     <description>The cosmic time at which this radiation object was set.</description>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: isDefined                =>Radiation_Is_Defined
     procedure :: define                   =>Radiation_Define
     procedure :: set                      =>Radiation_Set
     procedure :: temperature              =>Radiation_Temperature
     procedure :: flux                     =>Radiation_Flux
     procedure :: integrateOverCrossSection=>Radiation_Integrate_Over_Cross_Section
     procedure :: time                     =>Radiation_Time
  end type radiationStructure

  ! Module global variables for use in integrand routines.
  type     (radiationStructure)          :: radiationGlobal
  procedure(double precision  ), pointer :: crossSectionFunctionGlobal
  !$omp threadprivate(radiationGlobal,crossSectionFunctionGlobal)
contains

  logical function Radiation_Is_Defined(radiation)
    !% Return true if the radiation object has been defined, false otherwise.
    implicit none
    class(radiationStructure), intent(in   ) :: radiation

    Radiation_Is_Defined=allocated(radiation%radiationType)
    return
  end function Radiation_Is_Defined

  subroutine Radiation_Define(radiation,radiationTypes)
    !% Define which radiation fields are active in this {\tt radiation} object.
    use Memory_Management
    implicit none
    class  (radiationStructure)              , intent(inout) :: radiation
    integer                    , dimension(:), intent(in   ) :: radiationTypes

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
    use Galacticus_Nodes
    implicit none
    class  (radiationStructure), intent(inout)          :: radiation
    type   (treeNode          ), intent(inout), pointer :: thisNode
    class  (nodeComponentBasic)               , pointer :: thisBasicComponent
    integer                                             :: iComponent

    ! For an unallocated radiation object, return immediately.
    if (.not.allocated(radiation%radiationType)) return

    ! Set the time.
    thisBasicComponent => thisNode%basic()
    radiation%timeValue=thisBasicComponent%time()

    ! Loop over all radiation components.
    do iComponent=1,size(radiation%radiationType)
       ! Call the appropriate routine to set the component.
       !# <include directive="radiationSet" type="functionCall" functionType="void">
       !#  <functionArgs>radiation%radiationType(iComponent)==radiationType#label,thisNode,radiation%components(iComponent)%properties</functionArgs>
       include 'objects.radiation.set.inc'
       !# </include>
    end do
    return
  end subroutine Radiation_Set

  double precision function Radiation_Time(radiation)
    !% Return the time of the {\tt radiation} object.
    implicit none
    class(radiationStructure), intent(in   ) :: radiation

    Radiation_Time=radiation%timeValue
    return
  end function Radiation_Time

  double precision function Radiation_Temperature(radiation,radiationType)
    !% Return the temperature of the {\tt radiation} object.
    !# <include directive="radiationTemperature" type="moduleUse">
    include 'objects.radiation.temperature.modules.inc'
    !# </include>
    implicit none
    class  (radiationStructure)              , intent(in   )           :: radiation
    integer                    , dimension(:), intent(in   ), optional :: radiationType
    integer                                                            :: iComponent

    ! Loop over all radiation components.
    Radiation_Temperature=0.0d0
    do iComponent=1,size(radiation%radiationType)
       ! Call the appropriate routine to get the temperature.
       !# <include directive="radiationTemperature" type="functionCall" functionType="void">
       !#  <functionArgs>radiation%radiationType(iComponent),radiationType#label,radiation%components(iComponent)%properties,Radiation_Temperature,radiationType</functionArgs>
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
    class           (radiationStructure)              , intent(in   )           :: radiation
    double precision                                  , intent(in   )           :: wavelength
    integer                             , dimension(:), intent(in   ), optional :: radiationType
    integer                                                                     :: iComponent

    ! Loop over all radiation components.
    Radiation_Flux=0.0d0
    do iComponent=1,size(radiation%radiationType)
       ! Call the appropriate routine to get the flux.
       !# <include directive="radiationFlux" type="functionCall" functionType="void">
       !#  <functionArgs>radiation%radiationType(iComponent),radiationType#label,radiation%components(iComponent)%properties,wavelength,Radiation_Flux,radiationType</functionArgs>
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
    use Numerical_Constants_Units
    use Numerical_Constants_Physical
    implicit none
    class           (radiationStructure        )              , intent(in   ) :: radiation
    double precision                            , dimension(2), intent(in   ) :: wavelengthRange
    double precision                            , external                    :: crossSectionFunction
    type            (c_ptr                     )                              :: parameterPointer
    type            (fgsl_function             )                              :: integrandFunction
    type            (fgsl_integration_workspace)                              :: integrationWorkspace

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
    real(kind=c_double)        :: Cross_Section_Integrand
    real(kind=c_double), value :: wavelength
    type(c_ptr        ), value :: parameterPointer

    if (wavelength > 0.0d0) then
       Cross_Section_Integrand=crossSectionFunctionGlobal(wavelength)*Radiation_Flux(radiationGlobal,wavelength)/wavelength
    else
       Cross_Section_Integrand=0.0d0
    end if
    return
  end function Cross_Section_Integrand

end module Radiation_Structure
