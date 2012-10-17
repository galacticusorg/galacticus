!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements an intergalatic background (excluding the CMB) radiation component.

module Radiation_Intergalactic_Background
  !% Implements an intergalatic background (excluding the CMB) radiation component.
  use FGSL
  use Tree_Nodes
  implicit none
  private
  public :: Radiation_Set_Intergalactic_Background, Radiation_Temperature_Intergalactic_Background, Radiation_Flux_Intergalactic_Background

  ! Pointer to the functions that actually do the calculation.
  procedure(Radiation_Set_Template ), pointer :: Radiation_Set_Intergalactic_Background_Do  => null()
  procedure(Radiation_Flux_Template), pointer :: Radiation_Flux_Intergalactic_Background_Do => null()
  abstract interface
     subroutine Radiation_Set_Template(thisNode,radiationProperties)
       import treeNode
       type(treeNode),   intent(inout), pointer                   :: thisNode
       double precision, intent(inout), allocatable, dimension(:) :: radiationProperties
     end subroutine Radiation_Set_Template
  end interface
  abstract interface
     subroutine Radiation_Flux_Template(radiationProperties,wavelength,radiationFlux)
       import treeNode
       double precision, intent(in),   dimension(:) :: radiationProperties
       double precision, intent(in)                 :: wavelength
       double precision, intent(inout)              :: radiationFlux
     end subroutine Radiation_Flux_Template
  end interface
  
  ! Flag indicating if module has been initialized.
  logical :: moduleInitialized=.false.

contains

  ! Specify the label for this radiation component.
  !# <radiationLabel>
  !#  <label>IGB</label>
  !# </radiationLabel>

  subroutine Radiation_Initialize_Intergalactic_Background
    !% Initialize the intergalatic background radiation component module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    !# <include directive="radiationIntergalacticBackgroundMethod" type="moduleUse">
    include 'radiation.intergalactic_background.modules.inc'
    !# </include>
    implicit none
    type(varying_string) :: radiationIntergalacticBackgroundMethod

    if (.not.moduleInitialized) then
       !$omp critical(Radiation_Initialize_Intergalactic_Background)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>radiationIntergalacticBackgroundMethod</name>
          !@   <defaultValue>file</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of the intergalatic background radiation field.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('radiationIntergalacticBackgroundMethod',radiationIntergalacticBackgroundMethod,defaultValue='file')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="radiationIntergalacticBackgroundMethod" type="code" action="subroutine">
          !#  <subroutineArgs>radiationIntergalacticBackgroundMethod,Radiation_Set_Intergalactic_Background_Do,Radiation_Flux_Intergalactic_Background_Do</subroutineArgs>
          include 'radiation.intergalactic_background.inc'
          !# </include>
          if (.not.(associated(Radiation_Set_Intergalactic_Background_Do).and.associated(Radiation_Flux_Intergalactic_Background_Do))) call&
               & Galacticus_Error_Report('Radiation_Initialize_Intergalactic_Background','method ' //char(radiationIntergalacticBackgroundMethod)//' is unrecognized')
          moduleInitialized=.true.
       end if
       !$omp end critical(Radiation_Initialize_Intergalactic_Background)
    end if
    return
  end subroutine Radiation_Initialize_Intergalactic_Background

  !# <radiationSet>
  !#  <unitName>Radiation_Set_Intergalactic_Background</unitName>
  !#  <label>IGB</label>
  !# </radiationSet>
  subroutine Radiation_Set_Intergalactic_Background(componentMatched,thisNode,radiationProperties)
    !% Property setting routine for the radiation component from file method.
    use Memory_Management
    implicit none
    logical,          intent(in)                               :: componentMatched
    type(treeNode),   intent(inout), pointer                   :: thisNode
    double precision, intent(inout), allocatable, dimension(:) :: radiationProperties

    ! Return immediately if this component was not matched.
    if (.not.componentMatched) return

    ! Ensure that the module is initialized.
    call Radiation_Initialize_Intergalactic_Background

    ! Call the routine to do the calculation.
    call Radiation_Set_Intergalactic_Background_Do(thisNode,radiationProperties)

    return
  end subroutine Radiation_Set_Intergalactic_Background

  !# <radiationTemperature>
  !#  <unitName>Radiation_Temperature_Intergalactic_Background</unitName>
  !#  <label>IGB</label>
  !# </radiationTemperature>
  subroutine Radiation_Temperature_Intergalactic_Background(requestedType,ourType,radiationProperties,radiationTemperature,radiationType)
    !% Returns the temperature for the radiation component from file method.
    implicit none
    integer,          intent(in)                               :: requestedType,ourType
    double precision, intent(in),    allocatable, dimension(:) :: radiationProperties
    double precision, intent(inout)                            :: radiationTemperature
    integer,          intent(in),    optional,    dimension(:) :: radiationType

    ! Leave temperature undefined.

    return
  end subroutine Radiation_Temperature_Intergalactic_Background

  !# <radiationFlux>
  !#  <unitName>Radiation_Flux_Intergalactic_Background</unitName>
  !#  <label>IGB</label>
  !# </radiationFlux>
  subroutine Radiation_Flux_Intergalactic_Background(requestedType,ourType,radiationProperties,wavelength,radiationFlux,radiationType)
    !% Flux method for the radiation component from file method.
    implicit none
    integer,          intent(in)                                :: requestedType,ourType
    double precision, intent(in)                                :: wavelength
    double precision, intent(in),   allocatable, dimension(:)   :: radiationProperties
    double precision, intent(inout)                             :: radiationFlux
    integer,          intent(in),   optional,    dimension(:)   :: radiationType
 
    ! Return immediately if this component was not matched.
    if (requestedType /= ourType) return

    ! Return immediately if the radiation object is not initialized.
    if (.not.allocated(radiationProperties)) return

    ! If specific radiation types were requested, check to see if they match our type.
    if (present(radiationType)) then
       if (all(radiationType /= ourType)) return
    end if

    ! Ensure that the module is initialized.
    call Radiation_Initialize_Intergalactic_Background

    ! Call the routine to do the calculation.
    call Radiation_Flux_Intergalactic_Background_Do(radiationProperties,wavelength,radiationFlux)

    return
  end subroutine Radiation_Flux_Intergalactic_Background

end module Radiation_Intergalactic_Background
