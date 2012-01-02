!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements a cosmic microwave background radiation component.

module Radiation_CMB
  !% Implements a cosmic microwave background radiation component.
  implicit none
  private
  public :: Radiation_Set_CMB, Radiation_Temperature_CMB, Radiation_Flux_CMB

contains

  ! Specify the label for this radiation component.
  !# <radiationLabel>
  !#  <label>CMB</label>
  !# </radiationLabel>

  !# <radiationSet>
  !#  <unitName>Radiation_Set_CMB</unitName>
  !#  <label>CMB</label>
  !# </radiationSet>
  subroutine Radiation_Set_CMB(componentMatched,thisNode,radiationProperties)
    !% Property setting routine for the cosmic microwave background radiation component.
    use Tree_Nodes
    use Memory_Management
    use Cosmology_Functions
    implicit none
    logical,          intent(in)                               :: componentMatched
    type(treeNode),   intent(inout), pointer                   :: thisNode
    double precision, intent(inout), allocatable, dimension(:) :: radiationProperties

    ! Return immediately if this component was not matched.
    if (.not.componentMatched) return

    ! Ensure that the properties array is allocated.
    if (.not.allocated(radiationProperties)) call Alloc_Array(radiationProperties,[1])

    ! Set the CMB temperature.
    radiationProperties(1)=CMB_Temperature(Tree_Node_Time(thisNode))

    return
  end subroutine Radiation_Set_CMB

  !# <radiationTemperature>
  !#  <unitName>Radiation_Temperature_CMB</unitName>
  !#  <label>CMB</label>
  !# </radiationTemperature>
  subroutine Radiation_Temperature_CMB(requestedType,ourType,radiationProperties,radiationTemperature,radiationType)
    !% Returns the temperature for the cosmic microwave background radiation component.
    implicit none 
    integer,          intent(in)                               :: requestedType,ourType
    double precision, intent(in),    allocatable, dimension(:) :: radiationProperties
    double precision, intent(inout)                            :: radiationTemperature
    integer,          intent(in),    optional,    dimension(:) :: radiationType

    ! Return immediately if this component was not matched.
    if (requestedType /= ourType) return

    ! Return immediately if the radiation object is not initialized.
    if (.not.allocated(radiationProperties)) return

    ! If specific radiation types were requested, check to see if they match our type.
    if (present(radiationType)) then
       if (all(radiationType /= ourType)) return
    end if

    ! Set the temperature.
    radiationTemperature=radiationProperties(1)

    return
  end subroutine Radiation_Temperature_CMB

  !# <radiationFlux>
  !#  <unitName>Radiation_Flux_CMB</unitName>
  !#  <label>CMB</label>
  !# </radiationFlux>
  subroutine Radiation_Flux_CMB(requestedType,ourType,radiationProperties,wavelength,radiationFlux,radiationType)
    !% Flux method for the CMB radiation component.
    use Thermodynamics_Radiation
    use Numerical_Constants_Units
    use Numerical_Constants_Prefixes
    implicit none
    integer,          intent(in)                              :: requestedType,ourType
    double precision, intent(in)                              :: wavelength
    double precision, intent(in),   allocatable, dimension(:) :: radiationProperties
    double precision, intent(inout)                           :: radiationFlux
    integer,          intent(in),   optional,    dimension(:) :: radiationType

    ! Return immediately if this component was not matched.
    if (requestedType /= ourType) return

    ! Return immediately if the radiation object is not initialized.
    if (.not.allocated(radiationProperties)) return

    ! If specific radiation types were requested, check to see if they match our type.
    if (present(radiationType)) then
       if (all(radiationType /= ourType)) return
    end if

    ! Set the flux.
    radiationFlux=radiationFlux+(centi**2)*Blackbody_Emission(wavelength,radiationProperties(1),radianceTypeFrequency)/ergs

    return
  end subroutine Radiation_Flux_CMB

end module Radiation_CMB
