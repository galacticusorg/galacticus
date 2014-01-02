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
    use Galacticus_Nodes
    use Memory_Management
    use Cosmology_Functions
    implicit none
    logical                                                             , intent(in   )          :: componentMatched
    type            (treeNode               )                           , intent(inout), pointer :: thisNode
    double precision                         , allocatable, dimension(:), intent(inout)          :: radiationProperties
    class           (nodeComponentBasic     )             , pointer                              :: thisBasicComponent
    class           (cosmologyFunctionsClass)             , pointer                              :: cosmologyFunctionsDefault

    ! Return immediately if this component was not matched.
    if (.not.componentMatched) return

    ! Ensure that the properties array is allocated.
    if (.not.allocated(radiationProperties)) call Alloc_Array(radiationProperties,[1])

    ! Get the default cosmology functions object.
    cosmologyFunctionsDefault => cosmologyFunctions()
    ! Set the CMB temperature.
    thisBasicComponent => thisNode%basic()
    radiationProperties(1)=cosmologyFunctionsDefault%temperatureCMBEpochal(thisBasicComponent%time())

    return
  end subroutine Radiation_Set_CMB

  !# <radiationTemperature>
  !#  <unitName>Radiation_Temperature_CMB</unitName>
  !#  <label>CMB</label>
  !# </radiationTemperature>
  subroutine Radiation_Temperature_CMB(requestedType,ourType,radiationProperties,radiationTemperature,radiationType)
    !% Returns the temperature for the cosmic microwave background radiation component.
    implicit none
    integer                                    , intent(in   )           :: ourType             , requestedType
    double precision, allocatable, dimension(:), intent(in   )           :: radiationProperties
    double precision                           , intent(inout)           :: radiationTemperature
    integer                      , dimension(:), intent(in   ), optional :: radiationType

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
    implicit none
    integer                                    , intent(in   )           :: ourType            , requestedType
    double precision                           , intent(in   )           :: wavelength
    double precision, allocatable, dimension(:), intent(in   )           :: radiationProperties
    double precision                           , intent(inout)           :: radiationFlux
    integer                      , dimension(:), intent(in   ), optional :: radiationType

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
