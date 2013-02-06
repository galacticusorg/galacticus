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

module Radiation_Null 
  implicit none
  private
  public :: Radiation_Set_Null, Radiation_Temperature_Null, Radiation_Flux_Null

contains

  !# <radiationLabel>
  !#  <label>Null</label>
  !# </radiationLabel>

  !# <radiationSet>
  !#  <unitName>Radiation_Set_Null</unitName>
  !#  <label>Null</label>
  !# </radiationSet>
  subroutine Radiation_Set_Null(componentMatched,thisNode,radiationProperties)
    !% Property setting routine for null radiation component.
    use Galacticus_Nodes
    implicit none
    logical,          intent(in)                               :: componentMatched
    type(treeNode),   intent(inout), pointer                   :: thisNode
    double precision, intent(inout), allocatable, dimension(:) :: radiationProperties

    return
  end subroutine Radiation_Set_Null

  !# <radiationTemperature>
  !#  <unitName>Radiation_Temperature_Null</unitName>
  !#  <label>Null</label>
  !# </radiationTemperature>
  subroutine Radiation_Temperature_Null(requestedType,ourType,radiationProperties,radiationTemperature,radiationType)
    !% Temperature method for the null radiation component.
    implicit none
    integer,          intent(in)                               :: requestedType,ourType
    double precision, intent(in),    allocatable, dimension(:) :: radiationProperties
    double precision, intent(inout)                            :: radiationTemperature
    integer,          intent(in),    optional,    dimension(:) :: radiationType

    return
  end subroutine Radiation_Temperature_Null

  !# <radiationFlux>
  !#  <unitName>Radiation_Flux_Null</unitName>
  !#  <label>Null</label>
  !# </radiationFlux>
  subroutine Radiation_Flux_Null(requestedType,ourType,radiationProperties,wavelength,radiationFlux,radiationType)
    !% Flux method for the null radiation component.
    implicit none
    integer,          intent(in)                               :: requestedType,ourType
    double precision, intent(in)                               :: wavelength
    double precision, intent(in),    allocatable, dimension(:) :: radiationProperties
    double precision, intent(inout)                            :: radiationFlux
    integer,          intent(in),    optional,    dimension(:) :: radiationType

    return
  end subroutine Radiation_Flux_Null

end module Radiation_Null
