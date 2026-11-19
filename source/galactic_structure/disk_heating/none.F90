!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implementation of a galactic disk heating class which applies no heating.
  !!}

  !![
  <galacticDiskHeating name="galacticDiskHeatingNone" docformat="rst">
    <description>
    A galactic disk heating class which applies no heating: the mass distributions of galactic disks are unaffected by tidal
    heating.
    </description>
  </galacticDiskHeating>
  !!]
  type, extends(galacticDiskHeatingClass) :: galacticDiskHeatingNone
     !!{RST
     A galactic disk heating class which applies no heating.
     !!}
     private
   contains
     procedure :: isActive => noneIsActive
     procedure :: heating  => noneHeating
  end type galacticDiskHeatingNone

  interface galacticDiskHeatingNone
     !!{RST
     Constructors for the :galacticus-class:`galacticDiskHeatingNone` galactic disk heating class.
     !!}
     module procedure noneConstructorParameters
  end interface galacticDiskHeatingNone

contains

  function noneConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticDiskHeatingNone` galactic disk heating class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticDiskHeatingNone)                :: self
    type(inputParameters        ), intent(inout) :: parameters

    self=galacticDiskHeatingNone()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function noneConstructorParameters

  logical function noneIsActive(self) result(isActive)
    !!{RST
    Return false, as this class applies no heating.
    !!}
    implicit none
    class(galacticDiskHeatingNone), intent(inout) :: self
    !$GLC attributes unused :: self

    isActive=.false.
    return
  end function noneIsActive

  double precision function noneHeating(self,tidalHeatingNormalized,radius,velocity) result(heating)
    !!{RST
    Return zero heating.
    !!}
    implicit none
    class           (galacticDiskHeatingNone), intent(inout) :: self
    double precision                         , intent(in   ) :: tidalHeatingNormalized, radius, &
         &                                                      velocity
    !$GLC attributes unused :: self, tidalHeatingNormalized, radius, velocity

    heating=0.0d0
    return
  end function noneHeating
