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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a galactic inclination class which supplies no inclination.
  !!}

  !![
  <galacticInclination name="galacticInclinationNull" docformat="rst">
   <description>
   A galactic inclination class which supplies no inclination, and reports itself unavailable. This is the default:
   most models have no need of orientation, and should pay nothing for it.

   Consumers which depend on orientation must test ``isAvailable`` at construction and report a clear error. For dust
   attenuation the alternative to an error is to average over orientation instead, which
   :galacticus-class:`dustAttenuationInclinationAveraged` does.
   </description>
  </galacticInclination>
  !!]
  type, extends(galacticInclinationClass) :: galacticInclinationNull
     !!{RST
     A galactic inclination class which supplies no inclination.
     !!}
     private
   contains
     procedure :: inclination => nullInclination
     procedure :: isAvailable => nullIsAvailable
  end type galacticInclinationNull

  interface galacticInclinationNull
     !!{RST
     Constructors for the :galacticus-class:`galacticInclinationNull` galactic inclination class.
     !!}
     module procedure nullConstructorParameters
  end interface galacticInclinationNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticInclinationNull` galactic inclination class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticInclinationNull)                :: self
    type(inputParameters        ), intent(inout) :: parameters

    self=galacticInclinationNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  double precision function nullInclination(self,node) result(inclination)
    !!{RST
    Return no inclination. Calling this is an error: a consumer must test ``isAvailable`` first.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(galacticInclinationNull), intent(inout)         :: self
    type (treeNode               ), intent(inout), target :: node
    !$GLC attributes unused :: self, node

    inclination=0.0d0
    call Error_Report('no inclination is available: this class supplies none, as its `isAvailable` method reports'//{introspection:location})
    return
  end function nullInclination

  logical function nullIsAvailable(self) result(isAvailable)
    !!{RST
    Return false: this class supplies no inclination.
    !!}
    implicit none
    class(galacticInclinationNull), intent(inout) :: self
    !$GLC attributes unused :: self

    isAvailable=.false.
    return
  end function nullIsAvailable
