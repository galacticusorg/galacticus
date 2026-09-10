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
  Implements a galactic inclination class in which every galaxy has the same inclination.
  !!}

  !![
  <galacticInclination name="galacticInclinationFixed" docformat="rst">
   <description>
   A galactic inclination class in which every galaxy is given the same inclination, set by the ``inclination``
   parameter (in degrees).

   This exists for testing, and for face-on versus edge-on comparisons: running the same model at
   ``inclination``:math:`=0` and :math:`=90` isolates the effect of orientation on an anisotropic process without
   any of the scatter that a distribution of orientations introduces.
   </description>
  </galacticInclination>
  !!]
  type, extends(galacticInclinationClass) :: galacticInclinationFixed
     !!{RST
     A galactic inclination class in which every galaxy has the same inclination.
     !!}
     private
     ! Stored in degrees, as the parameter is given, so that a descriptor round-trips in the same units.
     double precision :: inclination_
   contains
     procedure :: inclination => fixedInclination
  end type galacticInclinationFixed

  interface galacticInclinationFixed
     !!{RST
     Constructors for the :galacticus-class:`galacticInclinationFixed` galactic inclination class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface galacticInclinationFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticInclinationFixed` galactic inclination class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticInclinationFixed)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    double precision                                          :: inclination

    !![
    <inputParameter docformat="rst">
      <name>inclination</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The inclination, in degrees, assigned to every galaxy. Zero is face-on and 90 edge-on.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=galacticInclinationFixed(inclination)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(inclination_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`galacticInclinationFixed` galactic inclination class. The
    inclination is given in degrees.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (galacticInclinationFixed)                :: self
    double precision                          , intent(in   ) :: inclination_
    !![
    <constructorAssign variables="inclination_"/>
    !!]

    if (self%inclination_ < 0.0d0 .or. self%inclination_ > 90.0d0) call Error_Report('inclination must be in the range [0,90] degrees'//{introspection:location})
    return
  end function fixedConstructorInternal

  double precision function fixedInclination(self,node) result(inclination)
    !!{RST
    Return the fixed inclination, converted from the degrees in which it is stored to the radians in which the class
    interface reports it.
    !!}
    use :: Numerical_Constants_Astronomical, only : degreesToRadians
    implicit none
    class(galacticInclinationFixed), intent(inout)         :: self
    type (treeNode                ), intent(inout), target :: node
    !$GLC attributes unused :: node

    inclination=+self%inclination_    &
         &      *     degreesToRadians
    return
  end function fixedInclination
