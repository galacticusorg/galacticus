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
  Implements a node property extractor for the inclination of a galaxy.
  !!}

  use :: Galactic_Inclinations, only : galacticInclinationClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorInclination" docformat="rst">
   <description>
   A node property extractor which reports the inclination of a galaxy to the line of sight, in degrees, with zero
   being face-on and 90 edge-on.

   Whatever consumes inclination---dust attenuation above all---does so through the same
   :galacticus-class:`galacticInclinationClass` object, so this reports the angle that was actually applied, and
   makes it reproducible rather than hidden inside the attenuation.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorInclination
     !!{RST
     A node property extractor for the inclination of a galaxy.
     !!}
     private
     class(galacticInclinationClass), pointer :: galacticInclination_ => null()
   contains
     final     ::                inclinationDestructor
     procedure :: extract     => inclinationExtract
     procedure :: name        => inclinationName
     procedure :: description => inclinationDescription
     procedure :: unitsInSI   => inclinationUnitsInSI
     procedure :: units       => inclinationUnits
  end type nodePropertyExtractorInclination

  interface nodePropertyExtractorInclination
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorInclination` property extractor class.
     !!}
     module procedure inclinationConstructorParameters
     module procedure inclinationConstructorInternal
  end interface nodePropertyExtractorInclination

contains

  function inclinationConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorInclination` property extractor class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorInclination)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(galacticInclinationClass        ), pointer       :: galacticInclination_

    !![
    <objectBuilder class="galacticInclination" name="galacticInclination_" source="parameters"/>
    !!]
    self=nodePropertyExtractorInclination(galacticInclination_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticInclination_"/>
    !!]
    return
  end function inclinationConstructorParameters

  function inclinationConstructorInternal(galacticInclination_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorInclination` property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (nodePropertyExtractorInclination)                        :: self
    class(galacticInclinationClass        ), intent(in   ), target :: galacticInclination_
    !![
    <constructorAssign variables="*galacticInclination_"/>
    !!]

    if (.not.self%galacticInclination_%isAvailable())                                                               &
         & call Error_Report(                                                                                       &
         &                   'no inclination is available to extract: set `galacticInclination` to a class which'// &
         &                   ' supplies one'                                                                     // &
         &                   {introspection:location}                                                               &
         &                  )
    return
  end function inclinationConstructorInternal

  subroutine inclinationDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorInclination` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorInclination), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticInclination_"/>
    !!]
    return
  end subroutine inclinationDestructor

  double precision function inclinationExtract(self,node,instance)
    !!{RST
    Return the inclination of the galaxy, in degrees.
    !!}
    use :: Numerical_Constants_Astronomical, only : degreesToRadians
    implicit none
    class(nodePropertyExtractorInclination), intent(inout), target   :: self
    type (treeNode                        ), intent(inout), target   :: node
    type (multiCounter                    ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    inclinationExtract=+self%galacticInclination_%inclination     (node) &
         &             /                          degreesToRadians
    return
  end function inclinationExtract

  function inclinationName(self)
    !!{RST
    Return the name of the inclination property.
    !!}
    implicit none
    type (varying_string                  )                :: inclinationName
    class(nodePropertyExtractorInclination), intent(inout) :: self
    !$GLC attributes unused :: self

    inclinationName=var_str('inclination')
    return
  end function inclinationName

  function inclinationDescription(self)
    !!{RST
    Return a description of the inclination property.
    !!}
    implicit none
    type (varying_string                  )                :: inclinationDescription
    class(nodePropertyExtractorInclination), intent(inout) :: self
    !$GLC attributes unused :: self

    inclinationDescription=var_str('Inclination of the galaxy to the line of sight [degrees; 0 is face-on].')
    return
  end function inclinationDescription

  double precision function inclinationUnitsInSI(self)
    !!{RST
    Return the units of the inclination property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : degreesToRadians
    implicit none
    class(nodePropertyExtractorInclination), intent(inout) :: self
    !$GLC attributes unused :: self

    ! Radians are dimensionless in SI, so the conversion factor is simply that from degrees to radians.
    inclinationUnitsInSI=degreesToRadians
    return
  end function inclinationUnitsInSI

  function inclinationUnits(self) result(units)
    !!{RST
    Return the units of the inclination property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                        )                :: units
    class(nodePropertyExtractorInclination), intent(inout) :: self

    units=unitType(self%unitsInSI(),description='degrees',quantity='angle')
    return
  end function inclinationUnits
