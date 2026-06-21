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

!!{RST
Implements a dark matter profile scale radius output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDarkMatterProfileShapeParameter" docformat="rst">
   <description>
   Extracts the shape parameter of the dark matter halo density profile (e.g., the Einasto index or inner logarithmic slope), which characterizes the curvature of the profile and distinguishes between different dark matter density models such as NFW, Einasto, or cored profiles.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorDarkMatterProfileShapeParameter
     !!{RST
     A dark matter profile scale radius output property extractor class.
     !!}
     private
   contains
     procedure :: extract     => darkMatterProfileShapeParameterExtract
     procedure :: name        => darkMatterProfileShapeParameterName
     procedure :: description => darkMatterProfileShapeParameterDescription
     procedure :: unitsInSI   => darkMatterProfileShapeParameterUnitsInSI
     procedure :: units       => darkMatterProfileShapeParameterUnits
  end type nodePropertyExtractorDarkMatterProfileShapeParameter

  interface nodePropertyExtractorDarkMatterProfileShapeParameter
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorDarkMatterProfileShapeParameter` property extractor class.
     !!}
     module procedure darkMatterProfileShapeParameterConstructorParameters
  end interface nodePropertyExtractorDarkMatterProfileShapeParameter

contains

  function darkMatterProfileShapeParameterConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorDarkMatterProfileShapeParameter` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorDarkMatterProfileShapeParameter)                :: self
    type (inputParameters                                  ), intent(inout) :: parameters
    
    self=nodePropertyExtractorDarkMatterProfileShapeParameter()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkMatterProfileShapeParameterConstructorParameters

  double precision function darkMatterProfileShapeParameterExtract(self,node,instance)
    !!{RST
    Implement a ``darkMatterProfileShapeParameter`` output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class(nodePropertyExtractorDarkMatterProfileShapeParameter), intent(inout), target   :: self
    type (treeNode                                            ), intent(inout), target   :: node
    type (multiCounter                                        ), intent(inout), optional :: instance
    class(nodeComponentDarkMatterProfile                      ), pointer                 :: darkMatterProfile
    !$GLC attributes unused :: self, instance

    darkMatterProfile                      => node             %darkMatterProfile()
    darkMatterProfileShapeParameterExtract =  darkMatterProfile%shape            ()
    return
  end function darkMatterProfileShapeParameterExtract


  function darkMatterProfileShapeParameterName(self)
    !!{RST
    Return the name of the ``darkMatterProfileShapeParameter`` property.
    !!}
    implicit none
    type (varying_string                                      )                :: darkMatterProfileShapeParameterName
    class(nodePropertyExtractorDarkMatterProfileShapeParameter), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileShapeParameterName=var_str('darkMatterProfileShapeParameter')
    return
  end function darkMatterProfileShapeParameterName

  function darkMatterProfileShapeParameterDescription(self)
    !!{RST
    Return a description of the ``darkMatterProfileShapeParameter`` property.
    !!}
    implicit none
    type (varying_string                                      )                :: darkMatterProfileShapeParameterDescription
    class(nodePropertyExtractorDarkMatterProfileShapeParameter), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileShapeParameterDescription=var_str('The shape parameter of the dark matter profile.')
    return
  end function darkMatterProfileShapeParameterDescription

  double precision function darkMatterProfileShapeParameterUnitsInSI(self)
    !!{RST
    Return the units of the darkMatterProfileShapeParameter property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorDarkMatterProfileShapeParameter), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileShapeParameterUnitsInSI=1.0d0
    return
  end function darkMatterProfileShapeParameterUnitsInSI

  function darkMatterProfileShapeParameterUnits(self) result(units)
    !!{RST
    Return the units of the darkMatterProfileShapeParameter property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                                            )                :: units
    class(nodePropertyExtractorDarkMatterProfileShapeParameter), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(1.0d0)
    return
  end function darkMatterProfileShapeParameterUnits
