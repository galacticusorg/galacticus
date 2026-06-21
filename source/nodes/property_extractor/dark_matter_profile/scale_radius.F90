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
  <nodePropertyExtractor name="nodePropertyExtractorDarkMatterProfileScaleRadius" docformat="rst">
   <description>
   Extracts the scale radius of the dark matter halo density profile (e.g., the NFW scale radius :math:`r_\mathrm{s}` where the logarithmic slope equals :math:`-2`), a key structural parameter relating halo mass to its spatial extent and concentration.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorDarkMatterProfileScaleRadius
     !!{RST
     A dark matter profile scale radius output property extractor class.
     !!}
     private
   contains
     procedure :: extract     => darkMatterProfileScaleRadiusExtract
     procedure :: name        => darkMatterProfileScaleRadiusName
     procedure :: description => darkMatterProfileScaleRadiusDescription
     procedure :: unitsInSI   => darkMatterProfileScaleRadiusUnitsInSI
     procedure :: units       => darkMatterProfileScaleRadiusUnits
  end type nodePropertyExtractorDarkMatterProfileScaleRadius

  interface nodePropertyExtractorDarkMatterProfileScaleRadius
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorDarkMatterProfileScaleRadius` property extractor class.
     !!}
     module procedure darkMatterProfileScaleRadiusConstructorParameters
  end interface nodePropertyExtractorDarkMatterProfileScaleRadius

contains

  function darkMatterProfileScaleRadiusConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorDarkMatterProfileScaleRadius` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorDarkMatterProfileScaleRadius)                :: self
    type (inputParameters                                  ), intent(inout) :: parameters
    
    self=nodePropertyExtractorDarkMatterProfileScaleRadius()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkMatterProfileScaleRadiusConstructorParameters

  double precision function darkMatterProfileScaleRadiusExtract(self,node,instance)
    !!{RST
    Implement a ``darkMatterProfileScaleRadius`` output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class(nodePropertyExtractorDarkMatterProfileScaleRadius), intent(inout), target   :: self
    type (treeNode                                         ), intent(inout), target   :: node
    type (multiCounter                                     ), intent(inout), optional :: instance
    class(nodeComponentDarkMatterProfile                   ), pointer                 :: darkMatterProfile
    !$GLC attributes unused :: self, instance

    darkMatterProfile                   => node             %darkMatterProfile()
    darkMatterProfileScaleRadiusExtract =  darkMatterProfile%scale            ()
    return
  end function darkMatterProfileScaleRadiusExtract


  function darkMatterProfileScaleRadiusName(self)
    !!{RST
    Return the name of the ``darkMatterProfileScaleRadius`` property.
    !!}
    implicit none
    type (varying_string                                   )                :: darkMatterProfileScaleRadiusName
    class(nodePropertyExtractorDarkMatterProfileScaleRadius), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileScaleRadiusName=var_str('darkMatterProfileScaleRadius')
    return
  end function darkMatterProfileScaleRadiusName

  function darkMatterProfileScaleRadiusDescription(self)
    !!{RST
    Return a description of the ``darkMatterProfileScaleRadius`` property.
    !!}
    implicit none
    type (varying_string                                   )                :: darkMatterProfileScaleRadiusDescription
    class(nodePropertyExtractorDarkMatterProfileScaleRadius), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileScaleRadiusDescription=var_str('The scale radius of the dark matter profile.')
    return
  end function darkMatterProfileScaleRadiusDescription

  double precision function darkMatterProfileScaleRadiusUnitsInSI(self)
    !!{RST
    Return the units of the ``darkMatterProfileScaleRadius`` property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorDarkMatterProfileScaleRadius), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileScaleRadiusUnitsInSI=megaParsec
    return
  end function darkMatterProfileScaleRadiusUnitsInSI

  function darkMatterProfileScaleRadiusUnits(self) result(units)
    !!{RST
    Return the units of the darkMatterProfileScaleRadius property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                                         )                :: units
    class(nodePropertyExtractorDarkMatterProfileScaleRadius), intent(inout) :: self

    units=unitType(self%unitsInSI(),description='Mpc',quantity='Mpc')
    return
  end function darkMatterProfileScaleRadiusUnits
