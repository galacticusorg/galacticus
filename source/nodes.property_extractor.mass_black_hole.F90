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
Implements a black hole mass property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassBlackHole" docformat="rst">
   <description>
   Extracts the mass of the central supermassive black hole in a galaxy node, accounting for growth via gas accretion and black hole mergers, for comparison with :math:`M_\bullet`--:math:`\sigma` and :math:`M_\bullet`--:math:`M_\mathrm{bulge}` scaling relations.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassBlackHole
     !!{RST
     A black hole mass property extractor class.
     !!}
     private
   contains
     procedure :: extract     => massBlackHoleExtract
     procedure :: quantity    => massBlackHoleQuantity
     procedure :: name        => massBlackHoleName
     procedure :: description => massBlackHoleDescription
     procedure :: unitsInSI   => massBlackHoleUnitsInSI
     procedure :: units       => massBlackHoleUnits
  end type nodePropertyExtractorMassBlackHole

  interface nodePropertyExtractorMassBlackHole
     !!{RST
     Constructors for the ``nodePropertyExtractorMassBlackHole`` property extractor class.
     !!}
     module procedure massBlackHoleConstructorParameters
  end interface nodePropertyExtractorMassBlackHole

contains

  function massBlackHoleConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``nodePropertyExtractorMassBlackHole`` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorMassBlackHole)                :: self
    type(inputParameters                   ), intent(inout) :: parameters

    self=nodePropertyExtractorMassBlackHole()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massBlackHoleConstructorParameters

  double precision function massBlackHoleExtract(self,node,instance)
    !!{RST
    Implement a massBlackHole output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, treeNode
    implicit none
    class(nodePropertyExtractorMassBlackHole), intent(inout), target   :: self
    type (treeNode                          ), intent(inout), target   :: node
    type (multiCounter                      ), intent(inout), optional :: instance
    class(nodeComponentBlackHole            ), pointer                 :: blackHole
    !$GLC attributes unused :: self, instance

    blackHole            => node     %blackHole()
    massBlackHoleExtract =  blackHole%mass     ()
    return
  end function massBlackHoleExtract


  function massBlackHoleQuantity(self)
    !!{RST
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityMass
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: massBlackHoleQuantity
    class(nodePropertyExtractorMassBlackHole           ), intent(inout) :: self
    !$GLC attributes unused :: self

    massBlackHoleQuantity=outputAnalysisPropertyQuantityMass
    return
  end function massBlackHoleQuantity

  function massBlackHoleName(self)
    !!{RST
    Return the name of the massBlackHole property.
    !!}
    implicit none
    type (varying_string                    )                :: massBlackHoleName
    class(nodePropertyExtractorMassBlackHole), intent(inout) :: self
    !$GLC attributes unused :: self

    massBlackHoleName=var_str('massBlackHole')
    return
  end function massBlackHoleName

  function massBlackHoleDescription(self)
    !!{RST
    Return a description of the massBlackHole property.
    !!}
    implicit none
    type (varying_string                    )                :: massBlackHoleDescription
    class(nodePropertyExtractorMassBlackHole), intent(inout) :: self
    !$GLC attributes unused :: self

    massBlackHoleDescription=var_str('The mass of the central black hole in each galaxy.')
    return
  end function massBlackHoleDescription

  double precision function massBlackHoleUnitsInSI(self)
    !!{RST
    Return the units of the black hole mass property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassBlackHole), intent(inout) :: self
    !$GLC attributes unused :: self

    massBlackHoleUnitsInSI=massSolar
    return
  end function massBlackHoleUnitsInSI

  function massBlackHoleUnits(self) result(units)
    !!{RST
    Return the units of the black hole mass property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                          )                :: units
    class(nodePropertyExtractorMassBlackHole), intent(inout) :: self

    units=unitType(self%unitsInSI(),description='Solar masses',quantity='solMass')
    return
  end function massBlackHoleUnits
