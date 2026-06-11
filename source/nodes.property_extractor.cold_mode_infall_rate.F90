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
Implements a cold mode infall rate property extractor class.
!!}

  use :: Cooling_Cold_Mode_Infall_Rates, only : coldModeInfallRate, coldModeInfallRateClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRateInfallColdMode" docformat="rst">
   <description>
   Extracts the infall rate of cold-mode gas accreting onto a galaxy, capturing the filamentary cold-stream accretion channel that bypasses virial shock heating and directly feeds star formation. Particularly relevant for high-redshift galaxies in massive halos.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRateInfallColdMode
     !!{RST
     A cold mode infall rate property extractor class.
     !!}
     private
     class(coldModeInfallRateClass), pointer :: coldModeInfallRate_ => null()
   contains
     final     ::                rateInfallColdModeDestructor
     procedure :: extract     => rateInfallColdModeExtract
     procedure :: name        => rateInfallColdModeName
     procedure :: description => rateInfallColdModeDescription
     procedure :: unitsInSI   => rateInfallColdModeUnitsInSI
     procedure :: units       => rateInfallColdModeUnits
  end type nodePropertyExtractorRateInfallColdMode

  interface nodePropertyExtractorRateInfallColdMode
     !!{RST
     Constructors for the ``nodePropertyExtractorRateInfallColdMode`` property extractor class.
     !!}
     module procedure rateInfallColdModeConstructorParameters
     module procedure rateInfallColdModeConstructorInternal
  end interface nodePropertyExtractorRateInfallColdMode

contains

  function rateInfallColdModeConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``nodePropertyExtractorRateInfallColdMode`` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRateInfallColdMode)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(coldModeInfallRateClass                ), pointer       :: coldModeInfallRate_

    !![
    <objectBuilder class="coldModeInfallRate" name="coldModeInfallRate_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRateInfallColdMode(coldModeInfallRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coldModeInfallRate_"/>
    !!]
    return
  end function rateInfallColdModeConstructorParameters

  function rateInfallColdModeConstructorInternal(coldModeInfallRate_) result(self)
    !!{RST
    Internal constructor for the ``nodePropertyExtractorRateInfallColdMode`` property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRateInfallColdMode)                        :: self
    class(coldModeInfallRateClass                ), intent(in   ), target :: coldModeInfallRate_
    !![
    <constructorAssign variables="*coldModeInfallRate_"/>
    !!]

    return
  end function rateInfallColdModeConstructorInternal

  subroutine rateInfallColdModeDestructor(self)
    !!{RST
    Destructor for the ``nodePropertyExtractorRateInfallColdMode`` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRateInfallColdMode), intent(inout) :: self

    !![
    <objectDestructor name="self%coldModeInfallRate_"/>
    !!]
    return
  end subroutine rateInfallColdModeDestructor

  double precision function rateInfallColdModeExtract(self,node,instance)
    !!{RST
    Implement a ``rateInfallColdMode`` property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nodePropertyExtractorRateInfallColdMode), intent(inout), target   :: self
    type (treeNode                               ), intent(inout), target   :: node
    type (multiCounter                           ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    rateInfallColdModeExtract=self%coldModeInfallRate_%infallRate(node)
    return
  end function rateInfallColdModeExtract

  function rateInfallColdModeName(self)
    !!{RST
    Return the name of the ``rateInfallColdMode`` property.
    !!}
    implicit none
    type (varying_string                         )                :: rateInfallColdModeName
    class(nodePropertyExtractorRateInfallColdMode), intent(inout) :: self
    !$GLC attributes unused :: self

    rateInfallColdModeName=var_str('coldModeInfallRate')
    return
  end function rateInfallColdModeName

  function rateInfallColdModeDescription(self)
    !!{RST
    Return a description of the ``rateInfallColdMode`` property.
    !!}
    implicit none
    type (varying_string                         )                :: rateInfallColdModeDescription
    class(nodePropertyExtractorRateInfallColdMode), intent(inout) :: self
    !$GLC attributes unused :: self

    rateInfallColdModeDescription=var_str('Rate of infall of cold mode material onto the galaxy [M☉/Gyr].')
    return
  end function rateInfallColdModeDescription

  double precision function rateInfallColdModeUnitsInSI(self)
    !!{RST
    Return the units of the ``rateInfallColdMode`` property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, massSolar
    implicit none
    class(nodePropertyExtractorRateInfallColdMode), intent(inout) :: self
    !$GLC attributes unused :: self

    rateInfallColdModeUnitsInSI=massSolar/gigaYear
    return
  end function rateInfallColdModeUnitsInSI

  function rateInfallColdModeUnits(self) result(units)
    !!{RST
    Return the units of the rateInfallColdMode property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                               )                :: units
    class(nodePropertyExtractorRateInfallColdMode), intent(inout) :: self

    units=unitType(self%unitsInSI(),description='M☉/Gyr',quantity='solMass/Gyr')
    return
  end function rateInfallColdModeUnits
