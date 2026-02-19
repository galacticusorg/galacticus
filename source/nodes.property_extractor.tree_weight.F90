!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
Implements a merger tree weight property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorTreeWeight">
   <description>A merger tree weight property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorTreeWeight
     !!{
     A merger tree weight property extractor class.
     !!}
     private
   contains
     procedure :: extract     => treeWeightExtract
     procedure :: name        => treeWeightName
     procedure :: description => treeWeightDescription
     procedure :: unitsInSI   => treeWeightUnitsInSI
  end type nodePropertyExtractorTreeWeight

  interface nodePropertyExtractorTreeWeight
     !!{
     Constructors for the \refClass{nodePropertyExtractorTreeWeight} output analysis class.
     !!}
     module procedure treeWeightConstructorParameters
  end interface nodePropertyExtractorTreeWeight

contains

  function treeWeightConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorTreeWeight} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorTreeWeight)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=nodePropertyExtractorTreeWeight()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function treeWeightConstructorParameters

  double precision function treeWeightExtract(self,node,instance)
    !!{
    Implement a last isolated redshift output analysis.
    !!}
    implicit none
    class(nodePropertyExtractorTreeWeight), intent(inout), target   :: self
    type (treeNode                       ), intent(inout), target   :: node
    type (multiCounter                   ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance

    treeWeightExtract=node%hostTree%volumeWeight
    return
  end function treeWeightExtract

  function treeWeightName(self)
    !!{
    Return the name of the last isolated redshift property.
    !!}
    implicit none
    type (varying_string                 )                :: treeWeightName
    class(nodePropertyExtractorTreeWeight), intent(inout) :: self
    !$GLC attributes unused :: self

    treeWeightName=var_str('mergerTreeWeight')
    return
  end function treeWeightName

  function treeWeightDescription(self)
    !!{
    Return a description of the treeWeight property.
    !!}
    implicit none
    type (varying_string                 )                :: treeWeightDescription
    class(nodePropertyExtractorTreeWeight), intent(inout) :: self
    !$GLC attributes unused :: self

    treeWeightDescription=var_str('The weight assigned to this tree - typically the number of such trees per unit volume required to make a cosmologically-representative sample [Mpc⁻³].')
    return
  end function treeWeightDescription

  double precision function treeWeightUnitsInSI(self)
    !!{
    Return the units of the last isolated redshift property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorTreeWeight), intent(inout) :: self
    !$GLC attributes unused :: self

    treeWeightUnitsInSI=1.0d0/megaParsec**3
    return
  end function treeWeightUnitsInSI

