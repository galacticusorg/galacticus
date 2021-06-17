!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a half-stellar mass radius output analysis property extractor class.
!!}


  !![
  <nodePropertyExtractor name="nodePropertyExtractorHalfMassRadius">
   <description>A half-(stellar) mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorHalfMassRadius
     !!{
     A half-(stellar) mass property extractor output analysis class.
     !!}
     private
   contains
     procedure :: extract     => halfMassRadiusExtract
     procedure :: type        => halfMassRadiusType
     procedure :: name        => halfMassRadiusName
     procedure :: description => halfMassRadiusDescription
     procedure :: unitsInSI   => halfMassRadiusUnitsInSI
  end type nodePropertyExtractorHalfMassRadius

  interface nodePropertyExtractorHalfMassRadius
     !!{
     Constructors for the ``halfMassRadius'' output analysis class.
     !!}
     module procedure halfMassRadiusConstructorParameters
  end interface nodePropertyExtractorHalfMassRadius

contains

  function halfMassRadiusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``halfMassRadius'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorHalfMassRadius)                :: self
    type(inputParameters                    ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    ! Build the object.
    self=nodePropertyExtractorHalfMassRadius()
    return
  end function halfMassRadiusConstructorParameters

  double precision function halfMassRadiusExtract(self,node,instance)
    !!{
    Implement a half-mass output analysis.
    !!}
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Radius_Enclosing_Mass
    use :: Galactic_Structure_Options        , only : massTypeStellar
    implicit none
    class(nodePropertyExtractorHalfMassRadius), intent(inout)           :: self
    type (treeNode                           ), intent(inout), target   :: node
    type (multiCounter                       ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance

    halfMassRadiusExtract=Galactic_Structure_Radius_Enclosing_Mass(node,fractionalMass=0.5d0,massType=massTypeStellar)
    return
  end function halfMassRadiusExtract

  integer function halfMassRadiusType(self)
    !!{
    Return the type of the half-mass radius property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorHalfMassRadius), intent(inout) :: self
    !$GLC attributes unused :: self

    halfMassRadiusType=outputAnalysisPropertyTypeLinear
    return
  end function halfMassRadiusType

  function halfMassRadiusName(self)
    !!{
    Return the name of the halfMassRadius property.
    !!}
    implicit none
    type (varying_string                     )                :: halfMassRadiusName
    class(nodePropertyExtractorHalfMassRadius), intent(inout) :: self
    !$GLC attributes unused :: self

    halfMassRadiusName=var_str('radiusHalfMass')
    return
  end function halfMassRadiusName

  function halfMassRadiusDescription(self)
    !!{
    Return a description of the halfMassRadius property.
    !!}
    implicit none
    type (varying_string                     )                :: halfMassRadiusDescription
    class(nodePropertyExtractorHalfMassRadius), intent(inout) :: self
    !$GLC attributes unused :: self

    halfMassRadiusDescription=var_str('The stellar half-mass radius.')
    return
  end function halfMassRadiusDescription

  double precision function halfMassRadiusUnitsInSI(self)
    !!{
    Return the units of the halfMassRadius property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorHalfMassRadius), intent(inout) :: self
    !$GLC attributes unused :: self

    halfMassRadiusUnitsInSI=megaParsec
    return
  end function halfMassRadiusUnitsInSI
