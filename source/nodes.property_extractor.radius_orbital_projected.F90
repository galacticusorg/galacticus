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
  Implements a projected orbital radus property extractor class.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusOrbitalProjected">
   <description>A property extractor that extracts the projected orbital radius. Projection is always along the $z$-axis.</description>
   <deepCopy>
    <functionClass variables="nodePropertyExtractor_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="nodePropertyExtractor_"/>
   </stateStorable>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusOrbitalProjected
     !!{
     A property extractor that extracts the projected orbital radius. Projection is always along the $z$-axis.
     !!}
     private
     type(nodePropertyExtractorPositionOrbital), pointer :: nodePropertyExtractor_  => null()
   contains
     final     ::                radiusOrbitalProjectedDestructor
     procedure :: extract     => radiusOrbitalProjectedExtract
     procedure :: name        => radiusOrbitalProjectedName
     procedure :: description => radiusOrbitalProjectedDescription
     procedure :: unitsInSI   => radiusOrbitalProjectedUnitsInSI
  end type nodePropertyExtractorRadiusOrbitalProjected

  interface nodePropertyExtractorRadiusOrbitalProjected
     !!{
     Constructors for the {\normalfont \ttfamily radiusOrbitalProjected} output analysis class.
     !!}
     module procedure radiusOrbitalProjectedConstructorParameters
     module procedure radiusOrbitalProjectedConstructorInternal
  end interface nodePropertyExtractorRadiusOrbitalProjected

contains

  function radiusOrbitalProjectedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily radiusOrbitalProjected} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorRadiusOrbitalProjected)                :: self
    type(inputParameters                            ), intent(inout) :: parameters

    self=nodePropertyExtractorRadiusOrbitalProjected()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function radiusOrbitalProjectedConstructorParameters

  function radiusOrbitalProjectedConstructorInternal(useLastIsolatedTime,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily radiusOrbitalProjected} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusOrbitalProjected) :: self
   
    !![
    <referenceConstruct isResult="yes" owner="self" object="nodePropertyExtractor_" constructor="nodePropertyExtractorPositionOrbital()"/>
    !!]
    return
  end function radiusOrbitalProjectedConstructorInternal

  subroutine radiusOrbitalProjectedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily mass} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusOrbitalProjected), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine radiusOrbitalProjectedDestructor

  double precision function radiusOrbitalProjectedExtract(self,node,instance)
    !!{
    Implement a radiusOrbitalProjected output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodePropertyExtractorRadiusOrbitalProjected), intent(inout), target   :: self
    type            (treeNode                                   ), intent(inout), target   :: node
    type            (multiCounter                               ), intent(inout), optional :: instance
    class           (nodeComponentBasic                         ), pointer                 :: basic
    double precision                                             , dimension(3)            :: positionOrbital
    !$GLC attributes unused :: instance

    basic                         => node                       %basic (                           )
    positionOrbital               =  self%nodePropertyExtractor_%extract(node,basic%time(),instance)
    radiusOrbitalProjectedExtract =  sqrt(sum(positionOrbital(1:2)**2))
    return
  end function radiusOrbitalProjectedExtract

  function radiusOrbitalProjectedName(self)
    !!{
    Return the name of the radiusOrbitalProjected property.
    !!}
    implicit none
    type (varying_string                             )                :: radiusOrbitalProjectedName
    class(nodePropertyExtractorRadiusOrbitalProjected), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusOrbitalProjectedName=var_str('radiusOrbitalProjected')
    return
  end function radiusOrbitalProjectedName

  function radiusOrbitalProjectedDescription(self)
    !!{
    Return a description of the radiusOrbitalProjected property.
    !!}
    implicit none
    type (varying_string                             )                :: radiusOrbitalProjectedDescription
    class(nodePropertyExtractorRadiusOrbitalProjected), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusOrbitalProjectedDescription=var_str('The projected (along the z-axis) orbital radius of the halo.')
    return
  end function radiusOrbitalProjectedDescription

  double precision function radiusOrbitalProjectedUnitsInSI(self)
    !!{
    Return the units of the radiusOrbitalProjected property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusOrbitalProjected), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusOrbitalProjectedUnitsInSI=megaParsec
    return
  end function radiusOrbitalProjectedUnitsInSI
