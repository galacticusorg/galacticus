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
Implements a dark matter profile SIDM interaction radius property extractor class.
!!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO, darkMatterProfileDMOClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM">
   <description>A dark matter profile SIDM interaction radius property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM
     !!{
     A dark matter profile SIDM interaction radius property extractor class.
     !!}
     private
     class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
   contains
     final     ::                darkMatterProfileRadiusInteractionSIDMDestructor
     procedure :: extract     => darkMatterProfileRadiusInteractionSIDMExtract
     procedure :: name        => darkMatterProfileRadiusInteractionSIDMName
     procedure :: description => darkMatterProfileRadiusInteractionSIDMDescription
     procedure :: unitsInSI   => darkMatterProfileRadiusInteractionSIDMUnitsInSI
  end type nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM

  interface nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM
     !!{
     Constructors for the {\normalfont \ttfamily darkMatterProfileRadiusInteractionSIDM} output analysis class.
     !!}
     module procedure darkMatterProfileRadiusInteractionSIDMConstructorParameters
     module procedure darkMatterProfileRadiusInteractionSIDMConstructorInternal
  end interface nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM

contains

  function darkMatterProfileRadiusInteractionSIDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily darkMatterProfileRadiusInteractionSIDM} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM)                :: self
    type (inputParameters                                            ), intent(inout) :: parameters
    class(darkMatterProfileDMOClass                                  ), pointer       :: darkMatterProfileDMO_

    !![
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM(darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkMatterProfileRadiusInteractionSIDMConstructorParameters

  function darkMatterProfileRadiusInteractionSIDMConstructorInternal(darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily darkMatterProfileRadiusInteractionSIDM} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM)                        :: self
    class(darkMatterProfileDMOClass                                  ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="*darkMatterProfileDMO_"/>
    !!]

    return
  end function darkMatterProfileRadiusInteractionSIDMConstructorInternal

  subroutine darkMatterProfileRadiusInteractionSIDMDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily darkMatterProfileRadiusInteractionSIDM} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine darkMatterProfileRadiusInteractionSIDMDestructor

  double precision function darkMatterProfileRadiusInteractionSIDMExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily darkMatterProfileRadiusInteractionSIDM} output analysis.
    !!}
    use :: Mass_Distributions, only : massDistributionSphericalSIDM, massDistributionClass
    implicit none
    class(nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM), intent(inout), target   :: self
    type (treeNode                                                   ), intent(inout), target   :: node
    type (multiCounter                                               ), intent(inout), optional :: instance
    class(massDistributionClass                                      ), pointer                 :: massDistribution_
    !$GLC attributes unused :: instance

    massDistribution_ => self%darkMatterProfileDMO_%get(node)
    select type (massDistribution_)
    class is (massDistributionSphericalSIDM)
       darkMatterProfileRadiusInteractionSIDMExtract=massDistribution_%radiusInteraction()
    class default
       darkMatterProfileRadiusInteractionSIDMExtract=0.0d0
    end select
    return
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
  end function darkMatterProfileRadiusInteractionSIDMExtract

  function darkMatterProfileRadiusInteractionSIDMName(self)
    !!{
    Return the name of the {\normalfont \ttfamily darkMatterProfileRadiusInteractionSIDM} property.
    !!}
    implicit none
    type (varying_string                                             )                :: darkMatterProfileRadiusInteractionSIDMName
    class(nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileRadiusInteractionSIDMName=var_str('darkMatterProfileRadiusInteractionSIDM')
    return
  end function darkMatterProfileRadiusInteractionSIDMName

  function darkMatterProfileRadiusInteractionSIDMDescription(self)
    !!{
    Return a description of the {\normalfont \ttfamily darkMatterProfileRadiusInteractionSIDM} property.
    !!}
    implicit none
    type (varying_string                                             )                :: darkMatterProfileRadiusInteractionSIDMDescription
    class(nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileRadiusInteractionSIDMDescription=var_str('The SIDM interaction radius, r₁, in the dark matter profile.')
    return
  end function darkMatterProfileRadiusInteractionSIDMDescription

  double precision function darkMatterProfileRadiusInteractionSIDMUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily darkMatterProfileRadiusInteractionSIDM} property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorDarkMatterProfileRadiusInteractionSIDM), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileRadiusInteractionSIDMUnitsInSI=megaParsec
    return
  end function darkMatterProfileRadiusInteractionSIDMUnitsInSI
