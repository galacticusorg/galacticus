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
  Implements a property extractor class which extracts the radius at which the maximum velocity is reached in the dark
  matter-only density profile.
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO, darkMatterProfileDMOClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusVelocityMaximum">
   <description>A property extractor class which extracts the radius at which the maximum velocity is reached in the dark matter-only density profile.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusVelocityMaximum
     !!{
     A property extractor class which extracts the radius at which the maximum velocity is reached in the dark matter-only
     density profile.
     !!}
     private
     class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     ! Name of the property.
     type (varying_string           )          :: propertyName
   contains
     final     ::                radiusVelocityMaximumDestructor
     procedure :: extract     => radiusVelocityMaximumExtract
     procedure :: name        => radiusVelocityMaximumName
     procedure :: description => radiusVelocityMaximumDescription
     procedure :: unitsInSI   => radiusVelocityMaximumUnitsInSI
  end type nodePropertyExtractorRadiusVelocityMaximum

  interface nodePropertyExtractorRadiusVelocityMaximum
     !!{
     Constructors for the {\normalfont \ttfamily radiusVelocityMaximum} output analysis class.
     !!}
     module procedure radiusVelocityMaximumConstructorParameters
     module procedure radiusVelocityMaximumConstructorInternal
  end interface nodePropertyExtractorRadiusVelocityMaximum

contains

  function radiusVelocityMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily radiusVelocityMaximum} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRadiusVelocityMaximum)                :: self
    type (inputParameters                           ), intent(inout) :: parameters
    class(darkMatterProfileDMOClass                 ), pointer       :: darkMatterProfileDMO_
    type (varying_string                            )                :: propertyName

    !![
    <inputParameter>
      <name>propertyName</name>
      <defaultValue>var_str('RadiusVelocityMaximum')</defaultValue>
      <source>parameters</source>
      <description>Name of the property.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiusVelocityMaximum(propertyName,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function radiusVelocityMaximumConstructorParameters

  function radiusVelocityMaximumConstructorInternal(propertyName,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily radiusVelocityMaximum} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRadiusVelocityMaximum)                        :: self
    class(darkMatterProfileDMOClass                 ), intent(in   ), target :: darkMatterProfileDMO_
    type (varying_string                            ), intent(in   )         :: propertyName
    !![
    <constructorAssign variables="propertyName,*darkMatterProfileDMO_"/>
    !!]

    return
  end function radiusVelocityMaximumConstructorInternal

  subroutine radiusVelocityMaximumDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily radiusVelocityMaximum} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusVelocityMaximum), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine radiusVelocityMaximumDestructor

  double precision function radiusVelocityMaximumExtract(self,node,instance)
    !!{
    Implement a radius of maximum velocity output analysis.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic   , treeNode
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class(nodePropertyExtractorRadiusVelocityMaximum), intent(inout), target   :: self
    type (treeNode                                  ), intent(inout), target   :: node
    type (multiCounter                              ), intent(inout), optional :: instance
    class(massDistributionClass                     )               , pointer  :: massDistribution_
    !$GLC attributes unused :: instance

    massDistribution_            => self             %darkMatterProfileDMO_%get                       (node)
    radiusVelocityMaximumExtract =  massDistribution_                      %radiusRotationCurveMaximum(    )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function radiusVelocityMaximumExtract

  function radiusVelocityMaximumName(self)
    !!{
    Return the name of the radius of maximum velocity property.
    !!}
    implicit none
    type (varying_string                            )                :: radiusVelocityMaximumName
    class(nodePropertyExtractorRadiusVelocityMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusVelocityMaximumName=var_str('darkMatterProfileDMO'//char(self%propertyName))
    return
  end function radiusVelocityMaximumName

  function radiusVelocityMaximumDescription(self)
    !!{
    Return a description of the radius of maximum velocity property.
    !!}
    implicit none
    type (varying_string                            )                :: radiusVelocityMaximumDescription
    class(nodePropertyExtractorRadiusVelocityMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusVelocityMaximumDescription=var_str('Radius of the maximum rotation velocity of the dark matter profile [Mpc].')
    return
  end function radiusVelocityMaximumDescription

  double precision function radiusVelocityMaximumUnitsInSI(self)
    !!{
    Return the units of the radius of maximum velocity property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusVelocityMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusVelocityMaximumUnitsInSI=megaParsec
    return
  end function radiusVelocityMaximumUnitsInSI

