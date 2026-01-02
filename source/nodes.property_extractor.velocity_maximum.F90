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

!!{
Implements a cooling rate property extractor class.
!!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO, darkMatterProfileDMOClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorVelocityMaximum">
   <description>A cooling rate property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorVelocityMaximum
     !!{
     A velocityMaximum property extractor class.
     !!}
     private
     class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     ! Name of the property.
     type (varying_string           )          :: propertyName
   contains
     final     ::                velocityMaximumDestructor
     procedure :: extract     => velocityMaximumExtract
     procedure :: name        => velocityMaximumName
     procedure :: description => velocityMaximumDescription
     procedure :: unitsInSI   => velocityMaximumUnitsInSI
  end type nodePropertyExtractorVelocityMaximum

  interface nodePropertyExtractorVelocityMaximum
     !!{
     Constructors for the \refClass{nodePropertyExtractorVelocityMaximum} output analysis class.
     !!}
     module procedure velocityMaximumConstructorParameters
     module procedure velocityMaximumConstructorInternal
  end interface nodePropertyExtractorVelocityMaximum

contains

  function velocityMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorVelocityMaximum} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorVelocityMaximum)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(darkMatterProfileDMOClass           ), pointer       :: darkMatterProfileDMO_
    type (varying_string                      )                :: propertyName

    !![
    <inputParameter>
      <name>propertyName</name>
      <defaultValue>var_str('VelocityMaximum')</defaultValue>
      <source>parameters</source>
      <description>Name of the property.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=nodePropertyExtractorVelocityMaximum(propertyName,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function velocityMaximumConstructorParameters

  function velocityMaximumConstructorInternal(propertyName,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorVelocityMaximum} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorVelocityMaximum)                        :: self
    class(darkMatterProfileDMOClass           ), intent(in   ), target :: darkMatterProfileDMO_
    type (varying_string                      ), intent(in   )         :: propertyName
    !![
    <constructorAssign variables="propertyName,*darkMatterProfileDMO_"/>
    !!]

    return
  end function velocityMaximumConstructorInternal

  subroutine velocityMaximumDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorVelocityMaximum} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorVelocityMaximum), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine velocityMaximumDestructor

  double precision function velocityMaximumExtract(self,node,instance)
    !!{
    Implement a last isolated redshift output analysis.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class(nodePropertyExtractorVelocityMaximum), intent(inout), target   :: self
    type (treeNode                            ), intent(inout), target   :: node
    type (multiCounter                        ), intent(inout), optional :: instance
    class(massDistributionClass               )               , pointer  :: massDistribution_
    !$GLC attributes unused :: instance

    massDistribution_      => self             %darkMatterProfileDMO_%get                         (node)
    velocityMaximumExtract =  massDistribution_                      %velocityRotationCurveMaximum(    )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
   return
  end function velocityMaximumExtract

  function velocityMaximumName(self)
    !!{
    Return the name of the last isolated redshift property.
    !!}
    implicit none
    type (varying_string                      )                :: velocityMaximumName
    class(nodePropertyExtractorVelocityMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    velocityMaximumName=var_str('darkMatterProfileDMO'//char(self%propertyName))
    return
  end function velocityMaximumName

  function velocityMaximumDescription(self)
    !!{
    Return a description of the velocityMaximum property.
    !!}
    implicit none
    type (varying_string                      )                :: velocityMaximumDescription
    class(nodePropertyExtractorVelocityMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    velocityMaximumDescription=var_str('Maximum rotation velocity of the dark matter profile [km/s].')
    return
  end function velocityMaximumDescription

  double precision function velocityMaximumUnitsInSI(self)
    !!{
    Return the units of the last isolated redshift property in the SI system.
    !!}
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class(nodePropertyExtractorVelocityMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    velocityMaximumUnitsInSI=kilo
    return
  end function velocityMaximumUnitsInSI


