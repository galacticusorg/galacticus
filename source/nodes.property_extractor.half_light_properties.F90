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
Contains a module which implements a half-light radii property extractor class.
!!}

  use :: Galactic_Structure, only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiiHalfLightProperties">
   <description>
    A node property extractor which extracts half-light radii and the masses enclosed within them. The half-light radius in
    each specified luminosity band is extracted as {\normalfont \ttfamily [halfLightRadius\{luminosityID\}]} (in Mpc), where
    {\normalfont \ttfamily\{luminosityID\}} is the usual luminosity identifier suffix, and the total (dark + baryonic) mass
    within that radius is extracted as {\normalfont \ttfamily [halfLightMass\{luminosityID\}]} (in $M_\odot$).
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorRadiiHalfLightProperties
     !!{
     A half-light radii property extractor class.
     !!}
     private
     class(galacticStructureClass), pointer :: galacticStructure_ => null()
   contains
     final     ::                 radiiHalfLightPropertiesDestructor
     procedure :: elementCount => radiiHalfLightPropertiesElementCount
     procedure :: extract      => radiiHalfLightPropertiesExtract
     procedure :: names        => radiiHalfLightPropertiesNames
     procedure :: descriptions => radiiHalfLightPropertiesDescriptions
     procedure :: unitsInSI    => radiiHalfLightPropertiesUnitsInSI
     procedure :: type         => radiiHalfLightPropertiesType
  end type nodePropertyExtractorRadiiHalfLightProperties

  interface nodePropertyExtractorRadiiHalfLightProperties
     !!{
     Constructors for the ``radiiHalfLightProperties'' output analysis class.
     !!}
     module procedure radiiHalfLightPropertiesConstructorParameters
     module procedure radiiHalfLightPropertiesConstructorInternal
  end interface nodePropertyExtractorRadiiHalfLightProperties

contains

  function radiiHalfLightPropertiesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiiHalfLightProperties)                :: self
    type (inputParameters                              ), intent(inout) :: parameters
    class(galacticStructureClass                       ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiiHalfLightProperties(galacticStructure_)
     !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
   return
  end function radiiHalfLightPropertiesConstructorParameters

  function radiiHalfLightPropertiesConstructorInternal(galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiiHalfLightProperties)                        :: self
    class(galacticStructureClass                       ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*galacticStructure_"/>
    !!]

    return
  end function radiiHalfLightPropertiesConstructorInternal
  
  subroutine radiiHalfLightPropertiesDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiiHalfLightProperties), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine radiiHalfLightPropertiesDestructor

  integer function radiiHalfLightPropertiesElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    class           (nodePropertyExtractorRadiiHalfLightProperties), intent(inout) :: self
    double precision                                               , intent(in   ) :: time
    !$GLC attributes unused :: self

    radiiHalfLightPropertiesElementCount=2*unitStellarLuminosities%luminosityOutputCount(time)
    return
  end function radiiHalfLightPropertiesElementCount

  function radiiHalfLightPropertiesExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily radiiHalfLightProperties} property extractor.
    !!}
    use :: Galactic_Structure_Options    , only : componentTypeAll       , massTypeAll, massTypeStellar, weightByLuminosity
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    double precision                                               , dimension(:) , allocatable :: radiiHalfLightPropertiesExtract
    class           (nodePropertyExtractorRadiiHalfLightProperties), intent(inout), target      :: self
    type            (treeNode                                     ), intent(inout), target      :: node
    double precision                                               , intent(in   )              :: time
    type            (multiCounter                                 ), intent(inout), optional    :: instance
    integer                                                                                     :: i                              , j
    double precision                                                                            :: halfLightRadius                , massEnclosed
    !$GLC attributes unused :: self, instance

    allocate(radiiHalfLightPropertiesExtract(2*unitStellarLuminosities%luminosityOutputCount(time)))
    j=-1
    do i=1,unitStellarLuminosities%luminosityCount()
       if (unitStellarLuminosities%isOutput(i,time)) then
          halfLightRadius=self%galacticStructure_%radiusEnclosingMass(node,massFractional=0.5d0                                         ,massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=i)
          massEnclosed   =self%galacticStructure_%massEnclosed       (node,               halfLightRadius,componentType=componentTypeAll,massType=massTypeAll                                              )
          j=j+1
          radiiHalfLightPropertiesExtract(2*j+1:2*j+2)=[halfLightRadius,massEnclosed]
       end if
    end do
    return
  end function radiiHalfLightPropertiesExtract

  subroutine radiiHalfLightPropertiesNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily radiiHalfLightProperties} properties.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    class           (nodePropertyExtractorRadiiHalfLightProperties), intent(inout)                             :: self
    double precision                                               , intent(in   )                             :: time
    type            (varying_string                               ), intent(inout), dimension(:) , allocatable :: names
    integer                                                                                                    :: i                            , j
    !$GLC attributes unused :: self

    allocate(names(2*unitStellarLuminosities%luminosityOutputCount(time)))
    j=-1
    do i=1,unitStellarLuminosities%luminosityCount()
       if (unitStellarLuminosities%isOutput(i,time)) then
          j=j+1
          names(2*j+1)=var_str('halfLightRadius')//unitStellarLuminosities%name(i)
          names(2*j+2)=var_str('halfLightMass'  )//unitStellarLuminosities%name(i)
       end if
    end do
    return
  end subroutine radiiHalfLightPropertiesNames

  subroutine radiiHalfLightPropertiesDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    class           (nodePropertyExtractorRadiiHalfLightProperties), intent(inout)                             :: self
    double precision                                               , intent(in   )                             :: time
    type            (varying_string                               ), intent(inout), dimension(:) , allocatable :: descriptions
    integer                                                                                                    :: i
    !$GLC attributes unused :: self

    allocate(descriptions(2*unitStellarLuminosities%luminosityOutputCount(time)))
    do i=0,unitStellarLuminosities%luminosityOutputCount(time)-1
       descriptions(2*i+1)=var_str('Radius enclosing half the galaxy light [Mpc]'          )
       descriptions(2*i+2)=var_str('Mass enclosed within the galaxy half-light radius [M☉]')
    end do
    return
  end subroutine radiiHalfLightPropertiesDescriptions

  function radiiHalfLightPropertiesUnitsInSI(self,time)
    !!{
    Return the units of the last isolated redshift property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar              , megaParsec
    use :: Stellar_Luminosities_Structure  , only : unitStellarLuminosities
    implicit none
    double precision                                               , allocatable  , dimension(:) :: radiiHalfLightPropertiesUnitsInSI
    class           (nodePropertyExtractorRadiiHalfLightProperties), intent(inout)               :: self
    double precision                                               , intent(in   )               :: time
    integer                                                                                      :: i
    !$GLC attributes unused :: self

    allocate(radiiHalfLightPropertiesUnitsInSI(2*unitStellarLuminosities%luminosityOutputCount(time)))
    do i=0,unitStellarLuminosities%luminosityOutputCount(time)-1
       radiiHalfLightPropertiesUnitsInSI(2*i+1:2*i+2)=[megaParsec,massSolar]
    end do
    return
  end function radiiHalfLightPropertiesUnitsInSI

  integer function radiiHalfLightPropertiesType(self)
    !!{
    Return the type of the last isolated redshift property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorRadiiHalfLightProperties), intent(inout) :: self
    !$GLC attributes unused :: self

    radiiHalfLightPropertiesType=outputAnalysisPropertyTypeLinear
    return
  end function radiiHalfLightPropertiesType
