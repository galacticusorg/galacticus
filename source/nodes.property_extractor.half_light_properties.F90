!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements a half-light radii property extractor class.

  !# <nodePropertyExtractor name="nodePropertyExtractorRadiiHalfLightProperties">
  !#  <description>A half-light radii property extractor class.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorRadiiHalfLightProperties
     !% A half-light radii property extractor class.
     private
   contains
     procedure :: elementCount => radiiHalfLightPropertiesElementCount
     procedure :: extract      => radiiHalfLightPropertiesExtract
     procedure :: names        => radiiHalfLightPropertiesNames
     procedure :: descriptions => radiiHalfLightPropertiesDescriptions
     procedure :: unitsInSI    => radiiHalfLightPropertiesUnitsInSI
     procedure :: type         => radiiHalfLightPropertiesType
  end type nodePropertyExtractorRadiiHalfLightProperties

  interface nodePropertyExtractorRadiiHalfLightProperties
     !% Constructors for the ``radiiHalfLightProperties'' output analysis class.
     module procedure radiiHalfLightPropertiesConstructorParameters
  end interface nodePropertyExtractorRadiiHalfLightProperties

contains

  function radiiHalfLightPropertiesConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiiHalfLightProperties)                :: self
    type (inputParameters                              ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=nodePropertyExtractorRadiiHalfLightProperties()
    return
  end function radiiHalfLightPropertiesConstructorParameters

  integer function radiiHalfLightPropertiesElementCount(self,time)
    !% Return the number of elements in the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    class           (nodePropertyExtractorRadiiHalfLightProperties), intent(inout) :: self
    double precision                                               , intent(in   ) :: time
    !$GLC attributes unused :: self

    radiiHalfLightPropertiesElementCount=2*unitStellarLuminosities%luminosityOutputCount(time)
    return
  end function radiiHalfLightPropertiesElementCount

  function radiiHalfLightPropertiesExtract(self,node,time,instance)
    !% Implement a {\normalfont \ttfamily radiiHalfLightProperties} property extractor.
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass, Galactic_Structure_Radius_Enclosing_Mass
    use :: Galactic_Structure_Options        , only : componentTypeAll                , massTypeAll                             , massTypeStellar, weightByLuminosity
    use :: Stellar_Luminosities_Structure    , only : unitStellarLuminosities
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
          halfLightRadius=Galactic_Structure_Radius_Enclosing_Mass(node,fractionalMass=0.5d0                                         ,massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=i)
          massEnclosed   =Galactic_Structure_Enclosed_Mass        (node,               halfLightRadius,componentType=componentTypeAll,massType=massTypeAll                                              )
          j=j+1
          radiiHalfLightPropertiesExtract(2*j+1:2*j+2)=[halfLightRadius,massEnclosed]
       end if
    end do
    return
  end function radiiHalfLightPropertiesExtract

  function radiiHalfLightPropertiesNames(self,time)
    !% Return the names of the {\normalfont \ttfamily radiiHalfLightProperties} properties.
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    type            (varying_string                               ), dimension(:) , allocatable :: radiiHalfLightPropertiesNames
    class           (nodePropertyExtractorRadiiHalfLightProperties), intent(inout)              :: self
    double precision                                               , intent(in   )              :: time
    integer                                                                                     :: i
    !$GLC attributes unused :: self

    allocate(radiiHalfLightPropertiesNames(2*unitStellarLuminosities%luminosityOutputCount(time)))
    do i=0,unitStellarLuminosities%luminosityOutputCount(time)-1
       radiiHalfLightPropertiesNames(2*i+1:2*i+2)=[                                                               &
            &                                      var_str('halfLightRadius')//unitStellarLuminosities%name(i+1), &
            &                                      var_str('halfLightMass'  )//unitStellarLuminosities%name(i+1)  &
            &                                     ]
    end do
    return
  end function radiiHalfLightPropertiesNames

  function radiiHalfLightPropertiesDescriptions(self,time)
    !% Return descriptions of the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    type            (varying_string                               ), dimension(:) , allocatable :: radiiHalfLightPropertiesDescriptions
    class           (nodePropertyExtractorRadiiHalfLightProperties), intent(inout)              :: self
    double precision                                               , intent(in   )              :: time
    integer                                                                                     :: i
    !$GLC attributes unused :: self

    allocate(radiiHalfLightPropertiesDescriptions(2*unitStellarLuminosities%luminosityOutputCount(time)))
    do i=0,unitStellarLuminosities%luminosityOutputCount(time)-1
       radiiHalfLightPropertiesDescriptions(2*i+1:2*i+2)=[                                                                   &
            &                                             var_str('Radius enclosing half the galaxy light [Mpc]'          ), &
            &                                             var_str('Mass enclosed within the galaxy half-light radius [M☉]')  &
            &                                            ]
    end do
    return
  end function radiiHalfLightPropertiesDescriptions

  function radiiHalfLightPropertiesUnitsInSI(self,time)
    !% Return the units of the last isolated redshift property in the SI system.
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
    !% Return the type of the last isolated redshift property.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorRadiiHalfLightProperties), intent(inout) :: self
    !$GLC attributes unused :: self

    radiiHalfLightPropertiesType=outputAnalysisPropertyTypeLinear
    return
  end function radiiHalfLightPropertiesType
