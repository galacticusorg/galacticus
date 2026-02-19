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
Implements a node property extractor class for halo environment.
!!}

  use :: Cosmological_Density_Field, only : haloEnvironment, haloEnvironmentClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorHaloEnvironment">
   <description>A node property extractor class for halo environment.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorHaloEnvironment
     !!{
     A property extractor class for halo environment.
     !!}
     private
     class(haloEnvironmentClass), pointer :: haloEnvironment_ => null()
   contains
     final     ::                 haloEnvironmentDestructor
     procedure :: elementCount => haloEnvironmentElementCount
     procedure :: extract      => haloEnvironmentExtract
     procedure :: names        => haloEnvironmentNames
     procedure :: descriptions => haloEnvironmentDescriptions
     procedure :: unitsInSI    => haloEnvironmentUnitsInSI
  end type nodePropertyExtractorHaloEnvironment

  interface nodePropertyExtractorHaloEnvironment
     !!{
     Constructors for the \refClass{nodePropertyExtractorHaloEnvironment} output analysis class.
     !!}
     module procedure haloEnvironmentConstructorParameters
     module procedure haloEnvironmentConstructorInternal
  end interface nodePropertyExtractorHaloEnvironment

contains

  function haloEnvironmentConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorHaloEnvironment} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorHaloEnvironment)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(haloEnvironmentClass                ), pointer       :: haloEnvironment_

    !![
    <objectBuilder class="haloEnvironment" name="haloEnvironment_" source="parameters"/>
    !!]
    self=nodePropertyExtractorHaloEnvironment(haloEnvironment_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloEnvironment_"/>
    !!]
    return
  end function haloEnvironmentConstructorParameters

  function haloEnvironmentConstructorInternal(haloEnvironment_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorHaloEnvironment} output analysis property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorHaloEnvironment)                        :: self
    class(haloEnvironmentClass                ), intent(in   ), target :: haloEnvironment_
    !![
    <constructorAssign variables="*haloEnvironment_"/>
    !!]

    return
  end function haloEnvironmentConstructorInternal

  subroutine haloEnvironmentDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorHaloEnvironment} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorHaloEnvironment), intent(inout) :: self

    !![
    <objectDestructor name="self%haloEnvironment_"/>
    !!]
    return
  end subroutine haloEnvironmentDestructor

  integer function haloEnvironmentElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily haloEnvironment} property extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorHaloEnvironment), intent(inout) :: self
    double precision                                      , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    haloEnvironmentElementCount=2
    return
  end function haloEnvironmentElementCount

  function haloEnvironmentExtract(self,node,time,instance)
    !!{
    Implement extraction of halo environment properties.
    !!}
    implicit none
    double precision                                      , dimension(:) , allocatable :: haloEnvironmentExtract
    class           (nodePropertyExtractorHaloEnvironment), intent(inout), target      :: self
    type            (treeNode                            ), intent(inout), target      :: node
    double precision                                      , intent(in   )              :: time
    type            (multiCounter                        ), intent(inout), optional    :: instance
    !$GLC attributes unused :: time, instance

    allocate(haloEnvironmentExtract(2))
    haloEnvironmentExtract=[                                                  &
         &                  self%haloEnvironment_%overdensityLinear   (node), &
         &                  self%haloEnvironment_%overdensityNonLinear(node)  &
         &                 ]
    return
  end function haloEnvironmentExtract

  subroutine haloEnvironmentNames(self,time,names)
    !!{
    Return the name of the {\normalfont \ttfamily haloEnvironment} property.
    !!}
    implicit none
    class           (nodePropertyExtractorHaloEnvironment), intent(inout)                             :: self
    double precision                                      , intent(in   )                             :: time
    type            (varying_string                      ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(2))
    names(1)=var_str('haloEnvironmentOverdensityLinear'   )
    names(2)=var_str('haloEnvironmentOverdensityNonLinear')
    return
  end subroutine haloEnvironmentNames

  subroutine haloEnvironmentDescriptions(self,time,descriptions)
    !!{
    Return a description of the {\normalfont \ttfamily haloEnvironment} property.
    !!}
    implicit none
    class           (nodePropertyExtractorHaloEnvironment), intent(inout)                             :: self
    double precision                                      , intent(in   )                             :: time
    type            (varying_string                      ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(2))
    descriptions(1)=var_str('Environmental linear overdensity of the halo [].'    )
    descriptions(2)=var_str('Environmental non-linear overdensity of the halo [].')
    return
  end subroutine haloEnvironmentDescriptions

  function haloEnvironmentUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily haloEnvironment} property in the SI system.
    !!}
    implicit none
    double precision                                      , allocatable  , dimension(:) :: haloEnvironmentUnitsInSI
    class           (nodePropertyExtractorHaloEnvironment), intent(inout)               :: self
    double precision                                      , intent(in   )               :: time
    !$GLC attributes unused :: self, time

    allocate(haloEnvironmentUnitsInSI(2))
    haloEnvironmentUnitsInSI=0.0d0
    return
  end function haloEnvironmentUnitsInSI

