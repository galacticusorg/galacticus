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
Implements a node property extractor class for halo environment.
!!}

  use :: Cosmological_Density_Field, only : haloEnvironment, haloEnvironmentClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorHaloEnvironment" docformat="rst">
   <description>
   Extracts environmental metrics for dark matter halos, specifically the linear and non-linear local overdensity, characterizing the large-scale structure environment that influences halo formation rates, assembly bias, and galaxy evolution within the cosmological density field.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorHaloEnvironment
     !!{RST
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
     procedure :: units        => haloEnvironmentUnits
  end type nodePropertyExtractorHaloEnvironment

  interface nodePropertyExtractorHaloEnvironment
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorHaloEnvironment` property extractor class.
     !!}
     module procedure haloEnvironmentConstructorParameters
     module procedure haloEnvironmentConstructorInternal
  end interface nodePropertyExtractorHaloEnvironment

contains

  function haloEnvironmentConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorHaloEnvironment` property extractor class which takes a parameter set as input.
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
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorHaloEnvironment` property extractor class.
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
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorHaloEnvironment` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorHaloEnvironment), intent(inout) :: self

    !![
    <objectDestructor name="self%haloEnvironment_"/>
    !!]
    return
  end subroutine haloEnvironmentDestructor

  integer function haloEnvironmentElementCount(self,time)
    !!{RST
    Return the number of elements in the ``haloEnvironment`` property extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorHaloEnvironment), intent(inout) :: self
    double precision                                      , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    haloEnvironmentElementCount=2
    return
  end function haloEnvironmentElementCount

  function haloEnvironmentExtract(self,node,time,instance)
    !!{RST
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
    !!{RST
    Return the name of the ``haloEnvironment`` property.
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
    !!{RST
    Return a description of the ``haloEnvironment`` property.
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
    !!{RST
    Return the units of the ``haloEnvironment`` property in the SI system.
    !!}
    implicit none
    double precision                                      , allocatable  , dimension(:) :: haloEnvironmentUnitsInSI
    class           (nodePropertyExtractorHaloEnvironment), intent(inout)               :: self
    double precision                                      , intent(in   )               :: time
    !$GLC attributes unused :: self, time

    allocate(haloEnvironmentUnitsInSI(2))
    haloEnvironmentUnitsInSI=1.0d0
    return
  end function haloEnvironmentUnitsInSI

  function haloEnvironmentUnits(self,time) result(units)
    !!{RST
    Return the units of the haloEnvironment properties.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType                            ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorHaloEnvironment), intent(inout)             :: self
    double precision                                      , intent(in   )             :: time
    integer                                                                           :: i
    !$GLC attributes unused :: self, time

    allocate(units(2))
    do i=1,2
       units(i)=unitType(1.0d0)
    end do
    return
  end function haloEnvironmentUnits
