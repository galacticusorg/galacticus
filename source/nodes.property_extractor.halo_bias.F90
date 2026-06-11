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
Implements an ISM mass output analysis property extractor class.
!!}

  use :: Dark_Matter_Halo_Biases, only : darkMatterHaloBias, darkMatterHaloBiasClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorHaloBias" docformat="rst">
   <description>
   A node property extractor which extracts the large scale, linearly theory bias for each node. For satellite nodes, this corresponds to the bias of their host halo.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorHaloBias
     !!{RST
     A node property extractor for halo bias.
     !!}
     private
     class(darkMatterHaloBiasClass), pointer :: darkMatterHaloBias_ => null()
   contains
     final     ::                haloBiasDestructor
     procedure :: extract     => haloBiasExtract
     procedure :: name        => haloBiasName
     procedure :: description => haloBiasDescription
     procedure :: unitsInSI   => haloBiasUnitsInSI
     procedure :: units       => haloBiasUnits
  end type nodePropertyExtractorHaloBias

  interface nodePropertyExtractorHaloBias
     !!{RST
     Constructors for the ``nodePropertyExtractorHaloBias`` property extractor class.
     !!}
     module procedure haloBiasConstructorParameters
     module procedure haloBiasConstructorInternal
  end interface nodePropertyExtractorHaloBias

contains

  function haloBiasConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``nodePropertyExtractorHaloBias`` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorHaloBias)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(darkMatterHaloBiasClass      ), pointer       :: darkMatterHaloBias_

    !![
    <objectBuilder class="darkMatterHaloBias" name="darkMatterHaloBias_" source="parameters"/>
    !!]
    self=nodePropertyExtractorHaloBias(darkMatterHaloBias_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloBias_"/>
    !!]
    return
  end function haloBiasConstructorParameters

  function haloBiasConstructorInternal(darkMatterHaloBias_) result(self)
    !!{RST
    Internal constructor for the ``nodePropertyExtractorHaloBias`` property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorHaloBias)                        :: self
    class(darkMatterHaloBiasClass      ), intent(in   ), target :: darkMatterHaloBias_
    !![
    <constructorAssign variables="*darkMatterHaloBias_"/>
    !!]

    return
  end function haloBiasConstructorInternal

  subroutine haloBiasDestructor(self)
    !!{RST
    Destructor for the ``nodePropertyExtractorHaloBias`` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorHaloBias), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloBias_"/>
    !!]
    return
  end subroutine haloBiasDestructor

  double precision function haloBiasExtract(self,node,instance)
    !!{RST
    Implement a ``haloBias`` node property extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorHaloBias), intent(inout), target   :: self
    type            (treeNode                     ), intent(inout), target   :: node
    type            (multiCounter                 ), intent(inout), optional :: instance
    type            (treeNode                     )               , pointer  :: nodeIsolated
    !$GLC attributes unused :: instance

    nodeIsolated => node
    do while (nodeIsolated%isSatellite())
       nodeIsolated => nodeIsolated%parent
    end do
    haloBiasExtract=self%darkMatterHaloBias_%bias(nodeIsolated)
    return
  end function haloBiasExtract


  function haloBiasName(self)
    !!{RST
    Return the name of the haloBias property.
    !!}
    implicit none
    type (varying_string               )                :: haloBiasName
    class(nodePropertyExtractorHaloBias), intent(inout) :: self
    !$GLC attributes unused :: self

    haloBiasName=var_str('nodeBias')
    return
  end function haloBiasName

  function haloBiasDescription(self)
    !!{RST
    Return a description of the haloBias property.
    !!}
    implicit none
    type (varying_string               )                :: haloBiasDescription
    class(nodePropertyExtractorHaloBias), intent(inout) :: self
    !$GLC attributes unused :: self

    haloBiasDescription=var_str('The linear bias for this node.')
    return
  end function haloBiasDescription

  double precision function haloBiasUnitsInSI(self)
    !!{RST
    Return the units of the ``haloBias`` property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorHaloBias), intent(inout) :: self
    !$GLC attributes unused :: self

    haloBiasUnitsInSI=1.0d0
    return
  end function haloBiasUnitsInSI

  function haloBiasUnits(self) result(units)
    !!{RST
    Return the units of the haloBias property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                     )                :: units
    class(nodePropertyExtractorHaloBias), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(1.0d0)
    return
  end function haloBiasUnits
