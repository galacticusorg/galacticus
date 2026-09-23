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
Implements an orbital position output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorPositionOrbital" docformat="rst">
   <description>
   An orbital position output analysis property extractor class. Specifically, the orbital position is defined relative to the top-level halo in any sub-halo hierarchy. That is, relative to the host halo which is itself not a sub-halo of any other halo. If the position of a (sub)\ :math:`^i`-halo with respect to the center of its (sub)\ :math:`^{i-1}`-halo host is :math:`\mathbf{x}_i` then the orbital position computed by this class is

   .. math::

      \mathbf{x} = \sum_{i=1}^N \mathbf{x}_i,

   where :math:`N` is the depth of the node in the sub-halo hierarchy.

   If that top-level halo is not on the main branch of its tree---as when the orbits of halos prior to infall are tracked (see the ``trackPreInfallOrbit`` parameter of :galacticus-class:`nodeOperatorSatelliteOrbit`)---the sum continues from the halo, at the same time, on the branch with which the top-level halo will eventually merge, repeating until a halo on the main branch is reached. The orbital position is therefore always relative to a halo on the main branch.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorPositionOrbital
     !!{RST
     An orbital position property extractor output analysis class.
     !!}
     private
   contains
     procedure :: elementCount => positionOrbitalElementCount
     procedure :: extract      => positionOrbitalExtract
     procedure :: names        => positionOrbitalNames
     procedure :: descriptions => positionOrbitalDescriptions
     procedure :: unitsInSI    => positionOrbitalUnitsInSI
     procedure :: units       => positionOrbitalUnits
  end type nodePropertyExtractorPositionOrbital

  interface nodePropertyExtractorPositionOrbital
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorPositionOrbital` property extractor class.
     !!}
     module procedure positionOrbitalConstructorParameters
  end interface nodePropertyExtractorPositionOrbital

contains

  function positionOrbitalConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorPositionOrbital` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorPositionOrbital)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=nodePropertyExtractorPositionOrbital()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function positionOrbitalConstructorParameters

  integer function positionOrbitalElementCount(self,time)
    !!{RST
    Return the number of elements in the ``positionOrbital`` property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorPositionOrbital), intent(inout) :: self
    double precision                                      , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    positionOrbitalElementCount=3
    return
  end function positionOrbitalElementCount

  function positionOrbitalExtract(self,node,time,instance) result(position)
    !!{RST
    Implement a positionOrbital output analysis.
    !!}
    use :: Node_Orbital_Offsets, only : Node_Orbital_Offset
    implicit none
    double precision                                      , dimension(:) , allocatable :: position
    class           (nodePropertyExtractorPositionOrbital), intent(inout), target      :: self
    type            (treeNode                            ), intent(inout), target      :: node
    double precision                                      , intent(in   )              :: time
    type            (multiCounter                        ), intent(inout), optional    :: instance
    !$GLC attributes unused :: self, instance

    allocate(position(3))
    call Node_Orbital_Offset(node,time,position=position)
    return
  end function positionOrbitalExtract

  subroutine positionOrbitalNames(self,time,names)
    !!{RST
    Return the name of the positionOrbital property.
    !!}
    implicit none
    class(nodePropertyExtractorPositionOrbital), intent(inout)                             :: self
    double precision                           , intent(in   )                             :: time
    type (varying_string                      ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(3))
    names(1)=var_str('positionOrbitalX')
    names(2)=var_str('positionOrbitalY')
    names(3)=var_str('positionOrbitalZ')
    return
  end subroutine positionOrbitalNames

  subroutine positionOrbitalDescriptions(self,time,descriptions)
    !!{RST
    Return a description of the positionOrbital property.
    !!}
    implicit none
    class(nodePropertyExtractorPositionOrbital), intent(inout)                            :: self
    double precision                           , intent(in   )                            :: time
    type (varying_string                      ), intent(inout), dimension(:), allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)=var_str('The orbital x-position of the halo relative to the top-level host halo (i.e. the host which is not a sub-halo of any other halo).')
    descriptions(2)=var_str('The orbital y-position of the halo relative to the top-level host halo (i.e. the host which is not a sub-halo of any other halo).')
    descriptions(3)=var_str('The orbital z-position of the halo relative to the top-level host halo (i.e. the host which is not a sub-halo of any other halo).')
    return
  end subroutine positionOrbitalDescriptions

  function positionOrbitalUnitsInSI(self,time)
    !!{RST
    Return the units of the positionOrbital property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    double precision                                      , dimension(:) , allocatable :: positionOrbitalUnitsInSI
    class           (nodePropertyExtractorPositionOrbital), intent(inout)              :: self
    double precision                                      , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(positionOrbitalUnitsInSI(3))
    positionOrbitalUnitsInSI=megaParsec
    return
  end function positionOrbitalUnitsInSI

  function positionOrbitalUnits(self,time) result(units)
    !!{RST
    Return the units of the orbital position properties.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType                            ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorPositionOrbital), intent(inout)             :: self
    double precision                                      , intent(in   )             :: time
    double precision                                      , dimension(:), allocatable :: siValues
    integer                                                                           :: i

    siValues=self%unitsInSI(time)
    allocate(units(size(siValues)))
    do i=1,size(siValues)
       units(i)=unitType(siValues(i),description='Mpc',quantity='Mpc')
    end do
    return
  end function positionOrbitalUnits
