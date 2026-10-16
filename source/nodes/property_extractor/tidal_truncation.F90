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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !![
  <nodePropertyExtractor name="nodePropertyExtractorTidalTruncation" docformat="rst">
   <description>
   A node property extractor which extracts the degree of tidal truncation, :math:`x=-\log f`, where :math:`f` is the fraction of
   the mass of its untruncated mass distribution which the component retains, of the disk, spheroid, and nuclear star cluster
   components. These are accumulated by the tidal mass loss node operators (e.g.
   :galacticus-class:`nodeOperatorTidalMassLossSpheroids`). Zero is returned for any component which does not exist, or for
   which no tidal mass loss operator is active.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorTidalTruncation
     !!{RST
     A property extractor which extracts the degree of tidal truncation of disk, spheroid, and nuclear star cluster components.
     !!}
     private
     integer :: tidalTruncationDiskID, tidalTruncationSpheroidID, &
          &     tidalTruncationNSCID
   contains
     procedure :: elementCount => tidalTruncationElementCount
     procedure :: extract      => tidalTruncationExtract
     procedure :: names        => tidalTruncationNames
     procedure :: descriptions => tidalTruncationDescriptions
     procedure :: unitsInSI    => tidalTruncationUnitsInSI
     procedure :: units        => tidalTruncationUnits
  end type nodePropertyExtractorTidalTruncation

  interface nodePropertyExtractorTidalTruncation
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorTidalTruncation` property extractor class.
     !!}
     module procedure tidalTruncationConstructorParameters
     module procedure tidalTruncationConstructorInternal
  end interface nodePropertyExtractorTidalTruncation

contains

  function tidalTruncationConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorTidalTruncation` property extractor class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorTidalTruncation)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=nodePropertyExtractorTidalTruncation()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function tidalTruncationConstructorParameters

  function tidalTruncationConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorTidalTruncation` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorTidalTruncation) :: self

    !![
    <addMetaProperty component="disk"     name="tidalTruncation" id="self%tidalTruncationDiskID"     isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="spheroid" name="tidalTruncation" id="self%tidalTruncationSpheroidID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="NSC"      name="tidalTruncation" id="self%tidalTruncationNSCID"      isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function tidalTruncationConstructorInternal

  integer function tidalTruncationElementCount(self,time)
    !!{RST
    Return the number of elements in the ``tidalTruncation`` property extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorTidalTruncation), intent(inout) :: self
    double precision                                      , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    tidalTruncationElementCount=3
    return
  end function tidalTruncationElementCount

  function tidalTruncationExtract(self,node,time,instance)
    !!{RST
    Implement a ``tidalTruncation`` output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    implicit none
    double precision                                      , dimension(:) , allocatable :: tidalTruncationExtract
    class           (nodePropertyExtractorTidalTruncation), intent(inout), target      :: self
    type            (treeNode                            ), intent(inout), target      :: node
    double precision                                      , intent(in   )              :: time
    type            (multiCounter                        ), intent(inout), optional    :: instance
    class           (nodeComponentDisk                   )               , pointer     :: disk
    class           (nodeComponentSpheroid               )               , pointer     :: spheroid
    class           (nodeComponentNSC                    )               , pointer     :: nuclearStarCluster
    !$GLC attributes unused :: time, instance

    allocate(tidalTruncationExtract(3))
    tidalTruncationExtract=0.0d0
    disk               => node%disk    ()
    spheroid           => node%spheroid()
    nuclearStarCluster => node%NSC     ()
    select type (disk              )
    type is (nodeComponentDisk    )
       ! Disk does not yet exist.
    class default
       if (disk              %floatRank0MetaPropertyIsCreated(self%tidalTruncationDiskID    )) &
            & tidalTruncationExtract(1)=disk              %floatRank0MetaPropertyGet(self%tidalTruncationDiskID    )
    end select
    select type (spheroid          )
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist.
    class default
       if (spheroid          %floatRank0MetaPropertyIsCreated(self%tidalTruncationSpheroidID)) &
            & tidalTruncationExtract(2)=spheroid          %floatRank0MetaPropertyGet(self%tidalTruncationSpheroidID)
    end select
    select type (nuclearStarCluster)
    type is (nodeComponentNSC     )
       ! Nuclear star cluster does not yet exist.
    class default
       if (nuclearStarCluster%floatRank0MetaPropertyIsCreated(self%tidalTruncationNSCID     )) &
            & tidalTruncationExtract(3)=nuclearStarCluster%floatRank0MetaPropertyGet(self%tidalTruncationNSCID     )
    end select
    return
  end function tidalTruncationExtract

  subroutine tidalTruncationNames(self,time,names)
    !!{RST
    Return the names of the ``tidalTruncation`` properties.
    !!}
    implicit none
    class           (nodePropertyExtractorTidalTruncation), intent(inout)                             :: self
    double precision                                      , intent(in   )                             :: time
    type            (varying_string                      ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(3))
    names(1)=var_str('diskTidalTruncation'              )
    names(2)=var_str('spheroidTidalTruncation'          )
    names(3)=var_str('nuclearStarClusterTidalTruncation')
    return
  end subroutine tidalTruncationNames

  subroutine tidalTruncationDescriptions(self,time,descriptions)
    !!{RST
    Return the descriptions of the ``tidalTruncation`` properties.
    !!}
    implicit none
    class           (nodePropertyExtractorTidalTruncation), intent(inout)                             :: self
    double precision                                      , intent(in   )                             :: time
    type            (varying_string                      ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)=var_str('Degree of tidal truncation, -log(f), of the disk, where f is the fraction of its untruncated mass retained.'                )
    descriptions(2)=var_str('Degree of tidal truncation, -log(f), of the spheroid, where f is the fraction of its untruncated mass retained.'            )
    descriptions(3)=var_str('Degree of tidal truncation, -log(f), of the nuclear star cluster, where f is the fraction of its untruncated mass retained.')
    return
  end subroutine tidalTruncationDescriptions

  function tidalTruncationUnitsInSI(self,time)
    !!{RST
    Return the units of the ``tidalTruncation`` properties in the SI system.
    !!}
    implicit none
    double precision                                      , dimension(:) , allocatable :: tidalTruncationUnitsInSI
    class           (nodePropertyExtractorTidalTruncation), intent(inout)              :: self
    double precision                                      , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(tidalTruncationUnitsInSI(3))
    tidalTruncationUnitsInSI=1.0d0
    return
  end function tidalTruncationUnitsInSI

  function tidalTruncationUnits(self,time) result(units)
    !!{RST
    Return the units of the ``tidalTruncation`` properties.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType                            ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorTidalTruncation), intent(inout)             :: self
    double precision                                      , intent(in   )             :: time
    integer                                                                           :: i
    !$GLC attributes unused :: self, time

    allocate(units(3))
    do i=1,3
       units(i)=unitType(1.0d0)
    end do
    return
  end function tidalTruncationUnits
