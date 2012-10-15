!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements an object to store merger tree data for processing into \glc's preferred file format.

module Merger_Tree_Data_Structure
  !% Implements an object to store merger tree data for processing into \glc's preferred file format.
  use Kind_Numbers
  use ISO_Varying_String
  implicit none
  private
  public :: mergerTreeData

  ! Property labels.
  integer, parameter         :: propertyTypeNull                  = 1
  integer, parameter, public :: propertyTypeTreeIndex             = 2
  integer, parameter, public :: propertyTypeNodeIndex             = 3
  integer, parameter, public :: propertyTypeDescendentIndex       = 4
  integer, parameter, public :: propertyTypeHostIndex             = 5
  integer, parameter, public :: propertyTypeRedshift              = 6
  integer, parameter, public :: propertyTypeNodeMass              = 7
  integer, parameter, public :: propertyTypeParticleCount         = 8
  integer, parameter, public :: propertyTypePositionX             = 9
  integer, parameter, public :: propertyTypePositionY             =10
  integer, parameter, public :: propertyTypePositionZ             =11
  integer, parameter, public :: propertyTypeVelocityX             =12
  integer, parameter, public :: propertyTypeVelocityY             =13
  integer, parameter, public :: propertyTypeVelocityZ             =14
  integer, parameter, public :: propertyTypeSpinX                 =15
  integer, parameter, public :: propertyTypeSpinY                 =16
  integer, parameter, public :: propertyTypeSpinZ                 =17
  integer, parameter, public :: propertyTypeSpin                  =18
  integer, parameter, public :: propertyTypeAngularMomentumX      =19
  integer, parameter, public :: propertyTypeAngularMomentumY      =20
  integer, parameter, public :: propertyTypeAngularMomentumZ      =21
  integer, parameter, public :: propertyTypeAngularMomentum       =22
  integer, parameter, public :: propertyTypeHalfMassRadius        =23
  integer, parameter, public :: propertyTypeParticleIndex         =24
  integer, parameter, public :: propertyTypeMostBoundParticleIndex=25
  integer, parameter, public :: propertyTypeSnapshot              =26

  ! Property names.
  character(len=*), parameter :: propertyNames(26)=[ &
       & 'null                  ',                   &
       & 'treeIndex             ',                   &
       & 'nodeIndex             ',                   &
       & 'descendentIndex       ',                   &
       & 'hostIndex             ',                   &
       & 'redshift              ',                   &
       & 'nodeMass              ',                   &
       & 'particleCount         ',                   &
       & 'positionX             ',                   &
       & 'positionY             ',                   &
       & 'positionZ             ',                   &
       & 'velocityX             ',                   &
       & 'velocityY             ',                   &
       & 'velocityZ             ',                   &
       & 'spinX                 ',                   &
       & 'spinY                 ',                   &
       & 'spinZ                 ',                   &
       & 'spin                  ',                   &
       & 'angularMomentumX      ',                   &
       & 'angularMomentumY      ',                   &
       & 'angularMomentumZ      ',                   &
       & 'angularMomentum       ',                   &
       & 'halfMassRadius        ',                   &
       & 'particleIndex         ',                   &
       & 'mostBoundParticleIndex',                   &
       & 'snapshot              '                    &
       &                                           ]

  ! Names of 3-D datasets (i.e. those which give properties in 3-D space).
  character(len=*), parameter :: propertyNames3D(4)=[ &
       & 'position       ',                           &
       & 'velocity       ',                           &
       & 'spin           ',                           &
       & 'angularMomentum'                            &
       &                                            ]

  ! Units labels.
  integer, parameter         :: unitTypeCount=4
  integer, parameter, public :: unitsMass    =1
  integer, parameter, public :: unitsLength  =2
  integer, parameter, public :: unitsTime    =3
  integer, parameter, public :: unitsVelocity=4

  type unitsMetaData
     !% A structure that holds metadata on units used.
     double precision     :: unitsInSI
     integer              :: hubbleExponent,scaleFactorExponent
     type(varying_string) :: name
  end type unitsMetaData

  ! Metadata labels.
  integer, parameter         :: metaDataTypeCount  =6 
  integer, parameter, public :: metaDataGeneric    =1
  integer, parameter, public :: metaDataCosmology  =2
  integer, parameter, public :: metaDataSimulation =3
  integer, parameter, public :: metaDataGroupFinder=4
  integer, parameter, public :: metaDataTreeBuilder=5
  integer, parameter, public :: metaDataProvenance =6

  ! Data types for metadata.
  integer, parameter         :: dataTypeNull       =0
  integer, parameter         :: dataTypeInteger    =1
  integer, parameter         :: dataTypeDouble     =2
  integer, parameter         :: dataTypeText       =3

  type treeMetaData
     !% Structure that holds metadata for the trees.
     integer              :: metadataType
     type(varying_string) :: label
     integer              :: dataType
     integer              :: integerAttribute
     double precision     :: doubleAttribute
     type(varying_string) :: textAttribute
  end type treeMetaData

  ! Labels for operators used in pruning conditions.
  integer, parameter, public :: operatorEqual             =1
  integer, parameter, public :: operatorNotEqual          =2
  integer, parameter, public :: operatorGreaterThan       =3
  integer, parameter, public :: operatorLessThan          =4
  integer, parameter, public :: operatorGreaterThanOrEqual=5
  integer, parameter, public :: operatorLessThanOrEqual   =6

  type mergerTreeData
     !% A structure that holds raw merger tree data.
     private
     integer                                              :: treeCount,nodeCount,particlesCount
     double precision                                     :: particleMass=0.0d0
     integer,                 allocatable, dimension(:)   :: columnProperties,particleColumnProperties,treeBeginsAt,treeNodeCount
     integer(kind=kind_int8), allocatable, dimension(:)   :: nodeIndex,treeIndex,descendentIndex,hostIndex,particleCount,treeID&
          &,particleIndex,mostBoundParticleIndex,particleReferenceStart,particleReferenceCount,snapshot,particleSnapshot
     double precision,        allocatable, dimension(:)   :: redshift,nodeMass,halfMassRadius,spinMagnitude&
          &,angularMomentumMagnitude,particleRedshift
     double precision,        allocatable, dimension(:,:) :: position,velocity,spin,angularMomentum,particlePosition&
          &,particleVelocity
     logical                                              :: hasNodeIndex,hasTreeIndex,hasDescendentIndex,hasHostIndex &
          &,hasRedshift ,hasNodeMass,hasPositionX,hasPositionY,hasPositionZ,hasVelocityX,hasVelocityY,hasVelocityZ &
          &,hasParticleCount,hasSpinX,hasSpinY,hasSpinZ,hasSpinMagnitude,hasAngularMomentumX,hasAngularMomentumY &
          &,hasAngularMomentumZ,hasAngularMomentumMagnitude,hasHalfMassRadius,hasParticleRedshift,hasParticlePositionX &
          &,hasParticlePositionY,hasParticlePositionZ,hasParticleVelocityX,hasParticleVelocityY,hasParticleVelocityZ&
          &,hasParticleIndex,hasParticles=.false.,hasMostBoundParticleIndex,hasSnapshot,hasParticleSnapshot
     logical                                              :: areSelfContained=.true., includesHubbleFlow=.false.,&
          & includesSubhaloMasses=.false.,doMakeReferences=.true., isPeriodic=.false.
     type(unitsMetaData),     dimension(unitTypeCount)    :: units
     logical,                 dimension(unitTypeCount)    :: unitsSet=.false.
     integer                                              :: metaDataCount=0
     type(treeMetaData),      allocatable, dimension(:)   :: metaData
   contains
     procedure                                            :: reset                     => Merger_Tree_Data_Structure_Reset
     procedure                                            :: treeCountSet              => Merger_Tree_Data_Structure_Set_Tree_Count
     procedure                                            :: nodeCountSet              => Merger_Tree_Data_Structure_Set_Node_Count
     procedure                                            :: particleCountSet          => Merger_Tree_Data_Structure_Set_Particle_Count
     procedure                                            :: readASCII                 => Merger_Tree_Data_Structure_Read_ASCII
     procedure                                            :: readParticlesASCII        => Merger_Tree_Data_Structure_Read_Particles_ASCII
     procedure                                            ::                              Merger_Tree_Data_Structure_Set_Property_Integer8
     procedure                                            ::                              Merger_Tree_Data_Structure_Set_Property_Double
     generic                                              :: setProperty               => Merger_Tree_Data_Structure_Set_Property_Integer8, &
          &                                                                               Merger_Tree_Data_Structure_Set_Property_Double
     procedure                                            :: setPropertyColumn         => Merger_Tree_Data_Structure_Set_Property_Column
     procedure                                            :: setParticlePropertyColumn => Merger_Tree_Data_Structure_Set_Particle_Property_Column
     procedure                                            :: setParticleMass           => Merger_Tree_Data_Structure_Set_Particle_Mass
     procedure                                            :: setSelfContained          => Merger_Tree_Data_Structure_Set_Self_Contained
     procedure                                            :: setIncludesHubbleFlow     => Merger_Tree_Data_Structure_Set_Includes_Hubble_Flow
     procedure                                            :: setPositionsArePeriodic   => Merger_Tree_Data_Structure_Set_Is_Periodic
     procedure                                            :: setIncludesSubhaloMasses  => Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses
     procedure                                            :: setUnits                  => Merger_Tree_Data_Structure_Set_Units
     procedure                                            ::                              Merger_Tree_Data_Structure_Add_Metadata_Double
     procedure                                            ::                              Merger_Tree_Data_Structure_Add_Metadata_Integer
     procedure                                            ::                              Merger_Tree_Data_Structure_Add_Metadata_Text
     generic                                              :: addMetadata               => Merger_Tree_Data_Structure_Add_Metadata_Double , &
          &                                                                               Merger_Tree_Data_Structure_Add_Metadata_Integer, &
          &                                                                               Merger_Tree_Data_Structure_Add_Metadata_Text
     procedure                                            :: makeReferences            => Merger_Tree_Data_Structure_Make_References
     procedure                                            :: export                    => Merger_Tree_Data_Structure_Export
  end type mergerTreeData

contains

  subroutine Merger_Tree_Data_Structure_Reset(mergerTrees)
    !% Reset a merger tree data object.
    use Memory_Management
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees

    ! No properties.
    mergerTrees%hasTreeIndex               =.false.
    mergerTrees%hasNodeIndex               =.false.
    mergerTrees%hasDescendentIndex         =.false.
    mergerTrees%hasHostIndex               =.false.
    mergerTrees%hasRedshift                =.false.
    mergerTrees%hasNodeMass                =.false.
    mergerTrees%hasParticleCount           =.false.
    mergerTrees%hasPositionX               =.false.
    mergerTrees%hasPositionY               =.false.
    mergerTrees%hasPositionZ               =.false.
    mergerTrees%hasVelocityX               =.false.
    mergerTrees%hasVelocityY               =.false.
    mergerTrees%hasVelocityZ               =.false.
    mergerTrees%hasSpinX                   =.false.
    mergerTrees%hasSpinY                   =.false.
    mergerTrees%hasSpinZ                   =.false.
    mergerTrees%hasSpinMagnitude           =.false.
    mergerTrees%hasAngularMomentumX        =.false.
    mergerTrees%hasAngularMomentumY        =.false.
    mergerTrees%hasAngularMomentumZ        =.false.
    mergerTrees%hasAngularMomentumMagnitude=.false.
    mergerTrees%hasHalfMassRadius          =.false.
    mergerTrees%hasMostBoundParticleIndex  =.false.

    ! Deallocate any previous data.
    if (allocated(mergerTrees%treeBeginsAt          )) call Dealloc_Array(mergerTrees%treeBeginsAt          )
    if (allocated(mergerTrees%treeNodeCount         )) call Dealloc_Array(mergerTrees%treeNodeCount         )
    if (allocated(mergerTrees%treeID                )) call Dealloc_Array(mergerTrees%treeID                )
    if (allocated(mergerTrees%treeIndex             )) call Dealloc_Array(mergerTrees%treeIndex             )
    if (allocated(mergerTrees%nodeIndex             )) call Dealloc_Array(mergerTrees%nodeIndex             )
    if (allocated(mergerTrees%mostBoundParticleIndex)) call Dealloc_Array(mergerTrees%mostBoundParticleIndex)
    if (allocated(mergerTrees%descendentIndex       )) call Dealloc_Array(mergerTrees%descendentIndex       )
    if (allocated(mergerTrees%hostIndex             )) call Dealloc_Array(mergerTrees%hostIndex             )
    if (allocated(mergerTrees%redshift              )) call Dealloc_Array(mergerTrees%redshift              )
    if (allocated(mergerTrees%nodeMass              )) call Dealloc_Array(mergerTrees%nodeMass              )
    if (allocated(mergerTrees%particleCount         )) call Dealloc_Array(mergerTrees%particleCount         )
    if (allocated(mergerTrees%position              )) call Dealloc_Array(mergerTrees%position              )
    if (allocated(mergerTrees%velocity              )) call Dealloc_Array(mergerTrees%velocity              )
    if (allocated(mergerTrees%spin                  )) call Dealloc_Array(mergerTrees%spin                  )
    if (allocated(mergerTrees%halfMassRadius        )) call Dealloc_Array(mergerTrees%halfMassRadius        )
    return
  end subroutine Merger_Tree_Data_Structure_Reset

  subroutine Merger_Tree_Data_Structure_Make_References(mergerTrees,makeReferences)
    !% Specify whether or not to make merger tree dataset references.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    logical,               intent(in)    :: makeReferences

    mergerTrees%doMakeReferences=makeReferences
    return
  end subroutine Merger_Tree_Data_Structure_Make_References

  subroutine Merger_Tree_Data_Structure_Add_Metadata_Double(mergerTrees,metadataType,label,doubleValue)
    !% Add a double metadatum.
    use Memory_Management
    use Galacticus_Error
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    integer,               intent(in)    :: metadataType
    character(len=*),      intent(in)    :: label
    double precision,      intent(in)    :: doubleValue

    call Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,doubleValue=doubleValue)
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata_Double

  subroutine Merger_Tree_Data_Structure_Add_Metadata_Integer(mergerTrees,metadataType,label,integerValue)
    !% Add an integer metadatum.
    use Memory_Management
    use Galacticus_Error
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    integer,               intent(in)    :: metadataType
    character(len=*),      intent(in)    :: label
    integer,               intent(in)    :: integerValue

    call Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,integerValue=integerValue)
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata_Integer

  subroutine Merger_Tree_Data_Structure_Add_Metadata_Text(mergerTrees,metadataType,label,textValue)
    !% Add a double metadatum.
    use Memory_Management
    use Galacticus_Error
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    integer,               intent(in)    :: metadataType
    character(len=*),      intent(in)    :: label
    character(len=*),      intent(in)    :: textValue

    call Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,textValue=textValue)
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata_Text

  subroutine Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,integerValue,doubleValue,textValue)
    !% Add a metadatum.
    use Memory_Management
    use Galacticus_Error
    implicit none
    class(mergerTreeData), intent(inout)              :: mergerTrees
    integer,               intent(in)                 :: metadataType
    character(len=*),      intent(in)                 :: label
    integer,               intent(in),   optional     :: integerValue
    double precision,      intent(in),   optional     :: doubleValue
    character(len=*),      intent(in),   optional     :: textValue
    integer,               parameter                  :: metadataBlockSize=100
    type(treeMetaData),    allocatable,  dimension(:) :: metaDataTemporary

    ! Validate the metadata type.
    if (metadataType < 1 .or. metadataType > metaDataTypeCount) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Add_Metadata','invalid metadata type')

    ! Ensure we have enough space in the metadata properties array.
    if (mergerTrees%metaDataCount == 0) then
       allocate(mergerTrees%metaData(metadataBlockSize))
       call Memory_Usage_Record(sizeof(mergerTrees%metaData))
    else if (mergerTrees%metaDataCount == mergerTrees%metaDataCount) then
       call Move_Alloc(mergerTrees%metaData,metaDataTemporary)
       allocate(mergerTrees%metaData(size(metaDataTemporary)+metadataBlockSize))
       call Memory_Usage_Record(sizeof(metaDataTemporary(1)),addRemove=metadataBlockSize,blockCount=0)
       mergerTrees%metaData(1:size(metaDataTemporary))=metaDataTemporary
       deallocate(metaDataTemporary)
    end if

    ! Increment number of metadata.
    mergerTrees%metaDataCount=mergerTrees%metaDataCount+1

    ! Store the type.
    mergerTrees%metaData(mergerTrees%metaDataCount)%metadataType=metadataType

    ! Store the label.
    mergerTrees%metaData(mergerTrees%metaDataCount)%label       =label

    ! Store the data.
    mergerTrees%metaData(mergerTrees%metaDataCount)%dataType=dataTypeNull
    if (present(integerValue)) then
       if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType /= dataTypeNull) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Add_Metadata','only one data type can be specified')
       mergerTrees%metaData(mergerTrees%metaDataCount)%integerAttribute=integerValue
       mergerTrees%metaData(mergerTrees%metaDataCount)%dataType        =dataTypeInteger
    end if
    if (present(doubleValue)) then
       if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType /= dataTypeNull) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Add_Metadata','only one data type can be specified')
       mergerTrees%metaData(mergerTrees%metaDataCount)%doubleAttribute =doubleValue
       mergerTrees%metaData(mergerTrees%metaDataCount)%dataType        =dataTypeDouble
    end if
    if (present(textValue)) then
       if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType /= dataTypeNull) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Add_Metadata','only one data type can be specified')
       mergerTrees%metaData(mergerTrees%metaDataCount)%textAttribute   =textValue
       mergerTrees%metaData(mergerTrees%metaDataCount)%dataType        =dataTypeText
    end if
    if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType == dataTypeNull) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Add_Metadata','no data was given')
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata

  subroutine Merger_Tree_Data_Structure_Set_Tree_Count(mergerTrees,treeCount)
    !% Set the total number of trees in merger trees.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    integer,               intent(in)    :: treeCount

    ! Set the number of trees.
    mergerTrees%treeCount=treeCount

    return
  end subroutine Merger_Tree_Data_Structure_Set_Tree_Count

  subroutine Merger_Tree_Data_Structure_Set_Node_Count(mergerTrees,nodeCount)
    !% Set the total number of nodes in merger trees.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    integer,               intent(in)    :: nodeCount

    ! Set the number of nodes.
    mergerTrees%nodeCount=nodeCount

    return
  end subroutine Merger_Tree_Data_Structure_Set_Node_Count

  subroutine Merger_Tree_Data_Structure_Set_Particle_Count(mergerTrees,particleCount)
    !% Set the total number of particles in merger trees.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    integer,               intent(in)    :: particleCount

    ! Set the number of nodes.
    mergerTrees%particlesCount=particleCount

    return
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Count

  subroutine Merger_Tree_Data_Structure_Set_Particle_Mass(mergerTrees,particleMass)
    !% Set the particle mass used in the trees.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    double precision,      intent(in)    :: particleMass

    ! Set the particle mass.
    mergerTrees%particleMass=particleMass

    return
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Mass

  subroutine Merger_Tree_Data_Structure_Set_Self_Contained(mergerTrees,areSelfContained)
    !% Set the particle mass used in the trees.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    logical,               intent(in)    :: areSelfContained

    ! Set whether trees are self-contained.
    mergerTrees%areSelfContained=areSelfContained

    return
  end subroutine Merger_Tree_Data_Structure_Set_Self_Contained

  subroutine Merger_Tree_Data_Structure_Set_Includes_Hubble_Flow(mergerTrees,includesHubbleFlow)
    !% Set the particle mass used in the trees.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    logical,              intent(in)    :: includesHubbleFlow

    ! Set whether velocities include the Hubble flow.
    mergerTrees%includesHubbleFlow=includesHubbleFlow

    return
  end subroutine Merger_Tree_Data_Structure_Set_Includes_Hubble_Flow

  subroutine Merger_Tree_Data_Structure_Set_Is_Periodic(mergerTrees,isPeriodic)
    !% Set whether or not positions are periodic.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    logical,               intent(in)    :: isPeriodic

    ! Set whether positions are periodic.
    mergerTrees%isPeriodic=isPeriodic

    return
  end subroutine Merger_Tree_Data_Structure_Set_Is_Periodic

  subroutine Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses(mergerTrees,includesSubhaloMasses)
    !% Set the particle mass used in the trees.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    logical,               intent(in)    :: includesSubhaloMasses

    ! Set whether halo masses include subhalo contributions.
    mergerTrees%includesSubhaloMasses=includesSubhaloMasses

    return
  end subroutine Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses

  subroutine Merger_Tree_Data_Structure_Set_Units(mergerTrees,unitType,unitsInSI,hubbleExponent,scaleFactorExponent,name)
    !% Set the units system.
    use Galacticus_Error
    implicit none
    class(mergerTreeData), intent(inout)        :: mergerTrees
    integer,               intent(in)           :: unitType
    double precision,      intent(in)           :: unitsInSI
    integer,               intent(in), optional :: hubbleExponent,scaleFactorExponent
    character(len=*),      intent(in), optional :: name

    ! Ensure the unit type is valid.
    if (unitType < 1 .or. unitType > unitTypeCount) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Set_Units','invalid unit type')

    ! Flag the units as set.
    mergerTrees%unitsSet(unitType)=.true.

    ! Store the units in the SI system.
    mergerTrees%units(unitType)%unitsInSI=unitsInSI
    
    ! Store Hubble parameter exponent if given.
    if (present(hubbleExponent)) then
       mergerTrees%units(unitType)%hubbleExponent=hubbleExponent
    else
       ! No Hubble parameter exponent provided - assume no dependence.
       mergerTrees%units(unitType)%hubbleExponent=0
    end if

    ! Store scale factor exponent if given.
    if (present(scaleFactorExponent)) then
       mergerTrees%units(unitType)%scaleFactorExponent=scaleFactorExponent
    else
       ! No scale factor parameter exponent provided - assume no dependence.
       mergerTrees%units(unitType)%scaleFactorExponent=0
    end if

    ! Store the name if given.
    if (present(name)) then
       mergerTrees%units(unitType)%name=name
    else
       mergerTrees%units(unitType)%name=""
    end if

    return
  end subroutine Merger_Tree_Data_Structure_Set_Units

  subroutine Merger_Tree_Data_Structure_Set_Particle_Property_Column(mergerTrees,propertyType,columnNumber)
    !% Set column mapping from the input file.
    use Memory_Management
    implicit none
    class(mergerTreeData), intent(inout)              :: mergerTrees
    integer,               intent(in)                 :: propertyType,columnNumber
    integer,               allocatable,  dimension(:) :: columnPropertiesTemp

    ! Ensure the storage array is large enough.
    if (allocated(mergerTrees%particleColumnProperties)) then
       if (columnNumber > size(mergerTrees%particleColumnProperties)) then
          call Move_Alloc(mergerTrees%particleColumnProperties,columnPropertiesTemp)
          call Alloc_Array(mergerTrees%particleColumnProperties,[columnNumber])
          mergerTrees%particleColumnProperties(1                        :size(columnPropertiesTemp))=columnPropertiesTemp
          mergerTrees%particleColumnProperties(1+size(columnPropertiesTemp):columnNumber           )=propertyTypeNull
          call Dealloc_Array(columnPropertiesTemp)
       end if
    else
       call Alloc_Array(mergerTrees%particleColumnProperties,[columnNumber])
       mergerTrees%particleColumnProperties=propertyTypeNull
    end if
    ! Store the property type.
    mergerTrees%particleColumnProperties(columnNumber)=propertyType
    return
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Property_Column

  subroutine Merger_Tree_Data_Structure_Set_Property_Column(mergerTrees,propertyType,columnNumber)
    !% Set column mapping from the input file.
    use Memory_Management
    implicit none
    class(mergerTreeData), intent(inout)              :: mergerTrees
    integer,               intent(in)                 :: propertyType,columnNumber
    integer,               allocatable,  dimension(:) :: columnPropertiesTemp

    ! Ensure the storage array is large enough.
    if (allocated(mergerTrees%columnProperties)) then
       if (columnNumber > size(mergerTrees%columnProperties)) then
          call Move_Alloc(mergerTrees%columnProperties,columnPropertiesTemp)
          call Alloc_Array(mergerTrees%columnProperties,[columnNumber])
          mergerTrees%columnProperties(1                        :size(columnPropertiesTemp))=columnPropertiesTemp
          mergerTrees%columnProperties(1+size(columnPropertiesTemp):columnNumber           )=propertyTypeNull
          call Dealloc_Array(columnPropertiesTemp)
       end if
    else
       call Alloc_Array(mergerTrees%columnProperties,[columnNumber])
       mergerTrees%columnProperties=propertyTypeNull
    end if
    ! Store the property type.
    mergerTrees%columnProperties(columnNumber)=propertyType
    return
  end subroutine Merger_Tree_Data_Structure_Set_Property_Column

  subroutine Merger_Tree_Data_Structure_Set_Property_Integer8(mergerTrees,propertyType,property)
    !% Set a property in the merger trees.
    use Memory_Management
    use Galacticus_Error
    implicit none
    class  (mergerTreeData), intent(inout)               :: mergerTrees
    integer                , intent(in   )               :: propertyType
    integer(kind=kind_int8), intent(in   ), dimension(:) :: property

    ! Check the supplied arrays is of the correct size.
    if (size(property) /= mergerTrees%nodeCount) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Set_Property_Integer8','property array size is incorrect')

    ! Assign to the relevant property.
    select case (propertyType)
    case (propertyTypeTreeIndex      )
       mergerTrees%hasTreeIndex      =.true.
       if (allocated(mergerTrees%treeIndex      )) call Dealloc_Array(mergerTrees%treeIndex      )
       call Alloc_Array(mergerTrees%treeIndex      ,[size(property)])
       mergerTrees%treeIndex      =property
       call Merger_Tree_Data_Structure_Set_Tree_Indices(mergerTrees)
    case (propertyTypeNodeIndex      )
       mergerTrees%hasNodeIndex      =.true.
       if (allocated(mergerTrees%nodeIndex      )) call Dealloc_Array(mergerTrees%nodeIndex      )
       call Alloc_Array(mergerTrees%nodeIndex      ,[size(property)])
       mergerTrees%nodeIndex      =property
    case (propertyTypeHostIndex      )
       mergerTrees%hasHostIndex      =.true.
       if (allocated(mergerTrees%hostIndex      )) call Dealloc_Array(mergerTrees%hostIndex      )
       call Alloc_Array(mergerTrees%hostIndex      ,[size(property)])
       mergerTrees%hostIndex      =property
    case (propertyTypeDescendentIndex)
       mergerTrees%hasDescendentIndex=.true.
       if (allocated(mergerTrees%descendentIndex)) call Dealloc_Array(mergerTrees%descendentIndex)
       call Alloc_Array(mergerTrees%descendentIndex,[size(property)])
       mergerTrees%descendentIndex=property
    case (propertyTypeSnapshot       )
       mergerTrees%hasSnapshot       =.true.
       if (allocated(mergerTrees%snapshot       )) call Dealloc_Array(mergerTrees%snapshot       )
       call Alloc_Array(mergerTrees%snapshot       ,[size(property)])
       mergerTrees%snapshot=property
    case default
       call Galacticus_Error_Report('Merger_Tree_Data_Structure_Set_Property_Integer8','unrecognized integer property')  
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Set_Property_Integer8

  subroutine Merger_Tree_Data_Structure_Set_Property_Double(mergerTrees,propertyType,property)
    !% Set a property in the merger trees.
    use Memory_Management
    use Galacticus_Error
    implicit none
    class  (mergerTreeData), intent(inout)               :: mergerTrees
    integer                , intent(in   )               :: propertyType
    double precision       , intent(in   ), dimension(:) :: property

    ! Check the supplied arrays is of the correct size.
    if (size(property) /= mergerTrees%nodeCount) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Set_Property_Double','property array size is incorrect')

    ! Assign to the relevant property.
    select case (propertyType)
    case (propertyTypeRedshift       )
       mergerTrees%hasRedshift       =.true.
       if (allocated(mergerTrees%redshift       )) call Dealloc_Array(mergerTrees%redshift       )
       call Alloc_Array(mergerTrees%redshift       ,[size(property)])
       mergerTrees%redshift       =property
    case (propertyTypeNodeMass       )
       mergerTrees%hasNodeMass       =.true.
       if (allocated(mergerTrees%nodeMass       )) call Dealloc_Array(mergerTrees%nodeMass       )
       call Alloc_Array(mergerTrees%nodeMass       ,[size(property)])
       mergerTrees%nodeMass       =property
    case (propertyTypePositionX      )
       if (                                             &
            & allocated(mergerTrees%position)           &
            & .and..not.mergerTrees%hasPositionX        &
            & .and..not.mergerTrees%hasPositionY        &
            & .and..not.mergerTrees%hasPositionZ        &
            &) call Dealloc_Array(mergerTrees%position)
       mergerTrees%hasPositionX=.true.
       if (.not.allocated(mergerTrees%position)) call Alloc_Array(mergerTrees%position,[size(property),3])
       mergerTrees%position(:,1)=property
    case (propertyTypePositionY      )
       if (                                             &
            & allocated(mergerTrees%position)           &
            & .and..not.mergerTrees%hasPositionX        &
            & .and..not.mergerTrees%hasPositionY        &
            & .and..not.mergerTrees%hasPositionZ        &
            &) call Dealloc_Array(mergerTrees%position)
       mergerTrees%hasPositionY=.true.
       if (.not.allocated(mergerTrees%position)) call Alloc_Array(mergerTrees%position,[size(property),3])
       mergerTrees%position(:,2)=property
    case (propertyTypePositionZ      )
       if (                                             &
            & allocated(mergerTrees%position)           &
            & .and..not.mergerTrees%hasPositionX        &
            & .and..not.mergerTrees%hasPositionY        &
            & .and..not.mergerTrees%hasPositionZ        &
            &) call Dealloc_Array(mergerTrees%position)
       mergerTrees%hasPositionZ=.true.
       if (.not.allocated(mergerTrees%position)) call Alloc_Array(mergerTrees%position,[size(property),3])
       mergerTrees%position(:,3)=property
    case (propertyTypeVelocityX      )
       if (                                             &
            & allocated(mergerTrees%velocity)           &
            & .and..not.mergerTrees%hasVelocityX        &
            & .and..not.mergerTrees%hasVelocityY        &
            & .and..not.mergerTrees%hasVelocityZ        &
            &) call Dealloc_Array(mergerTrees%velocity)
       mergerTrees%hasVelocityX=.true.
       if (.not.allocated(mergerTrees%velocity)) call Alloc_Array(mergerTrees%velocity,[size(property),3])
       mergerTrees%velocity(:,1)=property
    case (propertyTypeVelocityY      )
       if (                                             &
            & allocated(mergerTrees%velocity)           &
            & .and..not.mergerTrees%hasVelocityX        &
            & .and..not.mergerTrees%hasVelocityY        &
            & .and..not.mergerTrees%hasVelocityZ        &
            &) call Dealloc_Array(mergerTrees%velocity)
       mergerTrees%hasVelocityY=.true.
       if (.not.allocated(mergerTrees%velocity)) call Alloc_Array(mergerTrees%velocity,[size(property),3])
       mergerTrees%velocity(:,2)=property
    case (propertyTypeVelocityZ      )
       if (                                             &
            & allocated(mergerTrees%velocity)           &
            & .and..not.mergerTrees%hasVelocityX        &
            & .and..not.mergerTrees%hasVelocityY        &
            & .and..not.mergerTrees%hasVelocityZ        &
            &) call Dealloc_Array(mergerTrees%velocity)
       mergerTrees%hasVelocityZ=.true.
       if (.not.allocated(mergerTrees%velocity)) call Alloc_Array(mergerTrees%velocity,[size(property),3])
       mergerTrees%velocity(:,3)=property
    case default
       call Galacticus_Error_Report('Merger_Tree_Data_Structure_Set_Property_Double','unrecognized double property')  
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Set_Property_Double
  
  subroutine Merger_Tree_Data_Structure_Read_ASCII(mergerTrees,inputFile,lineNumberStart,lineNumberStop,separator)
    !% Read in merger tree data from an ASCII file.
    use String_Handling
    use Memory_Management
    use Galacticus_Error
    use Galacticus_Display
    use File_Utilities
    implicit none
    class(mergerTreeData), intent(inout)              :: mergerTrees
    character(len=*),      intent(in)                 :: inputFile
    integer,               intent(in),   optional     :: lineNumberStart,lineNumberStop
    character(len=*),      intent(in),   optional     :: separator
    character(len=32),     allocatable,  dimension(:) :: inputColumns
    integer                                           :: lineNumberStartActual,lineNumberStopActual,columnsCount,lineNumber&
         &,fileUnit,iColumn,iNode
    logical                                           :: gotFirstDataLine
    character(len=1024)                               :: inputLine
    type(varying_string)                              :: message

    ! Get start and stop line numbers.
    if (present(lineNumberStart)) then
       lineNumberStartActual=lineNumberStart
    else
       lineNumberStartActual=1
    end if
    if (present(lineNumberStop )) then
       lineNumberStopActual =lineNumberStop
    else
       lineNumberStopActual =Count_Lines_In_File(inputFile)
    end if

    ! Determine number of nodes.
    call mergerTrees%nodeCountSet(lineNumberStopActual-lineNumberStartActual+1)

    ! Specify what properties these trees have.
    mergerTrees%hasTreeIndex               =any(mergerTrees%columnProperties == propertyTypeTreeIndex             )
    mergerTrees%hasNodeIndex               =any(mergerTrees%columnProperties == propertyTypeNodeIndex             )
    mergerTrees%hasDescendentIndex         =any(mergerTrees%columnProperties == propertyTypeDescendentIndex       )
    mergerTrees%hasHostIndex               =any(mergerTrees%columnProperties == propertyTypeHostIndex             )
    mergerTrees%hasRedshift                =any(mergerTrees%columnProperties == propertyTypeRedshift              )
    mergerTrees%hasNodeMass                =any(mergerTrees%columnProperties == propertyTypeNodeMass              )
    mergerTrees%hasParticleCount           =any(mergerTrees%columnProperties == propertyTypeParticleCount         )
    mergerTrees%hasPositionX               =any(mergerTrees%columnProperties == propertyTypePositionX             )
    mergerTrees%hasPositionY               =any(mergerTrees%columnProperties == propertyTypePositionY             )
    mergerTrees%hasPositionZ               =any(mergerTrees%columnProperties == propertyTypePositionZ             )
    mergerTrees%hasVelocityX               =any(mergerTrees%columnProperties == propertyTypeVelocityX             )
    mergerTrees%hasVelocityY               =any(mergerTrees%columnProperties == propertyTypeVelocityY             )
    mergerTrees%hasVelocityZ               =any(mergerTrees%columnProperties == propertyTypeVelocityZ             )
    mergerTrees%hasSpinX                   =any(mergerTrees%columnProperties == propertyTypeSpinX                 )
    mergerTrees%hasSpinY                   =any(mergerTrees%columnProperties == propertyTypeSpinY                 )
    mergerTrees%hasSpinZ                   =any(mergerTrees%columnProperties == propertyTypeSpinZ                 )
    mergerTrees%hasSpinMagnitude           =any(mergerTrees%columnProperties == propertyTypeSpin                  )
    mergerTrees%hasAngularMomentumX        =any(mergerTrees%columnProperties == propertyTypeAngularMomentumX      )
    mergerTrees%hasAngularMomentumY        =any(mergerTrees%columnProperties == propertyTypeAngularMomentumY      )
    mergerTrees%hasAngularMomentumZ        =any(mergerTrees%columnProperties == propertyTypeAngularMomentumZ      )
    mergerTrees%hasAngularMomentumMagnitude=any(mergerTrees%columnProperties == propertyTypeAngularMomentum       )
    mergerTrees%hasHalfMassRadius          =any(mergerTrees%columnProperties == propertyTypeHalfMassRadius        )
    mergerTrees%hasMostBoundParticleIndex  =any(mergerTrees%columnProperties == propertyTypeMostBoundParticleIndex)
    mergerTrees%hasSnapshot                =any(mergerTrees%columnProperties == propertyTypeSnapshot              )

    ! Validate 3-D datasets.
    if     (.not.((     mergerTrees%hasPositionX       .and.     mergerTrees%hasPositionY       .and.     mergerTrees%hasPositionZ       ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasPositionX       .and..not.mergerTrees%hasPositionY       .and..not.mergerTrees%hasPositionZ       ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Read_ASCII","all three axes or none must be supplied for position"        )
    if     (.not.((     mergerTrees%hasVelocityX       .and.     mergerTrees%hasVelocityY       .and.     mergerTrees%hasVelocityZ       ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasVelocityX       .and..not.mergerTrees%hasVelocityY       .and..not.mergerTrees%hasVelocityZ       ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Read_ASCII","all three axes or none must be supplied for velocity"        )
    if     (.not.((     mergerTrees%hasSpinX           .and.     mergerTrees%hasSpinY           .and.     mergerTrees%hasSpinZ           ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasSpinX           .and..not.mergerTrees%hasSpinY           .and..not.mergerTrees%hasSpinZ           ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Read_ASCII","all three axes or none must be supplied for spin"            )
    if     (.not.((     mergerTrees%hasAngularMomentumX.and.     mergerTrees%hasAngularMomentumY.and.     mergerTrees%hasAngularMomentumZ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasAngularMomentumX.and..not.mergerTrees%hasAngularMomentumY.and..not.mergerTrees%hasAngularMomentumZ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Read_ASCII","all three axes or none must be supplied for angular momentum")
    if (mergerTrees%hasSpinX           .and.mergerTrees%hasSpinMagnitude           ) call &
         & Galacticus_Error_Report("Merger_Tree_Data_Structure_Read_ASCII","can not specify both 3-D and scalar angular momentum") 
    if (mergerTrees%hasAngularMomentumX.and.mergerTrees%hasAngularMomentumMagnitude) call &
         & Galacticus_Error_Report("Merger_Tree_Data_Structure_Read_ASCII","can not specify both 3-D and scalar angular momentum") 

    ! Deallocate internal arrays.
    if (allocated(mergerTrees%treeIndex             )) call Dealloc_Array(mergerTrees%treeIndex             )
    if (allocated(mergerTrees%nodeIndex             )) call Dealloc_Array(mergerTrees%nodeIndex             )
    if (allocated(mergerTrees%mostBoundParticleIndex)) call Dealloc_Array(mergerTrees%mostBoundParticleIndex)
    if (allocated(mergerTrees%snapshot              )) call Dealloc_Array(mergerTrees%snapshot              )
    if (allocated(mergerTrees%descendentIndex       )) call Dealloc_Array(mergerTrees%descendentIndex       )
    if (allocated(mergerTrees%hostIndex             )) call Dealloc_Array(mergerTrees%hostIndex             )
    if (allocated(mergerTrees%redshift              )) call Dealloc_Array(mergerTrees%redshift              )
    if (allocated(mergerTrees%nodeMass              )) call Dealloc_Array(mergerTrees%nodeMass              )
    if (allocated(mergerTrees%particleCount         )) call Dealloc_Array(mergerTrees%particleCount         )
    if (allocated(mergerTrees%position              )) call Dealloc_Array(mergerTrees%position              )
    if (allocated(mergerTrees%velocity              )) call Dealloc_Array(mergerTrees%velocity              )
    if (allocated(mergerTrees%spin                  )) call Dealloc_Array(mergerTrees%spin                  )
    if (allocated(mergerTrees%halfMassRadius        )) call Dealloc_Array(mergerTrees%halfMassRadius        )

    ! Allocate internal arrays to correct size as needed.
    if (mergerTrees%hasTreeIndex               ) call Alloc_Array(mergerTrees%treeIndex               ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasNodeIndex               ) call Alloc_Array(mergerTrees%nodeIndex               ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasMostBoundParticleIndex  ) call Alloc_Array(mergerTrees%mostBoundParticleIndex  ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasSnapshot                ) call Alloc_Array(mergerTrees%snapshot                ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasDescendentIndex         ) call Alloc_Array(mergerTrees%descendentIndex         ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasHostIndex               ) call Alloc_Array(mergerTrees%hostIndex               ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasRedshift                ) call Alloc_Array(mergerTrees%redshift                ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasNodeMass                ) call Alloc_Array(mergerTrees%nodeMass                ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasParticleCount           ) call Alloc_Array(mergerTrees%particleCount           ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasPositionX               ) call Alloc_Array(mergerTrees%position                ,[3,mergerTrees%nodeCount])
    if (mergerTrees%hasVelocityX               ) call Alloc_Array(mergerTrees%velocity                ,[3,mergerTrees%nodeCount])
    if (mergerTrees%hasSpinX                   ) call Alloc_Array(mergerTrees%spin                    ,[3,mergerTrees%nodeCount])
    if (mergerTrees%hasAngularMomentumX        ) call Alloc_Array(mergerTrees%angularMomentum         ,[3,mergerTrees%nodeCount])
    if (mergerTrees%hasSpinMagnitude           ) call Alloc_Array(mergerTrees%spinMagnitude           ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasAngularMomentumMagnitude) call Alloc_Array(mergerTrees%angularMomentumMagnitude,[  mergerTrees%nodeCount])
    if (mergerTrees%hasHalfMassRadius          ) call Alloc_Array(mergerTrees%halfMassRadius          ,[  mergerTrees%nodeCount])

    ! Open the file and read lines.
    open(newunit=fileUnit,file=inputFile,status='old',form='formatted')
    lineNumber      =0
    iNode           =0
    gotFirstDataLine=.false.
    do while (lineNumber < lineNumberStopActual)
       ! Get the line.
       read (fileUnit,'(a)') inputLine
       ! Increment line number count.
       lineNumber=lineNumber+1
       ! Check if this is a data line.
       if (lineNumber >= lineNumberStartActual .and. lineNumber <= lineNumberStopActual) then
          ! Count nodes.
          iNode=iNode+1
          ! If this is the first data line, determine how many columns are present and allocate array to store them.
          if (.not.gotFirstDataLine) then
             columnsCount=String_Count_Words(inputLine,separator)
             call Alloc_Array(inputColumns,[columnsCount])
             gotFirstDataLine=.true.
          end if
          call String_Split_Words(inputColumns,inputLine,separator)
          do iColumn=1,min(columnsCount,size(mergerTrees%columnProperties))
             select case (mergerTrees%columnProperties(iColumn))
             case (propertyTypeNull                  )
                ! Ignore this column.
             case (propertyTypeTreeIndex             )
                ! Column is a tree index.
                read (inputColumns(iColumn),*) mergerTrees%treeIndex(iNode)
                if (iNode > 1) then
                   if (mergerTrees%treeIndex(iNode) < mergerTrees%treeIndex(iNode-1)) call&
                        & Galacticus_Error_Report('Merger_Tree_Data_Structure_Read_ASCII','tree indices must be in ascending &
                        & order')
                   if (mergerTrees%treeIndex(iNode) /= mergerTrees%treeIndex(iNode-1)) mergerTrees%treeCount&
                        &=mergerTrees%treeCount+1                      
                else
                   mergerTrees%treeCount=1
                end if
             case (propertyTypeNodeIndex             )
                ! Column is a node index.
                read (inputColumns(iColumn),*) mergerTrees%nodeIndex      (  iNode)
             case (propertyTypeDescendentIndex       )
                ! Column is a descendent node index.
                read (inputColumns(iColumn),*) mergerTrees%descendentIndex         (  iNode)
             case (propertyTypeHostIndex             )
                ! Column is a host index.
                read (inputColumns(iColumn),*) mergerTrees%hostIndex               (  iNode)
             case (propertyTypeRedshift              )
                ! Column is redshift.
                read (inputColumns(iColumn),*) mergerTrees%redshift                (  iNode)
             case (propertyTypeNodeMass              )
                ! Column is mass.
                read (inputColumns(iColumn),*) mergerTrees%nodeMass                (  iNode)
             case (propertyTypeParticleCount         )
                ! Column is particle count.
                read (inputColumns(iColumn),*) mergerTrees%particleCount           (  iNode)
             case (propertyTypePositionX             )
                ! Column is x position.
                read (inputColumns(iColumn),*) mergerTrees%position                (1,iNode)
             case (propertyTypePositionY             )
                ! Column is y position.
                read (inputColumns(iColumn),*) mergerTrees%position                (2,iNode)
             case (propertyTypePositionZ             )
                ! Column is z position.
                read (inputColumns(iColumn),*) mergerTrees%position                (3,iNode)
             case (propertyTypeVelocityX             )
                ! Column is x velocity.
                read (inputColumns(iColumn),*) mergerTrees%velocity                (1,iNode)
             case (propertyTypeVelocityY             )
                ! Column is y velocity.
                read (inputColumns(iColumn),*) mergerTrees%velocity                (2,iNode)
             case (propertyTypeVelocityZ             )
                ! Column is z velocity.
                read (inputColumns(iColumn),*) mergerTrees%velocity                (3,iNode)
             case (propertyTypeSpinX                 )
                ! Column is x spin.
                read (inputColumns(iColumn),*) mergerTrees%spin                    (1,iNode)
             case (propertyTypeSpinY                 )
                ! Column is y spin.
                read (inputColumns(iColumn),*) mergerTrees%spin                    (2,iNode)
             case (propertyTypeSpinZ                 )
                ! Column is z spin.
                read (inputColumns(iColumn),*) mergerTrees%spin                    (3,iNode)
             case (propertyTypeSpin                  )
                ! Column is scalar spin.
                read (inputColumns(iColumn),*) mergerTrees%spinMagnitude           (  iNode)
             case (propertyTypeAngularMomentumX      )
                ! Column is x angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%angularMomentum         (1,iNode)
             case (propertyTypeAngularMomentumY      )
                ! Column is y angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%angularMomentum         (2,iNode)
             case (propertyTypeAngularMomentumZ      )
                ! Column is z angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%angularMomentum         (3,iNode)
             case (propertyTypeAngularMomentum       )
                ! Column is scalar angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%angularMomentumMagnitude(  iNode)
             case (propertyTypeHalfMassRadius        )
                ! Column is half mass radius.
                read (inputColumns(iColumn),*) mergerTrees%halfMassRadius          (  iNode)
             case (propertyTypeMostBoundParticleIndex)
                ! Column is a node index.
                read (inputColumns(iColumn),*) mergerTrees%mostBoundParticleIndex  (  iNode)
             case (propertyTypeSnapshot              )
                ! Column is a snapshot index.
                read (inputColumns(iColumn),*) mergerTrees%snapshot                (  iNode)
             case default
                call Galacticus_Error_Report('Merger_Tree_Data_Structure_Read_ASCII','unknown column type')
             end select
          end do
       end if
    end do
    close(fileUnit)

    ! Report number of trees found.
    message='Found '
    message=message//mergerTrees%treeCount//' merger trees'
    call Galacticus_Display_Message(message)

    ! Deallocate workspace.
    if (allocated(inputColumns)) call Dealloc_Array(inputColumns)

    ! Set tree indices.
    call Merger_Tree_Data_Structure_Set_Tree_Indices(mergerTrees)
    return
  end subroutine Merger_Tree_Data_Structure_Read_ASCII

  subroutine Merger_Tree_Data_Structure_Set_Tree_Indices(mergerTrees)
    !% Set the merger tree index arrays.
    use Memory_Management
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    integer                              :: iTree,iNode

    ! Allocate arrays for tree start and stop indices and reference ID.
    if (allocated(mergerTrees%treeBeginsAt )) call Dealloc_Array(mergerTrees%treeBeginsAt )
    if (allocated(mergerTrees%treeNodeCount)) call Dealloc_Array(mergerTrees%treeNodeCount)
    if (allocated(mergerTrees%treeID       )) call Dealloc_Array(mergerTrees%treeID       )
    call Alloc_Array(mergerTrees%treeBeginsAt ,[mergerTrees%treeCount])
    call Alloc_Array(mergerTrees%treeNodeCount,[mergerTrees%treeCount])
    call Alloc_Array(mergerTrees%treeID       ,[mergerTrees%treeCount])

    ! Determine index in arrays where each tree begins.
    mergerTrees%treeBeginsAt(1)=0
    mergerTrees%treeID      (1)=mergerTrees%treeIndex(1)
    iTree                      =1
    do iNode=2,mergerTrees%nodeCount
       if (mergerTrees%treeIndex(iNode) /= mergerTrees%treeIndex(iNode-1)) then
          iTree=iTree+1
          mergerTrees%treeBeginsAt (iTree  )=iNode-1
          mergerTrees%treeNodeCount(iTree-1)=mergerTrees%treeBeginsAt(iTree)-mergerTrees%treeBeginsAt(iTree-1)
          mergerTrees%treeID       (iTree  )=mergerTrees%treeIndex(iNode)
       end if
    end do
    mergerTrees%treeNodeCount(mergerTrees%treeCount)=mergerTrees%nodeCount-mergerTrees%treeBeginsAt(mergerTrees%treeCount)
    return
  end subroutine Merger_Tree_Data_Structure_Set_Tree_Indices
  
  subroutine Merger_Tree_Data_Structure_Read_Particles_ASCII(mergerTrees,inputFile,lineNumberStart,lineNumberStop,separator)
    !% Read in particle data from an ASCII file.
    use String_Handling
    use Memory_Management
    use Galacticus_Error
    use Galacticus_Display
    use File_Utilities
    implicit none
    class(mergerTreeData), intent(inout)              :: mergerTrees
    character(len=*),      intent(in)                 :: inputFile
    integer,               intent(in),   optional     :: lineNumberStart,lineNumberStop
    character(len=*),      intent(in),   optional     :: separator
    character(len=32),     allocatable,  dimension(:) :: inputColumns
    integer                                           :: lineNumberStartActual,lineNumberStopActual,columnsCount,lineNumber&
         &,fileUnit,iColumn,iNode
    logical                                           :: gotFirstDataLine
    character(len=1024)                               :: inputLine

    ! Flag that these trees have particles.
    mergerTrees%hasParticles=.true.

    ! Get start and stop line numbers.
    if (present(lineNumberStart)) then
       lineNumberStartActual=lineNumberStart
    else
       lineNumberStartActual=1
    end if
    if (present(lineNumberStop )) then
       lineNumberStopActual =lineNumberStop
    else
       lineNumberStopActual =Count_Lines_In_File(inputFile)
    end if

    ! Determine number of particles.
    call mergerTrees%particleCountSet(lineNumberStopActual-lineNumberStartActual+1)

    ! Specify what properties these particles have.
    mergerTrees%hasParticleIndex    =any(mergerTrees%particleColumnProperties == propertyTypeParticleIndex)
    mergerTrees%hasParticleRedshift =any(mergerTrees%particleColumnProperties == propertyTypeRedshift     )
    mergerTrees%hasParticlePositionX=any(mergerTrees%particleColumnProperties == propertyTypePositionX    )
    mergerTrees%hasParticlePositionY=any(mergerTrees%particleColumnProperties == propertyTypePositionY    )
    mergerTrees%hasParticlePositionZ=any(mergerTrees%particleColumnProperties == propertyTypePositionZ    )
    mergerTrees%hasParticleVelocityX=any(mergerTrees%particleColumnProperties == propertyTypeVelocityX    )
    mergerTrees%hasParticleVelocityY=any(mergerTrees%particleColumnProperties == propertyTypeVelocityY    )
    mergerTrees%hasParticleVelocityZ=any(mergerTrees%particleColumnProperties == propertyTypeVelocityZ    )
    mergerTrees%hasParticleSnapshot =any(mergerTrees%particleColumnProperties == propertyTypeSnapshot     )

    ! Validate 3-D datasets.
    if     (.not.((     mergerTrees%hasParticlePositionX.and.     mergerTrees%hasParticlePositionY.and.     mergerTrees%hasParticlePositionZ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasParticlePositionX.and..not.mergerTrees%hasParticlePositionY.and..not.mergerTrees%hasParticlePositionZ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Read_Particle_ASCII","all three axes or none must be supplied for particle position")
    if     (.not.((     mergerTrees%hasParticleVelocityX.and.     mergerTrees%hasParticleVelocityY.and.     mergerTrees%hasParticleVelocityZ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasParticleVelocityX.and..not.mergerTrees%hasParticleVelocityY.and..not.mergerTrees%hasParticleVelocityZ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Read_Particle_ASCII","all three axes or none must be supplied for particle velocity")

    ! Ensure we have a redshift.
    if (.not.mergerTrees%hasParticleRedshift) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Read_Particle_ASCII","particle redshift must be supplied")

    ! Deallocate internal arrays.
    if (allocated(mergerTrees%particleIndex   )) call Dealloc_Array(mergerTrees%particleIndex   )
    if (allocated(mergerTrees%particleRedshift)) call Dealloc_Array(mergerTrees%particleRedshift)
    if (allocated(mergerTrees%particlePosition)) call Dealloc_Array(mergerTrees%particlePosition)
    if (allocated(mergerTrees%particleVelocity)) call Dealloc_Array(mergerTrees%particleVelocity)
    if (allocated(mergerTrees%particleSnapshot)) call Dealloc_Array(mergerTrees%particleSnapshot)
 
    ! Allocate internal arrays to correct size as needed.
    if (mergerTrees%hasParticleIndex    ) call Alloc_Array(mergerTrees%particleIndex   ,[  mergerTrees%particlesCount])
    if (mergerTrees%hasParticleRedshift ) call Alloc_Array(mergerTrees%particleRedshift,[  mergerTrees%particlesCount])
    if (mergerTrees%hasParticlePositionX) call Alloc_Array(mergerTrees%particlePosition,[3,mergerTrees%particlesCount])
    if (mergerTrees%hasParticleVelocityX) call Alloc_Array(mergerTrees%particleVelocity,[3,mergerTrees%particlesCount])
    if (mergerTrees%hasParticleSnapshot ) call Alloc_Array(mergerTrees%particleSnapshot,[  mergerTrees%particlesCount])

    ! Open the file and read lines.
    open(newunit=fileUnit,file=inputFile,status='old',form='formatted')
    lineNumber      =0
    iNode           =0
    gotFirstDataLine=.false.
    do while (lineNumber < lineNumberStopActual)
       ! Get the line.
       read (fileUnit,'(a)') inputLine
       ! Increment line number count.
       lineNumber=lineNumber+1
       ! Check if this is a data line.
       if (lineNumber >= lineNumberStartActual .and. lineNumber <= lineNumberStopActual) then
          ! Count nodes.
          iNode=iNode+1
          ! If this is the first data line, determine how many columns are present and allocate array to store them.
          if (.not.gotFirstDataLine) then
             columnsCount=String_Count_Words(inputLine,separator)
             call Alloc_Array(inputColumns,[columnsCount])
             gotFirstDataLine=.true.
          end if
          call String_Split_Words(inputColumns,inputLine,separator)
          do iColumn=1,min(columnsCount,size(mergerTrees%particleColumnProperties))
             select case (mergerTrees%particleColumnProperties(iColumn))
             case (propertyTypeNull         )
                ! Ignore this column.
             case (propertyTypeParticleIndex)
                ! Column is particle index.
                read (inputColumns(iColumn),*) mergerTrees%particleIndex   (  iNode)
             case (propertyTypeRedshift     )
                ! Column is redshift.
                read (inputColumns(iColumn),*) mergerTrees%particleRedshift(  iNode)
             case (propertyTypeSnapshot     )
                ! Column is snapshot.
                read (inputColumns(iColumn),*) mergerTrees%particleSnapshot(  iNode)
             case (propertyTypePositionX    )
                ! Column is x position.
                read (inputColumns(iColumn),*) mergerTrees%particlePosition(1,iNode)
             case (propertyTypePositionY    )
                ! Column is y position.
                read (inputColumns(iColumn),*) mergerTrees%particlePosition(2,iNode)
             case (propertyTypePositionZ    )
                ! Column is z position.
                read (inputColumns(iColumn),*) mergerTrees%particlePosition(3,iNode)
             case (propertyTypeVelocityX    )
                ! Column is x velocity.
                read (inputColumns(iColumn),*) mergerTrees%particleVelocity(1,iNode)
             case (propertyTypeVelocityY    )
                ! Column is y velocity.
                read (inputColumns(iColumn),*) mergerTrees%particleVelocity(2,iNode)
             case (propertyTypeVelocityZ    )
                ! Column is z velocity.
                read (inputColumns(iColumn),*) mergerTrees%particleVelocity(3,iNode)
             case default
                call Galacticus_Error_Report('Merger_Tree_Data_Structure_Read_Particle_ASCII','unknown column type')
             end select
          end do
       end if
    end do
    close(fileUnit)

    ! Deallocate workspace.
    if (allocated(inputColumns)) call Dealloc_Array(inputColumns)

    return
  end subroutine Merger_Tree_Data_Structure_Read_Particles_ASCII

  subroutine Merger_Tree_Data_Structure_Export(mergerTrees,outputFileName,outputFormat,hdfChunkSize,hdfCompressionLevel,append)
    !% Output a set of merger trees to an HDF5 file.
    use Galacticus_Error
    use String_Handling
    implicit none
    integer                  , intent(in   )           :: hdfChunkSize,hdfCompressionLevel
    class    (mergerTreeData), intent(inout)           :: mergerTrees
    character(len=*         ), intent(in   )           :: outputFileName,outputFormat
    logical                  , intent(in   ), optional :: append

    select case (String_Lower_Case(trim(outputFormat)))
    case ("galacticus")
       call Merger_Tree_Data_Structure_Export_Galacticus(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    case ("irate")
       call Merger_Tree_Data_Structure_Export_IRATE     (mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    case default
       call Galacticus_Error_Report('Merger_Tree_Data_Structure_Export','output format is not recognized')
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Export

  subroutine Merger_Tree_Data_Structure_Export_Galacticus(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    !% Output a set of merger trees to a Galacticus-format HDF5 file.
    use HDF5
    use IO_HDF5
    use Galacticus_Error
    use String_Handling
    use Memory_Management
    implicit none
    integer                  , intent(in   )              :: hdfChunkSize,hdfCompressionLevel
    class    (mergerTreeData), intent(inout)              :: mergerTrees
    character(len=*         ), intent(in   )              :: outputFileName
    logical                  , intent(in   ), optional    :: append
    integer  (kind=HSIZE_T  ),               dimension(2) :: hyperslabStart,hyperslabCount
    type     (hdf5Object    ), pointer                    :: attributeGroup
    type     (hdf5Object    ), target                     :: outputFile,haloTrees,treesGroup,treeGroup,treeDataset,treeIndexGroup&
         &,unitsGroup ,genericGroup,cosmologyGroup,simulationGroup,groupFinderGroup,treeBuilderGroup,provenanceGroup&
         &,particlesGroup
    integer                  , allocatable , dimension(:) :: firstNode,numberOfNodes
    integer                                               :: iTree,iProperty,integerAttribute,iAttribute
    type     (varying_string)                             :: groupName
    logical                                               :: appendActual

    ! Determine if we are to append to an existing file.
    appendActual=.false.
    if (present(append)) appendActual=append

    ! Validate the merger tree.
    call Merger_Tree_Data_Validate_Trees            (mergerTrees)

    ! If we have the particle mass, set the masses of any subhalos (which have zero mass by default) based on particle count.
    call Merger_Tree_Data_Set_Subhalo_Masses        (mergerTrees)

    ! If we have most-bound particle indices and particle data has been read, construct arrays giving position of particle data for each node.
    call Merger_Tree_Data_Construct_Particle_Indices(mergerTrees)

    ! Open the output file.
    !$omp critical (HDF5_Access)
    call outputFile%openFile(outputFileName,overWrite=.not.appendActual,chunkSize=hdfChunkSize,compressionLevel=hdfCompressionLevel)

    ! Create a group for the datasets.
    haloTrees=outputFile%openGroup("haloTrees","Stores all data for merger trees.")

    ! Write the data.
    if (mergerTrees%hasNodeIndex               ) call haloTrees%writeDataset(mergerTrees%nodeIndex               ,"nodeIndex"      ,"The index of each node."                ,appendTo=.true.)
    if (mergerTrees%hasDescendentIndex         ) call haloTrees%writeDataset(mergerTrees%descendentIndex         ,"descendentIndex","The index of each descendent node."     ,appendTo=.true.)
    if (mergerTrees%hasHostIndex               ) call haloTrees%writeDataset(mergerTrees%hostIndex               ,"hostIndex"      ,"The index of each host node."           ,appendTo=.true.)
    if (mergerTrees%hasNodeMass                ) call haloTrees%writeDataset(mergerTrees%nodeMass                ,"nodeMass"       ,"The mass of each node."                 ,appendTo=.true.)
    if (mergerTrees%hasRedshift                ) call haloTrees%writeDataset(mergerTrees%redshift                ,"redshift"       ,"The redshift of each node."             ,appendTo=.true.)
    if (mergerTrees%hasPositionX               ) call haloTrees%writeDataset(mergerTrees%position                ,"position"       ,"The position of each node."             ,appendTo=.true.)
    if (mergerTrees%hasVelocityX               ) call haloTrees%writeDataset(mergerTrees%velocity                ,"velocity"       ,"The velocity of each node."             ,appendTo=.true.)
    if (mergerTrees%hasSpinX                   ) call haloTrees%writeDataset(mergerTrees%spin                    ,"spin"           ,"The spin of each node."                 ,appendTo=.true.)
    if (mergerTrees%hasAngularMomentumX        ) call haloTrees%writeDataset(mergerTrees%angularMomentum         ,"angularMomentum","The angular momentum spin of each node.",appendTo=.true.)
    if (mergerTrees%hasSpinMagnitude           ) call haloTrees%writeDataset(mergerTrees%spinMagnitude           ,"spin"           ,"The spin of each node."                 ,appendTo=.true.)
    if (mergerTrees%hasAngularMomentumMagnitude) call haloTrees%writeDataset(mergerTrees%angularMomentumMagnitude,"angularMomentum","The angular momentum spin of each node.",appendTo=.true.)
    if (mergerTrees%hasHalfMassRadius          ) call haloTrees%writeDataset(mergerTrees%halfMassRadius          ,"halfMassRadius" ,"The half mass radius of each node."     ,appendTo=.true.)
    if (mergerTrees%hasMostBoundParticleIndex) then
       call haloTrees%writeDataset(mergerTrees%particleReferenceStart,"particleIndexStart","The starting index of particle data for each node.",appendTo=.true.)
       call haloTrees%writeDataset(mergerTrees%particleReferenceCount,"particleIndexCount","The number of particle data for each node."        ,appendTo=.true.)
    end if

    ! Begin creating individual merger tree datasets if requested.
    if (mergerTrees%doMakeReferences) then
       
       ! Create a containing group for individual trees.
       treesGroup=outputFile%openGroup("mergerTrees","Data for individual merger trees.")
       
       ! Create groups for trees and dataset references.
       do iTree=1,mergerTrees%treeCount
          groupName="mergerTree"
          groupName=groupName//iTree
          treeGroup=treesGroup%openGroup(char(groupName),"Data for a merger tree.")
          
          ! Standard datasets.
          hyperslabStart(1)=mergerTrees%treeBeginsAt (iTree)
          hyperslabCount(1)=mergerTrees%treeNodeCount(iTree)
          do iProperty=3,size(propertyNames)
             ! Skip cases where we have the corresponding 3-D dataset.
             if (trim(propertyNames(iProperty)) == "spin"            .and. .not.mergerTrees%hasSpinMagnitude           ) cycle
             if (trim(propertyNames(iProperty)) == "angularMomentum" .and. .not.mergerTrees%hasAngularMomentumMagnitude) cycle
             if (haloTrees%hasDataset(trim(propertyNames(iProperty)))) then
                treeDataset=haloTrees%openDataset(propertyNames(iProperty))
                call treeGroup%createReference1D(treeDataset,trim(propertyNames(iProperty)),hyperslabStart,hyperslabCount)
                call treeDataset%close()
             end if
          end do
          
          ! Datasets for properties defined in 3-D spaces.
          hyperslabStart(2)=hyperslabStart(1)
          hyperslabCount(2)=hyperslabCount(1)
          hyperslabStart(1)=int(1,kind=HSIZE_T)
          hyperslabCount(1)=int(3,kind=HSIZE_T)
          do iProperty=1,size(propertyNames3D)
             if (haloTrees%hasDataset(trim(propertyNames3D(iProperty)))) then
                treeDataset=haloTrees%openDataset(propertyNames3D(iProperty))
                call treeGroup%createReference2D(treeDataset,trim(propertyNames3D(iProperty)),hyperslabStart,hyperslabCount)
                call treeDataset%close()
             end if
          end do
          
          call treeGroup%close()
       end do
    
       ! Close the merger trees group.
       call treesGroup%close()
       
       ! Finished making individual merger tree datasets
    end if

    ! Write particle data if necessary.
    if (mergerTrees%hasParticles) then
       ! Open the particles group.
       particlesGroup=outputFile%openGroup("particles","Data for a particles.")

       ! Write datasets.
       if (mergerTrees%hasParticleIndex    ) call particlesGroup%writeDataset(mergerTrees%particleIndex   ,"particleID","The ID of each particle."      ,appendTo=.true.)
       if (mergerTrees%hasParticleRedshift ) call particlesGroup%writeDataset(mergerTrees%particleRedshift,"redshift"  ,"The redshift of each particle.",appendTo=.true.)
       if (mergerTrees%hasParticlePositionX) call particlesGroup%writeDataset(mergerTrees%particlePosition,"position"  ,"The position of each particle.",appendTo=.true.)
       if (mergerTrees%hasParticleVelocityX) call particlesGroup%writeDataset(mergerTrees%particleVelocity,"velocity"  ,"The velocity of each particle.",appendTo=.true.)

       ! Close the particles group.
       call particlesGroup%close()
    end if

    ! Create datasets giving positions of merger trees within the node arrays.
    treeIndexGroup=outputFile%openGroup("treeIndex","Locations of merger trees within the halo data arrays.")
    if (appendActual) then
       call treeIndexGroup%readDataset("firstNode"    ,firstNode    )
       call treeIndexGroup%readDataset("numberOfNodes",numberOfNodes)
       mergerTrees%treeBeginsAt=mergerTrees%treeBeginsAt+firstNode(size(firstNode))+numberOfNodes(size(numberOfNodes))
       call Dealloc_Array(firstNode    )
       call Dealloc_Array(numberOfNodes)
    end if
    call treeIndexGroup%writeDataset(mergerTrees%treeBeginsAt ,"firstNode"    ,"Position of the first node in each tree in the halo data arrays.",appendTo=.true.)
    call treeIndexGroup%writeDataset(mergerTrees%treeNodeCount,"numberOfNodes","Number of nodes in each tree."                                   ,appendTo=.true.)
    call treeIndexGroup%writeDataset(mergerTrees%treeID       ,"treeIndex"    ,"Unique index of tree."                                           ,appendTo=.true.)
    call treeIndexGroup%close()
       
    ! Only write remaining data if we are not appending to an existing file.
    if (.not.appendActual) then
       
       ! Set tree metadata.
       ! Determine if trees have subhalos.
       if (mergerTrees%hasHostIndex) then
          ! A host index is included. If any node is not its own host, then it's a subhalo.
          if (any(mergerTrees%nodeIndex /= mergerTrees%hostIndex)) then
             integerAttribute=1
          else
             integerAttribute=0
          end if
       else
          ! No host index is included - assume no nodes are subhalos.
          integerAttribute=0
       end if
       call haloTrees%writeAttribute(integerAttribute,"treesHaveSubhalos")
       ! Determine if trees are self-contained.
       if (mergerTrees%areSelfContained) then
          integerAttribute=1
       else
          integerAttribute=0
       end if
       call haloTrees%writeAttribute(integerAttribute,"treesAreSelfContained")
       ! Determine if velocities include the Hubble flow.
       if (mergerTrees%includesHubbleFlow) then
          integerAttribute=1
       else
          integerAttribute=0
       end if
       call haloTrees%writeAttribute(integerAttribute,"velocitiesIncludeHubbleFlow")
       ! Determine if positions are periodic.
       if (mergerTrees%isPeriodic) then
          integerAttribute=1
       else
          integerAttribute=0
       end if
       call haloTrees%writeAttribute(integerAttribute,"positionsArePeriodic")
       ! Determine if halo masses include subhalo contributions.
       if (mergerTrees%includesSubhaloMasses) then
          integerAttribute=1
       else
          integerAttribute=0
       end if
       call haloTrees%writeAttribute(integerAttribute,"haloMassesIncludeSubhalos")
       
       ! Store units.
       unitsGroup=outputFile%openGroup("units","The units system used.")
       if (mergerTrees%unitsSet(unitsMass    )) call Store_Unit_Attributes_Galacticus(unitsMass    ,"mass"    ,mergerTrees,unitsGroup)
       if (mergerTrees%unitsSet(unitsLength  )) call Store_Unit_Attributes_Galacticus(unitsLength  ,"length"  ,mergerTrees,unitsGroup)
       if (mergerTrees%unitsSet(unitsTime    )) call Store_Unit_Attributes_Galacticus(unitsTime    ,"time"    ,mergerTrees,unitsGroup)
       if (mergerTrees%unitsSet(unitsVelocity)) call Store_Unit_Attributes_Galacticus(unitsVelocity,"velocity",mergerTrees,unitsGroup)
       call unitsGroup%close()
 
       ! Create groups for attributes.
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataGeneric    )) genericGroup    =outputFile%openGroup("metaData"   ,"Generic metadata."                  )
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataCosmology  )) cosmologyGroup  =outputFile%openGroup("cosmology"  ,"Cosmological parameters."           )
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataSimulation )) simulationGroup =outputFile%openGroup("simulation" ,"Simulation parameters."             )
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataGroupFinder)) groupFinderGroup=outputFile%openGroup("groupFinder","Group finder parameters."           )
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTreeBuilder)) treeBuilderGroup=outputFile%openGroup("treeBuilder","Tree building algorithm parameters.")
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataProvenance )) provenanceGroup =outputFile%openGroup("provenance" ,"Data provenance."                   )
       
       ! Write attributes.
       do iAttribute=1,mergerTrees%metaDataCount
          
          ! Determine which group to write to.
          select case (mergerTrees%metaData(iAttribute)%metadataType)
          case (metaDataGeneric    )
             attributeGroup => genericGroup
          case (metaDataCosmology  )
             attributeGroup => cosmologyGroup
          case (metaDataSimulation )
             attributeGroup => simulationGroup
          case (metaDataGroupFinder)
             attributeGroup => groupFinderGroup
          case (metaDataTreeBuilder)
             attributeGroup => treeBuilderGroup
          case (metaDataProvenance )
             attributeGroup => provenanceGroup
          end select
          
          ! Determine what data type to write.
          select case (mergerTrees%metaData(iAttribute)%dataType)
          case (dataTypeInteger)
             call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%integerAttribute,char(mergerTrees%metaData(iAttribute)%label))
          case (dataTypeDouble )
             call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%doubleAttribute ,char(mergerTrees%metaData(iAttribute)%label))
          case (dataTypeText   )
             call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%textAttribute   ,char(mergerTrees%metaData(iAttribute)%label))
          end select
          
       end do
       
       ! Close attribute groups.
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataGeneric    )) call genericGroup    %close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataCosmology  )) call cosmologyGroup  %close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataSimulation )) call simulationGroup %close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataGroupFinder)) call groupFinderGroup%close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTreeBuilder)) call treeBuilderGroup%close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataProvenance )) call provenanceGroup %close()

    end if

    ! Close the group for datasets.
    call haloTrees%close()

    ! Close the output file.
    call outputFile%close()
    !$omp end critical (HDF5_Access)

    return
  end subroutine Merger_Tree_Data_Structure_Export_Galacticus

  subroutine Merger_Tree_Data_Structure_Export_IRATE(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    !% Output a set of merger trees to an IRATE-format HDF5 file.
    use HDF5
    use IO_HDF5
    use Galacticus_Error
    use String_Handling
    use Memory_Management
    use Array_Utilities
    implicit none
    integer,                 intent(in)                :: hdfChunkSize,hdfCompressionLevel
    class(mergerTreeData),   intent(inout)             :: mergerTrees
    character(len=*),        intent(in)                :: outputFileName
    logical,                 intent(in), optional      :: append
    type(hdf5Object),        pointer                   :: attributeGroup
    type(hdf5Object),        target                    :: outputFile,haloTrees,cosmologyGroup,simulationGroup,snapshotGroup&
         &,thisDataset,mergerTreesGroup,particlesGroup,darkParticlesGroup
    integer,                 dimension(:), allocatable :: thisSnapshotIndices,nodeSnapshotIndices
    integer(kind=kind_int8), dimension(:), allocatable :: descendentSnapshot
    double precision,        dimension(:), allocatable :: particleMass
    integer                                            :: iAttribute,nodesOnSnapshotCount,particlesOnSnapshotCount
    integer(kind=kind_int8)                            :: snapshotMinimum,snapshotMaximum,iSnapshot,iNode,iDescendent
    character(len=14)                                  :: snapshotGroupName
    logical                                            :: appendActual

    ! Determine if we are to append to an existing file.
    appendActual=.false.
    if (present(append)) appendActual=append

    ! Validate the merger tree.
    call Merger_Tree_Data_Validate_Trees            (mergerTrees)

    ! If we have the particle mass, set the masses of any subhalos (which have zero mass by default) based on particle count.
    call Merger_Tree_Data_Set_Subhalo_Masses        (mergerTrees)

    ! If we have most-bound particle indices and particle data has been read, construct arrays giving position of particle data for each node.
    call Merger_Tree_Data_Construct_Particle_Indices(mergerTrees)

    ! IRATE-specific validation.
    if (.not.mergerTrees%hasSnapshot       ) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Export_IRATE','snapshot indices are required for this format'  )
    if (.not.mergerTrees%hasPositionX      ) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Export_IRATE','halo positions are required for this format'    )
    if (.not.mergerTrees%hasNodeIndex      ) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Export_IRATE','halo indices are required for this format'      )
    if (.not.mergerTrees%hasDescendentIndex) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Export_IRATE','descendent indices are required for this format')

    ! Open the output file.
    !$omp critical (HDF5_Access)
    call outputFile%openFile(outputFileName,overWrite=.not.appendActual,chunkSize=hdfChunkSize,compressionLevel=hdfCompressionLevel)

    ! Write the IRATE version.
    if (.not.appendActual) call outputFile%writeAttribute(0,"IRATEVersion")

    ! Find the highest and lowest snapshot numbers in the trees.
    snapshotMinimum=minval(mergerTrees%snapshot)
    snapshotMaximum=maxval(mergerTrees%snapshot)

    ! Loop over snapshots to output.
    do iSnapshot=snapshotMinimum,snapshotMaximum
       
       ! Create a snapshot group.
       write (snapshotGroupName,'("Snapshot",i5.5)') iSnapshot
       snapshotGroup=outputFile%openGroup(trim(snapshotGroupName),"Stores all data for a snapshot.")
       
       ! Create a group for halo catalogs.
       haloTrees=snapshotGroup%openGroup("HaloCatalog","Stores all data for halo catalogs.")

       ! Find those nodes which exist at this snapshot.
       nodesOnSnapshotCount=count(mergerTrees%snapshot == iSnapshot)
       call Alloc_Array(thisSnapshotIndices,[nodesOnSnapshotCount])
       call Array_Which(mergerTrees%snapshot == iSnapshot,thisSnapshotIndices)

       ! Write redshift attribute.
       if (.not.appendActual) call snapshotGroup%writeAttribute(mergerTrees%redshift(thisSnapshotIndices(1)),"Redshift")

       ! Write the data.
       if (mergerTrees%hasNodeIndex               ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeIndex               ,thisSnapshotIndices),"Index"          ,"The index of each halo."                                            ,appendTo=.true.)
       end if
       if (mergerTrees%hasNodeMass                ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeMass                ,thisSnapshotIndices),"Mass"           ,"The mass of each halo."                 ,datasetReturned=thisDataset,appendTo=.true.)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasPositionX               ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%position                ,thisSnapshotIndices,indexOn=1),"Center"         ,"The position of each halo center."      ,datasetReturned=thisDataset,appendTo=.true.)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([          unitsLength              ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasVelocityX               ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%velocity                ,thisSnapshotIndices,indexOn=1),"Velocity"       ,"The velocity of each halo."             ,datasetReturned=thisDataset,appendTo=.true.)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsVelocity                      ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasSpinX                   ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%spin                    ,thisSnapshotIndices),"Spin"           ,"The spin of each halo."                                             ,appendTo=.true.)
       end if
       if (mergerTrees%hasAngularMomentumX        ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%angularMomentum         ,thisSnapshotIndices),"AngularMomentum","The angular momentum spin of each halo.",datasetReturned=thisDataset,appendTo=.true.)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass,unitsLength,unitsVelocity],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasSpinMagnitude           ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%spinMagnitude           ,thisSnapshotIndices),"Spin"           ,"The spin of each halo."                                             ,appendTo=.true.)
       end if
       if (mergerTrees%hasAngularMomentumMagnitude) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%angularMomentumMagnitude,thisSnapshotIndices),"AngularMomentum","The angular momentum spin of each halo.",datasetReturned=thisDataset,appendTo=.true.)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass,unitsLength,unitsVelocity],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasHalfMassRadius          ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%halfMassRadius          ,thisSnapshotIndices),"HalfMassRadius" ,"The half mass radius of each halo."     ,datasetReturned=thisDataset,appendTo=.true.)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       
       ! Destroy snapshot indices.
       call Dealloc_Array(thisSnapshotIndices)
       
       ! Close the group for halo catalogs
       call haloTrees%close()
       
       ! Close the snapshot group.
       call snapshotGroup%close()
       
    end do

    ! Create a group for merger trees.
    mergerTreesGroup=outputFile%openGroup("MergerTrees","Stores all data for merger trees.")

    ! Specify the name of the halo catalog group.
    if (.not.appendActual) call mergerTreesGroup%writeAttribute("HaloCatalog","HaloCatalogName")
    
    ! Build snapshot numbers for descendents.
    call Alloc_Array(descendentSnapshot,shape(mergerTrees%nodeIndex))
    descendentSnapshot=-1
    do iDescendent=1,size(mergerTrees%nodeIndex)
       if (mergerTrees%descendentIndex(iDescendent) >= 0) then
          iNode=0
          do while (iNode < size(mergerTrees%nodeIndex) .and. descendentSnapshot(iDescendent) < 0)
             iNode=iNode+1
             if (mergerTrees%nodeIndex(iNode) == mergerTrees%descendentIndex(iDescendent))&
                  & descendentSnapshot(iDescendent)=mergerTrees%snapshot(iNode)
          end do
       end if
    end do

    ! Output merger tree datasets.
    call                               mergerTreesGroup%writeDataset(mergerTrees%snapshot          ,"HaloSnapshot"      ,"The snapshot of each halo."           ,appendTo=.true.)
    call                               mergerTreesGroup%writeDataset(mergerTrees%nodeIndex         ,"HaloID"         ,"The index of each halo."              ,appendTo=.true.)
    call                               mergerTreesGroup%writeDataset(mergerTrees%descendentIndex   ,"DescendentID"   ,"The index of each descendent halo."   ,appendTo=.true.)
    call                               mergerTreesGroup%writeDataset(            descendentSnapshot,"DescendentSnapshot","The snapshot of each descendent halo.",appendTo=.true.)
    if (mergerTrees%hasHostIndex) call mergerTreesGroup%writeDataset(mergerTrees%hostIndex         ,"HostID"         ,"The index of each host halo."         ,appendTo=.true.)
    call                               mergerTreesGroup%writeDataset(mergerTrees%treeNodeCount     ,"HalosPerTree"      ,"Number of halos in each tree."        ,appendTo=.true.)
    call                               mergerTreesGroup%writeDataset(mergerTrees%treeID            ,"TreeID"         ,"Unique index of tree."                ,appendTo=.true.)
    call Dealloc_Array(descendentSnapshot)
    call mergerTreesGroup%close()

    if (mergerTrees%hasMostBoundParticleIndex) then
       ! Find the highest and lowest snapshot numbers in the particles.
       if (.not.mergerTrees%hasParticleSnapshot) call Galacticus_Error_Report('Merger_Tree_Data_Structure_Export_IRATE','particle snapshot numbers must be available for IRATE format export')
       snapshotMinimum=minval(mergerTrees%particleSnapshot)
       snapshotMaximum=maxval(mergerTrees%particleSnapshot)

       ! Loop over snapshots to output.
       do iSnapshot=snapshotMinimum,snapshotMaximum

          ! Find those particles which exist at this snapshot.
          particlesOnSnapshotCount=count(mergerTrees%particleSnapshot == iSnapshot)
          call Alloc_Array(thisSnapshotIndices,[particlesOnSnapshotCount])
          call Array_Which(mergerTrees%particleSnapshot == iSnapshot,thisSnapshotIndices)
          
          ! Find those nodes which exist at this snapshot.
          nodesOnSnapshotCount=count(mergerTrees%snapshot == iSnapshot)
          call Alloc_Array(nodeSnapshotIndices,[nodesOnSnapshotCount])
          call Array_Which(mergerTrees%snapshot == iSnapshot,nodeSnapshotIndices)
          
          ! Create a snapshot group.
          write (snapshotGroupName,'("Snapshot",i5.5)') iSnapshot
          snapshotGroup=outputFile%openGroup(trim(snapshotGroupName),"Stores all data for a snapshot.")

          ! Create a group for halo catalogs.
          haloTrees=snapshotGroup%openGroup("HaloCatalog","Stores all data for halo catalogs.")

          ! Write the particle indices.
          call haloTrees%writeDataset(Array_Index(mergerTrees%mostBoundParticleIndex,nodeSnapshotIndices),"MostBoundParticleID","The index of each particle.",appendTo=.true.)
          call haloTrees%close()

          ! Create a group for particles.
          particlesGroup=snapshotGroup%openGroup("ParticleData","Stores all data for particles.")

          ! Make a group for dark matter particles.
          darkParticlesGroup=particlesGroup%openGroup("Dark","Stores all data for dark matter particles.")

          ! Write redshift attribute.
          if (.not.snapshotGroup%hasAttribute("Redshift")) call snapshotGroup%writeAttribute(mergerTrees%particleRedshift(thisSnapshotIndices(1)),"Redshift")

          ! Write the data.
          call Alloc_Array(particleMass,[particlesOnSnapshotCount])
          particleMass=mergerTrees%particleMass
          call darkParticlesGroup%writeDataset(particleMass,"Mass","The mass of each particle.",datasetReturned=thisDataset,appendTo=.true.)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass],mergerTrees,thisDataset)
          call thisDataset%close()
          call Dealloc_Array(particleMass)
          if (mergerTrees%hasParticleIndex    ) then
             call darkParticlesGroup%writeDataset(Array_Index(mergerTrees%particleIndex   ,thisSnapshotIndices),"ID"      ,"The index of each particle."                              ,appendTo=.true.)
          end if
          if (mergerTrees%hasParticlePositionX) then
             call darkParticlesGroup%writeDataset(Array_Index(mergerTrees%particlePosition,thisSnapshotIndices),"Position","The position of each particle.",datasetReturned=thisDataset,appendTo=.true.)
             if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsLength  ],mergerTrees,thisDataset)
             call thisDataset%close()
          end if
          if (mergerTrees%hasParticleVelocityX) then
             call darkParticlesGroup%writeDataset(Array_Index(mergerTrees%particleVelocity,thisSnapshotIndices),"Velocity","The velocity of each particle.",datasetReturned=thisDataset,appendTo=.true.)
             if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsVelocity],mergerTrees,thisDataset)
             call thisDataset%close()
          end if
          ! Destroy the snapshot indices.
          call Dealloc_Array(thisSnapshotIndices)
          call Dealloc_Array(nodeSnapshotIndices)

          ! Close the groups.
          call darkParticlesGroup%close()
          call particlesGroup    %close()
          call snapshotGroup     %close()

       end do

    end if

    ! Create groups for attributes.
    if (.not.appendActual) then
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataCosmology  )) cosmologyGroup  =outputFile%openGroup("Cosmology"            ,"Cosmological parameters."           )
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataSimulation )) simulationGroup =outputFile%openGroup("SimulationProperties" ,"Simulation parameters."             )
       
       ! Write attributes.
       do iAttribute=1,mergerTrees%metaDataCount
          
          ! Determine which group to write to.
          attributeGroup => null()
          select case (mergerTrees%metaData(iAttribute)%metadataType)
          case (metaDataCosmology  )
             attributeGroup => cosmologyGroup
          case (metaDataSimulation )
             attributeGroup => simulationGroup
          end select
          
          ! Check if the group was recognized.
          if (associated(attributeGroup)) then
             
             ! Perform dictionary mapping from our internal names (which follow Galacticus format) to IRATE names.
             select case (mergerTrees%metaData(iAttribute)%metadataType)
             case (metaDataCosmology  )
                select case (char(mergerTrees%metaData(iAttribute)%label))
                case ("powerSpectrumIndex")
                   mergerTrees%metaData(iAttribute)%label="PowerSpectrumIndex"
                end select
             end select
             
             ! Determine what data type to write.
             select case (mergerTrees%metaData(iAttribute)%dataType)
             case (dataTypeInteger)
                call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%integerAttribute,char(mergerTrees%metaData(iAttribute)%label))
             case (dataTypeDouble )
                call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%doubleAttribute ,char(mergerTrees%metaData(iAttribute)%label))
             case (dataTypeText   )
                call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%textAttribute   ,char(mergerTrees%metaData(iAttribute)%label))
             end select
             
          end if
          
       end do
       
       ! Close attribute groups.
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataCosmology  )) call cosmologyGroup  %close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataSimulation )) call simulationGroup %close()
    end if

    ! Close the output file.
    call outputFile%close()
    !$omp end critical (HDF5_Access)

    return
  end subroutine Merger_Tree_Data_Structure_Export_IRATE

  subroutine Store_Unit_Attributes_Galacticus(unitType,unitLabel,mergerTrees,unitsGroup)
    !% Store attributes describing the unit system.
    use IO_HDF5
    implicit none
    integer,               intent(in)    :: unitType
    character(len=*),      intent(in)    :: unitLabel
    class(mergerTreeData), intent(in)    :: mergerTrees
    type(hdf5Object),      intent(inout) :: unitsGroup

    call unitsGroup%writeAttribute(mergerTrees%units(unitType)%unitsInSI          ,unitLabel//"UnitsInSI"          )
    call unitsGroup%writeAttribute(mergerTrees%units(unitType)%hubbleExponent     ,unitLabel//"HubbleExponent"     )
    call unitsGroup%writeAttribute(mergerTrees%units(unitType)%scaleFactorExponent,unitLabel//"ScaleFactorExponent")
    return
  end subroutine Store_Unit_Attributes_Galacticus

  subroutine Store_Unit_Attributes_IRATE(unitType,mergerTrees,thisDataset)
    !% Store unit attributes in IRATE format files.
    use IO_HDF5
    use Numerical_Constants_Prefixes
    implicit none
    integer,               intent(in),    dimension(:) :: unitType
    class(mergerTreeData), intent(in)                  :: mergerTrees
    type(hdf5Object),      intent(inout)               :: thisDataset
    integer                                            :: iUnit
    double precision                                   :: cgsUnits,hubbleExponent,scaleFactorExponent
    type(varying_string)                               :: unitName

    ! Get conversion factor to cgs units.
    cgsUnits           =1.0d0
    hubbleExponent     =0.0d0
    scaleFactorExponent=0.0d0
    unitName           =""
    do iUnit=1,size(unitType)
       cgsUnits=cgsUnits*mergerTrees%units(unitType(iUnit))%unitsInSI
       select case (unitType(iUnit))
       case (unitsMass    )
          cgsUnits=cgsUnits*kilo 
       case (unitsLength  )
          cgsUnits=cgsUnits*hecto
       case (unitsTime    )
          cgsUnits=cgsUnits*1.0d0
       case (unitsVelocity)
          cgsUnits=cgsUnits*hecto
       end select
       hubbleExponent     =hubbleExponent     +dble(mergerTrees%units(unitType(iUnit))%hubbleExponent     )
       scaleFactorExponent=scaleFactorExponent+dble(mergerTrees%units(unitType(iUnit))%scaleFactorExponent)
       if (iUnit > 1) unitName=unitName//" * "
       unitName=unitName//mergerTrees%units(unitType(iUnit))%name
    end do
    
    ! Write the attributes.
    call thisDataset%writeAttribute(                      &
         &                          [                     &
         &                           cgsUnits           , &
         &                           hubbleExponent     , &
         &                           scaleFactorExponent  &
         &                          ],                    &
         &                          "unitscgs"            &
         &                         ) 
    call thisDataset%writeAttribute(                      &
         &                           unitName,            &
         &                          "unitname"            &
         &                         )
    return
  end subroutine Store_Unit_Attributes_IRATE
  
  subroutine Merger_Tree_Data_Validate_Trees(mergerTrees)
    !% Validate the merger trees.
    use Galacticus_Error
    class(mergerTreeData), intent(in) :: mergerTrees
    
    if (.not.mergerTrees%hasTreeIndex      ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'treeIndex'"      )
    if (.not.mergerTrees%hasNodeIndex      ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'nodeIndex'"      )
    if (.not.mergerTrees%hasDescendentIndex) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'descendentIndex'")
    if (.not.mergerTrees%hasRedshift       ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'redshift'"       )
    if (.not.mergerTrees%hasNodeMass       ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'nodeMass'"       )
    return
  end subroutine Merger_Tree_Data_Validate_Trees

  subroutine Merger_Tree_Data_Set_Subhalo_Masses(mergerTrees)
    !% Set the masses of any subhalos (which have zero mass by default) based on particle count.
    class(mergerTreeData), intent(inout) :: mergerTrees
    
    if (mergerTrees%hasParticleCount) then
       where (mergerTrees%nodeMass <= 0.0d0)
          mergerTrees%nodeMass = dble(mergerTrees%particleCount)*mergerTrees%particleMass
       end where
    end if
    return
  end subroutine Merger_Tree_Data_Set_Subhalo_Masses

  subroutine Merger_Tree_Data_Construct_Particle_Indices(mergerTrees)
    !% If we have most-bound particle indices and particle data has been read, construct arrays giving position of particle data for each node.
    use Galacticus_Error
    use Memory_Management
    class(mergerTreeData), intent(inout) :: mergerTrees
    logical                              :: foundParticleData
    integer                              :: iNode,iParticle

    if (mergerTrees%hasMostBoundParticleIndex) then
       ! Insist on having particle data.
       if (.not.mergerTrees%hasParticles) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","most bound particle IDs provided, but no particle data was read")
       ! Allocate arrays for storing indices.
       if (allocated(mergerTrees%particleReferenceStart)) call Dealloc_Array(mergerTrees%particleReferenceStart)
       if (allocated(mergerTrees%particleReferenceCount)) call Dealloc_Array(mergerTrees%particleReferenceCount)
       call Alloc_Array(mergerTrees%particleReferenceStart,[mergerTrees%nodeCount])
       call Alloc_Array(mergerTrees%particleReferenceCount,[mergerTrees%nodeCount])
       mergerTrees%particleReferenceStart=-1
       mergerTrees%particleReferenceCount=-1
       ! Loop over nodes.
       do iNode=1,mergerTrees%nodeCount
          ! Flag that particle data has not yet been found.
          foundParticleData=.false.

          ! Loop over all particles.
          do iParticle=1,mergerTrees%particlesCount
             ! Check if this particle has the ID of the node's most bound particle.
             if     (     mergerTrees%particleIndex   (iParticle) == mergerTrees%mostBoundParticleIndex(iNode) &
                  & .and. mergerTrees%particleRedshift(iParticle) <  mergerTrees%redshift(iNode)               &
                  & ) then
                ! It does, so store the position at which the particle data starts (offset by -1 due to HDF5 array convention) if we haven't already done so.
                if (.not.foundParticleData) mergerTrees%particleReferenceStart(iNode)=iParticle-1
                ! Flag that data has been found for this particle.
                foundParticleData=.true.
             else
                ! Particle ID does not match. Have we already found data for this particle? If so, exit.
                if (foundParticleData) exit
             end if
          end do
          if (foundParticleData) mergerTrees%particleReferenceCount(iNode)=min(iParticle-1,mergerTrees%particlesCount)&
               &-mergerTrees%particleReferenceStart(iNode)
       end do
    end if
    return
  end subroutine Merger_Tree_Data_Construct_Particle_Indices
 
end module Merger_Tree_Data_Structure
