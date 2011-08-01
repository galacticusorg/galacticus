!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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

  ! Property names.
  character(len=*), parameter :: propertyNames(25)=[ &
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
       & 'mostBoundParticleIndex'                    &
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
     double precision :: unitsInSI
     integer          :: hubbleExponent,scaleFactorExponent
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
  implicit none
     private
     integer                                              :: treeCount,nodeCount,particlesCount
     double precision                                     :: particleMass=0.0d0
     integer,                 allocatable, dimension(:)   :: columnProperties,particleColumnProperties,treeBeginsAt,treeNodeCount
     integer(kind=kind_int8), allocatable, dimension(:)   :: nodeIndex,treeIndex,descendentIndex,hostIndex,particleCount,treeID&
          &,particleIndex,mostBoundParticleIndex,particleReferenceStart,particleReferenceCount
     double precision,        allocatable, dimension(:)   :: redshift,nodeMass,halfMassRadius,spinMagnitude&
          &,angularMomentumMagnitude,particleRedshift
     double precision,        allocatable, dimension(:,:) :: position,velocity,spin,angularMomentum,particlePosition&
          &,particleVelocity
     logical                                              :: hasNodeIndex,hasTreeIndex,hasDescendentIndex,hasHostIndex &
          &,hasRedshift ,hasNodeMass,hasPositionX,hasPositionY,hasPositionZ,hasVelocityX,hasVelocityY,hasVelocityZ &
          &,hasParticleCount,hasSpinX,hasSpinY,hasSpinZ,hasSpinMagnitude,hasAngularMomentumX,hasAngularMomentumY &
          &,hasAngularMomentumZ,hasAngularMomentumMagnitude,hasHalfMassRadius,hasParticleRedshift,hasParticlePositionX &
          &,hasParticlePositionY,hasParticlePositionZ,hasParticleVelocityX,hasParticleVelocityY,hasParticleVelocityZ&
          &,hasParticleIndex,hasParticles,hasMostBoundParticleIndex
     logical                                              :: areSelfContained=.true., includesHubbleFlow=.false.,&
          & includesSubhaloMasses=.false.,doMakeReferences=.true.
     type(unitsMetaData),     dimension(unitTypeCount)    :: units
     logical,                 dimension(unitTypeCount)    :: unitsSet=.false.
     integer                                              :: metaDataCount=0
     type(treeMetaData),      allocatable, dimension(:)   :: metaData
   contains
     procedure                                            :: nodeCountSet             => Merger_Tree_Data_Structure_Set_Node_Count
     procedure                                            :: particleCountSet         => Merger_Tree_Data_Structure_Set_Particle_Count
     procedure                                            :: readASCII                => Merger_Tree_Data_Structure_Read_ASCII
     procedure                                            :: readParticlesASCII       => Merger_Tree_Data_Structure_Read_Particles_ASCII
     procedure                                            :: setProperty              => Merger_Tree_Data_Structure_Set_Property
     procedure                                            :: setParticleProperty      => Merger_Tree_Data_Structure_Set_Particle_Property
     procedure                                            :: setParticleMass          => Merger_Tree_Data_Structure_Set_Particle_Mass
     procedure                                            :: setSelfContained         => Merger_Tree_Data_Structure_Set_Self_Contained
     procedure                                            :: setIncludesHubbleFlow    => Merger_Tree_Data_Structure_Set_Includes_Hubble_Flow
     procedure                                            :: setIncludesSubhaloMasses => Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses
     procedure                                            :: setUnits                 => Merger_Tree_Data_Structure_Set_Units
     procedure                                            ::                             Merger_Tree_Data_Structure_Add_Metadata_Double
     procedure                                            ::                             Merger_Tree_Data_Structure_Add_Metadata_Integer
     procedure                                            ::                             Merger_Tree_Data_Structure_Add_Metadata_Text
     generic                                              :: addMetadata              => Merger_Tree_Data_Structure_Add_Metadata_Double , &
          &                                                                              Merger_Tree_Data_Structure_Add_Metadata_Integer, &
          &                                                                              Merger_Tree_Data_Structure_Add_Metadata_Text
     procedure                                            :: makeReferences           => Merger_Tree_Data_Structure_Make_References
     procedure                                            :: export                   => Merger_Tree_Data_Structure_Export
  end type mergerTreeData

contains

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

  subroutine Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses(mergerTrees,includesSubhaloMasses)
    !% Set the particle mass used in the trees.
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    logical,               intent(in)    :: includesSubhaloMasses

    ! Set whether halo masses include subhalo contributions.
    mergerTrees%includesSubhaloMasses=includesSubhaloMasses

    return
  end subroutine Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses

  subroutine Merger_Tree_Data_Structure_Set_Units(mergerTrees,unitType,unitsInSI,hubbleExponent,scaleFactorExponent)
    !% Set the units system.
    use Galacticus_Error
    implicit none
    class(mergerTreeData), intent(inout)        :: mergerTrees
    integer,               intent(in)           :: unitType
    double precision,      intent(in)           :: unitsInSI
    integer,               intent(in), optional :: hubbleExponent,scaleFactorExponent

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
    return
  end subroutine Merger_Tree_Data_Structure_Set_Units

  subroutine Merger_Tree_Data_Structure_Set_Particle_Property(mergerTrees,propertyType,columnNumber)
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
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Property

  subroutine Merger_Tree_Data_Structure_Set_Property(mergerTrees,propertyType,columnNumber)
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
  end subroutine Merger_Tree_Data_Structure_Set_Property

  subroutine Merger_Tree_Data_Structure_Read_ASCII(mergerTrees,inputFile,lineNumberStart,lineNumberStop,separator)
    !% Read in merger tree data from an ASCII file.
    use File_Utilities
    use String_Handling
    use Memory_Management
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    class(mergerTreeData), intent(inout)              :: mergerTrees
    character(len=*),      intent(in)                 :: inputFile
    integer,               intent(in),   optional     :: lineNumberStart,lineNumberStop
    character(len=*),      intent(in),   optional     :: separator
    character(len=32),     allocatable,  dimension(:) :: inputColumns
    integer                                           :: lineNumberStartActual,lineNumberStopActual,columnsCount,lineNumber&
         &,fileUnit,iColumn,iNode,iTree
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
    fileUnit=File_Units_Get()
    open(fileUnit,file=inputFile,status='old',form='formatted')
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

    ! Allocate arrays for tree start and stop indices and reference ID.
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
  end subroutine Merger_Tree_Data_Structure_Read_ASCII

  subroutine Merger_Tree_Data_Structure_Read_Particles_ASCII(mergerTrees,inputFile,lineNumberStart,lineNumberStop,separator)
    !% Read in particle data from an ASCII file.
    use File_Utilities
    use String_Handling
    use Memory_Management
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    class(mergerTreeData), intent(inout)              :: mergerTrees
    character(len=*),      intent(in)                 :: inputFile
    integer,               intent(in),   optional     :: lineNumberStart,lineNumberStop
    character(len=*),      intent(in),   optional     :: separator
    character(len=32),     allocatable,  dimension(:) :: inputColumns
    integer                                           :: lineNumberStartActual,lineNumberStopActual,columnsCount,lineNumber&
         &,fileUnit,iColumn,iNode,iTree
    logical                                           :: gotFirstDataLine
    character(len=1024)                               :: inputLine
    type(varying_string)                              :: message

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
 
    ! Allocate internal arrays to correct size as needed.
    if (mergerTrees%hasParticleIndex    ) call Alloc_Array(mergerTrees%particleIndex   ,[  mergerTrees%particlesCount])
    if (mergerTrees%hasParticleRedshift ) call Alloc_Array(mergerTrees%particleRedshift,[  mergerTrees%particlesCount])
    if (mergerTrees%hasParticlePositionX) call Alloc_Array(mergerTrees%particlePosition,[3,mergerTrees%particlesCount])
    if (mergerTrees%hasParticleVelocityX) call Alloc_Array(mergerTrees%particleVelocity,[3,mergerTrees%particlesCount])

    ! Open the file and read lines.
    fileUnit=File_Units_Get()
    open(fileUnit,file=inputFile,status='old',form='formatted')
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

  subroutine Merger_Tree_Data_Structure_Export(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel)
    !% Output a set of merger trees to an HDF5 file.
    use HDF5
    use IO_HDF5
    use Galacticus_Error
    use String_Handling
    use Memory_Management
    implicit none
    integer,                 intent(in)    :: hdfChunkSize,hdfCompressionLevel
    class(mergerTreeData),   intent(inout) :: mergerTrees
    character(len=*),        intent(in)    :: outputFileName
    integer(kind=HSIZE_T) ,  dimension(2)  :: hyperslabStart,hyperslabCount
    type(hdf5Object),        pointer       :: attributeGroup
    type(hdf5Object),        target        :: outputFile,haloTrees,treesGroup,treeGroup,treeDataset,treeIndexGroup,unitsGroup&
         &,genericGroup,cosmologyGroup,simulationGroup,groupFinderGroup,treeBuilderGroup,provenanceGroup,particlesGroup
    integer                                :: iTree,iProperty,integerAttribute,iAttribute,iNode,iParticle
    logical                                :: foundParticleData
    type(varying_string)                   :: groupName

    ! Validate the merger tree.
    if (.not.mergerTrees%hasTreeIndex      ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'treeIndex'"      )
    if (.not.mergerTrees%hasNodeIndex      ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'nodeIndex'"      )
    if (.not.mergerTrees%hasDescendentIndex) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'descendentIndex'")
    if (.not.mergerTrees%hasRedshift       ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'redshift'"       )
    if (.not.mergerTrees%hasNodeMass       ) call Galacticus_Error_Report("Merger_Tree_Data_Structure_Export","merger trees do not have required property 'nodeMass'"       )

    ! If we have the particle mass, set the masses of any subhalos (which have zero mass by default) based on particle count.
    if (mergerTrees%hasParticleCount) then
       where (mergerTrees%nodeMass <= 0.0d0)
          mergerTrees%nodeMass = dble(mergerTrees%particleCount)*mergerTrees%particleMass
       end where
    end if

    ! If we have most-bound particle indices and particle data has been read, construct arrays giving position of particle data for each node.
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

    ! Open the output file.
    call outputFile%openFile(outputFileName,overWrite=.true.,chunkSize=hdfChunkSize,compressionLevel=hdfCompressionLevel)

    ! Create a group for the datasets.
    haloTrees=IO_HDF5_Open_Group(outputFile,"haloTrees","Stores all data for merger trees.")

    ! Write the data.
    if (mergerTrees%hasNodeIndex               ) call haloTrees%writeDataset(mergerTrees%nodeIndex               ,"nodeIndex"      ,"The index of each node."                )
    if (mergerTrees%hasDescendentIndex         ) call haloTrees%writeDataset(mergerTrees%descendentIndex         ,"descendentIndex","The index of each descendent node."     )
    if (mergerTrees%hasHostIndex               ) call haloTrees%writeDataset(mergerTrees%hostIndex               ,"hostIndex"      ,"The index of each host node."           )
    if (mergerTrees%hasNodeMass                ) call haloTrees%writeDataset(mergerTrees%nodeMass                ,"nodeMass"       ,"The mass of each node."                 )
    if (mergerTrees%hasRedshift                ) call haloTrees%writeDataset(mergerTrees%redshift                ,"redshift"       ,"The redshift of each node."             )
    if (mergerTrees%hasPositionX               ) call haloTrees%writeDataset(mergerTrees%position                ,"position"       ,"The position of each node."             )
    if (mergerTrees%hasVelocityX               ) call haloTrees%writeDataset(mergerTrees%velocity                ,"velocity"       ,"The velocity of each node."             )
    if (mergerTrees%hasSpinX                   ) call haloTrees%writeDataset(mergerTrees%spin                    ,"spin"           ,"The spin of each node."                 )
    if (mergerTrees%hasAngularMomentumX        ) call haloTrees%writeDataset(mergerTrees%angularMomentum         ,"angularMomentum","The angular momentum spin of each node.")
    if (mergerTrees%hasSpinMagnitude           ) call haloTrees%writeDataset(mergerTrees%spinMagnitude           ,"spin"           ,"The spin of each node."                 )
    if (mergerTrees%hasAngularMomentumMagnitude) call haloTrees%writeDataset(mergerTrees%angularMomentumMagnitude,"angularMomentum","The angular momentum spin of each node.")
    if (mergerTrees%hasHalfMassRadius          ) call haloTrees%writeDataset(mergerTrees%halfMassRadius          ,"halfMassRadius" ,"The half mass radius of each node."     )
    if (mergerTrees%hasMostBoundParticleIndex) then
       call haloTrees%writeDataset(mergerTrees%particleReferenceStart,"particleIndexStart","The starting index of particle data for each node.")
       call haloTrees%writeDataset(mergerTrees%particleReferenceCount,"particleIndexCount","The number of particle data for each node."        )
    end if

    ! Begin creating individual merger tree datasets if requested.
    if (mergerTrees%doMakeReferences) then
       
       ! Create a containing group for individual trees.
       treesGroup=IO_HDF5_Open_Group(outputFile,"mergerTrees","Data for individual merger trees.")
       
       ! Create groups for trees and dataset references.
       do iTree=1,mergerTrees%treeCount
          groupName="mergerTree"
          groupName=groupName//iTree
          treeGroup=IO_HDF5_Open_Group(treesGroup,char(groupName),"Data for a merger tree.")
          
          ! Standard datasets.
          hyperslabStart(1)=mergerTrees%treeBeginsAt (iTree)
          hyperslabCount(1)=mergerTrees%treeNodeCount(iTree)
          do iProperty=3,size(propertyNames)
             ! Skip cases where we have the corresponding 3-D dataset.
             if (trim(propertyNames(iProperty)) == "spin"            .and. .not.mergerTrees%hasSpinMagnitude           ) cycle
             if (trim(propertyNames(iProperty)) == "angularMomentum" .and. .not.mergerTrees%hasAngularMomentumMagnitude) cycle
             if (haloTrees%hasDataset(trim(propertyNames(iProperty)))) then
                treeDataset=IO_HDF5_Open_Dataset(haloTrees,propertyNames(iProperty))
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
                treeDataset=IO_HDF5_Open_Dataset(haloTrees,propertyNames3D(iProperty))
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
       particlesGroup=IO_HDF5_Open_Group(outputFile,"particles","Data for a particles.")

       ! Write datasets.
       if (mergerTrees%hasParticleIndex    ) call particlesGroup%writeDataset(mergerTrees%particleIndex   ,"particleID","The ID of each particle."      )
       if (mergerTrees%hasParticleRedshift ) call particlesGroup%writeDataset(mergerTrees%particleRedshift,"redshift"  ,"The redshift of each particle.")
       if (mergerTrees%hasParticlePositionX) call particlesGroup%writeDataset(mergerTrees%particlePosition,"position"  ,"The position of each particle.")
       if (mergerTrees%hasParticleVelocityX) call particlesGroup%writeDataset(mergerTrees%particleVelocity,"velocity"  ,"The velocity of each particle.")

       ! Close the particles group.
       call particlesGroup%close()
    end if

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
    ! Determine if halo masses include subhalo contributions.
    if (mergerTrees%includesSubhaloMasses) then
       integerAttribute=1
    else
       integerAttribute=0
    end if
    call haloTrees%writeAttribute(integerAttribute,"haloMassesIncludeSubhalos")

    ! Close the group for datasets.
    call haloTrees%close()

    ! Store units.
    unitsGroup=IO_HDF5_Open_Group(outputFile,"units","The units system used.")
    if (mergerTrees%unitsSet(unitsMass    )) call Store_Unit_Attributes(unitsMass    ,"mass"    ,mergerTrees,unitsGroup)
    if (mergerTrees%unitsSet(unitsLength  )) call Store_Unit_Attributes(unitsLength  ,"length"  ,mergerTrees,unitsGroup)
    if (mergerTrees%unitsSet(unitsTime    )) call Store_Unit_Attributes(unitsTime    ,"time"    ,mergerTrees,unitsGroup)
    if (mergerTrees%unitsSet(unitsVelocity)) call Store_Unit_Attributes(unitsVelocity,"velocity",mergerTrees,unitsGroup)

    ! Create datasets giving positions of merger trees within the node arrays.
    treeIndexGroup=IO_HDF5_Open_Group(outputFile,"treeIndex","Locations of merger trees within the halo data arrays.")
    call treeIndexGroup%writeDataset(mergerTrees%treeBeginsAt ,"firstNode"    ,"Position of the first node in each tree in the halo data arrays.")
    call treeIndexGroup%writeDataset(mergerTrees%treeNodeCount,"numberOfNodes","Number of nodes in each tree."                                   )
    call treeIndexGroup%writeDataset(mergerTrees%treeID       ,"treeIndex"    ,"Unique index of tree."                                           )
    call treeIndexGroup%close()

    ! Create groups for attributes.
    if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataGeneric    )) genericGroup    =IO_HDF5_Open_Group(outputFile,"metaData"   ,"Generic metadata."                  )
    if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataCosmology  )) cosmologyGroup  =IO_HDF5_Open_Group(outputFile,"cosmology"  ,"Cosmological parameters."           )
    if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataSimulation )) simulationGroup =IO_HDF5_Open_Group(outputFile,"simulation" ,"Simulation parameters."             )
    if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataGroupFinder)) groupFinderGroup=IO_HDF5_Open_Group(outputFile,"groupFinder","Group finder parameters."           )
    if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTreeBuilder)) treeBuilderGroup=IO_HDF5_Open_Group(outputFile,"treeBuilder","Tree building algorithm parameters.")
    if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataProvenance )) provenanceGroup =IO_HDF5_Open_Group(outputFile,"provenance" ,"Data provenance."                   )
    
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

    ! Close the output file.
    call outputFile%close()

    return
  end subroutine Merger_Tree_Data_Structure_Export

  subroutine Store_Unit_Attributes(unitType,unitLabel,mergerTrees,unitsGroup)
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
  end subroutine Store_Unit_Attributes

end module Merger_Tree_Data_Structure
