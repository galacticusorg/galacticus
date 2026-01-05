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

!+ Contributions to this file made by: Stephanie Dörschner.

!!{
Contains a module which implements an object to store merger tree data for processing into \glc's preferred file format.
!!}

module Merger_Tree_Data_Structure
  !!{
  Implements an object to store merger tree data for processing into \glc's preferred file format.
  !!}
  use, intrinsic :: ISO_C_Binding     , only : c_size_t
  use            :: ISO_Varying_String, only : varying_string
  use            :: Kind_Numbers      , only : kind_int8
  implicit none
  private
  public :: mergerTreeData

  ! Output formats.
  !![
  <enumeration>
   <name>mergerTreeFormat</name>
   <description>Used to specify which output format to use for merger tree data.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="galacticus" />
   <entry label="irate"      />
  </enumeration>
  !!]

  ! Property labels.
  !![
  <enumeration>
   <name>propertyType</name>
   <description>Used to specify properties in a {\normalfont \ttfamily mergerTreeData} structure.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <entry label="null"                    />
   <entry label="treeIndex"               />
   <entry label="nodeIndex"               />
   <entry label="descendantIndex"         />
   <entry label="hostIndex"               />
   <entry label="redshift"                />
   <entry label="scaleFactor"             />
   <entry label="nodeMass"                />
   <entry label="nodeMass200Mean"         />
   <entry label="nodeMass200Crit"         />
   <entry label="particleCount"           />
   <entry label="positionX"               />
   <entry label="positionY"               />
   <entry label="positionZ"               />
   <entry label="velocityX"               />
   <entry label="velocityY"               />
   <entry label="velocityZ"               />
   <entry label="spinX"                   />
   <entry label="spinY"                   />
   <entry label="spinZ"                   />
   <entry label="spin"                    />
   <entry label="angularMomentumX"        />
   <entry label="angularMomentumY"        />
   <entry label="angularMomentumZ"        />
   <entry label="angularMomentum"         />
   <entry label="specificAngularMomentumX"/>
   <entry label="specificAngularMomentumY"/>
   <entry label="specificAngularMomentumZ"/>
   <entry label="specificAngularMomentum" />
   <entry label="halfMassRadius"          />
   <entry label="scaleRadius"             />
   <entry label="particleIndex"           />
   <entry label="mostBoundParticleIndex"  />
   <entry label="snapshot"                />
   <entry label="treeWeight"              />
   <entry label="velocityMaximum"         />
   <entry label="velocityDispersion"      />
  </enumeration>
  !!]

  ! Names of 3-D datasets (i.e. those which give properties in 3-D space).
  character(len=*), parameter :: propertyNames3D(5)=[ &
       & 'position               ',                   &
       & 'velocity               ',                   &
       & 'spin                   ',                   &
       & 'angularMomentum        ',                   &
       & 'specificAngularMomentum'                    &
       &                                            ]

  ! Units labels.
  !![
  <enumeration>
   <name>units</name>
   <description>Used to specify the type of units being stored in a {\normalfont \ttfamily mergerTreeData} structure.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <entry label="mass"     />
   <entry label="length"   />
   <entry label="time"     />
   <entry label="velocity" />
  </enumeration>
  !!]

  type unitsMetaData
     !!{
     A structure that holds metadata on units used.
     !!}
     double precision                              :: unitsInSI
     integer                                       :: hubbleExponent, scaleFactorExponent
     type            (varying_string), allocatable :: name
  end type unitsMetaData

  ! Metadata labels.
  !![
  <enumeration>
   <name>metaDataType</name>
   <description>Used to specify the type of metadata being stored in a {\normalfont \ttfamily mergerTreeData} structure.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <entry label="generic"     />
   <entry label="cosmology"   />
   <entry label="simulation"  />
   <entry label="groupFinder" />
   <entry label="treeBuilder" />
   <entry label="provenance"  />
  </enumeration>
  !!]

  ! Data types for metadata.
  !![
  <enumeration>
   <name>dataType</name>
   <description>Used to specify the type of data being stored in a {\normalfont \ttfamily mergerTreeData} structure metadata entry.</description>
   <visibility>public</visibility>
   <entry label="null"    />
   <entry label="integer" />
   <entry label="double"  />
   <entry label="text"    />
  </enumeration>
  !!]

  type treeMetaData
     !!{
     Structure that holds metadata for the trees.
     !!}
     type            (enumerationMetaDataTypeType) :: metadataType
     type            (varying_string             ) :: label
     type            (enumerationDataTypeType    ) :: dataType
     integer                                       :: integerAttribute
     double precision                              :: doubleAttribute
     type            (varying_string             ) :: textAttribute
  end type treeMetaData

  type mergerTreeData
     !!{
     A structure that holds raw merger tree data.
     !!}
     private
     integer                                                                                                :: dummyHostId                        , nodeCount                        , &
          &                                                                                                    particlesCount                     , forestCount
     double precision                                                                                       :: particleMass               =0.0d0
     type            (enumerationPropertyTypeType), allocatable, dimension(:                              ) :: columnProperties                   , particleColumnProperties 
     integer                                      , allocatable, dimension(:                              ) :: treeBeginsAt                       , treeNodeCount
     integer         (kind=kind_int8             ), allocatable, dimension(:                              ) :: descendantIndex                    , hostIndex                        , &
          &                                                                                                    mostBoundParticleIndex             , nodeIndex                        , &
          &                                                                                                    particleIndex                      , forestID                         , &
          &                                                                                                    forestIndex
     integer         (c_size_t                   ), allocatable, dimension(:                              ) :: particleCount                      , snapshot                         , &
          &                                                                                                    particleReferenceCount             , particleReferenceStart           , &
          &                                                                                                    particleSnapshot
     double precision                             , allocatable, dimension(:                              ) :: angularMomentumMagnitude           , halfMassRadius                   , &
          &                                                                                                    nodeMass                           , particleRedshift                 , &
          &                                                                                                    redshift                           , scaleFactor                      , &
          &                                                                                                    scaleRadius                        , spinMagnitude                    , &
          &                                                                                                    treeWeight                         , forestWeightNode                 , &
          &                                                                                                    specificAngularMomentumMagnitude   , velocityMaximum                  , &
          &                                                                                                    velocityDispersion                 , nodeMass200Mean                  , &
          &                                                                                                    nodeMass200Crit
     double precision                             , allocatable, dimension(:,:                            ) :: angularMomentum                    , particlePosition                 , &
          &                                                                                                    particleVelocity                   , position                         , &
          &                                                                                                    spin                               , velocity                         , &
          &                                                                                                    specificAngularMomentum
     double precision                                          , dimension(propertyTypeMin:propertyTypeMax) :: convertProperty            =1.0d0
     logical                                                                                                :: hasAngularMomentumMagnitude        , hasAngularMomentumX              , &
          &                                                                                                    hasAngularMomentumY                , hasAngularMomentumZ              , &
          &                                                                                                    hasSpecificAngularMomentumMagnitude, hasSpecificAngularMomentumX      , &
          &                                                                                                    hasSpecificAngularMomentumY        , hasSpecificAngularMomentumZ      , &
          &                                                                                                    hasDescendantIndex                 , hasDummyHostId                   , &
          &                                                                                                    hasHalfMassRadius                                                     , &
          &                                                                                                    hasHostIndex                       , hasMostBoundParticleIndex        , &
          &                                                                                                    hasNodeIndex                       , hasNodeMass                      , &
          &                                                                                                    hasNodeMass200Mean                 , hasNodeMass200Crit               , &
          &                                                                                                    hasParticleCount                   , hasParticleIndex                 , &
          &                                                                                                    hasParticlePositionX               , hasParticlePositionY             , &
          &                                                                                                    hasParticlePositionZ               , hasParticleRedshift              , &
          &                                                                                                    hasParticleSnapshot                , hasParticleVelocityX             , &
          &                                                                                                    hasParticleVelocityY               , hasParticleVelocityZ             , &
          &                                                                                                    hasParticles               =.false., hasPositionX                     , &
          &                                                                                                    hasPositionY                       , hasPositionZ                     , &
          &                                                                                                    hasRedshift                        , hasScaleFactor                   , &
          &                                                                                                    hasScaleRadius                     , hasSnapshot                      , &
          &                                                                                                    hasSpinMagnitude                   , hasSpinX                         , &
          &                                                                                                    hasSpinY                           , hasSpinZ                         , &
          &                                                                                                    hasForestIndex                     , hasVelocityX                     , &
          &                                                                                                    hasVelocityY                       , hasVelocityZ                     , &
          &                                                                                                    hasVelocityMaximum                 , hasVelocityDispersion            , &
          &                                                                                                    hasForestWeight                                                       , &
          &                                                                                                    hasBoxSize
     logical                                                                                                :: areSelfContained           =.true. , doMakeReferences         =.true. , &
          &                                                                                                    includesHubbleFlow         =.false., includesSubhaloMasses    =.false., &
          &                                                                                                    isPeriodic                 =.false.
     type            (unitsMetaData              )             , dimension(unitsMin       :unitsMax       ) :: units
     logical                                                   , dimension(unitsMin       :unitsMax       ) :: unitsSet                   =.false.
     integer                                                                                                :: metaDataCount              =0
     type            (treeMetaData               ), allocatable, dimension(               :               ) :: metaData
   contains
     !![
     <methods>
       <method description="Reset the data structure." method="reset" />
       <method description="Set the total number of forests in the data structure." method="forestCountSet" />
       <method description="Set the total number of nodes in the data structure." method="nodeCountSet" />
       <method description="Set the total number of particles in the data structure." method="particleCountSet" />
       <method description="Read node data from an ASCII file into the data structure." method="readASCII" />
       <method description="Read particle data from an ASCII file into the data structure" method="readParticlesASCII" />
       <method description="Set a node property in the data structure." method="setProperty" />
       <method description="Set the column in an ASCII data file corresponding to a given node property." method="setPropertyColumn" />
       <method description="Set the column in an ASCII data file corresponding to a given particle property." method="setParticlePropertyColumn" />
       <method description="Set the mass of an N-body particle in the simulation from which the trees were derived." method="setParticleMass" />
       <method description="Specify if trees are self-contained (i.e. contain no cross-links to other trees)." method="setSelfContained" />
       <method description="Specify if velocities include the Hubble flow." method="setIncludesHubbleFlow" />
       <method description="Set if positions are periodic." method="setPositionsArePeriodic" />
       <method description="Set whether halo masses include the masses of the subhalos." method="setIncludesSubhaloMasses" />
       <method description="Set host ID for self-hosting halos if host ID is not node ID." method="setDummyHostId" />
       <method description="Set property type and conversion factor to adjust inconsistent units." method="setConversionFactor" />
       <method description="Set the units used." method="setUnits" />
       <method description="Add a metadatum to the tree data structure." method="addMetadata" />
       <method description="Specify whether or not merger tree dataset references should be made." method="makeReferences" />
       <method description="Export the tree data to an output file." method="export" />
     </methods>
     !!]
     procedure :: reset                                           =>Merger_Tree_Data_Structure_Reset
     procedure :: forestCountSet                                  =>Merger_Tree_Data_Structure_Set_Forest_Count
     procedure :: nodeCountSet                                    =>Merger_Tree_Data_Structure_Set_Node_Count
     procedure :: particleCountSet                                =>Merger_Tree_Data_Structure_Set_Particle_Count
     procedure :: readASCII                                       =>Merger_Tree_Data_Structure_Read_ASCII
     procedure :: readParticlesASCII                              =>Merger_Tree_Data_Structure_Read_Particles_ASCII
     procedure :: Merger_Tree_Data_Structure_Set_Property_Integer8
     procedure :: Merger_Tree_Data_Structure_Set_Property_Double
     generic                                              :: setProperty               => Merger_Tree_Data_Structure_Set_Property_Integer8, &
          &                                                                               Merger_Tree_Data_Structure_Set_Property_Double
     procedure :: setPropertyColumn                              =>Merger_Tree_Data_Structure_Set_Property_Column
     procedure :: setParticlePropertyColumn                      =>Merger_Tree_Data_Structure_Set_Particle_Property_Column
     procedure :: setParticleMass                                =>Merger_Tree_Data_Structure_Set_Particle_Mass
     procedure :: setSelfContained                               =>Merger_Tree_Data_Structure_Set_Self_Contained
     procedure :: setIncludesHubbleFlow                          =>Merger_Tree_Data_Structure_Set_Includes_Hubble_Flow
     procedure :: setPositionsArePeriodic                        =>Merger_Tree_Data_Structure_Set_Is_Periodic
     procedure :: setIncludesSubhaloMasses                       =>Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses
     procedure :: setDummyHostId                                 =>Merger_Tree_Data_Structure_Set_Self_Hosting_Halo_Id
     procedure :: setConversionFactor                            =>Merger_Tree_Data_Structure_Set_Conversion_Factor
     procedure :: setUnits                                       =>Merger_Tree_Data_Structure_Set_Units
     procedure :: Merger_Tree_Data_Structure_Add_Metadata_Double
     procedure :: Merger_Tree_Data_Structure_Add_Metadata_Integer
     procedure :: Merger_Tree_Data_Structure_Add_Metadata_Text
     generic                                              :: addMetadata               => Merger_Tree_Data_Structure_Add_Metadata_Double , &
          &                                                                               Merger_Tree_Data_Structure_Add_Metadata_Integer, &
          &                                                                               Merger_Tree_Data_Structure_Add_Metadata_Text
     procedure :: makeReferences=>Merger_Tree_Data_Structure_Make_References
     procedure :: export        =>Merger_Tree_Data_Structure_Export
  end type mergerTreeData

contains

  subroutine Merger_Tree_Data_Structure_Reset(mergerTrees)
    !!{
    Reset a merger tree data object.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                                :: i

    ! No properties.
    mergerTrees%hasForestIndex                     =.false.
    mergerTrees%hasForestWeight                    =.false.
    mergerTrees%hasBoxSize                         =.false.
    mergerTrees%hasNodeIndex                       =.false.
    mergerTrees%hasDescendantIndex                 =.false.
    mergerTrees%hasHostIndex                       =.false.
    mergerTrees%hasDummyHostID                     =.false.    
    mergerTrees%hasRedshift                        =.false.
    mergerTrees%hasScaleFactor                     =.false.
    mergerTrees%hasNodeMass                        =.false.
    mergerTrees%hasNodeMass200Mean                 =.false.
    mergerTrees%hasNodeMass200Crit                 =.false.
    mergerTrees%hasParticleCount                   =.false.
    mergerTrees%hasPositionX                       =.false.
    mergerTrees%hasPositionY                       =.false.
    mergerTrees%hasPositionZ                       =.false.
    mergerTrees%hasVelocityX                       =.false.
    mergerTrees%hasVelocityY                       =.false.
    mergerTrees%hasVelocityZ                       =.false.
    mergerTrees%hasSpinX                           =.false.
    mergerTrees%hasSpinY                           =.false.
    mergerTrees%hasSpinZ                           =.false.
    mergerTrees%hasSpinMagnitude                   =.false.
    mergerTrees%hasAngularMomentumX                =.false.
    mergerTrees%hasAngularMomentumY                =.false.
    mergerTrees%hasAngularMomentumZ                =.false.
    mergerTrees%hasAngularMomentumMagnitude        =.false.
    mergerTrees%hasSpecificAngularMomentumX        =.false.
    mergerTrees%hasSpecificAngularMomentumY        =.false.
    mergerTrees%hasSpecificAngularMomentumZ        =.false.
    mergerTrees%hasSpecificAngularMomentumMagnitude=.false.
    mergerTrees%hasHalfMassRadius                  =.false.
    mergerTrees%hasScaleRadius                     =.false.
    mergerTrees%hasMostBoundParticleIndex          =.false.
    mergerTrees%hasVelocityMaximum                 =.false.
    mergerTrees%hasVelocityDispersion              =.false.

    ! Deallocate any previous data.
    if (allocated(mergerTrees%treeBeginsAt          )) deallocate(mergerTrees%treeBeginsAt          )
    if (allocated(mergerTrees%treeNodeCount         )) deallocate(mergerTrees%treeNodeCount         )
    if (allocated(mergerTrees%forestID              )) deallocate(mergerTrees%forestID              )
    if (allocated(mergerTrees%treeWeight            )) deallocate(mergerTrees%treeWeight            )
    if (allocated(mergerTrees%forestWeightNode      )) deallocate(mergerTrees%forestWeightNode      )
    if (allocated(mergerTrees%forestIndex           )) deallocate(mergerTrees%forestIndex           )
    if (allocated(mergerTrees%nodeIndex             )) deallocate(mergerTrees%nodeIndex             )
    if (allocated(mergerTrees%mostBoundParticleIndex)) deallocate(mergerTrees%mostBoundParticleIndex)
    if (allocated(mergerTrees%descendantIndex       )) deallocate(mergerTrees%descendantIndex       )
    if (allocated(mergerTrees%hostIndex             )) deallocate(mergerTrees%hostIndex             )
    if (allocated(mergerTrees%redshift              )) deallocate(mergerTrees%redshift              )
    if (allocated(mergerTrees%scaleFactor           )) deallocate(mergerTrees%scaleFactor           )
    if (allocated(mergerTrees%nodeMass              )) deallocate(mergerTrees%nodeMass              )
    if (allocated(mergerTrees%nodeMass200Mean       )) deallocate(mergerTrees%nodeMass200Mean       )
    if (allocated(mergerTrees%nodeMass200Crit       )) deallocate(mergerTrees%nodeMass200Crit       )
    if (allocated(mergerTrees%particleCount         )) deallocate(mergerTrees%particleCount         )
    if (allocated(mergerTrees%position              )) deallocate(mergerTrees%position              )
    if (allocated(mergerTrees%velocity              )) deallocate(mergerTrees%velocity              )
    if (allocated(mergerTrees%spin                  )) deallocate(mergerTrees%spin                  )
    if (allocated(mergerTrees%halfMassRadius        )) deallocate(mergerTrees%halfMassRadius        )
    if (allocated(mergerTrees%scaleRadius           )) deallocate(mergerTrees%scaleRadius           )
    if (allocated(mergerTrees%velocityMaximum       )) deallocate(mergerTrees%velocityMaximum       )
    if (allocated(mergerTrees%velocityDispersion    )) deallocate(mergerTrees%velocityDispersion    )
    do i=unitsMin,unitsMax
       if (allocated(mergerTrees%units(i)%name)) deallocate(mergerTrees%units(i)%name)
    end do
    return
  end subroutine Merger_Tree_Data_Structure_Reset

  subroutine Merger_Tree_Data_Structure_Make_References(mergerTrees,makeReferences)
    !!{
    Specify whether or not to make merger tree dataset references.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: makeReferences

    mergerTrees%doMakeReferences=makeReferences
    return
  end subroutine Merger_Tree_Data_Structure_Make_References

  subroutine Merger_Tree_Data_Structure_Add_Metadata_Double(mergerTrees,metadataType,label,doubleValue)
    !!{
    Add a double metadatum.
    !!}
    implicit none
    class           (mergerTreeData             ), intent(inout) :: mergerTrees
    type            (enumerationMetaDataTypeType), intent(in   ) :: metadataType
    character       (len=*                      ), intent(in   ) :: label
    double precision                             , intent(in   ) :: doubleValue

    call Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,doubleValue=doubleValue)
    ! Check if this is box size.
    if (metadataType == metaDataTypeSimulation .and. trim(label) == "boxSize") mergerTrees%hasBoxSize=.true.
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata_Double

  subroutine Merger_Tree_Data_Structure_Add_Metadata_Integer(mergerTrees,metadataType,label,integerValue)
    !!{
    Add an integer metadatum.
    !!}
    implicit none
    class    (mergerTreeData             ), intent(inout) :: mergerTrees
    type     (enumerationMetaDataTypeType), intent(in   ) :: metadataType
    character(len=*                      ), intent(in   ) :: label
    integer                               , intent(in   ) :: integerValue

    call Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,integerValue=integerValue)
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata_Integer

  subroutine Merger_Tree_Data_Structure_Add_Metadata_Text(mergerTrees,metadataType,label,textValue)
    !!{
    Add a double metadatum.
    !!}
    implicit none
    class    (mergerTreeData             ), intent(inout) :: mergerTrees
    type     (enumerationMetaDataTypeType), intent(in   ) :: metadataType
    character(len=*                      ), intent(in   ) :: label
    character(len=*                      ), intent(in   ) :: textValue

    call Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,textValue=textValue)
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata_Text

  subroutine Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,integerValue,doubleValue,textValue)
    !!{
    Add a metadatum.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class           (mergerTreeData             ), intent(inout)               :: mergerTrees
    type            (enumerationMetaDataTypeType), intent(in   )               :: metadataType
    character       (len=*                      ), intent(in   )               :: label
    integer                                      , intent(in   ), optional     :: integerValue
    double precision                             , intent(in   ), optional     :: doubleValue
    character       (len=*                      ), intent(in   ), optional     :: textValue
    integer                                      , parameter                   :: metadataBlockSize=100
    type            (treeMetaData               ), allocatable  , dimension(:) :: metaDataTemporary

    ! Validate the metadata type.
    if (.not.enumerationMetadataTypeIsValid(metadataType)) call Error_Report('invalid metadata type'//{introspection:location})

    ! Ensure we have enough space in the metadata properties array.
    if (mergerTrees%metaDataCount == 0) then
       allocate(mergerTrees%metaData(metadataBlockSize))
    else if (mergerTrees%metaDataCount == mergerTrees%metaDataCount) then
       call Move_Alloc(mergerTrees%metaData,metaDataTemporary)
       allocate(mergerTrees%metaData(size(metaDataTemporary)+metadataBlockSize))
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
       if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType /= dataTypeNull) call Error_Report('only one data type can be specified'//{introspection:location})
       mergerTrees%metaData(mergerTrees%metaDataCount)%integerAttribute=integerValue
       mergerTrees%metaData(mergerTrees%metaDataCount)%dataType        =dataTypeInteger
    end if
    if (present(doubleValue)) then
       if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType /= dataTypeNull) call Error_Report('only one data type can be specified'//{introspection:location})
       mergerTrees%metaData(mergerTrees%metaDataCount)%doubleAttribute =doubleValue
       mergerTrees%metaData(mergerTrees%metaDataCount)%dataType        =dataTypeDouble
    end if
    if (present(textValue)) then
       if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType /= dataTypeNull) call Error_Report('only one data type can be specified'//{introspection:location})
       mergerTrees%metaData(mergerTrees%metaDataCount)%textAttribute   =textValue
       mergerTrees%metaData(mergerTrees%metaDataCount)%dataType        =dataTypeText
    else
       mergerTrees%metaData(mergerTrees%metaDataCount)%textAttribute   =""
    end if
    if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType == dataTypeNull) call Error_Report('no data was given'//{introspection:location})
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata

  subroutine Merger_Tree_Data_Structure_Set_Forest_Count(mergerTrees,forestCount)
    !!{
    Set the total number of trees in merger trees.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                , intent(in   ) :: forestCount

    ! Set the number of trees.
    mergerTrees%forestCount=forestCount

    return
  end subroutine Merger_Tree_Data_Structure_Set_Forest_Count

  subroutine Merger_Tree_Data_Structure_Set_Node_Count(mergerTrees,nodeCount)
    !!{
    Set the total number of nodes in merger trees.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                , intent(in   ) :: nodeCount

    ! Set the number of nodes.
    mergerTrees%nodeCount=nodeCount

    return
  end subroutine Merger_Tree_Data_Structure_Set_Node_Count

  subroutine Merger_Tree_Data_Structure_Set_Particle_Count(mergerTrees,particleCount)
    !!{
    Set the total number of particles in merger trees.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                , intent(in   ) :: particleCount

    ! Set the number of nodes.
    mergerTrees%particlesCount=particleCount

    return
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Count

  subroutine Merger_Tree_Data_Structure_Set_Particle_Mass(mergerTrees,particleMass)
    !!{
    Set the particle mass used in the trees.
    !!}
    implicit none
    class           (mergerTreeData), intent(inout) :: mergerTrees
    double precision                , intent(in   ) :: particleMass

    ! Set the particle mass.
    mergerTrees%particleMass=particleMass

    return
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Mass

  subroutine Merger_Tree_Data_Structure_Set_Self_Contained(mergerTrees,areSelfContained)
    !!{
    Set the particle mass used in the trees.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: areSelfContained

    ! Set whether trees are self-contained.
    mergerTrees%areSelfContained=areSelfContained

    return
  end subroutine Merger_Tree_Data_Structure_Set_Self_Contained

  subroutine Merger_Tree_Data_Structure_Set_Includes_Hubble_Flow(mergerTrees,includesHubbleFlow)
    !!{
    Set the particle mass used in the trees.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: includesHubbleFlow

    ! Set whether velocities include the Hubble flow.
    mergerTrees%includesHubbleFlow=includesHubbleFlow

    return
  end subroutine Merger_Tree_Data_Structure_Set_Includes_Hubble_Flow

  subroutine Merger_Tree_Data_Structure_Set_Is_Periodic(mergerTrees,isPeriodic)
    !!{
    Set whether or not positions are periodic.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: isPeriodic

    ! Set whether positions are periodic.
    mergerTrees%isPeriodic=isPeriodic

    return
  end subroutine Merger_Tree_Data_Structure_Set_Is_Periodic

  subroutine Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses(mergerTrees,includesSubhaloMasses)
    !!{
    Set the particle mass used in the trees.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: includesSubhaloMasses

    ! Set whether halo masses include subhalo contributions.
    mergerTrees%includesSubhaloMasses=includesSubhaloMasses

    return
  end subroutine Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses

  subroutine Merger_Tree_Data_Structure_Set_Self_Hosting_Halo_Id(mergerTrees,dummyHostId)
    !!{
    Set the host ID in case of self-hosting halos. Default is host ID = node ID.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                , intent(in   ) :: dummyHostId

    ! Set the value of the dummy-variable for self hosting halos.
    mergerTrees%dummyHostId   =dummyHostId
    mergerTrees%hasDummyHostId=.true.
    return
  end subroutine Merger_Tree_Data_Structure_Set_Self_Hosting_Halo_Id

  subroutine Merger_Tree_Data_Structure_Set_Conversion_Factor(mergerTrees,propertyType,conversionFactor)
    !!{
    Set Conversion factor for property type with inconsistent unit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (mergerTreeData             ), intent(inout) :: mergerTrees
    type            (enumerationPropertyTypeType), intent(in   ) :: propertyType
    double precision                             , intent(in   ) :: conversionFactor

    ! Ensure the property type is valid.
    if (.not.enumerationPropertyTypeIsValid(propertyType)) call Error_Report('invalid property type'//{introspection:location})

    ! Store conversion factor into array.
    mergerTrees%convertProperty(propertyType%ID)=conversionFactor
    return
  end subroutine Merger_Tree_Data_Structure_Set_Conversion_Factor

  subroutine Merger_Tree_Data_Structure_Set_Units(mergerTrees,unitType,unitsInSI,hubbleExponent,scaleFactorExponent,name)
    !!{
    Set the units system.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class           (mergerTreeData      ), intent(inout)           :: mergerTrees
    type            (enumerationUnitsType), intent(in   )           :: unitType
    double precision                      , intent(in   )           :: unitsInSI
    integer                               , intent(in   ), optional :: hubbleExponent, scaleFactorExponent
    character       (len=*               ), intent(in   ), optional :: name

    ! Ensure the unit type is valid.
    if (.not.enumerationUnitsIsValid(unitType)) call Error_Report('invalid unit type'//{introspection:location})

    ! Flag the units as set.
    mergerTrees%unitsSet(unitType%ID)=.true.

    ! Store the units in the SI system.
    mergerTrees%units(unitType%ID)%unitsInSI=unitsInSI

    ! Store Hubble parameter exponent if given.
    if (present(hubbleExponent)) then
       mergerTrees%units(unitType%ID)%hubbleExponent=hubbleExponent
    else
       ! No Hubble parameter exponent provided - assume no dependence.
       mergerTrees%units(unitType%ID)%hubbleExponent=0
    end if

    ! Store scale factor exponent if given.
    if (present(scaleFactorExponent)) then
       mergerTrees%units(unitType%ID)%scaleFactorExponent=scaleFactorExponent
    else
       ! No scale factor parameter exponent provided - assume no dependence.
       mergerTrees%units(unitType%ID)%scaleFactorExponent=0
    end if

    ! Store the name if given.
    allocate(mergerTrees%units(unitType%ID)%name)
    if (present(name)) then
       mergerTrees%units(unitType%ID)%name=name
    else
       mergerTrees%units(unitType%ID)%name=""
    end if
    return
  end subroutine Merger_Tree_Data_Structure_Set_Units

  subroutine Merger_Tree_Data_Structure_Set_Particle_Property_Column(mergerTrees,propertyType,columnNumber)
    !!{
    Set column mapping from the input file.
    !!}
    implicit none
    class  (mergerTreeData             ), intent(inout)               :: mergerTrees
    integer                             , intent(in   )               :: columnNumber
    type   (enumerationPropertyTypeType), intent(in   )               :: propertyType
    type   (enumerationPropertyTypeType), allocatable  , dimension(:) :: columnPropertiesTemp

    ! Ensure the storage array is large enough.
    if (allocated(mergerTrees%particleColumnProperties)) then
       if (columnNumber > size(mergerTrees%particleColumnProperties)) then
          call Move_Alloc(mergerTrees%particleColumnProperties,columnPropertiesTemp)
          allocate(mergerTrees%particleColumnProperties(columnNumber))
          mergerTrees%particleColumnProperties(1                           :size(columnPropertiesTemp))=columnPropertiesTemp
          mergerTrees%particleColumnProperties(1+size(columnPropertiesTemp):columnNumber              )=propertyTypeNull
          deallocate(columnPropertiesTemp)
       end if
    else
       allocate(mergerTrees%particleColumnProperties(columnNumber))
       mergerTrees%particleColumnProperties=propertyTypeNull
    end if
    ! Store the property type.
    mergerTrees%particleColumnProperties(columnNumber)=propertyType
    return
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Property_Column

  subroutine Merger_Tree_Data_Structure_Set_Property_Column(mergerTrees,propertyType,columnNumber)
    !!{
    Set column mapping from the input file.
    !!}
    implicit none
    class  (mergerTreeData             ), intent(inout)               :: mergerTrees
    type   (enumerationPropertyTypeType), intent(in   )               :: propertyType
    integer                             , intent(in   )               :: columnNumber
    type   (enumerationPropertyTypeType), allocatable  , dimension(:) :: columnPropertiesTemp

    ! Ensure the storage array is large enough.
    if (allocated(mergerTrees%columnProperties)) then
       if (columnNumber > size(mergerTrees%columnProperties)) then
          call Move_Alloc(mergerTrees%columnProperties,columnPropertiesTemp)
          allocate(mergerTrees%columnProperties(columnNumber))
          mergerTrees%columnProperties(1                           :size(columnPropertiesTemp))=columnPropertiesTemp
          mergerTrees%columnProperties(1+size(columnPropertiesTemp):columnNumber              )=propertyTypeNull
          deallocate(columnPropertiesTemp)
       end if
    else
       allocate(mergerTrees%columnProperties(columnNumber))
       mergerTrees%columnProperties=propertyTypeNull
    end if
    ! Store the property type.
    mergerTrees%columnProperties(columnNumber)=propertyType
    return
  end subroutine Merger_Tree_Data_Structure_Set_Property_Column

  subroutine Merger_Tree_Data_Structure_Set_Property_Integer8(mergerTrees,propertyType,property)
    !!{
    Set a property in the merger trees.
    !!}
    use :: Error            , only : Error_Report
    implicit none
    class  (mergerTreeData             )              , intent(inout) :: mergerTrees
    type   (enumerationPropertyTypeType)              , intent(in   ) :: propertyType
    integer(kind=kind_int8             ), dimension(:), intent(in   ) :: property

    ! Check the supplied arrays is of the correct size.
    if (size(property) /= mergerTrees%nodeCount) call Error_Report('property array size is incorrect'//{introspection:location})

    ! Assign to the relevant property.
    select case (propertyType%ID)
    case (propertyTypeTreeIndex      %ID)
       mergerTrees%hasForestIndex    =.true.
       if (allocated(mergerTrees%forestIndex    )) deallocate(mergerTrees%forestIndex    )
       allocate(mergerTrees%forestIndex    (size(property)))
       mergerTrees%forestIndex    =property
       call Merger_Tree_Data_Structure_Set_Tree_Indices(mergerTrees)
    case (propertyTypeNodeIndex      %ID)
       mergerTrees%hasNodeIndex      =.true.
       if (allocated(mergerTrees%nodeIndex      )) deallocate(mergerTrees%nodeIndex      )
       allocate(mergerTrees%nodeIndex      (size(property)))
       mergerTrees%nodeIndex      =property
    case (propertyTypeHostIndex      %ID)
       mergerTrees%hasHostIndex      =.true.
       if (allocated(mergerTrees%hostIndex      )) deallocate(mergerTrees%hostIndex      )
       allocate(mergerTrees%hostIndex      (size(property)))
       mergerTrees%hostIndex      =property
    case (propertyTypeDescendantIndex%ID)
       mergerTrees%hasDescendantIndex=.true.
       if (allocated(mergerTrees%descendantIndex)) deallocate(mergerTrees%descendantIndex)
       allocate(mergerTrees%descendantIndex(size(property)))
       mergerTrees%descendantIndex=property
    case (propertyTypeSnapshot       %ID)
       mergerTrees%hasSnapshot       =.true.
       if (allocated(mergerTrees%snapshot       )) deallocate(mergerTrees%snapshot       )
       allocate(mergerTrees%snapshot       (size(property)))
       mergerTrees%snapshot=property
    case default
       call Error_Report('unrecognized integer property'//{introspection:location})
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Set_Property_Integer8

  subroutine Merger_Tree_Data_Structure_Set_Property_Double(mergerTrees,propertyType,property)
    !!{
    Set a property in the merger trees.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (mergerTreeData             )              , intent(inout) :: mergerTrees
    type            (enumerationPropertyTypeType)              , intent(in   ) :: propertyType
    double precision                             , dimension(:), intent(in   ) :: property

    ! Check the supplied arrays is of the correct size.
    if (size(property) /= mergerTrees%nodeCount) call Error_Report('property array size is incorrect'//{introspection:location})

    ! Assign to the relevant property.
    select case (propertyType%ID)
    case (propertyTypeTreeWeight     %ID)
       mergerTrees%hasForestWeight           =.true.
       if (allocated(mergerTrees%forestWeightNode        )) deallocate(mergerTrees%forestWeightNode)
       allocate(mergerTrees%forestWeightNode        (size(property)))
       mergerTrees%forestWeightNode =property
    case (propertyTypeRedshift       %ID)
       mergerTrees%hasRedshift                =.true.
       if (allocated(mergerTrees%redshift                )) deallocate(mergerTrees%redshift        )
       allocate(mergerTrees%redshift                (size(property)))
       mergerTrees%redshift       =property
    case (propertyTypeNodeMass       %ID)
       mergerTrees%hasNodeMass                =.true.
       if (allocated(mergerTrees%nodeMass                )) deallocate(mergerTrees%nodeMass        )
       allocate(mergerTrees%nodeMass                (size(property)))
       mergerTrees%nodeMass       =property
    case (propertyTypeNodeMass200Mean%ID)
       mergerTrees%hasNodeMass200Mean         =.true.
       if (allocated(mergerTrees%nodeMass200Mean         )) deallocate(mergerTrees%nodeMass200Mean )
       allocate(mergerTrees%nodeMass200Mean         (size(property)))
       mergerTrees%nodeMass200Mean=property
    case (propertyTypeNodeMass200Crit%ID)
       mergerTrees%hasNodeMass200Crit         =.true.
       if (allocated(mergerTrees%nodeMass200Crit         )) deallocate(mergerTrees%nodeMass200Crit )
       allocate(mergerTrees%nodeMass200Crit         (size(property)))
       mergerTrees%nodeMass200Crit=property
    case (propertyTypeScaleRadius    %ID)
       mergerTrees%hasScaleRadius             =.true.
       if (allocated(mergerTrees%scaleRadius             )) deallocate(mergerTrees%scaleRadius )
       allocate(mergerTrees%scaleRadius             (size(property)))
       mergerTrees%scaleRadius=property
    case (propertyTypeAngularMomentum%ID)
       mergerTrees%hasAngularMomentumMagnitude=.true.
       if (allocated(mergerTrees%angularMomentumMagnitude)) deallocate(mergerTrees%angularMomentumMagnitude )
       allocate(mergerTrees%angularMomentumMagnitude(size(property)))
       mergerTrees%angularMomentumMagnitude=property
    case (propertyTypePositionX      %ID)
       if (                                             &
            & allocated(mergerTrees%position)           &
            & .and..not.mergerTrees%hasPositionX        &
            & .and..not.mergerTrees%hasPositionY        &
            & .and..not.mergerTrees%hasPositionZ        &
            &) deallocate(mergerTrees%position)
       mergerTrees%hasPositionX       =.true.
       if (.not.allocated(mergerTrees%position       )) allocate(mergerTrees%position       (3,size(property)))
       mergerTrees%position       (1,:)=property
    case (propertyTypePositionY      %ID)
       if (                                             &
            & allocated(mergerTrees%position)           &
            & .and..not.mergerTrees%hasPositionX        &
            & .and..not.mergerTrees%hasPositionY        &
            & .and..not.mergerTrees%hasPositionZ        &
            &) deallocate(mergerTrees%position)
       mergerTrees%hasPositionY       =.true.
       if (.not.allocated(mergerTrees%position       )) allocate(mergerTrees%position       (3,size(property)))
       mergerTrees%position       (2,:)=property
    case (propertyTypePositionZ      %ID)
       if (                                             &
            & allocated(mergerTrees%position)           &
            & .and..not.mergerTrees%hasPositionX        &
            & .and..not.mergerTrees%hasPositionY        &
            & .and..not.mergerTrees%hasPositionZ        &
            &) deallocate(mergerTrees%position)
       mergerTrees%hasPositionZ       =.true.
       if (.not.allocated(mergerTrees%position       )) allocate(mergerTrees%position       (3,size(property)))
       mergerTrees%position       (3,:)=property
    case (propertyTypeVelocityX      %ID)
       if (                                             &
            & allocated(mergerTrees%velocity)           &
            & .and..not.mergerTrees%hasVelocityX        &
            & .and..not.mergerTrees%hasVelocityY        &
            & .and..not.mergerTrees%hasVelocityZ        &
            &) deallocate(mergerTrees%velocity)
       mergerTrees%hasVelocityX       =.true.
       if (.not.allocated(mergerTrees%velocity       )) allocate(mergerTrees%velocity       (3,size(property)))
       mergerTrees%velocity       (1,:)=property
    case (propertyTypeVelocityY      %ID)
       if (                                             &
            & allocated(mergerTrees%velocity)           &
            & .and..not.mergerTrees%hasVelocityX        &
            & .and..not.mergerTrees%hasVelocityY        &
            & .and..not.mergerTrees%hasVelocityZ        &
            &) deallocate(mergerTrees%velocity)
       mergerTrees%hasVelocityY       =.true.
       if (.not.allocated(mergerTrees%velocity       )) allocate(mergerTrees%velocity       (3,size(property)))
       mergerTrees%velocity       (2,:)=property
    case (propertyTypeVelocityZ            %ID)
       if (                                             &
            & allocated(mergerTrees%velocity)           &
            & .and..not.mergerTrees%hasVelocityX        &
            & .and..not.mergerTrees%hasVelocityY        &
            & .and..not.mergerTrees%hasVelocityZ        &
            &) deallocate(mergerTrees%velocity)
       mergerTrees%hasVelocityZ       =.true.
       if (.not.allocated(mergerTrees%velocity       )) allocate(mergerTrees%velocity       (3,size(property)))
       mergerTrees%velocity       (3,:)=property
    case (propertyTypeAngularMomentumX      %ID)
       if (                                             &
            & allocated(mergerTrees%angularMomentum)    &
            & .and..not.mergerTrees%hasAngularMomentumX &
            & .and..not.mergerTrees%hasAngularMomentumY &
            & .and..not.mergerTrees%hasAngularMomentumZ &
            &) deallocate(mergerTrees%angularMomentum)
       mergerTrees%hasAngularMomentumX=.true.
       if (.not.allocated(mergerTrees%angularMomentum)) allocate(mergerTrees%angularMomentum(3,size(property)))
       mergerTrees%angularMomentum(1,:)=property
    case (propertyTypeAngularMomentumY      %ID)
       if (                                             &
            & allocated(mergerTrees%angularMomentum)    &
            & .and..not.mergerTrees%hasAngularMomentumX &
            & .and..not.mergerTrees%hasAngularMomentumY &
            & .and..not.mergerTrees%hasAngularMomentumZ &
            &) deallocate(mergerTrees%angularMomentum)
       mergerTrees%hasAngularMomentumY=.true.
       if (.not.allocated(mergerTrees%angularMomentum)) allocate(mergerTrees%angularMomentum(3,size(property)))
       mergerTrees%angularMomentum(2,:)=property
    case (propertyTypeAngularMomentumZ      %ID)
       if (                                             &
            & allocated(mergerTrees%angularMomentum)    &
            & .and..not.mergerTrees%hasAngularMomentumX &
            & .and..not.mergerTrees%hasAngularMomentumY &
            & .and..not.mergerTrees%hasAngularMomentumZ &
            &) deallocate(mergerTrees%angularMomentum)
       mergerTrees%hasAngularMomentumZ=.true.
       if (.not.allocated(mergerTrees%angularMomentum)) allocate(mergerTrees%angularMomentum(3,size(property)))
       mergerTrees%angularMomentum(3,:)=property
    case default
       call Error_Report('unrecognized double property'//{introspection:location})
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Set_Property_Double

  subroutine Merger_Tree_Data_Structure_Read_ASCII(mergerTrees,inputFile,columnHeaders,commentCharacter,separator,maximumRedshift)
    !!{
    Read in merger tree data from an ASCII file.
    !!}
    use :: Display           , only : displayMessage
    use :: File_Utilities    , only : Count_Lines_In_File
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=)      , operator(//)
    use :: String_Handling   , only : String_Count_Words , String_Split_Words, operator(//)
    implicit none
    class    (mergerTreeData), intent(inout)               :: mergerTrees
    character(len=*         ), intent(in   )               :: inputFile
    character(len=*         ), intent(in   ), optional     :: separator
    character(len=1         ), intent(in   ), optional     :: commentCharacter
    double precision         , intent(in   ), optional     :: maximumRedshift
    logical                  , intent(in   ), optional     :: columnHeaders
    character(len=32        ), allocatable  , dimension(:) :: inputColumns
    double precision         , parameter                   :: maximumRedshiftDefault=800.0d0
    integer                                                :: columnsCount                  , fileUnit           , &
         &                                                    iColumn                       , iNode              , &
         &                                                    i                             , nodeCount
    double precision                                       :: maximumRedshiftActual
    logical                                                :: gotFirstDataLine              , gotColumnHeaderLine
    character(len=1024      )                              :: inputLine
    type     (varying_string)                              :: message

    ! Determine number of nodes.
    nodeCount=Count_Lines_In_File(inputFile,commentCharacter)
    if (present(columnHeaders).and.columnHeaders) nodeCount=nodeCount-1
    call mergerTrees%nodeCountSet(nodeCount)

    ! Specify what properties these trees have.
    mergerTrees%hasForestIndex                     =any(mergerTrees%columnProperties == propertyTypeTreeIndex               )
    mergerTrees%hasForestWeight                    =any(mergerTrees%columnProperties == propertyTypeTreeWeight              )
    mergerTrees%hasNodeIndex                       =any(mergerTrees%columnProperties == propertyTypeNodeIndex               )
    mergerTrees%hasDescendantIndex                 =any(mergerTrees%columnProperties == propertyTypeDescendantIndex         )
    mergerTrees%hasHostIndex                       =any(mergerTrees%columnProperties == propertyTypeHostIndex               )
    mergerTrees%hasRedshift                        =any(mergerTrees%columnProperties == propertyTypeRedshift                )
    mergerTrees%hasScaleFactor                     =any(mergerTrees%columnProperties == propertyTypeScaleFactor             )
    mergerTrees%hasNodeMass                        =any(mergerTrees%columnProperties == propertyTypeNodeMass                )
    mergerTrees%hasNodeMass200Mean                 =any(mergerTrees%columnProperties == propertyTypeNodeMass200Mean         )
    mergerTrees%hasNodeMass200Crit                 =any(mergerTrees%columnProperties == propertyTypeNodeMass200Crit         )
    mergerTrees%hasParticleCount                   =any(mergerTrees%columnProperties == propertyTypeParticleCount           )
    mergerTrees%hasPositionX                       =any(mergerTrees%columnProperties == propertyTypePositionX               )
    mergerTrees%hasPositionY                       =any(mergerTrees%columnProperties == propertyTypePositionY               )
    mergerTrees%hasPositionZ                       =any(mergerTrees%columnProperties == propertyTypePositionZ               )
    mergerTrees%hasVelocityX                       =any(mergerTrees%columnProperties == propertyTypeVelocityX               )
    mergerTrees%hasVelocityY                       =any(mergerTrees%columnProperties == propertyTypeVelocityY               )
    mergerTrees%hasVelocityZ                       =any(mergerTrees%columnProperties == propertyTypeVelocityZ               )
    mergerTrees%hasSpinX                           =any(mergerTrees%columnProperties == propertyTypeSpinX                   )
    mergerTrees%hasSpinY                           =any(mergerTrees%columnProperties == propertyTypeSpinY                   )
    mergerTrees%hasSpinZ                           =any(mergerTrees%columnProperties == propertyTypeSpinZ                   )
    mergerTrees%hasSpinMagnitude                   =any(mergerTrees%columnProperties == propertyTypeSpin                    )
    mergerTrees%hasAngularMomentumX                =any(mergerTrees%columnProperties == propertyTypeAngularMomentumX        )
    mergerTrees%hasAngularMomentumY                =any(mergerTrees%columnProperties == propertyTypeAngularMomentumY        )
    mergerTrees%hasAngularMomentumZ                =any(mergerTrees%columnProperties == propertyTypeAngularMomentumZ        )
    mergerTrees%hasAngularMomentumMagnitude        =any(mergerTrees%columnProperties == propertyTypeAngularMomentum         )
    mergerTrees%hasSpecificAngularMomentumX        =any(mergerTrees%columnProperties == propertyTypeSpecificAngularMomentumX)
    mergerTrees%hasSpecificAngularMomentumY        =any(mergerTrees%columnProperties == propertyTypeSpecificAngularMomentumY)
    mergerTrees%hasSpecificAngularMomentumZ        =any(mergerTrees%columnProperties == propertyTypeSpecificAngularMomentumZ)
    mergerTrees%hasSpecificAngularMomentumMagnitude=any(mergerTrees%columnProperties == propertyTypeSpecificAngularMomentum )
    mergerTrees%hasHalfMassRadius                  =any(mergerTrees%columnProperties == propertyTypeHalfMassRadius          )
    mergerTrees%hasScaleRadius                     =any(mergerTrees%columnProperties == propertyTypeScaleRadius             )
    mergerTrees%hasMostBoundParticleIndex          =any(mergerTrees%columnProperties == propertyTypeMostBoundParticleIndex  )
    mergerTrees%hasSnapshot                        =any(mergerTrees%columnProperties == propertyTypeSnapshot                )
    mergerTrees%hasVelocityMaximum                 =any(mergerTrees%columnProperties == propertyTypeVelocityMaximum         )
    mergerTrees%hasVelocityDispersion              =any(mergerTrees%columnProperties == propertyTypeVelocityDispersion      )

    ! Validate 3-D datasets.
    if     (.not.((     mergerTrees%hasPositionX       .and.     mergerTrees%hasPositionY       .and.     mergerTrees%hasPositionZ       ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasPositionX       .and..not.mergerTrees%hasPositionY       .and..not.mergerTrees%hasPositionZ       ) &
         &        )                                                                                                   &
         & ) call Error_Report("all three axes or none must be supplied for position"        //{introspection:location})
    if     (.not.((     mergerTrees%hasVelocityX       .and.     mergerTrees%hasVelocityY       .and.     mergerTrees%hasVelocityZ       ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasVelocityX       .and..not.mergerTrees%hasVelocityY       .and..not.mergerTrees%hasVelocityZ       ) &
         &        )                                                                                                   &
         & ) call Error_Report("all three axes or none must be supplied for velocity"        //{introspection:location})
    if     (.not.((     mergerTrees%hasSpinX           .and.     mergerTrees%hasSpinY           .and.     mergerTrees%hasSpinZ           ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasSpinX           .and..not.mergerTrees%hasSpinY           .and..not.mergerTrees%hasSpinZ           ) &
         &        )                                                                                                   &
         & ) call Error_Report("all three axes or none must be supplied for spin"            //{introspection:location})
    if     (.not.((     mergerTrees%hasAngularMomentumX        .and.     mergerTrees%hasAngularMomentumY        .and.     mergerTrees%hasAngularMomentumZ        ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasAngularMomentumX        .and..not.mergerTrees%hasAngularMomentumY        .and..not.mergerTrees%hasAngularMomentumZ        ) &
         &        )                                                                                                   &
         & ) call Error_Report("all three axes or none must be supplied for angular momentum"         //{introspection:location})
    if     (.not.((     mergerTrees%hasSpecificAngularMomentumX.and.     mergerTrees%hasSpecificAngularMomentumY.and.     mergerTrees%hasSpecificAngularMomentumZ) &
         &        .or.                                                                                                                                             &
         &        (.not.mergerTrees%hasSpecificAngularMomentumX.and..not.mergerTrees%hasSpecificAngularMomentumY.and..not.mergerTrees%hasSpecificAngularMomentumZ) &
         &        )                                                                                                                                                &
         & ) call Error_Report("all three axes or none must be supplied for specific angular momentum"//{introspection:location})
    if (mergerTrees%hasSpinX                   .and.mergerTrees%hasSpinMagnitude                   ) &
         & call Error_Report("can not specify both 3-D and scalar spin"                     //{introspection:location})
    if (mergerTrees%hasAngularMomentumX        .and.mergerTrees%hasAngularMomentumMagnitude        ) &
         & call Error_Report("can not specify both 3-D and scalar angular momentum"         //{introspection:location})
    if (mergerTrees%hasSpecificAngularMomentumX.and.mergerTrees%hasSpecificAngularMomentumMagnitude) &
         & call Error_Report("can not specify both 3-D and scalar specific angular momentum"//{introspection:location})

    ! Validate that either redshift or scale factor is given.
    if (.not.(mergerTrees%hasRedshift .or. mergerTrees%hasScaleFactor)) &
         & call Error_Report("either redshift or scale factor has to be given"//{introspection:location})

    ! Deallocate internal arrays.
    if (allocated(mergerTrees%forestIndex           )) deallocate(mergerTrees%forestIndex           )
    if (allocated(mergerTrees%forestWeightNode      )) deallocate(mergerTrees%forestWeightNode      )
    if (allocated(mergerTrees%nodeIndex             )) deallocate(mergerTrees%nodeIndex             )
    if (allocated(mergerTrees%mostBoundParticleIndex)) deallocate(mergerTrees%mostBoundParticleIndex)
    if (allocated(mergerTrees%snapshot              )) deallocate(mergerTrees%snapshot              )
    if (allocated(mergerTrees%descendantIndex       )) deallocate(mergerTrees%descendantIndex       )
    if (allocated(mergerTrees%hostIndex             )) deallocate(mergerTrees%hostIndex             )
    if (allocated(mergerTrees%redshift              )) deallocate(mergerTrees%redshift              )
    if (allocated(mergerTrees%scaleFactor           )) deallocate(mergerTrees%scaleFactor           )
    if (allocated(mergerTrees%nodeMass              )) deallocate(mergerTrees%nodeMass              )
    if (allocated(mergerTrees%nodeMass200Mean       )) deallocate(mergerTrees%nodeMass200Mean       )
    if (allocated(mergerTrees%nodeMass200Crit       )) deallocate(mergerTrees%nodeMass200Crit       )
    if (allocated(mergerTrees%particleCount         )) deallocate(mergerTrees%particleCount         )
    if (allocated(mergerTrees%position              )) deallocate(mergerTrees%position              )
    if (allocated(mergerTrees%velocity              )) deallocate(mergerTrees%velocity              )
    if (allocated(mergerTrees%spin                  )) deallocate(mergerTrees%spin                  )
    if (allocated(mergerTrees%halfMassRadius        )) deallocate(mergerTrees%halfMassRadius        )
    if (allocated(mergerTrees%scaleRadius           )) deallocate(mergerTrees%scaleRadius           )
    if (allocated(mergerTrees%velocityMaximum       )) deallocate(mergerTrees%velocityMaximum       )
    if (allocated(mergerTrees%velocityDispersion    )) deallocate(mergerTrees%velocityDispersion    )

    ! Allocate internal arrays to correct size as needed.
    if (mergerTrees%hasForestIndex                     ) allocate(mergerTrees%forestIndex                     (mergerTrees%nodeCount))
    if (mergerTrees%hasForestWeight                    ) allocate(mergerTrees%forestWeightNode                (mergerTrees%nodeCount))
    if (mergerTrees%hasNodeIndex                       ) allocate(mergerTrees%nodeIndex                       (mergerTrees%nodeCount))
    if (mergerTrees%hasMostBoundParticleIndex          ) allocate(mergerTrees%mostBoundParticleIndex          (mergerTrees%nodeCount))
    if (mergerTrees%hasSnapshot                        ) allocate(mergerTrees%snapshot                        (mergerTrees%nodeCount))
    if (mergerTrees%hasDescendantIndex                 ) allocate(mergerTrees%descendantIndex                 (mergerTrees%nodeCount))
    if (mergerTrees%hasHostIndex                       ) allocate(mergerTrees%hostIndex                       (mergerTrees%nodeCount))
    if (mergerTrees%hasRedshift                        ) allocate(mergerTrees%redshift                        (mergerTrees%nodeCount))
    if (mergerTrees%hasScaleFactor                     ) allocate(mergerTrees%scaleFactor                     (mergerTrees%nodeCount))
    if (mergerTrees%hasNodeMass                        ) allocate(mergerTrees%nodeMass                        (mergerTrees%nodeCount))
    if (mergerTrees%hasNodeMass200Mean                 ) allocate(mergerTrees%nodeMass200Mean                 (mergerTrees%nodeCount))
    if (mergerTrees%hasNodeMass200Crit                 ) allocate(mergerTrees%nodeMass200Crit                 (mergerTrees%nodeCount))
    if (mergerTrees%hasParticleCount                   ) allocate(mergerTrees%particleCount                   (mergerTrees%nodeCount))
    if (mergerTrees%hasPositionX                       ) allocate(mergerTrees%position                        (3,mergerTrees%nodeCount))
    if (mergerTrees%hasVelocityX                       ) allocate(mergerTrees%velocity                        (3,mergerTrees%nodeCount))
    if (mergerTrees%hasSpinX                           ) allocate(mergerTrees%spin                            (3,mergerTrees%nodeCount))
    if (mergerTrees%hasAngularMomentumX                ) allocate(mergerTrees%angularMomentum                 (3,mergerTrees%nodeCount))
    if (mergerTrees%hasSpecificAngularMomentumX        ) allocate(mergerTrees%specificAngularMomentum         (3,mergerTrees%nodeCount))
    if (mergerTrees%hasSpinMagnitude                   ) allocate(mergerTrees%spinMagnitude                   (mergerTrees%nodeCount))
    if (mergerTrees%hasAngularMomentumMagnitude        ) allocate(mergerTrees%angularMomentumMagnitude        (mergerTrees%nodeCount))
    if (mergerTrees%hasSpecificAngularMomentumMagnitude) allocate(mergerTrees%specificAngularMomentumMagnitude(mergerTrees%nodeCount))
    if (mergerTrees%hasHalfMassRadius                  ) allocate(mergerTrees%halfMassRadius                  (mergerTrees%nodeCount))
    if (mergerTrees%hasScaleRadius                     ) allocate(mergerTrees%scaleRadius                     (mergerTrees%nodeCount))
    if (mergerTrees%hasVelocityMaximum                 ) allocate(mergerTrees%velocityMaximum                 (mergerTrees%nodeCount))
    if (mergerTrees%hasVelocityDispersion              ) allocate(mergerTrees%velocityDispersion              (mergerTrees%nodeCount))

    ! Open the file and read lines.
    open(newunit=fileUnit,file=inputFile,status='old',form='formatted')
    iNode              =0
    columnsCount       =0
    gotFirstDataLine   =.false.
    gotColumnHeaderLine=.not.present(columnHeaders).or..not.columnHeaders
    do while (iNode < nodeCount)
       ! Get the line.
       read (fileUnit,'(a)') inputLine
       ! Check if this is a data line.
       if (.not.present(commentCharacter) .or. inputLine(1:1) /= commentCharacter) then
          ! Skip header line.
          if (.not.gotColumnHeaderLine) then
             gotColumnHeaderLine=.true.
             cycle
          end if
          ! If this is the first data line, determine how many columns are present and allocate array to store them.
          if (.not.gotFirstDataLine) then
             columnsCount=String_Count_Words(inputLine,separator)
             allocate(inputColumns(columnsCount))
             gotFirstDataLine=.true.
          end if
          ! Count nodes.
          iNode=iNode+1
          call String_Split_Words(inputColumns,inputLine,separator)
          do iColumn=1,min(columnsCount,size(mergerTrees%columnProperties))
             select case (mergerTrees%columnProperties(iColumn)%ID)
             case (propertyTypeNull                  %ID)
                ! Ignore this column.
             case (propertyTypeTreeIndex             %ID)
                ! Column is a tree index.
                read (inputColumns(iColumn),*) mergerTrees%forestIndex(iNode)
                if (iNode > 1) then
                   if (mergerTrees%forestIndex(iNode) < mergerTrees%forestIndex(iNode-1)) &
                        & call Error_Report('tree indices must be in ascending order'//{introspection:location})
                   if (mergerTrees%forestIndex(iNode) /= mergerTrees%forestIndex(iNode-1)) mergerTrees%forestCount=mergerTrees%forestCount+1
                else
                   mergerTrees%forestCount=1
                end if
             case (propertyTypeTreeWeight              %ID)
                ! Column is a tree weight.
                read (inputColumns(iColumn),*) mergerTrees%forestWeightNode                (  iNode)
             case (propertyTypeNodeIndex               %ID)
                ! Column is a node index.
                read (inputColumns(iColumn),*) mergerTrees%nodeIndex                       (  iNode)
             case (propertyTypeDescendantIndex         %ID)
                ! Column is a descendant node index.
                read (inputColumns(iColumn),*) mergerTrees%descendantIndex                 (  iNode)
             case (propertyTypeHostIndex               %ID)
                ! Column is a host index.
                read (inputColumns(iColumn),*) mergerTrees%hostIndex                       (  iNode)
             case (propertyTypeRedshift                %ID)
                ! Column is redshift.
                read (inputColumns(iColumn),*) mergerTrees%redshift                        (  iNode)
              case (propertyTypeScaleFactor            %ID)
                ! Column is scale factor.
                read (inputColumns(iColumn),*) mergerTrees%scaleFactor                     (  iNode)
             case (propertyTypeNodeMass                %ID)
                ! Column is mass.
                read (inputColumns(iColumn),*) mergerTrees%nodeMass                        (  iNode)
             case (propertyTypeNodeMass200Mean         %ID)
                ! Column is mass.
                read (inputColumns(iColumn),*) mergerTrees%nodeMass200Mean                 (  iNode)
             case (propertyTypeNodeMass200Crit         %ID)
                ! Column is mass.
                read (inputColumns(iColumn),*) mergerTrees%nodeMass200Crit                 (  iNode)
             case (propertyTypeParticleCount           %ID)
                ! Column is particle count.
                read (inputColumns(iColumn),*) mergerTrees%particleCount                   (  iNode)
             case (propertyTypePositionX               %ID)
                ! Column is x position.
                read (inputColumns(iColumn),*) mergerTrees%position                        (1,iNode)
             case (propertyTypePositionY               %ID)
                ! Column is y position.
                read (inputColumns(iColumn),*) mergerTrees%position                        (2,iNode)
             case (propertyTypePositionZ               %ID)
                ! Column is z position.
                read (inputColumns(iColumn),*) mergerTrees%position                        (3,iNode)
             case (propertyTypeVelocityX               %ID)
                ! Column is x velocity.
                read (inputColumns(iColumn),*) mergerTrees%velocity                        (1,iNode)
             case (propertyTypeVelocityY               %ID)
                ! Column is y velocity.
                read (inputColumns(iColumn),*) mergerTrees%velocity                        (2,iNode)
             case (propertyTypeVelocityZ               %ID)
                ! Column is z velocity.
                read (inputColumns(iColumn),*) mergerTrees%velocity                        (3,iNode)
             case (propertyTypeSpinX                   %ID)
                ! Column is x spin.
                read (inputColumns(iColumn),*) mergerTrees%spin                            (1,iNode)
             case (propertyTypeSpinY                   %ID)
                ! Column is y spin.
                read (inputColumns(iColumn),*) mergerTrees%spin                            (2,iNode)
             case (propertyTypeSpinZ                   %ID)
                ! Column is z spin.
                read (inputColumns(iColumn),*) mergerTrees%spin                            (3,iNode)
             case (propertyTypeSpin                    %ID)
                ! Column is scalar spin.
                read (inputColumns(iColumn),*) mergerTrees%spinMagnitude                   (  iNode)
             case (propertyTypeAngularMomentumX        %ID)
                ! Column is x angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%angularMomentum                 (1,iNode)
             case (propertyTypeAngularMomentumY        %ID)
                ! Column is y angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%angularMomentum                 (2,iNode)
             case (propertyTypeAngularMomentumZ        %ID)
                ! Column is z angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%angularMomentum                 (3,iNode)
             case (propertyTypeAngularMomentum         %ID)
                ! Column is scalar angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%angularMomentumMagnitude        (  iNode)
             case (propertyTypeSpecificAngularMomentumX%ID)
                ! Column is x specific angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%specificAngularMomentum         (1,iNode)
             case (propertyTypeSpecificAngularMomentumY%ID)
                ! Column is y specific angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%specificAngularMomentum         (2,iNode)
             case (propertyTypeSpecificAngularMomentumZ%ID)
                ! Column is z specific angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%specificAngularMomentum         (3,iNode)
             case (propertyTypeSpecificAngularMomentum %ID)
                ! Column is scalar specific angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%specificAngularMomentumMagnitude(  iNode)
             case (propertyTypeHalfMassRadius          %ID)
                ! Column is half mass radius.
                read (inputColumns(iColumn),*) mergerTrees%halfMassRadius                  (  iNode)
             case (propertyTypeScaleRadius             %ID)
                ! Column is scale radius.
                read (inputColumns(iColumn),*) mergerTrees%scaleRadius                     (  iNode)
             case (propertyTypeMostBoundParticleIndex  %ID)
                ! Column is a most bound particle index.
                read (inputColumns(iColumn),*) mergerTrees%mostBoundParticleIndex          (  iNode)
             case (propertyTypeSnapshot                %ID)
                ! Column is a snapshot index.
                read (inputColumns(iColumn),*) mergerTrees%snapshot                        (  iNode)
             case (propertyTypeVelocityMaximum         %ID)
                ! Column is a maximum velocity.
                read (inputColumns(iColumn),*) mergerTrees%velocityMaximum                 (  iNode)
             case (propertyTypeVelocityDispersion      %ID)
                ! Column is a velocity dispersion.
                read (inputColumns(iColumn),*) mergerTrees%velocityDispersion              (  iNode)
             case default
                call Error_Report('unknown column type'//{introspection:location})
             end select
          end do
       end if
    end do
    close(fileUnit)

    ! Report number of forests found.
    message='Found '
    message=message//mergerTrees%forestCount//' forests'
    call displayMessage(message)

    ! Deallocate workspace.
    if (allocated(inputColumns)) deallocate(inputColumns)

    ! If we have the particle mass, set the masses of any subhalos (which have zero mass by default) based on particle count.
    call Merger_Tree_Data_Set_Subhalo_Masses(mergerTrees)

    ! Convert specific angular momenta as needed.
    if (mergerTrees%hasSpecificAngularMomentumMagnitude.and..not.mergerTrees%hasAngularMomentumMagnitude) then
       allocate(mergerTrees%angularMomentumMagnitude        (mergerTrees%nodeCount))
       mergerTrees%angularMomentumMagnitude=mergerTrees%specificAngularMomentumMagnitude*mergerTrees%nodeMass
       deallocate(mergerTrees%specificAngularMomentumMagnitude                          )
       mergerTrees%hasSpecificAngularMomentumMagnitude=.false.
       mergerTrees%hasAngularMomentumMagnitude        =.true.
    end if
    if (mergerTrees%hasSpecificAngularMomentumX.and..not.mergerTrees%hasAngularMomentumX) then
       allocate(mergerTrees%angularMomentum        (3,mergerTrees%nodeCount))
       forall(i=1:3)
          mergerTrees%angularMomentum(i,:)=mergerTrees%specificAngularMomentum(i,:)*mergerTrees%nodeMass
       end forall
       deallocate(mergerTrees%specificAngularMomentum                          )
       mergerTrees%hasSpecificAngularMomentumX=.false.
       mergerTrees%hasSpecificAngularMomentumY=.false.
       mergerTrees%hasSpecificAngularMomentumZ=.false.
       mergerTrees%hasAngularMomentumX        =.true.
       mergerTrees%hasAngularMomentumY        =.true.
       mergerTrees%hasAngularMomentumZ        =.true.
    end if
    
    ! If needed convert host IDs of self-hosting halos.
    if (mergerTrees%hasHostIndex .and. mergerTrees%hasDummyHostId) then
       where (mergerTrees%hostIndex == mergerTrees%dummyHostId)
          mergerTrees%hostIndex=mergerTrees%nodeIndex
       end where
    end if

    ! If no redshift is given convert scale factor to redshift.
    if(mergerTrees%hasScaleFactor.and..not.mergerTrees%hasRedshift) then
        allocate(mergerTrees%redshift(mergerTrees%nodeCount))
        if(present(maximumRedshift)) then
           maximumRedshiftActual=maximumRedshift
        else
           maximumRedshiftActual=maximumRedshiftDefault
        end if
        do iNode=1,mergerTrees%nodeCount
           if (mergerTrees%scaleFactor(iNode) > 0.0d0) then
              mergerTrees%redshift(iNode)=min((1.0d0/mergerTrees%scaleFactor(iNode))-1.0d0,maximumRedshiftActual)
           else
              mergerTrees%redshift(iNode)=                                                 maximumRedshiftActual
           end if
        end do
        deallocate(mergerTrees%scaleFactor)
        mergerTrees%hasScaleFactor=.false.
        mergerTrees%hasRedshift   =.true.
    end if

    ! Convert properties with inconsistent units.
    do i=propertyTypeMin,propertyTypeMax
       if (mergerTrees%convertProperty(i) /= 1.0d0) then
           call Merger_Tree_Data_Structure_Convert_Property_Units(mergerTrees,i,mergerTrees%convertProperty(i))
       end if
    end do

    ! Set tree indices.
    call Merger_Tree_Data_Structure_Set_Tree_Indices(mergerTrees)
    return
  end subroutine Merger_Tree_Data_Structure_Read_ASCII

  subroutine Merger_Tree_Data_Structure_Convert_Property_Units(mergerTrees,propertyType,conversionFactor)
    !!{
    Convert the property with inconsistent units.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (mergerTreeData), intent(inout) :: mergerTrees
    integer                         , intent(in   ) :: propertyType
    double precision                , intent(in   ) :: conversionFactor
    integer                                         :: j

    select case (propertyType)
    case (propertyTypeNodeMass          %ID)
       ! Property is mass.
       mergerTrees%nodeMass                    =mergerTrees%nodeMass                     *conversionFactor
    case (propertyTypeNodeMass200Mean   %ID)
       ! Property is mass.
       mergerTrees%nodeMass200Mean             =mergerTrees%nodeMass200Mean              *conversionFactor
    case (propertyTypeNodeMass200Crit   %ID)
       ! Property is mass.
       mergerTrees%nodeMass200Crit             =mergerTrees%nodeMass200Crit              *conversionFactor
    case (propertyTypePositionX         %ID)
       ! Property is position.
       forall(j=1:3)
          mergerTrees%position            (j,:)=mergerTrees%position                (j,:)*conversionFactor
       end forall
    case (propertyTypeVelocityX         %ID)
       ! Property is velocity.
       forall(j=1:3)
          mergerTrees%velocity            (j,:)=mergerTrees%velocity                (j,:)*conversionFactor
       end forall
    case (propertyTypeAngularMomentumX  %ID)
       ! Property is angular momentum vector.
       forall(j=1:3)
          mergerTrees%angularMomentum     (j,:)=mergerTrees%angularMomentum         (j,:)*conversionFactor
       end forall
    case (propertyTypeAngularMomentum   %ID)
       ! Property is scalar angular momentum.
       mergerTrees%angularMomentumMagnitude    =mergerTrees%angularMomentumMagnitude     *conversionFactor
    case (propertyTypeHalfMassRadius    %ID)
       ! Property is half mass radius.
       mergerTrees%halfMassRadius              =mergerTrees%halfMassRadius               *conversionFactor
    case (propertyTypeScaleRadius       %ID)
       ! Property is scale radius.
       mergerTrees%scaleRadius                 =mergerTrees%scaleRadius                  *conversionFactor
    case (propertyTypeVelocityMaximum   %ID)
       ! Property is maximum velocity.
       mergerTrees%velocityMaximum             =mergerTrees%velocityMaximum              *conversionFactor
    case (propertyTypeVelocityDispersion%ID)
       ! Property is velocity dispersion.
       mergerTrees%velocityDispersion          =mergerTrees%velocityDispersion           *conversionFactor
    case default
       ! Property has no units.
       call Error_Report('property has no units to convert.'//{introspection:location})
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Convert_Property_Units

  subroutine Merger_Tree_Data_Structure_Set_Tree_Indices(mergerTrees)
    !!{
    Set the merger tree index arrays.
    !!}
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                                :: iNode      , iTree

    ! Allocate arrays for tree start and stop indices and reference ID.
    if (allocated(mergerTrees%treeBeginsAt )) deallocate(mergerTrees%treeBeginsAt )
    if (allocated(mergerTrees%treeNodeCount)) deallocate(mergerTrees%treeNodeCount)
    if (allocated(mergerTrees%forestID     )) deallocate(mergerTrees%forestID     )
    if (allocated(mergerTrees%treeWeight   )) deallocate(mergerTrees%treeWeight   )
    allocate(mergerTrees%treeBeginsAt (mergerTrees%forestCount))
    allocate(mergerTrees%treeNodeCount(mergerTrees%forestCount))
    allocate(mergerTrees%forestID     (mergerTrees%forestCount))
    allocate(mergerTrees%treeWeight   (mergerTrees%forestCount))

    ! Determine index in arrays where each tree begins.
    mergerTrees%treeBeginsAt (1)=0
    mergerTrees%forestID     (1)=mergerTrees%forestIndex     (1)
    if (mergerTrees%hasForestWeight) then
       mergerTrees%treeWeight(1)=mergerTrees%forestWeightNode(1)
    else
       mergerTrees%treeWeight(1)=1.0d0
    end if
    iTree                      =1
    do iNode=2,mergerTrees%nodeCount
       if (mergerTrees%forestIndex(iNode) /= mergerTrees%forestIndex(iNode-1)) then
          iTree=iTree+1
          mergerTrees%treeBeginsAt (iTree  )=iNode-1
          mergerTrees%treeNodeCount(iTree-1)=mergerTrees%treeBeginsAt    (iTree)-mergerTrees%treeBeginsAt(iTree-1)
          mergerTrees%forestID     (iTree  )=mergerTrees%forestIndex     (iNode)
          if (mergerTrees%hasForestWeight) then
             mergerTrees%treeWeight(iTree  )=mergerTrees%forestWeightNode(iNode)
          else
             mergerTrees%treeWeight(iTree  )=1.0d0
          end if
       end if
    end do
    mergerTrees%treeNodeCount(mergerTrees%forestCount)=mergerTrees%nodeCount-mergerTrees%treeBeginsAt(mergerTrees%forestCount)
    return
  end subroutine Merger_Tree_Data_Structure_Set_Tree_Indices

  subroutine Merger_Tree_Data_Structure_Read_Particles_ASCII(mergerTrees,inputFile,columnHeaders,commentCharacter,separator)
    !!{
    Read in particle data from an ASCII file.
    !!}
    use :: File_Utilities   , only : Count_Lines_In_File
    use :: Error            , only : Error_Report
    use :: String_Handling  , only : String_Count_Words , String_Split_Words
    implicit none
    class    (mergerTreeData), intent(inout)               :: mergerTrees
    character(len=*         ), intent(in   )               :: inputFile
    character(len=1         ), intent(in   ), optional     :: commentCharacter
    character(len=*         ), intent(in   ), optional     :: separator
    logical                  , intent(in   ), optional     :: columnHeaders
    character(len=32        ), allocatable  , dimension(:) :: inputColumns
    integer                                                :: columnsCount    , fileUnit           , &
         &                                                    iNode           , nodeCount          , &
         &                                                    iColumn
    logical                                                :: gotFirstDataLine, gotColumnHeaderLine
    character(len=1024      )                              :: inputLine

    ! Flag that these trees have particles.
    mergerTrees%hasParticles=.true.

    ! Determine number of particles.
    nodeCount=Count_Lines_In_File(inputFile,commentCharacter)
    if (present(columnHeaders).and.columnHeaders) nodeCount=nodeCount-1
    call mergerTrees%particleCountSet(nodeCount)

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
         & ) call Error_Report("all three axes or none must be supplied for particle position"//{introspection:location})
    if     (.not.((     mergerTrees%hasParticleVelocityX.and.     mergerTrees%hasParticleVelocityY.and.     mergerTrees%hasParticleVelocityZ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasParticleVelocityX.and..not.mergerTrees%hasParticleVelocityY.and..not.mergerTrees%hasParticleVelocityZ) &
         &        )                                                                                                   &
         & ) call Error_Report("all three axes or none must be supplied for particle velocity"//{introspection:location})

    ! Ensure we have a redshift.
    if (.not.mergerTrees%hasParticleRedshift) call Error_Report("particle redshift must be supplied"//{introspection:location})

    ! Deallocate internal arrays.
    if (allocated(mergerTrees%particleIndex   )) deallocate(mergerTrees%particleIndex   )
    if (allocated(mergerTrees%particleRedshift)) deallocate(mergerTrees%particleRedshift)
    if (allocated(mergerTrees%particlePosition)) deallocate(mergerTrees%particlePosition)
    if (allocated(mergerTrees%particleVelocity)) deallocate(mergerTrees%particleVelocity)
    if (allocated(mergerTrees%particleSnapshot)) deallocate(mergerTrees%particleSnapshot)

    ! Allocate internal arrays to correct size as needed.
    if (mergerTrees%hasParticleIndex    ) allocate(mergerTrees%particleIndex   (mergerTrees%particlesCount))
    if (mergerTrees%hasParticleRedshift ) allocate(mergerTrees%particleRedshift(mergerTrees%particlesCount))
    if (mergerTrees%hasParticlePositionX) allocate(mergerTrees%particlePosition(3,mergerTrees%particlesCount))
    if (mergerTrees%hasParticleVelocityX) allocate(mergerTrees%particleVelocity(3,mergerTrees%particlesCount))
    if (mergerTrees%hasParticleSnapshot ) allocate(mergerTrees%particleSnapshot(mergerTrees%particlesCount))

    ! Open the file and read lines.
    open(newunit=fileUnit,file=inputFile,status='old',form='formatted')
    iNode              =0
    columnsCount       =0
    gotFirstDataLine   =.false.
    gotColumnHeaderLine=.not.present(columnHeaders).or..not.columnHeaders
    do while (iNode < nodeCount)
       ! Get the line.
       read (fileUnit,'(a)') inputLine
       ! Check if this is a data line.
       if (.not.present(commentCharacter) .or. inputLine(1:1) /= commentCharacter) then
          ! Skip header line.
          if (.not.gotColumnHeaderLine) then
             gotColumnHeaderLine=.true.
             cycle
          end if
          ! If this is the first data line, determine how many columns are present and allocate array to store them.
          if (.not.gotFirstDataLine) then
             columnsCount=String_Count_Words(inputLine,separator)
             allocate(inputColumns(columnsCount))
             gotFirstDataLine=.true.
          end if
          ! Count nodes.
          iNode=iNode+1
          call String_Split_Words(inputColumns,inputLine,separator)
          do iColumn=1,min(columnsCount,size(mergerTrees%particleColumnProperties))
             select case (mergerTrees%particleColumnProperties(iColumn)%ID)
             case (propertyTypeNull         %ID)
                ! Ignore this column.
             case (propertyTypeParticleIndex%ID)
                ! Column is particle index.
                read (inputColumns(iColumn),*) mergerTrees%particleIndex   (  iNode)
             case (propertyTypeRedshift     %ID)
                ! Column is redshift.
                read (inputColumns(iColumn),*) mergerTrees%particleRedshift(  iNode)
             case (propertyTypeSnapshot     %ID)
                ! Column is snapshot.
                read (inputColumns(iColumn),*) mergerTrees%particleSnapshot(  iNode)
             case (propertyTypePositionX    %ID)
                ! Column is x position.
                read (inputColumns(iColumn),*) mergerTrees%particlePosition(1,iNode)
             case (propertyTypePositionY    %ID)
                ! Column is y position.
                read (inputColumns(iColumn),*) mergerTrees%particlePosition(2,iNode)
             case (propertyTypePositionZ    %ID)
                ! Column is z position.
                read (inputColumns(iColumn),*) mergerTrees%particlePosition(3,iNode)
             case (propertyTypeVelocityX    %ID)
                ! Column is x velocity.
                read (inputColumns(iColumn),*) mergerTrees%particleVelocity(1,iNode)
             case (propertyTypeVelocityY    %ID)
                ! Column is y velocity.
                read (inputColumns(iColumn),*) mergerTrees%particleVelocity(2,iNode)
             case (propertyTypeVelocityZ    %ID)
                ! Column is z velocity.
                read (inputColumns(iColumn),*) mergerTrees%particleVelocity(3,iNode)
             case default
                call Error_Report('unknown column type'//{introspection:location})
             end select
          end do
       end if
    end do
    close(fileUnit)

    ! Deallocate workspace.
    if (allocated(inputColumns)) deallocate(inputColumns)

    return
  end subroutine Merger_Tree_Data_Structure_Read_Particles_ASCII

  subroutine Merger_Tree_Data_Structure_Export(mergerTrees,outputFileName,outputFormat,hdfChunkSize,hdfCompressionLevel,append)
    !!{
    Output a set of merger trees to an HDF5 file.
    !!}
    use :: Error         , only : Error_Report
    use :: HDF5          , only : hsize_t
    use :: File_Utilities, only : lockDescriptor, File_Lock, File_Unlock
    implicit none
    integer  (kind=hsize_t                   ), intent(in   )           :: hdfChunkSize
    integer                                   , intent(in   )           :: hdfCompressionLevel
    type     (enumerationMergerTreeFormatType), intent(in   )           :: outputFormat
    class    (mergerTreeData                 ), intent(inout)           :: mergerTrees
    character(len=*                          ), intent(in   )           :: outputFileName
    logical                                   , intent(in   ), optional :: append
    type     (lockDescriptor                 )                          :: fileLock

    ! Validate the merger tree.
    call Merger_Tree_Data_Validate_Trees            (mergerTrees)

    ! If we have most-bound particle indices and particle data has been read, construct arrays giving position of particle data for each node.
    call Merger_Tree_Data_Construct_Particle_Indices(mergerTrees)

    call File_Lock(outputFileName,fileLock,lockIsShared=.false.)
    select case (outputFormat%ID)
    case (mergerTreeFormatGalacticus%ID)
       call Merger_Tree_Data_Structure_Export_Galacticus(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    case (mergerTreeFormatIrate     %ID)
       call Merger_Tree_Data_Structure_Export_IRATE     (mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    case default
       call Error_Report('output format is not recognized'//{introspection:location})
    end select
    call File_Unlock(fileLock)
    return
  end subroutine Merger_Tree_Data_Structure_Export

  subroutine Merger_Tree_Data_Structure_Export_Galacticus(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    !!{
    Output a set of merger trees to a Galacticus-format HDF5 file.
    !!}
    use :: File_Utilities    , only : File_Exists
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HSIZE_T        , hsize_t
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : assignment(=)  , char
    use :: String_Handling   , only : operator(//)
    implicit none
    integer  (kind=hsize_t )                            , intent(in   ) :: hdfChunkSize
    integer                                             , intent(in   ) :: hdfCompressionLevel
    class    (mergerTreeData)                           , intent(inout) :: mergerTrees
    character(len=*         )                           , intent(in   ) :: outputFileName
    logical                               , optional    , intent(in   ) :: append
    integer  (kind=HSIZE_T  )             , dimension(2)                :: hyperslabCount        , hyperslabStart
    type     (hdf5Object    ), pointer                                  :: attributeGroup
    type     (hdf5Object    ), target                                   :: cosmologyGroup        , genericGroup       , groupFinderGroup, &
         &                                                                 forestHalos           , outputFile         , particlesGroup  , &
         &                                                                 provenanceGroup       , simulationGroup    , treeBuilderGroup, &
         &                                                                 treeDataset           , treeGroup          , forestIndexGroup, &
         &                                                                 forestsGroup          , unitsGroup
    integer                  , allocatable, dimension(:)                :: firstNode             , numberOfNodes
    integer                                                             :: iAttribute            , iProperty          , iTree           , &
         &                                                                 integerAttribute      , completeCount
    type     (varying_string)                                           :: groupName
    logical                                                             :: appendActual          , fileExists

    ! Determine if we are to append to an existing file.
    appendActual=.false.
    if (present(append)) appendActual=append

    ! Determine if file exists.
    fileExists=appendActual.and.File_Exists(outputFileName)

    ! Open the output file.
    !$ call hdf5Access%set()
    outputFile=hdf5Object(outputFileName,overWrite=.not.appendActual,objectsOverwritable=.true.,chunkSize=hdfChunkSize,compressionLevel=hdfCompressionLevel)

    ! Write a format version attribute.
    if (.not.fileExists) call outputFile%writeAttribute(2,"formatVersion")

    ! Create a group for the datasets.
    forestHalos=outputFile%openGroup("forestHalos","Stores all data for merger trees.")

    ! Write the data.
    if (mergerTrees%hasNodeIndex               ) call forestHalos%writeDataset(mergerTrees%nodeIndex               ,"nodeIndex"          ,"The index of each node."                            ,appendTo=appendActual                  )
    if (mergerTrees%hasDescendantIndex         ) call forestHalos%writeDataset(mergerTrees%descendantIndex         ,"descendantIndex"    ,"The index of each descendant node."                 ,appendTo=appendActual                  )
    if (mergerTrees%hasHostIndex               ) call forestHalos%writeDataset(mergerTrees%hostIndex               ,"hostIndex"          ,"The index of each host node."                       ,appendTo=appendActual                  )
    if (mergerTrees%hasNodeMass                ) call forestHalos%writeDataset(mergerTrees%nodeMass                ,"nodeMass"           ,"The mass of each node."                             ,appendTo=appendActual                  )
    if (mergerTrees%hasRedshift                ) call forestHalos%writeDataset(mergerTrees%redshift                ,"redshift"           ,"The redshift of each node."                         ,appendTo=appendActual                  )
    if (mergerTrees%hasNodeMass200Mean         ) call forestHalos%writeDataset(mergerTrees%nodeMass200Mean         ,"nodeMass200Mean"    ,"The M200 mass of each node (200 * mean density)."   ,appendTo=appendActual                  )
    if (mergerTrees%hasNodeMass200Crit         ) call forestHalos%writeDataset(mergerTrees%nodeMass200Crit         ,"nodeMass200Crit"    ,"The M200 mass of each node (200 * crit density)."   ,appendTo=appendActual                  )
    if (mergerTrees%hasScaleFactor             ) call forestHalos%writeDataset(mergerTrees%scaleFactor             ,"scaleFactor"        ,"The scale factor of each node."                     ,appendTo=appendActual                  )
    if (mergerTrees%hasPositionX               ) call forestHalos%writeDataset(mergerTrees%position                ,"position"           ,"The position of each node."                         ,appendTo=appendActual,appendDimension=2)
    if (mergerTrees%hasVelocityX               ) call forestHalos%writeDataset(mergerTrees%velocity                ,"velocity"           ,"The velocity of each node."                         ,appendTo=appendActual,appendDimension=2)
    if (mergerTrees%hasSpinX                   ) call forestHalos%writeDataset(mergerTrees%spin                    ,"spin"               ,"The spin of each node."                             ,appendTo=appendActual                  )
    if (mergerTrees%hasAngularMomentumX        ) call forestHalos%writeDataset(mergerTrees%angularMomentum         ,"angularMomentum"    ,"The angular momentum spin of each node."            ,appendTo=appendActual                  )
    if (mergerTrees%hasSpinMagnitude           ) call forestHalos%writeDataset(mergerTrees%spinMagnitude           ,"spin"               ,"The spin of each node."                             ,appendTo=appendActual                  )
    if (mergerTrees%hasAngularMomentumMagnitude) call forestHalos%writeDataset(mergerTrees%angularMomentumMagnitude,"angularMomentum"    ,"The angular momentum spin of each node."            ,appendTo=appendActual                  )
    if (mergerTrees%hasHalfMassRadius          ) call forestHalos%writeDataset(mergerTrees%halfMassRadius          ,"halfMassRadius"     ,"The half mass radius of each node."                 ,appendTo=appendActual                  )
    if (mergerTrees%hasScaleRadius             ) call forestHalos%writeDataset(mergerTrees%scaleRadius             ,"scaleRadius"        ,"The scale radius of each node."                     ,appendTo=appendActual                  )
    if (mergerTrees%hasVelocityMaximum         ) call forestHalos%writeDataset(mergerTrees%velocityMaximum         ,"velocityMaximum"    ,"The maximum velocity of each node's rotation curve.",appendTo=appendActual                  )
    if (mergerTrees%hasVelocityDispersion      ) call forestHalos%writeDataset(mergerTrees%velocityDispersion      ,"velocityDispersion" ,"The velocity dispersion of each node."              ,appendTo=appendActual                  )
    if (mergerTrees%hasParticleCount           ) call forestHalos%writeDataset(mergerTrees%particleCount           ,"particleCount"      ,"The number of particles in each node."              ,appendTo=appendActual                  )
    if (mergerTrees%hasMostBoundParticleIndex) then
       call forestHalos%writeDataset(mergerTrees%particleReferenceStart,"particleIndexStart","The starting index of particle data for each node.",appendTo=appendActual)
       call forestHalos%writeDataset(mergerTrees%particleReferenceCount,"particleIndexCount","The number of particle data for each node."        ,appendTo=appendActual)
    end if

    ! Begin creating individual merger tree datasets if requested.
    if (mergerTrees%doMakeReferences) then

       ! Create a containing group for individual trees.
       forestsGroup=outputFile%openGroup("forests","Data for individual forests.")

       ! Create groups for trees and dataset references.
       do iTree=1,mergerTrees%forestCount
          groupName="forest"
          groupName=groupName//iTree
          treeGroup=forestsGroup%openGroup(char(groupName),"Data for a forest.")

          ! Standard datasets.
          hyperslabStart(1)=mergerTrees%treeBeginsAt (iTree)
          hyperslabCount(1)=mergerTrees%treeNodeCount(iTree)
          do iProperty=propertyTypeMin,propertyTypeMax
             ! Skip null and tree index cases.
             if     (                                       &
                  &   iProperty == propertyTypeNull     %ID &
                  &  .or.                                   &
                  &   iProperty == propertyTypeTreeIndex%ID &
                  & ) cycle
             ! Skip cases where we have the corresponding 3-D dataset.
             if (iProperty == propertyTypeSpin           %ID .and. .not.mergerTrees%hasSpinMagnitude           ) cycle
             if (iProperty == propertyTypeAngularMomentum%ID .and. .not.mergerTrees%hasAngularMomentumMagnitude) cycle
             if (forestHalos%hasDataset(char(enumerationPropertyTypeDecode(iProperty)))) then
                treeDataset=forestHalos%openDataset(char(enumerationPropertyTypeDecode(iProperty)))
                call treeGroup%createReference1D(treeDataset,char(enumerationPropertyTypeDecode(iProperty)),hyperslabStart,hyperslabCount)
             end if
          end do

          ! Datasets for properties defined in 3-D spaces.
          hyperslabStart(2)=hyperslabStart(1)
          hyperslabCount(2)=hyperslabCount(1)
          hyperslabStart(1)=int(1,kind=HSIZE_T)
          hyperslabCount(1)=int(3,kind=HSIZE_T)
          do iProperty=1,size(propertyNames3D)
             if (forestHalos%hasDataset(trim(propertyNames3D(iProperty)))) then
                treeDataset=forestHalos%openDataset(propertyNames3D(iProperty))
                call treeGroup%createReference2D(treeDataset,trim(propertyNames3D(iProperty)),hyperslabStart,hyperslabCount)
             end if
          end do
       end do
       ! Finished making individual merger tree datasets
    end if

    ! Write particle data if necessary.
    if (mergerTrees%hasParticles) then
       ! Open the particles group.
       particlesGroup=outputFile%openGroup("particles","Data for a particles.")

       ! Write datasets.
       if (mergerTrees%hasParticleIndex    ) call particlesGroup%writeDataset(mergerTrees%particleIndex   ,"particleID","The ID of each particle."      ,appendTo=appendActual)
       if (mergerTrees%hasParticleRedshift ) call particlesGroup%writeDataset(mergerTrees%particleRedshift,"redshift"  ,"The redshift of each particle.",appendTo=appendActual)
       if (mergerTrees%hasParticlePositionX) call particlesGroup%writeDataset(mergerTrees%particlePosition,"position"  ,"The position of each particle.",appendTo=appendActual)
       if (mergerTrees%hasParticleVelocityX) call particlesGroup%writeDataset(mergerTrees%particleVelocity,"velocity"  ,"The velocity of each particle.",appendTo=appendActual)
    end if

    ! Create datasets giving positions of merger trees within the node arrays.
    forestIndexGroup=outputFile%openGroup("forestIndex","Locations of forests within the halo data arrays.")
    if (fileExists) then
       call forestIndexGroup%readDataset("firstNode"    ,firstNode    )
       call forestIndexGroup%readDataset("numberOfNodes",numberOfNodes)
       mergerTrees%treeBeginsAt=mergerTrees%treeBeginsAt+firstNode(size(firstNode))+numberOfNodes(size(numberOfNodes))
       deallocate(firstNode    )
       deallocate(numberOfNodes)
    end if
    call        forestIndexGroup%writeDataset(mergerTrees%treeBeginsAt ,"firstNode"    ,"Position of the first node in each forest in the halo data arrays.",appendTo=appendActual)
    call        forestIndexGroup%writeDataset(mergerTrees%treeNodeCount,"numberOfNodes","Number of nodes in each forest."                                   ,appendTo=appendActual)
    call        forestIndexGroup%writeDataset(mergerTrees%forestID     ,"forestIndex"  ,"Unique index of forest."                                           ,appendTo=appendActual)
    if (mergerTrees%hasForestWeight.or..not.mergerTrees%hasBoxSize)                                                                                                               &
         & call forestIndexGroup%writeDataset(mergerTrees%treeWeight   ,"forestWeight" ,"Weight of forest."                                                 ,appendTo=appendActual)

    ! Only write remaining data if we are not appending to an existing file.
    if (.not.fileExists) then

       ! Set tree metadata.
       ! Determine if trees have subhalos.
       if (mergerTrees%hasHostIndex) then
          ! A host index is included, so potentially there could be subhalos. (A given tree may not actually
          ! have any subhalos, but we take the presence of the host index to indicate that there could be
          ! subhalos, in principle.)
          integerAttribute=1
       else
          ! No host index is included - assume no nodes are subhalos.
          integerAttribute=0
       end if
       call forestHalos%writeAttribute(integerAttribute,"treesHaveSubhalos")
       ! Determine if trees are self-contained.
       if (mergerTrees%areSelfContained) then
          integerAttribute=1
       else
          integerAttribute=0
       end if
       call forestHalos%writeAttribute(integerAttribute,"forestsAreSelfContained")
       ! Determine if velocities include the Hubble flow.
       if (mergerTrees%includesHubbleFlow) then
          integerAttribute=1
       else
          integerAttribute=0
       end if
       call forestHalos%writeAttribute(integerAttribute,"velocitiesIncludeHubbleFlow")
       ! Determine if positions are periodic.
       if (mergerTrees%isPeriodic) then
          integerAttribute=1
       else
          integerAttribute=0
       end if
       call forestHalos%writeAttribute(integerAttribute,"positionsArePeriodic")
       ! Determine if halo masses include subhalo contributions.
       if (mergerTrees%includesSubhaloMasses) then
          integerAttribute=1
       else
          integerAttribute=0
       end if
       call forestHalos%writeAttribute(integerAttribute,"haloMassesIncludeSubhalos")

       ! Store units.
       unitsGroup=outputFile%openGroup("units","The units system used.")
       if (mergerTrees%unitsSet(unitsMass    %ID)) call Store_Unit_Attributes_Galacticus(unitsMass    ,"mass"    ,mergerTrees,unitsGroup)
       if (mergerTrees%unitsSet(unitsLength  %ID)) call Store_Unit_Attributes_Galacticus(unitsLength  ,"length"  ,mergerTrees,unitsGroup)
       if (mergerTrees%unitsSet(unitsTime    %ID)) call Store_Unit_Attributes_Galacticus(unitsTime    ,"time"    ,mergerTrees,unitsGroup)
       if (mergerTrees%unitsSet(unitsVelocity%ID)) call Store_Unit_Attributes_Galacticus(unitsVelocity,"velocity",mergerTrees,unitsGroup)

       ! Create groups for attributes.
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeGeneric    )) genericGroup    =outputFile%openGroup("metaData"   ,"Generic metadata."                  )
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeCosmology  )) cosmologyGroup  =outputFile%openGroup("cosmology"  ,"Cosmological parameters."           )
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeSimulation )) simulationGroup =outputFile%openGroup("simulation" ,"Simulation parameters."             )
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeGroupFinder)) groupFinderGroup=outputFile%openGroup("groupFinder","Group finder parameters."           )
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeTreeBuilder)) treeBuilderGroup=outputFile%openGroup("treeBuilder","Tree building algorithm parameters.")
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeProvenance )) provenanceGroup =outputFile%openGroup("provenance" ,"Data provenance."                   )

       ! Write attributes.
       do iAttribute=1,mergerTrees%metaDataCount

          ! Determine which group to write to.
          select case (mergerTrees%metaData(iAttribute)%metadataType%ID)
          case (metaDataTypeGeneric    %ID)
             attributeGroup => genericGroup
          case (metaDataTypeCosmology  %ID)
             attributeGroup => cosmologyGroup
          case (metaDataTypeSimulation %ID)
             attributeGroup => simulationGroup
          case (metaDataTypeGroupFinder%ID)
             attributeGroup => groupFinderGroup
          case (metaDataTypeTreeBuilder%ID)
             attributeGroup => treeBuilderGroup
          case (metaDataTypeProvenance %ID)
             attributeGroup => provenanceGroup
          case default
             attributeGroup => null()
             call Error_Report('unknown meta-data group'//{introspection:location})
          end select

          ! Determine what data type to write.
          select case (mergerTrees%metaData(iAttribute)%dataType%ID)
          case (dataTypeInteger%ID)
             call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%integerAttribute,char(mergerTrees%metaData(iAttribute)%label))
          case (dataTypeDouble %ID)
             call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%doubleAttribute ,char(mergerTrees%metaData(iAttribute)%label))
          case (dataTypeText   %ID)
             call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%textAttribute   ,char(mergerTrees%metaData(iAttribute)%label))
          end select
       end do
    end if

    ! Add a flag to indicate successfully completed writing merger tree information.
    if (fileExists) then
       call outputFile%readAttribute('fileCompleteFlag',completeCount)
       completeCount=completeCount+1
    else
       completeCount=             +1
    end if
    call outputFile%writeAttribute(completeCount,"fileCompleteFlag")
    !$ call hdf5Access%unset()

    return
  end subroutine Merger_Tree_Data_Structure_Export_Galacticus

  subroutine Merger_Tree_Data_Structure_Export_IRATE(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    !!{
    Output a set of merger trees to an IRATE-format HDF5 file.
    !!}
    use :: Array_Utilities   , only : Array_Index  , Array_Which
    use :: File_Utilities    , only : File_Exists
    use :: Error             , only : Error_Report
    use :: HDF5              , only : hsize_t
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : assignment(=), char
    implicit none
    integer         (kind=hsize_t  )                           , intent(in   ) ::        hdfChunkSize
    integer                                                    , intent(in   ) ::        hdfCompressionLevel
    class           (mergerTreeData)                           , intent(inout) ::        mergerTrees
    character       (len=*         )                           , intent(in   ) ::        outputFileName
    logical                                                    , intent(in   ) , optional::                      append
    type            (hdf5Object    ), pointer                                  ::        attributeGroup
    type            (hdf5Object    ), target                                   ::        cosmologyGroup                , darkParticlesGroup  , &
         &                                                                               haloTrees                     , mergerTreesGroup    , &
         &                                                                               outputFile                    , particlesGroup      , &
         &                                                                               simulationGroup               , snapshotGroup       , &
         &                                                                               dataset
    integer                         , allocatable, dimension(:)                ::        nodeSnapshotIndices           , snapshotIndices
    integer         (c_size_t      ), allocatable, dimension(:)                ::        descendantSnapshot
    double precision                , allocatable, dimension(:)                ::        particleMass
    integer                                                                    ::        iAttribute                    , nodesOnSnapshotCount, &
         &                                                                               particlesOnSnapshotCount
    integer         (c_size_t      )                                           ::        iDescendant                   , iNode               , &
         &                                                                               iSnapshot                     , snapshotMaximum     , &
         &                                                                               snapshotMinimum
    character       (len=14        )                                           ::        snapshotGroupName
    logical                                                                    ::        appendActual                  , fileExists

    ! Determine if we are to append to an existing file.
    appendActual=.false.
    if (present(append)) appendActual=append

    ! Determine if file exists.
    fileExists=appendActual.and.File_Exists(outputFileName)

    ! IRATE-specific validation.
    if (.not.mergerTrees%hasSnapshot       ) call Error_Report('snapshot indices are required for this format'  //{introspection:location})
    if (.not.mergerTrees%hasPositionX      ) call Error_Report('halo positions are required for this format'    //{introspection:location})
    if (.not.mergerTrees%hasNodeIndex      ) call Error_Report('halo indices are required for this format'      //{introspection:location})
    if (.not.mergerTrees%hasDescendantIndex) call Error_Report('descendant indices are required for this format'//{introspection:location})

    ! Open the output file.
    !$ call hdf5Access%set()
    outputFile=hdf5Object(outputFileName,overWrite=.not.appendActual,chunkSize=hdfChunkSize,compressionLevel=hdfCompressionLevel)

    ! Write the IRATE version.
    if (.not.fileExists) call outputFile%writeAttribute(0,"IRATEVersion")

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
       allocate(snapshotIndices(nodesOnSnapshotCount))
       call Array_Which(mergerTrees%snapshot == iSnapshot,snapshotIndices)

       ! Write redshift attribute.
       if (.not.fileExists) call snapshotGroup%writeAttribute(mergerTrees%redshift(snapshotIndices(1)),"Redshift")

       ! Write the data.
       if (mergerTrees%hasNodeIndex               ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeIndex               ,snapshotIndices),"Index"          ,"The index of each halo."                                                     ,appendTo=appendActual                  )
       end if
       if (mergerTrees%hasNodeMass                ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeMass                ,snapshotIndices),"Mass"           ,"The mass of each halo."                          ,datasetReturned=dataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,dataset)
       end if
       if (mergerTrees%hasNodeMass200Mean         ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeMass200Mean         ,snapshotIndices),"Mass200Mean"    ,"The M200 mass of each halo (200 * mean density).",datasetReturned=dataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,dataset)
       end if
       if (mergerTrees%hasNodeMass200Crit         ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeMass200Crit         ,snapshotIndices),"Mass200Crit"    ,"The M200 mass of each halo (200 * mean density).",datasetReturned=dataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,dataset)
       end if
       if (mergerTrees%hasPositionX               ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%position    ,snapshotIndices,indexOn=2),"Center"         ,"The position of each halo center."                 ,datasetReturned=dataset,appendTo=appendActual,appendDimension=2)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([          unitsLength              ],mergerTrees,dataset)
       end if
       if (mergerTrees%hasVelocityX               ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%velocity   ,snapshotIndices,indexOn=2),"Velocity"       ,"The velocity of each halo."                         ,datasetReturned=dataset,appendTo=appendActual,appendDimension=2)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsVelocity                      ],mergerTrees,dataset)
       end if
       if (mergerTrees%hasSpinX                   ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%spin                    ,snapshotIndices),"Spin"           ,"The spin of each halo."                                                      ,appendTo=appendActual                  )
       end if
       if (mergerTrees%hasAngularMomentumX        ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%angularMomentum         ,snapshotIndices),"AngularMomentum","The angular momentum spin of each halo."         ,datasetReturned=dataset,appendTo=appendActual,appendDimension=2)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass,unitsLength,unitsVelocity],mergerTrees,dataset)
       end if
       if (mergerTrees%hasSpinMagnitude           ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%spinMagnitude           ,snapshotIndices),"Spin"           ,"The spin of each halo."                                                      ,appendTo=appendActual                  )
       end if
       if (mergerTrees%hasAngularMomentumMagnitude) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%angularMomentumMagnitude,snapshotIndices),"AngularMomentum","The angular momentum spin of each halo."         ,datasetReturned=dataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass,unitsLength,unitsVelocity],mergerTrees,dataset)
       end if
       if (mergerTrees%hasHalfMassRadius          ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%halfMassRadius          ,snapshotIndices),"HalfMassRadius" ,"The half mass radius of each halo."              ,datasetReturned=dataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,dataset)
       end if
       ! Destroy snapshot indices.
       deallocate(snapshotIndices)
    end do

    ! Create a group for merger trees.
    mergerTreesGroup=outputFile%openGroup("MergerTrees","Stores all data for merger trees.")

    ! Specify the name of the halo catalog group.
    if (.not.fileExists) call mergerTreesGroup%writeAttribute("HaloCatalog","HaloCatalogName")

    ! Build snapshot numbers for descendants.
    allocate(descendantSnapshot(size(mergerTrees%nodeIndex)))
    descendantSnapshot=-1
    do iDescendant=1,size(mergerTrees%nodeIndex)
       if (mergerTrees%descendantIndex(iDescendant) >= 0) then
          iNode=0
          do while (iNode < size(mergerTrees%nodeIndex) .and. descendantSnapshot(iDescendant) < 0)
             iNode=iNode+1
             if (mergerTrees%nodeIndex(iNode) == mergerTrees%descendantIndex(iDescendant))&
                  & descendantSnapshot(iDescendant)=mergerTrees%snapshot(iNode)
          end do
       end if
    end do

    ! Output merger tree datasets.
    call                               mergerTreesGroup%writeDataset(mergerTrees%snapshot          ,"HaloSnapshot"      ,"The snapshot of each halo."           ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(mergerTrees%nodeIndex         ,"HaloID"            ,"The index of each halo."              ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(mergerTrees%descendantIndex   ,"DescendentID"      ,"The index of each descendant halo."   ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(            descendantSnapshot,"DescendentSnapshot","The snapshot of each descendant halo.",appendTo=appendActual)
    if (mergerTrees%hasHostIndex) call mergerTreesGroup%writeDataset(mergerTrees%hostIndex         ,"HostID"            ,"The index of each host halo."         ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(mergerTrees%treeNodeCount     ,"HalosPerTree"      ,"Number of halos in each tree."        ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(mergerTrees%forestID          ,"TreeID"            ,"Unique index of tree."                ,appendTo=appendActual)
    deallocate(descendantSnapshot)

    if (mergerTrees%hasMostBoundParticleIndex) then
       ! Find the highest and lowest snapshot numbers in the particles.
       if (.not.mergerTrees%hasParticleSnapshot) call Error_Report('particle snapshot numbers must be available for IRATE format export'//{introspection:location})
       snapshotMinimum=minval(mergerTrees%particleSnapshot)
       snapshotMaximum=maxval(mergerTrees%particleSnapshot)

       ! Loop over snapshots to output.
       do iSnapshot=snapshotMinimum,snapshotMaximum

          ! Find those particles which exist at this snapshot.
          particlesOnSnapshotCount=count(mergerTrees%particleSnapshot == iSnapshot)
          allocate(snapshotIndices(particlesOnSnapshotCount))
          call Array_Which(mergerTrees%particleSnapshot == iSnapshot,snapshotIndices)

          ! Find those nodes which exist at this snapshot.
          nodesOnSnapshotCount=count(mergerTrees%snapshot == iSnapshot)
          allocate(nodeSnapshotIndices(nodesOnSnapshotCount))
          call Array_Which(mergerTrees%snapshot == iSnapshot,nodeSnapshotIndices)

          ! Create a snapshot group.
          write (snapshotGroupName,'("Snapshot",i5.5)') iSnapshot
          snapshotGroup=outputFile%openGroup(trim(snapshotGroupName),"Stores all data for a snapshot.")

          ! Create a group for halo catalogs.
          haloTrees=snapshotGroup%openGroup("HaloCatalog","Stores all data for halo catalogs.")

          ! Write the particle indices.
          call haloTrees%writeDataset(Array_Index(mergerTrees%mostBoundParticleIndex,nodeSnapshotIndices),"MostBoundParticleID","The index of each particle.",appendTo=appendActual)

          ! Create a group for particles.
          particlesGroup=snapshotGroup%openGroup("ParticleData","Stores all data for particles.")

          ! Make a group for dark matter particles.
          darkParticlesGroup=particlesGroup%openGroup("Dark","Stores all data for dark matter particles.")

          ! Write redshift attribute.
          if (.not.snapshotGroup%hasAttribute("Redshift")) call snapshotGroup%writeAttribute(mergerTrees%particleRedshift(snapshotIndices(1)),"Redshift")

          ! Write the data.
          allocate(particleMass(particlesOnSnapshotCount))
          particleMass=mergerTrees%particleMass
          call darkParticlesGroup%writeDataset(particleMass,"Mass","The mass of each particle.",datasetReturned=dataset,appendTo=appendActual)
          if (.not.fileExists) call Store_Unit_Attributes_IRATE([unitsMass],mergerTrees,dataset)
          deallocate(particleMass)
          if (mergerTrees%hasParticleIndex    ) then
             call darkParticlesGroup%writeDataset(Array_Index(mergerTrees%particleIndex   ,snapshotIndices),"ID"      ,"The index of each particle."                               ,appendTo=appendActual)
          end if
          if (mergerTrees%hasParticlePositionX) then
             call darkParticlesGroup%writeDataset(Array_Index(mergerTrees%particlePosition,snapshotIndices),"Position","The position of each particle.",datasetReturned=dataset,appendTo=appendActual)
             if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsLength  ],mergerTrees,dataset)
          end if
          if (mergerTrees%hasParticleVelocityX) then
             call darkParticlesGroup%writeDataset(Array_Index(mergerTrees%particleVelocity,snapshotIndices),"Velocity","The velocity of each particle.",datasetReturned=dataset,appendTo=appendActual)
             if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsVelocity],mergerTrees,dataset)
          end if
          ! Destroy the snapshot indices.
          deallocate(snapshotIndices    )
          deallocate(nodeSnapshotIndices)
       end do

    end if

    ! Create groups for attributes.
    if (.not.fileExists) then
       cosmologyGroup  =outputFile%openGroup("Cosmology"            ,"Cosmological parameters."           )
       simulationGroup =outputFile%openGroup("SimulationProperties" ,"Simulation parameters."             )

       ! Write attributes.
       do iAttribute=1,mergerTrees%metaDataCount

          ! Determine which group to write to.
          attributeGroup => null()
          select case (mergerTrees%metaData(iAttribute)%metadataType%ID)
          case (metaDataTypeCosmology %ID)
             attributeGroup => cosmologyGroup
          case (metaDataTypeSimulation%ID)
             attributeGroup => simulationGroup
          end select

          ! Check if the group was recognized.
          if (associated(attributeGroup)) then

             ! Perform dictionary mapping from our internal names (which follow Galacticus format) to IRATE names.
             select case (mergerTrees%metaData(iAttribute)%metadataType%ID)
             case (metaDataTypeCosmology%ID)
                select case (char(mergerTrees%metaData(iAttribute)%label))
                case ("powerSpectrumIndex")
                   mergerTrees%metaData(iAttribute)%label="PowerSpectrumIndex"
                end select
             end select

             ! Determine what data type to write.
             select case (mergerTrees%metaData(iAttribute)%dataType%ID)
             case (dataTypeInteger%ID)
                call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%integerAttribute,char(mergerTrees%metaData(iAttribute)%label))
             case (dataTypeDouble %ID)
                call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%doubleAttribute ,char(mergerTrees%metaData(iAttribute)%label))
             case (dataTypeText   %ID)
                call attributeGroup%writeAttribute(mergerTrees%metaData(iAttribute)%textAttribute   ,char(mergerTrees%metaData(iAttribute)%label))
             end select

          end if

       end do
    end if
    !$ call hdf5Access%unset()
    return
  end subroutine Merger_Tree_Data_Structure_Export_IRATE

  subroutine Store_Unit_Attributes_Galacticus(unitType,unitLabel,mergerTrees,unitsGroup)
    !!{
    Store attributes describing the unit system.
    !!}
    use :: IO_HDF5, only : hdf5Object
    implicit none
    type     (enumerationUnitsType), intent(in   ) :: unitType
    character(len=*               ), intent(in   ) :: unitLabel
    class    (mergerTreeData      ), intent(in   ) :: mergerTrees
    type     (hdf5Object          ), intent(inout) :: unitsGroup

    call unitsGroup%writeAttribute(mergerTrees%units(unitType%ID)%unitsInSI          ,unitLabel//"UnitsInSI"          )
    call unitsGroup%writeAttribute(mergerTrees%units(unitType%ID)%hubbleExponent     ,unitLabel//"HubbleExponent"     )
    call unitsGroup%writeAttribute(mergerTrees%units(unitType%ID)%scaleFactorExponent,unitLabel//"ScaleFactorExponent")
    return
  end subroutine Store_Unit_Attributes_Galacticus

  subroutine Store_Unit_Attributes_IRATE(unitType,mergerTrees,dataset)
    !!{
    Store unit attributes in IRATE format files.
    !!}
    use :: IO_HDF5                     , only : hdf5Object
    use :: ISO_Varying_String          , only : assignment(=), operator(//)
    use :: Numerical_Constants_Prefixes, only : hecto        , kilo
    implicit none
    type            (enumerationUnitsType), dimension(:), intent(in   ) :: unitType
    class           (mergerTreeData      )              , intent(in   ) :: mergerTrees
    type            (hdf5Object          )              , intent(inout) :: dataset
    integer                                                             :: iUnit
    double precision                                                    :: cgsUnits   , hubbleExponent, scaleFactorExponent
    type            (varying_string      )                              :: unitName

    ! Get conversion factor to cgs units.
    cgsUnits           =1.0d0
    hubbleExponent     =0.0d0
    scaleFactorExponent=0.0d0
    unitName           =""
    do iUnit=1,size(unitType)
       cgsUnits=cgsUnits*mergerTrees%units(unitType(iUnit)%ID)%unitsInSI
       select case (unitType(iUnit)%ID)
       case (unitsMass    %ID)
          cgsUnits=cgsUnits*kilo
       case (unitsLength  %ID)
          cgsUnits=cgsUnits*hecto
       case (unitsTime    %ID)
          cgsUnits=cgsUnits*1.0d0
       case (unitsVelocity%ID)
          cgsUnits=cgsUnits*hecto
       end select
       hubbleExponent     =hubbleExponent     +dble(mergerTrees%units(unitType(iUnit)%ID)%hubbleExponent     )
       scaleFactorExponent=scaleFactorExponent+dble(mergerTrees%units(unitType(iUnit)%ID)%scaleFactorExponent)
       if (iUnit > 1) unitName=unitName//" * "
       unitName=unitName//mergerTrees%units(unitType(iUnit)%ID)%name
    end do

    ! Write the attributes.
    call dataset%writeAttribute(                      &
         &                      [                     &
         &                       cgsUnits           , &
         &                       hubbleExponent     , &
         &                       scaleFactorExponent  &
         &                      ],                    &
         &                      "unitscgs"            &
         &                     )
    call dataset%writeAttribute(                      &
         &                       unitName,            &
         &                      "unitname"            &
         &                     )
    return
  end subroutine Store_Unit_Attributes_IRATE

  subroutine Merger_Tree_Data_Validate_Trees(mergerTrees)
    !!{
    Validate the merger trees.
    !!}
    use :: Error, only : Error_Report
    class(mergerTreeData), intent(in   ) :: mergerTrees

    if (.not.mergerTrees%hasForestIndex    ) call Error_Report("merger trees do not have required property 'forestIndex'"    //{introspection:location})
    if (.not.mergerTrees%hasNodeIndex      ) call Error_Report("merger trees do not have required property 'nodeIndex'"      //{introspection:location})
    if (.not.mergerTrees%hasDescendantIndex) call Error_Report("merger trees do not have required property 'descendantIndex'"//{introspection:location})
    if (.not.mergerTrees%hasRedshift       ) call Error_Report("merger trees do not have required property 'redshift'"       //{introspection:location})
    if (.not.mergerTrees%hasNodeMass       ) call Error_Report("merger trees do not have required property 'nodeMass'"       //{introspection:location})
    return
  end subroutine Merger_Tree_Data_Validate_Trees

  subroutine Merger_Tree_Data_Set_Subhalo_Masses(mergerTrees)
    !!{
    Set the masses of any subhalos (which have zero mass by default) based on particle count.
    !!}
    use :: Display, only : displayMagenta, displayReset
    use :: Error  , only : Warn
    class(mergerTreeData), intent(inout) :: mergerTrees

    if (mergerTrees%hasParticleCount) then
       where (mergerTrees%nodeMass <= 0.0d0)
          mergerTrees%nodeMass=dble(mergerTrees%particleCount)*mergerTrees%particleMass
       end where
    end if
    if (any(mergerTrees%nodeMass <= 0.0d0)) call Warn(displayMagenta()//"WARNING:"//displayReset()//" some nodes have non-positive mass"//{introspection:location})
    return
  end subroutine Merger_Tree_Data_Set_Subhalo_Masses

  subroutine Merger_Tree_Data_Construct_Particle_Indices(mergerTrees)
    !!{
    If we have most-bound particle indices and particle data has been read, construct arrays giving position of particle data for each node.
    !!}
    use :: Error            , only : Error_Report
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                                :: foundParticleData
    integer                                :: iNode            , iParticle

    if (mergerTrees%hasMostBoundParticleIndex) then
       ! Insist on having particle data.
       if (.not.mergerTrees%hasParticles) call Error_Report("most bound particle IDs provided, but no particle data was read"//{introspection:location})
       ! Allocate arrays for storing indices.
       if (allocated(mergerTrees%particleReferenceStart)) deallocate(mergerTrees%particleReferenceStart)
       if (allocated(mergerTrees%particleReferenceCount)) deallocate(mergerTrees%particleReferenceCount)
       allocate(mergerTrees%particleReferenceStart(mergerTrees%nodeCount))
       allocate(mergerTrees%particleReferenceCount(mergerTrees%nodeCount))
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
