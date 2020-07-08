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

!+ Contributions to this file made by: Stephanie Dörschner.

!% Contains a module which implements an object to store merger tree data for processing into \glc's preferred file format.

module Merger_Tree_Data_Structure
  !% Implements an object to store merger tree data for processing into \glc's preferred file format.
  use, intrinsic :: ISO_C_Binding     , only : c_size_t
  use            :: ISO_Varying_String, only : varying_string
  use            :: Kind_Numbers      , only : kind_int8
  implicit none
  private
  public :: mergerTreeData

  ! Output formats.
  !# <enumeration>
  !#  <name>mergerTreeFormat</name>
  !#  <description>Used to specify which output format to use for merger tree data.</description>
  !#  <visibility>public</visibility>
  !#  <validator>yes</validator>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <entry label="galacticus" />
  !#  <entry label="irate"      />
  !# </enumeration>

  ! Property labels.
  !# <enumeration>
  !#  <name>propertyType</name>
  !#  <description>Used to specify properties in a {\normalfont \ttfamily mergerTreeData} structure.</description>
  !#  <visibility>public</visibility>
  !#  <validator>yes</validator>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <entry label="null"                    />
  !#  <entry label="treeIndex"               />
  !#  <entry label="nodeIndex"               />
  !#  <entry label="descendentIndex"         />
  !#  <entry label="hostIndex"               />
  !#  <entry label="redshift"                />
  !#  <entry label="scaleFactor"             />
  !#  <entry label="nodeMass"                />
  !#  <entry label="nodeMass200Mean"         />
  !#  <entry label="nodeMass200Crit"         />
  !#  <entry label="particleCount"           />
  !#  <entry label="positionX"               />
  !#  <entry label="positionY"               />
  !#  <entry label="positionZ"               />
  !#  <entry label="velocityX"               />
  !#  <entry label="velocityY"               />
  !#  <entry label="velocityZ"               />
  !#  <entry label="spinX"                   />
  !#  <entry label="spinY"                   />
  !#  <entry label="spinZ"                   />
  !#  <entry label="spin"                    />
  !#  <entry label="angularMomentumX"        />
  !#  <entry label="angularMomentumY"        />
  !#  <entry label="angularMomentumZ"        />
  !#  <entry label="angularMomentum"         />
  !#  <entry label="specificAngularMomentumX"/>
  !#  <entry label="specificAngularMomentumY"/>
  !#  <entry label="specificAngularMomentumZ"/>
  !#  <entry label="specificAngularMomentum" />
  !#  <entry label="halfMassRadius"          />
  !#  <entry label="scaleRadius"             />
  !#  <entry label="particleIndex"           />
  !#  <entry label="mostBoundParticleIndex"  />
  !#  <entry label="snapshot"                />
  !#  <entry label="treeWeight"              />
  !#  <entry label="velocityMaximum"         />
  !#  <entry label="velocityDispersion"      />
  !# </enumeration>

  ! Names of 3-D datasets (i.e. those which give properties in 3-D space).
  character(len=*), parameter :: propertyNames3D(5)=[ &
       & 'position               ',                   &
       & 'velocity               ',                   &
       & 'spin                   ',                   &
       & 'angularMomentum        ',                   &
       & 'specificAngularMomentum'                    &
       &                                            ]

  ! Units labels.
  !# <enumeration>
  !#  <name>units</name>
  !#  <description>Used to specify the type of units being stored in a {\normalfont \ttfamily mergerTreeData} structure.</description>
  !#  <visibility>public</visibility>
  !#  <validator>yes</validator>
  !#  <entry label="mass"     />
  !#  <entry label="length"   />
  !#  <entry label="time"     />
  !#  <entry label="velocity" />
  !# </enumeration>

  type unitsMetaData
     !% A structure that holds metadata on units used.
     double precision                              :: unitsInSI
     integer                                       :: hubbleExponent, scaleFactorExponent
     type            (varying_string), allocatable :: name
  end type unitsMetaData

  ! Metadata labels.
  !# <enumeration>
  !#  <name>metaDataType</name>
  !#  <description>Used to specify the type of metadata being stored in a {\normalfont \ttfamily mergerTreeData} structure.</description>
  !#  <visibility>public</visibility>
  !#  <validator>yes</validator>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <entry label="generic"     />
  !#  <entry label="cosmology"   />
  !#  <entry label="simulation"  />
  !#  <entry label="groupFinder" />
  !#  <entry label="treeBuilder" />
  !#  <entry label="provenance"  />
  !# </enumeration>

  ! Data types for metadata.
  !# <enumeration>
  !#  <name>dataType</name>
  !#  <description>Used to specify the type of data being stored in a {\normalfont \ttfamily mergerTreeData} structure metadata entry.</description>
  !#  <visibility>public</visibility>
  !#  <entry label="null"    />
  !#  <entry label="integer" />
  !#  <entry label="double"  />
  !#  <entry label="text"    />
  !# </enumeration>

  type treeMetaData
     !% Structure that holds metadata for the trees.
     integer                          :: metadataType
     type            (varying_string) :: label
     integer                          :: dataType
     integer                          :: integerAttribute
     double precision                 :: doubleAttribute
     type            (varying_string) :: textAttribute
  end type treeMetaData

  type mergerTreeData
     !% A structure that holds raw merger tree data.
     private
     integer                                                                 :: dummyHostId                        , nodeCount                        , &
          &                                                                     particlesCount                     , forestCount
     double precision                                                        :: particleMass               =0.0d0
     integer                         , allocatable, dimension(:)             :: columnProperties                   , particleColumnProperties         , &
          &                                                                     treeBeginsAt                       , treeNodeCount
     integer         (kind=kind_int8), allocatable, dimension(:)             :: descendentIndex                    , hostIndex                        , &
          &                                                                     mostBoundParticleIndex             , nodeIndex                        , &
          &                                                                     particleIndex                      , forestID                         , &
          &                                                                     forestIndex
     integer         (c_size_t      ), allocatable, dimension(:)             :: particleCount                      , snapshot                         , &
          &                                                                     particleReferenceCount             , particleReferenceStart           , &
          &                                                                     particleSnapshot
     double precision                , allocatable, dimension(:)             :: angularMomentumMagnitude           , halfMassRadius                   , &
          &                                                                     nodeMass                           , particleRedshift                 , &
          &                                                                     redshift                           , scaleFactor                      , &
          &                                                                     scaleRadius                        , spinMagnitude                    , &
          &                                                                     treeWeight                         , forestWeightNode                 , &
          &                                                                     specificAngularMomentumMagnitude   , velocityMaximum                  , &
          &                                                                     velocityDispersion                 , nodeMass200Mean                  , &
          &                                                                     nodeMass200Crit
     double precision                , allocatable, dimension(:,:)           :: angularMomentum                    , particlePosition                 , &
          &                                                                     particleVelocity                   , position                         , &
          &                                                                     spin                               , velocity                         , &
          &                                                                     specificAngularMomentum
     double precision                , dimension(propertyTypeMin:propertyTypeMax)          :: convertProperty            =1.0d0
     logical                                                                 :: hasAngularMomentumMagnitude        , hasAngularMomentumX              , &
          &                                                                     hasAngularMomentumY                , hasAngularMomentumZ              , &
          &                                                                     hasSpecificAngularMomentumMagnitude, hasSpecificAngularMomentumX      , &
          &                                                                     hasSpecificAngularMomentumY        , hasSpecificAngularMomentumZ      , &
          &                                                                     hasDescendentIndex                 , hasDummyHostId                   , &
          &                                                                     hasHalfMassRadius                                                     , &
          &                                                                     hasHostIndex                       , hasMostBoundParticleIndex        , &
          &                                                                     hasNodeIndex                       , hasNodeMass                      , &
          &                                                                     hasNodeMass200Mean                 , hasNodeMass200Crit               , &
          &                                                                     hasParticleCount                   , hasParticleIndex                 , &
          &                                                                     hasParticlePositionX               , hasParticlePositionY             , &
          &                                                                     hasParticlePositionZ               , hasParticleRedshift              , &
          &                                                                     hasParticleSnapshot                , hasParticleVelocityX             , &
          &                                                                     hasParticleVelocityY               , hasParticleVelocityZ             , &
          &                                                                     hasParticles               =.false., hasPositionX                     , &
          &                                                                     hasPositionY                       , hasPositionZ                     , &
          &                                                                     hasRedshift                        , hasScaleFactor                   , &
          &                                                                     hasScaleRadius                     , hasSnapshot                      , &
          &                                                                     hasSpinMagnitude                   , hasSpinX                         , &
          &                                                                     hasSpinY                           , hasSpinZ                         , &
          &                                                                     hasForestIndex                     , hasVelocityX                     , &
          &                                                                     hasVelocityY                       , hasVelocityZ                     , &
          &                                                                     hasVelocityMaximum                 , hasVelocityDispersion            , &
          &                                                                     hasForestWeight                                                       , &
          &                                                                     hasBoxSize
     logical                                                                 :: areSelfContained           =.true. , doMakeReferences         =.true. , &
          &                                                                     includesHubbleFlow         =.false., includesSubhaloMasses    =.false., &
          &                                                                     isPeriodic                 =.false.
     type            (unitsMetaData )             , dimension(unitsMin:unitsMax) :: units
     logical                                      , dimension(unitsMin:unitsMax) :: unitsSet                   =.false.
     integer                                                                 :: metaDataCount              =0
     type            (treeMetaData  ), allocatable, dimension(:)             :: metaData
   contains
     !@ <objectMethods>
     !@   <object>mergerTreeData</object>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <description>Reset the data structure.</description>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>forestCountSet</method>
     !@     <description>Set the total number of forests in the data structure.</description>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ forestCount\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>nodeCountSet</method>
     !@     <description>Set the total number of nodes in the data structure.</description>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ nodeCount\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>particleCountSet</method>
     !@     <description>Set the total number of particles in the data structure.</description>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ particleCount\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>readASCII</method>
     !@     <description>Read node data from an ASCII file into the data structure.</description>
     !@     <type>\void</type>
     !@     <arguments>inputFile\argin, [commentCharacter]\argin, [separator]\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>readParticlesASCII</method>
     !@     <description>Read particle data from an ASCII file into the data structure</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} inputFile\argin, \textcolor{red}{\textless character(len=1)\textgreater} [commentCharacter]\argin, \textcolor{red}{\textless character(len=*)\textgreater} [separator]\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setProperty</method>
     !@     <description>Set a node property in the data structure.</description>
     !@     <type>\void</type>
     !@     <arguments>\enumPropertyType\ propertyType\argin, \textcolor{red}{\textless integer(kind=kind\_int8)(:)|\doubleone\textgreater} property\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setPropertyColumn</method>
     !@     <description>Set the column in an ASCII data file corresponding to a given node property.</description>
     !@     <type>\void</type>
     !@     <arguments>\enumPropertyType\ propertyType\argin, columnNumber\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setParticlePropertyColumn</method>
     !@     <description>Set the column in an ASCII data file corresponding to a given particle property.</description>
     !@     <type>\void</type>
     !@     <arguments>\enumPropertyType\ propertyType\argin, \intzero\ columnNumber\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setParticleMass</method>
     !@     <description>Set the mass of an N-body particle in the simulation from which the trees were derived.</description>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ particleMass\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setSelfContained</method>
     !@     <description>Specify if trees are self-contained (i.e. contain no cross-links to other trees).</description>
     !@     <type>\void</type>
     !@     <arguments>\logicalzero\ areSelfContained\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setIncludesHubbleFlow</method>
     !@     <description>Specify if velocities include the Hubble flow.</description>
     !@     <type>\void</type>
     !@     <arguments>\logicalzero\ includesHubbleFlow\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setPositionsArePeriodic</method>
     !@     <description>Set if positions are periodic.</description>
     !@     <type>\void</type>
     !@     <arguments>\logicalzero\ isPeriodic\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setIncludesSubhaloMasses</method>
     !@     <description>Set whether halo masses include the masses of the subhalos.</description>
     !@     <type>\void</type>
     !@     <arguments>\logicalzero\ includesSubhaloMasses\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setDummyHostId</method>
     !@     <description>Set host ID for self-hosting halos if host ID is not node ID.</description>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ dummyHostId\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setConversionFactor</method>
     !@     <description>Set property type and conversion factor to adjust inconsistent units.</description>
     !@     <type>\void</type>
     !@     <arguments>\enumPropertyType\ propertyType\argin, \doublezero\ conversionFactor\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setUnits</method>
     !@     <description>Set the units used.</description>
     !@     <type>\void</type>
     !@     <arguments>\enumMetaDataType\ unitType\argin, \doublezero\ unitsInSI\argin, \intzero\ [hubbleExponent]\argin, \intzero\ [scaleFactorExponent]\argin, \textcolor{red}{\textless character(len=*)\textgreater} [name]\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>addMetadata</method>
     !@     <description>Add a metadatum to the tree data structure.</description>
     !@     <type>\void</type>
     !@     <arguments>\enumMetaDataType\ metadataType\argin, \textcolor{red}{\textless character(len=*)\textgreater} label\argin, \textcolor{red}{\textless (\doublezero,\intzero,character(len=*))\textgreater} value\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>makeReferences</method>
     !@     <description>Specify whether or not merger tree dataset references should be made.</description>
     !@     <type>\void</type>
     !@     <arguments>\logicalzero\ makeReferences\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>export</method>
     !@     <description>Export the tree data to an output file.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} outputFileName\argin, \enumMergerTreeFormat\ outputFormat\argin, \intzero\ hdfChunkSize\argin, \intzero\ hdfCompressionLevel\argin, \logicalzero\ [append]\argin</arguments>
     !@   </objectMethod>
     !@ </objectMethods>
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
    !% Reset a merger tree data object.
    use :: Memory_Management, only : deallocateArray
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees

    ! No properties.
    mergerTrees%hasForestIndex                     =.false.
    mergerTrees%hasForestWeight                    =.false.
    mergerTrees%hasBoxSize                         =.false.
    mergerTrees%hasNodeIndex                       =.false.
    mergerTrees%hasDescendentIndex                 =.false.
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
    if (allocated(mergerTrees%treeBeginsAt          )) call deallocateArray(mergerTrees%treeBeginsAt          )
    if (allocated(mergerTrees%treeNodeCount         )) call deallocateArray(mergerTrees%treeNodeCount         )
    if (allocated(mergerTrees%forestID              )) call deallocateArray(mergerTrees%forestID              )
    if (allocated(mergerTrees%treeWeight            )) call deallocateArray(mergerTrees%treeWeight            )
    if (allocated(mergerTrees%forestWeightNode      )) call deallocateArray(mergerTrees%forestWeightNode      )
    if (allocated(mergerTrees%forestIndex           )) call deallocateArray(mergerTrees%forestIndex           )
    if (allocated(mergerTrees%nodeIndex             )) call deallocateArray(mergerTrees%nodeIndex             )
    if (allocated(mergerTrees%mostBoundParticleIndex)) call deallocateArray(mergerTrees%mostBoundParticleIndex)
    if (allocated(mergerTrees%descendentIndex       )) call deallocateArray(mergerTrees%descendentIndex       )
    if (allocated(mergerTrees%hostIndex             )) call deallocateArray(mergerTrees%hostIndex             )
    if (allocated(mergerTrees%redshift              )) call deallocateArray(mergerTrees%redshift              )
    if (allocated(mergerTrees%scaleFactor           )) call deallocateArray(mergerTrees%scaleFactor           )
    if (allocated(mergerTrees%nodeMass              )) call deallocateArray(mergerTrees%nodeMass              )
    if (allocated(mergerTrees%nodeMass200Mean       )) call deallocateArray(mergerTrees%nodeMass200Mean       )
    if (allocated(mergerTrees%nodeMass200Crit       )) call deallocateArray(mergerTrees%nodeMass200Crit       )
    if (allocated(mergerTrees%particleCount         )) call deallocateArray(mergerTrees%particleCount         )
    if (allocated(mergerTrees%position              )) call deallocateArray(mergerTrees%position              )
    if (allocated(mergerTrees%velocity              )) call deallocateArray(mergerTrees%velocity              )
    if (allocated(mergerTrees%spin                  )) call deallocateArray(mergerTrees%spin                  )
    if (allocated(mergerTrees%halfMassRadius        )) call deallocateArray(mergerTrees%halfMassRadius        )
    if (allocated(mergerTrees%scaleRadius           )) call deallocateArray(mergerTrees%scaleRadius           )
    if (allocated(mergerTrees%velocityMaximum       )) call deallocateArray(mergerTrees%velocityMaximum       )
    if (allocated(mergerTrees%velocityDispersion    )) call deallocateArray(mergerTrees%velocityDispersion    )
    return
  end subroutine Merger_Tree_Data_Structure_Reset

  subroutine Merger_Tree_Data_Structure_Make_References(mergerTrees,makeReferences)
    !% Specify whether or not to make merger tree dataset references.
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: makeReferences

    mergerTrees%doMakeReferences=makeReferences
    return
  end subroutine Merger_Tree_Data_Structure_Make_References

  subroutine Merger_Tree_Data_Structure_Add_Metadata_Double(mergerTrees,metadataType,label,doubleValue)
    !% Add a double metadatum.
    implicit none
    class           (mergerTreeData), intent(inout) :: mergerTrees
    integer                         , intent(in   ) :: metadataType
    character       (len=*         ), intent(in   ) :: label
    double precision                , intent(in   ) :: doubleValue

    call Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,doubleValue=doubleValue)
    ! Check if this is box size.
    if (metadataType == metaDataTypeSimulation .and. trim(label) == "boxSize") mergerTrees%hasBoxSize=.true.
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata_Double

  subroutine Merger_Tree_Data_Structure_Add_Metadata_Integer(mergerTrees,metadataType,label,integerValue)
    !% Add an integer metadatum.
    implicit none
    class    (mergerTreeData), intent(inout) :: mergerTrees
    integer                  , intent(in   ) :: metadataType
    character(len=*         ), intent(in   ) :: label
    integer                  , intent(in   ) :: integerValue

    call Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,integerValue=integerValue)
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata_Integer

  subroutine Merger_Tree_Data_Structure_Add_Metadata_Text(mergerTrees,metadataType,label,textValue)
    !% Add a double metadatum.
    implicit none
    class    (mergerTreeData), intent(inout) :: mergerTrees
    integer                  , intent(in   ) :: metadataType
    character(len=*         ), intent(in   ) :: label
    character(len=*         ), intent(in   ) :: textValue

    call Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,textValue=textValue)
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata_Text

  subroutine Merger_Tree_Data_Structure_Add_Metadata(mergerTrees,metadataType,label,integerValue,doubleValue,textValue)
    !% Add a metadatum.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    use :: Memory_Management , only : Memory_Usage_Record
    implicit none
    class           (mergerTreeData), intent(inout)               :: mergerTrees
    integer                         , intent(in   )               :: metadataType
    character       (len=*         ), intent(in   )               :: label
    integer                         , intent(in   ), optional     :: integerValue
    double precision                , intent(in   ), optional     :: doubleValue
    character       (len=*         ), intent(in   ), optional     :: textValue
    integer                         , parameter                   :: metadataBlockSize=100
    type            (treeMetaData  ), allocatable  , dimension(:) :: metaDataTemporary

    ! Validate the metadata type.
    if (.not.enumerationMetadataTypeIsValid(metadataType)) call Galacticus_Error_Report('invalid metadata type'//{introspection:location})

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
       if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType /= dataTypeNull) call Galacticus_Error_Report('only one data type can be specified'//{introspection:location})
       mergerTrees%metaData(mergerTrees%metaDataCount)%integerAttribute=integerValue
       mergerTrees%metaData(mergerTrees%metaDataCount)%dataType        =dataTypeInteger
    end if
    if (present(doubleValue)) then
       if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType /= dataTypeNull) call Galacticus_Error_Report('only one data type can be specified'//{introspection:location})
       mergerTrees%metaData(mergerTrees%metaDataCount)%doubleAttribute =doubleValue
       mergerTrees%metaData(mergerTrees%metaDataCount)%dataType        =dataTypeDouble
    end if
    if (present(textValue)) then
       if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType /= dataTypeNull) call Galacticus_Error_Report('only one data type can be specified'//{introspection:location})
       mergerTrees%metaData(mergerTrees%metaDataCount)%textAttribute   =textValue
       mergerTrees%metaData(mergerTrees%metaDataCount)%dataType        =dataTypeText
    else
       mergerTrees%metaData(mergerTrees%metaDataCount)%textAttribute   =""
    end if
    if (mergerTrees%metaData(mergerTrees%metaDataCount)%dataType == dataTypeNull) call Galacticus_Error_Report('no data was given'//{introspection:location})
    return
  end subroutine Merger_Tree_Data_Structure_Add_Metadata

  subroutine Merger_Tree_Data_Structure_Set_Forest_Count(mergerTrees,forestCount)
    !% Set the total number of trees in merger trees.
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                , intent(in   ) :: forestCount

    ! Set the number of trees.
    mergerTrees%forestCount=forestCount

    return
  end subroutine Merger_Tree_Data_Structure_Set_Forest_Count

  subroutine Merger_Tree_Data_Structure_Set_Node_Count(mergerTrees,nodeCount)
    !% Set the total number of nodes in merger trees.
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                , intent(in   ) :: nodeCount

    ! Set the number of nodes.
    mergerTrees%nodeCount=nodeCount

    return
  end subroutine Merger_Tree_Data_Structure_Set_Node_Count

  subroutine Merger_Tree_Data_Structure_Set_Particle_Count(mergerTrees,particleCount)
    !% Set the total number of particles in merger trees.
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                , intent(in   ) :: particleCount

    ! Set the number of nodes.
    mergerTrees%particlesCount=particleCount

    return
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Count

  subroutine Merger_Tree_Data_Structure_Set_Particle_Mass(mergerTrees,particleMass)
    !% Set the particle mass used in the trees.
    implicit none
    class           (mergerTreeData), intent(inout) :: mergerTrees
    double precision                , intent(in   ) :: particleMass

    ! Set the particle mass.
    mergerTrees%particleMass=particleMass

    return
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Mass

  subroutine Merger_Tree_Data_Structure_Set_Self_Contained(mergerTrees,areSelfContained)
    !% Set the particle mass used in the trees.
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: areSelfContained

    ! Set whether trees are self-contained.
    mergerTrees%areSelfContained=areSelfContained

    return
  end subroutine Merger_Tree_Data_Structure_Set_Self_Contained

  subroutine Merger_Tree_Data_Structure_Set_Includes_Hubble_Flow(mergerTrees,includesHubbleFlow)
    !% Set the particle mass used in the trees.
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: includesHubbleFlow

    ! Set whether velocities include the Hubble flow.
    mergerTrees%includesHubbleFlow=includesHubbleFlow

    return
  end subroutine Merger_Tree_Data_Structure_Set_Includes_Hubble_Flow

  subroutine Merger_Tree_Data_Structure_Set_Is_Periodic(mergerTrees,isPeriodic)
    !% Set whether or not positions are periodic.
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: isPeriodic

    ! Set whether positions are periodic.
    mergerTrees%isPeriodic=isPeriodic

    return
  end subroutine Merger_Tree_Data_Structure_Set_Is_Periodic

  subroutine Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses(mergerTrees,includesSubhaloMasses)
    !% Set the particle mass used in the trees.
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                , intent(in   ) :: includesSubhaloMasses

    ! Set whether halo masses include subhalo contributions.
    mergerTrees%includesSubhaloMasses=includesSubhaloMasses

    return
  end subroutine Merger_Tree_Data_Structure_Set_Includes_Subhalo_Masses

  subroutine Merger_Tree_Data_Structure_Set_Self_Hosting_Halo_Id(mergerTrees,dummyHostId)
    !% Set the host ID in case of self-hosting halos. Default is host ID = node ID.
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                , intent(in   ) :: dummyHostId

    ! Set the value of the dummy-variable for self hosting halos.
    mergerTrees%dummyHostId   =dummyHostId
    mergerTrees%hasDummyHostId=.true.
    return
  end subroutine Merger_Tree_Data_Structure_Set_Self_Hosting_Halo_Id

  subroutine Merger_Tree_Data_Structure_Set_Conversion_Factor(mergerTrees,propertyType,conversionFactor)
    !% Set Conversion factor for property type with inconsistent unit.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                , intent(in   ) :: propertyType
    double precision       , intent(in   ) :: conversionFactor

    ! Ensure the property type is valid.
    if (.not.enumerationPropertyTypeIsValid(propertyType)) call Galacticus_Error_Report('invalid property type'//{introspection:location})

    ! Store conversion factor into array.
    mergerTrees%convertProperty(propertyType)=conversionFactor
    return
  end subroutine Merger_Tree_Data_Structure_Set_Conversion_Factor

  subroutine Merger_Tree_Data_Structure_Set_Units(mergerTrees,unitType,unitsInSI,hubbleExponent,scaleFactorExponent,name)
    !% Set the units system.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class           (mergerTreeData), intent(inout)           :: mergerTrees
    integer                         , intent(in   )           :: unitType
    double precision                , intent(in   )           :: unitsInSI
    integer                         , intent(in   ), optional :: hubbleExponent, scaleFactorExponent
    character       (len=*         ), intent(in   ), optional :: name

    ! Ensure the unit type is valid.
    if (.not.enumerationUnitsIsValid(unitType)) call Galacticus_Error_Report('invalid unit type'//{introspection:location})

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
    allocate(mergerTrees%units(unitType)%name)
    if (present(name)) then
       mergerTrees%units(unitType)%name=name
    else
       mergerTrees%units(unitType)%name=""
    end if

    return
  end subroutine Merger_Tree_Data_Structure_Set_Units

  subroutine Merger_Tree_Data_Structure_Set_Particle_Property_Column(mergerTrees,propertyType,columnNumber)
    !% Set column mapping from the input file.
    use :: Memory_Management, only : allocateArray, deallocateArray
    implicit none
    class  (mergerTreeData), intent(inout)               :: mergerTrees
    integer                , intent(in   )               :: columnNumber        , propertyType
    integer                , allocatable  , dimension(:) :: columnPropertiesTemp

    ! Ensure the storage array is large enough.
    if (allocated(mergerTrees%particleColumnProperties)) then
       if (columnNumber > size(mergerTrees%particleColumnProperties)) then
          call Move_Alloc(mergerTrees%particleColumnProperties,columnPropertiesTemp)
          call allocateArray(mergerTrees%particleColumnProperties,[columnNumber])
          mergerTrees%particleColumnProperties(1                        :size(columnPropertiesTemp))=columnPropertiesTemp
          mergerTrees%particleColumnProperties(1+size(columnPropertiesTemp):columnNumber           )=propertyTypeNull
          call deallocateArray(columnPropertiesTemp)
       end if
    else
       call allocateArray(mergerTrees%particleColumnProperties,[columnNumber])
       mergerTrees%particleColumnProperties=propertyTypeNull
    end if
    ! Store the property type.
    mergerTrees%particleColumnProperties(columnNumber)=propertyType
    return
  end subroutine Merger_Tree_Data_Structure_Set_Particle_Property_Column

  subroutine Merger_Tree_Data_Structure_Set_Property_Column(mergerTrees,propertyType,columnNumber)
    !% Set column mapping from the input file.
    use :: Memory_Management, only : allocateArray, deallocateArray
    implicit none
    class  (mergerTreeData), intent(inout)               :: mergerTrees
    integer                , intent(in   )               :: columnNumber        , propertyType
    integer                , allocatable  , dimension(:) :: columnPropertiesTemp

    ! Ensure the storage array is large enough.
    if (allocated(mergerTrees%columnProperties)) then
       if (columnNumber > size(mergerTrees%columnProperties)) then
          call Move_Alloc(mergerTrees%columnProperties,columnPropertiesTemp)
          call allocateArray(mergerTrees%columnProperties,[columnNumber])
          mergerTrees%columnProperties(1                        :size(columnPropertiesTemp))=columnPropertiesTemp
          mergerTrees%columnProperties(1+size(columnPropertiesTemp):columnNumber           )=propertyTypeNull
          call deallocateArray(columnPropertiesTemp)
       end if
    else
       call allocateArray(mergerTrees%columnProperties,[columnNumber])
       mergerTrees%columnProperties=propertyTypeNull
    end if
    ! Store the property type.
    mergerTrees%columnProperties(columnNumber)=propertyType
    return
  end subroutine Merger_Tree_Data_Structure_Set_Property_Column

  subroutine Merger_Tree_Data_Structure_Set_Property_Integer8(mergerTrees,propertyType,property)
    !% Set a property in the merger trees.
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Memory_Management, only : allocateArray          , deallocateArray
    implicit none
    class  (mergerTreeData)              , intent(inout) :: mergerTrees
    integer                              , intent(in   ) :: propertyType
    integer(kind=kind_int8), dimension(:), intent(in   ) :: property

    ! Check the supplied arrays is of the correct size.
    if (size(property) /= mergerTrees%nodeCount) call Galacticus_Error_Report('property array size is incorrect'//{introspection:location})

    ! Assign to the relevant property.
    select case (propertyType)
    case (propertyTypeTreeIndex      )
       mergerTrees%hasForestIndex    =.true.
       if (allocated(mergerTrees%forestIndex    )) call deallocateArray(mergerTrees%forestIndex    )
       call allocateArray(mergerTrees%forestIndex    ,[size(property)])
       mergerTrees%forestIndex    =property
       call Merger_Tree_Data_Structure_Set_Tree_Indices(mergerTrees)
    case (propertyTypeNodeIndex      )
       mergerTrees%hasNodeIndex      =.true.
       if (allocated(mergerTrees%nodeIndex      )) call deallocateArray(mergerTrees%nodeIndex      )
       call allocateArray(mergerTrees%nodeIndex      ,[size(property)])
       mergerTrees%nodeIndex      =property
    case (propertyTypeHostIndex      )
       mergerTrees%hasHostIndex      =.true.
       if (allocated(mergerTrees%hostIndex      )) call deallocateArray(mergerTrees%hostIndex      )
       call allocateArray(mergerTrees%hostIndex      ,[size(property)])
       mergerTrees%hostIndex      =property
    case (propertyTypeDescendentIndex)
       mergerTrees%hasDescendentIndex=.true.
       if (allocated(mergerTrees%descendentIndex)) call deallocateArray(mergerTrees%descendentIndex)
       call allocateArray(mergerTrees%descendentIndex,[size(property)])
       mergerTrees%descendentIndex=property
    case (propertyTypeSnapshot       )
       mergerTrees%hasSnapshot       =.true.
       if (allocated(mergerTrees%snapshot       )) call deallocateArray(mergerTrees%snapshot       )
       call allocateArray(mergerTrees%snapshot       ,[size(property)])
       mergerTrees%snapshot=property
    case default
       call Galacticus_Error_Report('unrecognized integer property'//{introspection:location})
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Set_Property_Integer8

  subroutine Merger_Tree_Data_Structure_Set_Property_Double(mergerTrees,propertyType,property)
    !% Set a property in the merger trees.
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Memory_Management, only : allocateArray          , deallocateArray
    implicit none
    class           (mergerTreeData)              , intent(inout) :: mergerTrees
    integer                                       , intent(in   ) :: propertyType
    double precision                , dimension(:), intent(in   ) :: property

    ! Check the supplied arrays is of the correct size.
    if (size(property) /= mergerTrees%nodeCount) call Galacticus_Error_Report('property array size is incorrect'//{introspection:location})

    ! Assign to the relevant property.
    select case (propertyType)
    case (propertyTypeTreeWeight     )
       mergerTrees%hasForestWeight  =.true.
       if (allocated(mergerTrees%forestWeightNode)) call deallocateArray(mergerTrees%forestWeightNode)
       call allocateArray(mergerTrees%forestWeightNode ,[size(property)])
       mergerTrees%forestWeightNode =property
    case (propertyTypeRedshift       )
       mergerTrees%hasRedshift       =.true.
       if (allocated(mergerTrees%redshift        )) call deallocateArray(mergerTrees%redshift        )
       call allocateArray(mergerTrees%redshift       ,[size(property)])
       mergerTrees%redshift       =property
    case (propertyTypeNodeMass       )
       mergerTrees%hasNodeMass       =.true.
       if (allocated(mergerTrees%nodeMass        )) call deallocateArray(mergerTrees%nodeMass        )
       call allocateArray(mergerTrees%nodeMass       ,[size(property)])
       mergerTrees%nodeMass       =property
    case (propertyTypeNodeMass200Mean)
       mergerTrees%hasNodeMass200Mean=.true.
       if (allocated(mergerTrees%nodeMass200Mean )) call deallocateArray(mergerTrees%nodeMass200Mean )
       call allocateArray(mergerTrees%nodeMass200Mean,[size(property)])
       mergerTrees%nodeMass200Mean=property
    case (propertyTypeNodeMass200Crit)
       mergerTrees%hasNodeMass200Crit=.true.
       if (allocated(mergerTrees%nodeMass200Crit )) call deallocateArray(mergerTrees%nodeMass200Crit )
       call allocateArray(mergerTrees%nodeMass200Crit,[size(property)])
       mergerTrees%nodeMass200Crit=property
    case (propertyTypePositionX      )
       if (                                             &
            & allocated(mergerTrees%position)           &
            & .and..not.mergerTrees%hasPositionX        &
            & .and..not.mergerTrees%hasPositionY        &
            & .and..not.mergerTrees%hasPositionZ        &
            &) call deallocateArray(mergerTrees%position)
       mergerTrees%hasPositionX=.true.
       if (.not.allocated(mergerTrees%position)) call allocateArray(mergerTrees%position,[3,size(property)])
       mergerTrees%position(1,:)=property
    case (propertyTypePositionY      )
       if (                                             &
            & allocated(mergerTrees%position)           &
            & .and..not.mergerTrees%hasPositionX        &
            & .and..not.mergerTrees%hasPositionY        &
            & .and..not.mergerTrees%hasPositionZ        &
            &) call deallocateArray(mergerTrees%position)
       mergerTrees%hasPositionY=.true.
       if (.not.allocated(mergerTrees%position)) call allocateArray(mergerTrees%position,[3,size(property)])
       mergerTrees%position(2,:)=property
    case (propertyTypePositionZ      )
       if (                                             &
            & allocated(mergerTrees%position)           &
            & .and..not.mergerTrees%hasPositionX        &
            & .and..not.mergerTrees%hasPositionY        &
            & .and..not.mergerTrees%hasPositionZ        &
            &) call deallocateArray(mergerTrees%position)
       mergerTrees%hasPositionZ=.true.
       if (.not.allocated(mergerTrees%position)) call allocateArray(mergerTrees%position,[3,size(property)])
       mergerTrees%position(3,:)=property
    case (propertyTypeVelocityX      )
       if (                                             &
            & allocated(mergerTrees%velocity)           &
            & .and..not.mergerTrees%hasVelocityX        &
            & .and..not.mergerTrees%hasVelocityY        &
            & .and..not.mergerTrees%hasVelocityZ        &
            &) call deallocateArray(mergerTrees%velocity)
       mergerTrees%hasVelocityX=.true.
       if (.not.allocated(mergerTrees%velocity)) call allocateArray(mergerTrees%velocity,[3,size(property)])
       mergerTrees%velocity(1,:)=property
    case (propertyTypeVelocityY      )
       if (                                             &
            & allocated(mergerTrees%velocity)           &
            & .and..not.mergerTrees%hasVelocityX        &
            & .and..not.mergerTrees%hasVelocityY        &
            & .and..not.mergerTrees%hasVelocityZ        &
            &) call deallocateArray(mergerTrees%velocity)
       mergerTrees%hasVelocityY=.true.
       if (.not.allocated(mergerTrees%velocity)) call allocateArray(mergerTrees%velocity,[3,size(property)])
       mergerTrees%velocity(2,:)=property
    case (propertyTypeVelocityZ      )
       if (                                             &
            & allocated(mergerTrees%velocity)           &
            & .and..not.mergerTrees%hasVelocityX        &
            & .and..not.mergerTrees%hasVelocityY        &
            & .and..not.mergerTrees%hasVelocityZ        &
            &) call deallocateArray(mergerTrees%velocity)
       mergerTrees%hasVelocityZ=.true.
       if (.not.allocated(mergerTrees%velocity)) call allocateArray(mergerTrees%velocity,[3,size(property)])
       mergerTrees%velocity(3,:)=property
    case default
       call Galacticus_Error_Report('unrecognized double property'//{introspection:location})
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Set_Property_Double

  subroutine Merger_Tree_Data_Structure_Read_ASCII(mergerTrees,inputFile,columnHeaders,commentCharacter,separator,maximumRedshift)
    !% Read in merger tree data from an ASCII file.
    use :: File_Utilities    , only : Count_Lines_In_File
    use :: Galacticus_Display, only : Galacticus_Display_Message
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)             , operator(//)
    use :: Memory_Management , only : allocateArray             , deallocateArray
    use :: String_Handling   , only : String_Count_Words        , String_Split_Words, operator(//)
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
    mergerTrees%hasDescendentIndex                 =any(mergerTrees%columnProperties == propertyTypeDescendentIndex         )
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
         & ) call Galacticus_Error_Report("all three axes or none must be supplied for position"        //{introspection:location})
    if     (.not.((     mergerTrees%hasVelocityX       .and.     mergerTrees%hasVelocityY       .and.     mergerTrees%hasVelocityZ       ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasVelocityX       .and..not.mergerTrees%hasVelocityY       .and..not.mergerTrees%hasVelocityZ       ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("all three axes or none must be supplied for velocity"        //{introspection:location})
    if     (.not.((     mergerTrees%hasSpinX           .and.     mergerTrees%hasSpinY           .and.     mergerTrees%hasSpinZ           ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasSpinX           .and..not.mergerTrees%hasSpinY           .and..not.mergerTrees%hasSpinZ           ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("all three axes or none must be supplied for spin"            //{introspection:location})
    if     (.not.((     mergerTrees%hasAngularMomentumX        .and.     mergerTrees%hasAngularMomentumY        .and.     mergerTrees%hasAngularMomentumZ        ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasAngularMomentumX        .and..not.mergerTrees%hasAngularMomentumY        .and..not.mergerTrees%hasAngularMomentumZ        ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("all three axes or none must be supplied for angular momentum"         //{introspection:location})
    if     (.not.((     mergerTrees%hasSpecificAngularMomentumX.and.     mergerTrees%hasSpecificAngularMomentumY.and.     mergerTrees%hasSpecificAngularMomentumZ) &
         &        .or.                                                                                                                                             &
         &        (.not.mergerTrees%hasSpecificAngularMomentumX.and..not.mergerTrees%hasSpecificAngularMomentumY.and..not.mergerTrees%hasSpecificAngularMomentumZ) &
         &        )                                                                                                                                                &
         & ) call Galacticus_Error_Report("all three axes or none must be supplied for specific angular momentum"//{introspection:location})
    if (mergerTrees%hasSpinX                   .and.mergerTrees%hasSpinMagnitude                   ) &
         & call Galacticus_Error_Report("can not specify both 3-D and scalar spin"                     //{introspection:location})
    if (mergerTrees%hasAngularMomentumX        .and.mergerTrees%hasAngularMomentumMagnitude        ) &
         & call Galacticus_Error_Report("can not specify both 3-D and scalar angular momentum"         //{introspection:location})
    if (mergerTrees%hasSpecificAngularMomentumX.and.mergerTrees%hasSpecificAngularMomentumMagnitude) &
         & call Galacticus_Error_Report("can not specify both 3-D and scalar specific angular momentum"//{introspection:location})

    ! Validate that either redshift or scale factor is given.
    if (.not.(mergerTrees%hasRedshift .or. mergerTrees%hasScaleFactor)) &
         & call Galacticus_Error_Report("either redshift or scale factor has to be given"//{introspection:location})

    ! Deallocate internal arrays.
    if (allocated(mergerTrees%forestIndex           )) call deallocateArray(mergerTrees%forestIndex           )
    if (allocated(mergerTrees%forestWeightNode      )) call deallocateArray(mergerTrees%forestWeightNode      )
    if (allocated(mergerTrees%nodeIndex             )) call deallocateArray(mergerTrees%nodeIndex             )
    if (allocated(mergerTrees%mostBoundParticleIndex)) call deallocateArray(mergerTrees%mostBoundParticleIndex)
    if (allocated(mergerTrees%snapshot              )) call deallocateArray(mergerTrees%snapshot              )
    if (allocated(mergerTrees%descendentIndex       )) call deallocateArray(mergerTrees%descendentIndex       )
    if (allocated(mergerTrees%hostIndex             )) call deallocateArray(mergerTrees%hostIndex             )
    if (allocated(mergerTrees%redshift              )) call deallocateArray(mergerTrees%redshift              )
    if (allocated(mergerTrees%scaleFactor           )) call deallocateArray(mergerTrees%scaleFactor           )
    if (allocated(mergerTrees%nodeMass              )) call deallocateArray(mergerTrees%nodeMass              )
    if (allocated(mergerTrees%nodeMass200Mean       )) call deallocateArray(mergerTrees%nodeMass200Mean       )
    if (allocated(mergerTrees%nodeMass200Crit       )) call deallocateArray(mergerTrees%nodeMass200Crit       )
    if (allocated(mergerTrees%particleCount         )) call deallocateArray(mergerTrees%particleCount         )
    if (allocated(mergerTrees%position              )) call deallocateArray(mergerTrees%position              )
    if (allocated(mergerTrees%velocity              )) call deallocateArray(mergerTrees%velocity              )
    if (allocated(mergerTrees%spin                  )) call deallocateArray(mergerTrees%spin                  )
    if (allocated(mergerTrees%halfMassRadius        )) call deallocateArray(mergerTrees%halfMassRadius        )
    if (allocated(mergerTrees%scaleRadius           )) call deallocateArray(mergerTrees%scaleRadius           )
    if (allocated(mergerTrees%velocityMaximum       )) call deallocateArray(mergerTrees%velocityMaximum       )
    if (allocated(mergerTrees%velocityDispersion    )) call deallocateArray(mergerTrees%velocityDispersion    )

    ! Allocate internal arrays to correct size as needed.
    if (mergerTrees%hasForestIndex                     ) call allocateArray(mergerTrees%forestIndex                     ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasForestWeight                    ) call allocateArray(mergerTrees%forestWeightNode                ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasNodeIndex                       ) call allocateArray(mergerTrees%nodeIndex                       ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasMostBoundParticleIndex          ) call allocateArray(mergerTrees%mostBoundParticleIndex          ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasSnapshot                        ) call allocateArray(mergerTrees%snapshot                        ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasDescendentIndex                 ) call allocateArray(mergerTrees%descendentIndex                 ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasHostIndex                       ) call allocateArray(mergerTrees%hostIndex                       ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasRedshift                        ) call allocateArray(mergerTrees%redshift                        ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasScaleFactor                     ) call allocateArray(mergerTrees%scaleFactor                     ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasNodeMass                        ) call allocateArray(mergerTrees%nodeMass                        ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasNodeMass200Mean                 ) call allocateArray(mergerTrees%nodeMass200Mean                 ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasNodeMass200Crit                 ) call allocateArray(mergerTrees%nodeMass200Crit                 ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasParticleCount                   ) call allocateArray(mergerTrees%particleCount                   ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasPositionX                       ) call allocateArray(mergerTrees%position                        ,[3,mergerTrees%nodeCount])
    if (mergerTrees%hasVelocityX                       ) call allocateArray(mergerTrees%velocity                        ,[3,mergerTrees%nodeCount])
    if (mergerTrees%hasSpinX                           ) call allocateArray(mergerTrees%spin                            ,[3,mergerTrees%nodeCount])
    if (mergerTrees%hasAngularMomentumX                ) call allocateArray(mergerTrees%angularMomentum                 ,[3,mergerTrees%nodeCount])
    if (mergerTrees%hasSpecificAngularMomentumX        ) call allocateArray(mergerTrees%specificAngularMomentum         ,[3,mergerTrees%nodeCount])
    if (mergerTrees%hasSpinMagnitude                   ) call allocateArray(mergerTrees%spinMagnitude                   ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasAngularMomentumMagnitude        ) call allocateArray(mergerTrees%angularMomentumMagnitude        ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasSpecificAngularMomentumMagnitude) call allocateArray(mergerTrees%specificAngularMomentumMagnitude,[  mergerTrees%nodeCount])
    if (mergerTrees%hasHalfMassRadius                  ) call allocateArray(mergerTrees%halfMassRadius                  ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasScaleRadius                     ) call allocateArray(mergerTrees%scaleRadius                     ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasVelocityMaximum                 ) call allocateArray(mergerTrees%velocityMaximum                 ,[  mergerTrees%nodeCount])
    if (mergerTrees%hasVelocityDispersion              ) call allocateArray(mergerTrees%velocityDispersion              ,[  mergerTrees%nodeCount])

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
             call allocateArray(inputColumns,[columnsCount])
             gotFirstDataLine=.true.
          end if
          ! Count nodes.
          iNode=iNode+1
          call String_Split_Words(inputColumns,inputLine,separator)
          do iColumn=1,min(columnsCount,size(mergerTrees%columnProperties))
             select case (mergerTrees%columnProperties(iColumn))
             case (propertyTypeNull                  )
                ! Ignore this column.
             case (propertyTypeTreeIndex             )
                ! Column is a tree index.
                read (inputColumns(iColumn),*) mergerTrees%forestIndex(iNode)
                if (iNode > 1) then
                   if (mergerTrees%forestIndex(iNode) < mergerTrees%forestIndex(iNode-1)) &
                        & call Galacticus_Error_Report('tree indices must be in ascending order'//{introspection:location})
                   if (mergerTrees%forestIndex(iNode) /= mergerTrees%forestIndex(iNode-1)) mergerTrees%forestCount=mergerTrees%forestCount+1
                else
                   mergerTrees%forestCount=1
                end if
             case (propertyTypeTreeWeight            )
                ! Column is a tree weight.
                read (inputColumns(iColumn),*) mergerTrees%forestWeightNode        (  iNode)
             case (propertyTypeNodeIndex             )
                ! Column is a node index.
                read (inputColumns(iColumn),*) mergerTrees%nodeIndex               (  iNode)
             case (propertyTypeDescendentIndex       )
                ! Column is a descendent node index.
                read (inputColumns(iColumn),*) mergerTrees%descendentIndex         (  iNode)
             case (propertyTypeHostIndex             )
                ! Column is a host index.
                read (inputColumns(iColumn),*) mergerTrees%hostIndex               (  iNode)
             case (propertyTypeRedshift              )
                ! Column is redshift.
                read (inputColumns(iColumn),*) mergerTrees%redshift                (  iNode)
              case (propertyTypeScaleFactor          )
                ! Column is scale factor.
                read (inputColumns(iColumn),*) mergerTrees%scaleFactor             (  iNode)
             case (propertyTypeNodeMass              )
                ! Column is mass.
                read (inputColumns(iColumn),*) mergerTrees%nodeMass                (  iNode)
             case (propertyTypeNodeMass200Mean       )
                ! Column is mass.
                read (inputColumns(iColumn),*) mergerTrees%nodeMass200Mean         (  iNode)
             case (propertyTypeNodeMass200Crit       )
                ! Column is mass.
                read (inputColumns(iColumn),*) mergerTrees%nodeMass200Crit         (  iNode)
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
                read (inputColumns(iColumn),*) mergerTrees%angularMomentumMagnitude        (  iNode)
             case (propertyTypeSpecificAngularMomentumX)
                ! Column is x specific angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%specificAngularMomentum         (1,iNode)
             case (propertyTypeSpecificAngularMomentumY)
                ! Column is y specific angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%specificAngularMomentum         (2,iNode)
             case (propertyTypeSpecificAngularMomentumZ)
                ! Column is z specific angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%specificAngularMomentum         (3,iNode)
             case (propertyTypeSpecificAngularMomentum )
                ! Column is scalar specific angular momentum.
                read (inputColumns(iColumn),*) mergerTrees%specificAngularMomentumMagnitude(  iNode)
             case (propertyTypeHalfMassRadius        )
                ! Column is half mass radius.
                read (inputColumns(iColumn),*) mergerTrees%halfMassRadius          (  iNode)
             case (propertyTypeScaleRadius           )
                ! Column is scale radius.
                read (inputColumns(iColumn),*) mergerTrees%scaleRadius             (  iNode)
             case (propertyTypeMostBoundParticleIndex  )
                ! Column is a most bound particle index.
                read (inputColumns(iColumn),*) mergerTrees%mostBoundParticleIndex  (  iNode)
             case (propertyTypeSnapshot              )
                ! Column is a snapshot index.
                read (inputColumns(iColumn),*) mergerTrees%snapshot                (  iNode)
             case (propertyTypeVelocityMaximum         )
                ! Column is a maximum velocity.
                read (inputColumns(iColumn),*) mergerTrees%velocityMaximum                 (  iNode)
             case (propertyTypeVelocityDispersion      )
                ! Column is a velocity dispersion.
                read (inputColumns(iColumn),*) mergerTrees%velocityDispersion              (  iNode)
             case default
                call Galacticus_Error_Report('unknown column type'//{introspection:location})
             end select
          end do
       end if
    end do
    close(fileUnit)

    ! Report number of forests found.
    message='Found '
    message=message//mergerTrees%forestCount//' forests'
    call Galacticus_Display_Message(message)

    ! Deallocate workspace.
    if (allocated(inputColumns)) call deallocateArray(inputColumns)

    ! If we have the particle mass, set the masses of any subhalos (which have zero mass by default) based on particle count.
    call Merger_Tree_Data_Set_Subhalo_Masses(mergerTrees)

    ! Convert specific angular momenta as needed.
    if (mergerTrees%hasSpecificAngularMomentumMagnitude.and..not.mergerTrees%hasAngularMomentumMagnitude) then
       call allocateArray  (mergerTrees%angularMomentumMagnitude        ,[  mergerTrees%nodeCount])
       mergerTrees%angularMomentumMagnitude=mergerTrees%specificAngularMomentumMagnitude*mergerTrees%nodeMass
       call deallocateArray(mergerTrees%specificAngularMomentumMagnitude                          )
       mergerTrees%hasSpecificAngularMomentumMagnitude=.false.
       mergerTrees%hasAngularMomentumMagnitude        =.true.
    end if
    if (mergerTrees%hasSpecificAngularMomentumX.and..not.mergerTrees%hasAngularMomentumX) then
       call allocateArray  (mergerTrees%angularMomentum        ,[3,mergerTrees%nodeCount])
       forall(i=1:3)
          mergerTrees%angularMomentum(i,:)=mergerTrees%specificAngularMomentum(i,:)*mergerTrees%nodeMass
       end forall
       call deallocateArray(mergerTrees%specificAngularMomentum                          )
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
        call allocateArray(mergerTrees%redshift,[ mergerTrees%nodeCount])
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
        call deallocateArray(mergerTrees%scaleFactor)
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

  subroutine Merger_Tree_Data_Structure_Convert_Property_Units(mergerTrees, propertyType, conversionFactor)
    !% Convert the property with inconsistent units.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(mergerTreeData), intent(inout) :: mergerTrees
    integer              , intent(in   ) :: propertyType
    double precision     , intent(in   ) :: conversionFactor
    integer                              :: j

    select case (propertyType)
    case (propertyTypeNodeMass          )
       ! Property is mass.
       mergerTrees%nodeMass                    =mergerTrees%nodeMass                     *conversionFactor
    case (propertyTypeNodeMass200Mean   )
       ! Property is mass.
       mergerTrees%nodeMass200Mean             =mergerTrees%nodeMass200Mean              *conversionFactor
    case (propertyTypeNodeMass200Crit   )
       ! Property is mass.
       mergerTrees%nodeMass200Crit             =mergerTrees%nodeMass200Crit              *conversionFactor
    case (propertyTypePositionX         )
       ! Property is position.
       forall(j=1:3)
          mergerTrees%position            (j,:)=mergerTrees%position                (j,:)*conversionFactor
       end forall
    case (propertyTypeVelocityX         )
       ! Property is velocity.
       forall(j=1:3)
          mergerTrees%velocity            (j,:)=mergerTrees%velocity                (j,:)*conversionFactor
       end forall
    case (propertyTypeAngularMomentumX  )
       ! Property is angular momentum vector.
       forall(j=1:3)
          mergerTrees%angularMomentum     (j,:)=mergerTrees%angularMomentum         (j,:)*conversionFactor
       end forall
    case (propertyTypeAngularMomentum   )
       ! Property is scalar angular momentum.
       mergerTrees%angularMomentumMagnitude    =mergerTrees%angularMomentumMagnitude     *conversionFactor
    case (propertyTypeHalfMassRadius    )
       ! Property is half mass radius.
       mergerTrees%halfMassRadius              =mergerTrees%halfMassRadius               *conversionFactor
    case (propertyTypeScaleRadius       )
       ! Property is scale radius.
       mergerTrees%scaleRadius                 =mergerTrees%scaleRadius                  *conversionFactor
    case (propertyTypeVelocityMaximum   )
       ! Property is maximum velocity.
       mergerTrees%velocityMaximum             =mergerTrees%velocityMaximum              *conversionFactor
    case (propertyTypeVelocityDispersion)
       ! Property is velocity dispersion.
       mergerTrees%velocityDispersion          =mergerTrees%velocityDispersion           *conversionFactor
    case default
       ! Property has no units.
       call Galacticus_Error_Report('property has no units to convert.'//{introspection:location})
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Convert_Property_Units

  subroutine Merger_Tree_Data_Structure_Set_Tree_Indices(mergerTrees)
    !% Set the merger tree index arrays.
    use :: Memory_Management, only : allocateArray, deallocateArray
    implicit none
    class  (mergerTreeData), intent(inout) :: mergerTrees
    integer                                :: iNode      , iTree

    ! Allocate arrays for tree start and stop indices and reference ID.
    if (allocated(mergerTrees%treeBeginsAt )) call deallocateArray(mergerTrees%treeBeginsAt )
    if (allocated(mergerTrees%treeNodeCount)) call deallocateArray(mergerTrees%treeNodeCount)
    if (allocated(mergerTrees%forestID     )) call deallocateArray(mergerTrees%forestID     )
    if (allocated(mergerTrees%treeWeight   )) call deallocateArray(mergerTrees%treeWeight   )
    call allocateArray(mergerTrees%treeBeginsAt ,[mergerTrees%forestCount])
    call allocateArray(mergerTrees%treeNodeCount,[mergerTrees%forestCount])
    call allocateArray(mergerTrees%forestID     ,[mergerTrees%forestCount])
    call allocateArray(mergerTrees%treeWeight   ,[mergerTrees%forestCount])

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
    !% Read in particle data from an ASCII file.
    use :: File_Utilities   , only : Count_Lines_In_File
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Memory_Management, only : allocateArray          , deallocateArray
    use :: String_Handling  , only : String_Count_Words     , String_Split_Words
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
         & ) call Galacticus_Error_Report("all three axes or none must be supplied for particle position"//{introspection:location})
    if     (.not.((     mergerTrees%hasParticleVelocityX.and.     mergerTrees%hasParticleVelocityY.and.     mergerTrees%hasParticleVelocityZ) &
         &        .or.                                                                                                &
         &        (.not.mergerTrees%hasParticleVelocityX.and..not.mergerTrees%hasParticleVelocityY.and..not.mergerTrees%hasParticleVelocityZ) &
         &        )                                                                                                   &
         & ) call Galacticus_Error_Report("all three axes or none must be supplied for particle velocity"//{introspection:location})

    ! Ensure we have a redshift.
    if (.not.mergerTrees%hasParticleRedshift) call Galacticus_Error_Report("particle redshift must be supplied"//{introspection:location})

    ! Deallocate internal arrays.
    if (allocated(mergerTrees%particleIndex   )) call deallocateArray(mergerTrees%particleIndex   )
    if (allocated(mergerTrees%particleRedshift)) call deallocateArray(mergerTrees%particleRedshift)
    if (allocated(mergerTrees%particlePosition)) call deallocateArray(mergerTrees%particlePosition)
    if (allocated(mergerTrees%particleVelocity)) call deallocateArray(mergerTrees%particleVelocity)
    if (allocated(mergerTrees%particleSnapshot)) call deallocateArray(mergerTrees%particleSnapshot)

    ! Allocate internal arrays to correct size as needed.
    if (mergerTrees%hasParticleIndex    ) call allocateArray(mergerTrees%particleIndex   ,[  mergerTrees%particlesCount])
    if (mergerTrees%hasParticleRedshift ) call allocateArray(mergerTrees%particleRedshift,[  mergerTrees%particlesCount])
    if (mergerTrees%hasParticlePositionX) call allocateArray(mergerTrees%particlePosition,[3,mergerTrees%particlesCount])
    if (mergerTrees%hasParticleVelocityX) call allocateArray(mergerTrees%particleVelocity,[3,mergerTrees%particlesCount])
    if (mergerTrees%hasParticleSnapshot ) call allocateArray(mergerTrees%particleSnapshot,[  mergerTrees%particlesCount])

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
             call allocateArray(inputColumns,[columnsCount])
             gotFirstDataLine=.true.
          end if
          ! Count nodes.
          iNode=iNode+1
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
                call Galacticus_Error_Report('unknown column type'//{introspection:location})
             end select
          end do
       end if
    end do
    close(fileUnit)

    ! Deallocate workspace.
    if (allocated(inputColumns)) call deallocateArray(inputColumns)

    return
  end subroutine Merger_Tree_Data_Structure_Read_Particles_ASCII

  subroutine Merger_Tree_Data_Structure_Export(mergerTrees,outputFileName,outputFormat,hdfChunkSize,hdfCompressionLevel,append)
    !% Output a set of merger trees to an HDF5 file.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: HDF5            , only : hsize_t
    implicit none
    integer  (kind=hsize_t  ), intent(in   )           :: hdfChunkSize
    integer                  , intent(in   )           :: hdfCompressionLevel, outputFormat
    class    (mergerTreeData), intent(inout)           :: mergerTrees
    character(len=*         ), intent(in   )           :: outputFileName
    logical                  , intent(in   ), optional :: append

    ! Validate the merger tree.
    call Merger_Tree_Data_Validate_Trees            (mergerTrees)

    ! If we have most-bound particle indices and particle data has been read, construct arrays giving position of particle data for each node.
    call Merger_Tree_Data_Construct_Particle_Indices(mergerTrees)

    select case (outputFormat)
    case (mergerTreeFormatGalacticus)
       call Merger_Tree_Data_Structure_Export_Galacticus(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    case (mergerTreeFormatIrate     )
       call Merger_Tree_Data_Structure_Export_IRATE     (mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    case default
       call Galacticus_Error_Report('output format is not recognized'//{introspection:location})
    end select
    return
  end subroutine Merger_Tree_Data_Structure_Export

  subroutine Merger_Tree_Data_Structure_Export_Galacticus(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    !% Output a set of merger trees to a Galacticus-format HDF5 file.
    use :: File_Utilities    , only : File_Exists
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: HDF5              , only : HSIZE_T                , hsize_t
    use :: IO_HDF5           , only : hdf5Access             , hdf5Object
    use :: ISO_Varying_String, only : assignment(=)          , char
    use :: Memory_Management , only : deallocateArray
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
    call hdf5Access%set()
    call outputFile%openFile(outputFileName,overWrite=.not.appendActual,objectsOverwritable=.true.,chunkSize=hdfChunkSize,compressionLevel=hdfCompressionLevel)


    ! Write a format version attribute.
    if (.not.fileExists) call outputFile%writeAttribute(2,"formatVersion")

    ! Create a group for the datasets.
    forestHalos=outputFile%openGroup("forestHalos","Stores all data for merger trees.")

    ! Write the data.
    if (mergerTrees%hasNodeIndex               ) call forestHalos%writeDataset(mergerTrees%nodeIndex               ,"nodeIndex"          ,"The index of each node."                            ,appendTo=appendActual                  )
    if (mergerTrees%hasDescendentIndex         ) call forestHalos%writeDataset(mergerTrees%descendentIndex         ,"descendentIndex"    ,"The index of each descendent node."                 ,appendTo=appendActual                  )
    if (mergerTrees%hasHostIndex               ) call forestHalos%writeDataset(mergerTrees%hostIndex               ,"hostIndex"          ,"The index of each host node."                       ,appendTo=appendActual                  )
    if (mergerTrees%hasNodeMass                ) call forestHalos%writeDataset(mergerTrees%nodeMass                ,"nodeMass"           ,"The mass of each node."                             ,appendTo=appendActual                  )
    if (mergerTrees%hasRedshift                ) call forestHalos%writeDataset(mergerTrees%redshift                ,"redshift"           ,"The redshift of each node."                         ,appendTo=appendActual                  )
    if (mergerTrees%hasNodeMass200Mean         ) call forestHalos%writeDataset(mergerTrees%nodeMass200Mean         ,"nodeMass200Mean"    ,"The M200 mass of each node (200 * mean density)."   ,appendTo=appendActual                  )
    if (mergerTrees%hasNodeMass200Crit         ) call forestHalos%writeDataset(mergerTrees%nodeMass200Crit         ,"nodeMass200Crit"    ,"The M200 mass of each node (200 * crit density(."   ,appendTo=appendActual                  )
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
             ! Skip null and tree indxe cases.
             if     (                                    &
                  &   iProperty == propertyTypeNull      &
                  &  .or.                                &
                  &   iProperty == propertyTypeTreeIndex &
                  & ) cycle
             ! Skip cases where we have the corresponding 3-D dataset.
             if (iProperty == propertyTypeSpin            .and. .not.mergerTrees%hasSpinMagnitude           ) cycle
             if (iProperty == propertyTypeAngularMomentum .and. .not.mergerTrees%hasAngularMomentumMagnitude) cycle
             if (forestHalos%hasDataset(char(enumerationPropertyTypeDecode(iProperty)))) then
                treeDataset=forestHalos%openDataset(char(enumerationPropertyTypeDecode(iProperty)))
                call treeGroup%createReference1D(treeDataset,char(enumerationPropertyTypeDecode(iProperty)),hyperslabStart,hyperslabCount)
                call treeDataset%close()
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
                call treeDataset%close()
             end if
          end do

          call treeGroup%close()
       end do

       ! Close the merger trees group.
       call forestsGroup%close()

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

       ! Close the particles group.
       call particlesGroup%close()
    end if

    ! Create datasets giving positions of merger trees within the node arrays.
    forestIndexGroup=outputFile%openGroup("forestIndex","Locations of forests within the halo data arrays.")
    if (fileExists) then
       call forestIndexGroup%readDataset("firstNode"    ,firstNode    )
       call forestIndexGroup%readDataset("numberOfNodes",numberOfNodes)
       mergerTrees%treeBeginsAt=mergerTrees%treeBeginsAt+firstNode(size(firstNode))+numberOfNodes(size(numberOfNodes))
       call deallocateArray(firstNode    )
       call deallocateArray(numberOfNodes)
    end if
    call        forestIndexGroup%writeDataset(mergerTrees%treeBeginsAt ,"firstNode"    ,"Position of the first node in each forest in the halo data arrays.",appendTo=appendActual)
    call        forestIndexGroup%writeDataset(mergerTrees%treeNodeCount,"numberOfNodes","Number of nodes in each forest."                                   ,appendTo=appendActual)
    call        forestIndexGroup%writeDataset(mergerTrees%forestID     ,"forestIndex"  ,"Unique index of forest."                                           ,appendTo=appendActual)
    if (mergerTrees%hasForestWeight.or..not.mergerTrees%hasBoxSize)                                                                                                               &
         & call forestIndexGroup%writeDataset(mergerTrees%treeWeight   ,"forestWeight" ,"Weight of forest."                                                 ,appendTo=appendActual)
    call forestIndexGroup%close()

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
       if (mergerTrees%unitsSet(unitsMass    )) call Store_Unit_Attributes_Galacticus(unitsMass    ,"mass"    ,mergerTrees,unitsGroup)
       if (mergerTrees%unitsSet(unitsLength  )) call Store_Unit_Attributes_Galacticus(unitsLength  ,"length"  ,mergerTrees,unitsGroup)
       if (mergerTrees%unitsSet(unitsTime    )) call Store_Unit_Attributes_Galacticus(unitsTime    ,"time"    ,mergerTrees,unitsGroup)
       if (mergerTrees%unitsSet(unitsVelocity)) call Store_Unit_Attributes_Galacticus(unitsVelocity,"velocity",mergerTrees,unitsGroup)
       call unitsGroup%close()

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
          select case (mergerTrees%metaData(iAttribute)%metadataType)
          case (metaDataTypeGeneric    )
             attributeGroup => genericGroup
          case (metaDataTypeCosmology  )
             attributeGroup => cosmologyGroup
          case (metaDataTypeSimulation )
             attributeGroup => simulationGroup
          case (metaDataTypeGroupFinder)
             attributeGroup => groupFinderGroup
          case (metaDataTypeTreeBuilder)
             attributeGroup => treeBuilderGroup
          case (metaDataTypeProvenance )
             attributeGroup => provenanceGroup
          case default
             attributeGroup => null()
             call Galacticus_Error_Report('unknown meta-data group'//{introspection:location})
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
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeGeneric    )) call genericGroup    %close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeCosmology  )) call cosmologyGroup  %close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeSimulation )) call simulationGroup %close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeGroupFinder)) call groupFinderGroup%close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeTreeBuilder)) call treeBuilderGroup%close()
       if (any(mergerTrees%metaData(1:mergerTrees%metaDataCount)%metadataType == metaDataTypeProvenance )) call provenanceGroup %close()

    end if

    ! Close the group for datasets.
    call forestHalos%close()

    ! Add a flag to indicate successfully completed writing merger tree information.
    if (fileExists) then
       call outputFile%readAttribute('fileCompleteFlag',completeCount)
       completeCount=completeCount+1
    else
       completeCount=             +1
    end if
    call outputFile%writeAttribute(completeCount,"fileCompleteFlag")

    ! Close the output file.
    call outputFile%close()
    call hdf5Access%unset()

    return
  end subroutine Merger_Tree_Data_Structure_Export_Galacticus

  subroutine Merger_Tree_Data_Structure_Export_IRATE(mergerTrees,outputFileName,hdfChunkSize,hdfCompressionLevel,append)
    !% Output a set of merger trees to an IRATE-format HDF5 file.
    use :: Array_Utilities   , only : Array_Index            , Array_Which
    use :: File_Utilities    , only : File_Exists
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: HDF5              , only : hsize_t
    use :: IO_HDF5           , only : hdf5Access             , hdf5Object
    use :: ISO_Varying_String, only : assignment(=)          , char
    use :: Memory_Management , only : allocateArray          , deallocateArray
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
         &                                                                               thisDataset
    integer                         , allocatable, dimension(:)                ::        nodeSnapshotIndices           , thisSnapshotIndices
    integer         (c_size_t      ), allocatable, dimension(:)                ::        descendentSnapshot
    double precision                , allocatable, dimension(:)                ::        particleMass
    integer                                                                    ::        iAttribute                    , nodesOnSnapshotCount, &
         &                                                                               particlesOnSnapshotCount
    integer         (c_size_t      )                                           ::        iDescendent                   , iNode               , &
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
    if (.not.mergerTrees%hasSnapshot       ) call Galacticus_Error_Report('snapshot indices are required for this format'  //{introspection:location})
    if (.not.mergerTrees%hasPositionX      ) call Galacticus_Error_Report('halo positions are required for this format'    //{introspection:location})
    if (.not.mergerTrees%hasNodeIndex      ) call Galacticus_Error_Report('halo indices are required for this format'      //{introspection:location})
    if (.not.mergerTrees%hasDescendentIndex) call Galacticus_Error_Report('descendent indices are required for this format'//{introspection:location})

    ! Open the output file.
    call hdf5Access%set()
    call outputFile%openFile(outputFileName,overWrite=.not.appendActual,chunkSize=hdfChunkSize,compressionLevel=hdfCompressionLevel)

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
       call allocateArray(thisSnapshotIndices,[nodesOnSnapshotCount])
       call Array_Which(mergerTrees%snapshot == iSnapshot,thisSnapshotIndices)

       ! Write redshift attribute.
       if (.not.fileExists) call snapshotGroup%writeAttribute(mergerTrees%redshift(thisSnapshotIndices(1)),"Redshift")

       ! Write the data.
       if (mergerTrees%hasNodeIndex               ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeIndex               ,thisSnapshotIndices),"Index"          ,"The index of each halo."                                                     ,appendTo=appendActual                  )
       end if
       if (mergerTrees%hasNodeMass                ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeMass                ,thisSnapshotIndices),"Mass"           ,"The mass of each halo."                          ,datasetReturned=thisDataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasNodeMass200Mean         ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeMass200Mean         ,thisSnapshotIndices),"Mass200Mean"    ,"The M200 mass of each halo (200 * mean density).",datasetReturned=thisDataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasNodeMass200Crit         ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%nodeMass200Crit         ,thisSnapshotIndices),"Mass200Crit"    ,"The M200 mass of each halo (200 * mean density).",datasetReturned=thisDataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasPositionX               ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%position    ,thisSnapshotIndices,indexOn=2),"Center"         ,"The position of each halo center."                 ,datasetReturned=thisDataset,appendTo=appendActual,appendDimension=2)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([          unitsLength              ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasVelocityX               ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%velocity   ,thisSnapshotIndices,indexOn=2),"Velocity"       ,"The velocity of each halo."                         ,datasetReturned=thisDataset,appendTo=appendActual,appendDimension=2)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsVelocity                      ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasSpinX                   ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%spin                    ,thisSnapshotIndices),"Spin"           ,"The spin of each halo."                                                      ,appendTo=appendActual                  )
       end if
       if (mergerTrees%hasAngularMomentumX        ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%angularMomentum         ,thisSnapshotIndices),"AngularMomentum","The angular momentum spin of each halo."         ,datasetReturned=thisDataset,appendTo=appendActual,appendDimension=2)
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass,unitsLength,unitsVelocity],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasSpinMagnitude           ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%spinMagnitude           ,thisSnapshotIndices),"Spin"           ,"The spin of each halo."                                                      ,appendTo=appendActual                  )
       end if
       if (mergerTrees%hasAngularMomentumMagnitude) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%angularMomentumMagnitude,thisSnapshotIndices),"AngularMomentum","The angular momentum spin of each halo."         ,datasetReturned=thisDataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass,unitsLength,unitsVelocity],mergerTrees,thisDataset)
          call thisDataset%close()
       end if
       if (mergerTrees%hasHalfMassRadius          ) then
          call haloTrees%writeDataset(Array_Index(mergerTrees%halfMassRadius          ,thisSnapshotIndices),"HalfMassRadius" ,"The half mass radius of each halo."              ,datasetReturned=thisDataset,appendTo=appendActual                  )
          if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsMass                          ],mergerTrees,thisDataset)
          call thisDataset%close()
       end if

       ! Destroy snapshot indices.
       call deallocateArray(thisSnapshotIndices)

       ! Close the group for halo catalogs
       call haloTrees%close()

       ! Close the snapshot group.
       call snapshotGroup%close()

    end do

    ! Create a group for merger trees.
    mergerTreesGroup=outputFile%openGroup("MergerTrees","Stores all data for merger trees.")

    ! Specify the name of the halo catalog group.
    if (.not.fileExists) call mergerTreesGroup%writeAttribute("HaloCatalog","HaloCatalogName")

    ! Build snapshot numbers for descendents.
    call allocateArray(descendentSnapshot,shape(mergerTrees%nodeIndex))
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
    call                               mergerTreesGroup%writeDataset(mergerTrees%snapshot          ,"HaloSnapshot"      ,"The snapshot of each halo."           ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(mergerTrees%nodeIndex         ,"HaloID"            ,"The index of each halo."              ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(mergerTrees%descendentIndex   ,"DescendentID"      ,"The index of each descendent halo."   ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(            descendentSnapshot,"DescendentSnapshot","The snapshot of each descendent halo.",appendTo=appendActual)
    if (mergerTrees%hasHostIndex) call mergerTreesGroup%writeDataset(mergerTrees%hostIndex         ,"HostID"            ,"The index of each host halo."         ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(mergerTrees%treeNodeCount     ,"HalosPerTree"      ,"Number of halos in each tree."        ,appendTo=appendActual)
    call                               mergerTreesGroup%writeDataset(mergerTrees%forestID          ,"TreeID"            ,"Unique index of tree."                ,appendTo=appendActual)
    call deallocateArray(descendentSnapshot)
    call mergerTreesGroup%close()

    if (mergerTrees%hasMostBoundParticleIndex) then
       ! Find the highest and lowest snapshot numbers in the particles.
       if (.not.mergerTrees%hasParticleSnapshot) call Galacticus_Error_Report('particle snapshot numbers must be available for IRATE format export'//{introspection:location})
       snapshotMinimum=minval(mergerTrees%particleSnapshot)
       snapshotMaximum=maxval(mergerTrees%particleSnapshot)

       ! Loop over snapshots to output.
       do iSnapshot=snapshotMinimum,snapshotMaximum

          ! Find those particles which exist at this snapshot.
          particlesOnSnapshotCount=count(mergerTrees%particleSnapshot == iSnapshot)
          call allocateArray(thisSnapshotIndices,[particlesOnSnapshotCount])
          call Array_Which(mergerTrees%particleSnapshot == iSnapshot,thisSnapshotIndices)

          ! Find those nodes which exist at this snapshot.
          nodesOnSnapshotCount=count(mergerTrees%snapshot == iSnapshot)
          call allocateArray(nodeSnapshotIndices,[nodesOnSnapshotCount])
          call Array_Which(mergerTrees%snapshot == iSnapshot,nodeSnapshotIndices)

          ! Create a snapshot group.
          write (snapshotGroupName,'("Snapshot",i5.5)') iSnapshot
          snapshotGroup=outputFile%openGroup(trim(snapshotGroupName),"Stores all data for a snapshot.")

          ! Create a group for halo catalogs.
          haloTrees=snapshotGroup%openGroup("HaloCatalog","Stores all data for halo catalogs.")

          ! Write the particle indices.
          call haloTrees%writeDataset(Array_Index(mergerTrees%mostBoundParticleIndex,nodeSnapshotIndices),"MostBoundParticleID","The index of each particle.",appendTo=appendActual)
          call haloTrees%close()

          ! Create a group for particles.
          particlesGroup=snapshotGroup%openGroup("ParticleData","Stores all data for particles.")

          ! Make a group for dark matter particles.
          darkParticlesGroup=particlesGroup%openGroup("Dark","Stores all data for dark matter particles.")

          ! Write redshift attribute.
          if (.not.snapshotGroup%hasAttribute("Redshift")) call snapshotGroup%writeAttribute(mergerTrees%particleRedshift(thisSnapshotIndices(1)),"Redshift")

          ! Write the data.
          call allocateArray(particleMass,[particlesOnSnapshotCount])
          particleMass=mergerTrees%particleMass
          call darkParticlesGroup%writeDataset(particleMass,"Mass","The mass of each particle.",datasetReturned=thisDataset,appendTo=appendActual)
          if (.not.fileExists) call Store_Unit_Attributes_IRATE([unitsMass],mergerTrees,thisDataset)
          call thisDataset%close()
          call deallocateArray(particleMass)
          if (mergerTrees%hasParticleIndex    ) then
             call darkParticlesGroup%writeDataset(Array_Index(mergerTrees%particleIndex   ,thisSnapshotIndices),"ID"      ,"The index of each particle."                               ,appendTo=appendActual)
          end if
          if (mergerTrees%hasParticlePositionX) then
             call darkParticlesGroup%writeDataset(Array_Index(mergerTrees%particlePosition,thisSnapshotIndices),"Position","The position of each particle.",datasetReturned=thisDataset,appendTo=appendActual)
             if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsLength  ],mergerTrees,thisDataset)
             call thisDataset%close()
          end if
          if (mergerTrees%hasParticleVelocityX) then
             call darkParticlesGroup%writeDataset(Array_Index(mergerTrees%particleVelocity,thisSnapshotIndices),"Velocity","The velocity of each particle.",datasetReturned=thisDataset,appendTo=appendActual)
             if (.not.appendActual) call Store_Unit_Attributes_IRATE([unitsVelocity],mergerTrees,thisDataset)
             call thisDataset%close()
          end if
          ! Destroy the snapshot indices.
          call deallocateArray(thisSnapshotIndices)
          call deallocateArray(nodeSnapshotIndices)

          ! Close the groups.
          call darkParticlesGroup%close()
          call particlesGroup    %close()
          call snapshotGroup     %close()

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
          select case (mergerTrees%metaData(iAttribute)%metadataType)
          case (metaDataTypeCosmology  )
             attributeGroup => cosmologyGroup
          case (metaDataTypeSimulation )
             attributeGroup => simulationGroup
          end select

          ! Check if the group was recognized.
          if (associated(attributeGroup)) then

             ! Perform dictionary mapping from our internal names (which follow Galacticus format) to IRATE names.
             select case (mergerTrees%metaData(iAttribute)%metadataType)
             case (metaDataTypeCosmology  )
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
       call cosmologyGroup  %close()
       call simulationGroup %close()
    end if

    ! Close the output file.
    call outputFile%close()
    call hdf5Access%unset()

    return
  end subroutine Merger_Tree_Data_Structure_Export_IRATE

  subroutine Store_Unit_Attributes_Galacticus(unitType,unitLabel,mergerTrees,unitsGroup)
    !% Store attributes describing the unit system.
    use :: IO_HDF5, only : hdf5Object
    implicit none
    integer                  , intent(in   ) :: unitType
    character(len=*         ), intent(in   ) :: unitLabel
    class    (mergerTreeData), intent(in   ) :: mergerTrees
    type     (hdf5Object    ), intent(inout) :: unitsGroup

    call unitsGroup%writeAttribute(mergerTrees%units(unitType)%unitsInSI          ,unitLabel//"UnitsInSI"          )
    call unitsGroup%writeAttribute(mergerTrees%units(unitType)%hubbleExponent     ,unitLabel//"HubbleExponent"     )
    call unitsGroup%writeAttribute(mergerTrees%units(unitType)%scaleFactorExponent,unitLabel//"ScaleFactorExponent")
    return
  end subroutine Store_Unit_Attributes_Galacticus

  subroutine Store_Unit_Attributes_IRATE(unitType,mergerTrees,thisDataset)
    !% Store unit attributes in IRATE format files.
    use :: IO_HDF5                     , only : hdf5Object
    use :: ISO_Varying_String          , only : assignment(=), operator(//)
    use :: Numerical_Constants_Prefixes, only : hecto        , kilo
    implicit none
    integer                         , dimension(:), intent(in   ) :: unitType
    class           (mergerTreeData)              , intent(in   ) :: mergerTrees
    type            (hdf5Object    )              , intent(inout) :: thisDataset
    integer                                                       :: iUnit
    double precision                                              :: cgsUnits   , hubbleExponent, scaleFactorExponent
    type            (varying_string)                              :: unitName

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
    use :: Galacticus_Error, only : Galacticus_Error_Report
    class(mergerTreeData), intent(in   ) :: mergerTrees

    if (.not.mergerTrees%hasForestIndex    ) call Galacticus_Error_Report("merger trees do not have required property 'forestIndex'"    //{introspection:location})
    if (.not.mergerTrees%hasNodeIndex      ) call Galacticus_Error_Report("merger trees do not have required property 'nodeIndex'"      //{introspection:location})
    if (.not.mergerTrees%hasDescendentIndex) call Galacticus_Error_Report("merger trees do not have required property 'descendentIndex'"//{introspection:location})
    if (.not.mergerTrees%hasRedshift       ) call Galacticus_Error_Report("merger trees do not have required property 'redshift'"       //{introspection:location})
    if (.not.mergerTrees%hasNodeMass       ) call Galacticus_Error_Report("merger trees do not have required property 'nodeMass'"       //{introspection:location})
    return
  end subroutine Merger_Tree_Data_Validate_Trees

  subroutine Merger_Tree_Data_Set_Subhalo_Masses(mergerTrees)
    !% Set the masses of any subhalos (which have zero mass by default) based on particle count.
    use :: Galacticus_Error, only : Galacticus_Warn
    class(mergerTreeData), intent(inout) :: mergerTrees

    if (mergerTrees%hasParticleCount) then
       where (mergerTrees%nodeMass <= 0.0d0)
          mergerTrees%nodeMass=dble(mergerTrees%particleCount)*mergerTrees%particleMass
       end where
    end if
    if (any(mergerTrees%nodeMass <= 0.0d0)) call Galacticus_Warn("WARNING: some nodes have non-positive mass"//{introspection:location})
    return
  end subroutine Merger_Tree_Data_Set_Subhalo_Masses

  subroutine Merger_Tree_Data_Construct_Particle_Indices(mergerTrees)
    !% If we have most-bound particle indices and particle data has been read, construct arrays giving position of particle data for each node.
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Memory_Management, only : allocateArray          , deallocateArray
    class  (mergerTreeData), intent(inout) :: mergerTrees
    logical                                :: foundParticleData
    integer                                :: iNode            , iParticle

    if (mergerTrees%hasMostBoundParticleIndex) then
       ! Insist on having particle data.
       if (.not.mergerTrees%hasParticles) call Galacticus_Error_Report("most bound particle IDs provided, but no particle data was read"//{introspection:location})
       ! Allocate arrays for storing indices.
       if (allocated(mergerTrees%particleReferenceStart)) call deallocateArray(mergerTrees%particleReferenceStart)
       if (allocated(mergerTrees%particleReferenceCount)) call deallocateArray(mergerTrees%particleReferenceCount)
       call allocateArray(mergerTrees%particleReferenceStart,[mergerTrees%nodeCount])
       call allocateArray(mergerTrees%particleReferenceCount,[mergerTrees%nodeCount])
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
