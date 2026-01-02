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

  !!{
  An implementation of the merger tree importer class for ``Sussing Merger Trees'' format merger tree files.
  !!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: ISO_Varying_String      , only : varying_string
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  ! Enumeration of bad value test options
  !![
  <enumeration>
   <name>sussingBadValueTest</name>
   <description>Bad value test options.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="lessThan"   />
   <entry label="greaterThan"/>
  </enumeration>
  !!]

  ! Enumeration of halo mass definitions.
  !![
  <enumeration>
   <name>sussingMassOption</name>
   <description>Halo mass definitions.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="default"/>
   <entry label="FoF"    />
   <entry label="200Mean"/>
   <entry label="200Crit"/>
   <entry label="topHat" />
  </enumeration>
  !!]

  !![
  <mergerTreeImporter name="mergerTreeImporterSussing" abstract="yes">
   <description>Importer for ``Sussing Merger Trees'' format merger tree files \citep{srisawat_sussing_2013}.</description>
  </mergerTreeImporter>
  !!]
  type, extends(mergerTreeImporterClass) :: mergerTreeImporterSussing
     !!{
     A merger tree importer class for ``Sussing Merger Trees'' format merger tree files \citep{srisawat_sussing_2013}.
     !!}
     private
     class           (cosmologyParametersClass        ), pointer                    :: cosmologyParameters_     => null()
     class           (cosmologyFunctionsClass         ), pointer                    :: cosmologyFunctions_      => null()
     class           (randomNumberGeneratorClass      ), pointer                    :: randomNumberGenerator_   => null()
     logical                                                                        :: fatalMismatches                   , treeIndicesRead    , &
          &                                                                            scaleRadiiAvailableValue          , spinsAvailableValue, &
          &                                                                            fatalNonTreeNode
     integer                                                                        :: treesCount
     type            (enumerationSussingMassOptionType)                             :: massOption
     double precision                                                               :: boxLength
     type            (importerUnits                    )                            :: boxLengthUnits
     type            (varying_string                   )                            :: mergerTreeFile
     type            (varying_string                   ), allocatable, dimension(:) :: snapshotFileName
     integer         (kind_int8                        ), allocatable, dimension(:) :: treeIndices
     integer         (c_size_t                         ), allocatable, dimension(:) :: treeIndexRanks
     integer         (c_size_t                         ), allocatable, dimension(:) :: treeSizes                         , treeBegins
     type            (nodeData                         ), allocatable, dimension(:) :: nodes
     double precision                                   , allocatable, dimension(:) :: snapshotTimes
     integer                                                                        :: subvolumeCount
     integer                                                         , dimension(3) :: subvolumeIndex
     double precision                                                               :: subvolumeBuffer
     double precision                                                               :: badValue
     type            (enumerationSussingBadValueTestType)                           :: badValueTest
     double precision                                                               :: treeSampleRate
   contains
     !![
     <methods>
       <method description="Load the halo data." method="load" />
       <method description="Return true if the given {\normalfont \ttfamily x,y,z} position lies within the current subvolume (plus the buffer region if {\normalfont \ttfamily buffered} is true." method="inSubvolume" />
       <method description="Return true if the given {\normalfont \ttfamily x} position lies within the {\normalfont \ttfamily iSubvolume}$^\mathrm{th}$ subvolume (plus the buffer region if {\normalfont \ttfamily buffered} is true." method="inSubvolume1D" />
       <method description="Return true if the given {\normalfont \ttfamily x} value is bad." method="valueIsBad" />
     </methods>
     !!]
     final     ::                                  sussingDestructor
     procedure :: load                          => sussingLoad
     procedure :: close                         => sussingClose
     procedure :: canReadSubsets                => sussingCanReadSubsets
     procedure :: treesHaveSubhalos             => sussingTreesHaveSubhalos
     procedure :: massesIncludeSubhalos         => sussingMassesIncludeSubhalos
     procedure :: angularMomentaIncludeSubhalos => sussingAngularMomentaIncludeSubhalos
     procedure :: treesAreSelfContained         => sussingTreesAreSelfContained
     procedure :: velocitiesIncludeHubbleFlow   => sussingVelocitiesIncludeHubbleFlow
     procedure :: positionsArePeriodic          => sussingPositionsArePeriodic
     procedure :: cubeLength                    => sussingCubeLength
     procedure :: treeWeight                    => sussingTreeWeight
     procedure :: treeCount                     => sussingTreeCount
     procedure :: treeIndex                     => sussingTreeIndex
     procedure :: nodeCount                     => sussingNodeCount
     procedure :: positionsAvailable            => sussingPositionsAvailable
     procedure :: scaleRadiiAvailable           => sussingScaleRadiiAvailable
     procedure :: particleCountAvailable        => sussingParticleCountAvailable
     procedure :: velocityMaximumAvailable      => sussingVelocityMaximumAvailable
     procedure :: velocityDispersionAvailable   => sussingVelocityDispersionAvailable
     procedure :: angularMomentaAvailable       => sussingAngularMomentaAvailable
     procedure :: angularMomenta3DAvailable     => sussingAngularMomenta3DAvailable
     procedure :: spinAvailable                 => sussingSpinAvailable
     procedure :: spin3DAvailable               => sussingSpin3DAvailable
     procedure :: import                        => sussingImport
     procedure :: subhaloTrace                  => sussingSubhaloTrace
     procedure :: subhaloTraceCount             => sussingSubhaloTraceCount
     procedure :: inSubvolume                   => sussingInSubvolume
     procedure :: inSubvolume1D                 => sussingInSubvolume1D
     procedure :: valueIsBad                    => sussingValueIsBad
  end type mergerTreeImporterSussing

  interface mergerTreeImporterSussing
     module procedure sussingConstructorParameters
     module procedure sussingConstructorInternal
  end interface mergerTreeImporterSussing

contains

  function sussingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeImporterSussing} format \citep{srisawat_sussing_2013} merger tree importer which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeImporterSussing )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    integer                                     , dimension(3)  :: subvolumeIndex
    logical                                                     :: fatalMismatches       , fatalNonTreeNode
    integer                                                     :: subvolumeCount
    double precision                                            :: subvolumeBuffer       , badValue        , &
         &                                                         treeSampleRate
    type            (varying_string            )                :: badValueTestText      , massOptionText

    !![
    <inputParameter>
      <name>fatalMismatches</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether mismatches in cosmological parameter values between \glc\ and ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree files should be considered fatal.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fatalNonTreeNode</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether nodes in snapshot files but not in the merger tree file should be considered fatal when importing from the ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>subvolumeCount</name>
      <defaultValue>1</defaultValue>
      <description>Specifies the number of subvolumes \emph{along each axis} into which a ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree files should be split for processing through \glc.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>subvolumeBuffer</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Specifies the buffer region (in units of Mpc$/h$ to follow the format convention) around subvolumes of a ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree file which should be read in to ensure that no halos are missed from trees.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>subvolumeIndex</name>
      <defaultValue>[0,0,0]</defaultValue>
      <description>Specifies the index (in each dimension) of the subvolume of a ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree file to process. Indices range from 0 to {\normalfont \ttfamily [subvolumeCount]}$-1$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>badValue</name>
      <defaultValue>-0.5d0</defaultValue>
      <description>Use for bad value detection in ``Sussing'' merger trees. Values for scale radius and halo spin which exceed this threshold are assumed to be bad.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>badValueTest</name>
      <defaultValue>var_str('lessThan')</defaultValue>
      <description>Use for bad value detection in ``Sussing'' merger trees. Values which exceed the threshold in ths specified direction are assumed to be bad.</description>
      <source>parameters</source>
      <variable>badValueTestText</variable>
    </inputParameter>
    <inputParameter>
      <name>treeSampleRate</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specify the probability that any given tree should processed (to permit subsampling).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massOptions</name>
      <defaultValue>var_str('default')</defaultValue>
      <description>Mass option for Sussing merger trees.</description>
      <source>parameters</source>
      <variable>massOptionText</variable>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"    name="cosmologyParameters_"    source="parameters"/>
    <objectBuilder class="cosmologyFunctions"     name="cosmologyFunctions_"     source="parameters"/>
    <objectBuilder class="randomNumberGenerator"  name="randomNumberGenerator_"  source="parameters"/>
    !!]
    self=mergerTreeImporterSussing(                                                                                     &
         &                         fatalMismatches                                                                    , &
         &                         fatalNonTreeNode                                                                   , &
         &                         subvolumeCount                                                                     , &
         &                         subvolumeBuffer                                                                    , &
         &                         subvolumeIndex                                                                     , &
         &                         badValue                                                                           , &
         &                         enumerationSussingBadValueTestEncode(char(badValueTestText),includesPrefix=.false.), &
         &                         treeSampleRate                                                                     , &
         &                         enumerationSussingMassOptionEncode  (char(massOptionText  ),includesPrefix=.false.), &
         &                         cosmologyParameters_                                                               , &
         &                         cosmologyFunctions_                                                                , &
         &                         randomNumberGenerator_                                                               &
         &                        )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function sussingConstructorParameters

  function sussingConstructorInternal(fatalMismatches,fatalNonTreeNode,subvolumeCount,subvolumeBuffer,subvolumeIndex,badValue,badValueTest,treeSampleRate,massOption,cosmologyParameters_,cosmologyFunctions_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeImporterSussing} format \citep{srisawat_sussing_2013} merger tree importer class.
    !!}
    implicit none
    type            (mergerTreeImporterSussing         )                              :: self
    integer                                             , intent(in   ), dimension(3) :: subvolumeIndex
    logical                                             , intent(in   )               :: fatalMismatches       , fatalNonTreeNode
    integer                                             , intent(in   )               :: subvolumeCount
    type            (enumerationSussingBadValueTestType), intent(in   )               :: badValueTest
    type            (enumerationSussingMassOptionType  ), intent(in   )               :: massOption
    double precision                                    , intent(in   )               :: subvolumeBuffer       , badValue        , &
         &                                                                               treeSampleRate
    class           (cosmologyParametersClass          ), intent(in   ), target       :: cosmologyParameters_
    class           (cosmologyFunctionsClass           ), intent(in   ), target       :: cosmologyFunctions_
    class           (randomNumberGeneratorClass        ), intent(in   ), target       :: randomNumberGenerator_
    !![
    <constructorAssign variables="fatalMismatches,fatalNonTreeNode,subvolumeCount,subvolumeBuffer,subvolumeIndex,badValue,badValueTest,treeSampleRate,massOption,*cosmologyParameters_,*cosmologyFunctions_, *randomNumberGenerator_"/>
    !!]

    self%treeIndicesRead         =.false.
    self%scaleRadiiAvailableValue=.true.
    return
  end function sussingConstructorInternal

  subroutine sussingDestructor(self)
    implicit none
    type(mergerTreeImporterSussing), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine sussingDestructor

  subroutine sussingClose(self)
    !!{
    Close a {\normalfont \ttfamily sussing} format merger tree file.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine sussingClose

  logical function sussingCanReadSubsets(self)
    !!{
    Return false since this format does not permit reading of arbitrary subsets of halos from a forest.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingCanReadSubsets=.false.
    return
  end function sussingCanReadSubsets

  integer function sussingTreesHaveSubhalos(self)
    !!{
    Return a Boolean integer specifying whether or not the trees have subhalos.
    !!}
    use :: Numerical_Constants_Boolean, only : booleanTrue
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingTreesHaveSubhalos=booleanTrue
    return
  end function sussingTreesHaveSubhalos

  logical function sussingMassesIncludeSubhalos(self)
    !!{
    Return a Boolean specifying whether or not the halo masses include the contribution from subhalos.
    !!}
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingMassesIncludeSubhalos=.true.
    return
  end function sussingMassesIncludeSubhalos

  logical function sussingAngularMomentaIncludeSubhalos(self)
    !!{
    Return a Boolean specifying whether or not the halo angular momenta include the contribution from subhalos.
    !!}
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingAngularMomentaIncludeSubhalos=.true.
    return
  end function sussingAngularMomentaIncludeSubhalos

  integer function sussingTreesAreSelfContained(self)
    !!{
    Return a Boolean integer specifying whether or not the trees are self-contained.
    !!}
    use :: Numerical_Constants_Boolean, only : booleanTrue
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingTreesAreSelfContained=booleanTrue
    return
  end function sussingTreesAreSelfContained

  integer function sussingVelocitiesIncludeHubbleFlow(self)
    !!{
    Return a Boolean integer specifying whether or not velocities include the Hubble flow.
    !!}
    use :: Numerical_Constants_Boolean, only : booleanFalse
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingVelocitiesIncludeHubbleFlow=booleanFalse
    return
  end function sussingVelocitiesIncludeHubbleFlow

  integer function sussingPositionsArePeriodic(self)
    !!{
    Return a Boolean integer specifying whether or not positions are periodic.
    !!}
    use :: Numerical_Constants_Boolean, only : booleanTrue
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingPositionsArePeriodic=booleanTrue
    return
  end function sussingPositionsArePeriodic

  double precision function sussingCubeLength(self,time,status)
    !!{
    Return the length of the simulation cube.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Boolean     , only : booleanTrue
    implicit none
    class           (mergerTreeImporterSussing), intent(inout)           :: self
    double precision                           , intent(in   )           :: time
    integer                                    , intent(  out), optional :: status

    sussingCubeLength=importerUnitConvert(self%boxLength,time,self%boxLengthUnits,megaParsec,self%cosmologyParameters_,self%cosmologyFunctions_)
    if (present(status)) status=booleanTrue
    return
  end function sussingCubeLength

  integer(kind=c_size_t) function sussingTreeCount(self)
    !!{
    Return a count of the number of trees available.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    call sussingTreeIndicesRead(self)
    sussingTreeCount=self%treesCount
    return
  end function sussingTreeCount

  integer(kind=c_size_t) function sussingTreeIndex(self,i)
    !!{
    Return the index of the $i^\mathrm{th}$ tree.
    !!}
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call sussingTreeIndicesRead(self)
    sussingTreeIndex=self%treeIndices(i)
    return
  end function sussingTreeIndex

  function sussingNodeCount(self,i)
    !!{
    Return a count of the number of nodes in the $i^\mathrm{th}$ tree.
    !!}
    implicit none
    integer(c_size_t                 )                :: sussingNodeCount
    class  (mergerTreeImporterSussing), intent(inout) :: self
    integer                           , intent(in   ) :: i

    call sussingTreeIndicesRead(self)
    sussingNodeCount=self%treeSizes(i)
    return
  end function sussingNodeCount

  subroutine sussingTreeIndicesRead(self)
    !!{
    Read the tree indices.
    !!}
    use            :: Arrays_Search                   , only : searchIndexed
    use            :: Display                         , only : displayCounter         , displayCounterClear, displayIndent        , displayMessage, &
          &                                                    displayUnindent        , displayVerbosity   , verbosityLevelWorking
    use            :: Error                           , only : Error_Report
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Kind_Numbers                    , only : kind_int8
    use            :: Numerical_Constants_Astronomical, only : massSolar              , megaParsec
    use            :: Numerical_Constants_Prefixes    , only : kilo
    use            :: Sorting                         , only : sortIndex
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class    (mergerTreeImporterSussing), intent(inout)                 :: self
    integer  (kind=kind_int8           ), allocatable  , dimension(:  ) :: nodeSelfIndices            , nodeTreeIndices
    integer  (kind=c_size_t            ), allocatable  , dimension(:  ) :: nodeIndexRanks             , nodeDescendantLocations
    logical                             , allocatable  , dimension(:  ) :: nodeIncomplete
    integer                             , parameter                     :: stepCountMaximum   =1000000
    integer                                                             :: descendantStepCount        , hostStepCount
    logical                                                             :: treeIndicesAssigned        , branchJumpCheckRequired
    integer  (kind=c_size_t            )                                :: i                          , j                      , &
         &                                                                 k                          , l                      , &
         &                                                                 iStart                     , nodeCountTrees
    character(len=32                   )                                :: label
    type     (importerUnits            )                                :: massUnits                  , lengthUnits            , &
         &                                                                 velocityUnits
    type     (varying_string           )                                :: message
    integer  (kind=kind_int8           )                                :: treeIndexPrevious          , treeIndexCurrent       , &
         &                                                                 treeIndexTo                , treeIndexFrom

    ! Return if indices have been read previously.
    if (self%treeIndicesRead) return
    self%treeIndicesRead=.true.
    ! Load the data.
    call self%load(                         &
         &         nodeSelfIndices        , &
         &         nodeIndexRanks         , &
         &         nodeDescendantLocations, &
         &         nodeIncomplete         , &
         &         nodeCountTrees         , &
         &         nodeTreeIndices        , &
         &         treeIndicesAssigned    , &
         &         branchJumpCheckRequired, &
         &         massUnits              , &
         &         lengthUnits            , &
         &         velocityUnits            &
         &        )
    ! Search for and resolve hosting loops.
    call displayIndent('Resolving hosting loops',verbosityLevelWorking)
    do i=1,nodeCountTrees
       if (self%nodes(i)%hostIndex /= self%nodes(i)%nodeIndex) then
          l=searchIndexed(nodeSelfIndices,nodeIndexRanks,self%nodes(i)%hostIndex)
          ! Detect missing host.
          if (l < 1 .or. l > nodeCountTrees) cycle
          l=nodeIndexRanks(l)
          if (self%nodes(l)%nodeIndex /= self%nodes(i)%hostIndex) cycle
          ! Detect loops.
          if (self%nodes(l)%hostIndex == self%nodes(i)%nodeIndex) then
             ! Hosting loop detected. Reset the more massive halo to be self-hosting.
             if (self%nodes(l)%nodeMass > self%nodes(i)%nodeMass) then
                self%nodes(l)%hostIndex=self%nodes(l)%nodeIndex
             else
                self%nodes(i)%hostIndex=self%nodes(i)%nodeIndex
             end if
          end if
       end if
    end do
    call displayUnindent('done',verbosityLevelWorking)
    ! Search for deep hosting hierarchies and reset to single level hierarchies.
    call displayIndent('Resolving deep hosting',verbosityLevelWorking)
    do i=1,nodeCountTrees
       if (self%nodes(i)%hostIndex /= self%nodes(i)%nodeIndex) then
          ! Find the host.
          l=searchIndexed(nodeSelfIndices,nodeIndexRanks,self%nodes(i)%hostIndex)
          ! Detect missing host.
          if (l < 1 .or. l > nodeCountTrees) cycle
          l=nodeIndexRanks(l)
          if (self%nodes(l)%nodeIndex /= self%nodes(i)%hostIndex) cycle
          ! Detect hosted host.
          hostStepCount=0
          do while (self%nodes(l)%hostIndex /= self%nodes(l)%nodeIndex)
             ! Find the host.
             j=searchIndexed(nodeSelfIndices,nodeIndexRanks,self%nodes(l)%hostIndex)
             ! Detect missing host.
             if (j < 1 .or. j > nodeCountTrees) exit
             j=nodeIndexRanks(j)
             if (self%nodes(j)%nodeIndex /= self%nodes(l)%hostIndex) exit
             ! Detect infinite loops.
             hostStepCount=hostStepCount+1
             if (hostStepCount > stepCountMaximum) then
                message='infinite (or at least, very large) loop detect in halo hosting [#2]'
                message=message//char(10)//' hosted node index: '//self%nodes(l)%nodeIndex
                message=message//char(10)//'   host node index: '//self%nodes(j)%nodeIndex
                write (label,'(f7.4)') self%nodes(l)%nodeTime
                message=message//char(10)//'  hosted node time: '//label
                write (label,'(f7.4)') self%nodes(j)%nodeTime
                message=message//char(10)//'    host node time: '//label
                call Error_Report(message//{introspection:location})
             end if
             ! Move to the new host.
             l=j
          end do
          self%nodes(i)%hostIndex=self%nodes(l)%nodeIndex
       end if
    end do
    call displayUnindent('done',verbosityLevelWorking)
    ! Assign tree indices.
    call displayMessage('Assigning tree indices',verbosityLevelWorking)
    if (.not.treeIndicesAssigned) then
       nodeTreeIndices=nodeSelfIndices
       do i=1,nodeCountTrees
          call displayCounter(int(100.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityLevelWorking)
          l=i
          hostStepCount=0
          do while (l /= -1 .and. self%nodes(l)%hostIndex /= self%nodes(l)%nodeIndex)
             ! Find the host halo.
             k=searchIndexed(nodeSelfIndices,nodeIndexRanks,self%nodes(l)%hostIndex)
             ! Check for missing hosts.
             if (k < 1 .or. k > nodeCountTrees .or. self%nodes(l)%hostIndex /= self%nodes(nodeIndexRanks(k))%nodeIndex) then
                ! No host can be found (it must be outside of the buffered subvolume). Assign this node its own index as a tree
                ! index, and skip looking for hosts.
                nodeTreeIndices(i)=self%nodes(i)%nodeIndex
                exit
             end if
             ! Perform sanity checks.
             if (k >= 0 .and. self%nodes(nodeIndexRanks(k))%nodeTime /= self%nodes(l)%nodeTime) then
                message='host exists at different time from hosted node'
                message=message//char(10)//' hosted node index: '//self%nodes(               l )%nodeIndex
                message=message//char(10)//'   host node index: '//self%nodes(nodeIndexRanks(k))%nodeIndex
                write (label,'(f7.4)') self%nodes(               l )%nodeTime
                message=message//char(10)//'  hosted node time: '//label
                write (label,'(f7.4)') self%nodes(nodeIndexRanks(k))%nodeTime
                message=message//char(10)//'    host node time: '//label
                call Error_Report(message//{introspection:location})
             end if
             hostStepCount=hostStepCount+1
             if (hostStepCount > stepCountMaximum) then
                message='infinite (or at least, very large) loop detect in halo hosting'
                message=message//char(10)//' hosted node index: '//self%nodes(               l )%nodeIndex
                message=message//char(10)//'   host node index: '//self%nodes(nodeIndexRanks(k))%nodeIndex
                write (label,'(f7.4)') self%nodes(               l )%nodeTime
                message=message//char(10)//'  hosted node time: '//label
                write (label,'(f7.4)') self%nodes(nodeIndexRanks(k))%nodeTime
                message=message//char(10)//'    host node time: '//label
                call Error_Report(message//{introspection:location})
             end if
             l=nodeIndexRanks(k)
          end do
          descendantStepCount=0
          do while (l /= -1)
             nodeTreeIndices(i)=nodeTreeIndices(l)
             k=nodeDescendantLocations(l)
             ! Check for missing descendants.
             if (k < 0 .or. self%nodes(l)%descendantIndex /= self%nodes(k)%nodeIndex) exit
             ! Perform sanity checks.
             if (k >= 0 .and. self%nodes(k)%nodeTime <= self%nodes(l)%nodeTime) then
                message='descendant exists before progenitor node'
                message=message//char(10)//' progenitor node index: '//self%nodes(l)%nodeIndex
                message=message//char(10)//' descendant node index: '//self%nodes(k)%nodeIndex
                write (label,'(f7.4)') self%nodes(l)%nodeTime
                message=message//char(10)//'  progenitor node time: '//label
                write (label,'(f7.4)') self%nodes(k)%nodeTime
                message=message//char(10)//'  descendant node time: '//label
                call Error_Report(message//{introspection:location})
             end if
             descendantStepCount=descendantStepCount+1
             if (descendantStepCount > stepCountMaximum) then
                message='infinite (or at least, very large) loop detect in halo descent'
                message=message//char(10)//' progenitor node index: '//self%nodes(l)%nodeIndex
                message=message//char(10)//' descendant node index: '//self%nodes(k)%nodeIndex
                write (label,'(f7.4)') self%nodes(l)%nodeTime
                message=message//char(10)//'  progenitor node time: '//label
                write (label,'(f7.4)') self%nodes(k)%nodeTime
                message=message//char(10)//'  descendant node time: '//label
                call Error_Report(message//{introspection:location})
             end if
             l=k
             hostStepCount=0
             do while (l /= -1 .and. self%nodes(l)%hostIndex /= self%nodes(l)%nodeIndex)
                k=searchIndexed(nodeSelfIndices,nodeIndexRanks,self%nodes(l)%hostIndex)
                ! Check for missing hosts.
                if (self%nodes(l)%hostIndex /= self%nodes(nodeIndexRanks(k))%nodeIndex) exit
                ! Perform sanity checks.
                if (k >= 0 .and. self%nodes(nodeIndexRanks(k))%nodeTime /= self%nodes(l)%nodeTime) then
                   message='host exists at different time from hosted node'
                   message=message//char(10)//' hosted node index: '//self%nodes(               l )%nodeIndex
                   message=message//char(10)//'   host node index: '//self%nodes(nodeIndexRanks(k))%nodeIndex
                   write (label,'(f7.4)') self%nodes(               l )%nodeTime
                   message=message//char(10)//'  hosted node time: '//label
                   write (label,'(f7.4)') self%nodes(nodeIndexRanks(k))%nodeTime
                   message=message//char(10)//'    host node time: '//label
                   call Error_Report(message//{introspection:location})
                end if
                hostStepCount=hostStepCount+1
                if (hostStepCount > stepCountMaximum) then
                   message='infinite (or at least, very large) loop detect in halo hosting'
                   message=message//char(10)//' hosted node index: '//self%nodes(               l )%nodeIndex
                   message=message//char(10)//'   host node index: '//self%nodes(nodeIndexRanks(k))%nodeIndex
                   write (label,'(f7.4)') self%nodes(               l )%nodeTime
                   message=message//char(10)//'  hosted node time: '//label
                   write (label,'(f7.4)') self%nodes(nodeIndexRanks(k))%nodeTime
                   message=message//char(10)//'    host node time: '//label
                   call Error_Report(message//{introspection:location})
                end if
                l=nodeIndexRanks(k)
             end do
          end do
       end do
    end if
    call displayCounterClear(verbosityLevelWorking)
    ! Check for nodes jumping between trees and join any such trees.
    if (branchJumpCheckRequired) then
       call displayMessage('Checking for branch jumps',verbosityLevelWorking)
       do i=1,nodeCountTrees
          call displayCounter(int(100.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityLevelWorking)
          l=nodeDescendantLocations(i)
          if (l /= -1) then
             if (nodeTreeIndices(i) /= nodeTreeIndices(l)) then
                ! Merge the trees by assigning the higher tree index to all nodes in the other tree.
                if (nodeTreeIndices(i) > nodeTreeIndices(l)) then
                   treeIndexFrom=nodeTreeIndices(l)
                   treeIndexTo  =nodeTreeIndices(i)
                else
                   treeIndexFrom=nodeTreeIndices(i)
                   treeIndexTo  =nodeTreeIndices(l)
                end if
                !$omp parallel do
                do j=1,nodeCountTrees
                   if (nodeTreeIndices(j) == treeIndexFrom) nodeTreeIndices(j)=treeIndexTo
                end do
                !$omp end parallel do
             end if
          end if
       end do
       do i=1,nodeCountTrees
          l=nodeDescendantLocations(i)
          if (l /= -1) then
             if (nodeTreeIndices(i) /= nodeTreeIndices(l))                                    &
                  & call Error_Report('failed to cross-link trees'//{introspection:location})
          end if
       end do
       call displayCounterClear(verbosityLevelWorking)
    end if
    ! Generate an index into nodes sorted by tree index.
    self%treeIndexRanks=sortIndex(nodeTreeIndices)
    ! Identify trees which contain incomplete nodes.
    call displayMessage('Checking for incomplete trees',verbosityLevelWorking)
    i               =0
    iStart          =0
    treeIndexCurrent=-1
    do while (i < nodeCountTrees)
       i     =i+1
       ! Check if we started a new tree, and record the index at which it began.
       if (treeIndexCurrent /= nodeTreeIndices(self%treeIndexRanks(i))) iStart=i
       treeIndexCurrent=nodeTreeIndices(self%treeIndexRanks(i))
       if (nodeIncomplete(self%treeIndexRanks(i))) then
          if (displayVerbosity() >= verbosityLevelWorking) then
             j=self%treeIndexRanks(i)
             message='Marking tree '
             message=message//treeIndexCurrent//' as incomplete due to node '//nodeSelfIndices(j)//' at position:'
             write (label,'(e12.6)') self%nodes(j)%position(1)
             message=message//char(10)//'  x: '//trim(label)
             write (label,'(e12.6)') self%nodes(j)%position(2)
             message=message//char(10)//'  y: '//trim(label)
             write (label,'(e12.6)') self%nodes(j)%position(3)
             message=message//char(10)//'  z: '//trim(label)
             call displayMessage(message,verbosityLevelWorking)
          end if
          ! Reset index to the start of this tree.
          i=iStart-1
          ! Mark all nodes in the tree as incomplete.
          do while (i < nodeCountTrees)
             i=i+1
             if (nodeTreeIndices(self%treeIndexRanks(i)) == treeIndexCurrent) then
                nodeIncomplete(self%treeIndexRanks(i))=.true.
             else
                i=i-1
                exit
             end if
          end do
       end if
       call displayCounter(int(50.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityLevelWorking)
    end do
    ! Check for zero mass halos.
    call displayMessage('Checking for zero mass halos',verbosityLevelWorking)
    do i=1,nodeCountTrees
       call displayCounter(int(100.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityLevelWorking)
       if (self%nodes(i)%nodeMass == 0.0d0) then
          k=i
          do while (self%nodes(i)%nodeMass == 0.0d0)
             if (self%nodes(k)%nodeMass > 0.0d0) then
                self%nodes(i)%nodeMass   =self%nodes(k)%nodeMass
                self%nodes(i)%scaleRadius=self%nodes(k)%scaleRadius
             end if
             k=nodeDescendantLocations(k)
             if (k < 0 .or. self%nodes(i)%descendantIndex /= self%nodes(k)%nodeIndex) then
                message='zero mass halo ['
                message=message//self%nodes(i)%nodeIndex//'] has no non-zero mass descendants'
                call Error_Report(message//{introspection:location})
             end if
          end do
       end if
    end do
    call displayCounterClear(verbosityLevelWorking)
    ! Check for incomplete trees not in the buffer zone.
    i               = 0
    iStart          = 0
    treeIndexCurrent=-1
    do while (i<nodeCountTrees)
       i=i+1
       j=self%treeIndexRanks(i)
       ! Check if we started a new tree, and record the index at which it began.
       if (treeIndexCurrent /= nodeTreeIndices(j)) iStart=i
       treeIndexCurrent=nodeTreeIndices(j)
       if (nodeTreeIndices(j) == nodeSelfIndices(j)) then
          if (self%inSubVolume(self%nodes(j)%position(1),self%nodes(j)%position(2),self%nodes(j)%position(3),buffered=.false.)) then
             ! Tree root is in subvolume. Check if the tree is incomplete.
             if (nodeIncomplete(j)) then
                message='tree in subvolume is incomplete - try increasing buffer size'
                message=message//char(10)//' tree index: '//treeIndexCurrent
                write (label,'(e12.6)') self%nodes(j)%position(1)
                message=message//char(10)//'          x: '//label
                write (label,'(e12.6)') self%nodes(j)%position(2)
                message=message//char(10)//'          y: '//label
                write (label,'(e12.6)') self%nodes(j)%position(3)
                message=message//char(10)//'          z: '//label
                call Error_Report(message//{introspection:location})
             end if
          else
             ! Reset index to the start of this tree.
             i=iStart-1
             ! Tree is not in subvolume. Set index to -1 to indicate that we should ignore it.
             do while (i < nodeCountTrees)
                i=i+1
                j=self%treeIndexRanks(i)
                if (nodeTreeIndices(j) == treeIndexCurrent) then
                   nodeTreeIndices(j)=-1
                else
                   i=i-1
                   exit
                end if
             end do
          end if
       end if
       call displayCounter(int(50.0d0+50.0d0*dble(i)/dble(nodeCountTrees)),.false.,verbosityLevelWorking)
    end do
    call displayCounterClear(verbosityLevelWorking)
    ! Generate an index into nodes sorted by tree index.
    self%treeIndexRanks=sortIndex(nodeTreeIndices)
    ! Create a list of tree indices, sizes, and start locations.
    call displayMessage('Generating tree list',verbosityLevelWorking)
    self%treesCount=0
    treeIndexPrevious=-1
    do i=1,nodeCountTrees
       call displayCounter(int(100.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityLevelWorking)
       if (nodeTreeIndices(self%treeIndexRanks(i)) /= treeIndexPrevious) then
          treeIndexPrevious=nodeTreeIndices(self%treeIndexRanks(i))
          self%treesCount=self%treesCount+1
       end if
    end do
    call displayCounterClear(verbosityLevelWorking)
    message='Found '
    message=message//self%treesCount//' trees'
    call displayMessage(message,verbosityLevelWorking)
    allocate(self%treeIndices(self%treesCount))
    allocate(self%treeSizes  (self%treesCount))
    allocate(self%treeBegins (self%treesCount))
    treeIndexPrevious=-1
    j                = 0
    self%treeSizes   = 0
    do i=1,nodeCountTrees
       if (nodeTreeIndices(self%treeIndexRanks(i)) /= treeIndexPrevious) then
          treeIndexPrevious=nodeTreeIndices(self%treeIndexRanks(i))
          j=j+1
          self%treeIndices(j)=nodeTreeIndices(self%treeIndexRanks(i))
          self%treeBegins (j)=i
       end if
       if (j > 0) self%treeSizes(j)=self%treeSizes(j)+1
    end do
    ! Clean up display.
    call displayCounterClear(verbosityLevelWorking)
    ! Do unit conversion.
    self   %nodes%nodeMass             =importerUnitConvert(self%nodes%nodeMass             ,self%nodes%nodeTime,massUnits    ,massSolar ,self%cosmologyParameters_,self%cosmologyFunctions_)
    self   %nodes%scaleRadius          =importerUnitConvert(self%nodes%scaleRadius          ,self%nodes%nodeTime,lengthUnits  ,megaParsec,self%cosmologyParameters_,self%cosmologyFunctions_)
    self   %nodes%velocityMaximum      =importerUnitConvert(self%nodes%velocityMaximum      ,self%nodes%nodeTime,velocityUnits,kilo      ,self%cosmologyParameters_,self%cosmologyFunctions_)
    self   %nodes%velocityDispersion   =importerUnitConvert(self%nodes%velocityDispersion   ,self%nodes%nodeTime,velocityUnits,kilo      ,self%cosmologyParameters_,self%cosmologyFunctions_)
    do i=1,3
       self%nodes%position          (i)=importerUnitConvert(self%nodes%position          (i),self%nodes%nodeTime,lengthUnits  ,megaParsec,self%cosmologyParameters_,self%cosmologyFunctions_)
       self%nodes%velocity          (i)=importerUnitConvert(self%nodes%velocity          (i),self%nodes%nodeTime,velocityUnits,kilo      ,self%cosmologyParameters_,self%cosmologyFunctions_)
    end do
    ! Destroy temporary workspace.
    deallocate(nodeSelfIndices        )
    deallocate(nodeTreeIndices        )
    deallocate(nodeDescendantLocations)
    ! Write completion message.
    call displayUnindent('done',verbosityLevelWorking)
   return
  end subroutine sussingTreeIndicesRead

  double precision function sussingTreeWeight(self,i)
    !!{
    Return the weight to assign to trees.
    !!}
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    integer                           , intent(in   ) :: i
    !$GLC attributes unused :: i

    ! Compute the inverse of the cube volume.
    sussingTreeWeight=1.0d0/self%cubeLength(self%cosmologyFunctions_%cosmicTime(1.0d0))**3/self%treeSampleRate
    return
  end function sussingTreeWeight

  logical function sussingPositionsAvailable(self,positions,velocities)
    !!{
    Return true if positions and/or velocities are available.
    !!}
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    logical                           , intent(in   ) :: positions, velocities
    !$GLC attributes unused :: self, positions, velocities

    sussingPositionsAvailable=.true.
    return
  end function sussingPositionsAvailable

  logical function sussingScaleRadiiAvailable(self)
    !!{
    Return true if scale radii are available.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    call sussingTreeIndicesRead(self)
    sussingScaleRadiiAvailable=self%scaleRadiiAvailableValue
    return
  end function sussingScaleRadiiAvailable

  logical function sussingParticleCountAvailable(self)
    !!{
    Return true if particle counts are available.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingParticleCountAvailable=.false.
    return
  end function sussingParticleCountAvailable

  logical function sussingVelocityMaximumAvailable(self)
    !!{
    Return true if halo rotation curve velocity maxima are available.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingVelocityMaximumAvailable=.true.
    return
  end function sussingVelocityMaximumAvailable

  logical function sussingVelocityDispersionAvailable(self)
    !!{
    Return true if halo velocity dispersions are available.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingVelocityDispersionAvailable=.true.
    return
  end function sussingVelocityDispersionAvailable

  logical function sussingAngularMomentaAvailable(self)
    !!{
    Return true if angular momenta are available.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingAngularMomentaAvailable=.false.
    return
  end function sussingAngularMomentaAvailable

  logical function sussingAngularMomenta3DAvailable(self)
    !!{
    Return true if angular momenta vectors are available.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    !$GLC attributes unused :: self

    sussingAngularMomenta3DAvailable=.false.
    return
  end function sussingAngularMomenta3DAvailable

  logical function sussingSpinAvailable(self)
    !!{
    Return true if spins are available.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    call sussingTreeIndicesRead(self)
    sussingSpinAvailable=self%spinsAvailableValue
    return
  end function sussingSpinAvailable

  logical function sussingSpin3DAvailable(self)
    !!{
    Return true if spins vectors are available.
    !!}
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    call sussingTreeIndicesRead(self)
    sussingSpin3DAvailable=self%spinsAvailableValue
    return
  end function sussingSpin3DAvailable

  subroutine sussingSubhaloTrace(self,node,time,position,velocity)
    !!{
    Returns a trace of subhalo position/velocity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (mergerTreeImporterSussing), intent(inout)                 :: self
    class           (nodeData                 ), intent(in   )                 :: node
    double precision                           , intent(  out), dimension(:  ) :: time
    double precision                           , intent(  out), dimension(:,:) :: position, velocity
    !$GLC attributes unused :: self, node, time, position, velocity

    call Error_Report('subhalo traces are not available'//{introspection:location})
    return
  end subroutine sussingSubhaloTrace

  function sussingSubhaloTraceCount(self,node)
    !!{
    Returns the length of a subhalo trace.
    !!}
    implicit none
    integer(c_size_t                 )                :: sussingSubhaloTraceCount
    class  (mergerTreeImporterSussing), intent(inout) :: self
    class  (nodeData                 ), intent(in   ) :: node
    !$GLC attributes unused :: self, node

    ! No particle data is available.
    sussingSubhaloTraceCount=0
    return
  end function sussingSubhaloTraceCount

  subroutine sussingImport(self,i,nodes,nodeSubset,requireScaleRadii,requireAngularMomenta,requireAngularMomenta3D,requireSpin,requireSpin3D,requirePositions,structureOnly,requireNamedReals,requireNamedIntegers)
    !!{
    Import the $i^\mathrm{th}$ merger tree.
    !!}
    use :: Error            , only : Error_Report
    implicit none
    class           (mergerTreeImporterSussing), intent(inout)                            :: self
    integer                                    , intent(in   )                            :: i
    class           (nodeDataMinimal          ), intent(  out), allocatable, dimension(:) :: nodes
    integer         (c_size_t                 ), intent(in   ), optional   , dimension(:) :: nodeSubset
    logical                                    , intent(in   ), optional                  :: requireScaleRadii         , requireAngularMomenta, &
         &                                                                                   requireAngularMomenta3D   , requirePositions     , &
         &                                                                                   structureOnly             , requireSpin          , &
         &                                                                                   requireSpin3D
    type            (varying_string           ), intent(in   ), optional   , dimension(:) :: requireNamedReals         , requireNamedIntegers
    integer         (c_size_t                 )                                           :: j
    !$GLC attributes unused :: requireAngularMomenta, requireAngularMomenta3D, requireScaleRadii, requirePositions, requireSpin, requireSpin3D, requireNamedReals, requireNamedIntegers

    ! Decide if this tree should be included.
    if (self%randomNumberGenerator_%uniformSample() <= self%treeSampleRate) then
       ! Validate arguments.
       if (present(structureOnly).and.    structureOnly                ) call Error_Report('import of structure only is not supported'//{introspection:location})
       if (present(nodeSubset   ).and.any(nodeSubset    /= -1_c_size_t)) call Error_Report('import of subsets is not supported'       //{introspection:location})
       ! Allocate the nodes array.
       allocate(nodeData :: nodes(self%treeSizes(i)))
       ! Copy data to nodes.
       select type (nodes)
       type is (nodeData)
          do j=1,self%treeSizes(i)
             nodes(j)=self%nodes(self%treeIndexRanks(self%treeBegins(i)+j-1))
          end do
       end select
    end if
    return
  end subroutine sussingImport

  logical function sussingInSubvolume(self,x,y,z,buffered)
    !!{
    Determine if a point lies within a subvolume of the simulation box (possibly with some buffering).
    !!}
    implicit none
    class           (mergerTreeImporterSussing), intent(inout) :: self
    double precision                           , intent(in   ) :: x,y,z
    logical                                    , intent(in   ) :: buffered

    sussingInSubvolume=                                                       &
         &              self%inSubvolume1D(x,self%subvolumeIndex(1),buffered) &
         &             .and.                                                  &
         &              self%inSubvolume1D(y,self%subvolumeIndex(2),buffered) &
         &             .and.                                                  &
         &              self%inSubvolume1D(z,self%subvolumeIndex(3),buffered)
    return
  end function sussingInSubvolume

  logical function sussingInSubvolume1D(self,x,iSubvolume,buffered)
    !!{
    Determine if a point lies within the 1-D range of a subvolume of the simulation box (possibly with some buffering).
    !!}
    use :: Numerical_Constants_Astronomical, only : kiloParsec, megaParsec
    implicit none
    class           (mergerTreeImporterSussing), intent(inout) :: self
    double precision                           , intent(in   ) :: x
    integer                                    , intent(in   ) :: iSubvolume
    logical                                    , intent(in   ) :: buffered
    double precision                                           :: subvolumeCenter, boxLength, &
         &                                                        buffer

    if (self%subvolumeCount <= 1) then
       sussingInSubvolume1D=.true.
    else
       boxLength=self%boxLength*megaParsec/kiloParsec
       if (buffered) then
          buffer=self%subvolumeBuffer*megaParsec/kiloParsec
       else
          buffer=0.0d0
       end if
       subvolumeCenter=boxLength*(dble(iSubvolume)+0.5d0)/dble(self%subvolumeCount)
       sussingInSubvolume1D=                                                             &
            &                abs(sussingPeriodicSeparation(x,subvolumeCenter,boxLength)) &
            &               <=                                                           &
            &                boxLength/2.0d0/dble(self%subvolumeCount)+buffer
    end if
    return
  end function sussingInSubvolume1D

  double precision function sussingPeriodicSeparation(x1,x2,periodicLength)
    !!{
    Determine the separation between two points in a periodic cube.
    !!}
    implicit none
    double precision, intent(in   ) :: x1,x2,periodicLength

    sussingPeriodicSeparation=x1-x2
    do while (sussingPeriodicSeparation > 0.5d0*periodicLength)
       sussingPeriodicSeparation=sussingPeriodicSeparation-periodicLength
    end do
    do while (sussingPeriodicSeparation < -0.5d0*periodicLength)
       sussingPeriodicSeparation=sussingPeriodicSeparation+periodicLength
    end do
    return
  end function sussingPeriodicSeparation

  logical function sussingValueIsBad(self,x)
    !!{
    Determine if a value in a ``Sussing'' merger tree file is bad
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (mergerTreeImporterSussing), intent(inout) :: self
    double precision                           , intent(in   ) :: x

    select case (self%badValueTest%ID)
    case (sussingBadValueTestLessThan   %ID)
       sussingValueIsBad=(x < self%badValue)
    case (sussingBadValueTestGreaterThan%ID)
       sussingValueIsBad=(x > self%badValue)
    case default
       sussingValueIsBad=.false.
       call Error_Report('unknown badness test'//{introspection:location})
    end select
    return
  end function sussingValueIsBad

  subroutine sussingLoad(self,nodeSelfIndices,nodeIndexRanks,nodeDescendantLocations,nodeIncomplete,nodeCountTrees,nodeTreeIndices,treeIndicesAssigned,branchJumpCheckRequired,massUnits,lengthUnits,velocityUnits)
    !!{
    Stub function for the {\normalfont \ttfamily load} method of the {\normalfont \ttfamily sussing} merger tree importer.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (mergerTreeImporterSussing), intent(inout)                            :: self
    integer(kind_int8                ), intent(  out), dimension(:), allocatable :: nodeSelfIndices    , nodeTreeIndices
    integer(c_size_t                 ), intent(  out), dimension(:), allocatable :: nodeIndexRanks     , nodeDescendantLocations
    logical                           , intent(  out), dimension(:), allocatable :: nodeIncomplete
    integer(kind=c_size_t            ), intent(  out)                            :: nodeCountTrees
    logical                           , intent(  out)                            :: treeIndicesAssigned, branchJumpCheckRequired
    type   (importerUnits            ), intent(  out)                            :: massUnits          , lengthUnits            , &
         &                                                                          velocityUnits
    !$GLC attributes unused :: self,nodeSelfIndices,nodeIndexRanks,nodeDescendantLocations,nodeIncomplete,nodeCountTrees,nodeTreeIndices,treeIndicesAssigned,branchJumpCheckRequired,massUnits,lengthUnits,velocityUnits

    call Error_Report('this function should not be accessed'//{introspection:location})
    return
  end subroutine sussingLoad

