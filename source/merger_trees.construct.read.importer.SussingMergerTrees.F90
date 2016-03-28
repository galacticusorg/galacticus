!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% An implementation of the merger tree importer class for ``Sussing Merger Trees'' format merger tree files.

  !# <mergerTreeImporter name="mergerTreeImporterSussing" abstract="yes" description="Importer for ``Sussing Merger Trees'' format merger tree files \citep{srisawat_sussing_2013}." />
  use Stateful_Types
  use ISO_Varying_String

  type, abstract, extends(mergerTreeImporterClass) :: mergerTreeImporterSussing
     !% A merger tree importer class for ``Sussing Merger Trees'' format merger tree files \citep{srisawat_sussing_2013}.
     private
     logical                                                     :: fatalMismatches         , treeIndicesRead    , &
          &                                                         scaleRadiiAvailableValue, spinsAvailableValue
     integer                                                     :: treesCount
     double precision                                            :: boxLength
     type            (importerUnits )                            :: boxLengthUnits
     type            (varying_string)                            :: mergerTreeFile
     type            (varying_string), allocatable, dimension(:) :: snapshotFileName
     integer         (kind=kind_int8), allocatable, dimension(:) :: treeIndices
     integer         (kind=c_size_t ), allocatable, dimension(:) :: treeIndexRanks
     integer         (c_size_t      ), allocatable, dimension(:) :: treeSizes               , treeBegins
     type            (nodeData      ), allocatable, dimension(:) :: nodes
     double precision                , allocatable, dimension(:) :: snapshotTimes
     integer                                                     :: subvolumeCount
     integer                                      , dimension(3) :: subvolumeIndex
     double precision                                            :: subvolumeBuffer
     logical                                                     :: convertToBinary         , binaryFormatOld
     double precision                                            :: badValue
     integer                                                     :: badTest
   contains
     !@ <objectMethods>
     !@   <object>mergerTreeImporterSussing</object>
     !@   <objectMethod>
     !@     <method>load</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>{\textcolor{red}{\textless integer(kind\_int8,:)\textgreater}}\ nodeSelfIndices\argout, {\textcolor{red}{\textless integer(c\_size\_t,:)\textgreater}}\ nodeIndexRanks\argout, {\textcolor{red}{\textless integer(c\_size\_t,:)\textgreater}}\ nodeDescendentLocations\argout, \logicalzero\ nodeIncomplete\argout,  {\textcolor{red}{\textless integer(c\_size\_t)\textgreater}}\ nodeCountTrees\argout, {\textcolor{red}{\textless integer(kind\_int8,:)\textgreater}}\ nodeTreeIndices\argout, \logicalzero\ treeIndicesAssigned\argout, \logicalzero\ branchJumpCheckRequired\argout, {\textcolor{red}{\textless type(importerUnits)\textgreater}}\ massUnits\argout, {\textcolor{red}{\textless type(importerUnits)\textgreater}}\ lengthUnits\argout, {\textcolor{red}{\textless type(importerUnits)\textgreater}}\ velocityUnits\argout</arguments>
     !@     <description>Load the halo data.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>inSubvolume</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\doublezero\ x\argin, \doublezero\ y\argin, \doublezero\ z\argin, \logicalzero\ [buffered]\argin</arguments>
     !@     <description>Return true if the given {\normalfont \ttfamily x,y,z} position lies within the current subvolume (plus the buffer region if {\normalfont \ttfamily buffered} is true.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>inSubvolume1D</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\doublezero\ x\argin, \intzero\ iSubvolume\argin, \logicalzero\ [buffered]\argin</arguments>
     !@     <description>Return true if the given {\normalfont \ttfamily x} position lies within the {\normalfont \ttfamily iSubvolume}$^{\mathrm th}$ subvolume (plus the buffer region if {\normalfont \ttfamily buffered} is true.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>valueIsBad</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\doublezero\ x\argin</arguments>
     !@     <description>Return true if the given {\normalfont \ttfamily x} value is bad.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(sussingLoad), deferred :: load
     procedure                        :: close                         => sussingClose
     procedure                        :: treesHaveSubhalos             => sussingTreesHaveSubhalos
     procedure                        :: massesIncludeSubhalos         => sussingMassesIncludeSubhalos
     procedure                        :: angularMomentaIncludeSubhalos => sussingAngularMomentaIncludeSubhalos
     procedure                        :: treesAreSelfContained         => sussingTreesAreSelfContained
     procedure                        :: velocitiesIncludeHubbleFlow   => sussingVelocitiesIncludeHubbleFlow
     procedure                        :: positionsArePeriodic          => sussingPositionsArePeriodic
     procedure                        :: cubeLength                    => sussingCubeLength
     procedure                        :: treeWeight                    => sussingTreeWeight
     procedure                        :: treeCount                     => sussingTreeCount
     procedure                        :: treeIndex                     => sussingTreeIndex
     procedure                        :: nodeCount                     => sussingNodeCount
     procedure                        :: positionsAvailable            => sussingPositionsAvailable
     procedure                        :: scaleRadiiAvailable           => sussingScaleRadiiAvailable
     procedure                        :: particleCountAvailable        => sussingParticleCountAvailable
     procedure                        :: velocityMaximumAvailable      => sussingVelocityMaximumAvailable
     procedure                        :: velocityDispersionAvailable   => sussingVelocityDispersionAvailable
     procedure                        :: angularMomentaAvailable       => sussingAngularMomentaAvailable
     procedure                        :: angularMomenta3DAvailable     => sussingAngularMomenta3DAvailable
     procedure                        :: spinAvailable                 => sussingSpinAvailable
     procedure                        :: spin3DAvailable               => sussingSpin3DAvailable
     procedure                        :: import                        => sussingImport
     procedure                        :: subhaloTrace                  => sussingSubhaloTrace
     procedure                        :: subhaloTraceCount             => sussingSubhaloTraceCount
     procedure                        :: inSubvolume                   => sussingInSubvolume
     procedure                        :: inSubvolume1D                 => sussingInSubvolume1D
     procedure                        :: valueIsBad                    => sussingValueIsBad
  end type mergerTreeImporterSussing

  abstract interface
     !% Interface for the {\normalfont \ttfamily load} method of the {\normalfont \ttfamily sussing} merger tree importer.
     subroutine sussingLoad(self,nodeSelfIndices,nodeIndexRanks,nodeDescendentLocations,nodeIncomplete,nodeCountTrees,nodeTreeIndices,treeIndicesAssigned,branchJumpCheckRequired,massUnits,lengthUnits,velocityUnits)
       import kind_int8, c_size_t, mergerTreeImporterSussing, importerUnits
       class  (mergerTreeImporterSussing), intent(inout)                            :: self
       integer(kind_int8                ), intent(  out), dimension(:), allocatable :: nodeSelfIndices    , nodeTreeIndices
       integer(c_size_t                 ), intent(  out), dimension(:), allocatable :: nodeIndexRanks     , nodeDescendentLocations
       logical                           , intent(  out), dimension(:), allocatable :: nodeIncomplete
       integer(kind=c_size_t            ), intent(  out)                            :: nodeCountTrees
       logical                           , intent(  out)                            :: treeIndicesAssigned, branchJumpCheckRequired
       type   (importerUnits            ), intent(  out)                            :: massUnits          , lengthUnits            , &
            &                                                                          velocityUnits
     end subroutine sussingLoad
  end interface
  
  ! Record of implementation initialization state.
  logical                        :: sussingInitialized                       =.false.

  ! Default settings.
  logical                        :: mergerTreeImportSussingMismatchIsFatal
  logical                        :: mergerTreeImportSussingNonTreeNodeIsFatal
  integer                        :: mergerTreeImportSussingSubvolumeCount
  double precision               :: mergerTreeImportSussingSubvolumeBuffer
  integer         , dimension(3) :: mergerTreeImportSussingSubvolumeIndex
  double precision               :: mergerTreeImportSussingBadValue
  integer                        :: mergerTreeImportSussingBadValueTest
  integer                        :: mergerTreeImportSussingMassOption

  ! Bad value detection limits.
  integer         , parameter    :: sussingBadValueLessThan                  =-1
  integer         , parameter    :: sussingBadValueGreaterThan               =+1

  ! Mass options.
  integer         , parameter    :: sussingMassDefault                       = 1
  integer         , parameter    :: sussingMassFoF                           = 2
  integer         , parameter    :: sussingMass200Mean                       = 3
  integer         , parameter    :: sussingMass200Crit                       = 4
  integer         , parameter    :: sussingMassTopHat                        = 5
 
contains

  subroutine sussingInitialize(self)
    !% Initializor for the ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree importer.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    type (varying_string           )                :: mergerTreeImportSussingBadValueTestText, mergerTreeImportSussingMassOptionText

    if (.not.sussingInitialized) then
       !$omp critical (mergerTreeImporterSussingInitialize)
       if (.not.sussingInitialized) then
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingMismatchIsFatal</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>true</defaultValue>
          !@   <description>
          !@     Specifies whether mismatches in cosmological parameter values between \glc\ and ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree files should be considered fatal.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingMismatchIsFatal',mergerTreeImportSussingMismatchIsFatal,defaultValue=.true.)
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingNonTreeNodeIsFatal</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>true</defaultValue>
          !@   <description>
          !@     Specifies whether nodes in snapshot files but not in the merger tree file should be considered fatal when importing from the ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013}.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingNonTreeNodeIsFatal',mergerTreeImportSussingNonTreeNodeIsFatal,defaultValue=.true.)
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingSubvolumeCount</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>1</defaultValue>
          !@   <description>
          !@    Specifies the number of subvolumes \emph{along each axis} into which a ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree files should be split for processing through \glc.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingSubvolumeCount',mergerTreeImportSussingSubvolumeCount,defaultValue=1)
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingSubvolumeBuffer</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>$0.0$</defaultValue>
          !@   <description>
          !@     Specifies the buffer region (in units of Mpc$/h$ to follow the format convention) around subvolumes of a ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree file which should be read in to ensure that no halos are missed from trees.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingSubvolumeBuffer',mergerTreeImportSussingSubvolumeBuffer,defaultValue=0.0d0)
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingSubvolumeIndex</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>[0,0,0]</defaultValue>
          !@   <description>
          !@     Specifies the index (in each dimension) of the subvolume of a ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree file to process. Indices range from 0 to {\normalfont \ttfamily [mergerTreeImportSussingSubvolumeCount]}$-1$.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingSubvolumeIndex',mergerTreeImportSussingSubvolumeIndex,defaultValue=[0,0,0])
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingConvertToBinary</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>true</defaultValue>
          !@   <description>
          !@     Specifies whether halo and tree files in the ``Sussing'' format should be converted to binary the first time they are read and stored to file. This allows rapid re-reading in future.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingConvertToBinary',mergerTreeImportSussingConvertToBinary,defaultValue=.true.)
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingBinaryFormatOld</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>true</defaultValue>
          !@   <description>
          !@     Specifies whether the old binary format is to be used (for reading only).
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingBinaryFormatOld',mergerTreeImportSussingBinaryFormatOld,defaultValue=.false.)
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingBadValue</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>-0.5</defaultValue>
          !@   <description>
          !@     Use for bad value detection in ``Sussing'' merger trees. Values for scale radius and halo spin which exceed this threshold are assumed to be bad.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingBadValue',mergerTreeImportSussingBadValue,defaultValue=-0.5d0)
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingBadValueTest</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>lessThan</defaultValue>
          !@   <description>
          !@     Use for bad value detection in ``Sussing'' merger trees. Values which exceed the threshold in ths specified direction are assumed to be bad.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingBadValueTest',mergerTreeImportSussingBadValueTestText,defaultValue="lessThan")
          select case (char(mergerTreeImportSussingBadValueTestText))
          case ("lessThan"   )
             mergerTreeImportSussingBadValueTest=sussingBadValueLessThan
          case ("greaterThan")
             mergerTreeImportSussingBadValueTest=sussingBadValueGreaterThan
          case default
             call Galacticus_Error_Report('sussingDefaultConstructor','[mergerTreeImportSussingBadValueTest must be either "lessThan" or "greaterThan"]')
          end select
          ! Check for a forest file.
          mergerTreeImportSussingUseForestFile=Input_Parameter_Is_Present('mergerTreeImportSussingForestFile')
          if (mergerTreeImportSussingUseForestFile) then
             !@ <inputParameter>
             !@   <name>mergerTreeImportSussingForestFile</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Name of file containing data on number of halos in each forest.
             !@   </description>
             !@   <type>string</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeImportSussingForestFile',mergerTreeImportSussingForestFile)
             !@ <inputParameter>
             !@   <name>mergerTreeImportSussingForestFirst</name>
             !@   <attachedTo>module</attachedTo>
             !@   <defaultValue>1</defaultValue>
             !@   <description>
             !@     Index of first forest to include.
             !@   </description>
             !@   <type>integer</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeImportSussingForestFirst',mergerTreeImportSussingForestFirst,defaultValue=1)
             !@ <inputParameter>
             !@   <name>mergerTreeImportSussingForestLast</name>
             !@   <attachedTo>module</attachedTo>
             !@   <defaultValue>-1</defaultValue>
             !@   <description>
             !@     Index of last forest to include.
             !@   </description>
             !@   <type>integer</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeImportSussingForestLast',mergerTreeImportSussingForestLast,defaultValue=-1)
             !@ <inputParameter>
             !@   <name>mergerTreeImportSussingForestReverseSnapshotOrder</name>
             !@   <attachedTo>module</attachedTo>
             !@   <defaultValue>false</defaultValue>
             !@   <description>
             !@     If true, the order of forest snapshots will be reversed after being read. This may be necessary to cause them to match the order of snapshot files.
             !@   </description>
             !@   <type>integer</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeImportSussingForestReverseSnapshotOrder',mergerTreeImportSussingForestReverseSnapshotOrder,defaultValue=.false.)
          end if
          !@ <inputParameter>
          !@   <name>mergerTreeImportSussingMassOption</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>default</defaultValue>
          !@   <description>
          !@     Mass option for Sussing merger trees.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportSussingMassOption',mergerTreeImportSussingMassOptionText,defaultValue="default")
          select case (char(mergerTreeImportSussingMassOptionText))
          case ("default")
             mergerTreeImportSussingMassOption=sussingMassDefault
          case ("FoF"    )
             mergerTreeImportSussingMassOption=sussingMassFoF
          case ("200Mean")
             mergerTreeImportSussingMassOption=sussingMass200Mean
          case ("200Crit")
             mergerTreeImportSussingMassOption=sussingMass200Crit
          case ("TopHat" )
             mergerTreeImportSussingMassOption=sussingMassTopHat
          case default
             call Galacticus_Error_Report('sussingDefaultConstructor','unrecognized mass option')
          end select
          sussingInitialized=.true.
       end if
       !$omp end critical (mergerTreeImporterSussingInitialize)
    end if
    self%treeIndicesRead         =.false.
    self%fatalMismatches         =mergerTreeImportSussingMismatchIsFatal
    self%subvolumeCount          =mergerTreeImportSussingSubvolumeCount
    self%subvolumeBuffer         =mergerTreeImportSussingSubvolumeBuffer
    self%subvolumeIndex          =mergerTreeImportSussingSubvolumeIndex
    self%convertToBinary         =mergerTreeImportSussingConvertToBinary
    self%binaryFormatOld         =mergerTreeImportSussingBinaryFormatOld
    self%badValue                =mergerTreeImportSussingBadValue
    self%badTest                 =mergerTreeImportSussingBadValueTest
    self%spinsAvailableValue     =.true.
    self%scaleRadiiAvailableValue=.true.
    return
  end subroutine sussingInitialize

  subroutine sussingDestroy(self)
    !% Destructor for the {\normalfont \ttfamily sussing} format merger tree importer class.
    use Memory_Management
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    if (allocated(self%treeIndices)) call Dealloc_Array(self%treeIndices)
    if (allocated(self%treeSizes  )) call Dealloc_Array(self%treeSizes  )
    return
  end subroutine sussingDestroy

  subroutine sussingClose(self)
    !% Close a {\normalfont \ttfamily sussing} format merger tree file.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    
    return
  end subroutine sussingClose

  integer function sussingTreesHaveSubhalos(self)
    !% Return a Boolean integer specifying whether or not the trees have subhalos.
    use Numerical_Constants_Boolean
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    sussingTreesHaveSubhalos=booleanTrue
    return
  end function sussingTreesHaveSubhalos

  logical function sussingMassesIncludeSubhalos(self)
    !% Return a Boolean specifying whether or not the halo masses include the contribution from subhalos.
    use Galacticus_Error
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self

    sussingMassesIncludeSubhalos=.true.
    return
  end function sussingMassesIncludeSubhalos

  logical function sussingAngularMomentaIncludeSubhalos(self)
    !% Return a Boolean specifying whether or not the halo angular momenta include the contribution from subhalos.
    use Galacticus_Error
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self

    sussingAngularMomentaIncludeSubhalos=.true.
    return
  end function sussingAngularMomentaIncludeSubhalos

  integer function sussingTreesAreSelfContained(self)
    !% Return a Boolean integer specifying whether or not the trees are self-contained.
    use Numerical_Constants_Boolean
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self

    sussingTreesAreSelfContained=booleanTrue
    return
  end function sussingTreesAreSelfContained

  integer function sussingVelocitiesIncludeHubbleFlow(self)
    !% Return a Boolean integer specifying whether or not velocities include the Hubble flow.
    use Numerical_Constants_Boolean
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self

    sussingVelocitiesIncludeHubbleFlow=booleanFalse
    return
  end function sussingVelocitiesIncludeHubbleFlow

  integer function sussingPositionsArePeriodic(self)
    !% Return a Boolean integer specifying whether or not positions are periodic.
    use Numerical_Constants_Boolean
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self

    sussingPositionsArePeriodic=booleanTrue
    return
  end function sussingPositionsArePeriodic

  double precision function sussingCubeLength(self,time,status)
    !% Return the length of the simulation cube.
    use Numerical_Constants_Boolean
    use Numerical_Constants_Astronomical
    implicit none
    class           (mergerTreeImporterSussing), intent(inout)           :: self
    double precision                           , intent(in   )           :: time
    integer                                    , intent(  out), optional :: status
 
    sussingCubeLength=importerUnitConvert(self%boxLength,time,self%boxLengthUnits,megaParsec)
    if (present(status)) status=booleanTrue
    return
  end function sussingCubeLength

  integer(kind=c_size_t) function sussingTreeCount(self)
    !% Return a count of the number of trees available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    call sussingTreeIndicesRead(self)
    sussingTreeCount=self%treesCount
    return
  end function sussingTreeCount

  integer(kind=c_size_t) function sussingTreeIndex(self,i)
    !% Return the index of the $i^{\mathrm th}$ tree.
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call sussingTreeIndicesRead(self)
    sussingTreeIndex=self%treeIndices(i)
    return
  end function sussingTreeIndex

  function sussingNodeCount(self,i)
    !% Return a count of the number of nodes in the $i^{\mathrm th}$ tree.
    implicit none
    integer(c_size_t                 )                :: sussingNodeCount
    class  (mergerTreeImporterSussing), intent(inout) :: self
    integer                           , intent(in   ) :: i

    call sussingTreeIndicesRead(self)
    sussingNodeCount=self%treeSizes(i)
    return
  end function sussingNodeCount

  subroutine sussingTreeIndicesRead(self)
    !% Read the tree indices.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Display
    use Galacticus_Error
    use Kind_Numbers
    use String_Handling
    use Sort
    use Memory_Management
    use Arrays_Search
    use Array_Utilities
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    implicit none
    class    (mergerTreeImporterSussing), intent(inout)                 :: self
    integer  (kind=kind_int8           ), allocatable  , dimension(:  ) :: nodeSelfIndices            , nodeTreeIndices
    integer  (kind=c_size_t            ), allocatable  , dimension(:  ) :: nodeIndexRanks             , nodeDescendentLocations
    logical                             , allocatable  , dimension(:  ) :: nodeIncomplete
    integer                             , parameter                     :: stepCountMaximum   =1000000
    integer                                                             :: descendentStepCount        , hostStepCount
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
         &         nodeDescendentLocations, &
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
    call Galacticus_Display_Indent('Resolving hosting loops',verbosityWorking)
    do i=1,nodeCountTrees
       if (self%nodes(i)%hostIndex /= self%nodes(i)%nodeIndex) then
          l=Search_Indexed(nodeSelfIndices,nodeIndexRanks,self%nodes(i)%hostIndex)
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
    call Galacticus_Display_Unindent('done',verbosityWorking)
    ! Search for deep hosting hierarchies and reset to single level hierarchies.
    call Galacticus_Display_Indent('Resolving deep hosting',verbosityWorking)
    do i=1,nodeCountTrees
       if (self%nodes(i)%hostIndex /= self%nodes(i)%nodeIndex) then
          ! Find the host.
          l=Search_Indexed(nodeSelfIndices,nodeIndexRanks,self%nodes(i)%hostIndex)
          ! Detect missing host.          
          if (l < 1 .or. l > nodeCountTrees) cycle
          l=nodeIndexRanks(l)
          if (self%nodes(l)%nodeIndex /= self%nodes(i)%hostIndex) cycle
          ! Detect hosted host.
          hostStepCount=0
          do while (self%nodes(l)%hostIndex /= self%nodes(l)%nodeIndex)
             ! Find the host.
             j=Search_Indexed(nodeSelfIndices,nodeIndexRanks,self%nodes(l)%hostIndex)
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
                call Galacticus_Error_Report('sussingTreeIndicesRead',message)
             end if
             ! Move to the new host.
             l=j
          end do
          self%nodes(i)%hostIndex=self%nodes(l)%nodeIndex
       end if
    end do
    call Galacticus_Display_Unindent('done',verbosityWorking)
    ! Assign tree indices.
    call Galacticus_Display_Message('Assigning tree indices',verbosityWorking)
    if (.not.treeIndicesAssigned) then
       nodeTreeIndices=nodeSelfIndices
       do i=1,nodeCountTrees
          call Galacticus_Display_Counter(int(100.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityWorking)
          l=i
          hostStepCount=0
          do while (l /= -1 .and. self%nodes(l)%hostIndex /= self%nodes(l)%nodeIndex) 
             ! Find the host halo.
             k=Search_Indexed(nodeSelfIndices,nodeIndexRanks,self%nodes(l)%hostIndex)
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
                call Galacticus_Error_Report('sussingTreeIndicesRead',message)
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
                call Galacticus_Error_Report('sussingTreeIndicesRead',message)
             end if
             l=nodeIndexRanks(k)
          end do
          descendentStepCount=0
          do while (l /= -1)
             nodeTreeIndices(i)=nodeTreeIndices(l)
             k=nodeDescendentLocations(l)
             ! Check for missing descendents.
             if (k < 0 .or. self%nodes(l)%descendentIndex /= self%nodes(k)%nodeIndex) exit
             ! Perform sanity checks.
             if (k >= 0 .and. self%nodes(k)%nodeTime <= self%nodes(l)%nodeTime) then
                message='descendent exists before progenitor node'
                message=message//char(10)//' progenitor node index: '//self%nodes(l)%nodeIndex
                message=message//char(10)//' descendent node index: '//self%nodes(k)%nodeIndex
                write (label,'(f7.4)') self%nodes(l)%nodeTime
                message=message//char(10)//'  progenitor node time: '//label
                write (label,'(f7.4)') self%nodes(k)%nodeTime
                message=message//char(10)//'  descendent node time: '//label
                call Galacticus_Error_Report('sussingTreeIndicesRead',message)
             end if
             descendentStepCount=descendentStepCount+1
             if (descendentStepCount > stepCountMaximum) then
                message='infinite (or at least, very large) loop detect in halo descent'
                message=message//char(10)//' progenitor node index: '//self%nodes(l)%nodeIndex
                message=message//char(10)//' descendent node index: '//self%nodes(k)%nodeIndex
                write (label,'(f7.4)') self%nodes(l)%nodeTime
                message=message//char(10)//'  progenitor node time: '//label
                write (label,'(f7.4)') self%nodes(k)%nodeTime
                message=message//char(10)//'  descendent node time: '//label
                call Galacticus_Error_Report('sussingTreeIndicesRead',message)
             end if
             l=k
             hostStepCount=0
             do while (l /= -1 .and. self%nodes(l)%hostIndex /= self%nodes(l)%nodeIndex) 
                k=Search_Indexed(nodeSelfIndices,nodeIndexRanks,self%nodes(l)%hostIndex)
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
                   call Galacticus_Error_Report('sussingTreeIndicesRead',message)
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
                   call Galacticus_Error_Report('sussingTreeIndicesRead',message)
                end if
                l=nodeIndexRanks(k)
             end do
          end do
       end do
    end if
    call Galacticus_Display_Counter_Clear(verbosityWorking)
    ! Check for nodes jumping between trees and join any such trees.
    if (branchJumpCheckRequired) then
       call Galacticus_Display_Message('Checking for branch jumps',verbosityWorking)
       do i=1,nodeCountTrees
          call Galacticus_Display_Counter(int(100.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityWorking)
          l=nodeDescendentLocations(i)
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
          l=nodeDescendentLocations(i)
          if (l /= -1) then
             if (nodeTreeIndices(i) /= nodeTreeIndices(l))                                              &
                  & call Galacticus_Error_Report('sussingTreeIndicesRead','failed to cross-link trees')
          end if
       end do
       call Galacticus_Display_Counter_Clear(verbosityWorking)
    end if
    ! Generate an index into nodes sorted by tree index.
    self%treeIndexRanks=Sort_Index_Do(nodeTreeIndices)
    ! Identify trees which contain incomplete nodes.
    call Galacticus_Display_Message('Checking for incomplete trees',verbosityWorking)
    i               =0
    iStart          =0
    treeIndexCurrent=-1
    do while (i < nodeCountTrees)
       i     =i+1
       ! Check if we started a new tree, and record the index at which it began.
       if (treeIndexCurrent /= nodeTreeIndices(self%treeIndexRanks(i))) iStart=i
       treeIndexCurrent=nodeTreeIndices(self%treeIndexRanks(i))
       if (nodeIncomplete(self%treeIndexRanks(i))) then
          if (Galacticus_Verbosity_Level() >= verbosityWorking) then
             j=self%treeIndexRanks(i)
             message='Marking tree '
             message=message//treeIndexCurrent//' as incomplete due to node '//nodeSelfIndices(j)//' at position:'
             write (label,'(e12.6)') self%nodes(j)%position(1)
             message=message//char(10)//'  x: '//trim(label)
             write (label,'(e12.6)') self%nodes(j)%position(2)
             message=message//char(10)//'  y: '//trim(label)
             write (label,'(e12.6)') self%nodes(j)%position(3)
             message=message//char(10)//'  z: '//trim(label)
             call Galacticus_Display_Message(message,verbosityWorking)
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
       call Galacticus_Display_Counter(int(50.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityWorking)
    end do
    ! Check for zero mass halos.
    call Galacticus_Display_Message('Checking for zero mass halos',verbosityWorking)
    do i=1,nodeCountTrees
       call Galacticus_Display_Counter(int(100.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityWorking)
       if (self%nodes(i)%nodeMass == 0.0d0) then
          k=i
          do while (self%nodes(i)%nodeMass == 0.0d0)
             if (self%nodes(k)%nodeMass > 0.0d0) then
                self%nodes(i)%nodeMass   =self%nodes(k)%nodeMass
                self%nodes(i)%scaleRadius=self%nodes(k)%scaleRadius
             end if
             k=nodeDescendentLocations(k)
             if (k < 0 .or. self%nodes(i)%descendentIndex /= self%nodes(k)%nodeIndex) then
                message='zero mass halo ['
                message=message//self%nodes(i)%nodeIndex//'] has no non-zero mass descendents'
                call Galacticus_Error_Report('sussingTreeIndicesRead',message)
             end if
          end do
       end if
    end do
    call Galacticus_Display_Counter_Clear(verbosityWorking)
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
                call Galacticus_Error_Report('sussingTreeIndicesRead',message)
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
       call Galacticus_Display_Counter(int(50.0d0+50.0d0*dble(i)/dble(nodeCountTrees)),.false.,verbosityWorking)
    end do
    call Galacticus_Display_Counter_Clear(verbosityWorking)
    ! Generate an index into nodes sorted by tree index.
    self%treeIndexRanks=Sort_Index_Do(nodeTreeIndices)
    ! Create a list of tree indices, sizes, and start locations.
    call Galacticus_Display_Message('Generating tree list',verbosityWorking)
    self%treesCount=0
    treeIndexPrevious=-1
    do i=1,nodeCountTrees
       call Galacticus_Display_Counter(int(100.0d0*dble(i)/dble(nodeCountTrees)),i==1,verbosityWorking)
       if (nodeTreeIndices(self%treeIndexRanks(i)) /= treeIndexPrevious) then
          treeIndexPrevious=nodeTreeIndices(self%treeIndexRanks(i))
          self%treesCount=self%treesCount+1
       end if
    end do
    call Galacticus_Display_Counter_Clear(verbosityWorking)
    message='Found '
    message=message//self%treesCount//' trees'
    call Galacticus_Display_Message(message,verbosityWorking)
    call Alloc_Array(self%treeIndices,[self%treesCount])
    call Alloc_Array(self%treeSizes  ,[self%treesCount])
    call Alloc_Array(self%treeBegins ,[self%treesCount])
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
    call Galacticus_Display_Counter_Clear(verbosityWorking)
    ! Do unit conversion.
    self   %nodes%nodeMass             =importerUnitConvert(self%nodes%nodeMass             ,self%nodes%nodeTime,massUnits    ,massSolar )
    self   %nodes%scaleRadius          =importerUnitConvert(self%nodes%scaleRadius          ,self%nodes%nodeTime,lengthUnits  ,megaParsec)
    self   %nodes%velocityMaximum      =importerUnitConvert(self%nodes%velocityMaximum      ,self%nodes%nodeTime,velocityUnits,kilo      )
    self   %nodes%velocityDispersion   =importerUnitConvert(self%nodes%velocityDispersion   ,self%nodes%nodeTime,velocityUnits,kilo      )
    do i=1,3
       self%nodes%position          (i)=importerUnitConvert(self%nodes%position          (i),self%nodes%nodeTime,lengthUnits  ,megaParsec)
       self%nodes%velocity          (i)=importerUnitConvert(self%nodes%velocity          (i),self%nodes%nodeTime,velocityUnits,kilo      )
    end do
    ! Destroy temporary workspace.
    call Dealloc_Array(nodeSelfIndices        )
    call Dealloc_Array(nodeTreeIndices        )
    call Dealloc_Array(nodeDescendentLocations)
    ! Write completion message.
    call Galacticus_Display_Unindent('done',verbosityWorking)
   return
  end subroutine sussingTreeIndicesRead

  double precision function sussingTreeWeight(self,i)
    !% Return the weight to assign to trees.
    use Numerical_Constants_Boolean
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    use Cosmology_Functions
    implicit none
    class           (mergerTreeImporterSussing), intent(inout) :: self
    integer                                    , intent(in   ) :: i
    class           (cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_

    ! Compute the inverse of the cube volume.
    cosmologyFunctions_ => cosmologyFunctions()
    sussingTreeWeight   =  1.0d0/self%cubeLength(cosmologyFunctions_%cosmicTime(1.0d0))**3
    return
  end function sussingTreeWeight

  logical function sussingPositionsAvailable(self,positions,velocities)
    !% Return true if positions and/or velocities are available.
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    logical                              , intent(in   ) :: positions, velocities

    sussingPositionsAvailable=.true.
    return
  end function sussingPositionsAvailable

  logical function sussingScaleRadiiAvailable(self)
    !% Return true if scale radii are available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    call sussingTreeIndicesRead(self)
    sussingScaleRadiiAvailable=self%scaleRadiiAvailableValue
    return
  end function sussingScaleRadiiAvailable

  logical function sussingParticleCountAvailable(self)
    !% Return true if particle counts are available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    sussingParticleCountAvailable=.false.
    return
  end function sussingParticleCountAvailable

  logical function sussingVelocityMaximumAvailable(self)
    !% Return true if halo rotation curve velocity maxima are available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    sussingVelocityMaximumAvailable=.true.
    return
  end function sussingVelocityMaximumAvailable

  logical function sussingVelocityDispersionAvailable(self)
    !% Return true if halo velocity dispersions are available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    sussingVelocityDispersionAvailable=.true.
    return
  end function sussingVelocityDispersionAvailable

  logical function sussingAngularMomentaAvailable(self)
    !% Return true if angular momenta are available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    
    sussingAngularMomentaAvailable=.false.
    return
  end function sussingAngularMomentaAvailable

  logical function sussingAngularMomenta3DAvailable(self)
    !% Return true if angular momenta vectors are available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    sussingAngularMomenta3DAvailable=.false.
    return
  end function sussingAngularMomenta3DAvailable

  logical function sussingSpinAvailable(self)
    !% Return true if spins are available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    
    call sussingTreeIndicesRead(self)
    sussingSpinAvailable=self%spinsAvailableValue
    return
  end function sussingSpinAvailable

  logical function sussingSpin3DAvailable(self)
    !% Return true if spins vectors are available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    call sussingTreeIndicesRead(self)
    sussingSpin3DAvailable=self%spinsAvailableValue
    return
  end function sussingSpin3DAvailable

  subroutine sussingSubhaloTrace(self,node,time,position,velocity)
    !% Returns a trace of subhalo position/velocity.
    use Galacticus_Error
    implicit none
    class           (mergerTreeImporterSussing), intent(inout)                 :: self
    class           (nodeData                 ), intent(in   )                 :: node
    double precision                           , intent(  out), dimension(:  ) :: time
    double precision                           , intent(  out), dimension(:,:) :: position, velocity

    call Galacticus_Error_Report('sussingSubhaloTrace','subhalo traces are not available')
    return
  end subroutine sussingSubhaloTrace

  function sussingSubhaloTraceCount(self,node)
    !% Returns the length of a subhalo trace.
    implicit none
    integer(c_size_t                 )                :: sussingSubhaloTraceCount
    class  (mergerTreeImporterSussing), intent(inout) :: self
    class  (nodeData                 ), intent(in   ) :: node

    ! No particle data is available.
    sussingSubhaloTraceCount=0
    return
  end function sussingSubhaloTraceCount

  subroutine sussingImport(self,i,nodes,requireScaleRadii,requireAngularMomenta,requireAngularMomenta3D,requireSpin,requireSpin3D,requirePositions,requireParticleCounts,requireVelocityMaxima,requireVelocityDispersions)
    !% Import the $i^{\mathrm th}$ merger tree.
    use Memory_Management
    implicit none
    class           (mergerTreeImporterSussing), intent(inout)                              :: self
    integer                                    , intent(in   )                              :: i
    class           (nodeData                 ), intent(  out), allocatable, dimension(:  ) :: nodes
    logical                                    , intent(in   ), optional                    :: requireScaleRadii         , requireAngularMomenta, &
         &                                                                                     requireAngularMomenta3D   , requirePositions     , &
         &                                                                                     requireParticleCounts     , requireVelocityMaxima, &
         &                                                                                     requireVelocityDispersions, requireSpin         , &
         &                                                                                     requireSpin3D
    integer         (c_size_t                 )                                             :: j

    ! Allocate the nodes array.
    allocate(nodeData :: nodes(self%treeSizes(i)))
    !# <workaround type="gfortran" PR="65889" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=65889">
    select type (nodes)
    type is (nodeData)
       call Memory_Usage_Record(sizeof(nodes))
    end select
    !# </workaround>
    ! Copy data to nodes.
    select type (nodes)
    type is (nodeData)
       do j=1,self%treeSizes(i)
          nodes(j)=self%nodes(self%treeIndexRanks(self%treeBegins(i)+j-1))
       end do
    end select
    return
  end subroutine sussingImport

  logical function sussingInSubvolume(self,x,y,z,buffered)
    !% Determine if a point lies within a subvolume of the simulation box (possibly with some buffering).
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
    !% Determine if a point lies within the 1-D range of a subvolume of the simulation box (possibly with some buffering).
    use Numerical_Constants_Astronomical
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
    !% Determine the separation between two points in a periodic cube.
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
    !% Determine if a value in a ``Sussing'' merger tree file is bad
    implicit none
    class           (mergerTreeImporterSussing), intent(inout) :: self
    double precision                           , intent(in   ) :: x
  
    select case (self%badTest)
    case (sussingBadValueLessThan   )
       sussingValueIsBad=(x < self%badValue)
    case (sussingBadValueGreaterThan)
       sussingValueIsBad=(x > self%badValue)
    end select
    return
  end function sussingValueIsBad
