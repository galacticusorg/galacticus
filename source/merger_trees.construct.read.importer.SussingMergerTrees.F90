!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of the merger tree importer class for \glc\ format merger tree files.

  !# <mergerTreeImporter name="mergerTreeImporterSussing" description="Importer for ``Sussing Merger Trees'' format merger tree files \citep{srisawat_sussing_2013}." />
  use IO_HDF5
  use Stateful_Types
  use ISO_Varying_String

  type, public, extends(nodeData) :: nodeDataSussing
     !% Extension of the {\tt nodeData} class for ``Sussing Merger Trees'' format merger trees \citep{srisawat_sussing_2013}.
  end type nodeDataSussing

  type, extends(mergerTreeImporterClass) :: mergerTreeImporterSussing
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
     integer                         , allocatable, dimension(:) :: treeSizes               , treeBegins
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
     !@     <method>inSubvolume</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\doublezero\ x\argin, \doublezero\ y\argin, \doublezero\ z\argin, \logicalzero\ [buffered]\argin</arguments>
     !@     <description>Return true if the given {\tt x,y,z} position lies within the current subvolume (plus the buffer region if {\tt buffered} is true.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>inSubvolume1D</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\doublezero\ x\argin, \intzero\ iSubvolume\argin, \logicalzero\ [buffered]\argin</arguments>
     !@     <description>Return true if the given {\tt x} position lies within the {\tt iSubvolume}$^{\rm th}$ subvolume (plus the buffer region if {\tt buffered} is true.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>valueIsBad</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\doublezero\ x\argin</arguments>
     !@     <description>Return true if the given {\tt x} value is bad.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     !# <workaround type="gfortran" PR="58471 58470" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58471 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58470">
     !# final     :: sussingDestructor
     !# </workaround>
     procedure :: open                          => sussingOpen
     procedure :: close                         => sussingClose
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
     !% Constructors for the \glc\ format merger tree importer class.
     module procedure sussingDefaultConstructor
  end interface mergerTreeImporterSussing

  ! Record of implementation initialization state.
  logical                                        :: sussingInitialized                       =.false.

  ! Default settings.
  logical                                        :: mergerTreeImportSussingMismatchIsFatal
  logical                                        :: mergerTreeImportSussingNonTreeNodeIsFatal
  integer                                        :: mergerTreeImportSussingSubvolumeCount
  double precision                               :: mergerTreeImportSussingSubvolumeBuffer
  integer                         , dimension(3) :: mergerTreeImportSussingSubvolumeIndex
  logical                                        :: mergerTreeImportSussingConvertToBinary
  logical                                        :: mergerTreeImportSussingBinaryFormatOld
  double precision                               :: mergerTreeImportSussingBadValue
  integer                                        :: mergerTreeImportSussingBadValueTest
  logical                                        :: mergerTreeImportSussingUseForestFile
  type            (varying_string)               :: mergerTreeImportSussingForestFile
  integer                                        :: mergerTreeImportSussingForestFirst
  integer                                        :: mergerTreeImportSussingForestLast
  integer                                        :: mergerTreeImportSussingMassOption

  ! Bad value detection limits.
  integer                         , parameter    :: sussingBadValueLessThan                  =-1
  integer                         , parameter    :: sussingBadValueGreaterThan               =+1
 
  ! File format identifiers.
  integer                         , parameter    :: sussingHaloFormatOld                     = 1
  integer                         , parameter    :: sussingHaloFormatNew                     = 2
  integer                         , parameter    :: sussingHaloFormatAll                     = 3
 
  ! Mass options.
  integer                         , parameter    :: sussingMassDefault                       = 1
  integer                         , parameter    :: sussingMassFoF                           = 2
  integer                         , parameter    :: sussingMass200Mean                       = 3
  integer                         , parameter    :: sussingMass200Crit                       = 4
  integer                         , parameter    :: sussingMassTopHat                        = 5
 
contains

  function sussingDefaultConstructor()
    !% Default constructor for the ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree importer.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(mergerTreeImporterSussing), target :: sussingDefaultConstructor
    type(varying_string           )         :: mergerTreeImportSussingBadValueTestText, mergerTreeImportSussingMassOptionText

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
          !@     Specifies the index (in each dimension) of the subvolume of a ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree file to process. Indices range from 0 to {\tt [mergerTreeImportSussingSubvolumeCount]}$-1$.
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
          end if
          sussingInitialized=.true.
       end if
       !$omp end critical (mergerTreeImporterSussingInitialize)
    end if
    sussingDefaultConstructor%treeIndicesRead         =.false.
    sussingDefaultConstructor%fatalMismatches         =mergerTreeImportSussingMismatchIsFatal
    sussingDefaultConstructor%subvolumeCount          =mergerTreeImportSussingSubvolumeCount
    sussingDefaultConstructor%subvolumeBuffer         =mergerTreeImportSussingSubvolumeBuffer
    sussingDefaultConstructor%subvolumeIndex          =mergerTreeImportSussingSubvolumeIndex
    sussingDefaultConstructor%convertToBinary         =mergerTreeImportSussingConvertToBinary
    sussingDefaultConstructor%binaryFormatOld         =mergerTreeImportSussingBinaryFormatOld
    sussingDefaultConstructor%badValue                =mergerTreeImportSussingBadValue
    sussingDefaultConstructor%badTest                 =mergerTreeImportSussingBadValueTest
    sussingDefaultConstructor%spinsAvailableValue     =.true.
    sussingDefaultConstructor%scaleRadiiAvailableValue=.true.
   return
  end function sussingDefaultConstructor

  subroutine sussingDestructor(self)
    !% Destructor for the \glc\ format merger tree importer class.
    implicit none
    type(mergerTreeImporterSussing), intent(inout) :: self

    if (allocated(self%treeIndices)) call Dealloc_Array(self%treeIndices)
    if (allocated(self%treeSizes  )) call Dealloc_Array(self%treeSizes  )
    return
  end subroutine sussingDestructor

  subroutine sussingOpen(self,fileName)
    !% Validate a \glc\ format merger tree file.
    use Numerical_Comparison
    use Numerical_Constants_Astronomical
    use Galacticus_Display
    use Galacticus_Error
    use Cosmology_Parameters
    use Cosmology_Functions
    use Regular_Expressions
    use String_Handling
    use File_Utilities
    use Memory_Management
    implicit none
    class           (mergerTreeImporterSussing   ), intent(inout) :: self
    type            (varying_string              ), intent(in   ) :: fileName
    class           (cosmologyParametersClass    ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    type            (varying_string              )                :: message                 , baseDirectory        , &
         &                                                           simulationDefinitionFile, snapshotTimesFile
    character       (len=1024                    )                :: line                    , parameterName        , &
         &                                                           parameterValue
    character       (len=14                      )                :: valueString             , unitString
    type            (regEx                       )                :: parameterRegEx
    type            (importerUnits               )                :: boxLengthUnits
    integer                                                       :: simulationDefintion     , fileUnit             , &
         &                                                           ioStat                  , snapshotFileCount    , &
         &                                                           i                       , snapshotNumber
    double precision                                              :: localLittleH0           , localOmegaMatter     , &
         &                                                           localOmegaDE            , cosmologicalParameter, &
         &                                                           expansionFactor         , redshift             , &
         &                                                           timeNormalized          , time

    ! Get the default cosmology.
    cosmologyParameters_ => cosmologyParameters()
    cosmologyFunctions_  => cosmologyFunctions ()
    ! Get cosmological parameters. We do this in advance to avoid HDF5 thread conflicts.
    localLittleH0   =cosmologyParameters_%HubbleConstant (unitsLittleH)
    localOmegaMatter=cosmologyParameters_%OmegaMatter    (            )
    localOmegaDE    =cosmologyParameters_%OmegaDarkEnergy(            )
    ! Determine the base location for files.
    if (index(fileName,'/',back=.true.) == 0) then
       baseDirectory=""
    else
       baseDirectory=extract(fileName,1,index(fileName,'/',back=.true.))
    end if        
    ! Read the file definition file.
    snapshotFileCount=Count_Lines_in_File(fileName)-3
    allocate(self%snapshotFileName(snapshotFileCount))
    open(newUnit=fileUnit,file=char(fileName),status='old',form='formatted')
    read (fileUnit,'(a)') line
    simulationDefinitionFile=baseDirectory//trim(line)
    read (fileUnit,'(a)') line
    self%mergerTreeFile     =baseDirectory//trim(line)
    read (fileUnit,'(a)') line
    snapshotTimesFile       =baseDirectory//trim(line)
    do i=1,snapshotFileCount
       read (fileUnit,'(a)') line
       self%snapshotFileName(i)=baseDirectory//trim(line)
    end do
    close(fileUnit)
    ! Read the snapshots times file.
    call Alloc_Array(self%snapshotTimes,[snapshotFileCount])
    open(newUnit=fileUnit,file=char(snapshotTimesFile),status='old',form='formatted',ioStat=ioStat)
    read (fileUnit,*)
    do i=1,snapshotFileCount
       read (fileUnit,*) snapshotNumber,expansionFactor,redshift,timeNormalized,time
       self%snapshotTimes(i)=                                            &
            & cosmologyFunctions_ %cosmicTime                 (          &
            &  cosmologyFunctions_%expansionFactorFromRedshift (         &
            &                                                   redshift &
            &                                                  )         &
            &                                                 )
    end do
    close(fileUnit)
    ! Read the simulation definition file.
    parameterRegEx=regEx('^([a-zA-Z0-9]+)\s*=\s*(.*)')
    open(newUnit=fileUnit,file=char(simulationDefinitionFile),status='old',form='formatted',ioStat=ioStat)
    do while (ioStat == 0)
       read (fileUnit,'(a)',ioStat=ioStat) line
       if (parameterRegEx%matches(line)) then
          parameterName =String_Strip(line(1:index(line,'=')-1                          ))
          parameterValue=String_Strip(line(  index(line,'=')+1:len(line)-index(line,'=')))
          select case (trim(parameterName))
          case ('h')
             read (parameterValue,*) cosmologicalParameter
             if (Values_Differ(cosmologicalParameter,localLittleH0,absTol=0.00001d0)) then
                message='H_0 in merger tree file ['
                write (valueString,'(e14.8)') cosmologicalParameter
                message=message//trim(valueString)//'] differs from the internal value ['
                write (valueString,'(e14.8)') localLittleH0
                message=message//trim(valueString)//']'
                if (self%fatalMismatches) then
                   call Galacticus_Error_Report('sussingOpen',message)
                else
                   call Galacticus_Display_Message(message,verbosityWarn)
                end if
             end if
          case ('Omega0' )
             read (parameterValue,*) cosmologicalParameter
             if (Values_Differ(cosmologicalParameter,localOmegaMatter,absTol=0.0001d0)) then
                message='Omega_Matter in merger tree file ['
                write (valueString,'(e14.8)') cosmologicalParameter
                message=message//trim(valueString)//'] differs from the internal value ['
                write (valueString,'(e14.8)') localOmegaMatter
                message=message//trim(valueString)//']'
                if (self%fatalMismatches) then
                   call Galacticus_Error_Report('sussingOpen',message)
                else
                   call Galacticus_Display_Message(message,verbosityWarn)
                end if
             end if
          case ('OmegaDE')
             read (parameterValue,*) cosmologicalParameter
             if (Values_Differ(cosmologicalParameter,localOmegaDE,absTol=0.0001d0)) then
                message='Omega_DE in merger tree file ['
                write (valueString,'(e14.8)') cosmologicalParameter
                message=message//trim(valueString)//'] differs from the internal value ['
                write (valueString,'(e14.8)') localOmegaDE
                message=message//trim(valueString)//']'
                if (self%fatalMismatches) then
                   call Galacticus_Error_Report('sussingOpen',message)
                else
                   call Galacticus_Display_Message(message,verbosityWarn)
                end if
             end if
          case ('B')
             read (parameterValue,*) self%boxLength
             unitString=String_Strip(parameterValue(index(parameterValue,' '):len(parameterValue)-index(parameterValue,' ')+1))
             if (String_Strip(unitString) /= 'Mpc/h') call Galacticus_Error_Report('sussingOpen','box length should be reported in units of Mpc/h')
             self%boxLengthUnits=importerUnits(.true.,megaParsec,-1,0)
          end select
       end if
    end do
    close(fileUnit)
    call parameterRegEx%destroy()
    return
  end subroutine sussingOpen

  subroutine sussingClose(self)
    !% Close a \glc\ format merger tree file.
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

  integer function sussingTreeCount(self)
    !% Return a count of the number of trees available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    call sussingTreeIndicesRead(self)
    sussingTreeCount=self%treesCount
    return
  end function sussingTreeCount

  integer(kind=kind_int8) function sussingTreeIndex(self,i)
    !% Return the index of the $i^{\rm th}$ tree.
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call sussingTreeIndicesRead(self)
    sussingTreeIndex=self%treeIndices(i)
    return
  end function sussingTreeIndex

  integer function sussingNodeCount(self,i)
    !% Return a count of the number of nodes in the $i^{\rm th}$ tree.
    implicit none
    class  (mergerTreeImporterSussing), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call sussingTreeIndicesRead(self)
    sussingNodeCount=self%treeSizes(i)
    return
  end function sussingNodeCount

  subroutine sussingTreeIndicesRead(self)
    !% Read the tree indices.
    use Galacticus_Display
    use Galacticus_Error
    use Kind_Numbers
    use String_Handling
    use Sort
    use File_Utilities
    use Memory_Management
    use Arrays_Search
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use, intrinsic :: ISO_C_Binding
    implicit none
    class    (mergerTreeImporterSussing), intent(inout)               :: self
    integer  (kind=c_size_t            ), allocatable  , dimension(:) :: nodeIndexRanks
    integer  (kind=kind_int8           ), allocatable  , dimension(:) :: nodeSelfIndices           , nodeTreeIndices       , &
         &                                                               nodeDescendentLocations   , nodesInSubvolume      , &
         &                                                               nodesTmp                  , hostsInSubvolume
    logical                             , allocatable  , dimension(:) :: nodeIncomplete            , nodeIncompleteTmp
    integer  (kind=kind_int8           ), allocatable  , dimension(:) :: forestSnapshotHaloCount   , forestSnapshotHaloCountFirst, &
         &                                                               forestSnapshotHaloCountLast
    integer                             , parameter                   :: fileFormatVersionCurrent=1
    integer                             , parameter                   :: stepCountMaximum        =1000000
    logical                                                           :: nodeIsActive              , doBinaryConversion    , &
         &                                                               readBinary                , mergerTreeFileIsBinary, &
         &                                                               mergerTreeFileConvert     , processHalo
    integer                                                           :: fileUnit                  , progenitorCount       , &
         &                                                               ioStat                    , lineStat              , &
         &                                                               nodeCountSubvolume        , nodeCount             , &
         &                                                               fileFormatVersion         , iProgenitor           , &
         &                                                               snapshotUnit              , snapshotOutUnit       , &
         &                                                               fileUnitOut               , nodeCountTrees
    character(len=1024                 )                              :: line
    character(len=32                   )                              :: label
    integer  (kind=c_size_t            )                              :: l
    integer  (kind=kind_int8           )                              :: nodeIndex                 , treeIndexPrevious     , &
         &                                                               treeIndexCurrent          , i                     , &
         &                                                               j                         , k                     , &
         &                                                               iNode                     , iStart                , &
         &                                                               jNode                     , iCount
    type     (varying_string           )                              :: message
    integer  (kind=kind_int8           )                              :: ID                        , hostHalo              , &
         &                                                               treeIndexFrom             , treeIndexTo           , &
         &                                                               progenitorIndex           , forestCount           , &
         &                                                               forestFirst               , forestLast            , &
         &                                                               forestHaloCountLast       , forestHaloCountFirst  , &
         &                                                               forestHaloCount
    integer                                                           :: numSubStruct              , npart                 , &
         &                                                               descendentStepCount       , hostStepCount         , &
         &                                                               haloFormat
    double precision                                                  :: Mvir                      , Xc                    , &
               &                                                         Yc                        , Zc                    , &
               &                                                         VXc                       , Vyc                   , &
               &                                                         VZc                       , Rvir                  , &
               &                                                         Rmax                      , r2                    , &
               &                                                         mbp_offset                , com_offset            , &
               &                                                         Vmax                      , v_esc                 , &
               &                                                         sigV                      , lambda                , &
               &                                                         lambdaE                   , Lx                    , &
               &                                                         Ly                        , Lz                    , &
               &                                                         b                         , c                     , &
               &                                                         Eax                       , Eay                   , &
               &                                                         Eaz                       , Ebx                   , &
               &                                                         Eby                       , Ebz                   , &
               &                                                         Ecx                       , Ecy                   , &
               &                                                         Ecz                       , ovdens                , &
               &                                                         fMhires                   , Ekin                  , &
               &                                                         Epot                      , SurfP                 , &
               &                                                         Phi0                      , cNFW                  , &
               &                                                         nbins                     , FoFMass               , &
               &                                                         M_200Mean                 , M_200Crit             , &
               &                                                         M_TopHat                  , R_200Mean             , &
               &                                                         R_200Crit                 , R_TopHat              , &
               &                                                         HalfMassRadius            , sigV_200Mean          , &
               &                                                         sigV_200Crit              , sigV_TopHat           , &
               &                                                         Xcm                       , Ycm                   , &
               &                                                         Zcm                       , Xgroup                , &
               &                                                         Ygroup                    , Zgroup
    type            (importerUnits     )                              :: massUnits                 , lengthUnits           , &
         &                                                               velocityUnits
     
    ! Return if indices have been read previously.
    if (self%treeIndicesRead) return
    self%treeIndicesRead=.true.
    ! Display counter.
    call Galacticus_Display_Indent ('Parsing "Sussing Merger Trees" format merger tree file',verbosityWorking)
    ! If a forest field is provided, scan it now to find the ranges to read from subsequent files.
    if (mergerTreeImportSussingUseForestFile) then
       forestCount=Count_Lines_in_File(mergerTreeImportSussingForestFile,'#')
       call Alloc_Array(forestSnapshotHaloCount     ,shape(self%snapshotFileName))
       call Alloc_Array(forestSnapshotHaloCountFirst,shape(self%snapshotFileName))
       call Alloc_Array(forestSnapshotHaloCountLast ,shape(self%snapshotFileName))
       forestHaloCountFirst        =0
       forestHaloCountLast         =0
       forestSnapshotHaloCountFirst=0
       forestSnapshotHaloCountLast =0
       forestFirst=mergerTreeImportSussingForestFirst
       forestLast =mergerTreeImportSussingForestLast
       if (forestLast < 0) forestLast=forestCount
       open(newUnit=fileUnit,file=char(mergerTreeImportSussingForestFile),status='old',form='formatted')
       read (fileUnit,'(a)') line
       do i=1,forestLast
          read (fileUnit,*) ID,forestHaloCount,forestSnapshotHaloCount
          forestHaloCountLast        =forestHaloCountLast        +forestHaloCount
          forestSnapshotHaloCountLast=forestSnapshotHaloCountLast+forestSnapshotHaloCount
          if (i < forestFirst) then
             forestHaloCountFirst        =forestHaloCountLast        +1
             forestSnapshotHaloCountFirst=forestSnapshotHaloCountLast+1
          end if
       end do
       close(fileUnit)
       call Dealloc_Array(forestSnapshotHaloCount)
    end if
    ! Open the merger tree file.
    mergerTreeFileIsBinary=File_Exists(char(self%mergerTreeFile//".bin"))
    mergerTreeFileConvert =.false.
    if (mergerTreeFileIsBinary) then
       open   (newUnit=fileUnit   ,file=char(self%mergerTreeFile//".bin"),status='old'    ,form='unformatted',ioStat=ioStat)
    else
       open   (newUnit=fileUnit   ,file=char(self%mergerTreeFile        ),status='old'    ,form='formatted'  ,ioStat=ioStat)
       if (self%convertToBinary) then
          mergerTreeFileConvert=.true.
          open(newUnit=fileUnitOut,file=char(self%mergerTreeFile//".bin"),status='unknown',form='unformatted'              )
       end if
    end if
    ! Read header information.
    call Galacticus_Display_Message('Reading header',verbosityWorking)
    if (mergerTreeFileIsBinary) then
       read     (fileUnit     ,ioStat=ioStat) fileFormatVersion
       read     (fileUnit     ,ioStat=ioStat) nodeCount
    else
       read     (fileUnit   ,*,ioStat=ioStat) fileFormatVersion
       read     (fileUnit   ,*,ioStat=ioStat) line
       read     (fileUnit   ,*,ioStat=ioStat) nodeCount
       if (mergerTreeFileConvert) then
          write (fileUnitOut                ) fileFormatVersion
          write (fileUnitOut                ) nodeCount
       end if
    end if
    ! Validate file format version/
    if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('sussingTreeIndicesRead','incorrect file format version')
    ! Allocate storage for list of nodes in subvolume.
    nodeCountSubVolume=int(dble(nodeCount)/dble(self%subvolumeCount)**3)+1
    call Alloc_Array(nodesInSubvolume,[nodeCountSubVolume])
    call Alloc_Array(hostsInSubvolume,[nodeCountSubVolume])
    ! Read snapshot halo catalogs.
    call Galacticus_Display_Indent('Finding halos in subvolume from AHF format snapshot halo catalogs',verbosityWorking)
    j                 =0
    nodeCountSubVolume=0
    do i=1,size(self%snapshotFileName)
       call Galacticus_Display_Message(self%snapshotFileName(i),verbosityWorking)
       doBinaryConversion=.false.
       readBinary        =.false.
       if (File_Exists(self%snapshotFileName(i)//".bin")) then
          ! A binary version of this file exists, use it.
          open   (newUnit=snapshotUnit   ,file=char(self%snapshotFileName(i)//".bin"),status='old'    ,form='unformatted',ioStat=ioStat)
          readBinary=.true.
       else
          ! No binary version of this file exists, use the ASCII version.
          open   (newUnit=snapshotUnit   ,file=char(self%snapshotFileName(i)        ),status='old'    ,form='formatted'  ,ioStat=ioStat)
          if (self%convertToBinary.and..not.mergerTreeImportSussingUseForestFile) then
             ! Open a binary file to write the converted halo data to.
             open(newUnit=snapshotOutUnit,file=char(self%snapshotFileName(i)//".bin"),status='unknown',form='unformatted'              )
             doBinaryConversion=.true.
          end if
          read (snapshotUnit,'(a)',ioStat=ioStat) line
          ! Test for format of halo files here. If "cNFW(25)" appears, it's a "new" format file, if "cNFW(43)" appears it's old
          ! format. Otherwise, we don't recognize it
          if      (index(line,"cNFW(43)"   ) /= 0) then
             ! Old format file.
             haloFormat=sussingHaloFormatOld
          else if (index(line,"cNFW(25)"   ) /= 0) then
             ! New format file.
             haloFormat=sussingHaloFormatNew
          else if (index(line,"FoFMass(44)") /= 0) then
             ! New format file.
             haloFormat=sussingHaloFormatAll
          else
             ! Unrecognized format.
             call Galacticus_Error_Report('sussingTreeIndicesRead','unrecognized format for halo files')
          end if         
       end if
       iCount=0
       do while (ioStat == 0)
          ! Increment count of number of halos read.
          iCount=iCount+1
          processHalo=                                             &
               &       (                                           &
               &            mergerTreeImportSussingUseForestFile   &
               &        .and.                                      &
               &         iCount >= forestSnapshotHaloCountFirst(i) &
               &        .and.                                      &
               &         iCount <= forestSnapshotHaloCountLast (i) &
               &       )                                           &
               &      .or.                                         &
               &       .not.mergerTreeImportSussingUseForestFile
          if (readBinary) then
             if (processHalo) then
                if (self%binaryFormatOld) then
                   read          (snapshotUnit     ,ioStat=ioStat)    &
                        &   ID            ,                           &
                        &   hostHalo      ,                           &
                        &   numSubStruct  ,                           &
                        &   Mvir          ,                           &
                        &   npart         ,                           &
                        &   Xc            ,                           &
                        &   Yc            ,                           &
                        &   Zc            ,                           &
                        &   VXc           ,                           &
                        &   Vyc           ,                           &
                        &   VZc           ,                           &
                        &   Rvir          ,                           &
                        &   Rmax          ,                           &
                        &   r2            ,                           &
                        &   mbp_offset    ,                           &
                        &   com_offset    ,                           &
                        &   Vmax          ,                           &
                        &   v_esc         ,                           &
                        &   sigV          ,                           &
                        &   lambda        ,                           &
                        &   lambdaE       ,                           &
                        &   Lx            ,                           &
                        &   Ly            ,                           &
                        &   Lz            ,                           &
                        &   b             ,                           &
                        &   c             ,                           &
                        &   Eax           ,                           &
                        &   Eay           ,                           &
                        &   Eaz           ,                           &
                        &   Ebx           ,                           &
                        &   Eby           ,                           &
                        &   Ebz           ,                           &
                        &   Ecx           ,                           &
                        &   Ecy           ,                           &
                        &   Ecz           ,                           &
                        &   ovdens        ,                           &
                        &   nbins         ,                           &
                        &   fMhires       ,                           &
                        &   Ekin          ,                           &
                        &   Epot          ,                           &
                        &   SurfP         ,                           &
                        &   Phi0          ,                           &
                        &   cNFW
                else
                   read          (snapshotUnit     ,ioStat=ioStat) &
                        &   ID            ,                           &
                        &   hostHalo      ,                           &
                        &   numSubStruct  ,                           &
                        &   Mvir          ,                           &
                        &   npart         ,                           &
                        &   Xc            ,                           &
                        &   Yc            ,                           &
                        &   Zc            ,                           &
                        &   VXc           ,                           &
                        &   Vyc           ,                           &
                        &   VZc           ,                           &
                        &   Rvir          ,                           &
                        &   Rmax          ,                           &
                        &   r2            ,                           &
                        &   mbp_offset    ,                           &
                        &   com_offset    ,                           &
                        &   Vmax          ,                           &
                        &   v_esc         ,                           &
                        &   sigV          ,                           &
                        &   lambda        ,                           &
                        &   lambdaE       ,                           &
                        &   Lx            ,                           &
                        &   Ly            ,                           &
                        &   Lz            ,                           &
                        &   b             ,                           &
                        &   c             ,                           &
                        &   Eax           ,                           &
                        &   Eay           ,                           &
                        &   Eaz           ,                           &
                        &   Ebx           ,                           &
                        &   Eby           ,                           &
                        &   Ebz           ,                           &
                        &   Ecx           ,                           &
                        &   Ecy           ,                           &
                        &   Ecz           ,                           &
                        &   ovdens        ,                           &
                        &   nbins         ,                           &
                        &   fMhires       ,                           &
                        &   Ekin          ,                           &
                        &   Epot          ,                           &
                        &   SurfP         ,                           &
                        &   Phi0          ,                           &
                        &   cNFW          ,                           &
                        &   FoFMass       ,                           &
                        &   M_200Mean     ,                           &
                        &   M_200Crit     ,                           &
                        &   M_TopHat      ,                           &
                        &   R_200Mean     ,                           &
                        &   R_200Crit     ,                           &
                        &   R_TopHat      ,                           &
                        &   HalfMassRadius,                           &
                        &   sigV_200Mean  ,                           &
                        &   sigV_200Crit  ,                           &
                        &   sigV_TopHat   ,                           &
                        &   Xcm           ,                           &
                        &   Ycm           ,                           &
                        &   Zcm           ,                           &
                        &   Xgroup        ,                           &
                        &   Ygroup        ,                           &
                        &   Zgroup
                end if
             else
                read (snapshotUnit,ioStat=ioStat)
             end if
          else
             call sussingReadHaloASCII(                         &
                  &   haloFormat    ,                           &
                  &   snapshotUnit  ,                           &
                  &   ioStat        ,                           &
                  &   ID            ,                           &
                  &   hostHalo      ,                           &
                  &   numSubStruct  ,                           &
                  &   Mvir          ,                           &
                  &   npart         ,                           &
                  &   Xc            ,                           &
                  &   Yc            ,                           &
                  &   Zc            ,                           &
                  &   VXc           ,                           &
                  &   Vyc           ,                           &
                  &   VZc           ,                           &
                  &   Rvir          ,                           &
                  &   Rmax          ,                           &
                  &   r2            ,                           &
                  &   mbp_offset    ,                           &
                  &   com_offset    ,                           &
                  &   Vmax          ,                           &
                  &   v_esc         ,                           &
                  &   sigV          ,                           &
                  &   lambda        ,                           &
                  &   lambdaE       ,                           &
                  &   Lx            ,                           &
                  &   Ly            ,                           &
                  &   Lz            ,                           &
                  &   b             ,                           &
                  &   c             ,                           &
                  &   Eax           ,                           &
                  &   Eay           ,                           &
                  &   Eaz           ,                           &
                  &   Ebx           ,                           &
                  &   Eby           ,                           &
                  &   Ebz           ,                           &
                  &   Ecx           ,                           &
                  &   Ecy           ,                           &
                  &   Ecz           ,                           &
                  &   ovdens        ,                           &
                  &   nbins         ,                           &
                  &   fMhires       ,                           &
                  &   Ekin          ,                           &
                  &   Epot          ,                           &
                  &   SurfP         ,                           &
                  &   Phi0          ,                           &
                  &   cNFW          ,                           &
                  &   FoFMass       ,                           &
                  &   M_200Mean     ,                           &
                  &   M_200Crit     ,                           &
                  &   M_TopHat      ,                           &
                  &   R_200Mean     ,                           &
                  &   R_200Crit     ,                           &
                  &   R_TopHat      ,                           &
                  &   HalfMassRadius,                           &
                  &   sigV_200Mean  ,                           &
                  &   sigV_200Crit  ,                           &
                  &   sigV_TopHat   ,                           &
                  &   Xcm           ,                           &
                  &   Ycm           ,                           &
                  &   Zcm           ,                           &
                  &   Xgroup        ,                           &
                  &   Ygroup        ,                           &
                  &   Zgroup        ,                           &
                  &   quickRead=.not.processHalo                &
                  & )
          end if
          if (ioStat /= 0) exit
          ! Write back in binary.
          if (doBinaryConversion)                            &
               &  write (snapshotOutUnit                   ) &
               &   ID            ,                           &
               &   hostHalo      ,                           &
               &   numSubStruct  ,                           &
               &   Mvir          ,                           &
               &   npart         ,                           &
               &   Xc            ,                           &
               &   Yc            ,                           &
               &   Zc            ,                           &
               &   VXc           ,                           &
               &   Vyc           ,                           &
               &   VZc           ,                           &
               &   Rvir          ,                           &
               &   Rmax          ,                           &
               &   r2            ,                           &
               &   mbp_offset    ,                           &
               &   com_offset    ,                           &
               &   Vmax          ,                           &
               &   v_esc         ,                           &
               &   sigV          ,                           &
               &   lambda        ,                           &
               &   lambdaE       ,                           &
               &   Lx            ,                           &
               &   Ly            ,                           &
               &   Lz            ,                           &
               &   b             ,                           &
               &   c             ,                           &
               &   Eax           ,                           &
               &   Eay           ,                           &
               &   Eaz           ,                           &
               &   Ebx           ,                           &
               &   Eby           ,                           &
               &   Ebz           ,                           &
               &   Ecx           ,                           &
               &   Ecy           ,                           &
               &   Ecz           ,                           &
               &   ovdens        ,                           &
               &   nbins         ,                           &
               &   fMhires       ,                           &
               &   Ekin          ,                           &
               &   Epot          ,                           &
               &   SurfP         ,                           &
               &   Phi0          ,                           &
               &   cNFW          ,                           &
               &   FoFMass       ,                           &
               &   M_200Mean     ,                           &
               &   M_200Crit     ,                           &
               &   M_TopHat      ,                           &
               &   R_200Mean     ,                           &
               &   R_200Crit     ,                           &
               &   R_TopHat      ,                           &
               &   HalfMassRadius,                           &
               &   sigV_200Mean  ,                           &
               &   sigV_200Crit  ,                           &
               &   sigV_TopHat   ,                           &
               &   Xcm           ,                           &
               &   Ycm           ,                           &
               &   Zcm           ,                           &
               &   Xgroup        ,                           &
               &   Ygroup        ,                           &
               &   Zgroup
           ! Check if halo is in our subvolume.
          if     (                                            &
               &   processHalo                                &
               &  .and.                                       &
               &   self%inSubVolume(Xc,Yc,Zc,buffered=.true.) &
               & ) then
             nodeCountSubvolume=nodeCountSubvolume+1
             if (nodeCountSubvolume > size(nodesInSubvolume)) then
                call Move_Alloc(nodesInSubvolume,nodesTmp)
                call Alloc_Array(nodesInSubvolume,int(dble(shape(nodesTmp))*1.4d0)+1)
                nodesInSubvolume(1:size(nodesTmp))=nodesTmp
                call Dealloc_Array(nodesTmp)
                call Move_Alloc(hostsInSubvolume,nodesTmp)
                call Alloc_Array(hostsInSubvolume,int(dble(shape(nodesTmp))*1.4d0)+1)
                hostsInSubvolume(1:size(nodesTmp))=nodesTmp
                call Dealloc_Array(nodesTmp)
             end if
             nodesInSubvolume(nodeCountSubvolume)=ID
             if (hostHalo <= 0) then
                hostsInSubvolume(nodeCountSubvolume)=ID
             else
                hostsInSubvolume(nodeCountSubvolume)=hostHalo
             end if
          end if
          ! Update the counter.
          j=j+1
          call Galacticus_Display_Counter(int(100.0d0*dble(j)/dble(nodeCount)),j == 1,verbosityWorking)
          ! If all required forests are processed, exit.
          if (mergerTreeImportSussingUseForestFile .and. iCount == forestSnapshotHaloCountLast(i)) exit
       end do
       close                        (snapshotUnit   )
       if (doBinaryConversion) close(snapshotOutUnit)
    end do
    call Galacticus_Display_Counter_Clear(verbosityWorking)
    call Sort_Do(nodesInSubvolume(1:nodeCountSubvolume),hostsInSubvolume(1:nodeCountSubvolume))
    message='Found '
    message=message//nodeCountSubvolume//' nodes in subvolume [from '//nodeCount//' total nodes]'
    call Galacticus_Display_Message(message,verbosityWorking)
    ! Allocate workspaces for merger trees.
    call Alloc_Array(nodeSelfIndices,[nodeCountSubvolume])
    ! Read node indices.
    call Galacticus_Display_Message("Reading node indices",verbosityWorking)
    i     =0
    iCount=0
    ioStat=0
    call Galacticus_Display_Counter(0,.true.,verbosityWorking)
    do while (ioStat == 0)
       if (mergerTreeFileIsBinary) then
          read         (fileUnit         ,ioStat=ioStat) nodeIndex,progenitorCount
          if (ioStat /= 0                         ) exit
       else
          read         (fileUnit   ,'(a)',ioStat=ioStat) line
          if (ioStat /= 0 .or. trim(line) == "END") exit
          read         (line       ,  *                ) nodeIndex,progenitorCount
          if (mergerTreeFileConvert)                                               &
               & write (fileUnitOut                    ) nodeIndex,progenitorCount
      end if
       ! Skip progenitors.
       if (mergerTreeFileIsBinary) then
          do j=1,progenitorCount
             read         (fileUnit     ,ioStat=ioStat) progenitorIndex
          end do
       else
          do j=1,progenitorCount
             read         (fileUnit   ,*,ioStat=ioStat) progenitorIndex
             if (mergerTreeFileConvert) &
                  & write (fileUnitOut                ) progenitorIndex
         end do
       end if
       ! Locate this node in the list of nodes in our subvolume.
       iNode=Search_Array(nodesInSubvolume(1:nodeCountSubvolume),nodeIndex)
       if (iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == nodeIndex) then
          ! This line represents a node in the tree.
          i                 =i+1
          nodeSelfIndices(i)=nodeIndex
       end if
       call Galacticus_Display_Counter(int(50.0d0*dble(i)/dble(nodeCountSubvolume)),.false.,verbosityWorking)
       if (i      == nodeCount          ) exit
       iCount=iCount+1
       if (mergerTreeImportSussingUseForestFile .and. iCount == forestHaloCountLast) exit
    end do
    close(fileUnit)
    if (mergerTreeFileConvert) close(fileUnitOut)
    ! Some halos in our subvolume might not be part of any tree. Adjust number of halos accordingly.
    nodeCountTrees=i
    if (nodeCountTrees < nodeCountSubvolume) then
       message='Found '
       message=message//nodeCountTrees//' nodes in subvolume trees [from '//nodeCountSubvolume//' total nodes in subvolume]'
       call Galacticus_Display_Message(message,verbosityWorking)
       call Move_Alloc(nodeSelfIndices,nodesTmp)       
       call Alloc_Array(nodeSelfIndices,[nodeCountTrees])
       nodeSelfIndices(1:nodeCountTrees)=nodesTmp(1:nodeCountTrees)
       call Dealloc_Array(nodesTmp)
    end if
    ! Allocate workspaces for merger trees.
    call Alloc_Array(nodeDescendentLocations,[nodeCountTrees])
    call Alloc_Array(nodeIncomplete         ,[nodeCountTrees])
    ! Get a sorted index into the list of nodes.
    call Galacticus_Display_Message('Building node index',verbosityWorking)
    nodeIndexRanks=Sort_Index_Do(nodeSelfIndices)
    ! Re-open the merger tree file. 
    mergerTreeFileIsBinary=File_Exists(char(self%mergerTreeFile//".bin"))
    if (mergerTreeFileIsBinary) then
       open(newUnit=fileUnit,file=char(self%mergerTreeFile//".bin"),status='old',form='unformatted',ioStat=ioStat)
    else
       open(newUnit=fileUnit,file=char(self%mergerTreeFile        ),status='old',form='formatted'  ,ioStat=ioStat)
    end if
    ! Read progenitor indices and make links.
    call Galacticus_Display_Message('Reading trees',verbosityWorking)
    if (mergerTreeFileIsBinary) then
       read (fileUnit  ,ioStat=ioStat) fileFormatVersion
       read (fileUnit  ,ioStat=ioStat) nodeCount
    else
       read (fileUnit,*,ioStat=ioStat) fileFormatVersion
       read (fileUnit,*,ioStat=ioStat) line
       read (fileUnit,*,ioStat=ioStat) nodeCount
    end if
    i                      = 0
    nodeDescendentLocations=-1
    nodeIncomplete         =.false.
    nodeIsActive           =.false.
    line                   =""
    do while (ioStat == 0)
       if (mergerTreeFileIsBinary) then
          read (fileUnit     ,ioStat=ioStat) nodeIndex,progenitorCount
          if (ioStat /= 0                         ) exit
       else
          read (fileUnit,'(a)',ioStat=ioStat) line
          if (ioStat /= 0 .or. trim(line) == "END") exit
          read (line,*) nodeIndex,progenitorCount
       end if
       ! Locate this node in the list of nodes in our subvolume.
       iNode=Search_Array(nodesInSubvolume(1:nodeCountSubvolume),nodeIndex)
       nodeIsActive=(iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == nodeIndex)
       if (nodeIsActive) then
          ! This line represents a node in the tree.
          i                 =i+1
          nodeSelfIndices(i)=nodeIndex
          call Galacticus_Display_Counter(int(50.0d0+50.0d0*dble(i)/dble(nodeCountTrees)),.false.,verbosityWorking)
       end if
       do j=1,progenitorCount
          if (mergerTreeFileIsBinary) then
             read (fileUnit  ,ioStat=ioStat) nodeIndex
          else
             read (fileUnit,*,ioStat=ioStat) nodeIndex
          end if          
          if (nodeIsActive) then
             ! This line represents a progenitor. Locate the progenitor in the list of halos.
             iProgenitor=Search_Indexed(nodeSelfIndices,nodeIndexRanks,nodeIndex)
             ! Does this progenitor exist within our subvolume?

if (nodeIndex==269000001798532_kind_int8) then
write (0,*) "CHECK ",nodeIndex,nodeSelfIndices(i),i,j,iprogenitor,nodeCountTrees,nodeIndexRanks(iProgenitor),size(nodeSelfIndices),nodeSelfIndices(nodeIndexRanks(iProgenitor))
end if

             if (iProgenitor <= 0 .or. iProgenitor > nodeCountTrees .or. nodeSelfIndices(nodeIndexRanks(iProgenitor)) /= nodeIndex) then
                nodeIncomplete(i)=.true.
             else
                if (nodeDescendentLocations(nodeIndexRanks(iProgenitor)) /= -1) then
                   message="multiple descendent trees not allowed"
                   message=message//char(10)//" first descendent: "//nodeSelfIndices(nodeDescendentLocations(nodeIndexRanks(iProgenitor)))
                   message=message//char(10)//"   new descendent: "//nodeSelfIndices(i)
                   message=message//char(10)//" progenitor index: "//nodeIndex
                   call Galacticus_Error_Report('sussingTreeIndicesRead',message)
                end if
                nodeDescendentLocations(nodeIndexRanks(iProgenitor))=i
                ! Find the progenitor node in the list of halos in the subvolume.
                iNode=Search_Array(nodesInSubvolume(1:nodeCountSubvolume),nodeIndex)
                if (iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == nodeIndex) then
                   ! Find hosted halos.
                   hostHalo=hostsInSubvolume(iNode)
                   if (hostHalo /= nodeIndex) then
                      ! Check if the host halo is in the subvolume.
                      jNode=Search_Array(nodesInSubvolume(1:nodeCountSubvolume),hostHalo)
                      if (jNode > 0 .and. jNode <= nodeCountSubvolume .and. nodesInSubvolume(jNode) == hostHalo) then
                         ! Check if the host halo is in the trees.
                         jNode=Search_Indexed(nodeSelfIndices,nodeIndexRanks,hostHalo)
                         if (.not.(jNode > 0 .and. jNode <= nodeCountTrees .and. nodeSelfIndices(nodeIndexRanks(jNode)) == hostHalo)) then
                            ! Host halo is in subvolume, but not in trees. Add it to the trees now.
                            message='host halo ['
                            message=message//hostHalo//'] in subvolume but not in trees - adding it now'
                            call Galacticus_Display_Message(message,verbosityWorking)
                            ! Expand arrays.
                            call Move_Alloc   (nodeSelfIndices        ,nodesTmp          )
                            call Alloc_Array  (nodeSelfIndices        ,[nodeCountTrees+1])
                            nodeSelfIndices        (1:nodeCountTrees)=nodesTmp
                            call Dealloc_Array(nodesTmp                                  )
                            call Move_Alloc   (nodeDescendentLocations,nodesTmp          )
                            call Alloc_Array  (nodeDescendentLocations,[nodeCountTrees+1])
                            nodeDescendentLocations(1:nodeCountTrees)=nodesTmp
                            call Dealloc_Array(nodesTmp                                  )
                            call Move_Alloc   (nodeIncomplete         ,nodeIncompleteTmp )
                            call Alloc_Array  (nodeIncomplete         ,[nodeCountTrees+1])
                            nodeIncomplete         (1:nodeCountTrees)=nodeIncompleteTmp
                            call Dealloc_Array(nodeIncompleteTmp                         )
                            ! Increment the number of halos in trees.
                            nodeCountTrees=nodeCountTrees+1
                            ! Insert the new halo, assigning the same descendent as its hosted halo.
                            nodeSelfIndices        (nodeCountTrees)=hostHalo
                            nodeDescendentLocations(nodeCountTrees)=i
                            nodeIncomplete         (nodeCountTrees)=.false.
                            ! Recompute the sort index into the node self indices.
                            deallocate(nodeIndexRanks)
                            nodeIndexRanks=Sort_Index_Do(nodeSelfIndices)
                         end if
                      end if
                   end if
                else
                   message='can not find halo ['
                   message=message//nodeIndex//'] in subvolume'
                   call Galacticus_Error_Report('sussingTreeIndicesRead',message)
                end if
             end if
          end if
       end do
    end do
    call Dealloc_Array(hostsInSubvolume)
    ! Close the merger tree file.
    close(fileUnit)
    ! Clear counter.
    call Galacticus_Display_Counter_Clear(       verbosityWorking)
    call Galacticus_Display_Unindent     ('done',verbosityWorking)
    ! Transfer tree structure to nodes array.
    allocate(self%nodes(nodeCountTrees))
    do i=1,nodeCountTrees
       self   %nodes(i)%nodeIndex      =nodeSelfIndices(                        i )
       if (nodeDescendentLocations(i) >= 0) then
          self%nodes(i)%descendentIndex=nodeSelfIndices(nodeDescendentLocations(i))
       else
          self%nodes(i)%descendentIndex=-1
       end if
    end do
    ! Read snapshot halo catalogs.
    call Galacticus_Display_Indent('Parsing AHF format snapshot halo catalogs',verbosityWorking)
    j=0
    do i=1,size(self%snapshotFileName)
       call Galacticus_Display_Message(self%snapshotFileName(i),verbosityWorking)
       readBinary=.false.
       if (File_Exists(self%snapshotFileName(i)//".bin")) then
          ! A binary version of this file exists, use it.
          open(newUnit=snapshotUnit,file=char(self%snapshotFileName(i)//".bin"),status='old',form='unformatted',ioStat=ioStat)
          readBinary=.true.
       else
          ! No binary version of this file exists, use the ASCII version.
          open(newUnit=snapshotUnit,file=char(self%snapshotFileName(i)        ),status='old',form='formatted'  ,ioStat=ioStat)
          read (snapshotUnit,*,ioStat=ioStat) line
       end if
       iCount=0
       do while (ioStat == 0)
          ! Increment count of number of halos read.
          iCount=iCount+1
          processHalo=                                             &
               &       (                                           &
               &            mergerTreeImportSussingUseForestFile   &
               &        .and.                                      &
               &         iCount >= forestSnapshotHaloCountFirst(i) &
               &        .and.                                      &
               &         iCount <= forestSnapshotHaloCountLast (i) &
               &       )                                           &
               &      .or.                                         &
               &       .not.mergerTreeImportSussingUseForestFile
          if (readBinary) then
             if (processHalo) then
                if (self%binaryFormatOld) then
                   read (snapshotUnit  ,ioStat=ioStat)                &
                        &   ID            ,                           &
                        &   hostHalo      ,                           &
                        &   numSubStruct  ,                           &
                        &   Mvir          ,                           &
                        &   npart         ,                           &
                        &   Xc            ,                           &
                        &   Yc            ,                           &
                        &   Zc            ,                           &
                        &   VXc           ,                           &
                        &   Vyc           ,                           &
                        &   VZc           ,                           &
                        &   Rvir          ,                           &
                        &   Rmax          ,                           &
                        &   r2            ,                           &
                        &   mbp_offset    ,                           &
                        &   com_offset    ,                           &
                        &   Vmax          ,                           &
                        &   v_esc         ,                           &
                        &   sigV          ,                           &
                        &   lambda        ,                           &
                        &   lambdaE       ,                           &
                        &   Lx            ,                           &
                        &   Ly            ,                           &
                        &   Lz            ,                           &
                        &   b             ,                           &
                        &   c             ,                           &
                        &   Eax           ,                           &
                        &   Eay           ,                           &
                        &   Eaz           ,                           &
                        &   Ebx           ,                           &
                        &   Eby           ,                           &
                        &   Ebz           ,                           &
                        &   Ecx           ,                           &
                        &   Ecy           ,                           &
                        &   Ecz           ,                           &
                        &   ovdens        ,                           &
                        &   nbins         ,                           &
                        &   fMhires       ,                           &
                        &   Ekin          ,                           &
                        &   Epot          ,                           &
                        &   SurfP         ,                           &
                        &   Phi0          ,                           &
                        &   cNFW
                else
                   read (snapshotUnit  ,ioStat=ioStat)                &
                        &   ID            ,                           &
                        &   hostHalo      ,                           &
                        &   numSubStruct  ,                           &
                        &   Mvir          ,                           &
                        &   npart         ,                           &
                        &   Xc            ,                           &
                        &   Yc            ,                           &
                        &   Zc            ,                           &
                        &   VXc           ,                           &
                        &   Vyc           ,                           &
                        &   VZc           ,                           &
                        &   Rvir          ,                           &
                        &   Rmax          ,                           &
                        &   r2            ,                           &
                        &   mbp_offset    ,                           &
                        &   com_offset    ,                           &
                        &   Vmax          ,                           &
                        &   v_esc         ,                           &
                        &   sigV          ,                           &
                        &   lambda        ,                           &
                        &   lambdaE       ,                           &
                        &   Lx            ,                           &
                        &   Ly            ,                           &
                        &   Lz            ,                           &
                        &   b             ,                           &
                        &   c             ,                           &
                        &   Eax           ,                           &
                        &   Eay           ,                           &
                        &   Eaz           ,                           &
                        &   Ebx           ,                           &
                        &   Eby           ,                           &
                        &   Ebz           ,                           &
                        &   Ecx           ,                           &
                        &   Ecy           ,                           &
                        &   Ecz           ,                           &
                        &   ovdens        ,                           &
                        &   nbins         ,                           &
                        &   fMhires       ,                           &
                        &   Ekin          ,                           &
                        &   Epot          ,                           &
                        &   SurfP         ,                           &
                        &   Phi0          ,                           &
                        &   cNFW          ,                           &
                        &   FoFMass       ,                           &
                        &   M_200Mean     ,                           &
                        &   M_200Crit     ,                           &
                        &   M_TopHat      ,                           &
                        &   R_200Mean     ,                           &
                        &   R_200Crit     ,                           &
                        &   R_TopHat      ,                           &
                        &   HalfMassRadius,                           &
                        &   sigV_200Mean  ,                           &
                        &   sigV_200Crit  ,                           &
                        &   sigV_TopHat   ,                           &
                        &   Xcm           ,                           &
                        &   Ycm           ,                           &
                        &   Zcm           ,                           &
                        &   Xgroup        ,                           &
                        &   Ygroup        ,                           &
                        &   Zgroup
                end if
             else
                read (snapshotUnit,ioStat=ioStat)
             end if
          else
             call sussingReadHaloASCII(                         &
                  &   haloFormat    ,                           &
                  &   snapshotUnit  ,                           &
                  &   ioStat        ,                           &
                  &   ID            ,                           &
                  &   hostHalo      ,                           &
                  &   numSubStruct  ,                           &
                  &   Mvir          ,                           &
                  &   npart         ,                           &
                  &   Xc            ,                           &
                  &   Yc            ,                           &
                  &   Zc            ,                           &
                  &   VXc           ,                           &
                  &   Vyc           ,                           &
                  &   VZc           ,                           &
                  &   Rvir          ,                           &
                  &   Rmax          ,                           &
                  &   r2            ,                           &
                  &   mbp_offset    ,                           &
                  &   com_offset    ,                           &
                  &   Vmax          ,                           &
                  &   v_esc         ,                           &
                  &   sigV          ,                           &
                  &   lambda        ,                           &
                  &   lambdaE       ,                           &
                  &   Lx            ,                           &
                  &   Ly            ,                           &
                  &   Lz            ,                           &
                  &   b             ,                           &
                  &   c             ,                           &
                  &   Eax           ,                           &
                  &   Eay           ,                           &
                  &   Eaz           ,                           &
                  &   Ebx           ,                           &
                  &   Eby           ,                           &
                  &   Ebz           ,                           &
                  &   Ecx           ,                           &
                  &   Ecy           ,                           &
                  &   Ecz           ,                           &
                  &   ovdens        ,                           &
                  &   nbins         ,                           &
                  &   fMhires       ,                           &
                  &   Ekin          ,                           &
                  &   Epot          ,                           &
                  &   SurfP         ,                           &
                  &   Phi0          ,                           &
                  &   cNFW          ,                           &
                  &   FoFMass       ,                           &
                  &   M_200Mean     ,                           &
                  &   M_200Crit     ,                           &
                  &   M_TopHat      ,                           &
                  &   R_200Mean     ,                           &
                  &   R_200Crit     ,                           &
                  &   R_TopHat      ,                           &
                  &   HalfMassRadius,                           &
                  &   sigV_200Mean  ,                           &
                  &   sigV_200Crit  ,                           &
                  &   sigV_TopHat   ,                           &
                  &   Xcm           ,                           &
                  &   Ycm           ,                           &
                  &   Zcm           ,                           &
                  &   Xgroup        ,                           &
                  &   Ygroup        ,                           &
                  &   Zgroup        ,                           &
                  &   quickRead=.not.processHalo                &
                  & )
          end if
          if (ioStat /= 0) exit
          ! Check if halo is to be processed.
          if (processHalo) then 
             ! Locate this node in the list of nodes in our subvolume.
             iNode=Search_Array(nodesInSubvolume(1:nodeCountSubvolume),ID)
             if (iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == ID) then
                ! Locate this node in the node list.
                l=Search_Indexed(nodeSelfIndices,nodeIndexRanks,ID)
                l=nodeIndexRanks(l)
                if (ID /= nodeSelfIndices(l)) then
                   if (mergerTreeImportSussingNonTreeNodeIsFatal) then
                      ! Node cannot be found.
                      message="node indexing failure"
                      message=message//char(10)//"     node index: "//ID
                      message=message//char(10)//"    found index: "//nodeSelfIndices(l)
                      message=message//char(10)//" found location: "//l
                      call Galacticus_Error_Report('sussingTreeIndicesRead',message)
                   else
                      ! Just skip this node.
                      cycle
                   end if
                end if
                ! Store properties to node array.
                if (hostHalo <= 0) then
                   self%nodes(l)%hostIndex         =ID
                else
                   self%nodes(l)%hostIndex         =hostHalo                
                   ! Check that the host halo is in the subvolume.
                   iNode=Search_Array(nodesInSubvolume(1:nodeCountSubvolume),hostHalo)
                   if (.not.(iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == hostHalo)) nodeIncomplete(l)=.true.
                end if
                self   %nodes(l)%particleCount     =npart
                ! Select a mass to use.
                select case (mergerTreeImportSussingMassOption)
                case (sussingMassDefault)
                   self%nodes(l)%nodeMass          =Mvir
                case (sussingMassFoF    )
                   self%nodes(l)%nodeMass          =FoFMass
                case (sussingMass200Mean)
                   self%nodes(l)%nodeMass          =M_200Mean
                case (sussingMass200Crit)
                   self%nodes(l)%nodeMass          =M_200Crit
                case (sussingMassTopHat )
                   self%nodes(l)%nodeMass          =M_TopHat
                end select
                if (self%nodes(l)%nodeMass == 0.0d0 .or. self%valueIsBad(self%nodes(l)%nodeMass)) self%nodes(l)%nodeMass=Mvir
                self   %nodes(l)%nodeTime          =self%snapshotTimes(i)
                if (.not.self%valueIsBad(cNFW)) then
                   if (cNFW > 0.0d0) then
                      self%nodes(l)%scaleRadius    =Rvir/cNFW
                   else
                      self%nodes(l)%scaleRadius    =0.0d0
                   end if
                else
                   self%scaleRadiiAvailableValue   =.false.
                   self%nodes(l)%scaleRadius       =-1.0d0
                end if
                self   %nodes(l)%halfMassRadius    =-1.0d0
                self   %nodes(l)%velocityMaximum   =Vmax
                self   %nodes(l)%velocityDispersion=sigV
                if (.not.self%valueIsBad(lambdaE)) then
                   self%nodes(l)%spin              =              lambdaE
                   self%nodes(l)%spin3D            =[Lx ,Ly ,Lz ]*lambdaE
                else
                   self%spinsAvailableValue        =.false.
                   self%nodes(l)%spin              =-1.0d0
                   self%nodes(l)%spin3D            =-1.0d0
                end if
                self   %nodes(l)%position          =[ Xc, Yc, Zc]
                self   %nodes(l)%velocity          =[VXc,VYc,VZc]
                ! Update the counter.
                j=j+1
                call Galacticus_Display_Counter(int(100.0d0*dble(j)/dble(nodeCountTrees)),j == 1,verbosityWorking)
             end if
          end if
          ! If all required forests are processed, exit.
          if (mergerTreeImportSussingUseForestFile .and. iCount == forestSnapshotHaloCountLast(i)) exit
       end do
       close(snapshotUnit)
    end do
    call Dealloc_Array(nodesInSubvolume)
    ! Clean up display.
    call Galacticus_Display_Counter_Clear(       verbosityWorking)
    call Galacticus_Display_Unindent     ('done',verbosityWorking)
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
    call Alloc_Array(nodeTreeIndices,[nodeCountTrees])
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
    call Galacticus_Display_Counter_Clear(       verbosityWorking)
    ! Check for nodes jumping between trees and join any such trees.
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
             where (nodeTreeIndices == treeIndexFrom)
                nodeTreeIndices=treeIndexTo
             end where
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
    massUnits                          =importerUnits(.true.,massSolar ,-1, 0)
    lengthUnits                        =importerUnits(.true.,kiloParsec,-1,+1)
    velocityUnits                      =importerUnits(.true.,kilo      , 0, 0)
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

  integer function sussingSubhaloTraceCount(self,node)
    !% Returns the length of a subhalo trace.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self
    class(nodeData                 ), intent(in   ) :: node

    ! No particle data is available.
    sussingSubhaloTraceCount=0
    return
  end function sussingSubhaloTraceCount

  subroutine sussingImport(self,i,nodes,requireScaleRadii,requireAngularMomenta,requireAngularMomenta3D,requireSpin,requireSpin3D,requirePositions,requireParticleCounts,requireVelocityMaxima,requireVelocityDispersions)
    !% Import the $i^{\rm th}$ merger tree.
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
    integer                                                                                 :: j

    ! Allocate the nodes array.
    allocate(nodeData :: nodes(self%treeSizes(i)))
    call Memory_Usage_Record(sizeof(nodes))
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
  
  subroutine sussingReadHaloASCII                 &
       &  (                                       &
       &   haloFormat    ,                        &
       &   snapshotUnit  ,                        &
       &   ioStat        ,                        &
       &   ID            ,                        &
       &   hostHalo      ,                        &
       &   numSubStruct  ,                        &
       &   Mvir          ,                        &
       &   npart         ,                        &
       &   Xc            ,                        &
       &   Yc            ,                        &
       &   Zc            ,                        &
       &   VXc           ,                        &
       &   Vyc           ,                        &
       &   VZc           ,                        &
       &   Rvir          ,                        &
       &   Rmax          ,                        &
       &   r2            ,                        &
       &   mbp_offset    ,                        &
       &   com_offset    ,                        &
       &   Vmax          ,                        &
       &   v_esc         ,                        &
       &   sigV          ,                        &
       &   lambda        ,                        &
       &   lambdaE       ,                        &
       &   Lx            ,                        &
       &   Ly            ,                        &
       &   Lz            ,                        &
       &   b             ,                        &
       &   c             ,                        &
       &   Eax           ,                        &
       &   Eay           ,                        &
       &   Eaz           ,                        &
       &   Ebx           ,                        &
       &   Eby           ,                        &
       &   Ebz           ,                        &
       &   Ecx           ,                        &
       &   Ecy           ,                        &
       &   Ecz           ,                        &
       &   ovdens        ,                        &
       &   nbins         ,                        &
       &   fMhires       ,                        &
       &   Ekin          ,                        &
       &   Epot          ,                        &
       &   SurfP         ,                        &
       &   Phi0          ,                        &
       &   cNFW          ,                        & 
       &   FoFMass       ,                        &
       &   M_200Mean     ,                        &
       &   M_200Crit     ,                        &
       &   M_TopHat      ,                        &
       &   R_200Mean     ,                        &
       &   R_200Crit     ,                        &
       &   R_TopHat      ,                        &
       &   HalfMassRadius,                        &
       &   sigV_200Mean  ,                        &
       &   sigV_200Crit  ,                        &
       &   sigV_TopHat   ,                        &
       &   Xcm           ,                        &
       &   Ycm           ,                        &
       &   Zcm           ,                        &
       &   Xgroup        ,                        &
       &   Ygroup        ,                        &
       &   Zgroup        ,                        &
       &   quickRead                              &
       & )
    !% Read an ASCII halo definition.
    use Galacticus_Error
    implicit none
    integer                         , intent(in   )           :: haloFormat    , snapshotUnit
    double precision                                          :: Mvir          , Xc          , &
         &                                                       Yc            , Zc          , &
         &                                                       VXc           , Vyc         , &
         &                                                       VZc           , Rvir        , &
         &                                                       Rmax          , r2          , &
         &                                                       mbp_offset    , com_offset  , &
         &                                                       Vmax          , v_esc       , &
         &                                                       sigV          , lambda      , &
         &                                                       lambdaE       , Lx          , &
         &                                                       Ly            , Lz          , &
         &                                                       b             , c           , &
         &                                                       Eax           , Eay         , &
         &                                                       Eaz           , Ebx         , &
         &                                                       Eby           , Ebz         , &
         &                                                       Ecx           , Ecy         , &
         &                                                       Ecz           , ovdens      , &
         &                                                       fMhires       , Ekin        , &
         &                                                       Epot          , SurfP       , &
         &                                                       Phi0          , cNFW        , &
         &                                                       nbins         , FoFMass     , &
         &                                                       M_200Mean     , M_200Crit   , &
         &                                                       M_TopHat      , R_200Mean   , &
         &                                                       R_200Crit     , R_TopHat    , &
         &                                                       HalfMassRadius, sigV_200Mean, &
         &                                                       sigV_200Crit  , sigV_TopHat , &
         &                                                       Xcm           , Ycm         , &
         &                                                       Zcm           , Xgroup      , &
         &                                                       Ygroup        , Zgroup
    integer         (kind=kind_int8), intent(  out)           :: ID            , hostHalo
    integer                         , intent(  out)           :: numSubStruct  , npart       , &
         &                                                       ioStat
    logical                         , intent(in   ), optional :: quickRead
    
    if (present(quickRead).and.quickRead) then
       read (snapshotUnit,*,ioStat=ioStat)
    else
       if (haloFormat == sussingHaloFormatOld) then
          read          (snapshotUnit   ,*,ioStat=ioStat) &
               &   ID          ,                          &
               &   hostHalo    ,                          &
               &   numSubStruct,                          &
               &   Mvir        ,                          &
               &   npart       ,                          &
               &   Xc          ,                          &
               &   Yc          ,                          &
               &   Zc          ,                          &
               &   VXc         ,                          &
               &   Vyc         ,                          &
               &   VZc         ,                          &
               &   Rvir        ,                          &
               &   Rmax        ,                          &
               &   r2          ,                          &
               &   mbp_offset  ,                          &
               &   com_offset  ,                          &
               &   Vmax        ,                          &
               &   v_esc       ,                          &
               &   sigV        ,                          &
               &   lambda      ,                          &
               &   lambdaE     ,                          &
               &   Lx          ,                          &
               &   Ly          ,                          &
               &   Lz          ,                          &
               &   b           ,                          &
               &   c           ,                          &
               &   Eax         ,                          &
               &   Eay         ,                          &
               &   Eaz         ,                          &
               &   Ebx         ,                          &
               &   Eby         ,                          &
               &   Ebz         ,                          &
               &   Ecx         ,                          &
               &   Ecy         ,                          &
               &   Ecz         ,                          &
               &   ovdens      ,                          &
               &   nbins       ,                          &
               &   fMhires     ,                          &
               &   Ekin        ,                          &
               &   Epot        ,                          &
               &   SurfP       ,                          &
               &   Phi0        ,                          &
               &   cNFW
          FoFMass       =0.0d0
          M_200Mean     =0.0d0
          M_200Crit     =0.0d0
          M_TopHat      =0.0d0
          R_200Mean     =0.0d0
          R_200Crit     =0.0d0
          R_TopHat      =0.0d0
          HalfMassRadius=0.0d0
          sigV_200Mean  =0.0d0
          sigV_200Crit  =0.0d0
          sigV_TopHat   =0.0d0
          Xcm           =0.0d0
          Ycm           =0.0d0
          Zcm           =0.0d0
          Xgroup        =0.0d0
          Ygroup        =0.0d0
          Zgroup        =0.0d0
        else if (haloFormat == sussingHaloFormatNew) then
          read          (snapshotUnit   ,*,ioStat=ioStat) &
               &   ID          ,                          &
               &   hostHalo    ,                          &
               &   numSubStruct,                          &
               &   Mvir        ,                          &
               &   npart       ,                          &
               &   Xc          ,                          &
               &   Yc          ,                          &
               &   Zc          ,                          &
               &   VXc         ,                          &
               &   Vyc         ,                          &
               &   VZc         ,                          &
               &   Rvir        ,                          &
               &   Rmax        ,                          &
               &   r2          ,                          &
               &   mbp_offset  ,                          &
               &   com_offset  ,                          &
               &   Vmax        ,                          &
               &   v_esc       ,                          &
               &   sigV        ,                          &
               &   lambda      ,                          &
               &   lambdaE     ,                          &
               &   Lx          ,                          &
               &   Ly          ,                          &
               &   Lz          ,                          &
               &   cNFW
          b             =0.0d0
          c             =0.0d0
          Eax           =0.0d0
          Eay           =0.0d0
          Eaz           =0.0d0
          Ebx           =0.0d0
          Eby           =0.0d0
          Ebz           =0.0d0
          Ecx           =0.0d0
          Ecy           =0.0d0
          Ecz           =0.0d0
          ovdens        =0.0d0
          nbins         =0.0d0
          fMhires       =0.0d0
          Ekin          =0.0d0
          Epot          =0.0d0
          SurfP         =0.0d0
          Phi0          =0.0d0
          FoFMass       =0.0d0
          M_200Mean     =0.0d0
          M_200Crit     =0.0d0
          M_TopHat      =0.0d0
          R_200Mean     =0.0d0
          R_200Crit     =0.0d0
          R_TopHat      =0.0d0
          HalfMassRadius=0.0d0
          sigV_200Mean  =0.0d0
          sigV_200Crit  =0.0d0
          sigV_TopHat   =0.0d0
          Xcm           =0.0d0
          Ycm           =0.0d0
          Zcm           =0.0d0
          Xgroup        =0.0d0
          Ygroup        =0.0d0
          Zgroup        =0.0d0
      else if (haloFormat == sussingHaloFormatAll) then
          read          (snapshotUnit   ,*,ioStat=ioStat) &
               &   ID            ,                        &
               &   hostHalo      ,                        &
               &   numSubStruct  ,                        &
               &   Mvir          ,                        &
               &   npart         ,                        &
               &   Xc            ,                        &
               &   Yc            ,                        &
               &   Zc            ,                        &
               &   VXc           ,                        &
               &   Vyc           ,                        &
               &   VZc           ,                        &
               &   Rvir          ,                        &
               &   Rmax          ,                        &
               &   r2            ,                        &
               &   mbp_offset    ,                        &
               &   com_offset    ,                        &
               &   Vmax          ,                        &
               &   v_esc         ,                        &
               &   sigV          ,                        &
               &   lambda        ,                        &
               &   lambdaE       ,                        &
               &   Lx            ,                        &
               &   Ly            ,                        &
               &   Lz            ,                        &
               &   b             ,                        &
               &   c             ,                        &
               &   Eax           ,                        &
               &   Eay           ,                        &
               &   Eaz           ,                        &
               &   Ebx           ,                        &
               &   Eby           ,                        &
               &   Ebz           ,                        &
               &   Ecx           ,                        &
               &   Ecy           ,                        &
               &   Ecz           ,                        &
               &   ovdens        ,                        &
               &   nbins         ,                        &
               &   fMhires       ,                        &
               &   Ekin          ,                        &
               &   Epot          ,                        &
               &   SurfP         ,                        &
               &   Phi0          ,                        &
               &   cNFW          ,                        &
               &   FoFMass       ,                        &
               &   M_200Mean     ,                        &
               &   M_200Crit     ,                        &
               &   M_TopHat      ,                        &
               &   R_200Mean     ,                        &
               &   R_200Crit     ,                        &
               &   R_TopHat      ,                        &
               &   HalfMassRadius,                        &
               &   sigV_200Mean  ,                        &
               &   sigV_200Crit  ,                        &
               &   sigV_TopHat   ,                        &
               &   Xcm           ,                        &
               &   Ycm           ,                        &
               &   Zcm           ,                        &
               &   Xgroup        ,                        &
               &   Ygroup        ,                        &
               &   Zgroup
    else
          call Galacticus_Error_Report('sussingReadHaloASCII','unknown halo file format')
       end if
    end if
    return
  end subroutine sussingReadHaloASCII
