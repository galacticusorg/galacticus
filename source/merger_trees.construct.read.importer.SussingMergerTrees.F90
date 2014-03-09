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
     logical                                                     :: fatalMismatches, treeIndicesRead
     integer                                                     :: treesCount
     double precision                                            :: boxLength
     type            (importerUnits )                            :: boxLengthUnits
     type            (varying_string)                            :: mergerTreeFile
     type            (varying_string), allocatable, dimension(:) :: snapshotFileName
     integer         (kind=kind_int8), allocatable, dimension(:) :: treeIndices
     integer         (kind=c_size_t ), allocatable, dimension(:) :: treeIndexRanks
     integer                         , allocatable, dimension(:) :: treeSizes      , treeBegins
     type            (nodeData      ), allocatable, dimension(:) :: nodes
     double precision                , allocatable, dimension(:) :: snapshotTimes
    contains
     !# <workaround type="gfortran" PR="58471 58470" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58471 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58470">
     !# final     :: sussingDestructor
     !# </workaround>
     procedure :: open                        => sussingOpen
     procedure :: close                       => sussingClose
     procedure :: treesHaveSubhalos           => sussingTreesHaveSubhalos
     procedure :: massesIncludeSubhalos       => sussingMassesIncludeSubhalos
     procedure :: treesAreSelfContained       => sussingTreesAreSelfContained
     procedure :: velocitiesIncludeHubbleFlow => sussingVelocitiesIncludeHubbleFlow
     procedure :: positionsArePeriodic        => sussingPositionsArePeriodic
     procedure :: cubeLength                  => sussingCubeLength
     procedure :: treeWeight                  => sussingTreeWeight
     procedure :: treeCount                   => sussingTreeCount
     procedure :: treeIndex                   => sussingTreeIndex
     procedure :: nodeCount                   => sussingNodeCount
     procedure :: positionsAvailable          => sussingPositionsAvailable
     procedure :: scaleRadiiAvailable         => sussingScaleRadiiAvailable
     procedure :: particleCountAvailable      => sussingParticleCountAvailable
     procedure :: velocityMaximumAvailable    => sussingVelocityMaximumAvailable
     procedure :: velocityDispersionAvailable => sussingVelocityDispersionAvailable
     procedure :: angularMomentaAvailable     => sussingAngularMomentaAvailable
     procedure :: angularMomenta3DAvailable   => sussingAngularMomenta3DAvailable
     procedure :: spinAvailable               => sussingSpinAvailable
     procedure :: spin3DAvailable             => sussingSpin3DAvailable
     procedure :: import                      => sussingImport
     procedure :: subhaloTrace                => sussingSubhaloTrace
     procedure :: subhaloTraceCount           => sussingSubhaloTraceCount
  end type mergerTreeImporterSussing

  interface mergerTreeImporterSussing
     !% Constructors for the \glc\ format merger tree importer class.
     module procedure sussingDefaultConstructor
  end interface mergerTreeImporterSussing

  ! Record of implementation initialization state.
  logical :: sussingInitialized                     =.false.

  ! Default settings.
  logical :: mergerTreeImportSussingMismatchIsFatal

contains

  function sussingDefaultConstructor()
    !% Default constructor for the ``Sussing Merger Trees'' format \citep{srisawat_sussing_2013} merger tree importer.
    use Input_Parameters
    implicit none
    type (mergerTreeImporterSussing), target :: sussingDefaultConstructor
        
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
          sussingInitialized=.true.
       end if
       !$omp end critical (mergerTreeImporterSussingInitialize)
    end if
    sussingDefaultConstructor%treeIndicesRead=.false.
    sussingDefaultConstructor%fatalMismatches=mergerTreeImportSussingMismatchIsFatal
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
    call Alloc_Array(self%snapshotTimes,[snapshotFileCount],lowerBounds=[0])
    open(newUnit=fileUnit,file=char(snapshotTimesFile),status='old',form='formatted',ioStat=ioStat)
    read (fileUnit,*)
    do i=1,snapshotFileCount
       read (fileUnit,*) snapshotNumber,expansionFactor,redshift,timeNormalized,time
       self%snapshotTimes(snapshotNumber)=                               &
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
    use Memory_Management
    use Arrays_Search
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use, intrinsic :: ISO_C_Binding
    implicit none
    class    (mergerTreeImporterSussing), intent(inout)               :: self
    integer  (kind=c_size_t            ), allocatable  , dimension(:) :: nodeIndexRanks            , nodeNodeLocations
    integer  (kind=kind_int8           ), allocatable  , dimension(:) :: nodeSelfIndices           , nodeTreeIndices  , &
         &                                                               nodeDescendentLocations
    integer                             , parameter                   :: fileFormatVersionCurrent=1
    integer                                                           :: fileUnit                  , progenitorCount  , &
         &                                                               ioStat                    , lineStat         , &
         &                                                               i                         , nodeCount        , &
         &                                                               fileFormatVersion         , iProgenitor      , &
         &                                                               l                         , j                , &
         &                                                               snapshotUnit              , k
    character(len=32                   )                              :: line
    integer  (kind=kind_int8           )                              :: nodeIndex                 , treeIndexPrevious
    type     (varying_string           )                              :: message
    integer  (kind=kind_int8           )                              :: ID                        , hostHalo         , &
         &                                                               treeIndexFrom             , treeIndexTo
    integer                                                           :: numSubStruct              , npart            , &
         &                                                               nbins
    double precision                                                  :: Mvir                      , Xc               , &
               &                                                         Yc                        , Zc               , &
               &                                                         VXc                       , Vyc              , &
               &                                                         VZc                       , Rvir             , &
               &                                                         Rmax                      , r2               , &
               &                                                         mbp_offset                , com_offset       , &
               &                                                         Vmax                      , v_esc            , &
               &                                                         sigV                      , lambda           , &
               &                                                         lambdaE                   , Lx               , &
               &                                                         Ly                        , Lz               , &
               &                                                         b                         , c                , &
               &                                                         Eax                       , Eay              , &
               &                                                         Eaz                       , Ebx              , &
               &                                                         Eby                       , Ebz              , &
               &                                                         Ecx                       , Ecy              , &
               &                                                         Ecz                       , ovdens           , &
               &                                                         fMhires                   , Ekin             , &
               &                                                         Epot                      , SurfP            , &
               &                                                         Phi0                      , cNFW
    type            (importerUnits     )                              :: massUnits                 , lengthUnits      , &
         &                                                               velocityUnits
     
    ! Return if indices have been read previously.
    if (self%treeIndicesRead) return
    self%treeIndicesRead=.true.
    ! Display counter.
    call Galacticus_Display_Indent ('Parsing "Sussing Merger Trees" format merger tree file',verbosityWorking)
    call Galacticus_Display_Counter(0,.true.,verbosityWorking)
    ! Open the merger tree file.
    open(newUnit=fileUnit,file=char(self%mergerTreeFile),status='old',form='formatted',ioStat=ioStat)
    ! Read header information.
    call Galacticus_Display_Message('Reading nodes',verbosityWorking)
    read (fileUnit,*,ioStat=ioStat) fileFormatVersion
    read (fileUnit,*,ioStat=ioStat) line
    read (fileUnit,*,ioStat=ioStat) nodeCount
    ! Validate file format version/
    if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('sussingTreeIndicesRead','incorrect file format version')
    ! Allocate workspaces for merger trees.
    call Alloc_Array(nodeSelfIndices        ,[nodeCount])
    call Alloc_Array(nodeTreeIndices        ,[nodeCount])
    call Alloc_Array(nodeNodeLocations      ,[nodeCount])
    call Alloc_Array(nodeDescendentLocations,[nodeCount])
    message='  -> Found '
    message=message//nodeCount//' total nodes'
    call Galacticus_Display_Message(message,verbosityWorking)
    ! Read node indices.
    i     =0
    ioStat=0
    do while (ioStat == 0)
       read (fileUnit,'(a)',ioStat=ioStat) line
       if (ioStat /= 0) exit
       read (line,*,ioStat=lineStat) nodeIndex,progenitorCount
       if (lineStat == 0) then
          ! This line represents a node in the tree.
          i                 =i+1
          nodeSelfIndices(i)=nodeIndex
          call Galacticus_Display_Counter(int(50.0d0*dble(i)/dble(nodeCount)),.false.,verbosityWorking)
       end if
    end do
    rewind(fileUnit)
    ! Get a sorted index into the list of nodes.
    call Galacticus_Display_Message('Building node index',verbosityWorking)
    nodeIndexRanks=Sort_Index_Do(nodeSelfIndices)
    ! Read progenitor indices and make links.
    call Galacticus_Display_Message('Reading trees',verbosityWorking)
    read (fileUnit,*,ioStat=ioStat) fileFormatVersion
    read (fileUnit,*,ioStat=ioStat) line
    read (fileUnit,*,ioStat=ioStat) nodeCount
    i                      = 0
    nodeDescendentLocations=-1
    do while (ioStat == 0)
       read (fileUnit,'(a)',ioStat=ioStat) line
       if (ioStat /= 0 .or. trim(line) == "END") exit
       read (line,*,ioStat=lineStat) nodeIndex,progenitorCount
       if (lineStat == 0) then
          ! This line represents a node in the tree.
          i                 =i+1
          nodeSelfIndices(i)=nodeIndex
          call Galacticus_Display_Counter(int(50.0d0+50.0d0*dble(i)/dble(nodeCount)),.false.,verbosityWorking)
       else
          ! This line represents a progenitor. Locate the progenitor in the list of halos.
          iProgenitor=Search_Indexed(nodeSelfIndices,nodeIndexRanks,nodeIndex)
          if (nodeDescendentLocations(nodeIndexRanks(iProgenitor)) /= -1) call Galacticus_Error_Report('sussingTreeIndicesRead','multiple descendent trees not allowed')
          nodeDescendentLocations(nodeIndexRanks(iProgenitor))=i
       end if
    end do
    ! Clear counter.
    call Galacticus_Display_Counter_Clear(       verbosityWorking)
    call Galacticus_Display_Unindent     ('done',verbosityWorking)
    ! Transfer tree structure to nodes array.
    allocate(self%nodes(nodeCount))
    do i=1,nodeCount
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
       open(newUnit=snapshotUnit,file=char(self%snapshotFileName(i)),status='old',form='formatted',ioStat=ioStat)
       read (snapshotUnit,*,ioStat=ioStat) line
       do while (ioStat == 0)
          read (snapshotUnit,*,ioStat=ioStat) &
               & ID          ,                &
               & hostHalo    ,                &
               & numSubStruct,                &
               & Mvir        ,                &
               & npart       ,                &
               & Xc          ,                &
               & Yc          ,                &
               & Zc          ,                &
               & VXc         ,                &
               & Vyc         ,                &
               & VZc         ,                &
               & Rvir        ,                &
               & Rmax        ,                &
               & r2          ,                &
               & mbp_offset  ,                &
               & com_offset  ,                &
               & Vmax        ,                &
               & v_esc       ,                &
               & sigV        ,                &
               & lambda      ,                &
               & lambdaE     ,                &
               & Lx          ,                &
               & Ly          ,                &
               & Lz          ,                &
               & b           ,                &
               & c           ,                &
               & Eax         ,                &
               & Eay         ,                &
               & Eaz         ,                &
               & Ebx         ,                &
               & Eby         ,                &
               & Ebz         ,                &
               & Ecx         ,                &
               & Ecy         ,                &
               & Ecz         ,                &
               & ovdens      ,                &
               & nbins       ,                &
               & fMhires     ,                &
               & Ekin        ,                &
               & Epot        ,                &
               & SurfP       ,                &
               & Phi0        ,                &
               & cNFW
          if (ioStat /= 0) exit
          ! Locate this node in the node list.
          l=Search_Indexed(nodeSelfIndices,nodeIndexRanks,ID)
          l=nodeIndexRanks(l)
          if (ID /= nodeSelfIndices(l)) call Galacticus_Error_Report('sussingTreeIndicesRead','node indexing failure')
          ! Store properties to node array.
          if (hostHalo <= 0) then
             self%nodes(l)%hostIndex         =ID
          else
             self%nodes(l)%hostIndex         =hostHalo
          end if
          self   %nodes(l)%particleCount     =npart
          self   %nodes(l)%nodeMass          =Mvir
          self   %nodes(l)%nodeTime          =self%snapshotTimes(i-1)
          self   %nodes(l)%scaleRadius       =Rvir/cNFW
          self   %nodes(l)%halfMassRadius    =-1.0d0
          self   %nodes(l)%velocityMaximum   =Vmax
          self   %nodes(l)%velocityDispersion=sigV
          self   %nodes(l)%spin              =              lambdaE
          self   %nodes(l)%spin3D            =[Lx ,Ly ,Lz ]*lambdaE
          self   %nodes(l)%position          =[ Xc, Yc, Zc]
          self   %nodes(l)%velocity          =[VXc,VYc,VZc]
          ! Update the counter.
          j=j+1
          call Galacticus_Display_Counter(int(100.0d0*dble(j)/dble(nodeCount)),j == 1,verbosityWorking)
       end do
       close(fileUnit)
    end do
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
    ! Clean up display.
    call Galacticus_Display_Counter_Clear(       verbosityWorking)
    call Galacticus_Display_Unindent     ('done',verbosityWorking)
    ! Assign tree indices.
    call Galacticus_Display_Message('Assigning tree indices',verbosityWorking)
    nodeTreeIndices=-1
    do i=1,nodeCount
       call Galacticus_Display_Counter(int(50.0d0*dble(i)/dble(nodeCount)),i==1,verbosityWorking)
       l=nodeDescendentLocations(i)
       if (l == -1 .and. self%nodes(i)%hostIndex == self%nodes(i)%nodeIndex) nodeTreeIndices(i)=nodeSelfIndices(i)
    end do
    do i=1,nodeCount
       call Galacticus_Display_Counter(int(50.0d0+50.0d0*dble(i)/dble(nodeCount)),.false.,verbosityWorking)
       l=i
       do while (l /= -1 .and. self%nodes(l)%hostIndex /= self%nodes(l)%nodeIndex) 
          k=Search_Indexed(nodeSelfIndices,nodeIndexRanks,self%nodes(l)%hostIndex)
          l=nodeIndexRanks(k)
       end do
       do while (l /= -1)
          nodeTreeIndices(i)=nodeTreeIndices(l)
          l=nodeDescendentLocations(l)
          do while (l /= -1 .and. self%nodes(l)%hostIndex /= self%nodes(l)%nodeIndex) 
             k=Search_Indexed(nodeSelfIndices,nodeIndexRanks,self%nodes(l)%hostIndex)
             l=nodeIndexRanks(k)
          end do
       end do
    end do
    ! Check for nodes jumping between trees and join any such trees.
    do i=1,nodeCount
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
    do i=1,nodeCount
       l=nodeDescendentLocations(i)
       if (l /= -1) then
          if (nodeTreeIndices(i) /= nodeTreeIndices(l))                                              &
               & call Galacticus_Error_Report('sussingTreeIndicesRead','failed to cross-link trees')
       end if
    end do
    call Galacticus_Display_Counter_Clear(verbosityWorking)
    ! Close the merger tree file.
    close(fileUnit)
    ! Generate an index into nodes sorted by tree index.
    self%treeIndexRanks=Sort_Index_Do(nodeTreeIndices)
    ! Create a list of tree indices, sizes, and start locations.
    self%treesCount=0
    treeIndexPrevious=-1
    do i=1,nodeCount
       if (nodeTreeIndices(self%treeIndexRanks(i)) /= treeIndexPrevious) then
          treeIndexPrevious=nodeTreeIndices(self%treeIndexRanks(i))
          self%treesCount=self%treesCount+1
       end if
    end do
    message='Found '
    message=message//self%treesCount//' trees'
    call Galacticus_Display_Message(message,verbosityWorking)
    call Alloc_Array(self%treeIndices,[self%treesCount])
    call Alloc_Array(self%treeSizes  ,[self%treesCount])
    call Alloc_Array(self%treeBegins ,[self%treesCount])
    treeIndexPrevious=-1
    j                = 0
    self%treeSizes   = 0
    do i=1,nodeCount
       if (nodeTreeIndices(self%treeIndexRanks(i)) /= treeIndexPrevious) then
          treeIndexPrevious=nodeTreeIndices(self%treeIndexRanks(i))
          j=j+1
          self%treeIndices(j)=nodeTreeIndices(self%treeIndexRanks(i))
          self%treeBegins (j)=i
       end if
       self%treeSizes(j)=self%treeSizes(j)+1
    end do
    ! Clean up display.
    call Galacticus_Display_Counter_Clear(verbosityWorking)
    ! Destroy temporary workspace.
    call Dealloc_Array(nodeSelfIndices        )
    call Dealloc_Array(nodeTreeIndices        )
    call Dealloc_Array(nodeNodeLocations      )
    call Dealloc_Array(nodeDescendentLocations)
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

    sussingScaleRadiiAvailable=.true.
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
    
    sussingSpinAvailable=.true.
    return
  end function sussingSpinAvailable

  logical function sussingSpin3DAvailable(self)
    !% Return true if spins vectors are available.
    implicit none
    class(mergerTreeImporterSussing), intent(inout) :: self

    sussingSpin3DAvailable=.true.
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
