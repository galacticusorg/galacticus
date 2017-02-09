!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !# <mergerTreeImporter name="mergerTreeImporterSussingASCII" description="Importer for ``Sussing Merger Trees'' ASCII format merger tree files \citep{srisawat_sussing_2013}." />

  type, extends(mergerTreeImporterSussing) :: mergerTreeImporterSussingASCII
     !% A merger tree importer class for ``Sussing Merger Trees'' ASCII format merger tree files \citep{srisawat_sussing_2013}.
     private
   contains
     final     ::         sussingASCIIDestructor
     procedure :: open => sussingASCIIOpen
     procedure :: load => sussingASCIILoad
  end type mergerTreeImporterSussingASCII

  interface mergerTreeImporterSussingASCII
     !% Constructors for the {\normalfont \ttfamily sussing} ASCII format merger tree importer class.
     module procedure sussingASCIIDefaultConstructor
  end interface mergerTreeImporterSussingASCII

  ! Record of implementation initialization state.
  logical :: sussingASCIIInitialized=.false.

  ! Default settings.
  logical                            :: mergerTreeImportSussingUseForestFile
  type   (varying_string)            :: mergerTreeImportSussingForestFile
  integer                            :: mergerTreeImportSussingForestFirst
  integer                            :: mergerTreeImportSussingForestLast
  logical                            :: mergerTreeImportSussingConvertToBinary
  logical                            :: mergerTreeImportSussingBinaryFormatOld
  logical                            :: mergerTreeImportSussingForestReverseSnapshotOrder
 
  ! File format identifiers.
  integer                , parameter :: sussingHaloFormatOld                             =1
  integer                , parameter :: sussingHaloFormatNew                             =2
  integer                , parameter :: sussingHaloFormatAll                             =3
 
contains

  function sussingASCIIDefaultConstructor()
    !% Default constructor for the ``Sussing Merger Trees'' ASCII format \citep{srisawat_sussing_2013} merger tree importer.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(mergerTreeImporterSussingASCII), target :: sussingASCIIDefaultConstructor

    if (.not.sussingASCIIInitialized) then
       !$omp critical (mergerTreeImporterSussingASCIIInitialize)
       if (.not.sussingASCIIInitialized) then
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
          sussingASCIIInitialized=.true.
       end if
       !$omp end critical (mergerTreeImporterSussingASCIIInitialize)
    end if
    call sussingInitialize(sussingASCIIDefaultConstructor)
    return
  end function sussingASCIIDefaultConstructor

  subroutine sussingASCIIDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sussing} ASCII format merger tree importer class.
    use Memory_Management
    implicit none
    type(mergerTreeImporterSussingASCII), intent(inout) :: self

    call sussingDestroy(self)
    return
  end subroutine sussingASCIIDestructor

  subroutine sussingASCIIOpen(self,fileName)
    !% Validate a {\normalfont \ttfamily sussing} ASCII format merger tree file.
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
    class           (mergerTreeImporterSussingASCII), intent(inout) :: self
    type            (varying_string                ), intent(in   ) :: fileName
    class           (cosmologyParametersClass      ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    type            (varying_string                )                :: message                 , baseDirectory        , &
         &                                                             simulationDefinitionFile, snapshotTimesFile
    character       (len=1024                      )                :: line                    , parameterName        , &
         &                                                             parameterValue
    character       (len=14                        )                :: valueString             , unitString
    type            (regEx                         )                :: parameterRegEx
    integer                                                         :: fileUnit                , snapshotNumber       , &
         &                                                             ioStat                  , snapshotFileCount    , &
         &                                                             i
    double precision                                                :: localLittleH0           , localOmegaMatter     , &
         &                                                             localOmegaDE            , cosmologicalParameter, &
         &                                                             expansionFactor         , redshift             , &
         &                                                             timeNormalized          , time

    ! Get the default cosmology.
    cosmologyParameters_ => cosmologyParameters()
    cosmologyFunctions_  => cosmologyFunctions ()
    ! Get cosmological parameters. We do this in advance to avoid HDF5 thread conflicts.
    localLittleH0   =cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH)
    localOmegaMatter=cosmologyParameters_%OmegaMatter    (                  )
    localOmegaDE    =cosmologyParameters_%OmegaDarkEnergy(                  )
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
    call allocateArray(self%snapshotTimes,[snapshotFileCount])
    open(newUnit=fileUnit,file=char(snapshotTimesFile),status='old',form='formatted',ioStat=ioStat)
    if (ioStat /= 0) call Galacticus_Error_Report('sussingASCIIOpen','can not open file "'//char(snapshotTimesFile)//'"')
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
             if (String_Strip(unitString) /= 'Mpc/h') call Galacticus_Error_Report('sussingASCIIOpen','box length should be reported in units of Mpc/h')
             self%boxLengthUnits=importerUnits(.true.,megaParsec,-1,0)
          end select
       end if
    end do
    close(fileUnit)
    call parameterRegEx%destroy()
    return
  end subroutine sussingASCIIOpen

  subroutine sussingASCIILoad(self,nodeSelfIndices,nodeIndexRanks,nodeDescendentLocations,nodeIncomplete,nodeCountTrees,nodeTreeIndices,treeIndicesAssigned,branchJumpCheckRequired,massUnits,lengthUnits,velocityUnits)
    !% Load a {\normalfont \ttfamily sussing} ASCII format merger tree data.
    use Galacticus_Display
    use Galacticus_Error
    use Kind_Numbers
    use String_Handling
    use Sort
    use File_Utilities
    use Memory_Management
    use Arrays_Search
    use Array_Utilities
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use, intrinsic :: ISO_C_Binding
    implicit none
    class    (mergerTreeImporterSussingASCII), intent(inout)                              :: self
    integer  (kind_int8                     ), intent(  out), dimension(:  ), allocatable :: nodeSelfIndices              , nodeTreeIndices
    integer  (c_size_t                      ), intent(  out), dimension(:  ), allocatable :: nodeIndexRanks               , nodeDescendentLocations
    logical                                  , intent(  out), dimension(:  ), allocatable :: nodeIncomplete
    integer  (kind=c_size_t                 ), intent(  out)                              :: nodeCountTrees
    logical                                  , intent(  out)                              :: treeIndicesAssigned          , branchJumpCheckRequired
    type     (importerUnits                 ), intent(  out)                              :: massUnits                    , lengthUnits            , &
         &                                                                                   velocityUnits
    integer  (kind=kind_int8                )               , dimension(:  ), allocatable :: nodesInSubvolume             , nodesTmp                    , &
         &                                                                                   hostsInSubvolume
    logical                                                 , dimension(:  ), allocatable :: nodeIncompleteTmp
    integer  (c_size_t                      )               , dimension(:  ), allocatable :: forestSnapshotHaloCount      , forestSnapshotHaloCountFirst, &
         &                                                                                   forestSnapshotHaloCountLast  , forestID
    integer  (c_size_t                      )               , dimension(:,:), allocatable :: forestSnapshotHaloCounts
    integer                                  , parameter                                  :: fileFormatVersionCurrent   =1
    logical                                                                               :: nodeIsActive                 , doBinaryConversion          , &
         &                                                                                   readBinary                   , mergerTreeFileIsBinary      , &
         &                                                                                   mergerTreeFileConvert        , processHalo
    integer                                                                               :: fileUnit                     , progenitorCount             , &
         &                                                                                   fileFormatVersion            , fileUnitOut                 , &
         &                                                                                   snapshotUnit                 , snapshotOutUnit             , &
         &                                                                                   ioStat
    character(len=1024                      )                                             :: line
    integer  (kind=c_size_t                 )                                             :: l                            , i                           , &
         &                                                                                   j                                                          , &
         &                                                                                   iNode                                                      , &
         &                                                                                   jNode                        , iCount                      , &
         &                                                                                   nodeCount                    , nodeCountSubvolume          , &
         &                                                                                   iProgenitor                  , jCount                      , &
         &                                                                                   jForest
    integer  (kind=kind_int8                )                                             :: nodeIndex
    type     (varying_string                )                                             :: message
    integer  (kind=kind_int8                )                                             :: ID                           , hostHalo                    , &
         &                                                                                   progenitorIndex 
    integer  (c_size_t                      )                                             :: forestCount                  , forestHaloCount             , &
         &                                                                                   forestFirst                  , forestLast                  , &
         &                                                                                   forestHaloCountLast          , forestHaloCountFirst
    integer                                                                               :: numSubStruct                 , npart                       , &
         &                                                                                   haloFormat
    double precision                                                                      :: Mvir                         , Xc                          , &
               &                                                                             Yc                           , Zc                          , &
               &                                                                             VXc                          , Vyc                         , &
               &                                                                             VZc                          , Rvir                        , &
               &                                                                             Rmax                         , r2                          , &
               &                                                                             mbp_offset                   , com_offset                  , &
               &                                                                             Vmax                         , v_esc                       , &
               &                                                                             sigV                         , lambda                      , &
               &                                                                             lambdaE                      , Lx                          , &
               &                                                                             Ly                           , Lz                          , &
               &                                                                             b                            , c                           , &
               &                                                                             Eax                          , Eay                         , &
               &                                                                             Eaz                          , Ebx                         , &
               &                                                                             Eby                          , Ebz                         , &
               &                                                                             Ecx                          , Ecy                         , &
               &                                                                             Ecz                          , ovdens                      , &
               &                                                                             fMhires                      , Ekin                        , &
               &                                                                             Epot                         , SurfP                       , &
               &                                                                             Phi0                         , cNFW                        , &
               &                                                                             nbins                        , FoFMass                     , &
               &                                                                             M_200Mean                    , M_200Crit                   , &
               &                                                                             M_TopHat                     , R_200Mean                   , &
               &                                                                             R_200Crit                    , R_TopHat                    , &
               &                                                                             HalfMassRadius               , sigV_200Mean                , &
               &                                                                             sigV_200Crit                 , sigV_TopHat                 , &
               &                                                                             Xcm                          , Ycm                         , &
               &                                                                             Zcm                          , Xgroup                      , &
               &                                                                             Ygroup                       , Zgroup

    ! Display counter.
    call Galacticus_Display_Indent ('Parsing "Sussing Merger Trees" format merger tree file',verbosityWorking)
    ! If a forest field is provided, scan it now to find the ranges to read from subsequent files.
    forestHaloCountFirst=0
    forestHaloCountLast =0
    if (mergerTreeImportSussingUseForestFile) then
       forestCount=Count_Lines_in_File(mergerTreeImportSussingForestFile,'#')
       call allocateArray(forestSnapshotHaloCount     ,                                                                        shape(self%snapshotFileName) )
       call allocateArray(forestSnapshotHaloCountFirst,                                                                        shape(self%snapshotFileName) )
       call allocateArray(forestSnapshotHaloCountLast ,                                                                        shape(self%snapshotFileName) )
       call allocateArray(forestSnapshotHaloCounts    ,[mergerTreeImportSussingForestLast-mergerTreeImportSussingForestFirst+1,size (self%snapshotFileName)])
       call allocateArray(forestID                    ,[mergerTreeImportSussingForestLast-mergerTreeImportSussingForestFirst+1                             ])
       forestFirst                 =mergerTreeImportSussingForestFirst
       forestLast                  =mergerTreeImportSussingForestLast
       forestSnapshotHaloCountFirst=0
       forestSnapshotHaloCountLast =0
       j                           =0
       if (forestLast < 0) forestLast=forestCount
       open(newUnit=fileUnit,file=char(mergerTreeImportSussingForestFile),status='old',form='formatted')
       read (fileUnit,'(a)') line
       do i=1,forestLast
          read (fileUnit,*) ID,forestHaloCount,forestSnapshotHaloCount
          forestHaloCountLast        =forestHaloCountLast        +forestHaloCount
          forestSnapshotHaloCountLast=forestSnapshotHaloCountLast+forestSnapshotHaloCount
          if (i < forestFirst) then
             forestHaloCountFirst         =forestHaloCountLast        +1
             forestSnapshotHaloCountFirst =forestSnapshotHaloCountLast+1
          else
             j                            =j                          +1
             forestID                (j  )=ID
             forestSnapshotHaloCounts(j,:)=forestSnapshotHaloCount
          end if
       end do
       close(fileUnit)
       call deallocateArray(forestSnapshotHaloCount)
       ! Reverse order of forest snapshots to match the order of snapshot files if necessary.
       if (mergerTreeImportSussingForestReverseSnapshotOrder) then
          forestSnapshotHaloCountFirst=Array_Reverse(forestSnapshotHaloCountFirst)
          forestSnapshotHaloCountLast =Array_Reverse(forestSnapshotHaloCountLast )
          do j=1,size(forestSnapshotHaloCounts,dim=1)
             forestSnapshotHaloCounts(j,:)=Array_Reverse(forestSnapshotHaloCounts(j,:))
          end do
       end if
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
       read     (fileUnit         ,ioStat=ioStat) fileFormatVersion
       read     (fileUnit         ,ioStat=ioStat) nodeCount
    else
       read     (fileUnit   ,*    ,ioStat=ioStat) fileFormatVersion
       read     (fileUnit   ,'(a)',ioStat=ioStat) line
       read     (fileUnit   ,*    ,ioStat=ioStat) nodeCount
       if (mergerTreeFileConvert) then
          write (fileUnitOut                ) fileFormatVersion
          write (fileUnitOut                ) nodeCount
       end if
    end if
    ! Validate file format version.
    if (fileFormatVersion /= fileFormatVersionCurrent) then
       message='incorrect file format version [found '
       message=message//fileFormatVersion//'; expected '//fileFormatVersionCurrent//';]'
       call Galacticus_Error_Report('sussingTreeIndicesRead',message)
    end if
    ! Allocate storage for list of nodes in subvolume.
    nodeCountSubVolume=int(dble(nodeCount)/dble(self%subvolumeCount)**3,kind=c_size_t)+1
    call allocateArray(nodesInSubvolume,[nodeCountSubVolume])
    call allocateArray(hostsInSubvolume,[nodeCountSubVolume])
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
             call sussingASCIIReadHalo(                         &
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
                call allocateArray(nodesInSubvolume,int(dble(shape(nodesTmp))*1.4d0)+1)
                nodesInSubvolume(1:size(nodesTmp))=nodesTmp
                call deallocateArray(nodesTmp)
                call Move_Alloc(hostsInSubvolume,nodesTmp)
                call allocateArray(hostsInSubvolume,int(dble(shape(nodesTmp))*1.4d0)+1)
                hostsInSubvolume(1:size(nodesTmp))=nodesTmp
                call deallocateArray(nodesTmp)
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
    call allocateArray(nodeSelfIndices,[nodeCountSubvolume])
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
       call allocateArray(nodeSelfIndices,[nodeCountTrees])
       nodeSelfIndices(1:nodeCountTrees)=nodesTmp(1:nodeCountTrees)
       call deallocateArray(nodesTmp)
    end if
    ! Allocate workspaces for merger trees.
    call allocateArray(nodeDescendentLocations,[nodeCountTrees])
    call allocateArray(nodeIncomplete         ,[nodeCountTrees])
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
    if (ioStat /= 0) call Galacticus_Error_Report('sussingTreeIndicesRead','failed to open merger tree file "'//char(self%mergerTreeFile)//'"')
    ! Read progenitor indices and make links.
    call Galacticus_Display_Message('Reading trees',verbosityWorking)
    if (mergerTreeFileIsBinary) then
       read (fileUnit  ,ioStat=ioStat) fileFormatVersion
       read (fileUnit  ,ioStat=ioStat) nodeCount
    else
       read (fileUnit,*,ioStat=ioStat) fileFormatVersion
       if (ioStat /= 0) call Galacticus_Error_Report('sussingTreeIndicesRead','failed to read merger tree file "'//char(self%mergerTreeFile)//'" header line 1')
       read (fileUnit,'(a)',ioStat=ioStat) line
       if (ioStat /= 0) call Galacticus_Error_Report('sussingTreeIndicesRead','failed to read merger tree file "'//char(self%mergerTreeFile)//'" header line 2')
       read (fileUnit,*,ioStat=ioStat) nodeCount
       if (ioStat /= 0) call Galacticus_Error_Report('sussingTreeIndicesRead','failed to read merger tree file "'//char(self%mergerTreeFile)//'" header line 3')
    end if
    i                      = 0
    iCount                 = 0
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
                            call allocateArray  (nodeSelfIndices        ,[nodeCountTrees+1])
                            nodeSelfIndices        (1:nodeCountTrees)=nodesTmp
                            call deallocateArray(nodesTmp                                  )
                            call Move_Alloc   (nodeDescendentLocations,nodesTmp          )
                            call allocateArray  (nodeDescendentLocations,[nodeCountTrees+1])
                            nodeDescendentLocations(1:nodeCountTrees)=nodesTmp
                            call deallocateArray(nodesTmp                                  )
                            call Move_Alloc   (nodeIncomplete         ,nodeIncompleteTmp )
                            call allocateArray  (nodeIncomplete         ,[nodeCountTrees+1])
                            nodeIncomplete         (1:nodeCountTrees)=nodeIncompleteTmp
                            call deallocateArray(nodeIncompleteTmp                         )
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
       iCount=iCount+1
       if (mergerTreeImportSussingUseForestFile .and. iCount == forestHaloCountLast) exit
    end do
    call deallocateArray(hostsInSubvolume)
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
    call allocateArray(nodeTreeIndices,[nodeCountTrees])
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
          read (snapshotUnit,'(a)',ioStat=ioStat) line
       end if
       iCount =0
       jCount =0
       jForest=1
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
             call sussingASCIIReadHalo(                         &
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
              ! Update forest ID.
             if (mergerTreeImportSussingUseForestFile) then
                jCount=jCount+1
                do while (jCount > forestSnapshotHaloCounts(jForest,i))
                   jForest=jForest+1
                   jCount =        1
                end do
             end if
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
                if (mergerTreeImportSussingUseForestFile) nodeTreeIndices(l)=forestID(jForest)
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
                case default
                   call Galacticus_Error_Report('sussingTreeIndicesRead','unrecognized mass option')
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
    call deallocateArray(nodesInSubvolume)
    ! Record whether tree indices were assigned.
    treeIndicesAssigned    =     mergerTreeImportSussingUseForestFile
    ! Record whether branch jump checks are required.
    branchJumpCheckRequired=.not.mergerTreeImportSussingUseForestFile
    ! Specify units.
    massUnits    =importerUnits(.true.,massSolar ,-1, 0)
    lengthUnits  =importerUnits(.true.,kiloParsec,-1,+1)
    velocityUnits=importerUnits(.true.,kilo      , 0, 0)
    ! Clean up display.
    call Galacticus_Display_Counter_Clear(       verbosityWorking)
    call Galacticus_Display_Unindent     ('done',verbosityWorking)
    return
  end subroutine sussingASCIILoad

  subroutine sussingASCIIReadHalo                 &
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
          call Galacticus_Error_Report('sussingASCIIReadHalo','unknown halo file format')
       end if
    end if
    return
  end subroutine sussingASCIIReadHalo
