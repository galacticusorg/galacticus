!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !![
  <mergerTreeImporter name="mergerTreeImporterSussingASCII">
   <description>
    A merger tree importer class for ``Sussing Merger Trees'' ASCII format merger tree files \citep{srisawat_sussing_2013},
    along with \gls{ahf} format halo catalogs. A descriptor file must be specified via the {\normalfont \ttfamily
    [mergerTreeReadFileName]} parameter. This descriptor file should have the following format:
    \begin{verbatim}
    simulation.txt
    MergerTree+AHF.txt
    snapidzred.txt
    AHF/62.5_dm_000.z50.000.AHF_halos
    AHF/62.5_dm_001.z30.000.AHF_halos
    AHF/62.5_dm_002.z19.916.AHF_halos
    .
    .
    .
    AHF/62.5_dm_061.z0.000.AHF_halos
    \end{verbatim}
    in which each line specifies a file to be read (by default path names are relative to the location of the descriptor
    file---fully-qualified path names can also be given).
  
    The first line identifies a file which specifies properties of the simulation. This file should look like:
    \begin{verbatim}
    WMAP7 cosmology:
    ----------------
    Omega0          =       0.272
    OmegaLambda0    =       0.728
    h               =       0.704
    
    simulation:
    -----------
    B               =       62.5 Mpc/h
    N               =       270^3 particles
    \end{verbatim}
    Currently only the cosmological parameter and box length are read from this file.
    
    The second line identifies the merger tree file which must be in the format specified by \cite{srisawat_sussing_2013}.
    
    The third line of the descriptor file specifies a snapshot file which should have the following format:
    \begin{verbatim}
    #    snapnum       a             z           t(t0)      t(year)
               0    0.0196080      49.9996   0.00354284  4.87485e+07
               1    0.0322580      30.0001   0.00747572  1.02864e+08
               2    0.0478110      19.9157    0.0134888  1.85602e+08
               3    0.0519650      18.2437    0.0152842  2.10306e+08
               4    0.0564190      16.7245    0.0172905  2.37912e+08
               5    0.0611880      15.3431    0.0195280  2.68700e+08
               6    0.0662870      14.0859    0.0220186  3.02969e+08
                 .
                 .
                 .
    \end{verbatim}
    This file must contain one line for each snapshot of the simulation, giving the snapshot number, expansion factor,
    redshift, fractional time (relative to present day), and age of the universe (in years).
    
    Subsequent lines identify the \gls{ahf} halo files for each snapshot (files can be listed in any order).
    
    Merger tree files of this type can be split into subvolumes before processing. This is useful if the file is too large to
    read into memory in one go. The number of subvolumes to use (in each of the three dimensions of the simulation cube) is
    specified by the {\normalfont \ttfamily [subvolumeCount]} parameter. The specific subvolume to process is specified by the
    {\normalfont \ttfamily [subvolumeIndex]} parameter, which should give the index (running from $0$ to {\normalfont \ttfamily
    [subvolumeCount]}$-1$) in each dimension (whitespace separated). To ensure that no halos are missed from trees near the
    edge of the subvolume, a buffer region around the subvolume is also read. The width of this buffer (in units of Mpc$/h$ to
    follow the format convention) is specified via the {\normalfont \ttfamily [subvolumeBuffer]} parameter.
   </description>
  </mergerTreeImporter>
  !!]
  type, extends(mergerTreeImporterSussing) :: mergerTreeImporterSussingASCII
     !!{
     A merger tree importer class for ``Sussing Merger Trees'' ASCII format merger tree files \citep{srisawat_sussing_2013}.
     !!}
     private
     logical                 :: convertToBinary           , binaryFormatOld, &
          &                     forestReverseSnapshotOrder, useForestFile
     integer                 :: forestFirst               , forestLast
     type   (varying_string) :: forestFile
   contains
     !![
     <methods>
       <method method="initialize" description="Initialize the object after construction."/>
     </methods>
     !!]
     procedure :: initialize => sussingASCIIInitialize
     procedure :: open       => sussingASCIIOpen
     procedure :: load       => sussingASCIILoad
  end type mergerTreeImporterSussingASCII

  interface mergerTreeImporterSussingASCII
     !!{
     Constructors for the \refClass{mergerTreeImporterSussingASCII} ASCII format merger tree importer class.
     !!}
     module procedure sussingASCIIConstructorParameters
     module procedure sussingASCIIConstructorInternal
  end interface mergerTreeImporterSussingASCII

  ! Enumeration of file formats.
  !![
  <enumeration>
   <name>sussingHaloFormat</name>
   <description>Halo file formats.</description>
   <entry label="old"/>
   <entry label="new"/>
   <entry label="all"/>
  </enumeration>
  !!]

contains

  function sussingASCIIConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeImporterSussingASCII} ASCII format \citep{srisawat_sussing_2013} merger tree importer which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(mergerTreeImporterSussingASCII)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    !![
    <inputParameter>
      <name>convertToBinary</name>
      <variable>self%convertToBinary</variable>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether halo and tree files in the ``Sussing'' format should be converted to binary the first time they are read and stored to file. This allows rapid re-reading in future.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>binaryFormatOld</name>
      <variable>self%binaryFormatOld</variable>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether the old binary format is to be used (for reading only).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>forestFile</name>
      <variable>self%forestFile</variable>
      <defaultValue>var_str('none')</defaultValue>
      <description>Name of file containing data on number of halos in each forest.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>forestFirst</name>
      <variable>self%forestFirst</variable>
      <defaultValue>1</defaultValue>
      <description>Index of first forest to include.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>forestLast</name>
      <variable>self%forestLast</variable>
      <defaultValue>-1</defaultValue>
      <description>Index of last forest to include.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>forestReverseSnapshotOrder</name>
      <variable>self%forestReverseSnapshotOrder</variable>
      <defaultValue>.false.</defaultValue>
      <description>If true, the order of forest snapshots will be reversed after being read. This may be necessary to cause them to match the order of snapshot files.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self%mergerTreeImporterSussing=mergerTreeImporterSussing(parameters)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    call self%initialize()
    return
  end function sussingASCIIConstructorParameters

  function sussingASCIIConstructorInternal(fatalMismatches,fatalNonTreeNode,subvolumeCount,subvolumeBuffer,subvolumeIndex,badValue,badValueTest,treeSampleRate,massOption,convertToBinary,binaryFormatOld,forestFile,forestFirst,forestLast,forestReverseSnapshotOrder,cosmologyParameters_,cosmologyFunctions_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeImporterSussingASCII} ASCII format \citep{srisawat_sussing_2013} merger tree importer.
    !!}
    implicit none
    type            (mergerTreeImporterSussingASCII    )                              :: self
    class           (cosmologyParametersClass          ), intent(in   ), target       :: cosmologyParameters_
    class           (cosmologyFunctionsClass           ), intent(in   ), target       :: cosmologyFunctions_
    class           (randomNumberGeneratorClass        ), intent(in   ), target       :: randomNumberGenerator_
    integer                                             , intent(in   ), dimension(3) :: subvolumeIndex
    logical                                             , intent(in   )               :: fatalMismatches           , fatalNonTreeNode, &
         &                                                                               convertToBinary           , binaryFormatOld , &
         &                                                                               forestReverseSnapshotOrder
    integer                                             , intent(in   )               :: subvolumeCount            , forestFirst     , &
         &                                                                               forestLast
    type            (enumerationSussingBadValueTestType), intent(in   )               :: badValueTest
    type            (enumerationSussingMassOptionType  ), intent(in   )               :: massOption
    double precision                                    , intent(in   )               :: subvolumeBuffer           , badValue        , &
         &                                                                               treeSampleRate
    type            (varying_string                    ), intent(in   )               :: forestFile
    !![
    <constructorAssign variables="convertToBinary,binaryFormatOld,forestFile,forestFirst,forestLast,forestReverseSnapshotOrder,*cosmologyParameters_,*randomNumberGenerator_"/>
    !!]

    self%mergerTreeImporterSussing=mergerTreeImporterSussing(fatalMismatches,fatalNonTreeNode,subvolumeCount,subvolumeBuffer,subvolumeIndex,badValue,badValueTest,treeSampleRate,massOption,cosmologyParameters_,cosmologyFunctions_,randomNumberGenerator_)
    call self%initialize()
    return
  end function sussingASCIIConstructorInternal

  subroutine sussingASCIIInitialize(self)
    !!{
    Initialize the object after construction.
    !!}
    implicit none
    class(mergerTreeImporterSussingASCII), intent(inout) :: self

    self%useForestFile=self%forestFile /= "none"
    return
  end subroutine sussingASCIIInitialize

  subroutine sussingASCIIOpen(self,fileName)
    !!{
    Validate a {\normalfont \ttfamily sussing} ASCII format merger tree file.
    !!}
    use :: Cosmology_Parameters            , only : hubbleUnitsLittleH
    use :: Display                         , only : displayMessage         , verbosityLevelWarn
    use :: File_Utilities                  , only : Count_Lines_in_File
    use :: Error                           , only : Error_Report
    use :: Numerical_Comparison            , only : Values_Differ
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Regular_Expressions             , only : regEx
    use :: String_Handling                 , only : String_Strip           , operator(//)
    implicit none
    class           (mergerTreeImporterSussingASCII), intent(inout) :: self
    type            (varying_string                ), intent(in   ) :: fileName
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

    ! Get cosmological parameters.
    localLittleH0   =self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH)
    localOmegaMatter=self%cosmologyParameters_%OmegaMatter    (                  )
    localOmegaDE    =self%cosmologyParameters_%OmegaDarkEnergy(                  )
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
    allocate(self%snapshotTimes(snapshotFileCount))
    open(newUnit=fileUnit,file=char(snapshotTimesFile),status='old',form='formatted',ioStat=ioStat)
    if (ioStat /= 0) call Error_Report('can not open file "'//char(snapshotTimesFile)//'"'//{introspection:location})
    read (fileUnit,*)
    do i=1,snapshotFileCount
       read (fileUnit,*) snapshotNumber,expansionFactor,redshift,timeNormalized,time
       self%snapshotTimes(i)=                                                 &
            & self%cosmologyFunctions_ %cosmicTime                 (          &
            &  self%cosmologyFunctions_%expansionFactorFromRedshift (         &
            &                                                        redshift &
            &                                                       )         &
            &                                                      )
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
                   call Error_Report(message//{introspection:location})
                else
                   call displayMessage(message,verbosityLevelWarn)
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
                   call Error_Report(message//{introspection:location})
                else
                   call displayMessage(message,verbosityLevelWarn)
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
                   call Error_Report(message//{introspection:location})
                else
                   call displayMessage(message,verbosityLevelWarn)
                end if
             end if
          case ('B')
             read (parameterValue,*) self%boxLength
             unitString=String_Strip(parameterValue(index(parameterValue,' '):len(parameterValue)-index(parameterValue,' ')+1))
             if (String_Strip(unitString) /= 'Mpc/h') call Error_Report('box length should be reported in units of Mpc/h'//{introspection:location})
             self%boxLengthUnits=importerUnits(.true.,megaParsec,-1,0)
          end select
       end if
    end do
    close(fileUnit)
    return
  end subroutine sussingASCIIOpen

  subroutine sussingASCIILoad(self,nodeSelfIndices,nodeIndexRanks,nodeDescendantLocations,nodeIncomplete,nodeCountTrees,nodeTreeIndices,treeIndicesAssigned,branchJumpCheckRequired,massUnits,lengthUnits,velocityUnits)
    !!{
    Load a {\normalfont \ttfamily sussing} ASCII format merger tree data.
    !!}
    use            :: Array_Utilities                 , only : Array_Reverse
    use            :: Arrays_Search                   , only : searchArray            , searchIndexed
    use            :: Display                         , only : displayCounter         , displayCounterClear  , displayIndent, displayMessage, &
          &                                                    displayUnindent        , verbosityLevelWorking
    use            :: File_Utilities                  , only : Count_Lines_in_File    , File_Exists
    use            :: Error                           , only : Error_Report
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Kind_Numbers                    , only : kind_int8
    use            :: Numerical_Constants_Astronomical, only : kiloParsec             , massSolar
    use            :: Numerical_Constants_Prefixes    , only : kilo
    use            :: Sorting                         , only : sort                   , sortIndex
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class           (mergerTreeImporterSussingASCII), intent(inout)                              :: self
    integer         (kind_int8                     ), intent(  out), dimension(:  ), allocatable :: nodeSelfIndices              , nodeTreeIndices
    integer         (c_size_t                      ), intent(  out), dimension(:  ), allocatable :: nodeIndexRanks               , nodeDescendantLocations
    logical                                         , intent(  out), dimension(:  ), allocatable :: nodeIncomplete
    integer         (kind=c_size_t                 ), intent(  out)                              :: nodeCountTrees
    logical                                         , intent(  out)                              :: treeIndicesAssigned          , branchJumpCheckRequired
    type            (importerUnits                 ), intent(  out)                              :: massUnits                    , lengthUnits            , &
         &                                                                                          velocityUnits
    integer         (kind=kind_int8                )               , dimension(:  ), allocatable :: nodesInSubvolume             , nodesTmp                    , &
         &                                                                                   hostsInSubvolume
    logical                                                        , dimension(:  ), allocatable :: nodeIncompleteTmp
    integer         (c_size_t                      )               , dimension(:  ), allocatable :: forestSnapshotHaloCount      , forestSnapshotHaloCountFirst, &
         &                                                                                          forestSnapshotHaloCountLast  , forestID
    integer         (c_size_t                      )               , dimension(:,:), allocatable :: forestSnapshotHaloCounts
    integer                                         , parameter                                  :: fileFormatVersionCurrent   =1
    logical                                                                                      :: nodeIsActive                 , doBinaryConversion          , &
         &                                                                                          readBinary                   , mergerTreeFileIsBinary      , &
         &                                                                                          mergerTreeFileConvert        , processHalo
    integer                                                                                      :: fileUnit                     , progenitorCount             , &
         &                                                                                          fileFormatVersion            , fileUnitOut                 , &
         &                                                                                          snapshotUnit                 , snapshotOutUnit             , &
         &                                                                                          ioStat
    character       (len=1024                      )                                             :: line
    integer         (kind=c_size_t                 )                                             :: l                            , i                           , &
         &                                                                                          j                                                          , &
         &                                                                                          iNode                                                      , &
         &                                                                                          jNode                        , iCount                      , &
         &                                                                                          nodeCount                    , nodeCountSubvolume          , &
         &                                                                                          iProgenitor                  , jCount                      , &
         &                                                                                          jForest
    integer         (kind=kind_int8                )                                             :: nodeIndex
    type            (varying_string                )                                             :: message
    integer         (kind=kind_int8                )                                             :: ID                           , hostHalo                    , &
         &                                                                                          progenitorIndex
    integer         (c_size_t                      )                                             :: forestCount                  , forestHaloCount             , &
         &                                                                                          forestFirst                  , forestLast                  , &
         &                                                                                          forestHaloCountLast          , forestHaloCountFirst
    integer                                                                                      :: numSubStruct                 , npart
    type            (enumerationSussingHaloFormatType)                                           :: haloFormat
    double precision                                                                             :: Mvir                         , Xc                          , &
               &                                                                                    Yc                           , Zc                          , &
               &                                                                                    VXc                          , Vyc                         , &
               &                                                                                    VZc                          , Rvir                        , &
               &                                                                                    Rmax                         , r2                          , &
               &                                                                                    mbp_offset                   , com_offset                  , &
               &                                                                                    Vmax                         , v_esc                       , &
               &                                                                                    sigV                         , lambda                      , &
               &                                                                                    lambdaE                      , Lx                          , &
               &                                                                                    Ly                           , Lz                          , &
               &                                                                                    b                            , c                           , &
               &                                                                                    Eax                          , Eay                         , &
               &                                                                                    Eaz                          , Ebx                         , &
               &                                                                                    Eby                          , Ebz                         , &
               &                                                                                    Ecx                          , Ecy                         , &
               &                                                                                    Ecz                          , ovdens                      , &
               &                                                                                    fMhires                      , Ekin                        , &
               &                                                                                    Epot                         , SurfP                       , &
               &                                                                                    Phi0                         , cNFW                        , &
               &                                                                                    nbins                        , FoFMass                     , &
               &                                                                                    M_200Mean                    , M_200Crit                   , &
               &                                                                                    M_TopHat                     , R_200Mean                   , &
               &                                                                                    R_200Crit                    , R_TopHat                    , &
               &                                                                                    HalfMassRadius               , sigV_200Mean                , &
               &                                                                                    sigV_200Crit                 , sigV_TopHat                 , &
               &                                                                                    Xcm                          , Ycm                         , &
               &                                                                                    Zcm                          , Xgroup                      , &
               &                                                                                    Ygroup                       , Zgroup

    ! Display counter.
    call displayIndent ('Parsing "Sussing Merger Trees" format merger tree file',verbosityLevelWorking)
    ! If a forest field is provided, scan it now to find the ranges to read from subsequent files.
    forestHaloCountFirst=0
    forestHaloCountLast =0
    if (self%useForestFile) then
       forestCount=Count_Lines_in_File(self%forestFile,'#')
       allocate(forestSnapshotHaloCount     (                                   size(self%snapshotFileName)))
       allocate(forestSnapshotHaloCountFirst(                                   size(self%snapshotFileName)))
       allocate(forestSnapshotHaloCountLast (                                   size(self%snapshotFileName)))
       allocate(forestSnapshotHaloCounts    (self%forestLast-self%forestFirst+1,size(self%snapshotFileName)))
       allocate(forestID                    (self%forestLast-self%forestFirst+1                            ))
       forestFirst                 =self%forestFirst
       forestLast                  =self%forestLast
       forestSnapshotHaloCountFirst=0
       forestSnapshotHaloCountLast =0
       j                           =0
       if (forestLast < 0) forestLast=forestCount
       open(newUnit=fileUnit,file=char(self%forestFile),status='old',form='formatted')
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
       deallocate(forestSnapshotHaloCount)
       ! Reverse order of forest snapshots to match the order of snapshot files if necessary.
       if (self%forestReverseSnapshotOrder) then
          forestSnapshotHaloCountFirst=Array_Reverse(forestSnapshotHaloCountFirst)
          forestSnapshotHaloCountLast =Array_Reverse(forestSnapshotHaloCountLast )
          do j=1,size(forestSnapshotHaloCounts,dim=1)
             forestSnapshotHaloCounts(j,:)=Array_Reverse(forestSnapshotHaloCounts(j,:))
          end do
       end if
    else
       allocate(forestID                (0  ))
       allocate(forestSnapshotHaloCounts(0,0))
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
    call displayMessage('Reading header',verbosityLevelWorking)
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
       call Error_Report(message//{introspection:location})
    end if
    ! Allocate storage for list of nodes in subvolume.
    nodeCountSubVolume=int(dble(nodeCount)/dble(self%subvolumeCount)**3,kind=c_size_t)+1
    allocate(nodesInSubvolume(nodeCountSubVolume))
    allocate(hostsInSubvolume(nodeCountSubVolume))
    ! Read snapshot halo catalogs.
    call displayIndent('Finding halos in subvolume from AHF format snapshot halo catalogs',verbosityLevelWorking)
    j                 =0
    nodeCountSubVolume=0
    do i=1,size(self%snapshotFileName)
       call displayMessage(self%snapshotFileName(i),verbosityLevelWorking)
       doBinaryConversion=.false.
       readBinary        =.false.
       if (File_Exists(self%snapshotFileName(i)//".bin")) then
          ! A binary version of this file exists, use it.
          open   (newUnit=snapshotUnit   ,file=char(self%snapshotFileName(i)//".bin"),status='old'    ,form='unformatted',ioStat=ioStat)
          readBinary=.true.
       else
          ! No binary version of this file exists, use the ASCII version.
          open   (newUnit=snapshotUnit   ,file=char(self%snapshotFileName(i)        ),status='old'    ,form='formatted'  ,ioStat=ioStat)
          if (self%convertToBinary.and..not.self%useForestFile) then
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
             call Error_Report('unrecognized format for halo files'//{introspection:location})
          end if
       end if
       iCount=0
       do while (ioStat == 0)
          ! Increment count of number of halos read.
          iCount=iCount+1
          processHalo=                                             &
               &       (                                           &
               &            self%useForestFile   &
               &        .and.                                      &
               &         iCount >= forestSnapshotHaloCountFirst(i) &
               &        .and.                                      &
               &         iCount <= forestSnapshotHaloCountLast (i) &
               &       )                                           &
               &      .or.                                         &
               &       .not.self%useForestFile
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
                call move_alloc(nodesInSubvolume,              nodesTmp            )
                allocate       (nodesInSubvolume(int(dble(size(nodesTmp))*1.4d0)+1))
                nodesInSubvolume(1:size(nodesTmp))=nodesTmp
                deallocate(nodesTmp)
                call move_alloc(hostsInSubvolume,              nodesTmp            )
                allocate       (hostsInSubvolume(int(dble(size(nodesTmp))*1.4d0)+1))
                hostsInSubvolume(1:size(nodesTmp))=nodesTmp
                deallocate(nodesTmp)
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
          call displayCounter(int(100.0d0*dble(j)/dble(nodeCount)),j == 1,verbosityLevelWorking)
          ! If all required forests are processed, exit.
          if (self%useForestFile .and. iCount == forestSnapshotHaloCountLast(i)) exit
       end do
       close                        (snapshotUnit   )
       if (doBinaryConversion) close(snapshotOutUnit)
    end do
    call displayCounterClear(verbosityLevelWorking)
    call sort(nodesInSubvolume(1:nodeCountSubvolume),hostsInSubvolume(1:nodeCountSubvolume))
    message='Found '
    message=message//nodeCountSubvolume//' nodes in subvolume [from '//nodeCount//' total nodes]'
    call displayMessage(message,verbosityLevelWorking)
    ! Allocate workspaces for merger trees.
    allocate(nodeSelfIndices(nodeCountSubvolume))
    ! Read node indices.
    call displayMessage("Reading node indices",verbosityLevelWorking)
    i     =0
    iCount=0
    ioStat=0
    call displayCounter(0,.true.,verbosityLevelWorking)
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
       iNode=searchArray(nodesInSubvolume(1:nodeCountSubvolume),nodeIndex)
       if (iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == nodeIndex) then
          ! This line represents a node in the tree.
          i                 =i+1
          nodeSelfIndices(i)=nodeIndex
       end if
       call displayCounter(int(50.0d0*dble(i)/dble(nodeCountSubvolume)),.false.,verbosityLevelWorking)
       if (i      == nodeCount          ) exit
       iCount=iCount+1
       if (self%useForestFile .and. iCount == forestHaloCountLast) exit
    end do
    close(fileUnit)
    if (mergerTreeFileConvert) close(fileUnitOut)
    ! Some halos in our subvolume might not be part of any tree. Adjust number of halos accordingly.
    nodeCountTrees=i
    if (nodeCountTrees < nodeCountSubvolume) then
       message='Found '
       message=message//nodeCountTrees//' nodes in subvolume trees [from '//nodeCountSubvolume//' total nodes in subvolume]'
       call displayMessage(message,verbosityLevelWorking)
       call Move_Alloc(nodeSelfIndices,nodesTmp)
       allocate(nodeSelfIndices(nodeCountTrees))
       nodeSelfIndices(1:nodeCountTrees)=nodesTmp(1:nodeCountTrees)
       deallocate(nodesTmp)
    end if
    ! Allocate workspaces for merger trees.
    allocate(nodeDescendantLocations(nodeCountTrees))
    allocate(nodeIncomplete         (nodeCountTrees))
    ! Get a sorted index into the list of nodes.
    call displayMessage('Building node index',verbosityLevelWorking)
    nodeIndexRanks=sortIndex(nodeSelfIndices)
    ! Re-open the merger tree file.
    mergerTreeFileIsBinary=File_Exists(char(self%mergerTreeFile//".bin"))
    if (mergerTreeFileIsBinary) then
       open(newUnit=fileUnit,file=char(self%mergerTreeFile//".bin"),status='old',form='unformatted',ioStat=ioStat)
    else
       open(newUnit=fileUnit,file=char(self%mergerTreeFile        ),status='old',form='formatted'  ,ioStat=ioStat)
    end if
    if (ioStat /= 0) call Error_Report('failed to open merger tree file "'//char(self%mergerTreeFile)//'"'//{introspection:location})
    ! Read progenitor indices and make links.
    call displayMessage('Reading trees',verbosityLevelWorking)
    if (mergerTreeFileIsBinary) then
       read (fileUnit  ,ioStat=ioStat) fileFormatVersion
       read (fileUnit  ,ioStat=ioStat) nodeCount
    else
       read (fileUnit,*,ioStat=ioStat) fileFormatVersion
       if (ioStat /= 0) call Error_Report('failed to read merger tree file "'//char(self%mergerTreeFile)//'" header line 1'//{introspection:location})
       read (fileUnit,'(a)',ioStat=ioStat) line
       if (ioStat /= 0) call Error_Report('failed to read merger tree file "'//char(self%mergerTreeFile)//'" header line 2'//{introspection:location})
       read (fileUnit,*,ioStat=ioStat) nodeCount
       if (ioStat /= 0) call Error_Report('failed to read merger tree file "'//char(self%mergerTreeFile)//'" header line 3'//{introspection:location})
    end if
    i                      = 0
    iCount                 = 0
    nodeDescendantLocations=-1
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
       iNode=searchArray(nodesInSubvolume(1:nodeCountSubvolume),nodeIndex)
       nodeIsActive=(iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == nodeIndex)
       if (nodeIsActive) then
          ! This line represents a node in the tree.
          i                 =i+1
          nodeSelfIndices(i)=nodeIndex
          call displayCounter(int(50.0d0+50.0d0*dble(i)/dble(nodeCountTrees)),.false.,verbosityLevelWorking)
       end if
       do j=1,progenitorCount
          if (mergerTreeFileIsBinary) then
             read (fileUnit  ,ioStat=ioStat) nodeIndex
          else
             read (fileUnit,*,ioStat=ioStat) nodeIndex
          end if
          if (nodeIsActive) then
             ! This line represents a progenitor. Locate the progenitor in the list of halos.
             iProgenitor=searchIndexed(nodeSelfIndices,nodeIndexRanks,nodeIndex)
             ! Does this progenitor exist within our subvolume?
             if (iProgenitor <= 0 .or. iProgenitor > nodeCountTrees .or. nodeSelfIndices(nodeIndexRanks(iProgenitor)) /= nodeIndex) then
                nodeIncomplete(i)=.true.
             else
                if (nodeDescendantLocations(nodeIndexRanks(iProgenitor)) /= -1) then
                   message="multiple descendant trees not allowed"
                   message=message//char(10)//" first descendant: "//nodeSelfIndices(nodeDescendantLocations(nodeIndexRanks(iProgenitor)))
                   message=message//char(10)//"   new descendant: "//nodeSelfIndices(i)
                   message=message//char(10)//" progenitor index: "//nodeIndex
                   call Error_Report(message//{introspection:location})
                end if
                nodeDescendantLocations(nodeIndexRanks(iProgenitor))=i
                ! Find the progenitor node in the list of halos in the subvolume.
                iNode=searchArray(nodesInSubvolume(1:nodeCountSubvolume),nodeIndex)
                if (iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == nodeIndex) then
                   ! Find hosted halos.
                   hostHalo=hostsInSubvolume(iNode)
                   if (hostHalo /= nodeIndex) then
                      ! Check if the host halo is in the subvolume.
                      jNode=searchArray(nodesInSubvolume(1:nodeCountSubvolume),hostHalo)
                      if (jNode > 0 .and. jNode <= nodeCountSubvolume .and. nodesInSubvolume(jNode) == hostHalo) then
                         ! Check if the host halo is in the trees.
                         jNode=searchIndexed(nodeSelfIndices,nodeIndexRanks,hostHalo)
                         if (.not.(jNode > 0 .and. jNode <= nodeCountTrees .and. nodeSelfIndices(nodeIndexRanks(jNode)) == hostHalo)) then
                            ! Host halo is in subvolume, but not in trees. Add it to the trees now.
                            message='host halo ['
                            message=message//hostHalo//'] in subvolume but not in trees - adding it now'
                            call displayMessage(message,verbosityLevelWorking)
                            ! Expand arrays.
                            call Move_Alloc   (nodeSelfIndices        ,nodesTmp          )
                            allocate(nodeSelfIndices        (nodeCountTrees+1))
                            nodeSelfIndices        (1:nodeCountTrees)=nodesTmp
                            deallocate(nodesTmp                                  )
                            call Move_Alloc   (nodeDescendantLocations,nodesTmp          )
                            allocate(nodeDescendantLocations(nodeCountTrees+1))
                            nodeDescendantLocations(1:nodeCountTrees)=nodesTmp
                            deallocate(nodesTmp                                  )
                            call Move_Alloc   (nodeIncomplete         ,nodeIncompleteTmp )
                            allocate(nodeIncomplete         (nodeCountTrees+1))
                            nodeIncomplete         (1:nodeCountTrees)=nodeIncompleteTmp
                            deallocate(nodeIncompleteTmp                         )
                            ! Increment the number of halos in trees.
                            nodeCountTrees=nodeCountTrees+1
                            ! Insert the new halo, assigning the same descendant as its hosted halo.
                            nodeSelfIndices        (nodeCountTrees)=hostHalo
                            nodeDescendantLocations(nodeCountTrees)=i
                            nodeIncomplete         (nodeCountTrees)=.false.
                            ! Recompute the sort index into the node self indices.
                            deallocate(nodeIndexRanks)
                            nodeIndexRanks=sortIndex(nodeSelfIndices)
                         end if
                      end if
                   end if
                else
                   message='can not find halo ['
                   message=message//nodeIndex//'] in subvolume'
                   call Error_Report(message//{introspection:location})
                end if
             end if
          end if
       end do
       iCount=iCount+1
       if (self%useForestFile .and. iCount == forestHaloCountLast) exit
    end do
    deallocate(hostsInSubvolume)
    ! Close the merger tree file.
    close(fileUnit)
    ! Clear counter.
    call displayCounterClear(       verbosityLevelWorking)
    call displayUnindent     ('done',verbosityLevelWorking)
    ! Transfer tree structure to nodes array.
    allocate(self%nodes(nodeCountTrees))
    do i=1,nodeCountTrees
       self   %nodes(i)%nodeIndex      =nodeSelfIndices(                        i )
       if (nodeDescendantLocations(i) >= 0) then
          self%nodes(i)%descendantIndex=nodeSelfIndices(nodeDescendantLocations(i))
       else
          self%nodes(i)%descendantIndex=-1
       end if
    end do
    ! Read snapshot halo catalogs.
    call displayIndent('Parsing AHF format snapshot halo catalogs',verbosityLevelWorking)
    allocate(nodeTreeIndices(nodeCountTrees))
    j=0
    do i=1,size(self%snapshotFileName)
       call displayMessage(self%snapshotFileName(i),verbosityLevelWorking)
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
               &            self%useForestFile   &
               &        .and.                                      &
               &         iCount >= forestSnapshotHaloCountFirst(i) &
               &        .and.                                      &
               &         iCount <= forestSnapshotHaloCountLast (i) &
               &       )                                           &
               &      .or.                                         &
               &       .not.self%useForestFile
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
             if (self%useForestFile) then
                jCount=jCount+1
                do while (jCount > forestSnapshotHaloCounts(jForest,i))
                   jForest=jForest+1
                   jCount =        1
                end do
             end if
             ! Locate this node in the list of nodes in our subvolume.
             iNode=searchArray(nodesInSubvolume(1:nodeCountSubvolume),ID)
             if (iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == ID) then
                ! Locate this node in the node list.
                l=searchIndexed(nodeSelfIndices,nodeIndexRanks,ID)
                l=nodeIndexRanks(l)
                if (ID /= nodeSelfIndices(l)) then
                   if (self%fatalNonTreeNode) then
                      ! Node cannot be found.
                      message="node indexing failure"
                      message=message//char(10)//"     node index: "//ID
                      message=message//char(10)//"    found index: "//nodeSelfIndices(l)
                      message=message//char(10)//" found location: "//l
                      call Error_Report(message//{introspection:location})
                   else
                      ! Just skip this node.
                      cycle
                   end if
                end if
                ! Store properties to node array.
                if (self%useForestFile) nodeTreeIndices(l)=forestID(jForest)
                if (hostHalo <= 0) then
                   self%nodes(l)%hostIndex         =ID
                else
                   self%nodes(l)%hostIndex         =hostHalo
                   ! Check that the host halo is in the subvolume.
                   iNode=searchArray(nodesInSubvolume(1:nodeCountSubvolume),hostHalo)
                   if (.not.(iNode > 0 .and. iNode <= nodeCountSubvolume .and. nodesInSubvolume(iNode) == hostHalo)) nodeIncomplete(l)=.true.
                end if
                self   %nodes(l)%particleCount     =npart
                ! Select a mass to use.
                select case (self%massOption%ID)
                case (sussingMassOptionDefault%ID)
                   self%nodes(l)%nodeMass          =Mvir
                case (sussingMassOptionFoF    %ID)
                   self%nodes(l)%nodeMass          =FoFMass
                case (sussingMassOption200Mean%ID)
                   self%nodes(l)%nodeMass          =M_200Mean
                case (sussingMassOption200Crit%ID)
                   self%nodes(l)%nodeMass          =M_200Crit
                case (sussingMassOptionTopHat %ID)
                   self%nodes(l)%nodeMass          =M_TopHat
                case default
                   call Error_Report('unrecognized mass option'//{introspection:location})
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
                call displayCounter(int(100.0d0*dble(j)/dble(nodeCountTrees)),j == 1,verbosityLevelWorking)
             end if
          end if
          ! If all required forests are processed, exit.
          if (self%useForestFile .and. iCount == forestSnapshotHaloCountLast(i)) exit
       end do
       close(snapshotUnit)
    end do
    deallocate(nodesInSubvolume)
    ! Record whether tree indices were assigned.
    treeIndicesAssigned    =     self%useForestFile
    ! Record whether branch jump checks are required.
    branchJumpCheckRequired=.not.self%useForestFile
    ! Specify units.
    massUnits    =importerUnits(.true.,massSolar ,-1, 0)
    lengthUnits  =importerUnits(.true.,kiloParsec,-1,+1)
    velocityUnits=importerUnits(.true.,kilo      , 0, 0)
    ! Clean up display.
    call displayCounterClear(       verbosityLevelWorking)
    call displayUnindent     ('done',verbosityLevelWorking)
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
    !!{
    Read an ASCII halo definition.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (enumerationSussingHaloFormatType), intent(in   )           :: haloFormat
    integer                                           , intent(in   )           :: snapshotUnit
    double precision                                                            :: Mvir          , Xc          , &
         &                                                                         Yc            , Zc          , &
         &                                                                         VXc           , Vyc         , &
         &                                                                         VZc           , Rvir        , &
         &                                                                         Rmax          , r2          , &
         &                                                                         mbp_offset    , com_offset  , &
         &                                                                         Vmax          , v_esc       , &
         &                                                                         sigV          , lambda      , &
         &                                                                         lambdaE       , Lx          , &
         &                                                                         Ly            , Lz          , &
         &                                                                         b             , c           , &
         &                                                                         Eax           , Eay         , &
         &                                                                         Eaz           , Ebx         , &
         &                                                                         Eby           , Ebz         , &
         &                                                                         Ecx           , Ecy         , &
         &                                                                         Ecz           , ovdens      , &
         &                                                                         fMhires       , Ekin        , &
         &                                                                         Epot          , SurfP       , &
         &                                                                         Phi0          , cNFW        , &
         &                                                                         nbins         , FoFMass     , &
         &                                                                         M_200Mean     , M_200Crit   , &
         &                                                                         M_TopHat      , R_200Mean   , &
         &                                                                         R_200Crit     , R_TopHat    , &
         &                                                                         HalfMassRadius, sigV_200Mean, &
         &                                                                         sigV_200Crit  , sigV_TopHat , &
         &                                                                         Xcm           , Ycm         , &
         &                                                                         Zcm           , Xgroup      , &
         &                                                                         Ygroup        , Zgroup
    integer         (kind=kind_int8                  ), intent(  out)           :: ID            , hostHalo
    integer                                           , intent(  out)           :: numSubStruct  , npart       , &
         &                                                                         ioStat
    logical                                           , intent(in   ), optional :: quickRead

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
          call Error_Report('unknown halo file format'//{introspection:location})
       end if
    end if
    return
  end subroutine sussingASCIIReadHalo
