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
Contains a module which provides various interfaces to the \gls{class} code.
!!}

module Interfaces_CLASS
  !!{
  Provides various interfaces to the \gls{class} code.
  !!}
  use :: File_Utilities, only : lockDescriptor
  private
  public :: Interface_CLASS_Initialize, Interface_CLASS_Transfer_Function, Interface_CLASS_Normalization

  ! Current file format version for transfer function files. Note that this file format matches that used by the "file" transfer
  ! function class.
  integer                         , parameter :: classFormatVersionCurrent       =     2

  ! Default maximum wavenumber to tabulate.
  double precision                , parameter :: classLogWavenumberMaximumDefault=log(10.0d0)

  !![
  <enumeration>
   <name>classSpecies</name>
   <description>Particle species in CLASS.</description>
   <visibility>public</visibility>
   <indexing>1</indexing>
   <entry label="photons"   />
   <entry label="darkMatter"/>
   <entry label="baryons"   />
  </enumeration>
  !!]

  ! Generate a source digest.
  !![
  <sourceDigest name="classSourceDigest"/>
  !!]

contains

  subroutine Interface_CLASS_Initialize(classPath,classVersion,static)
    !!{
    Initialize the interface with CLASS, including downloading and compiling CLASS if necessary.
    !!}
    use :: Dependencies      , only : dependencyVersion
    use :: Display           , only : displayMessage   , verbosityLevelWorking
    use :: File_Utilities    , only : Directory_Make   , File_Exists          , File_Lock   , File_Unlock, &
          &                           lockDescriptor
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath        , pathTypeDataDynamic
    use :: ISO_Varying_String, only : assignment(=)    , char                 , operator(//), replace    , &
          &                           varying_string
    use :: String_Handling   , only : stringSubstitute
    use :: System_Command    , only : System_Command_Do
    use :: System_Download   , only : download
    use :: System_Compilers  , only : compiler         , compilerOptions      , languageC
    implicit none
    type   (varying_string), intent(  out)           :: classPath, classVersion
    logical                , intent(in   ), optional :: static
    integer                                          :: status
    type   (varying_string)                          :: command
    type   (lockDescriptor)                          :: fileLock
    !![
    <optionalArgument name="static" defaultsTo=".false." />
    !!]

    ! Set path and version
    classVersion=dependencyVersion("class")
    classPath   =inputPath(pathTypeDataDynamic)//"class_public-"//classVersion//"/"
    ! Build the CLASS code.
    if (.not.File_Exists(classPath//"class")) then
       call Directory_Make(     classPath                                        )
       call File_Lock     (char(classPath//"class"),fileLock,lockIsShared=.false.)
       ! Unpack the code.
       if (.not.File_Exists(classPath//"Makefile")) then
          ! Download CLASS if necessary.
          if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"class_public-"//char(classVersion)//".tar.gz")) then
             call displayMessage("downloading CLASS code....",verbosityLevelWorking)
             call download("https://github.com/lesgourg/class_public/archive/refs/tags/v"//char(classVersion)//".tar.gz",char(inputPath(pathTypeDataDynamic))//"class_public-"//char(classVersion)//".tar.gz",status=status,retries=5,retryWait=60)
             if (status /= 0 .or. .not.File_Exists(inputPath(pathTypeDataDynamic)//"class_public-"//char(classVersion)//".tar.gz")) call Error_Report("unable to download CLASS"//{introspection:location})
          end if
          call displayMessage("unpacking CLASS code....",verbosityLevelWorking)
          call System_Command_Do("tar -x -v -z -C "//inputPath(pathTypeDataDynamic)//" -f "//inputPath(pathTypeDataDynamic)//"class_public-"//char(classVersion)//".tar.gz",status)
          if (status /= 0 .or. .not.File_Exists(classPath)) call Error_Report('failed to unpack CLASS code'//{introspection:location})        
       end if
       call displayMessage("compiling CLASS code",verbosityLevelWorking)
       command='cd '//classPath//'; cp Makefile Makefile.tmp; '
       ! Include Galacticus compilation flags here.
       command=command//'sed -E -i~ s/"^CC[[:space:]]+=[[:space:]]+gcc"/"CC='//char(compiler(languageC))//'"/ Makefile.tmp; sed -E -i~ s/"^(CC|LD)FLAG = "/"\1FLAG = '//char(stringSubstitute(compilerOptions(languageC),"/","\/"))
       if (static_) command=command//' -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive'
       command=command//' "/ Makefile.tmp; make -f Makefile.tmp -j1 class'
       call System_Command_Do(char(command),status);
       if (status /= 0 .or. .not.File_Exists(classPath//"class")) call Error_Report("failed to build CLASS code"//{introspection:location})
       call File_Unlock(fileLock)
    end if
    return

  contains

    function flagsRetrieve(flagsLength)
      !!{
      Retrieve the compiler flags.
      !!}
      implicit none
      type     (varying_string )                :: flagsRetrieve
      integer                   , intent(in   ) :: flagsLength
      character(len=flagsLength)                :: flags

      call Get_Environment_Variable('GALACTICUS_FCFLAGS',value=flags)
      flagsRetrieve=replace(flags,"/","\/",every=.true.)
      return
    end function flagsRetrieve

  end subroutine Interface_CLASS_Initialize

  subroutine Interface_CLASS_Transfer_Function(cosmologyParameters_,redshifts,wavenumberRequired,wavenumberMaximum,countPerDecade,fileName,wavenumberMaximumReached,transferFunctionDarkMatter,transferFunctionBaryons)
    !!{
    Run CLASS as necessary to compute transfer functions.
    !!}
    use               :: Cosmology_Parameters            , only : cosmologyParametersClass    , hubbleUnitsLittleH
    use               :: File_Utilities                  , only : Count_Lines_In_File         , Directory_Make     , File_Exists   , File_Lock     , &
         &                                                        File_Path                   , File_Remove        , File_Unlock   , lockDescriptor
    use               :: Error                           , only : Error_Report
    use               :: Input_Paths                     , only : inputPath                   , pathTypeDataDynamic
    use               :: HDF5                            , only : hsize_t
    use               :: Hashes_Cryptographic            , only : Hash_MD5
    use               :: HDF5_Access                     , only : hdf5Access
    use               :: IO_HDF5                         , only : hdf5Object
    use   , intrinsic :: ISO_C_Binding                   , only : c_size_t
    use               :: ISO_Varying_String              , only : assignment(=)               , char               , extract       , len           , &
          &                                                       operator(//)                , operator(==)       , varying_string
    use               :: Input_Parameters                , only : inputParameters
#ifdef USEMPI
    use               :: MPI_Utilities                   , only : mpiSelf
#endif
    use               :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial
    use               :: Numerical_Interpolation         , only : GSL_Interp_cSpline
    !$ use            :: OMP_Lib                         , only : OMP_Get_Thread_Num
    use               :: Sorting                         , only : sortIndex
    use               :: String_Handling                 , only : String_C_To_Fortran         , operator(//)
    use               :: System_Command                  , only : System_Command_Do
    use               :: Table_Labels                    , only : extrapolationTypeExtrapolate
    use               :: Tables                          , only : table                       , table1DGeneric
    implicit none
    class           (cosmologyParametersClass), intent(inout)                   :: cosmologyParameters_
    double precision                          , intent(in   ), dimension(:    ) :: redshifts
    double precision                          , intent(in   )                   :: wavenumberRequired                      , wavenumberMaximum
    integer                                   , intent(in   ), optional         :: countPerDecade
    type            (varying_string          ), intent(  out), optional         :: fileName
    type            (table1DGeneric          ), intent(  out), optional         :: transferFunctionDarkMatter              , transferFunctionBaryons
    logical                                   , intent(inout), optional         :: wavenumberMaximumReached
    double precision                          , allocatable  , dimension(:    ) :: wavenumbers                             , wavenumbersLogarithmic  , &
         &                                                                         transferFunctionLogarithmic             , redshiftsCombined
    double precision                          , allocatable  , dimension(:,:,:) :: transferFunctions
    character       (len= 9                  ), allocatable  , dimension(:    ) :: redshiftLabels                          , redshiftLabelsCombined
    integer         (c_size_t                ), allocatable  , dimension(:    ) :: redshiftRanks                           , redshiftRanksCombined
    type            (varying_string          ), allocatable  , dimension(:    ) :: datasetNames
    integer         (hsize_t                 ), parameter                       :: chunkSize                   =100_hsize_t
    type            (lockDescriptor          )                                  :: fileLock
    integer                                                                     :: i                                       , j                       , &
         &                                                                         countRedshiftsUnique
    type            (hdf5Object              )                                  :: classOutput                             , parametersGroup         , &
         &                                                                         extrapolationWavenumberGroup            , extrapolationGroup      , &
         &                                                                         speciesGroup
    character       (len=32                  )                                  :: parameterLabel                          , datasetName             , &
         &                                                                         redshiftLabel
    type            (varying_string          )                                  :: uniqueLabel                             , fileName_
    type            (inputParameters         )                                  :: descriptor
    logical                                                                     :: allEpochsFound
    !![
    <optionalArgument name="countPerDecade" defaultsTo="0"/>
    !!]

    ! Build a sorted array of all redshift labels.
    allocate(redshiftRanks (size(redshifts)))
    allocate(redshiftLabels(size(redshifts)))
    redshiftRanks=sortIndex(redshifts)
    do i=1,size(redshifts)
       write (redshiftLabels(i),'(f9.4)') redshifts(redshiftRanks(i))
    end do
    ! Get a constructor descriptor for this object.
    descriptor=inputParameters()
    call cosmologyParameters_%descriptor(descriptor)
    ! Add primordial helium abundance to the descriptor.
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    call descriptor%addParameter("Y_He"          ,parameterLabel)
    ! Add wavenumber resolution to descriptor.
    write (parameterLabel,'(i4)'  ) countPerDecade_
    call descriptor%addParameter("countPerDecade",parameterLabel)
    ! Add the unique label string to the descriptor.
    uniqueLabel=descriptor%serializeToString()       // &
         &      "_sourceDigest:"                     // &
         &      String_C_To_Fortran(classSourceDigest)
    call descriptor%destroy()
    ! Build the file name.
    fileName_=char(inputPath(pathTypeDataDynamic))                        // &
         &                  'largeScaleStructure/transfer_function_CLASS_'// &
         &                  Hash_MD5(uniqueLabel)                         // &
         &                  '.hdf5'
    if (present(fileName)) fileName=fileName_
    ! Create the directory.
    call Directory_Make(File_Path(fileName_))
    ! If the file exists but has not yet been read, read it now.
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName_),fileLock)
    allEpochsFound=.false.
    if (File_Exists(fileName_)) then
       allEpochsFound=.true.
       !$ call hdf5Access%set()
       call    classOutput%openFile(char(fileName_))
       call    classOutput%readDataset           ('wavenumber',wavenumbers                                    )
       allocate(transferFunctions(size(wavenumbers),3,size(redshifts)))
       speciesGroup=classOutput%openGroup('darkMatter')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(redshiftRanks(i))))
          if (speciesGroup%hasDataset(datasetName)) then
             call speciesGroup%readDatasetStatic(datasetName,transferFunctions(:,classSpeciesDarkMatter%ID,i))
          else
             allEpochsFound=.false.
          end if
       end do
       call speciesGroup%close()
       speciesGroup=classOutput%openGroup('baryons')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(redshiftRanks(i))))
          if (speciesGroup%hasDataset(datasetName)) then
             call speciesGroup%readDatasetStatic(datasetName,transferFunctions(:,classSpeciesBaryons   %ID,i))
          else
             allEpochsFound=.false.
          end if
       end do
       call   speciesGroup%close()
       call   classOutput %close()
       !$ call hdf5Access %unset()
    end if
    if (.not.allocated(wavenumbers) .or. wavenumberRequired > wavenumbers(size(wavenumbers)) .or. .not.allEpochsFound) then
       ! If the wavenumber is out of range, or if not all requested epochs exist within the file, recompute the CLASS transfer function.
       ! Find all existing epochs in the file, create a union of these and the requested epochs.
       if (File_Exists(fileName_)) then
          !$ call hdf5Access%set       (               )
          call    classOutput%openFile  (char(fileName_))
          speciesGroup=classOutput%openGroup('darkMatter')
          call    speciesGroup%datasets(datasetNames   )
          call    speciesGroup%close   (               )
          call    classOutput  %close   (               )
          !$ call hdf5Access  %unset   (               )
       else
          allocate(datasetNames(0))
       end if
       allocate(redshiftsCombined(size(redshifts)+size(datasetNames)))
       redshiftsCombined(1:size(redshifts))=redshifts
       do i=1,size(datasetNames)
          if (extract(datasetNames(i),1,17) == 'transferFunctionZ') then
             redshiftLabel=extract(datasetNames(i),18,len(datasetNames(i)))
             read (redshiftLabel,*) redshiftsCombined(size(redshifts)+i)
          else
             call Error_Report('unknown dataset'//{introspection:location})
          end if
       end do
       allocate(redshiftRanksCombined (size(redshiftsCombined)))
       allocate(redshiftLabelsCombined(size(redshiftsCombined)))
       redshiftRanksCombined=sortIndex(redshiftsCombined)
       do i=1,size(redshiftsCombined)
          write (redshiftLabelsCombined(i),'(f9.4)') redshiftsCombined(redshiftRanksCombined(i))
       end do
       ! Remove duplicated redshifts.
       countRedshiftsUnique=size(redshiftLabelsCombined)
       i=2
       do while (i <= countRedshiftsUnique)
          if (redshiftLabelsCombined(i) == redshiftLabelsCombined(i-1)) then
             do j=i,countRedshiftsUnique-1
                redshiftLabelsCombined(j)=redshiftLabelsCombined(j+1)
             end do
             countRedshiftsUnique=countRedshiftsUnique-1
          else
             i=i+1
          end if
       end do
       ! Run CLASS to get transfer functions.
       call Interface_CLASS_Run(cosmologyParameters_,redshiftsCombined,wavenumberRequired,wavenumberMaximum,countPerDecade,wavenumberMaximumReached,wavenumbers,transferFunctions)
       ! Construct the output HDF5 file.
       !$ call hdf5Access  %set     (                                          )
       call    classOutput %openFile(char(fileName_),objectsOverwritable=.true.)
       call    classOutput %writeAttribute('Transfer functions created by CLASS.','description')
       call    classOutput %writeAttribute(classFormatVersionCurrent,'fileFormat')
       call    classOutput %writeDataset(wavenumbers ,'wavenumber'                               ,chunkSize=chunkSize,appendTo=.not. classOutput%hasDataset('wavenumber'))
       speciesGroup=classOutput%openGroup('darkMatter','Group containing transfer functions for dark matter.')
       do i=1,countRedshiftsUnique
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabelsCombined(i)))
          call speciesGroup%writeDataset(transferFunctions(:,classSpeciesDarkMatter%ID,i),datasetName,chunkSize=chunkSize,appendTo=.not.speciesGroup%hasDataset(datasetName ))
       end do
       call speciesGroup%close()
       speciesGroup=classOutput%openGroup('baryons'   ,'Group containing transfer functions for baryons.'    )
       do i=1,countRedshiftsUnique
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabelsCombined(i)))
          call speciesGroup%writeDataset(transferFunctions(:,classSpeciesBaryons   %ID,i),datasetName,chunkSize=chunkSize,appendTo=.not.speciesGroup%hasDataset(datasetName ))
       end do
       call speciesGroup%close()
       parametersGroup=classOutput%openGroup('parameters')
       call parametersGroup%writeAttribute(cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
       call parametersGroup%writeAttribute(cosmologyParameters_%OmegaBaryon    (),'OmegaBaryon'    )
       call parametersGroup%writeAttribute(cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
       call parametersGroup%writeAttribute(cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
       call parametersGroup%writeAttribute(cosmologyParameters_%temperatureCMB (),'temperatureCMB' )
       call parametersGroup%close()
       extrapolationGroup          =classOutput       %openGroup('extrapolation')
       extrapolationWavenumberGroup=extrapolationGroup%openGroup('wavenumber'   )
       call    extrapolationWavenumberGroup%writeAttribute('extrapolate','low' )
       call    extrapolationWavenumberGroup%writeAttribute('extrapolate','high')
       call    extrapolationWavenumberGroup%close()
       call    extrapolationGroup          %close()
       call    classOutput                 %close()
       !$ call hdf5Access                  %unset()
    end if
    ! If necessary, construct tables of transfer functions.
    if (present(transferFunctionDarkMatter)) then
       !$ call hdf5Access%set()
       call classOutput%openFile(char(fileName_))
       call classOutput%readDataset('wavenumber',wavenumbersLogarithmic)
       wavenumbersLogarithmic=log(wavenumbersLogarithmic)
       call transferFunctionDarkMatter%create(                                                 &
            &                                                   wavenumbersLogarithmic       , &
            &                                 tableCount       =size(redshifts)              , &
            &                                 extrapolationType=[                              &
            &                                                    extrapolationTypeExtrapolate, &
            &                                                    extrapolationTypeExtrapolate  &
            &                                                   ]                            , &
            &                                 interpolationType=GSL_Interp_cSpline             &
            &                                )
       deallocate(wavenumbersLogarithmic)
       speciesGroup=classOutput%openGroup('darkMatter')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(i)))
          call speciesGroup%readDataset(datasetName,transferFunctionLogarithmic)
          transferFunctionLogarithmic=log(transferFunctionLogarithmic)
          call transferFunctionDarkMatter%populate(transferFunctionLogarithmic,table=int(redshiftRanks(i)))
          deallocate(transferFunctionLogarithmic)
       end do
       call    speciesGroup%close()
       call    classOutput  %close()
       !$ call hdf5Access  %unset()
    end if
    if (present(transferFunctionBaryons)) then
       !$ call hdf5Access%set()
       call    classOutput%openFile(char(fileName_))
       call    classOutput%readDataset('wavenumber',wavenumbersLogarithmic)
       wavenumbersLogarithmic=log(wavenumbersLogarithmic)
       call transferFunctionBaryons   %create(                                                 &
            &                                                   wavenumbersLogarithmic       , &
            &                                 tableCount       =size(redshifts)              , &
            &                                 extrapolationType=[                              &
            &                                                    extrapolationTypeExtrapolate, &
            &                                                    extrapolationTypeExtrapolate  &
            &                                                   ]                            , &
            &                                 interpolationType=GSL_Interp_cSpline             &
            &                                )
       deallocate(wavenumbersLogarithmic)
       speciesGroup=classOutput%openGroup('baryons')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(i)))
          call speciesGroup%readDataset(datasetName,transferFunctionLogarithmic)
          transferFunctionLogarithmic=log(transferFunctionLogarithmic)
          call transferFunctionBaryons   %populate(transferFunctionLogarithmic,table=int(redshiftRanks(i)))
          deallocate(transferFunctionLogarithmic)
       end do
       call    speciesGroup%close()
       call    classOutput  %close()
       !$ call hdf5Access  %unset()
    end if
    ! Unlock the file.
    call File_Unlock(fileLock)
    return
  end subroutine Interface_CLASS_Transfer_Function

  double precision function Interface_CLASS_Normalization(cosmologyParameters_) result(normalization)
    !!{
    Run CLASS to compute the ratio $\sigma_8^2/A_s$.
    !!}
    use               :: Cosmology_Parameters            , only : cosmologyParametersClass
    use               :: File_Utilities                  , only : Directory_Make          , File_Exists       , File_Lock     , File_Unlock, &
         &                                                        File_Path               , lockDescriptor
    use               :: Input_Paths                     , only : inputPath               , pathTypeDataDynamic
    use               :: Hashes_Cryptographic            , only : Hash_MD5
    use               :: HDF5_Access                     , only : hdf5Access
    use               :: IO_HDF5                         , only : hdf5Object
    use               :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial
    use   , intrinsic :: ISO_C_Binding                   , only : c_size_t
    use               :: ISO_Varying_String              , only : varying_string          , char               , operator(//)
    use               :: Input_Parameters                , only : inputParameters
    use               :: String_Handling                 , only : String_C_To_Fortran     , operator(//)
    use               :: System_Command                  , only : System_Command_Do
    implicit none
    class    (cosmologyParametersClass), intent(inout)                   :: cosmologyParameters_
    type     (lockDescriptor          )                                  :: fileLock
    type     (hdf5Object              )                                  :: classOutput
    character(len=32                  )                                  :: parameterLabel
    type     (varying_string          )                                  :: uniqueLabel                , fileName
    type     (inputParameters         )                                  :: descriptor
  
    ! Get a constructor descriptor for this object.
    descriptor=inputParameters()
    call cosmologyParameters_%descriptor(descriptor)
    ! Add primordial helium abundance to the descriptor.
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    call descriptor%addParameter("Y_He",parameterLabel)
    ! Add the unique label string to the descriptor.
    uniqueLabel=descriptor%serializeToString()        // &
         &      "_sourceDigest:"                      // &
         &      String_C_To_Fortran(classSourceDigest)
    call descriptor%destroy()
    ! Build the file name.
    fileName=char(inputPath(pathTypeDataDynamic))                    // &
         &                 'largeScaleStructure/normalization_CLASS_'// &
         &                 Hash_MD5(uniqueLabel)                     // &
         &                 '.hdf5'
    ! Create the directory.
    call Directory_Make(File_Path(fileName))
    ! If the file exists but has not yet been read, read it now.
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName),fileLock)
    if (File_Exists(fileName)) then
       !$ call hdf5Access %set          (                                 )
       call    classOutput%openFile     (char(fileName     )              )
       call    classOutput%readAttribute(    'normalization',normalization)      
       call    classOutput%close        (                                 )
       !$ call hdf5Access %unset        (                                 )
    else
       ! Run CLASS to compute the normalization.
       call Interface_CLASS_Run(cosmologyParameters_,ratioSigma8SquaredAs=normalization)
       ! Construct the output HDF5 file.
       !$ call hdf5Access %set           (                                              )
       call    classOutput%openFile      (char(fileName     ),objectsOverwritable=.true.)
       call    classOutput%writeAttribute(     normalization ,'normalization'           )
       call    classOutput%close         (                                              )
       !$ call hdf5Access %unset         (                                              )
    end if
    ! Unlock the file.
    call File_Unlock(fileLock)
    return
  end function Interface_CLASS_Normalization

  subroutine Interface_CLASS_Run(cosmologyParameters_,redshifts,wavenumberRequired,wavenumberMaximum,countPerDecade,wavenumberMaximumReached,wavenumbers,transferFunctions,ratioSigma8SquaredAs)
    !!{
    Run CLASS and extract required information.
    !!}
#ifdef USEMPI
    use               :: MPI_Utilities                   , only : mpiSelf
#endif
    !$ use            :: OMP_Lib                         , only : OMP_Get_Thread_Num
    use               :: Cosmology_Parameters            , only : cosmologyParametersClass, hubbleUnitsLittleH
    use               :: Error                           , only : Error_Report
    use               :: File_Utilities                  , only : Directory_Make          , File_Remove        , Directory_Remove, Count_Lines_In_File
    use               :: Input_Paths                     , only : inputPath               , pathTypeDataDynamic
    use   , intrinsic :: ISO_C_Binding                   , only : c_size_t
    use               :: ISO_Varying_String              , only : varying_string          , char               , operator(//)
    use               :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial
    use               :: Sorting                         , only : sortIndex
    use               :: String_Handling                 , only : operator(//)
    use               :: System_Command                  , only : System_Command_Do
    implicit none
    class           (cosmologyParametersClass), intent(inout)                                          :: cosmologyParameters_
    double precision                          , intent(in   ), optional, dimension(:    )              :: redshifts
    double precision                          , intent(in   ), optional                                :: wavenumberRequired               , wavenumberMaximum
    integer                                   , intent(in   ), optional                                :: countPerDecade
    logical                                   , intent(inout), optional                                :: wavenumberMaximumReached
    double precision                          , intent(inout), optional, dimension(:    ), allocatable :: wavenumbers
    double precision                          , intent(inout), optional, dimension(:,:,:), allocatable :: transferFunctions
    double precision                          , intent(  out), optional                                :: ratioSigma8SquaredAs
    double precision                          , parameter                                              :: normalizationAs           =1.0d-9
    integer         (c_size_t                ), allocatable            , dimension(:    )              :: redshiftRanks
    character       (len= 9                  ), allocatable            , dimension(:    )              :: redshiftLabels
    character       (len=255                 )                                                         :: hostName                         , classTransferLine
    character       (len= 32                 )                                                         :: line                             , outputs
    integer         (c_size_t                )                                                         :: countWavenumber
    integer                                                                                            :: classParameterFile               , classLog               , &
         &                                                                                                offset                           , status                 , &
         &                                                                                                classTransferFile                , i                      , &
         &                                                                                                j
    type            (varying_string          )                                                         :: classVersion                     , parameterFile          , &
         &                                                                                                classPath                        , workPath               , &
         &                                                                                                transferFileName
    logical                                                                                            :: found                            , haveTransferFunctions  , &
         &                                                                                                haveNormalization
    double precision                                                                                   :: sigma8                           , wavenumberCLASS
    !![
    <optionalArgument name="countPerDecade" defaultsTo="0"/>
    !!]

    ! Ensure CLASS is initialized.
    call Interface_CLASS_Initialize(classPath,classVersion)
    ! Determine the types of output needed.
    outputs=''
    haveNormalization=present(ratioSigma8SquaredAs)
    if (haveNormalization) then
       outputs=trim(outputs)//' mPk'
    end if
    haveTransferFunctions=present(transferFunctions)
    if (haveTransferFunctions) then
       outputs=trim(outputs)//' mTk'
       if (.not.present(redshifts         )) call Error_Report('`redshifts` is required for transfer function calculation'         //{introspection:location})
       if (.not.present(wavenumberRequired)) call Error_Report('`wavenumberRequired` is required for transfer function calculation'//{introspection:location})
       if (.not.present(wavenumberMaximum )) call Error_Report('`wavenumberMaximum` is required for transfer function calculation' //{introspection:location})
       if (.not.present(wavenumbers       )) call Error_Report('`wavenumbers` is required for transfer function calculation'       //{introspection:location})
    end if
    if (trim(outputs) == 'no outputs added') call Error_Report(''//{introspection:location})
    ! Determine maximum wavenumber.
    if (haveTransferFunctions) then
       wavenumberCLASS=exp(max(log(wavenumberRequired)+1.0d0,classLogWavenumberMaximumDefault))
       if (wavenumberCLASS > wavenumberMaximum) then
          wavenumberCLASS=wavenumberMaximum
          if (present(wavenumberMaximumReached)) wavenumberMaximumReached=.true.
       end if
       if (allocated(wavenumbers)) wavenumberCLASS=max(wavenumberCLASS,wavenumbers(size(wavenumbers)))
    end if
    ! Generate redshift labels.
    if (haveTransferFunctions) then
       allocate(redshiftRanks (size(redshifts)))
       allocate(redshiftLabels(size(redshifts)))
       redshiftRanks=sortIndex(redshifts)
       do i=1,size(redshifts)
          write (redshiftLabels(i),'(f9.4)') redshifts(redshiftRanks(i))
       end do
    end if
    ! Construct input file for CLASS.
    call Get_Environment_Variable('HOSTNAME',hostName)
    workPath        =inputPath(pathTypeDataDynamic)//'largeScaleStructure/class_work_'//trim(hostName)//'_'//GetPID()
    !$ workPath     =workPath     //'_'//OMP_Get_Thread_Num()
    !$ parameterFile=parameterFile//'_'//OMP_Get_Thread_Num()
#ifdef USEMPI
    workPath        =workPath     //'_'//mpiSelf%rankLabel()
#endif
    parameterFile=workPath//'/parameters.ini'
    call Directory_Make(workPath)
    open(newunit=classParameterFile,file=char(parameterFile),status='unknown',form='formatted')
    write (classParameterFile,'(a,1x,"=",1x,a,a  )') 'root            ',char(workPath),'/output'
    write (classParameterFile,'(a,1x,"=",1x,a    )') 'overwrite_root  ','yes'
    write (classParameterFile,'(a,1x,"=",1x,a    )') 'headers         ','yes'
    write (classParameterFile,'(a,1x,"=",1x,a    )') 'format          ','class'
    write (classParameterFile,'(a,1x,"=",1x,a    )') 'output          ',trim(adjustl(outputs))
    write (classParameterFile,'(a,1x,"=",1x,a    )') 'modes           ','s'
    write (classParameterFile,'(a,1x,"=",1x,a    )') 'ic              ','ad'
    write (classParameterFile,'(a,1x,"=",1x,e12.6)') 'h               ',                                                                                     cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    write (classParameterFile,'(a,1x,"=",1x,e12.6)') 'T_cmb           ',       cosmologyParameters_%temperatureCMB()
    write (classParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_b         ',(      cosmologyParameters_%OmegaBaryon   ()                                       )*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
    write (classParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_cdm       ',(      cosmologyParameters_%OmegaMatter   ()-cosmologyParameters_%OmegaBaryon    ())*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
    write (classParameterFile,'(a,1x,"=",1x,e12.6)') 'Omega_k         ',(1.0d0-cosmologyParameters_%OmegaMatter   ()-cosmologyParameters_%OmegaDarkEnergy())*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
    write (classParameterFile,'(a,1x,"=",1x,e12.6)') 'N_ur            ',3.044d+0
    write (classParameterFile,'(a,1x,"=",1x,e12.6)') 'YHe             ',heliumByMassPrimordial
    write (classParameterFile,'(a,1x,"=",1x,e12.6)') 'A_s             ',normalizationAs
    write (classParameterFile,'(a,1x,"=",1x,i1   )') 'fourier_verbose ',1
    if (haveTransferFunctions) then
       write        (classParameterFile,'(a,1x,"=",1x,e12.6)') 'P_k_max_1/Mpc      ',wavenumberCLASS
       if (countPerDecade_ > 0) &
            & write (classParameterFile,'(a,1x,"=",1x,i3   )') 'k_per_decade_for_pk',countPerDecade_
       write        (classParameterFile,'(a,1x,"=",$       )') 'z_pk               '
       do i=size(redshifts),1,-1
          write     (classParameterFile,'(1x,a,$           )') trim(adjustl(redshiftLabels(i)))
          if (i > 1 ) write (classParameterFile,'(",",$)')
       end do
    end if
    write (classParameterFile,'(a)') ''
    close(classParameterFile)
    ! Run CLASS.
    call System_Command_Do(classPath//"class "//parameterFile//" > "//workPath//"/class.log")
    ! Extract the ratio σ₈²/Aₛ.
    if (haveNormalization) then
       found=.false.
       open(newUnit=classLog,file=char(workPath)//"/class.log",form="formatted",status="old",iostat=status)
       do while (status == 0)
          read (classLog,'(a)',iostat=status) line
          offset=index(line,"sigma8=")
          if (offset > 0) then
             found=.true.
             read (line(offset+7:),*) sigma8
          end if
       end do
       close(classLog)
       if (.not.found) call Error_Report("CLASS σ₈ value not found"//{introspection:location})
       ratioSigma8SquaredAs=+sigma8         **2 &
            &               /normalizationAs
    end if
    ! Read the CLASS transfer function file.
    if (haveTransferFunctions) then
       if (allocated(wavenumbers      )) deallocate(wavenumbers      )
       if (allocated(transferFunctions)) deallocate(transferFunctions)
       allocate(wavenumbers      (0    ))
       allocate(transferFunctions(0,0,0))
       countWavenumber=0
       do j=1,size(redshifts)
          if (size(redshifts) > 1) then
             transferFileName=workPath//'/output_z'//j//'_tk.dat'
          else
             transferFileName=workPath//'/output_tk.dat'
          end if
          if (j == 1) then
             countWavenumber=Count_Lines_In_File(transferFileName,"#")
             if (allocated(wavenumbers      )) deallocate(wavenumbers      )
             if (allocated(transferFunctions)) deallocate(transferFunctions)
             allocate(wavenumbers      (countWavenumber                  ))
             allocate(transferFunctions(countWavenumber,3,size(redshifts)))
          end if
          open(newunit=classTransferFile,file=char(transferFileName),status='old',form='formatted')
          i=0
          do while (i < countWavenumber)
             read (classTransferFile,'(a)',iostat=status) classTransferLine
             if (status == 0) then
                if (classTransferLine(1:1) /= "#") then
                   i=i+1
                   read (classTransferLine,*) wavenumbers(i),transferFunctions(i,classSpeciesPhotons%ID,j),transferFunctions(i,classSpeciesBaryons%ID,j),transferFunctions(i,classSpeciesDarkMatter%ID,j)
                end if
             else
                call Error_Report('unable to read CLASS transfer function file'//{introspection:location})
             end if
          end do
          close(classTransferFile)
       end do
       ! Convert from CLASS units to Galacticus units.
       wavenumbers=+wavenumbers                                                   &
            &      *cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)
       ! Convert transfer functions to standard form.
       do i=1,3
          do j=1,size(redshifts)
             transferFunctions(:,i,j)=-transferFunctions(:,i,j)/wavenumbers**2
          end do
       end do
    end if
    ! Remove temporary files.
    call File_Remove(parameterFile             )
    call File_Remove(workPath//"/class.log"    )
    if (haveNormalization) then
       call File_Remove(workPath//"/output_pk.dat")
    end if
    if (haveTransferFunctions) then
       do i=1,size(redshifts)
          if (size(redshifts) > 1) then
             transferFileName=workPath//'/output_z'//i//'_tk.dat'
          else
             transferFileName=workPath//'/output_tk.dat'
          end if
          call File_Remove(transferFileName)
       end do
    end if
    call Directory_Remove(workPath)
    return
  end subroutine Interface_CLASS_Run

end module Interfaces_CLASS
