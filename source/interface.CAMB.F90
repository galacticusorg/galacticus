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
Contains a module which provides various interfaces to the \gls{camb} code.
!!}

module Interfaces_CAMB
  !!{
  Provides various interfaces to the \gls{camb} code.
  !!}
  use :: File_Utilities, only : lockDescriptor
  private
  public :: Interface_CAMB_Initialize, Interface_CAMB_Transfer_Function

  ! Current file format version for transfer function files. Note that this file format matches that used by the "file" transfer
  ! function class.
  integer                         , parameter :: cambFormatVersionCurrent       =     2

  ! Default maximum wavenumber to tabulate.
  double precision                , parameter :: cambLogWavenumberMaximumDefault=log(10.0d0)

  !![
  <enumeration>
   <name>cambSpecies</name>
   <description>Particle species in CAMB.</description>
   <visibility>public</visibility>
   <indexing>1</indexing>
   <entry label="darkMatter"/>
   <entry label="baryons"   />
  </enumeration>
  !!]

  ! Generate a source digest.
  !![
  <sourceDigest name="cambSourceDigest"/>
  !!]

contains

  subroutine Interface_CAMB_Initialize(cambPath,cambVersion,static)
    !!{
    Initialize the interface with CAMB, including downloading and compiling CAMB if necessary.
    !!}
    use :: Dependencies      , only : dependencyVersion
    use :: Display           , only : displayMessage   , verbosityLevelWorking
    use :: File_Utilities    , only : Directory_Make   , File_Exists          , File_Lock      , File_Unlock, &
          &                           lockDescriptor
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath        , pathTypeDataDynamic
    use :: ISO_Varying_String, only : assignment(=)    , char                 , operator(//)   , replace    , &
          &                           varying_string
    use :: String_Handling   , only : stringSubstitute
    use :: System_Command    , only : System_Command_Do
    use :: System_Download   , only : download
    use :: System_Compilers  , only : compiler         , compilerOptions      , languageFortran
    implicit none
    type   (varying_string), intent(  out)           :: cambPath, cambVersion
    logical                , intent(in   ), optional :: static
    integer                                          :: status
    type   (varying_string)                          :: command , forutilsVersion
    type   (lockDescriptor)                          :: fileLock
    !![
    <optionalArgument name="static" defaultsTo=".false." />
    !!]

    ! Set path and version
    cambVersion    =dependencyVersion("camb"    )
    forutilsVersion=dependencyVersion("forutils")
    cambPath       =inputPath(pathTypeDataDynamic)//"CAMB-"//cambVersion//"/fortran/"
    ! Build the CAMB code.
    if (.not.File_Exists(cambPath//"camb")) then
       call Directory_Make(     cambPath                                       )
       call File_Lock     (char(cambPath//"camb"),fileLock,lockIsShared=.false.)
       ! Unpack the code.
       if (.not.File_Exists(cambPath//"Makefile")) then
          ! Download CAMB if necessary.
          if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"CAMB_"//char(cambVersion)//".tar.gz")) then
             call displayMessage("downloading CAMB code....",verbosityLevelWorking)
             call download("https://github.com/cmbant/CAMB/archive/refs/tags/"//char(cambVersion)//".tar.gz",char(inputPath(pathTypeDataDynamic))//"CAMB_"//char(cambVersion)//".tar.gz",status=status,retries=5,retryWait=60)
             if (status /= 0 .or. .not.File_Exists(inputPath(pathTypeDataDynamic)//"CAMB_"//char(cambVersion)//".tar.gz")) call Error_Report("unable to download CAMB"//{introspection:location})
          end if
          call displayMessage("unpacking CAMB code....",verbosityLevelWorking)
          call System_Command_Do("tar -x -v -z -C "//inputPath(pathTypeDataDynamic)//" -f "//inputPath(pathTypeDataDynamic)//"CAMB_"//char(cambVersion)//".tar.gz",status);
          if (status /= 0 .or. .not.File_Exists(cambPath)) call Error_Report('failed to unpack CAMB code'//{introspection:location})
          ! Download the "forutils" package if necessary.
          if (.not.File_Exists(cambPath//"../forutils/Makefile")) then
             if (.not.File_Exists(cambPath//"../forutils_"//char(forutilsVersion)//".tar.gz")) then
                call displayMessage("downloading forutils code....",verbosityLevelWorking)
                call download("https://github.com/cmbant/forutils/archive/refs/tags/"//char(forutilsVersion)//".tar.gz",char(cambPath)//"../forutils_"//char(forutilsVersion)//".tar.gz",status=status,retries=5,retryWait=60)
                if (status /= 0 .or. .not.File_Exists(cambPath//"../forutils_"//char(forutilsVersion)//".tar.gz")) call Error_Report("unable to download forutils"//{introspection:location})
             end if
             call displayMessage("unpacking forutils code....",verbosityLevelWorking)
             call System_Command_Do("tar -x -v -z -C "//cambPath//"../forutils -f "//cambPath//"../forutils_"//char(forutilsVersion)//".tar.gz --strip-components 1");          
             if (status /= 0 .or. .not.File_Exists(cambPath//"../forutils/Makefile")) call Error_Report('failed to unpack forutils code'//{introspection:location})
          end if
       end if       
       call displayMessage("compiling CAMB code",verbosityLevelWorking)
       command='cd '//cambPath//'; sed -E -i~ s/"ifortErr[[:space:]]*=.*"/"ifortErr = 1"/ Makefile; sed -E -i~ s/"gfortErr[[:space:]]*=.*"/"gfortErr = 0"/ Makefile; sed -E -i~ s/"gfortran"/"'//compiler(languageFortran)//'"/ Makefile; sed -E -i~ s/"gfortran"/"'//compiler(languageFortran)//'"/ ../forutils/Makefile_compiler; sed -E -i~ s/"^FFLAGS[[:space:]]*\+=[[:space:]]*\-march=native"/"FFLAGS+="/ Makefile; sed -E -i~ s/"^FFLAGS[[:space:]]*=[[:space:]]*.*"/"FFLAGS = -cpp -Ofast -fopenmp '//stringSubstitute(compilerOptions(languageFortran),"/","\/")
       if (static_) command=command//" -static -Wl,--whole-archive -lpthread -ldl -Wl,--no-whole-archive"
       command=command//'"/ Makefile'
       if (static_) command=command//'; cp $GALACTICUS_EXEC_PATH/source/utility.OpenMP.workaround.c '//cambPath//'; gcc -DSTATIC -c utility.OpenMP.workaround.c -o utility.OpenMP.workaround.o; sed -E -i~ s/"\-o camb$"/"utility\.OpenMP\.workaround\.o \-o camb"/ Makefile_main'
       command=command//'; find . -name "*.f90" | xargs sed -E -i~ s/"error stop"/"error stop "/; make -j1 camb'
       call System_Command_Do(char(command),status);
       if (status /= 0 .or. .not.File_Exists(cambPath//"camb")) call Error_Report("failed to build CAMB code"//{introspection:location})
       call File_Unlock(fileLock)
    end if
    return
  end subroutine Interface_CAMB_Initialize

  subroutine Interface_CAMB_Transfer_Function(cosmologyParameters_,redshifts,wavenumberRequired,wavenumberMaximum,countPerDecade,fileName,wavenumberMaximumReached,transferFunctionDarkMatter,transferFunctionBaryons)
    !!{
    Run CAMB as necessary to compute transfer functions.
    !!}
    use            :: Cosmology_Parameters            , only : cosmologyParametersClass    , hubbleUnitsLittleH
    use            :: File_Utilities                  , only : Count_Lines_In_File         , Directory_Make     , File_Exists   , File_Lock     , &
         &                                                     File_Path                   , File_Remove        , File_Unlock   , lockDescriptor, &
         &                                                     File_Name_Temporary
    use            :: Error                           , only : Error_Report
    use            :: Input_Paths                     , only : inputPath                   , pathTypeDataDynamic
    use            :: HDF5                            , only : hsize_t
    use            :: Hashes_Cryptographic            , only : Hash_MD5
    use            :: HDF5_Access                     , only : hdf5Access
    use            :: IO_HDF5                         , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: ISO_Varying_String              , only : assignment(=)               , char               , extract       , len           , &
          &                                                    operator(//)                , operator(==)       , varying_string
    use            :: Input_Parameters                , only : inputParameters
    use            :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial
    use            :: Numerical_Interpolation         , only : GSL_Interp_cSpline
    use            :: Sorting                         , only : sortIndex
    use            :: String_Handling                 , only : String_C_To_Fortran         , operator(//)
    use            :: System_Command                  , only : System_Command_Do
    use            :: Table_Labels                    , only : extrapolationTypeExtrapolate
    use            :: Tables                          , only : table                       , table1DGeneric
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
    character       (len=255                 )                                  :: cambTransferLine
    type            (varying_string          )                                  :: cambPath                                , cambVersion             , &
         &                                                                         parameterFile                           , outputRoot
    double precision                                                            :: wavenumberCAMB
    integer                                                                     :: status                                  , cambParameterFile       , &
         &                                                                         i                                       , cambTransferFile        , &
         &                                                                         j                                       , countRedshiftsUnique    , &
         &                                                                         iLock
    integer         (c_size_t                )                                  :: countWavenumber
    type            (hdf5Object              )                                  :: cambOutput                              , parametersGroup         , &
         &                                                                         extrapolationWavenumberGroup            , extrapolationGroup      , &
         &                                                                         speciesGroup
    character       (len=32                  )                                  :: parameterLabel                          , datasetName             , &
         &                                                                         redshiftLabel                           , indexLabel
    type            (varying_string          )                                  :: uniqueLabel                             , workPath                , &
         &                                                                         transferFileName                        , fileName_
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
         &      String_C_To_Fortran(cambSourceDigest)
    call descriptor%destroy()
    ! Build the file name.
    fileName_=char(inputPath(pathTypeDataDynamic))                       // &
         &                  'largeScaleStructure/transfer_function_CAMB_'// &
         &                  Hash_MD5(uniqueLabel)                        // &
         &                  '.hdf5'
    if (present(fileName)) fileName=fileName_
    ! Create the directory.
    call Directory_Make(File_Path(fileName_))
    ! If the file exists but has not yet been read, read it now.
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    do iLock=1,2
       ! Get a lock on the file - first a shared lock, but if we need to remake the file, use an exclusive lock on the second pass.
       call File_Lock(char(fileName_),fileLock,lockIsShared=iLock == 1)
       ! Clean up workspace.
       if (allocated(transferFunctions     )) deallocate(transferFunctions     )
       if (allocated(datasetNames          )) deallocate(datasetNames          )
       if (allocated(redshiftsCombined     )) deallocate(redshiftsCombined     )
       if (allocated(redshiftRanksCombined )) deallocate(redshiftRanksCombined )
       if (allocated(redshiftLabelsCombined)) deallocate(redshiftLabelsCombined)
       ! Search for a existing tabulation.
       allEpochsFound=.false.
       if (File_Exists(fileName_)) then
          allEpochsFound=.true.
          !$ call hdf5Access%set()
          call    cambOutput%openFile(char(fileName_),readOnly=.true.)
          call    cambOutput%readDataset           ('wavenumber'  ,wavenumbers                                 )
          allocate(transferFunctions(size(wavenumbers),2,size(redshifts)))
          speciesGroup=cambOutput%openGroup('darkMatter')
          do i=1,size(redshifts)
             datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(redshiftRanks(i))))
             if (speciesGroup%hasDataset(datasetName)) then
                call speciesGroup%readDatasetStatic(datasetName,transferFunctions(:,cambSpeciesDarkMatter%ID,i))
             else
                allEpochsFound=.false.
             end if
          end do
          call speciesGroup%close()
          speciesGroup=cambOutput%openGroup('baryons')
          do i=1,size(redshifts)
             datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(redshiftRanks(i))))
             if (speciesGroup%hasDataset(datasetName)) then
                call speciesGroup%readDatasetStatic(datasetName,transferFunctions(:,cambSpeciesBaryons   %ID,i))
             else
                allEpochsFound=.false.
             end if
          end do
          call   speciesGroup%close()
          call   cambOutput  %close()
          !$ call hdf5Access  %unset()
       end if
       if (.not.allocated(wavenumbers) .or. wavenumberRequired > wavenumbers(size(wavenumbers)) .or. .not.allEpochsFound) then
          ! The table is insufficient. If this is the first pass, cycle so that we get an exclusive lock, and then check again (in
          ! case the file has been made by another thread/process while we were waiting).
          if (iLock == 1) then
             ! Unlock then cycle so that we can request an exclusive lock.
             call File_Unlock(fileLock,sync=.false.)
             cycle
          else
             ! We now have an exclusive lock and the table is insufficient - remake it now.
             ! Find all existing epochs in the file, create a union of these and the requested epochs.
             if (File_Exists(fileName_)) then
                !$ call hdf5Access%set       (               )
                call    cambOutput%openFile  (char(fileName_))
                speciesGroup=cambOutput%openGroup('darkMatter')
                call    speciesGroup%datasets(datasetNames   )
                call    speciesGroup%close   (               )
                call    cambOutput  %close   (               )
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
             ! Ensure CAMB is initialized.
             call Interface_CAMB_Initialize(cambPath,cambVersion)
             ! Determine maximum wavenumber.
             wavenumberCAMB=exp(max(log(wavenumberRequired)+1.0d0,cambLogWavenumberMaximumDefault))
             if (wavenumberCAMB > wavenumberMaximum) then
                wavenumberCAMB=wavenumberMaximum
                if (present(wavenumberMaximumReached)) wavenumberMaximumReached=.true.
             end if
             if (allocated(wavenumbers)) wavenumberCAMB=max(wavenumberCAMB,wavenumbers(size(wavenumbers)))
             ! Construct input file for CAMB.
             workPath     =inputPath(pathTypeDataDynamic)//'largeScaleStructure/'
             parameterFile=File_Name_Temporary('transfer_function_parameters',char(workPath))//'.txt'
             outputRoot   =File_Name_Temporary('camb'                        ,char(workPath))
             call Directory_Make(workPath)
             open(newunit=cambParameterFile,file=char(parameterFile),status='unknown',form='formatted')
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'output_root                  ',char(outputRoot)
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'get_scalar_cls               ','F'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'get_vector_cls               ','F'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'get_tensor_cls               ','F'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'get_transfer                 ','T'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'do_lensing                   ','F'
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'do_nonlinear                 ',0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'l_max_scalar                 ',2200.0d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'l_max_tensor                 ',1500.0d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'k_eta_max_tensor             ',3000.0d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'ombh2                        ',(      cosmologyParameters_%OmegaBaryon   ()                                       )*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'omch2                        ',(      cosmologyParameters_%OmegaMatter   ()-cosmologyParameters_%OmegaBaryon    ())*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'omk                          ',(1.0d0-cosmologyParameters_%OmegaMatter   ()-cosmologyParameters_%OmegaDarkEnergy())*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'omnuh2                       ',(0.0d0                                                                             )*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'hubble                       ',                                                                                     cosmologyParameters_%HubbleConstant(                  )
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'w                            ',-1.0d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'cs2_lam                      ',1.0d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'temp_cmb                     ',      cosmologyParameters_%temperatureCMB()
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'helium_fraction              ',heliumByMassPrimordial
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'massless_neutrinos           ',2.046d0
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'nu_mass_eigenstates          ',1
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'massive_neutrinos            ',1
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'share_delta_neff             ','T'
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'nu_mass_fractions            ',1.0d0
             write (cambParameterFile,'(a,1x,"=",1x      )') 'nu_mass_degeneracies         '
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'initial_power_num            ',1
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'pivot_scalar                 ',0.05d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'pivot_tensor                 ',0.05d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'scalar_amp(1)                ',2.1d-9
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'scalar_spectral_index(1)     ',0.96d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'scalar_nrun(1)               ',0.0d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'tensor_spectral_index(1)     ',0.0d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'initial_ratio(1)             ',1.0d0
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'reionization                 ','T'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 're_use_optical_depth         ','T'
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 're_optical_depth             ',0.09d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 're_redshift                  ',11.0d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 're_delta_redshift            ',1.5d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 're_ionization_frac           ',-1.0d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'RECFAST_fudge                ',1.14d0
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'RECFAST_fudge_He             ',0.86d0
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'RECFAST_Heswitch             ',6
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'RECFAST_Hswitch              ','T'
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'initial_condition            ',1
             write (cambParameterFile,'(a,1x,"=",1x,5(i2))') 'initial_vector               ',-1,0,0,0,0
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'vector_mode                  ',0
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'COBE_normalize               ','F'
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'CMB_outputscale              ',7.42835025d12
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'transfer_high_precision      ','F'
             write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'transfer_kmax                ',wavenumberCAMB/cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)
             write (cambParameterFile,'(a,1x,"=",1x,i3   )') 'transfer_k_per_logint        ',countPerDecade_
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'transfer_num_redshifts       ',countRedshiftsUnique
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'transfer_interp_matterpower  ','T'
             do i=countRedshiftsUnique,1,-1
                write (indexLabel,'(i4)') countRedshiftsUnique+1-i
                write (cambParameterFile,'(a,a,a,1x,"=",1x,a      )') 'transfer_redshift('   ,trim(adjustl(indexLabel)),')'               ,trim(adjustl(redshiftLabelsCombined(i)))
                write (cambParameterFile,'(a,a,a,1x,"=",1x,a,a,a,a)') 'transfer_filename('   ,trim(adjustl(indexLabel)),')','transfer_'   ,trim(adjustl(redshiftLabelsCombined(i))),'.dat'
                write (cambParameterFile,'(a,a,a,1x,"=",1x,a,a,a,a)') 'transfer_matterpower(',trim(adjustl(indexLabel)),')','matterpower_',trim(adjustl(redshiftLabelsCombined(i))),'.dat'
             end do
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'scalar_output_file           ','scalCls.dat'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'vector_output_file           ','vecCls.dat'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'tensor_output_file           ','tensCls.dat'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'total_output_file            ','totCls.dat'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'lensed_output_file           ','lensedCls.dat'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'lensed_total_output_file     ','lensedtotCls.dat'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'lens_potential_output_file   ','lenspotentialCls.dat'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'FITS_filename                ','scalCls.fits'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'do_lensing_bispectrum        ','F'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'do_primordial_bispectrum     ','F'
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_nfields           ',1
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_slice_base_L      ',0
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_ndelta            ',3
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_delta(1)          ',0
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_delta(2)          ',2
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_delta(3)          ',4
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'bispectrum_do_fisher         ','F'
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_fisher_noise      ',0
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_fisher_noise_pol  ',0
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_fisher_fwhm_arcmin',7
             write (cambParameterFile,'(a,1x,"=",1x      )') 'bispectrum_full_output_file  '
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'bispectrum_full_output_sparse','F'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'bispectrum_export_alpha_beta ','F'
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'feedback_level               ',1
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'derived_parameters           ','T'
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'lensing_method               ',1
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'accurate_BB                  ','F'
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'massive_nu_approx            ',1
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'accurate_polarization        ','T'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'accurate_reionization        ','T'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'do_tensor_neutrinos          ','T'
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'do_late_rad_truncation       ','T'
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'number_of_threads            ',0
             write (cambParameterFile,'(a,1x,"=",1x,a    )') 'high_accuracy_default        ','T'
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'accuracy_boost               ',1
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'l_accuracy_boost             ',1
             write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'l_sample_boost               ',1
             close(cambParameterFile)
             ! Run CAMB.
             call System_Command_Do(cambPath//"camb "//parameterFile)
             ! Read the CAMB transfer function file.
             if (allocated(wavenumbers      )) deallocate(wavenumbers      )
             if (allocated(transferFunctions)) deallocate(transferFunctions)
             allocate(wavenumbers      (0    ))
             allocate(transferFunctions(0,0,0))
             do j=1,countRedshiftsUnique
                transferFileName=outputRoot//'_transfer_'//trim(adjustl(redshiftLabelsCombined(j)))//'.dat'
                if (j == 1) then
                   countWavenumber=Count_Lines_In_File(transferFileName,"#")
                   if (allocated(wavenumbers      )) deallocate(wavenumbers      )
                   if (allocated(transferFunctions)) deallocate(transferFunctions)
                   allocate(wavenumbers      (countWavenumber                       ))
                   allocate(transferFunctions(countWavenumber,2,countRedshiftsUnique))
                end if
                open(newunit=cambTransferFile,file=char(transferFileName),status='old',form='formatted')
                i=0
                do while (i < countWavenumber)
                   read (cambTransferFile,'(a)',iostat=status) cambTransferLine
                   if (status == 0) then
                      if (cambTransferLine(1:1) /= "#") then
                         i=i+1
                         read (cambTransferLine,*) wavenumbers(i),transferFunctions(i,cambSpeciesDarkMatter%ID,j),transferFunctions(i,cambSpeciesBaryons%ID,j)
                      end if
                   else
                      call Error_Report('unable to read CAMB transfer function file'//{introspection:location})
                   end if
                end do
                close(cambTransferFile)
             end do
             ! Remove temporary files.
             call File_Remove(parameterFile            )
             call File_Remove(outputRoot//'_params.ini')
             do i=1,countRedshiftsUnique
                call File_Remove(outputRoot//'_transfer_'   //trim(adjustl(redshiftLabelsCombined(i)))//'.dat')
                call File_Remove(outputRoot//'_matterpower_'//trim(adjustl(redshiftLabelsCombined(i)))//'.dat')
             end do
             ! Convert from CAMB units to Galacticus units.
             wavenumbers=+wavenumbers                                                   &
                  &      *cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)
             ! Construct the output HDF5 file.
             !$ call hdf5Access  %set     (                                          )
             call    cambOutput  %openFile(char(fileName_),objectsOverwritable=.true.)
             call    cambOutput  %writeAttribute('Transfer functions created by CAMB.','description')
             call    cambOutput  %writeAttribute(cambFormatVersionCurrent,'fileFormat')
             call    cambOutput  %writeDataset(wavenumbers    ,'wavenumber'                                  ,chunkSize=chunkSize,appendTo=.not.  cambOutput%hasDataset('wavenumber'))
             speciesGroup=cambOutput%openGroup('darkMatter','Group containing transfer functions for dark matter.')
             do i=1,countRedshiftsUnique
                datasetName='transferFunctionZ'//trim(adjustl(redshiftLabelsCombined(i)))
                call speciesGroup%writeDataset(transferFunctions(:,cambSpeciesDarkMatter%ID,i),datasetName,chunkSize=chunkSize,appendTo=.not.speciesGroup%hasDataset(datasetName ))
             end do
             call speciesGroup%close()
             speciesGroup=cambOutput%openGroup('baryons'   ,'Group containing transfer functions for baryons.'    )
             do i=1,countRedshiftsUnique
                datasetName='transferFunctionZ'//trim(adjustl(redshiftLabelsCombined(i)))
                call speciesGroup%writeDataset(transferFunctions(:,cambSpeciesBaryons   %ID,i),datasetName,chunkSize=chunkSize,appendTo=.not.speciesGroup%hasDataset(datasetName ))
             end do
             call speciesGroup%close()
             parametersGroup=cambOutput%openGroup('parameters')
             call parametersGroup%writeAttribute(cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
             call parametersGroup%writeAttribute(cosmologyParameters_%OmegaBaryon    (),'OmegaBaryon'    )
             call parametersGroup%writeAttribute(cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
             call parametersGroup%writeAttribute(cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
             call parametersGroup%writeAttribute(cosmologyParameters_%temperatureCMB (),'temperatureCMB' )
             call parametersGroup%close()
             extrapolationGroup          =cambOutput        %openGroup('extrapolation')
             extrapolationWavenumberGroup=extrapolationGroup%openGroup('wavenumber'   )
             call    extrapolationWavenumberGroup%writeAttribute('extrapolate','low' )
             call    extrapolationWavenumberGroup%writeAttribute('extrapolate','high')
             call    extrapolationWavenumberGroup%close()
             call    extrapolationGroup          %close()
             call    cambOutput                  %close()
             !$ call hdf5Access                  %unset()
          end if
          ! If necessary, construct tables of transfer functions.
          if (present(transferFunctionDarkMatter)) then
             !$ call hdf5Access%set()
             call cambOutput%openFile(char(fileName_))
             call cambOutput%readDataset('wavenumber',wavenumbersLogarithmic)
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
             speciesGroup=cambOutput%openGroup('darkMatter')
             do i=1,size(redshifts)
                datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(i)))
                call speciesGroup%readDataset(datasetName,transferFunctionLogarithmic)
                transferFunctionLogarithmic=log(transferFunctionLogarithmic)
                call transferFunctionDarkMatter%populate(transferFunctionLogarithmic,table=int(redshiftRanks(i)))
                deallocate(transferFunctionLogarithmic)
             end do
             call    speciesGroup%close()
             call    cambOutput  %close()
             !$ call hdf5Access  %unset()
          end if
          if (present(transferFunctionBaryons)) then
             !$ call hdf5Access%set()
             call    cambOutput%openFile(char(fileName_))
             call    cambOutput%readDataset('wavenumber',wavenumbersLogarithmic)
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
             speciesGroup=cambOutput%openGroup('baryons')
             do i=1,size(redshifts)
                datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(i)))
                call speciesGroup%readDataset(datasetName,transferFunctionLogarithmic)
                transferFunctionLogarithmic=log(transferFunctionLogarithmic)
                call transferFunctionBaryons   %populate(transferFunctionLogarithmic,table=int(redshiftRanks(i)))
                deallocate(transferFunctionLogarithmic)
             end do
             call    speciesGroup%close()
             call    cambOutput  %close()
             !$ call hdf5Access  %unset()
          end if
       end if
       ! Unlock the file.
       call File_Unlock(fileLock,sync=iLock == 2)
       ! If we have reached this point on the first pass, the table is sufficient and we are done.
       if (iLock == 1) exit
    end do
    return
  end subroutine Interface_CAMB_Transfer_Function

end module Interfaces_CAMB
