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

!+    Contributions to this file made by: Andrew Benson, Xiaolong Du.

!!{
Contains a module which provides various interfaces to the \gls{axioncamb} code.
!!}

module Interfaces_AxionCAMB
  !!{
  Provides various interfaces to the \gls{axioncamb} code.
  !!}
  use :: File_Utilities, only : lockDescriptor
  private
  public :: Interface_AxionCAMB_Initialize, Interface_AxionCAMB_Transfer_Function

  ! Current file format version for transfer function files. Note that this file format matches that used by the "file" transfer
  ! function class.
  integer                         , parameter :: axionCambFormatVersionCurrent       =     2

  ! Default maximum wavenumber to tabulate.
  double precision                , parameter :: axionCambLogWavenumberMaximumDefault=log(2500.0d0)

  !![
  <enumeration>
   <name>axionCambSpecies</name>
   <description>Particle species in AxionCAMB.</description>
   <visibility>public</visibility>
   <indexing>1</indexing>
   <entry label="darkMatter"     />
   <entry label="coldDarkMatter" />
   <entry label="fuzzyDarkMatter"/>
   <entry label="baryons"        />
  </enumeration>
  !!]

  ! Generate a source digest.
  !![
  <sourceDigest name="axionCambSourceDigest"/>
  !!]

contains

  subroutine Interface_AxionCAMB_Initialize(axionCambPath,axionCambVersion,static)
    !!{
    Initialize the interface with AxionCAMB, including downloading and compiling AxionCAMB if necessary.
    !!}
    use :: File_Utilities    , only : File_Exists      , File_Lock            , File_Unlock , lockDescriptor, &
         &                            Directory_Make
    use :: Display           , only : displayMessage   , verbosityLevelWorking
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath        , pathTypeDataDynamic
    use :: ISO_Varying_String, only : assignment(=)    , char                 , operator(//), replace       , &
          &                           varying_string
    use :: String_Handling   , only : stringSubstitute
    use :: System_Command    , only : System_Command_Do
    use :: System_Compilers  , only : compiler         , compilerOptions      , languageFortran
    implicit none
    type   (varying_string), intent(  out)           :: axionCambPath, axionCambVersion
    logical                , intent(in   ), optional :: static
    integer                                          :: status
    type   (varying_string)                          :: command
    type   (lockDescriptor)                          :: fileLock
    type   (varying_string)                          :: lockPath
    !![
    <optionalArgument name="static" defaultsTo=".false." />
    !!]

    ! Set path and version
    axionCambPath   =inputPath(pathTypeDataDynamic)//"AxionCAMB/"
    lockPath        =inputPath(pathTypeDataDynamic)//"axion_camb"
    axionCambVersion="?"
    call File_Lock(char(lockPath),fileLock,lockIsShared=.false.)
    ! Build the AxionCAMB code.
    if (.not.File_Exists(axionCambPath//"camb")) then
       if (.not.File_Exists(axionCambPath)) then
          ! Download AxionCAMB if necessary.
          call displayMessage("downloading AxionCAMB code....",verbosityLevelWorking)
          call System_Command_Do("git clone https://github.com/dgrin1/axionCAMB.git "//axionCambPath,status)
          if (status /= 0 .or. .not.File_Exists(axionCambPath)) call Error_Report("unable to download AxionCAMB"//{introspection:location})
       end if
       call displayMessage("compiling AxionCAMB code",verbosityLevelWorking)
       command='cd '//axionCambPath//'; sed -E -i~ s/"Ini_Read_Double\('//"'"//'omega_axion'//"'"//'\)\/\(P%H0\/100\)\*\*2"/"Ini_Read_Double\('//"'"//'omega_axion'//"'"//'\)"/ inidriver_axion.F90; sed -E -i~ s/"^F90C[[:space:]]*=[[:space:]]*[[:alpha:]]+"/"F90C = '//compiler(languageFortran)//'"/ Makefile; sed -E -i~ s/"^FFLAGS[[:space:]]*\+=[[:space:]]*\-march=native"/"FFLAGS+="/ Makefile; sed -E -i~ s/"^FFLAGS[[:space:]]*=[[:space:]]*.*"/"FFLAGS = -O3 '//stringSubstitute(compilerOptions(languageFortran),"/","\/")
       if (static_) command=command//" -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive"
       command=command//'"/ Makefile; find . -name "*.f90" | xargs sed -E -i~ s/"error stop"/"error stop "/; make -j1 camb'
       call System_Command_Do(char(command),status);
       if (status /= 0 .or. .not.File_Exists(axionCambPath//"camb")) call Error_Report("failed to build AxionCAMB code"//{introspection:location})
    end if
    call File_Unlock(fileLock)
    return
  end subroutine Interface_AxionCAMB_Initialize

  subroutine Interface_AxionCAMB_Transfer_Function(cosmologyParameters_,darkMatterParticle_,redshifts,wavenumberRequired,wavenumberMaximum,countPerDecade,fileName,wavenumberMaximumReached,transferFunctionDarkMatter,transferFunctionColdDarkMatter,transferFunctionFuzzyDarkMatter,transferFunctionBaryons)
    !!{
    Run AxionCAMB as necessary to compute transfer functions.
    !!}
    use               :: Cosmology_Parameters            , only : cosmologyParametersClass    , hubbleUnitsLittleH
    use               :: Dark_Matter_Particles           , only : darkMatterParticleClass     , darkMatterParticleFuzzyDarkMatter
    use               :: File_Utilities                  , only : Count_Lines_In_File         , Directory_Make                   , File_Exists , File_Lock     , &
          &                                                       File_Path                   , File_Remove                      , File_Unlock , lockDescriptor, &
          &                                                       File_Name_Temporary
    use               :: Error                           , only : Error_Report
    use               :: Input_Paths                     , only : inputPath                   , pathTypeDataDynamic
    use               :: HDF5                            , only : hsize_t
    use               :: Hashes_Cryptographic            , only : Hash_MD5
    use               :: HDF5_Access                     , only : hdf5Access
    use               :: IO_HDF5                         , only : hdf5Object
    use   , intrinsic :: ISO_C_Binding                   , only : c_size_t
    use               :: ISO_Varying_String              , only : assignment(=)               , char                             , extract     , len           , &
          &                                                       operator(==)                , varying_string                   , operator(//)
    use               :: Input_Parameters                , only : inputParameters
    use               :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial
    use               :: Numerical_Constants_Prefixes    , only : kilo
    use               :: Numerical_Interpolation         , only : GSL_Interp_cSpline
    !$ use            :: OMP_Lib                         , only : OMP_Get_Thread_Num
    use               :: Sorting                         , only : sortIndex
    use               :: String_Handling                 , only : operator(//)                , String_C_To_Fortran
    use               :: System_Command                  , only : System_Command_Do
    use               :: Table_Labels                    , only : extrapolationTypeExtrapolate, extrapolationTypeFix
    use               :: Tables                          , only : table                       , table1DGeneric
    implicit none
    class           (cosmologyParametersClass), intent(inout)                   :: cosmologyParameters_
    class           (darkMatterParticleClass ), intent(inout)                   :: darkMatterParticle_
    double precision                          , intent(in   ), dimension(:    ) :: redshifts
    double precision                          , intent(in   )                   :: wavenumberRequired                      , wavenumberMaximum
    integer                                   , intent(in   ), optional         :: countPerDecade
    type            (varying_string          ), intent(  out), optional         :: fileName
    type            (table1DGeneric          ), intent(  out), optional         :: transferFunctionDarkMatter              , transferFunctionColdDarkMatter, &
         &                                                                         transferFunctionFuzzyDarkMatter         , transferFunctionBaryons
    logical                                   , intent(inout), optional         :: wavenumberMaximumReached
    double precision                          , allocatable  , dimension(:    ) :: wavenumbers                             , wavenumbersLogarithmic        , &
         &                                                                         transferFunctionLogarithmic             , redshiftsCombined
    double precision                          , allocatable  , dimension(:,:,:) :: transferFunctions
    character       (len= 9                  ), allocatable  , dimension(:    ) :: redshiftLabels                          , redshiftLabelsCombined
    integer         (c_size_t                ), allocatable  , dimension(:    ) :: redshiftRanks                           , redshiftRanksCombined
    type            (varying_string          ), allocatable  , dimension(:    ) :: datasetNames
    integer         (hsize_t                 ), parameter                       :: chunkSize                   =100_hsize_t
    type            (lockDescriptor          )                                  :: fileLock
    character       (len=255                 )                                  :: axionCambTransferLine
    type            (varying_string          )                                  :: axionCambVersion                        , parameterFile                 , &
         &                                                                         axionCambPath                           , outputRoot
    double precision                                                            :: wavenumberAxionCAMB                     , coldDarkMatterDensityFraction , &
         &                                                                         fuzzyDarkMatterDensityFraction
    double precision                                                            :: transferFunctionUnused
    integer                                                                     :: status                                  , axionCambParameterFile        , &
         &                                                                         i                                       , axionCambTransferFile         , &
         &                                                                         j                                       , countRedshiftsUnique
    integer         (c_size_t                )                                  :: countWavenumber
    type            (hdf5Object              )                                  :: axionCambOutput                         , parametersGroup               , &
         &                                                                         extrapolationWavenumberGroup            , extrapolationGroup            , &
         &                                                                         speciesGroup
    character       (len=32                  )                                  :: parameterLabel                          , datasetName                   , &
         &                                                                         redshiftLabel                           , indexLabel
    type            (varying_string          )                                  :: uniqueLabel                             , workPath                      , &
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
    select type (darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       ! Add fuzzy dark matter to descriptor.
       write (parameterLabel,'(e12.6)') darkMatterParticle_%mass           ()
       call descriptor%addParameter("fuzzyDMMass"           ,parameterLabel)
       ! Add fuzzy dark matter density fraction to descriptor.
       write (parameterLabel,'(e12.6)') darkMatterParticle_%densityFraction()
       call descriptor%addParameter("fuzzyDMDensityFraction",parameterLabel)
    class default
       call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
    end select
    ! Add the unique label string to the descriptor.
    uniqueLabel=descriptor%serializeToString()            // &
         &      "_sourceDigest:"                          // &
         &      String_C_To_Fortran(axionCambSourceDigest)
    call descriptor%destroy()
    ! Build the file name.
    fileName_=char(inputPath(pathTypeDataDynamic))                            // &
         &                  'largeScaleStructure/transfer_function_AxionCAMB_'// &
         &                  Hash_MD5(uniqueLabel)                             // &
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
       call    axionCambOutput%openFile(char(fileName_))
       call    axionCambOutput%readDataset           ('wavenumber'  ,wavenumbers                                    )
       allocate(transferFunctions(size(wavenumbers),4,size(redshifts)))
       speciesGroup=axionCambOutput%openGroup('darkMatter')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(redshiftRanks(i))))
          if (speciesGroup%hasDataset(datasetName)) then
             call speciesGroup%readDatasetStatic(datasetName,transferFunctions(:,axionCambSpeciesDarkMatter     %ID,i))
          else
             allEpochsFound=.false.
          end if
       end do
       call speciesGroup%close()
       speciesGroup=axionCambOutput%openGroup('coldDarkMatter')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(redshiftRanks(i))))
          if (speciesGroup%hasDataset(datasetName)) then
             call speciesGroup%readDatasetStatic(datasetName,transferFunctions(:,axionCambSpeciesColdDarkMatter %ID,i))
          else
             allEpochsFound=.false.
          end if
       end do
       call speciesGroup%close()
       speciesGroup=axionCambOutput%openGroup('fuzzyDarkMatter')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(redshiftRanks(i))))
          if (speciesGroup%hasDataset(datasetName)) then
             call speciesGroup%readDatasetStatic(datasetName,transferFunctions(:,axionCambSpeciesFuzzyDarkMatter%ID,i))
          else
             allEpochsFound=.false.
          end if
       end do
       call speciesGroup%close()
       speciesGroup=axionCambOutput%openGroup('baryons')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(redshiftRanks(i))))
          if (speciesGroup%hasDataset(datasetName)) then
             call speciesGroup%readDatasetStatic(datasetName,transferFunctions(:,axionCambSpeciesBaryons        %ID,i))
          else
             allEpochsFound=.false.
          end if
       end do
       call   speciesGroup   %close()
       call   axionCambOutput%close()
       !$ call hdf5Access    %unset()
    end if
    ! Get density fractions of cold dark matter and fuzzy dark matter.
    coldDarkMatterDensityFraction =cosmologyParameters_%OmegaMatter()-cosmologyParameters_%OmegaBaryon()
    fuzzyDarkMatterDensityFraction=0.0d0
    select type (darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       coldDarkMatterDensityFraction =(1.0d0-darkMatterParticle_%densityFraction())*(cosmologyParameters_%OmegaMatter()-cosmologyParameters_%OmegaBaryon())
       fuzzyDarkMatterDensityFraction=       darkMatterParticle_%densityFraction() *(cosmologyParameters_%OmegaMatter()-cosmologyParameters_%OmegaBaryon())
    class default
       call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
    end select
    if (fuzzyDarkMatterDensityFraction == 0.0d0) call Error_Report('density fraction of fuzzy dark matter can not be exactally 0'//{introspection:location})
    if (.not.allocated(wavenumbers) .or. wavenumberRequired > wavenumbers(size(wavenumbers)) .or. .not.allEpochsFound) then
       ! If the wavenumber if out of range, or if not all requested epochs exist within the file, recompute the AxionCAMB transfer function.
       ! Find all existing epochs in the file, create a union of these and the requested epochs.
       if (File_Exists(fileName_)) then
          !$ call hdf5Access     %set       (               )
          call    axionCambOutput%openFile  (char(fileName_))
          speciesGroup=axionCambOutput%openGroup('darkMatter')
          call    speciesGroup   %datasets  (datasetNames   )
          call    speciesGroup   %close     (               )
          call    axionCambOutput%close     (               )
          !$ call hdf5Access     %unset     (               )
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
       ! Ensure AxionCAMB is initialized.
       call Interface_AxionCAMB_Initialize(axionCambPath,axionCambVersion)
       ! Determine maximum wavenumber.
       wavenumberAxionCAMB=exp(max(log(wavenumberRequired)+1.0d0,axionCambLogWavenumberMaximumDefault))
       if (wavenumberAxionCAMB > wavenumberMaximum) then
          wavenumberAxionCAMB=wavenumberMaximum
          if (present(wavenumberMaximumReached)) wavenumberMaximumReached=.true.
       end if
       if (allocated(wavenumbers)) wavenumberAxionCAMB=max(wavenumberAxionCAMB,wavenumbers(size(wavenumbers)))
       ! Construct input file for AxionCAMB.
       workPath     =inputPath(pathTypeDataDynamic)//'largeScaleStructure/'
       parameterFile=File_Name_Temporary('transfer_function_parameters',char(workPath))//'.txt'
       outputRoot   =File_Name_Temporary('axionCAMB'                   ,char(workPath))
       !$ parameterFile=parameterFile//'_'//OMP_Get_Thread_Num()
       parameterFile=parameterFile//'.txt'
       call Directory_Make(workPath)
       open(newunit=axionCambParameterFile,file=char(parameterFile),status='unknown',form='formatted')
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'output_root                  ',char(outputRoot)
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'get_scalar_cls               ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'get_vector_cls               ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'get_tensor_cls               ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'get_transfer                 ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'do_lensing                   ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'do_nonlinear                 ',0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'l_max_scalar                 ',2200.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'l_max_tensor                 ',1500.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'k_eta_max_tensor             ',3000.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'use_physical                 ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'use_axfrac                   ','F'
       select type (darkMatterParticle_)
       class is (darkMatterParticleFuzzyDarkMatter)
          write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'm_ax                      ',darkMatterParticle_ %mass           ()*kilo
          write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_axion               ',fuzzyDarkMatterDensityFraction
       class default
          call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
       end select
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_baryon                 ',cosmologyParameters_%OmegaBaryon    ()
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_cdm                    ',coldDarkMatterDensityFraction
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_lambda                 ',cosmologyParameters_%OmegaDarkEnergy()
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_neutrino               ',0.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'omk                          ',0.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'hubble                       ',cosmologyParameters_%HubbleConstant ()
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'w                            ',-1.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'cs2_lam                      ',1.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'temp_cmb                     ',cosmologyParameters_%temperatureCMB ()
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'helium_fraction              ',heliumByMassPrimordial
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'massless_neutrinos           ',2.046d0
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'nu_mass_eigenstates          ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'massive_neutrinos            ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'share_delta_neff             ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'nu_mass_fractions            ',1.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x      )') 'nu_mass_degeneracies         '
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'initial_power_num            ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'pivot_scalar                 ',0.05d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'pivot_tensor                 ',0.05d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'scalar_amp(1)                ',2.1d-9
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'scalar_spectral_index(1)     ',0.96d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'scalar_nrun(1)               ',0.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'tensor_spectral_index(1)     ',0.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'initial_ratio(1)             ',1.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'reionization                 ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 're_use_optical_depth         ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 're_optical_depth             ',0.09d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 're_redshift                  ',11.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 're_delta_redshift            ',1.5d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 're_ionization_frac           ',-1.0d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'RECFAST_fudge                ',1.14d0
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'RECFAST_fudge_He             ',0.86d0
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'RECFAST_Heswitch             ',6
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'RECFAST_Hswitch              ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'initial_condition            ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,5(i2))') 'initial_vector               ',-1,0,0,0,0
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'vector_mode                  ',0
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'COBE_normalize               ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'CMB_outputscale              ',7.42835025d12
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'transfer_high_precision      ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,e12.6)') 'transfer_kmax                ',wavenumberAxionCAMB/cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)
       write (axionCambParameterFile,'(a,1x,"=",1x,i5   )') 'transfer_k_per_logint        ',countPerDecade_
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'transfer_num_redshifts       ',countRedshiftsUnique
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'transfer_interp_matterpower  ','T'
       do i=countRedshiftsUnique,1,-1
          write (indexLabel,'(i4)') countRedshiftsUnique+1-i
          write (axionCambParameterFile,'(a,a,a,1x,"=",1x,a    )') 'transfer_redshift('   ,trim(adjustl(indexLabel)),')'               ,trim(adjustl(redshiftLabelsCombined(i)))
          write (axionCambParameterFile,'(a,a,a,1x,"=",1x,a,a,a)') 'transfer_filename('   ,trim(adjustl(indexLabel)),')','transfer_'   ,trim(adjustl(redshiftLabelsCombined(i))),'.dat'
          write (axionCambParameterFile,'(a,a,a,1x,"=",1x,a,a,a)') 'transfer_matterpower(',trim(adjustl(indexLabel)),')','matterpower_',trim(adjustl(redshiftLabelsCombined(i))),'.dat'
       end do
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'scalar_output_file           ','scalCls.dat'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'vector_output_file           ','vecCls.dat'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'tensor_output_file           ','tensCls.dat'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'total_output_file            ','totCls.dat'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'lensed_output_file           ','lensedCls.dat'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'lensed_total_output_file     ','lensedtotCls.dat'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'lens_potential_output_file   ','lenspotentialCls.dat'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'FITS_filename                ','scalCls.fits'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'do_lensing_bispectrum        ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'do_primordial_bispectrum     ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_nfields           ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_slice_base_L      ',0
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_ndelta            ',3
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_delta(1)          ',0
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_delta(2)          ',2
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_delta(3)          ',4
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'bispectrum_do_fisher         ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_fisher_noise      ',0
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_fisher_noise_pol  ',0
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'bispectrum_fisher_fwhm_arcmin',7
       write (axionCambParameterFile,'(a,1x,"=",1x      )') 'bispectrum_full_output_file  '
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'bispectrum_full_output_sparse','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'bispectrum_export_alpha_beta ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'feedback_level               ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'derived_parameters           ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'lensing_method               ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'accurate_BB                  ','F'
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'massive_nu_approx            ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'accurate_polarization        ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'accurate_reionization        ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'do_tensor_neutrinos          ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'do_late_rad_truncation       ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'number_of_threads            ',0
       write (axionCambParameterFile,'(a,1x,"=",1x,a    )') 'high_accuracy_default        ','T'
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'accuracy_boost               ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'l_accuracy_boost             ',1
       write (axionCambParameterFile,'(a,1x,"=",1x,i1   )') 'l_sample_boost               ',1
       close(axionCambParameterFile)
       ! Run AxionCAMB.
       call System_Command_Do(axionCambPath//"camb "//parameterFile)
       ! Read the AxionCAMB transfer function file.
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
             allocate(transferFunctions(countWavenumber,4,countRedshiftsUnique))
          end if
          open(newunit=axionCambTransferFile,file=char(transferFileName),status='old',form='formatted')
          i=0
          do while (i < countWavenumber)
             read (axionCambTransferFile,'(a)',iostat=status) axionCambTransferLine
             if (status == 0) then
                if (axionCambTransferLine(1:1) /= "#") then
                   i=i+1
                   read (axionCambTransferLine,*) wavenumbers           (i                                     ),transferFunctions(i,axionCambSpeciesColdDarkMatter%ID,j), &
                        &                         transferFunctions     (i,axionCambSpeciesBaryons        %ID,j),transferFunctionUnused                                  , &
                        &                         transferFunctionUnused                                        ,transferFunctionUnused                                  , &
                        &                         transferFunctions     (i,axionCambSpeciesFuzzyDarkMatter%ID,j)
                   ! Transfer function for total dark matter perturbations.
                   transferFunctions(i,axionCambSpeciesDarkMatter%ID,j)=(                                                                                          &
                        &                                                +transferFunctions(i,axionCambSpeciesColdDarkMatter %ID,j)*coldDarkMatterDensityFraction  &
                        &                                                +transferFunctions(i,axionCambSpeciesFuzzyDarkMatter%ID,j)*fuzzyDarkMatterDensityFraction &
                        &                                               )                                                                                          &
                        &                                               /(coldDarkMatterDensityFraction+fuzzyDarkMatterDensityFraction)
                end if
             else
                call Error_Report('unable to read AxionCAMB transfer function file'//{introspection:location})
             end if
          end do
          close(axionCambTransferFile)
       end do
       ! Make sure that the transfer function is non-negative.
       transferFunctions=abs(transferFunctions)
       ! Remove temporary files.
       call File_Remove(parameterFile            )
       call File_Remove(outputRoot//'_params.ini')
       do i=1,countRedshiftsUnique
          call File_Remove(outputRoot//'_transfer_'   //trim(adjustl(redshiftLabelsCombined(i)))//'.dat')
          call File_Remove(outputRoot//'_matterpower_'//trim(adjustl(redshiftLabelsCombined(i)))//'.dat')
       end do
       ! Convert from AxionCAMB units to Galacticus units.
       wavenumbers=+wavenumbers                                                   &
            &      *cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)
       ! Construct the output HDF5 file.
       !$ call hdf5Access       %set     (                                          )
       call    axionCambOutput  %openFile(char(fileName_),objectsOverwritable=.true.)
       call    axionCambOutput  %writeAttribute('Transfer functions created by AxionCAMB.','description')
       call    axionCambOutput  %writeAttribute(axionCambFormatVersionCurrent,'fileFormat')
       call    axionCambOutput  %writeDataset(wavenumbers    ,'wavenumber'                                 ,chunkSize=chunkSize,appendTo=.not.axionCambOutput%hasDataset('wavenumber'))
       speciesGroup=axionCambOutput%openGroup('darkMatter','Group containing transfer functions for dark matter.')
       do i=1,countRedshiftsUnique
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabelsCombined(i)))
          call speciesGroup%writeDataset(transferFunctions(:,axionCambSpeciesDarkMatter     %ID,i),datasetName,chunkSize=chunkSize,appendTo=.not.speciesGroup%hasDataset(datasetName ))
       end do
       call speciesGroup%close()
       speciesGroup=axionCambOutput%openGroup('coldDarkMatter' ,'Group containing transfer functions for cold dark matter.' )
       do i=1,countRedshiftsUnique
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabelsCombined(i)))
          call speciesGroup%writeDataset(transferFunctions(:,axionCambSpeciesColdDarkMatter %ID,i),datasetName,chunkSize=chunkSize,appendTo=.not.speciesGroup%hasDataset(datasetName ))
       end do
       call speciesGroup%close()
       speciesGroup=axionCambOutput%openGroup('fuzzyDarkMatter','Group containing transfer functions for fuzzy dark matter.')
       do i=1,countRedshiftsUnique
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabelsCombined(i)))
          call speciesGroup%writeDataset(transferFunctions(:,axionCambSpeciesFuzzyDarkMatter%ID,i),datasetName,chunkSize=chunkSize,appendTo=.not.speciesGroup%hasDataset(datasetName ))
       end do
       call speciesGroup%close()
       speciesGroup=axionCambOutput%openGroup('baryons'        ,'Group containing transfer functions for baryons.'          )
       do i=1,countRedshiftsUnique
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabelsCombined(i)))
          call speciesGroup%writeDataset(transferFunctions(:,axionCambSpeciesBaryons        %ID,i),datasetName,chunkSize=chunkSize,appendTo=.not.speciesGroup%hasDataset(datasetName ))
       end do
       call speciesGroup%close()
       parametersGroup=axionCambOutput%openGroup('parameters')
       call parametersGroup%writeAttribute(cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
       call parametersGroup%writeAttribute(cosmologyParameters_%OmegaBaryon    (),'OmegaBaryon'    )
       call parametersGroup%writeAttribute(cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
       call parametersGroup%writeAttribute(cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
       call parametersGroup%writeAttribute(cosmologyParameters_%temperatureCMB (),'temperatureCMB' )
       select type (darkMatterParticle_)
       class is (darkMatterParticleFuzzyDarkMatter)
          call parametersGroup%writeAttribute(darkMatterParticle_%mass           (),'fuzzyDMMass'           )
          call parametersGroup%writeAttribute(darkMatterParticle_%densityFraction(),'fuzzyDMDensityFraction')
       class default
          call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
       end select
       call parametersGroup%close()
       extrapolationGroup          =axionCambOutput   %openGroup('extrapolation')
       extrapolationWavenumberGroup=extrapolationGroup%openGroup('wavenumber'   )
       call    extrapolationWavenumberGroup%writeAttribute('extrapolate','low' )
       call    extrapolationWavenumberGroup%writeAttribute('fix'        ,'high')
       call    extrapolationWavenumberGroup%close()
       call    extrapolationGroup          %close()
       call    axionCambOutput             %close()
       !$ call hdf5Access                  %unset()
    end if
    ! If necessary, construct tables of transfer functions.
    if (present(transferFunctionDarkMatter)) then
       !$ call hdf5Access%set()
       call axionCambOutput%openFile(char(fileName_))
       call axionCambOutput%readDataset('wavenumber',wavenumbersLogarithmic)
       wavenumbersLogarithmic=log(wavenumbersLogarithmic)
       call transferFunctionDarkMatter     %create(                                                 &
            &                                                        wavenumbersLogarithmic       , &
            &                                      tableCount       =size(redshifts)              , &
            &                                      extrapolationType=[                              &
            &                                                         extrapolationTypeExtrapolate, &
            &                                                         extrapolationTypeFix          &
            &                                                        ]                            , &
            &                                      interpolationType=GSL_Interp_cSpline             &
            &                                     )
       deallocate(wavenumbersLogarithmic)
       speciesGroup=axionCambOutput%openGroup('darkMatter')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(i)))
          call speciesGroup%readDataset(datasetName,transferFunctionLogarithmic)
          transferFunctionLogarithmic=log(abs(transferFunctionLogarithmic))
          call transferFunctionDarkMatter     %populate(transferFunctionLogarithmic,table=int(redshiftRanks(i)))
          deallocate(transferFunctionLogarithmic)
       end do
       call    speciesGroup   %close()
       call    axionCambOutput%close()
       !$ call hdf5Access     %unset()
    end if
    if (present(transferFunctionColdDarkMatter)) then
       !$ call hdf5Access%set()
       call axionCambOutput%openFile(char(fileName_))
       call axionCambOutput%readDataset('wavenumber',wavenumbersLogarithmic)
       wavenumbersLogarithmic=log(wavenumbersLogarithmic)
       call transferFunctionColdDarkMatter %create(                                                 &
            &                                                        wavenumbersLogarithmic       , &
            &                                      tableCount       =size(redshifts)              , &
            &                                      extrapolationType=[                              &
            &                                                         extrapolationTypeExtrapolate, &
            &                                                         extrapolationTypeFix          &
            &                                                        ]                            , &
            &                                      interpolationType=GSL_Interp_cSpline             &
            &                                     )
       deallocate(wavenumbersLogarithmic)
       speciesGroup=axionCambOutput%openGroup('coldDarkMatter')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(i)))
          call speciesGroup%readDataset(datasetName,transferFunctionLogarithmic)
          transferFunctionLogarithmic=log(abs(transferFunctionLogarithmic))
          call transferFunctionColdDarkMatter %populate(transferFunctionLogarithmic,table=int(redshiftRanks(i)))
          deallocate(transferFunctionLogarithmic)
       end do
       call    speciesGroup   %close()
       call    axionCambOutput%close()
       !$ call hdf5Access     %unset()
    end if
    if (present(transferFunctionFuzzyDarkMatter)) then
       !$ call hdf5Access%set()
       call axionCambOutput%openFile(char(fileName_))
       call axionCambOutput%readDataset('wavenumber',wavenumbersLogarithmic)
       wavenumbersLogarithmic=log(wavenumbersLogarithmic)
       call transferFunctionFuzzyDarkMatter%create(                                                 &
            &                                                        wavenumbersLogarithmic       , &
            &                                      tableCount       =size(redshifts)              , &
            &                                      extrapolationType=[                              &
            &                                                         extrapolationTypeExtrapolate, &
            &                                                         extrapolationTypeFix          &
            &                                                        ]                            , &
            &                                      interpolationType=GSL_Interp_cSpline             &
            &                                     )
       deallocate(wavenumbersLogarithmic)
       speciesGroup=axionCambOutput%openGroup('fuzzyDarkMatter')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(i)))
          call speciesGroup%readDataset(datasetName,transferFunctionLogarithmic)
          transferFunctionLogarithmic=log(abs(transferFunctionLogarithmic))
          call transferFunctionFuzzyDarkMatter%populate(transferFunctionLogarithmic,table=int(redshiftRanks(i)))
          deallocate(transferFunctionLogarithmic)
       end do
       call    speciesGroup   %close()
       call    axionCambOutput%close()
       !$ call hdf5Access     %unset()
    end if
    if (present(transferFunctionBaryons)) then
       !$ call hdf5Access%set()
       call    axionCambOutput%openFile(char(fileName_))
       call    axionCambOutput%readDataset('wavenumber',wavenumbersLogarithmic)
       wavenumbersLogarithmic=log(wavenumbersLogarithmic)
       call transferFunctionBaryons        %create(                                                 &
            &                                                        wavenumbersLogarithmic       , &
            &                                      tableCount       =size(redshifts)              , &
            &                                      extrapolationType=[                              &
            &                                                         extrapolationTypeExtrapolate, &
            &                                                         extrapolationTypeFix          &
            &                                                        ]                            , &
            &                                      interpolationType=GSL_Interp_cSpline             &
            &                                     )
       deallocate(wavenumbersLogarithmic)
       speciesGroup=axionCambOutput%openGroup('baryons')
       do i=1,size(redshifts)
          datasetName='transferFunctionZ'//trim(adjustl(redshiftLabels(i)))
          call speciesGroup%readDataset(datasetName,transferFunctionLogarithmic)
          transferFunctionLogarithmic=log(abs(transferFunctionLogarithmic))
          call transferFunctionBaryons        %populate(transferFunctionLogarithmic,table=int(redshiftRanks(i)))
          deallocate(transferFunctionLogarithmic)
       end do
       call    speciesGroup   %close()
       call    axionCambOutput%close()
       !$ call hdf5Access     %unset()
    end if
    ! Unlock the file.
    call File_Unlock(fileLock)
    return
  end subroutine Interface_AxionCAMB_Transfer_Function

end module Interfaces_AxionCAMB
