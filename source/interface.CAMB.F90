!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which provides various interfaces to the \gls{camb} code.

module Interfaces_CAMB
  !% Provides various interfaces to the \gls{camb} code.
  use File_Utilities, only : lockDescriptor
  private
  public :: Interface_CAMB_Initialize, Interface_CAMB_Transfer_Function

  ! Global lock descriptor to be used in non-(recent Linux) cases
  type            (lockDescriptor)            :: cambFileLockGlobal
  logical                                     :: cambFileLockInitialized        =.false.

  ! Current file format version for transfer function files.
  integer                         , parameter :: cambFormatVersionCurrent       =     1

  ! Default maximum wavenumber to tabulate.
  double precision                , parameter :: cambLogWavenumberMaximumDefault=log(10.0d0)

  ! Generate a source digest.
  !# <sourceDigest name="cambSourceDigest"/>
  
contains

  subroutine Interface_CAMB_Initialize(cambPath,cambVersion,static)
    !% Initialize the interface with CAMB, including downloading and compiling CAMB if necessary.
    use ISO_Varying_String, only : varying_string            , replace            , operator(//), assignment(=), &
         &                         char
    use Galacticus_Paths  , only : galacticusPath            , pathTypeDataDynamic
    use File_Utilities    , only : File_Exists
    use System_Command    , only : System_Command_Do
    use Galacticus_Display, only : Galacticus_Display_Message, verbosityWorking
    use Galacticus_Error  , only : Galacticus_Error_Report
    implicit none
    type   (varying_string), intent(  out)           :: cambPath, cambVersion
    logical                , intent(in   ), optional :: static
    integer                                          :: status  , flagsLength
    type   (varying_string)                          :: command
    !# <optionalArgument name="static" defaultsTo=".false." />

    ! Set path and version
    cambPath   =galacticusPath(pathTypeDataDynamic)//"CAMB/"
    cambVersion="?"
    ! Build the CAMB code.
    if (.not.File_Exists(cambPath//"camb")) then
       ! Unpack the code.
       if (.not.File_Exists(cambPath)) then
          ! Download CAMB if necessary.
          if (.not.File_Exists(galacticusPath(pathTypeDataDynamic)//"CAMB.tar.gz")) then
             call Galacticus_Display_Message("downloading CAMB code....",verbosityWorking)
             call System_Command_Do("wget http://camb.info/CAMB.tar.gz -O "//galacticusPath(pathTypeDataDynamic)//"CAMB.tar.gz",status)
             if (status /= 0 .or. .not.File_Exists(galacticusPath(pathTypeDataDynamic)//"CAMB.tar.gz")) call Galacticus_Error_Report("unable to download CAMB"//{introspection:location})
          end if
          call Galacticus_Display_Message("unpacking CAMB code....",verbosityWorking)
          call System_Command_Do("tar -x -v -z -C "//galacticusPath(pathTypeDataDynamic)//" -f "//galacticusPath(pathTypeDataDynamic)//"CAMB.tar.gz");
          if (status /= 0 .or. .not.File_Exists(cambPath)) call Galacticus_Error_Report('failed to unpack CAMB code'//{introspection:location})
       end if
       call Galacticus_Display_Message("compiling CAMB code",verbosityWorking)
       command='cd '//cambPath//'; sed -r -i~ s/"ifortErr\s*=.*"/"ifortErr = 1"/ Makefile; sed -r -i~ s/"gfortErr\s*=.*"/"gfortErr = 0"/ Makefile; sed -r -i~ s/"^FFLAGS\s*\+=\s*\-march=native"/"FFLAGS+="/ Makefile; sed -r -i~ s/"^FFLAGS\s*=\s*.*"/"FFLAGS = -Ofast -fopenmp'
       if (static_) then
          ! Include Galacticus compilation flags here - may be necessary for static linking.
          call Get_Environment_Variable("GALACTICUS_FCFLAGS",length=flagsLength,status=status)
          if (status  == 0) command=command//" "//flagsRetrieve(flagsLength)
          command=command//" -static"
       end if
       command=command//'"/ Makefile; find . -name "*.f90" | xargs sed -r -i~ s/"error stop"/"error stop "/; make -j1 camb'
       call System_Command_Do(char(command),status);
       if (status /= 0 .or. .not.File_Exists(cambPath//"camb")) call Galacticus_Error_Report("failed to build CAMB code"//{introspection:location})
    end if
    return

  contains
    
    function flagsRetrieve(flagsLength)
      !% Retrieve the compiler flags.
      implicit none
      type     (varying_string )                :: flagsRetrieve
      integer                   , intent(in   ) :: flagsLength
      character(len=flagsLength)                :: flags

      call Get_Environment_Variable('GALACTICUS_FCFLAGS',value=flags)
      flagsRetrieve=replace(flags,"/","\/",every=.true.)
      return
    end function flagsRetrieve

  end subroutine Interface_CAMB_Initialize

  subroutine Interface_CAMB_Transfer_Function(cosmologyParameters_,wavenumberRequired,wavenumberMaximum,lockFileGlobally,fileName,wavenumberMaximumReached)
    !% Run CAMB as necessary to compute transfer functions.
    !$ use            :: OMP_Lib                         , only : OMP_Get_Thread_Num
    use   , intrinsic :: ISO_C_Binding                   , only : c_size_t
    use               :: Input_Parameters                , only : inputParameters
    use               :: ISO_Varying_String              , only : varying_string          , char
    use               :: IO_HDF5                         , only : hdf5Object
    use               :: File_Utilities                  , only : File_Lock_Initialize    , File_Lock          , File_Unlock, File_Exists, &
         &                                                        Count_Lines_In_File
    use               :: System_Command                  , only : System_Command_Do
    use               :: Galacticus_Error                , only : Galacticus_Error_Report
    use               :: Galacticus_Paths                , only : galacticusPath          , pathTypeDataDynamic
    use               :: String_Handling                 , only : operator(//)
    use               :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial
    use               :: Hashes_Cryptographic            , only : Hash_MD5
    use               :: Cosmology_Parameters            , only : cosmologyParametersClass, hubbleUnitsLittleH
    use               :: HDF5                            , only : hsize_t
    implicit none
    class           (cosmologyParametersClass), intent(inout)               :: cosmologyParameters_
    double precision                          , intent(in   )               :: wavenumberRequired                      , wavenumberMaximum
    logical                                   , intent(in   )               :: lockFileGlobally
    type            (varying_string          ), intent(  out)               :: fileName
    logical                                   , intent(inout)               :: wavenumberMaximumReached
    double precision                          , allocatable  , dimension(:) :: wavenumbers                             , transferFunctions
    integer         (hsize_t                 ), parameter                   :: chunkSize                   =100_hsize_t
    type            (lockDescriptor          )                              :: fileLock
    character       (len=32                  )                              :: wavenumberLabel
    character       (len=255                 )                              :: hostName                                , cambTransferLine
    type            (varying_string          )                              :: command                                 , parameterFile     , &
         &                                                                     cambPath                                , cambVersion
    double precision                                                        :: wavenumberCAMB
    integer                                                                 :: status                                  , cambParameterFile , &
         &                                                                     i                                       , cambTransferFile
    integer         (c_size_t                )                              :: countWavenumber
    type            (hdf5Object              )                              :: cambOutput                              , parametersGroup   , &
         &                                                                     extrapolationWavenumberGroup            , extrapolationGroup
    character       (len=32                  )                              :: parameterLabel
    type            (varying_string          )                              :: uniqueLabel                             , workPath
    type            (inputParameters         )                              :: descriptor
    logical                                                                 :: fileIsNew
    
    ! Get a constructor descriptor for this object.
    descriptor=inputParameters()
    call cosmologyParameters_%descriptor(descriptor)
    ! Add primordial helium abundance to the descriptor.
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    call descriptor%addParameter("Y_He",parameterLabel)
    ! Add the unique label string to the descriptor.
    uniqueLabel=descriptor%serializeToString()// &
         &      "_sourceDigest:"              // &
         &      cambSourceDigest
    call descriptor%destroy()
    ! Build the file name.
    fileName=char(galacticusPath(pathTypeDataDynamic))                       // &
         &                      'largeScaleStructure/transfer_function_CAMB_'// &
         &                      Hash_MD5(uniqueLabel)                        // &
         &                      '.hdf5'
    ! Create the directory.
    call System_Command_Do("mkdir -p `dirname "//fileName//"`")
    ! If the file exists but has not yet been read, read it now.
    if (lockFileGlobally) then
       if (.not.cambFileLockInitialized) then
          !$omp critical (cambFileLockInitialize)
          if (.not.cambFileLockInitialized) then
             call File_Lock_Initialize(cambFileLockGlobal)
             cambFileLockInitialized=.true.
          end if
          !$omp end critical (cambFileLockInitialize)
       end if
    else
       call File_Lock_Initialize(fileLock)
    end if
    if (File_Exists(fileName)) then
       if (lockFileGlobally) then
          call File_Lock(char(fileName),cambFileLockGlobal)
       else
          call File_Lock(char(fileName),fileLock          )
       end if
       call cambOutput%openFile(char(fileName))
       call cambOutput%readDataset('wavenumber'      ,wavenumbers      )
       call cambOutput%readDataset('transferFunction',transferFunctions)
       call cambOutput%close()
       if (lockFileGlobally) then
          call File_Unlock(cambFileLockGlobal)
       else
          call File_Unlock(fileLock          )
       end if
    end if
    if (.not.allocated(wavenumbers) .or. wavenumberRequired > wavenumbers(size(wavenumbers))) then
       ! If the wavenumber if out of range, recompute the CAMB transfer function.
       ! Get a lock on the relevant lock file.
       if (lockFileGlobally) then
          call File_Lock(char(fileName),cambFileLockGlobal)
       else
          call File_Lock(char(fileName),fileLock          )
       end if
       ! Ensure CAMB is initialized.
       call Interface_CAMB_Initialize(cambPath,cambVersion)
       ! Determine maximum wavenumber.
       wavenumberCAMB=exp(max(log(wavenumberRequired)+1.0d0,cambLogWavenumberMaximumDefault))
       if (wavenumberCAMB > wavenumberMaximum) then
          wavenumberCAMB=wavenumberMaximum
          wavenumberMaximumReached=.true.
       end if
       write (wavenumberLabel,'(e12.6)') wavenumberCAMB
       ! Construct input file for CAMB.
       call Get_Environment_Variable('HOSTNAME',hostName)
       workPath     =galacticusPath(pathTypeDataDynamic)//'largeScaleStructure/'
       parameterFile=workPath//'transfer_function_parameters'//'_'//trim(hostName)//'_'//GetPID()
       !$ parameterFile=parameterFile//'_'//OMP_Get_Thread_Num()
       parameterFile=parameterFile//'.txt'
       call System_Command_Do("mkdir -p "//char(workPath))
       open(newunit=cambParameterFile,file=char(parameterFile),status='unknown',form='formatted')
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'output_root                  ','camb'
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'get_scalar_cls               ','F'
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'get_vector_cls               ','F'
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'get_tensor_cls               ','F'
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'get_transfer                 ','T'
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'do_lensing                   ','F'
       write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'do_nonlinear                 ',0
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'l_max_scalar                 ',2200.0d0
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'l_max_tensor                 ',1500.0d0
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'k_eta_max_tensor             ',3000.0d0
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'use_physical                 ','F'
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_baryon                 ',cosmologyParameters_%OmegaBaryon    ()
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_cdm                    ',cosmologyParameters_%OmegaMatter    ()-cosmologyParameters_%OmegaBaryon()
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_lambda                 ',cosmologyParameters_%OmegaDarkEnergy()
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'omega_neutrino               ',0.0d0
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'omk                          ',0.0d0
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'hubble                       ',cosmologyParameters_%HubbleConstant ()
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'w                            ',-1.0d0
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'cs2_lam                      ',1.0d0
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'temp_cmb                     ',cosmologyParameters_%temperatureCMB ()
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
       write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'transfer_k_per_logint        ',0
       write (cambParameterFile,'(a,1x,"=",1x,i1   )') 'transfer_num_redshifts       ',1
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'transfer_interp_matterpower  ','T'
       write (cambParameterFile,'(a,1x,"=",1x,e12.6)') 'transfer_redshift(1)         ',0.0d0
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'transfer_filename(1)         ','transfer_out.dat'
       write (cambParameterFile,'(a,1x,"=",1x,a    )') 'transfer_matterpower(1)      ','matterpower.dat'
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
       countWavenumber=Count_Lines_In_File("camb_transfer_out.dat","#")
       if (allocated(wavenumbers      )) deallocate(wavenumbers      )
       if (allocated(transferFunctions)) deallocate(transferFunctions)
       allocate(wavenumbers      (countWavenumber))
       allocate(transferFunctions(countWavenumber))
       open(newunit=cambTransferFile,file="camb_transfer_out.dat",status='old',form='formatted')
       i=0
       do while (i < countWavenumber)
          read (cambTransferFile,'(a)',iostat=status) cambTransferLine
          if (status == 0) then
             if (cambTransferLine(1:1) /= "#") then
                i=i+1
                read (cambTransferLine,*) wavenumbers(i),transferFunctions(i)
             end if
          else
             call Galacticus_Error_Report('unable to read CAMB transfer function file'//{introspection:location})
          end if
       end do
       close(cambTransferFile)
       ! Remove temporary files.
       command="rm -f "                    // &
            &   parameterFile         //" "// &
            &  "camb_params.ini"      //" "// &
            &  "camb_transfer_out.dat"//" "// &
            &  "camb_matterpower.dat"
       call System_Command_Do(command)
       ! Convert from CAMB units to Galacticus units.
       wavenumbers=+wavenumbers                                                   &
            &      *cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)       
       ! Construct the output HDF5 file.
       fileIsNew=.not.File_Exists(fileName)
       call cambOutput%openFile(char(fileName),objectsOverwritable=.true.)
       call cambOutput%writeDataset(wavenumbers      ,'wavenumber'      ,chunkSize=chunkSize,appendTo=fileIsNew)
       call cambOutput%writeDataset(transferFunctions,'transferFunction',chunkSize=chunkSize,appendTo=fileIsNew)
       call cambOutput%writeAttribute('Cold dark matter transfer function created by CAMB.','description')
       call cambOutput%writeAttribute(cambFormatVersionCurrent,'fileFormat')
       parametersGroup=cambOutput%openGroup('parameters')
       call parametersGroup%writeAttribute(cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
       call parametersGroup%writeAttribute(cosmologyParameters_%OmegaBaryon    (),'OmegaBaryon'    )
       call parametersGroup%writeAttribute(cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
       call parametersGroup%writeAttribute(cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
       call parametersGroup%writeAttribute(cosmologyParameters_%temperatureCMB (),'temperatureCMB' )
       call parametersGroup%close()
       extrapolationGroup          =cambOutput        %openGroup('extrapolation')
       extrapolationWavenumberGroup=extrapolationGroup%openGroup('wavenumber'   )
       call extrapolationWavenumberGroup%writeAttribute('extrapolate','low' )
       call extrapolationWavenumberGroup%writeAttribute('extrapolate','high')
       call extrapolationWavenumberGroup%close()
       call extrapolationGroup          %close()
       call cambOutput                  %close()
       ! Unlock the lock file.
       if (lockFileGlobally) then
          call File_Unlock(cambFileLockGlobal)
       else
          call File_Unlock(fileLock          )
       end if      
    end if
    return
  end subroutine Interface_CAMB_Transfer_Function
  
end module Interfaces_CAMB
