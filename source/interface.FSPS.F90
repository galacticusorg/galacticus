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
Contains a module which provides various interfaces to the FSPS code \citep{conroy_propagation_2009}.
!!}

module Interfaces_FSPS
  !!{
  Provides various interfaces to the FSPS code \citep{conroy_propagation_2009}.
  !!}
  use :: File_Utilities, only : lockDescriptor
  private
  public :: Interface_FSPS_Initialize, Interface_FSPS_SSPs_Tabulate, Interface_FSPS_Version

  ! Lock object to prevent multiple threads/processes attempting to build the code simultaneously.
  type(lockDescriptor) :: fspsLock

  ! Lock object to prevent multiple threads/processes attempting to build the same IMF simultaneously.
  type(lockDescriptor) :: imfLock

contains

  subroutine Interface_FSPS_Version(fspsVersion)
    !!{
    Set the version of FSPS being used.
    !!}
    use :: Dependencies      , only : dependencyVersion
    use :: ISO_Varying_String, only : varying_string   , assignment(=)
    implicit none
    type(varying_string), intent(  out) :: fspsVersion
    
    fspsVersion=dependencyVersion("fsps")
    return
  end subroutine Interface_FSPS_Version
  
  subroutine Interface_FSPS_Initialize(fspsPath,fspsVersion,static)
    !!{
    Initialize the interface with FSPS, including downloading and compiling FSPS if necessary.
    !!}
    use :: Display           , only : displayMessage   , verbosityLevelWorking
    use :: File_Utilities    , only : File_Exists      , File_Lock            , File_Remove       , File_Unlock
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath        , pathTypeDataDynamic  , pathTypeDataStatic
    use :: ISO_Varying_String, only : assignment(=)    , char                 , operator(//)      , varying_string
    use :: String_Handling   , only : operator(//)
    use :: System_Command    , only : System_Command_Do
    use :: System_Download   , only : download
    use :: System_Compilers  , only : compiler         , compilerOptions      , languageFortran
    implicit none
    type     (varying_string), intent(  out)           :: fspsPath, fspsVersion
    logical                  , intent(in   ), optional :: static
    integer                                            :: status
    type     (varying_string)                          :: lockPath
    !![
    <optionalArgument name="static" defaultsTo=".false." />
    !!]
#include "os.inc"

    ! Specify source code path.
    call Interface_FSPS_Version(fspsVersion)
    fspsPath=inputPath(pathTypeDataDynamic)//"fsps-"//fspsVersion
    lockPath=inputPath(pathTypeDataDynamic)//"fsps" //fspsVersion
    call File_Lock(char(lockPath),fspsLock)
    !  Build the code if the executable does not exist.
    if (.not.File_Exists(fspsPath//"/src/autosps.exe")) then
       ! Download the code if not already done.
       if (.not.File_Exists(fspsPath)) then
          if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"FSPS_"//char(fspsVersion)//".tar.gz")) then
             call displayMessage("downloading FSPS source code....",verbosityLevelWorking)
             call download("https://github.com/cconroy20/fsps/archive/refs/tags/v"//char(fspsVersion)//".tar.gz",char(inputPath(pathTypeDataDynamic))//"FSPS_"//char(fspsVersion)//".tar.gz",status=status,retries=5,retryWait=60)
             if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"FSPS_"//char(fspsVersion)//".tar.gz") .or. status /= 0) call Error_Report("failed to download FSPS"//{introspection:location})
          end if
          call displayMessage("unpacking FSPS code....",verbosityLevelWorking)
          call System_Command_Do("tar -x -v -z -C "//inputPath(pathTypeDataDynamic)//" -f "//inputPath(pathTypeDataDynamic)//"FSPS_"//char(fspsVersion)//".tar.gz",status)          
          if (status /= 0 .or. .not.File_Exists(fspsPath)) call Error_Report('failed to unpack FSPS code'//{introspection:location})
       end if
       ! Patch the code if not already patched.
       if (.not.File_Exists(fspsPath//"/src/galacticus_IMF.f90")) then
          call System_Command_Do("cp "//inputPath(pathTypeDataStatic)//"patches/FSPS/galacticus_IMF.f90 "//fspsPath//"/src/"                                                    ,status)
          if (status /= 0) call Error_Report("failed to copy FSPS patch 'galacticus_IMF.f90'"//{introspection:location})
          call System_Command_Do("cp "//inputPath(pathTypeDataStatic)//"patches/FSPS/imf.f90.patch "     //fspsPath//"/src/; cd "//fspsPath//"/src/; patch < imf.f90.patch"     ,status)
          if (status /= 0) call Error_Report("failed to patch FSPS file 'imf.f90'"           //{introspection:location})
          call System_Command_Do("cp "//inputPath(pathTypeDataStatic)//"patches/FSPS/ssp_gen.f90.patch " //fspsPath//"/src/; cd "//fspsPath//"/src/; patch < ssp_gen.f90.patch" ,status)
          if (status /= 0) call Error_Report("failed to patch FSPS file 'ssp_gen.f90'"       //{introspection:location})
          call System_Command_Do("cp "//inputPath(pathTypeDataStatic)//"patches/FSPS/sps_vars.f90.patch "//fspsPath//"/src/; cd "//fspsPath//"/src/; patch < sps_vars.f90.patch",status)
          if (status /= 0) call Error_Report("failed to patch FSPS file 'sps_vars.f90'"      //{introspection:location})
          call System_Command_Do("cp "//inputPath(pathTypeDataStatic)//"patches/FSPS/autosps.f90.patch " //fspsPath//"/src/; cd "//fspsPath//"/src/; patch < autosps.f90.patch" ,status)
          if (status /= 0) call Error_Report("failed to patch FSPS file 'autosps.f90'"       //{introspection:location})
          call System_Command_Do("cp "//inputPath(pathTypeDataStatic)//"patches/FSPS/Makefile.patch "    //fspsPath//"/src/; cd "//fspsPath//"/src/; patch < Makefile.patch"    ,status)
          if (status /= 0) call Error_Report("failed to patch FSPS file 'Makefile'"          //{introspection:location})
          call File_Remove(fspsPath//"/src/autosps.exe")
       end if
       call displayMessage("compiling autosps.exe code",verbosityLevelWorking)
       if (static_) then
          call System_Command_Do(                                                                                &
               &                 "cd "//fspsPath//"/src; "                                                    // &
#ifndef __APPLE__
               &                 "grep -P '^F90FLAGS := ' Makefile && "                                       // &
#endif
               &                 "sed -i~ -E s/'^(F90FLAGS := [^#]*)'/'\1 \-static'/g Makefile"               ,  &
               &                 status                                                                          &
               &                )
       else
          call System_Command_Do(                                                                                &
               &                 "cd "//fspsPath//"/src; "                                                    // &
#ifndef __APPLE__
               &                 "grep -P '^F90FLAGS := ' Makefile && "                                       // &
#endif
               &                 "sed -i~ -E s/'^(F90FLAGS := .*)[[:space:]]*\-static(.*)'/'\1 \2'/g Makefile",  &
               &                 status                                                                          &
               &                )
       end if
       if (status /= 0) call Error_Report("failed to patch FSPS file 'Makefile' for static/dynamic build"//{introspection:location})
       call System_Command_Do(                                                                                                                                                &
            &                 "cd "//fspsPath//"/src; export SPS_HOME="//fspsPath//'; export F90FLAGS="'                                                                   // &
#ifndef __aarch64__
            &                 '-mcmodel=medium '                                                                                                                           // & ! Larger memory model required except on Arm64.
#endif
            &                 char(compilerOptions(languageFortran))//'"; sed -i~ -E s/"gfortran"/"'//char(compiler(languageFortran))//'"/ Makefile; make clean; make -j 1',  &
            &                 status                                                                                                                                          &
            &                )
       if (.not.File_Exists(fspsPath//"/src/autosps.exe") .or. status /= 0) call Error_Report("failed to build autosps.exe code"//{introspection:location})
    end if
    call File_Unlock(fspsLock)
    return
  end subroutine Interface_FSPS_Initialize

  subroutine Interface_FSPS_SSPs_Tabulate(imf,imfName,fileFormat,spectraFileName)
    !!{
    Tabulate simple stellar populations for the given \gls{imf} using FSPS.
    !!}
    use :: Dates_and_Times                 , only : Formatted_Date_and_Time
    use :: File_Utilities                  , only : Directory_Make         , File_Exists    , File_Name_Temporary, File_Path, &
          &                                         File_Remove            , File_Lock      , File_Unlock
    use :: Error                           , only : Error_Report
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: ISO_Varying_String              , only : char                   , operator(//)   , trim               , var_str  , &
          &                                         varying_string
    use :: Numerical_Constants_Astronomical, only : gigaYear               , luminositySolar, massSolar
    use :: Numerical_Constants_Units       , only : metersToAngstroms
    use :: String_Handling                 , only : operator(//)
    use :: System_Command                  , only : System_Command_Do
    use :: Tables                          , only : table1D
    implicit none
    class           (table1D       ), intent(inout)                              :: imf
    type            (varying_string), intent(in   )                              :: imfName             , spectraFileName
    integer                         , intent(in   )                              :: fileFormat
    integer                         , parameter                                  :: metallicityCount =22
    double precision                               , dimension(metallicityCount) :: metallicity
    double precision                , allocatable  , dimension(:               ) :: wavelength          , age
    double precision                , allocatable  , dimension(:,:,:           ) :: spectrum
    integer                         , parameter                                  :: fileFormatCurrent= 1
    type            (varying_string)                                             :: fspsVersion         , fspsPath         , &
         &                                                                          outputFileName      , fspsInputFileName, &
         &                                                                          imfFileName
    integer                                                                      :: iIMF                , outputFile       , &
         &                                                                          iMetallicity        , inputFile        , &
         &                                                                          iAge                , ageCount         , &
         &                                                                          wavelengthCount
    character       (len=256       )                                             :: line
    type            (hdf5Object    )                                             :: spectraFile         , imfGroup         , &
         &                                                                          dataset

    ! Validate file format.
    if (fileFormat /= fileFormatCurrent) call Error_Report(var_str("FSPS interface supports file format version ")//fileFormatCurrent//" but version "//fileFormat//" was requested"//{introspection:location})
    ! Ensure FSPS is available.
    call Interface_FSPS_Initialize(fspsPath,fspsVersion)
    ! Output the IMF file.
    imfFileName=File_Name_Temporary("fsps.imf")
    open(newUnit=outputFile,file=char(imfFileName),status='unknown',form='formatted')
    do iIMF=1,imf%size()
       write (outputFile,'(2(1x,e12.6))') imf%x(iIMF),imf%y(iIMF)
    end do
    close(outputFile)
    ! Iterate over metallicities.
    do iMetallicity=1,metallicityCount
       ! Construct output file name.
       outputFileName="imf"//trim(imfName)//".iZ"//iMetallicity
       ! Generate output file if necessary.
       if (.not.File_Exists(fspsPath//"/OUTPUTS/"//outputFileName//".spec")) then
          call File_Lock(char(fspsPath//"/OUTPUTS/"//outputFileName),imfLock)
          if (.not.File_Exists(fspsPath//"/OUTPUTS/"//outputFileName//".spec")) then
             ! Create parameter file for FSPS.
             fspsInputFileName=File_Name_Temporary("fsps.inp")
             open(newUnit=outputFile,file=char(fspsInputFileName),status='unknown',form='formatted')
             write (outputFile,'(i1)') 6                    ! IMF.
             write (outputFile,'(a)' ) char(   imfFileName) ! Specify IMF filename.
             write (outputFile,'(i1)') 0                    ! Generate SSP.
             write (outputFile,'(i2)') iMetallicity         ! Specify metallicity.
             write (outputFile,'(a)' ) "no"                 ! Do not include dust.
             write (outputFile,'(a)' ) char(outputFileName) ! Specify filename.
             close(outputFile)
             call System_Command_Do("export SPS_HOME="//fspsPath//"; "//fspsPath//"/src/autosps.exe < "//fspsInputFileName)
             call File_Remove(fspsInputFileName)
          end if
          call File_Unlock(imfLock)
       end if
       ! Parse the output file.
       open(newUnit=inputFile,file=char(fspsPath//"/OUTPUTS/"//outputFileName//".spec"),status='old',form='formatted')
       ! First line contains metallicity information.
       read (inputFile,'(a)') line
       read (line(index(line,":")+1:),*) metallicity(iMetallicity)
       do while (line(1:1) == "#")
          read (inputFile,'(a)') line
       end do
       read (line,*) ageCount,wavelengthCount
       if (.not.allocated(wavelength)) then
          allocate(age       (                ageCount                 ))
          allocate(wavelength(wavelengthCount                          ))
          allocate(spectrum  (wavelengthCount,ageCount,metallicityCount))
       end if
       read (inputFile,*) wavelength
       do iAge=1,ageCount
          read (inputFile,*) age     (  iAge             )
          read (inputFile,*) spectrum(:,iAge,iMetallicity)
       end do
       close(inputFile)
    end do
    ! Clean up.
    call File_Remove(imfFileName)
    ! Convert ages from logarithmic form.
    age=10.0d0**(age-9.0d0)
    ! Write output file.
    call Directory_Make(File_Path(spectraFileName))
    !$ call hdf5Access%set()
    call spectraFile%openFile(char(spectraFileName))
    ! Add metadata.
    call spectraFile%writeAttribute('Galacticus'                                                                           ,'createdBy'  )
    call spectraFile%writeAttribute(Formatted_Date_and_Time()                                                              ,'timestep'   )
    call spectraFile%writeAttribute(fileFormatCurrent                                                                      ,'fileFormat' )
    call spectraFile%writeAttribute(fspsVersion                                                                            ,'fspsVersion')
    call spectraFile%writeAttribute("Simple stellar population spectra from FSPS for a "//imfName//" initial mass function",'description')
    ! Add IMF.
    imfGroup=spectraFile%openGroup('initialMassFunction')
    call imfGroup%writeDataset  (imf%xs(),'mass'                  ,datasetReturned=dataset)
    call dataset %writeAttribute('M☉'                 ,'units'                           )
    call dataset %writeAttribute(      massSolar      ,'unitsInSI'                       )
    call dataset %close         (                                                         )
    call imfGroup%writeDataset  (imf%ys()   ,'initialMassFunction',datasetReturned=dataset)
    call dataset %writeAttribute('M☉⁻¹'               ,'units'                           )
    call dataset %writeAttribute(1.0d0/massSolar      ,'unitsInSI'                       )
    call dataset %close         (                                                        )
    call imfGroup%close         (                                                        )
    ! Write datasets.
    call spectraFile%writeDataset  (wavelength ,'wavelengths'        ,datasetReturned=dataset)
    call dataset    %writeAttribute('Å'                  ,'units'                            )
    call dataset    %writeAttribute(1.0d0/metersToAngstroms             ,'unitsInSI'         )
    call dataset    %close         (                                                         )
    call spectraFile%writeDataset  (age        ,'ages'         ,      datasetReturned=dataset)
    call dataset    %writeAttribute('Gyr'                ,'units'                            )
    call dataset    %writeAttribute(gigaYear             ,'unitsInSI'                        )
    call dataset    %close         (                                                         )
    call spectraFile%writeDataset  (metallicity,'metallicities'                              )
    call spectraFile%writeDataset  (spectrum   ,'spectra'            ,datasetReturned=dataset)
    call dataset    %writeAttribute('L☉ Hz⁻¹'            ,'units'                            )
    call dataset    %writeAttribute(luminositySolar      ,'unitsInSI'                        )
    call dataset    %close         (                                                         )
    call spectraFile%close()
    !$ call hdf5Access%unset()
  return
  end subroutine Interface_FSPS_SSPs_Tabulate

end module Interfaces_FSPS
