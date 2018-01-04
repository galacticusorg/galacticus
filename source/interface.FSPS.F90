!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which provides various interfaces to the FSPS code \citep{conroy_propagation_2009}.

module Interfaces_FSPS
  !% Provides various interfaces to the FSPS code \citep{conroy_propagation_2009}.
  private
  public :: Interface_FSPS_Initialize, Interface_FSPS_SSPs_Tabulate
  
contains

  subroutine Interface_FSPS_Initialize(fspsPath,fspsVersion)
    !% Initialize the interface with FSPS, including downloading and compiling FSPS if necessary.
    use ISO_Varying_String
    use Galacticus_Input_Paths
    use File_Utilities
    use System_Command
    use Galacticus_Display
    use Galacticus_Error
    use String_Handling
    implicit none
    type     (varying_string), intent(  out) :: fspsPath       , fspsVersion
    integer                                  :: status         , inputFile
    logical                                  :: upToDate
    character(len=40        )                :: currentRevision
    
    ! Specify source code path.
    fspsPath=Galacticus_Input_Path()//"/aux/FSPS_v2.5"
    ! Check out the code.
    if (.not.File_Exists(fspsPath)) then
       call Galacticus_Display_Message("downloading FSPS source code....",verbosityWorking)
       call System_Command_Do("git clone git://github.com/cconroy20/fsps.git/ "//fspsPath,status)
       if (.not.File_Exists(fspsPath) .or. status /= 0) call Galacticus_Error_Report("failed to clone FSPS git repository"//{introspection:location})
    end if
    ! Get the code revision number.
    call System_Command_Do("cd "//fspsPath//"; git rev-parse HEAD > currentRevision.txt",status)
    if (status /= 0) call Galacticus_Error_Report("unable to find FSPS revision"//{introspection:location})
    open(newUnit=inputFile,file=char(fspsPath)//"/currentRevision.txt",status='old',form='formatted')
    read (inputFile,'(a)') currentRevision
    close(inputFile)
    fspsVersion="v2.5; "//currentRevision
    ! Check for updates to the code.
    call System_Command_Do("cd "//fspsPath//"; git status | grep -q ""Your branch is up-to-date with 'origin/master'""",status)
    upToDate=(status == 0)
    if (.not.upToDate) then
       call Galacticus_Display_Message("updating FSPS source code",verbosityWorking)
       ! Update and remove the galacticus_IMF.f90 file to trigger re-patching of the code.
       call System_Command_Do("cd "//fspsPath//"; git checkout -- .; git pull; rm -f src/galacticus_IMF.f90")
    end if
    ! Patch the code.
    if (.not.File_Exists(fspsPath//"/src/galacticus_IMF.f90")) then
       call System_Command_Do("cp "//Galacticus_Input_Path()//"/aux/FSPS_v2.5_Galacticus_Modifications/galacticus_IMF.f90 "//fspsPath//"/src/"                                                   ,status)
       if (status /= 0) call Galacticus_Error_Report("failed to copy FSPS patch 'galacticus_IMF.f90'"//{introspection:location})
       call System_Command_Do("cp "//Galacticus_Input_Path()//"/aux/FSPS_v2.5_Galacticus_Modifications/imf.f90.patch "     //fspsPath//"/src/; cd "//fspsPath//"/src/; patch < imf.f90.patch"    ,status)
       if (status /= 0) call Galacticus_Error_Report("failed to patch FSPS file 'imf.f90'"           //{introspection:location})
       call System_Command_Do("cp "//Galacticus_Input_Path()//"/aux/FSPS_v2.5_Galacticus_Modifications/ssp_gen.f90.patch " //fspsPath//"/src/; cd "//fspsPath//"/src/; patch < ssp_gen.f90.patch",status)
       if (status /= 0) call Galacticus_Error_Report("failed to patch FSPS file 'ssp_gen.f90'"       //{introspection:location})
       call System_Command_Do("cp "//Galacticus_Input_Path()//"/aux/FSPS_v2.5_Galacticus_Modifications/autosps.f90.patch " //fspsPath//"/src/; cd "//fspsPath//"/src/; patch < autosps.f90.patch",status)
       if (status /= 0) call Galacticus_Error_Report("failed to patch FSPS file 'autosps.f90'"       //{introspection:location})
       call System_Command_Do("cp "//Galacticus_Input_Path()//"/aux/FSPS_v2.5_Galacticus_Modifications/Makefile.patch "    //fspsPath//"/src/; cd "//fspsPath//"/src/; patch < Makefile.patch"   ,status)
       if (status /= 0) call Galacticus_Error_Report("failed to patch FSPS file 'Makefile'"          //{introspection:location})
       call System_Command_Do("rm -f "//fspsPath//"/src/autosps.exe")
    end if
    !  Build the code.
    if (.not.File_Exists(fspsPath//"/src/autosps.exe")) then
       call Galacticus_Display_Message("compiling autosps.exe code",verbosityWorking)
       call System_Command_Do("cd "//fspsPath//"/src; export SPS_HOME="//fspsPath//"; make clean; make -j 1",status)
       if (.not.File_Exists(fspsPath//"/src/autosps.exe") .or. status /= 0) call Galacticus_Error_Report("failed to build autosps.exe code"//{introspection:location})
    end if
    return
  end subroutine Interface_FSPS_Initialize

  subroutine Interface_FSPS_SSPs_Tabulate(imf,imfName,fileFormat,spectraFileName)
    !% Tabulate simple stellar populations for the given \gls{imf} using FSPS.
    use ISO_Varying_String
    use Tables
    use System_Command
    use File_Utilities
    use String_Handling
    use IO_HDF5
    use Dates_and_Times
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    implicit none
    class           (table1D       ), intent(inout)                              :: imf
    type            (varying_string), intent(in   )                              :: imfName             , spectraFileName
    integer                         , intent(in   )                              :: fileFormat
    integer                         , parameter                                  :: metallicityCount =22
    double precision                               , dimension(metallicityCount) :: metallicity
    double precision                , allocatable  , dimension(:               ) :: wavelength          , age
    double precision                , allocatable  , dimension(:,:,:           ) :: spectrum
    integer                         , parameter                                  :: fileFormatCurrent= 1
    type            (varying_string)                                             :: fspsVersion         , fspsPath      , &
         &                                                                          outputFileName
    integer                                                                      :: iIMF                , outputFile    , &
         &                                                                          iMetallicity        , inputFile     , &
         &                                                                          iAge                , ageCount      , &
         &                                                                          wavelengthCount
    character       (len=256       )                                             :: line
    type            (hdf5Object    )                                             :: spectraFile         , imfGroup      , &
         &                                                                          dataset

    ! Validate file format.
    if (fileFormat /= fileFormatCurrent) call Galacticus_Error_Report(var_str("FSPS interface supports file format version ")//fileFormatCurrent//" but version "//fileFormat//" was requested"//{introspection:location})
    ! Ensure FSPS is available.
    call Interface_FSPS_Initialize(fspsPath,fspsVersion)
    ! Output the IMF file.
    open(newUnit=outputFile,file="galacticus.imf",status='unknown',form='formatted')
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
          ! Create parameter file for FSPS.
          open(newUnit=outputFile,file="fsps.inp",status='unknown',form='formatted')
          write (outputFile,'(i1)') 6                    ! IMF.
          write (outputFile,'(i1)') 0                    ! Generate SSP.
          write (outputFile,'(i2)') iMetallicity         ! Specify metallicity.
          write (outputFile,'(a)' ) "no"                 ! Do not include dust.
          write (outputFile,'(a)' ) char(outputFileName) ! Specify filename.
          close(outputFile)       
          call System_Command_Do("export SPS_HOME="//fspsPath//"; "//fspsPath//"/src/autosps.exe < fsps.inp")
          call System_Command_Do("rm -f fsps.inp")
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
    call System_Command_Do("rm -f galacticus.imf")
    ! Convert ages from loagrithmic form.
    age=10.0d0**(age-9.0d0)
    ! Write output file.
    !$omp critical (HDF5_Access)
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
    call dataset    %writeAttribute(1.0d0/angstromsPerMeter             ,'unitsInSI'                        )
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
    !$omp end critical (HDF5_Access)
  return
  end subroutine Interface_FSPS_SSPs_Tabulate
  
end module Interfaces_FSPS
