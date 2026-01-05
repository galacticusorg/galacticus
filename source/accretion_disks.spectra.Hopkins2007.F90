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

! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
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
  An implementation of the accretion disk spectra class using the model of \cite{hopkins_observational_2007}.
  !!}

  use :: File_Utilities, only : lockDescriptor

  !![
  <accretionDiskSpectra name="accretionDiskSpectraHopkins2007">
   <description>Accretion disk spectra using the model of \cite{hopkins_observational_2007}.</description>
  </accretionDiskSpectra>
  !!]
  type, extends(accretionDiskSpectraFile) :: accretionDiskSpectraHopkins2007
     !!{
     An accretion disk spectra class which uses the algorithm of \cite{hopkins_observational_2007}.
     !!}
     private
     type(lockDescriptor) :: fileLock
   contains
     !![
     <methods>
       <method description="Build the tabulation file containing AGN spectra." method="buildFile" />
     </methods>
     !!]
     procedure :: buildFile => hopkins2007BuildFile
  end type accretionDiskSpectraHopkins2007

  interface accretionDiskSpectraHopkins2007
     !!{
     Constructors for the \refClass{accretionDiskSpectraHopkins2007} accretion disk spectra class.
     !!}
     module procedure hopkins2007ConstructorParameters
     module procedure hopkins2007ConstructorInternal
  end interface accretionDiskSpectraHopkins2007

contains

  function hopkins2007ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{accretionDiskSpectraHopkins2007} accretion disk spectra class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (accretionDiskSpectraHopkins2007)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(blackHoleAccretionRateClass    ), pointer       :: blackHoleAccretionRate_
    class(accretionDisksClass            ), pointer       :: accretionDisks_

    !![
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    !!]
    self=accretionDiskSpectraHopkins2007(blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function hopkins2007ConstructorParameters

  function hopkins2007ConstructorInternal(blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Constructor for the \refClass{accretionDiskSpectraHopkins2007} accretion disk spectra class.
    !!}
    use :: File_Utilities, only : File_Lock, File_Unlock
    use :: Input_Paths   , only : inputPath, pathTypeDataStatic
    implicit none
    type (accretionDiskSpectraHopkins2007)                        :: self
    class(blackHoleAccretionRateClass    ), target, intent(in   ) :: blackHoleAccretionRate_
    class(accretionDisksClass            ), target, intent(in   ) :: accretionDisks_
    !![
    <constructorAssign variables="*blackHoleAccretionRate_, *accretionDisks_"/>
    !!]
    
    ! Set the file name.
    self%fileName=inputPath(pathTypeDataStatic)//"blackHoles/AGN_SEDs_Hopkins2007.hdf5"
    ! Build the file.
    call self%buildFile()
    ! Load the file.
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock    (char(self%fileName),self%fileLock,lockIsShared=.true.)
    call self%loadFile(char(self%fileName)                                  )
    call File_Unlock  (                    self%fileLock                    )
    return
  end function hopkins2007ConstructorInternal

  subroutine hopkins2007BuildFile(self)
    !!{
    Build a file containing a tabulation of the \cite{hopkins_observational_2007} model AGN spectra.
    !!}
    use            :: Dates_and_Times                 , only : Formatted_Date_and_Time
    use            :: Display                         , only : displayCounter         , displayCounterClear , displayIndent     , displayUnindent, &
          &                                                    verbosityLevelWorking
    use            :: File_Utilities                  , only : Count_Lines_in_File    , Directory_Make      , File_Exists       , File_Lock      , &
          &                                                    File_Unlock
    use            :: Error                           , only : Error_Report
    use            :: Input_Paths                     , only : inputPath              , pathTypeDataDynamic , pathTypeDataStatic
    use            :: HDF5_Access                     , only : hdf5Access
    use            :: IO_HDF5                         , only : hdf5Object
    use, intrinsic :: ISO_Fortran_Env
    use            :: Numerical_Constants_Astronomical, only : luminositySolar
    use            :: Numerical_Constants_Physical    , only : speedLight
    use            :: Numerical_Constants_Units       , only : metersToAngstroms
    use            :: Numerical_Ranges                , only : Make_Range             , rangeTypeLogarithmic
    use            :: String_Handling                 , only : operator(//)
    use            :: System_Command                  , only : System_Command_Do
    use            :: System_Download                 , only : download
    implicit none
    class           (accretionDiskSpectraHopkins2007), intent(inout)               :: self
    double precision                                 , dimension(:  ), allocatable :: wavelength                        , luminosityBolometric
    double precision                                 , dimension(:,:), allocatable :: SED
    double precision                                 , parameter                   :: luminosityBolometricMinimum=1.0d06
    double precision                                 , parameter                   :: luminosityBolometricMaximum=1.0d28
    integer                                          , parameter                   :: luminosityBolometricCount  =200
    logical                                                                        :: makeFile
    type            (hdf5Object                     )                              :: file                              , dataset
    integer                                                                        :: fileFormatCurrentFile             , sedUnit             , &
         &                                                                            i                                 , j                   , &
         &                                                                            status                            , wavelengthCount
    character       (len= 16                        )                              :: label
    character       (len=256                        )                              :: line
    double precision                                                               :: frequencyLogarithmic              , spectrumLogarithmic

    ! Determine if we need to make the file.
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(self%fileName),self%fileLock,lockIsShared=.true.)
    makeFile=.false.
    if (File_Exists(self%fileName)) then
       !$ call hdf5Access%set()
       file=hdf5Object(char(self%fileName),readOnly=.true.)
       if (file%hasAttribute('fileFormat')) then
          call file%readAttribute('fileFormat',fileFormatCurrentFile,allowPseudoScalar=.true.)
          makeFile=fileFormatCurrentFile /= fileFormatCurrent
       else
          makeFile=.true.
       end if
       !$ call hdf5Access%unset()
    else
       makeFile=.true.
    end if
    call File_Unlock(self%fileLock)
    ! Make the file if necessary.
    if (makeFile) then
       call displayIndent('Building file of tabulated AGN spectra for Hopkins2007 class',verbosity=verbosityLevelWorking)
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),self%fileLock,lockIsShared=.false.)
       ! Download the AGN SED code.
       if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"AGN_Spectrum/agn_spectrum.c")) then
          call Directory_Make(inputPath(pathTypeDataDynamic)//"/AGN_Spectrum")
          call download("http://www.tapir.caltech.edu/~phopkins/Site/qlf_files/agn_spectrum.c",char(inputPath(pathTypeDataStatic))//"aux/AGN_Spectrum/agn_spectrum.c",status=status)
          if (status /= 0 .or. .not.File_Exists(inputPath(pathTypeDataDynamic)//"AGN_Spectrum/agn_spectrum.c")) call Error_Report('failed to download agn_spectrum.c'//{introspection:location})
       end if
       ! Compile the AGN SED code.
       if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"AGN_Spectrum/agn_spectrum.x")) then
          call System_Command_Do("cd "//char(inputPath(pathTypeDataStatic))//"aux/AGN_Spectrum; gcc agn_spectrum.c -o agn_spectrum.x -lm");
          if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"AGN_Spectrum/agn_spectrum.x")) call Error_Report('failed to compile agn_spectrum.c'//{introspection:location})
       end if
       ! Generate a tabulation of AGN spectra over a sufficiently large range of AGN luminosity.
       call Directory_Make(inputPath(pathTypeDataStatic)//"blackHoles")
       allocate(luminosityBolometric(luminosityBolometricCount))
       luminosityBolometric=Make_Range(luminosityBolometricMinimum,luminosityBolometricMaximum,luminosityBolometricCount,rangeTypeLogarithmic)
       do i=1,luminosityBolometricCount
          call displayCounter(int(100.0*dble(i-1)/dble(luminosityBolometricCount)),isNew=i==1,verbosity=verbosityLevelWorking)
          write (label,'(e12.6)') log10(luminosityBolometric(i))
          call System_Command_Do(inputPath(pathTypeDataDynamic)//"AGN_Spectrum/agn_spectrum.x "//label//" > "//inputPath(pathTypeDataDynamic)//"AGN_Spectrum/SED.txt")
          wavelengthCount=Count_Lines_in_File(inputPath(pathTypeDataDynamic)//"AGN_Spectrum/SED.txt",";")-4
          if (allocated(wavelength)) then
             if (wavelengthCount /= size(wavelength)) call Error_Report('inconsistent number of wavelengths'//{introspection:location})
          else
             allocate(wavelength(wavelengthCount                          ))
             allocate(SED       (wavelengthCount,luminosityBolometricCount))
          end if
          open(newUnit=sedUnit,file=char(inputPath(pathTypeDataStatic))//"aux/AGN_Spectrum/SED.txt",status="old",form="formatted")
          j=wavelengthCount+1
          do
             read(sedUnit,'(a)',iostat=status) line
             if (status               == iostat_end) exit
             if (line(1:1)            == ";"       ) cycle
             read (line,*) frequencyLogarithmic,spectrumLogarithmic
             if (frequencyLogarithmic <  0.0d0     ) cycle
             j=j-1
             wavelength(j)=speedLight/10.0d0**frequencyLogarithmic*metersToAngstroms
             SED(j,i)=10.0d0**spectrumLogarithmic/10.0d0**frequencyLogarithmic
          end do
          close(sedUnit)
       end do
       call displayCounterClear(verbosity=verbosityLevelWorking)
       ! Store the data to file.
       !$ call hdf5Access%set()
       file=hdf5Object(char(self%fileName),overWrite=.true.)
       call file   %writeDataset  (wavelength               ,"wavelength"          ,datasetReturned=dataset)
       call dataset%writeAttribute("Angstroms (Å)"          ,"units"                                       )
       call dataset%writeAttribute(1.0d0/metersToAngstroms  ,"unitsInSI"                                   )
       call file   %writeDataset  (luminosityBolometric     ,"bolometricLuminosity",datasetReturned=dataset)
       call dataset%writeAttribute("Solar luminosities (L☉)","units"                                       )
       call dataset%writeAttribute(luminositySolar          ,"unitsInSI"                                   )
       call file   %writeDataset  (SED                      ,"SED"                 ,datasetReturned=dataset)
       call dataset%writeAttribute("L☉/Hz"                  ,"units"                                       )
       call dataset%writeAttribute(luminositySolar          ,"unitsInSI"                                   )
       ! Add some metadata.
       call file%writeAttribute("Computed using agn_spectrum.c downloaded from  http://www.tapir.caltech.edu/~phopkins/Site/qlf.html","source"      )
       call file%writeAttribute("http://adsabs.harvard.edu/abs/2007ApJ...654..731H"                                                  ,"URL"         )
       call file%writeAttribute("Hopkins et al. (2007)"                                                                              ,"reference"   )
       call file%writeAttribute(Formatted_Date_and_Time()                                                                            ,"creationTime")
       call file%writeAttribute(fileFormatCurrent                                                                                    ,"fileFormat"  )
       !$ call hdf5Access%unset()
       call File_Unlock(self%fileLock)
       call displayUnindent('done',verbosity=verbosityLevelWorking)
    end if
    return
  end subroutine hopkins2007BuildFile
