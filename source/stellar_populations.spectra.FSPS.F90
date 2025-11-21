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
Implements a stellar population spectra class which utilizes the FSPS package \citep{conroy_propagation_2009}.
!!}

  use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunction, initialMassFunctionClass

  !![
  <stellarPopulationSpectra name="stellarPopulationSpectraFSPS">
   <description>
    A stellar population spectra class utilizing the FSPS package \citep{conroy_propagation_2009}. If necessary, the
    \href{https://github.com/cconroy20/fsps}{\normalfont \ttfamily FSPS} code will be downloaded, patched and compiled and run
    to generate spectra. These tabulations are then stored to file for later re-use. The file name used is {\normalfont
    \ttfamily datasets/dynamic/stellarPopulations/simpleStellarPopulationsFSPS:v$&lt;$version$&gt;$\_$&lt;$descriptor$&gt;$.hdf5} where
    $&lt;${\normalfont \ttfamily version}$&gt;$ is the FSPS version used, and $&lt;${\normalfont \ttfamily descriptor}$&gt;$ is
    an MD5 hash descriptor of the selected stellar population.
   </description>
  </stellarPopulationSpectra>
  !!]
  type, extends(stellarPopulationSpectraFile) :: stellarPopulationSpectraFSPS
     !!{
     A stellar population spectra class which utilizes the FSPS package \citep{conroy_propagation_2009}.
     !!}
     private
     class(initialMassFunctionClass), pointer :: initialMassFunction_ => null()
   contains
     final     ::             fspsDestructor
     procedure :: readFile => fspsReadFile
  end type stellarPopulationSpectraFSPS

  interface stellarPopulationSpectraFSPS
     !!{
     Constructors for the FSPS stellar spectra class.
     !!}
     module procedure fspsConstructorParameters
     module procedure fspsConstructorInternal
  end interface stellarPopulationSpectraFSPS

contains

  function fspsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the FSPS stellar spectra class which takes a parameter set as input.
    !!}
    implicit none
    type   (stellarPopulationSpectraFSPS)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    class  (initialMassFunctionClass    ), pointer       :: initialMassFunction_
    logical                                              :: forceZeroMetallicity

    !![
    <inputParameter>
      <name>forceZeroMetallicity</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>Force the use of zero metallicity (or lowest metallicity available) for all stellar populations.</description>
    </inputParameter>
    <objectBuilder class="initialMassFunction" name="initialMassFunction_" source="parameters"/>
    !!]
    self=stellarPopulationSpectraFSPS(forceZeroMetallicity,initialMassFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="initialMassFunction_"/>
    !!]
    return
  end function fspsConstructorParameters

  function fspsConstructorInternal(forceZeroMetallicity,initialMassFunction_) result(self)
    !!{
    Internal constructor for the FSPS stellar spectra class.
    !!}
    use :: Input_Paths       , only : inputPath             , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string        , operator(//)
    use :: Interfaces_FSPS   , only : Interface_FSPS_Version
    implicit none
    type   (stellarPopulationSpectraFSPS)                        :: self
    logical                              , intent(in   )         :: forceZeroMetallicity
    class  (initialMassFunctionClass    ), intent(in   ), target :: initialMassFunction_
    type   (varying_string              )                        :: fspsVersion
    !![
    <constructorAssign variables="forceZeroMetallicity, *initialMassFunction_"/>
    !!]

    call Interface_FSPS_Version(fspsVersion)
    self%stellarPopulationSpectraFile=stellarPopulationSpectraFile(forceZeroMetallicity,char(inputPath(pathTypeDataDynamic)//'stellarPopulations/simpleStellarPopulationsFSPS:v'//fspsVersion//'_'//self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)//'.hdf5'))
    return
  end function fspsConstructorInternal

  subroutine fspsDestructor(self)
    !!{
    Destructor for the \refClass{stellarPopulationSpectraFSPS} stellar population class.
    !!}
    implicit none
    type(stellarPopulationSpectraFSPS), intent(inout) :: self

    !![
    <objectDestructor name="self%initialMassFunction_"/>
    !!]
    return
  end subroutine fspsDestructor

  subroutine fspsReadFile(self)
    !!{
    Ensure that the requested stellar population has been generated.
    !!}
    use :: File_Utilities , only : File_Exists                 , File_Lock     , File_Unlock, lockDescriptor, &
         &                         File_Path                   , Directory_Make
    use :: HDF5_Access    , only : hdf5Access
    use :: IO_HDF5        , only : hdf5Object
    use :: Interfaces_FSPS, only : Interface_FSPS_SSPs_Tabulate
    use :: Tables         , only : table1D
    implicit none
    class  (stellarPopulationSpectraFSPS), intent(inout) :: self
    class  (table1D                     ), allocatable   :: imf
    logical                                              :: remakeFile
    integer                                              :: fileFormatVersion
    type   (hdf5Object                  )                :: spectraFile
    type   (lockDescriptor              )                :: fileLock
    integer                                              :: i

    ! Decide if we need to read the file.
    if (.not.self%fileRead) then
       ! Check if the file exists and has the correct version.
       remakeFile=.false.
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call Directory_Make(File_Path(self%fileName))
       do i=1,2
          call File_Lock(self%fileName,fileLock,lockIsShared=i == 1)
          if (File_Exists(self%fileName)) then
             !$ call hdf5Access%set()
             call spectraFile%openFile     (self%fileName,readOnly         =.true.)
             call spectraFile%readAttribute('fileFormat' ,fileFormatVersion       )
             if (fileFormatVersion /= fileFormatVersionCurrent) remakeFile=.true.
             !$ call hdf5Access%unset()
          else
             remakeFile=.true.
          end if
          if (remakeFile .and. i == 1) then
             remakeFile=.false.
             call File_Unlock(fileLock,sync=.false.)
             cycle
          end if
          ! Generate the file if necessary.
          if (remakeFile) then
             ! Generate the IMF tabulation.
             call self%initialMassFunction_%tabulate(imf)
             ! Build the SSPs.
             call Interface_FSPS_SSPs_Tabulate(imf,self%initialMassFunction_%label(),fileFormatVersionCurrent,self%fileName)
          end if
          ! Call the parent file reader to complete reading.
          call self%stellarPopulationSpectraFile%readFile()
          call File_Unlock(fileLock)
       end do
    end if
    return
  end subroutine fspsReadFile
