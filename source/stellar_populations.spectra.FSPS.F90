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

!% Implements a stellar population spectra class which utilizes the FSPS package \citep{conroy_propagation_2009}.

  use Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionClass, initialMassFunction
  
  !# <stellarPopulationSpectra name="stellarPopulationSpectraFSPS">
  !#  <description>Provides stellar population spectra utilizing the FSPS package \citep{conroy_propagation_2009}. If necessary, the {\normalfont \ttfamily FSPS} code will be downloaded, patched and compiled and run to generate spectra. These tabulations are then stored to file for later re-use.</description>
  !# </stellarPopulationSpectra>
  type, extends(stellarPopulationSpectraFile) :: stellarPopulationSpectraFSPS
     !% A stellar population spectra class which utilizes the FSPS package \citep{conroy_propagation_2009}.
     private
     class(initialMassFunctionClass), pointer :: initialMassFunction_ => null()
   contains
     final     ::             fspsDestructor
     procedure :: readFile => fspsReadFile
  end type stellarPopulationSpectraFSPS

  interface stellarPopulationSpectraFSPS
     !% Constructors for the FSPS stellar spectra class.
     module procedure fspsConstructorParameters
     module procedure fspsConstructorInternal
  end interface stellarPopulationSpectraFSPS

contains
  
  function fspsConstructorParameters(parameters) result(self)
    !% Constructor for the FSPS stellar spectra class which takes a parameter set as input.
    implicit none
    type   (stellarPopulationSpectraFSPS)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    class  (initialMassFunctionClass    ), pointer       :: initialMassFunction_
    logical                                              :: forceZeroMetallicity

    !# <inputParameter>
    !#   <name>forceZeroMetallicity</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Force the use of zero metallicity (or lowest metallicity available) for all stellar populations.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="initialMassFunction" name="initialMassFunction_" source="parameters"/>
    self=stellarPopulationSpectraFSPS(forceZeroMetallicity,initialMassFunction_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="initialMassFunction_"/>
    return
  end function fspsConstructorParameters
  
  function fspsConstructorInternal(forceZeroMetallicity,initialMassFunction_) result(self)
    !% Internal constructor for the FSPS stellar spectra class.
    use Galacticus_Paths
    implicit none
    type   (stellarPopulationSpectraFSPS)                        :: self
    logical                              , intent(in   )         :: forceZeroMetallicity
    class  (initialMassFunctionClass    ), intent(in   ), target :: initialMassFunction_
    !# <constructorAssign variables="forceZeroMetallicity, *initialMassFunction_"/>
    
    self%stellarPopulationSpectraFile=stellarPopulationSpectraFile(forceZeroMetallicity,char(galacticusPath(pathTypeDataDynamic)//'stellarPopulations/simpleStellarPopulationsFSPS:v2.5_'//self%hashedDescriptor(includeSourceDigest=.true.)//'.hdf5'))
    return
  end function fspsConstructorInternal
  
  subroutine fspsDestructor(self)
    !% Destructor for the {\normalfont \ttfamily FSPS} stellar population class.
    implicit none
    type(stellarPopulationSpectraFSPS), intent(inout) :: self

    !# <objectDestructor name="self%initialMassFunction_"/>
    return
  end subroutine fspsDestructor

  subroutine fspsReadFile(self)
    !% Ensure that the requested stellar population has been generated.
    use IO_HDF5
    use File_Utilities
    use Tables
    use Interfaces_FSPS
    implicit none
    class  (stellarPopulationSpectraFSPS), intent(inout) :: self
    class  (table1D                     ), allocatable   :: imf
    logical                                              :: remakeFile
    integer                                              :: fileFormatVersion
    type   (hdf5Object                  )                :: spectraFile

    ! Decide if we need to read the file.
    if (.not.self%fileRead) then
       ! Check if the file exists and has the correct version.
       remakeFile=.false.
       if (File_Exists(char(self%fileName))) then
          !$ call hdf5Access%set()
          call spectraFile%openFile     (char(self%fileName),readOnly         =.true.)
          call spectraFile%readAttribute('fileFormat'       ,fileFormatVersion       )
          if (fileFormatVersion /= fileFormatVersionCurrent) remakeFile=.true.
          !$ call hdf5Access%unset()
       else
          remakeFile=.true.
       end if
       ! Generate the file if necessary.
       if (remakeFile) then
          ! Generate the IMF tabulation.
          call self%initialMassFunction_%tabulate(imf)
          ! Build the SSPs.
          !$omp critical (stellarPopulationsFSPS)
          call Interface_FSPS_SSPs_Tabulate(imf,self%initialMassFunction_%label(),fileFormatVersionCurrent,self%fileName)
          !$omp end critical (stellarPopulationsFSPS)
       end if
       ! Call the parent file reader to complete reading.
       call self%stellarPopulationSpectraFile%readFile()
    end if
    return
  end subroutine fspsReadFile
