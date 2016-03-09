!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !# <stellarPopulationSpectra name="stellarPopulationSpectraFSPS">
  !#  <description>Provides stellar population spectra utilizing the FSPS package \citep{conroy_propagation_2009}. If necessary, the {\normalfont \ttfamily FSPS} code will be downloaded, patched and compiled and run to generate spectra. These tabulations are then stored to file for later re-use.</description>
  !# </stellarPopulationSpectra>
  type, extends(stellarPopulationSpectraFile) :: stellarPopulationSpectraFSPS
     !% A stellar population spectra class which utilizes the FSPS package \citep{conroy_propagation_2009}.
     private
   contains
     final     ::                  fspsDestructor
     procedure :: imfInitialize => fspsIMFInitialize
     procedure :: descriptor    => fspsDescriptor
  end type stellarPopulationSpectraFSPS

  interface stellarPopulationSpectraFSPS
     !% Constructors for the FSPS stellar spectra class.
     module procedure fspsConstructorParameters
     module procedure fspsConstructorInternal
  end interface stellarPopulationSpectraFSPS

contains
  
  function fspsConstructorParameters(parameters)
    !% Constructor for the FSPS stellar spectra class which takes a parameter set as input.
    use Star_Formation_IMF
    implicit none
    type   (stellarPopulationSpectraFSPS)                :: fspsConstructorParameters
    type   (inputParameters             ), intent(inout) :: parameters
    logical                                              :: forceZeroMetallicity
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>forceZeroMetallicity</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Force the use of zero metallicity (or lowest metallicity available) for all stellar populations.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    ! Simply use the internal constructor.
    fspsConstructorParameters=fspsConstructorInternal(forceZeroMetallicity)
    return
  end function fspsConstructorParameters
  
  function fspsConstructorInternal(forceZeroMetallicity)
    !% Internal constructor for the FSPS stellar spectra class.
    use Star_Formation_IMF
    use Galacticus_Input_Paths
    implicit none
    type   (stellarPopulationSpectraFSPS)                :: fspsConstructorInternal
    logical                              , intent(in   ) :: forceZeroMetallicity
    integer                                              :: imfCount               , i

    ! Store options.
    fspsConstructorInternal%forceZeroMetallicity=forceZeroMetallicity
    ! Set file names for all possible IMFs.
    imfCount=IMF_Available_Count()
    allocate(fspsConstructorInternal%fileName(imfCount))
    do i=1,imfCount
       fspsConstructorInternal%fileName(i)=Galacticus_Input_Path()//'data/stellarPopulations/SSP_Spectra_Conroy-et-al_v2.5_imf'//IMF_Descriptor(i)//'.hdf5'
    end do
    return
  end function fspsConstructorInternal
  
  subroutine fspsDestructor(self)
    !% Destructor for the file stellar spectra class.
    implicit none
    type(stellarPopulationSpectraFSPS), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine fspsDestructor

  subroutine fspsIMFInitialize(self,imfIndex)
    !% Ensure that the requested IMF has been generated.
    use Galacticus_Input_Paths
    use IO_HDF5
    use Star_Formation_IMF
    use File_Utilities
    use String_Handling
    use Tables
    use System_Command
    implicit none
    class  (stellarPopulationSpectraFSPS), intent(inout) :: self
    integer                              , intent(in   ) :: imfIndex
    class  (table1D                     ), allocatable   :: imf
    logical                                              :: readFile         , remakeFile
    integer                                              :: fileFormatVersion, iIMF      , &
         &                                                  imfUnit
    type   (varying_string              )                :: command
    type   (hdf5Object                  )                :: spectraFile

    ! Decide if we need to read the file.
    readFile=.not.allocated(self%imfLookup)
    if (.not.readFile) then
       if (size(self%imfLookup) < imfIndex) then
          readFile=.true.
       else
          readFile=(self%imfLookup(imfIndex)==0)
       end if
    end if
    ! If the file must be read, check if it needs to be generated.
    if (readFile) then
       ! Check if the file exists and has the correct version.
       remakeFile=.false.
       if (File_Exists(char(self%fileName(imfIndex)))) then
          !$omp critical(HDF5_Access)
          call spectraFile%openFile     (char(self%fileName(imfIndex)),readOnly         =.true.)
          call spectraFile%readAttribute('fileFormat'                 ,fileFormatVersion       )
          if (fileFormatVersion /= fileFormatVersionCurrent) remakeFile=.true.
          !$omp end critical(HDF5_Access)
       else
          remakeFile=.true.
       end if
       ! Generate the file if necessary.
       if (remakeFile) then
          ! Generate the IMF tabulation.
          call IMF_Tabulate(imfIndex,imf)
          open(newunit=imfUnit,file="galacticus.imf",status="unknown",form="formatted")
          do iIMF=1,imf%size()
             write (imfUnit,'(2(1x,e12.6))') imf%x(iIMF),imf%y(iIMF)
          end do
          close(imfUnit)
          call imf%destroy()
          ! Call the driver script to generate this file.
          command=char(Galacticus_Input_Path())//'scripts/aux/Conroy_SPS_Driver.pl '//IMF_Descriptor(imfIndex)//' '//self%fileName(imfIndex)
          command=command//' '//fileFormatVersionCurrent
          call System_Command_Do(command)
       end if
       ! Call the parent IMF initializor to complete initialization.
       call self%stellarPopulationSpectraFile%imfInitialize(imfIndex)
    end if
    return
  end subroutine fspsIMFInitialize

  subroutine fspsDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    implicit none
    class(stellarPopulationSpectraFSPS), intent(inout) :: self
    type (inputParameters             ), intent(inout) :: descriptor
    
    call descriptor%addParameter("stellarPopulationSpectraMethod","fsps")
    return
  end subroutine fspsDescriptor
