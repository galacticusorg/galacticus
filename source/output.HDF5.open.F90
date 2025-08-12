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
Contains a module which handles opening of the \glc\ output file.
!!}

module Output_HDF5_Open
  !!{
  Handles opening of the \glc\ output file.
  !!}
  use :: ISO_Varying_String, only : varying_string
  use :: Error             , only : errorStatusSuccess
  implicit none
  private
  public :: Output_HDF5_Open_File, Output_HDF5_Close_File, Output_HDF5_Completion_Status, Output_HDF5_Set_Group

  ! Output file name.
  type   (varying_string) :: outputFileName                     , outputScratchFileName

  ! Completion status.
  integer                 :: statusCompletion=errorStatusSuccess
 
contains

  subroutine Output_HDF5_Open_File(parameters)
    !!{
    Open the file for \glc\ output.
    !!}
    use :: Output_HDF5       , only : hdf5SieveBufferSize , hdf5UseLatestFormat, hdf5CompressionLevel  , hdf5CacheElementsCount, &
         &                            outputFileIsOpen    , outputFile         , hdf5CacheSizeBytes    , hdf5ChunkSize         , &
         &                            outputGroup
    use :: HDF5              , only : hsize_t             , size_t
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : IO_HDF5_Set_Defaults, hdf5Object          , ioHDF5AccessInitialize
    use :: ISO_Varying_String, only : var_str             , char                , operator(//)          , extract              , &
         &                            len                 , operator(==)        , adjustl               , trim
    use :: Input_Parameters  , only : inputParameters     , inputParameter
#ifdef USEMPI
    use :: MPI_Utilities     , only : mpiSelf
#endif
    use :: String_Handling   , only : operator(//)
    implicit none
    type   (inputParameters), intent(inout) :: parameters
    integer(hsize_t        )                :: chunkSize
    integer                                 :: sieveBufferSize
    integer(size_t         )                :: cacheElementsCount, cacheSizeBytes
    type   (varying_string )                :: outputFileName_   , outputScratchFileName_
#ifdef USEMPI
    type   (varying_string )                :: fileNamePrefix
#endif
    
    if (.not.outputFileIsOpen) then
       !![
       <inputParameter>
         <name>outputFileName</name>
         <variable>outputFileName_</variable>
         <defaultValue>var_str('galacticus.hdf5')</defaultValue>
         <description>The name of the file to which \glc\ results will be written.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>outputScratchFileName</name>
         <variable>outputScratchFileName_</variable>
         <defaultValue>outputFileName_</defaultValue>
         <description>The name of the file to which \glc\ results will be written temporarily during runs.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>hdf5SieveBufferSize</name>
         <defaultValue>65536</defaultValue>
         <description>The size of the sieve buffer used by the HDF5 library to speed reading/writing of partial datasets.</description>
         <source>parameters</source>
         <variable>sieveBufferSize</variable>
       </inputParameter>
       !!]
       hdf5SieveBufferSize=sieveBufferSize
       !![
       <inputParameter>
         <name>hdf5UseLatestFormat</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether to use the latest HDF5 file format.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>hdf5CacheElementsCount</name>
         <defaultValue>521_size_t</defaultValue>
         <description>Number of elements in the raw data chunk cache.</description>
         <source>parameters</source>
         <variable>cacheElementsCount</variable>
       </inputParameter>
       !!]
       hdf5CacheElementsCount=cacheElementsCount
       !![
       <inputParameter>
         <name>hdf5CacheSizeBytes</name>
         <defaultValue>1048576_size_t</defaultValue>
         <description>Size of the raw data chunk cache in bytes.</description>
         <source>parameters</source>
         <variable>cacheSizeBytes</variable>
       </inputParameter>
       !!]
       hdf5CacheSizeBytes=cacheSizeBytes
       ! Remove leading and trailing spaces.
       outputFileName       =trim(adjustl(outputFileName_       ))
       outputScratchFileName=trim(adjustl(outputScratchFileName_))
       ! Modify the file name on a per-process basis if running under MPI.
#ifdef USEMPI
       if (extract(outputFileName       ,len(outputFileName       )-4,len(outputFileName       )) == ".hdf5") then
          fileNamePrefix                 =extract(outputFileName       ,1,len(outputFileName       )-5)
          outputFileName       =fileNamePrefix//':MPI'//mpiSelf%rankLabel()//'.hdf5'
       else
          fileNamePrefix                 =outputFileName
          outputFileName       =fileNamePrefix//':MPI'//mpiSelf%rankLabel()
       end if
       if (extract(outputScratchFileName,len(outputScratchFileName)-4,len(outputScratchFileName)) == ".hdf5") then
          fileNamePrefix                 =extract(outputScratchFileName,1,len(outputScratchFileName)-5)
          outputScratchFileName=fileNamePrefix//':MPI'//mpiSelf%rankLabel()//'.hdf5'
       else
          fileNamePrefix                 =outputScratchFileName
          outputScratchFileName=fileNamePrefix//':MPI'//mpiSelf%rankLabel()
       end if
#endif
       ! Open the file.
       call ioHDF5AccessInitialize()
       !$ call hdf5Access%set()
       outputFile=hdf5Object(                                                 &
            &                                    char(outputScratchFileName), &
            &                overWrite          =.true.                     , &
            &                objectsOverwritable=.true.                     , &
            &                sieveBufferSize    =hdf5SieveBufferSize        , &
            &                useLatestFormat    =hdf5UseLatestFormat        , &
            &                cacheElementsCount =hdf5CacheElementsCount     , &
            &                cacheSizeBytes     =hdf5CacheSizeBytes           &
            &               )
       call outputFile%deepCopy(outputGroup)
       !$ call hdf5Access%unset()
       ! Now that the parameter file is open, we can open an output group in it for parameters.
       call parameters%parametersGroupOpen(outputFile)
       ! Read parameters.
       !![
       <inputParameter>
         <name>hdf5ChunkSize</name>
         <defaultValue>1024_hsize_t</defaultValue>
         <description>The chunk size used for outputting HDF5 datasets.</description>
         <source>parameters</source>
         <variable>chunksize</variable>
       </inputParameter>
       !!]
       hdf5ChunkSize=chunksize
       !![
       <inputParameter>
         <name>hdf5CompressionLevel</name>
         <defaultValue>-1</defaultValue>
         <description>The compression level used for outputting HDF5 datasets.</description>
         <source>parameters</source>
       </inputParameter>
       !!]

       ! Set default chunking and compression levels.
       call IO_HDF5_Set_Defaults(hdf5ChunkSize,hdf5CompressionLevel)

       ! Call all functions that requested to output to the file on start up.
       !![
       <eventHookStatic name="outputFileOpen"/>
       !!]
       
       ! Flag that the file is now open.
       outputFileIsOpen=.true.
    end if
    return
  end subroutine Output_HDF5_Open_File

  subroutine Output_HDF5_Close_File
    !!{
    Close the \glc\ output file.
    !!}
    use :: Output_HDF5       , only : outputFileIsOpen, outputFile, outputGroup
    use :: File_Utilities    , only : File_Rename
    use :: HDF5_Access       , only : hdf5Access
    use :: ISO_Varying_String, only : operator(/=)
    implicit none

    ! Perform any final tasks prior to shutdown.
    if (outputFileIsOpen) then
       !$omp critical (Output_HDF5_Close_File)
       if (outputFileIsOpen) then
          !![
	  <eventHookStatic name="outputFileClose"/>
          <eventHook       name="outputFileClose"/>
          !!]
          ! Close the file.
          !$ call hdf5Access%set()
          call outputFile%writeAttribute(statusCompletion,"statusCompletion")
          !$ call hdf5Access%unset()
          ! Move the scratch file to the final file if necessary.
          if (outputFileName /= outputScratchFileName) call File_Rename(outputScratchFileName,outputFileName,overwrite=.true.)
          ! Record that the file is now closed.
          outputFileIsOpen=.false.
       end if
       !$omp end critical (Output_HDF5_Close_File)
    end if
    return
  end subroutine Output_HDF5_Close_File

  subroutine Output_HDF5_Completion_Status(status)
    !!{
    Set the completion status.
    !!}
    implicit none
    integer, intent(in   ) :: status

    statusCompletion=status
    return
  end subroutine Output_HDF5_Completion_Status
  
  subroutine Output_HDF5_Set_Group(nameGroup)
    !!{
    Set the name of the current output group.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: ISO_Varying_String, only : varying_string  , char
    use :: Output_HDF5       , only : outputFileIsOpen, outputFile, outputGroup
    implicit none
    type(varying_string), intent(in   ) :: nameGroup

    if (.not.outputFileIsOpen) call Error_Report('can not set the output group - file is not open'//{introspection:location})
    !$ call hdf5Access%set()
    outputGroup=outputFile%openGroup(char(nameGroup))
    !$ call hdf5Access%unset()
    return
  end subroutine Output_HDF5_Set_Group
  
end module Output_HDF5_Open
