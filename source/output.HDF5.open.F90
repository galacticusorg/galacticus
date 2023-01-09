!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  implicit none
  private
  public :: Output_HDF5_Open_File, Output_HDF5_Close_File

  ! Output file name.
  type(varying_string) :: outputFileName, outputScratchFileName

contains

  subroutine Output_HDF5_Open_File(parameters)
    !!{
    Open the file for \glc\ output.
    !!}
    use :: Output_HDF5       , only : hdf5SieveBufferSize , hdf5UseLatestFormat , hdf5CompressionLevel, hdf5CacheElementsCount, &
         &                            outputFileIsOpen    , outputFile, hdf5CacheSizeBytes  , hdf5ChunkSize
    use :: HDF5              , only : hsize_t             , size_t
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : IO_HDF5_Set_Defaults
    use :: ISO_Varying_String, only : var_str             , char                , operator(//)        , extract               , &
         &                            len                 , operator(==)
    use :: Input_Parameters  , only : inputParameters     , inputParameter
#ifdef USEMPI
    use :: MPI_Utilities     , only : mpiSelf
#endif
    use :: String_Handling   , only : operator(//)
    !![
    <include directive="outputFileOpenTask" type="moduleUse">
    !!]
    include 'output.open.modules.inc'
    !![
    </include>
    !!]
    implicit none
    type   (inputParameters), intent(inout) :: parameters
    integer(hsize_t        )                :: chunkSize
    integer                                 :: sieveBufferSize
    integer(size_t         )                :: cacheElementsCount, cacheSizeBytes
#ifdef USEMPI
    type   (varying_string )                :: fileNamePrefix
#endif
    
    if (.not.outputFileIsOpen) then
       !![
       <inputParameter>
         <name>outputFileName</name>
         <defaultValue>var_str('galacticus.hdf5')</defaultValue>
         <description>The name of the file to which \glc\ results will be written.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>outputScratchFileName</name>
         <defaultValue>outputFileName</defaultValue>
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
       !$ call hdf5Access%set()
       call outputFile%openFile(                                                 &
            &                                       char(outputScratchFileName), &
            &                   overWrite          =.true.                     , &
            &                   objectsOverwritable=.true.                     , &
            &                   sieveBufferSize    =hdf5SieveBufferSize        , &
            &                   useLatestFormat    =hdf5UseLatestFormat        , &
            &                   cacheElementsCount =hdf5CacheElementsCount     , &
            &                   cacheSizeBytes     =hdf5CacheSizeBytes           &
            &                  )
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

       ! Call all routines that requested to output to the file on start up.
       !![
       <include directive="outputFileOpenTask" type="functionCall" functionType="void">
       !!]
       include 'output.open.inc'
       !![
       </include>
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
    use :: Output_HDF5       , only : outputFileIsOpen, outputFile
    use :: File_Utilities    , only : File_Rename
    use :: HDF5_Access       , only : hdf5Access
    use :: ISO_Varying_String, only : operator(/=)
    !![
    <include directive="hdfPreCloseTask" type="moduleUse">
    !!]
    include 'output.HDF5.pre_close_tasks.moduleUse.inc'
    !![
    </include>
    !!]
    implicit none

    ! Perform any final tasks prior to shutdown.
    if (outputFileIsOpen) then
       !$omp critical (Output_HDF5_Close_File)
       if (outputFileIsOpen) then
          !![
          <include directive="hdfPreCloseTask" type="functionCall" functionType="void">
          !!]
          include 'output.HDF5.pre_close_tasks.inc'
          !![
          </include>
          <eventHook name="hdf5PreClose"/>
          !!]
          ! Close the file.
          !$ call hdf5Access%set()
          call outputFile%writeAttribute(1,"galacticusCompleted")
          call outputFile%close()
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

end module Output_HDF5_Open
