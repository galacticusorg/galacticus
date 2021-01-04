!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which handles opening of the \glc\ output file.

module Galacticus_Output_Open
  !% Handles opening of the \glc\ output file.
  use :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: Galacticus_Output_Open_File, Galacticus_Output_Close_File

  ! Output file name.
  type(varying_string) :: galacticusOutputFileName, galacticusOutputScratchFileName

contains

  subroutine Galacticus_Output_Open_File(parameters)
    !% Open the file for \glc\ output.
    use :: Galacticus_HDF5   , only : hdf5SieveBufferSize       , hdf5UseLatestFormat , hdf5CompressionLevel, hdf5CacheElementsCount, &
         &                            galacticusOutputFileIsOpen, galacticusOutputFile, hdf5CacheSizeBytes  , hdf5ChunkSize
    use :: HDF5              , only : hsize_t                   , size_t
    use :: IO_HDF5           , only : hdf5Access                , IO_HDF5_Set_Defaults
    use :: ISO_Varying_String, only : var_str                   , char                , operator(//)        , extract               , &
         &                            len                       , operator(==)
    use :: Input_Parameters  , only : inputParameters           , inputParameter
#ifdef USEMPI
    use :: MPI_Utilities     , only : mpiSelf
#endif
    use :: String_Handling   , only : operator(//)
    !# <include directive="outputFileOpenTask" type="moduleUse">
    include 'galacticus.output.open.modules.inc'
    !# </include>
    implicit none
    type   (inputParameters), intent(inout) :: parameters
    integer(hsize_t        )                :: chunkSize
    integer                                 :: sieveBufferSize
    integer(size_t         )                :: cacheElementsCount, cacheSizeBytes

    if (.not.galacticusOutputFileIsOpen) then
       !# <inputParameter>
       !#   <name>galacticusOutputFileName</name>
       !#   <defaultValue>var_str('galacticus.hdf5')</defaultValue>
       !#   <description>The name of the file to which \glc\ results will be written.</description>
       !#   <source>parameters</source>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>galacticusOutputScratchFileName</name>
       !#   <defaultValue>galacticusOutputFileName</defaultValue>
       !#   <description>The name of the file to which \glc\ results will be written temporarily during runs.</description>
       !#   <source>parameters</source>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>hdf5SieveBufferSize</name>
       !#   <defaultValue>65536</defaultValue>
       !#   <description>The size of the sieve buffer used by the HDF5 library to speed reading/writing of partial datasets.</description>
       !#   <source>parameters</source>
       !#   <variable>sieveBufferSize</variable>
       !# </inputParameter>
       hdf5SieveBufferSize=sieveBufferSize
       !# <inputParameter>
       !#   <name>hdf5UseLatestFormat</name>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>Specifies whether to use the latest HDF5 file format.</description>
       !#   <source>parameters</source>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>hdf5CacheElementsCount</name>
       !#   <defaultValue>521_size_t</defaultValue>
       !#   <description>Number of elements in the raw data chunk cache.</description>
       !#   <source>parameters</source>
       !#   <variable>cacheElementsCount</variable>
       !# </inputParameter>
       hdf5CacheElementsCount=cacheElementsCount
       !# <inputParameter>
       !#   <name>hdf5CacheSizeBytes</name>
       !#   <defaultValue>1048576_size_t</defaultValue>
       !#   <description>Size of the raw data chunk cache in bytes.</description>
       !#   <source>parameters</source>
       !#   <variable>cacheSizeBytes</variable>
       !# </inputParameter>
       hdf5CacheSizeBytes=cacheSizeBytes
       ! Modify the file name on a per-process basis if running under MPI.
#ifdef USEMPI
       if (extract(galacticusOutputFileName       ,len(galacticusOutputFileName       )-4,len(galacticusOutputFileName       )) == ".hdf5") then
          galacticusOutputFileName       =extract(galacticusOutputFileName       ,1,len(galacticusOutputFileName       )-5)//':MPI'//mpiSelf%rankLabel()//'.hdf5'
       else
          galacticusOutputFileName       =galacticusOutputFileName       //':MPI'//mpiSelf%rankLabel()
       end if
       if (extract(galacticusOutputScratchFileName,len(galacticusOutputScratchFileName)-4,len(galacticusOutputScratchFileName)) == ".hdf5") then
          galacticusOutputScratchFileName=extract(galacticusOutputScratchFileName,1,len(galacticusOutputScratchFileName)-5)//':MPI'//mpiSelf%rankLabel()//'.hdf5'
       else
          galacticusOutputScratchFileName=galacticusOutputScratchFileName//':MPI'//mpiSelf%rankLabel()
       end if
#endif
       ! Open the file.
       !$ call hdf5Access%set()
       call galacticusOutputFile%openFile(                                            &
            &                             char(galacticusOutputScratchFileName)     , &
            &                             overWrite          =.true.                , &
            &                             objectsOverwritable=.true.                , &
            &                             sieveBufferSize    =hdf5SieveBufferSize   , &
            &                             useLatestFormat    =hdf5UseLatestFormat   , &
            &                             cacheElementsCount =hdf5CacheElementsCount, &
            &                             cacheSizeBytes     =hdf5CacheSizeBytes      &
            &                            )
       !$ call hdf5Access%unset()
       ! Now that the parameter file is open, we can open an output group in it for parameters.
       call parameters%parametersGroupOpen(galacticusOutputFile)
       ! Read parameters.
       !# <inputParameter>
       !#   <name>hdf5ChunkSize</name>
       !#   <defaultValue>1024_hsize_t</defaultValue>
       !#   <description>The chunk size used for outputting HDF5 datasets.</description>
       !#   <source>parameters</source>
       !#   <variable>chunksize</variable>
       !# </inputParameter>
       hdf5ChunkSize=chunksize
       !# <inputParameter>
       !#   <name>hdf5CompressionLevel</name>
       !#   <defaultValue>-1</defaultValue>
       !#   <description>The compression level used for outputting HDF5 datasets.</description>
       !#   <source>parameters</source>
       !# </inputParameter>

       ! Set default chunking and compression levels.
       call IO_HDF5_Set_Defaults(hdf5ChunkSize,hdf5CompressionLevel)

       ! Call all routines that requested to output to the file on start up.
       !# <include directive="outputFileOpenTask" type="functionCall" functionType="void">
       include 'galacticus.output.open.inc'
       !# </include>

       ! Flag that the file is now open.
       galacticusOutputFileIsOpen=.true.
    end if
    return
  end subroutine Galacticus_Output_Open_File

  subroutine Galacticus_Output_Close_File
    !% Close the \glc\ output file.
    use :: Galacticus_HDF5   , only : galacticusOutputFileIsOpen, galacticusOutputFile
    use :: File_Utilities    , only : File_Rename
    use :: IO_HDF5           , only : hdf5Access
    use :: ISO_Varying_String, only : operator(/=)
    !# <include directive="hdfPreCloseTask" type="moduleUse">
    include 'galacticus.output.HDF5.pre_close_tasks.moduleUse.inc'
    !# </include>
    implicit none

    ! Perform any final tasks prior to shutdown.
    if (galacticusOutputFileIsOpen) then
       !$omp critical (Galacticus_Output_Close_File)
       if (galacticusOutputFileIsOpen) then
          !# <include directive="hdfPreCloseTask" type="functionCall" functionType="void">
          include 'galacticus.output.HDF5.pre_close_tasks.inc'
          !# </include>
          !# <eventHook name="hdf5PreClose"/>
          ! Close the file.
          !$ call hdf5Access%set()
          call galacticusOutputFile%writeAttribute(1,"galacticusCompleted")
          call galacticusOutputFile%close()
          !$ call hdf5Access%unset()
          ! Move the scratch file to the final file if necessary.
          if (galacticusOutputFileName /= galacticusOutputScratchFileName) call File_Rename(galacticusOutputScratchFileName,galacticusOutputFileName)
          ! Record that the file is now closed.
          galacticusOutputFileIsOpen=.false.
       end if
       !$omp end critical (Galacticus_Output_Close_File)
    end if
    return
  end subroutine Galacticus_Output_Close_File

end module Galacticus_Output_Open
