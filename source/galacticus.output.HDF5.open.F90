!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use ISO_Varying_String
  use Galacticus_HDF5
  use IO_HDF5
  use Galacticus_Error
  use HDF5
  implicit none
  private
  public :: Galacticus_Output_Open_File, Galacticus_Output_Close_File

  ! Output file name.
  type(varying_string) :: galacticusOutputFileName

contains

  subroutine Galacticus_Output_Open_File
    !% Open the file for \glc\ output.
    use Input_Parameters
    !# <include directive="outputFileOpenTask" type="moduleUse">
    include 'galacticus.output.open.modules.inc'
    !# </include>
    implicit none
    integer :: chunkSize
    
    if (.not.galacticusOutputFileIsOpen) then
       ! Get file name parameter.
       !@ <inputParameter>
       !@   <name>galacticusOutputFileName</name>
       !@   <defaultValue>galacticus.hdf5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the file to which \glc\ results will be written.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('galacticusOutputFileName',galacticusOutputFileName,defaultValue='galacticus.hdf5',writeOutput=.false.)
       ! Open the file.
       !$omp critical(HDF5_Access)
       call galacticusOutputFile%openFile(char(galacticusOutputFileName),overWrite=.true.,objectsOverwritable=.false.)
       !$omp end critical(HDF5_Access)
       ! Get file name parameter again and write it to the output file.
       call Get_Input_Parameter('galacticusOutputFileName',galacticusOutputFileName,defaultValue='galacticus.hdf5')
       ! Read parameters.
       !@ <inputParameter>
       !@   <name>hdf5ChunkSize</name>
       !@   <defaultValue>1024</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The chunk size used for outputting HDF5 datasets.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('hdf5ChunkSize',chunksize,defaultValue=1024)
       hdf5ChunkSize=chunksize
       !@ <inputParameter>
       !@   <name>hdf5CompressionLevel</name>
       !@   <defaultValue>-1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The compression level used for outputting HDF5 datasets.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('hdf5CompressionLevel',hdf5CompressionLevel,defaultValue=-1)

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
    !# <include directive="hdfPreCloseTask" type="moduleUse">
    include 'galacticus.output.HDF5.pre_close_tasks.moduleUse.inc'
    !# </include>
    implicit none

    ! Perform any final tasks prior to shutdown.
    !# <include directive="hdfPreCloseTask" type="functionCall" functionType="void">
    include 'galacticus.output.HDF5.pre_close_tasks.inc'
    !# </include>

    ! Close the file.
    !$omp critical(HDF5_Access)
    call galacticusOutputFile%close()
    !$omp end critical(HDF5_Access)
    return
  end subroutine Galacticus_Output_Close_File

end module Galacticus_Output_Open
