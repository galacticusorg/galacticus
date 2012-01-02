!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
       !@   <defaultValue>9</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The compression level used for outputting HDF5 datasets.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('hdf5CompressionLevel',hdf5CompressionLevel,defaultValue=9)

       ! Set default chunking and compression levels.
       call IO_HDF5_Set_Defaults(hdf5ChunkSize,hdf5CompressionLevel)

       ! Call all routines that requested to output to the file on start up.
       !# <include directive="outputFileOpenTask" type="code" action="subroutine">
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
    !# <include directive="hdfPreCloseTask" type="code" action="subroutine">
    include 'galacticus.output.HDF5.pre_close_tasks.inc'
    !# </include>

    ! Close the file.
    !$omp critical(HDF5_Access)
    call galacticusOutputFile%close()
    !$omp end critical(HDF5_Access)
    return
  end subroutine Galacticus_Output_Close_File

end module Galacticus_Output_Open
