!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
  private
  public :: Galacticus_Output_Open_File, Galacticus_Output_Close_File

  ! Output file name.
  type(varying_string) :: galacticusOutputFile

contains

  subroutine Galacticus_Output_Open_File
    !% Open the file for \glc\ output.
    use Input_Parameters
    !# <include directive="outputFileOpenTask" type="moduleUse">
    include 'galacticus.output.open.modules.inc'
    !# </include>
    implicit none
    integer :: errorCode
    
    if (galacticusOutputID < 0) then
       ! Ensure HDF5 system is initialized.
       call IO_HDF5_Initialize
       ! Get file name parameter.
       !@ <inputParameter>
       !@   <name>galacticusOutputFile</name>
       !@   <defaultValue>galacticus.hdf5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the file to which \glc\ results will be written.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('galacticusOutputFile',galacticusOutputFile,defaultValue='galacticus.hdf5',writeOutput=.false.)
       ! Open the file.
       call h5fcreate_f(char(galacticusOutputFile),H5F_ACC_TRUNC_F,galacticusOutputID,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Merger_Tree_Output_Open_File','failed to open output file')
       ! Get file name parameter again and write it to the output file.
       call Get_Input_Parameter('galacticusOutputFile',galacticusOutputFile,defaultValue='galacticus.hdf5')

       ! Call all routines that requested to output to the file on start up.
       !# <include directive="outputFileOpenTask" type="code" action="subroutine">
       include 'galacticus.output.open.inc'
       !# </include>

    end if
    return
  end subroutine Galacticus_Output_Open_File

  subroutine Galacticus_Output_Close_File
    !% Close the \glc\ output file.
    !# <include directive="hdfPreCloseTask" type="moduleUse">
    include 'galacticus.output.HDF5.pre_close_tasks.moduleUse.inc'
    !# </include>
    implicit none
    integer :: errorCode

    ! Perform any final tasks prior to shutdown.
    !# <include directive="hdfPreCloseTask" type="code" action="subroutine">
    include 'galacticus.output.HDF5.pre_close_tasks.inc'
    !# </include>

    call h5fclose_f(galacticusOutputID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Close_File','failed to close output file')
    call IO_HDF5_Uninitialize
    return
  end subroutine Galacticus_Output_Close_File

end module Galacticus_Output_Open
