!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements storage and recovery of the Galacticus internal state. Used for restoring random number
!% generator sequences for example.

module Galacticus_State
  !% Implements storage and recovery of the Galacticus internal state. Used for restoring random number
  !% generator sequences for example.
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_State_Snapshot, Galacticus_State_Store, Galacticus_State_Retrieve

  ! Flag indicating if we have retrieved the internal state already.
  logical              :: stateHasBeenRetrieved=.false.

  ! Root name for state files.
  type(varying_string) :: stateFileRoot,stateRetrieveFileRoot

  ! Flag indicating if module has been initialized.
  logical              :: stateInitialized=.false.

contains

  subroutine Galacticus_State_Snapshot
    !% Take a snapshot of the internal state.
    !# <include directive="galacticusStateSnapshotTask" type="moduleUse">
    include 'galacticus.state.snapshot.modules.inc'
    !# </include>
    implicit none

    ! Ensure that module is initialized.
    call State_Initialize

    !# <include directive="galacticusStateSnapshotTask" type="code" action="subroutine">
    include 'galacticus.state.snapshot.inc'
    !# </include>
    return
  end subroutine Galacticus_State_Snapshot

  subroutine Galacticus_State_Store(logMessage)
    !% Store the internal state.
    !$ use OMP_Lib
    use String_Handling
    use FGSL
    !# <include directive="galacticusStateStoreTask" type="moduleUse">
    include 'galacticus.state.store.modules.inc'
    !# </include>
    implicit none
    type(varying_string), intent(in), optional :: logMessage
    integer                                    :: stateUnit,iError
    type(fgsl_file)                            :: fgslStateFile
    type(varying_string)                       :: fileName,fileNameFGSL,fileNameLog

    ! Ensure that module is initialized.
    call State_Initialize

    ! Check if a file has been specified.
    if (stateFileRoot /= "none") then
       
       ! Open a file in which to store the state and an additional file for FGSL state.
       !$ if (omp_in_parallel()) then
       !$    fileName    =stateFileRoot//     '.state.'
       !$    fileNameFGSL=stateFileRoot//'.fgsl.state.'
       !$    fileNameLog =stateFileRoot//     '.state.log.'
       !$    fileName    =fileName    //omp_get_thread_num()
       !$    fileNameFGSL=fileNameFGSL//omp_get_thread_num()
       !$    fileNameLog =fileNameLog //omp_get_thread_num()
       !$ else
       fileName    =stateFileRoot//'.state'
       fileNameFGSL=stateFileRoot//'.fgsl.state'
       fileNameLog =stateFileRoot//'.state.log'
       !$ end if
       if (present(logMessage)) then
          open(newunit=stateUnit,file=char(fileNameLog),form='formatted',status='unknown',access='append')
          write (stateUnit,*) char(logMessage)
          close(stateUnit)
       end if

       open(newunit=stateUnit,file=char(fileName),form='unformatted',status='unknown')
       fgslStateFile=FGSL_Open(char(fileNameFGSL),'w')
       
       !# <include directive="galacticusStateStoreTask" type="code" action="subroutine">
       !#  <subroutineArgs>stateUnit,fgslStateFile</subroutineArgs>
       include 'galacticus.state.store.inc'
       !# </include>
       
       ! Close the state files.
       close(stateUnit)
       iError=FGSL_Close(fgslStateFile)

       ! Flush standard output to ensure that any output log has a record of where the code reached at the last state store.
       call Flush(0)

    end if
    return
  end subroutine Galacticus_State_Store
  
  subroutine Galacticus_State_Retrieve
    !% Retrieve the interal state.
    !$ use OMP_Lib
    use String_Handling
    use FGSL
    !# <include directive="galacticusStateRetrieveTask" type="moduleUse">
    include 'galacticus.state.retrieve.modules.inc'
    !# </include>
    implicit none
    integer              :: stateUnit,iError
    type(fgsl_file)      :: fgslStateFile
    type(varying_string) :: fileName,fileNameFGSL

    ! Check if we have already retrieved the internal state.
    if (.not.stateHasBeenRetrieved) then

       ! Ensure that module is initialized.
       call State_Initialize

       ! Check if a file has been specified.
       if (stateRetrieveFileRoot /= "none") then
          
          ! Open a file in which to retrieve the state and an additional file for FGSL state.
          !$ if (omp_in_parallel()) then
          !$    fileName    =stateRetrieveFileRoot//     '.state.'
          !$    fileNameFGSL=stateRetrieveFileRoot//'.fgsl.state.'
          !$    fileName    =fileName    //omp_get_thread_num()
          !$    fileNameFGSL=fileNameFGSL//omp_get_thread_num()
          !$ else
          fileName    =stateRetrieveFileRoot//'.state'
          fileNameFGSL=stateRetrieveFileRoot//'.fgsl.state'
          !$ end if
          open(newunit=stateUnit,file=char(fileName),form='unformatted',status='old')
          fgslStateFile=FGSL_Open(char(fileNameFGSL),'r')
          
          !# <include directive="galacticusStateRetrieveTask" type="code" action="subroutine">
          !#  <subroutineArgs>stateUnit,fgslStateFile</subroutineArgs>
          include 'galacticus.state.retrieve.inc'
          !# </include>
          
          ! Close the state files.
          close(stateUnit)
          iError=FGSL_Close(fgslStateFile)
          
       end if
       
       ! Flag that internal state has been retrieved
       stateHasBeenRetrieved=.true.
    end if
    
    return
  end subroutine Galacticus_State_Retrieve

  subroutine State_Initialize
    !% Initialize the state module by getting the name of the file to which states should be stored and whether or not we are to
    !% retrieve a state.
    use Input_Parameters
    implicit none
  
    if (.not.stateInitialized) then
       
       ! Get the base name of the state files.
       !@ <inputParameter>
       !@   <name>stateFileRoot</name>
       !@   <defaultValue>none</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The root name of files to which the internal state is written (to permit restarts).
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stateFileRoot'        ,stateFileRoot        ,defaultValue="none")
   
       ! Get the base name of the files to retrieve from.
       !@ <inputParameter>
       !@   <name>stateRetrieveFileRoot</name>
       !@   <defaultValue>none</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The root name of files to which the internal state is retrieved from (to restart).
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stateRetrieveFileRoot',stateRetrieveFileRoot,defaultValue="none")
       
       ! Flag that module is now initialized.
       stateInitialized=.true. 
       
    end if
    return
  end subroutine State_Initialize
  
end module Galacticus_State
