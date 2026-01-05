!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Contains a module which implements storage and recovery of the Galacticus internal state. Used for restoring random number
generator sequences for example.
!!}

module State
  !!{
  Implements storage and recovery of the Galacticus internal state. Used for restoring random number
  generator sequences for example.
  !!}
  use, intrinsic :: ISO_C_Binding     , only : c_size_t
  use            :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: State_Store, State_Retrieve, State_Initialize, State_Set

  ! Flag indicating if we have retrieved the internal state already.
  logical                 :: stateHasBeenRetrieved=.false.

  ! Root name for state files.
  type   (varying_string) :: stateFileRoot                  , stateRetrieveFileRoot

  ! Active status of store and retrieve.
  logical                 :: stateStoreActive               , stateRetrieveActive

  ! Counter which tracks state operators, used to ensure objects are stored to file only once per operation.
  integer(c_size_t      ) :: stateOperatorID      =0_c_size_t

contains

  !![
  <functionGlobal>
   <unitName>State_Store</unitName>
   <type>void</type>
   <module>ISO_Varying_String, only : varying_string</module>
   <arguments>type(varying_string) , intent(in   ), optional :: logMessage</arguments>
  </functionGlobal>
  !!]
  subroutine State_Store(logMessage)
    !!{
    Store the internal state.
    !!}
#ifdef USEMPI
    use            :: MPI_Utilities     , only : mpiSelf
#endif
    !$ use         :: OMP_Lib           , only : omp_get_thread_num, omp_in_parallel
    use            :: Interface_GSL     , only : gslFileOpen       , gslFileClose
    use, intrinsic :: ISO_C_Binding     , only : c_ptr
    use            :: ISO_Varying_String, only : operator(//)      , char
    use            :: String_Handling   , only : operator(//)
    !![
    <include directive="stateStoreTask" type="moduleUse">
    !!]
    include 'state.store.modules.inc'
    !![
    </include>
    !!]
    implicit none
    type   (varying_string), intent(in   ), optional :: logMessage
    integer                                          :: stateUnit
    integer(c_size_t      )                          :: stateOperatorID_
    type   (c_ptr         )                          :: gslStateFile
    type   (varying_string)                          :: fileName        , fileNameGSL , &
         &                                              fileNameLog     , fileName_   , &
         &                                              fileNameGSL_    , fileNameLog_

    ! Check if state store is active.
    if (stateStoreActive) then
       ! Open a file in which to store the state and an additional file for GSL state.
       fileName   =stateFileRoot//'.state'
       fileNameGSL=stateFileRoot//'.gsl.state'
       fileNameLog=stateFileRoot//'.state.log'
       !$ if (omp_in_parallel()) then
       !$    fileName_   =fileName   //':openMP'//omp_get_thread_num()
       !$    fileNameGSL_=fileNameGSL//':openMP'//omp_get_thread_num()
       !$    fileNameLog_=fileNameLog//':openMP'//omp_get_thread_num()
       !$    fileName    =fileName_
       !$    fileNameGSL =fileNameGSL_
       !$    fileNameLog =fileNameLog_
       !$ end if
#ifdef USEMPI
       fileName_   =fileName   //':MPI'//mpiSelf%rankLabel()
       fileNameGSL_=fileNameGSL//':MPI'//mpiSelf%rankLabel()
       fileNameLog_=fileNameLog//':MPI'//mpiSelf%rankLabel()
       fileName    =fileName_
       fileNameGSL =fileNameGSL_
       fileNameLog =fileNameLog_
#endif
       if (present(logMessage)) then
          open(newunit=stateUnit,file=char(fileNameLog),form='formatted',status='unknown',access='append')
          write (stateUnit,*) char(logMessage)
          close(stateUnit)
       end if

       open(newunit=stateUnit,file=char(fileName),form='unformatted',status='unknown')
       gslStateFile=gslFileOpen(char(fileNameGSL),'w')

       !$omp critical(stateOperationID)
       stateOperatorID =stateOperatorID+1_c_size_t
       stateOperatorID_=stateOperatorID
       !$omp end critical(stateOperationID)
       !![
       <include directive="stateStoreTask" type="functionCall" functionType="void">
        <functionArgs>stateUnit,gslStateFile,stateOperatorID_</functionArgs>
       !!]
       include 'state.store.inc'
       !![
       </include>
       <eventHook name="stateStore">
        <callWith>stateUnit,gslStateFile,stateOperatorID_</callWith>
       </eventHook>
       !!]

       ! Close the state files.
       close(stateUnit)
       call gslFileClose(gslStateFile)

       ! Flush standard output to ensure that any output log has a record of where the code reached at the last state store.
       call Flush(0)

    end if
    return
  end subroutine State_Store

  !![
  <functionGlobal>
   <unitName>State_Retrieve</unitName>
   <type>void</type>
  </functionGlobal>
  !!]
  subroutine State_Retrieve
    !!{
    Retrieve the internal state.
    !!}
    use            :: Interface_GSL     , only : gslFileOpen , gslFileClose
#ifdef USEMPI
    use            :: MPI_Utilities     , only : mpiSelf
#endif
    !$ use :: OMP_Lib           , only : omp_get_thread_num  , omp_in_parallel
    use, intrinsic :: ISO_C_Binding     , only : c_ptr
    use            :: ISO_Varying_String, only : operator(//), char
    use            :: String_Handling   , only : operator(//)
    !![
    <include directive="stateRetrieveTask" type="moduleUse">
    !!]
    include 'state.retrieve.modules.inc'
    !![
    </include>
    !!]
    implicit none
    integer                 :: stateUnit
    integer(c_size_t      ) :: stateOperatorID_
    type   (c_ptr         ) :: gslStateFile
    type   (varying_string) :: fileName        , fileNameGSL , &
         &                     fileName_       , fileNameGSL_

    ! Check if we have already retrieved the internal state.
    if (.not.stateHasBeenRetrieved) then
       !$omp critical (stateRetrieve)
       if (.not.stateHasBeenRetrieved) then
          
          ! Check if state retrieve is active.
          if (stateRetrieveActive) then
             
             ! Open a file in which to retrieve the state and an additional file for GSL state.
             fileName   =stateRetrieveFileRoot//'.state'
             fileNameGSL=stateRetrieveFileRoot//'.gsl.state'
             !$ if (omp_in_parallel()) then
             !$    fileName_   =fileName    //':openMP'//omp_get_thread_num()
             !$    fileNameGSL_=fileNameGSL //':openMP'//omp_get_thread_num()
             !$    fileName    =fileName_
             !$    fileNameGSL =fileNameGSL_
             !$ end if
#ifdef USEMPI
             fileName_   =fileName   //':MPI'//mpiSelf%rankLabel()
             fileNameGSL_=fileNameGSL//':MPI'//mpiSelf%rankLabel()
             fileName    =fileName_
             fileNameGSL =fileNameGSL_
#endif
             open(newunit=stateUnit,file=char(fileName),form='unformatted',status='old')
             gslStateFile=gslFileOpen(char(fileNameGSL),'r')
             
             !$omp critical(stateOperationID)
             stateOperatorID =stateOperatorID+1_c_size_t
             stateOperatorID_=stateOperatorID
             !$omp end critical(stateOperationID)
             !![
             <include directive="stateRetrieveTask" type="functionCall" functionType="void">
              <functionArgs>stateUnit,gslStateFile,stateOperatorID_</functionArgs>
             !!]
	     include 'state.retrieve.inc'
             !![
             </include>
             <eventHook name="stateRestore">
              <callWith>stateUnit,gslStateFile,stateOperatorID_</callWith>
             </eventHook>
             !!]
      
             ! Close the state files.
             close(stateUnit)
             call gslFileClose(gslStateFile)

          end if
          
          ! Flag that internal state has been retrieved
          stateHasBeenRetrieved=.true.
       end if
       !$omp end critical (stateRetrieve)
    end if

    return
  end subroutine State_Retrieve

  !![
  <nodeComponentInitializationTask>
   <unitName>State_Initialize</unitName>
   <useGlobal>yes</useGlobal>
  </nodeComponentInitializationTask>
  <functionGlobal>
   <unitName>State_Initialize</unitName>
   <type>void</type>
   <module>Input_Parameters  , only : inputParameters</module>
   <arguments>type(inputParameters), intent(inout) :: parameters_</arguments>
  </functionGlobal>
  !!]
  subroutine State_Initialize(parameters_)
    !!{
    Initialize the state module by getting the name of the file to which states should be stored and whether or not we are to
    retrieve a state.
    !!}
    use :: Input_Parameters  , only : inputParameters
    use :: ISO_Varying_String, only : var_str        , operator(/=)
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    ! Get the base name of the state files.
    !![
    <inputParameter>
      <name>stateFileRoot</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>The root name of files to which the internal state is written (to permit restarts).</description>
      <source>parameters_</source>
    </inputParameter>
    !!]
    ! Get the base name of the files to retrieve from.
    !![
    <inputParameter>
      <name>stateRetrieveFileRoot</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>The root name of files to which the internal state is retrieved from (to restart).</description>
      <source>parameters_</source>
    </inputParameter>
    !!]
    ! Record active status of store and retrieve.
    stateStoreActive   =(stateFileRoot         /= "none")
    stateRetrieveActive=(stateRetrieveFileRoot /= "none")
    return
  end subroutine State_Initialize

  !![
  <functionGlobal>
    <unitName>State_Set</unitName>
    <type>void</type>
    <module>ISO_Varying_String, only : varying_string</module>
    <arguments>type(varying_string) , intent(in   ) :: stateFileRoot_</arguments>
  </functionGlobal>
  !!]
  subroutine State_Set(stateFileRoot_)
    !!{
    Set the state system---can be used for force storing of state on prior to calls to state store functions for debugging purposes.
    !!}
    use :: ISO_Varying_String, only : operator(/=)
    implicit none
    type(varying_string), intent(in   ) :: stateFileRoot_

    stateFileRoot   = stateFileRoot_
    stateStoreActive=(stateFileRoot /= "none")
    return
  end subroutine State_Set

end module State
