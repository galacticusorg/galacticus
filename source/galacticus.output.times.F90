!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which provides output times.

module Galacticus_Output_Times
  !% Provides output times.
  implicit none
  private
  public :: Galacticus_Output_Time_Count, Galacticus_Output_Time         , &
       &    Galacticus_Next_Output_Time , Galacticus_Previous_Output_Time, &
       &    Galacticus_Output_Time_Index

  ! Flag to indicate if output times have been initialized.
  logical                                     :: outputsInitialized=.false.

  ! Array of output times.
  integer                                     :: outputCount
  double precision, allocatable, dimension(:) :: outputTimes

contains

  subroutine Output_Times_Initialize()
    !% Initialize the output times.
    use Input_Parameters
    use Sort
    use Memory_Management
    use Histories
    use Cosmology_Functions
    implicit none
    class           (cosmologyFunctionsClass), pointer :: cosmologyFunctionsDefault
    integer                                            :: iOutput
    double precision                                   :: aExpansion

    if (.not.outputsInitialized) then
       !$omp critical (Tasks_Evolve_Tree_Initialize)
       if (.not.outputsInitialized) then
          ! Get a list of output redshifts - stored temporarily in the outputTimes array.
          outputCount=max(Get_Input_Parameter_Array_Size('outputRedshifts'),1)
          call Alloc_Array(outputTimes,[outputCount])
          !@ <inputParameter>
          !@   <name>outputRedshifts</name>
          !@   <defaultValue>0</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     A list of (space-separated) redshifts at which \glc\ results should be output. Redshifts need not be in any particular order.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1..*</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          if (outputCount == 1) then
             ! If only one (or zero) output redshifts present, make redshift zero the default.
             call Get_Input_Parameter('outputRedshifts',outputTimes,defaultValue=[0.0d0])
          else
             call Get_Input_Parameter('outputRedshifts',outputTimes                     )
          end if
          ! Get the default cosmology functions object.
          cosmologyFunctionsDefault => cosmologyFunctions()
          ! Convert redshifts to times.
          do iOutput=1,outputCount
             aExpansion=cosmologyFunctionsDefault%expansionFactorFromRedshift(outputTimes(iOutput))
             outputTimes(iOutput)=cosmologyFunctionsDefault%cosmicTime(aExpansion)
          end do

          ! Sort the times.
          call Sort_Do(outputTimes)

          ! Set history ranges to include these times.
          call History_Set_Times(timeEarliest=outputTimes(1),timeLatest=outputTimes(outputCount))

          ! Flag that this module is now initialized.
          outputsInitialized=.true.
       end if
       !$omp end critical (Tasks_Evolve_Tree_Initialize)
    end if
    return
  end subroutine Output_Times_Initialize

  integer function Galacticus_Output_Time_Count()
    !% Return the number of outputs.
    implicit none

    ! Ensure the module is initialized.
    call Output_Times_Initialize()

    ! Return the number of outputs.
    Galacticus_Output_Time_Count=outputCount
    return
  end function Galacticus_Output_Time_Count

  double precision function Galacticus_Output_Time(iOutput)
    !% Returns the time of the output indexed by {\tt iOutput}.
    implicit none
    integer, intent(in   ) :: iOutput

    ! Ensure the module is initialized.
    call Output_Times_Initialize()

    ! Return the requested output time.
    if (iOutput >=1 .and. iOutput <= outputCount) then
       Galacticus_Output_Time=outputTimes(iOutput)
    else
       Galacticus_Output_Time=-1.0d0
    end if
    return
  end function Galacticus_Output_Time

  integer function Galacticus_Output_Time_Index(time)
    !% Returns the index of the output given the corresponding time.
    use Galacticus_Error
    use Numerical_Comparison
    use Arrays_Search
    implicit none
    double precision, intent(in   ) :: time

    ! Ensure the module is initialized.
    call Output_Times_Initialize()

    Galacticus_Output_Time_Index=Search_Array(outputTimes,time)
    if (Values_Differ(time,outputTimes(Galacticus_Output_Time_Index),relTol=1.0d-6)) &
         & call Galacticus_Error_Report('Galacticus_Output_Time_Index','time does not correspond to an output')
    return
  end function Galacticus_Output_Time_Index

  double precision function Galacticus_Next_Output_Time(currentTime,outputIndex)
    !% Returns the time of the next output after {\tt currentTime}.
    use Arrays_Search
    implicit none
    double precision, intent(in   )           :: currentTime
    integer         , intent(  out), optional :: outputIndex
    integer                                   :: i

    ! Ensure the module is initialized.
    call Output_Times_Initialize()

    ! If the current time exceeds the last output, return an unphysical value.
    if      (currentTime >= outputTimes(outputCount)) then
       Galacticus_Next_Output_Time=-1.0d0
       if (present(outputIndex)) outputIndex=-1
    else if (currentTime <  outputTimes(          1)) then
       Galacticus_Next_Output_Time=outputTimes(                                      1)
       if (present(outputIndex)) outputIndex=+1
    else
       i=Search_Array(outputTimes,currentTime)+1
       Galacticus_Next_Output_Time=outputTimes(i)
       if (present(outputIndex)) outputIndex=i
    end if
    return
  end function Galacticus_Next_Output_Time

  double precision function Galacticus_Previous_Output_Time(currentTime)
    !% Returns the time of the previous output prior to {\tt currentTime}.
    use Arrays_Search
    implicit none
    double precision, intent(in   ) :: currentTime

    ! Ensure the module is initialized.
    call Output_Times_Initialize()

    if      (currentTime > outputTimes(outputCount)) then
       ! If the current time exceeds the last output, return the last output.
       Galacticus_Previous_Output_Time=outputTimes(outputCount)
    else if (currentTime < outputTimes(          1)) then
       ! If the current time preceeds the first output, return and unphysical value.
       Galacticus_Previous_Output_Time=-1.0d0
    else
       Galacticus_Previous_Output_Time=outputTimes(Search_Array(outputTimes,currentTime))
    end if
    return
  end function Galacticus_Previous_Output_Time

end module Galacticus_Output_Times
