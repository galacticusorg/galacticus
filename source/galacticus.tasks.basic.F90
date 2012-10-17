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

module Galacticus_Tasks_Basic
  use Galacticus_Display
  implicit none
  private
  public :: Galacticus_Task_Start, Galacticus_Task_End
  
contains

  !# <galacticusTask>
  !#  <unitName>Galacticus_Task_Start</unitName>
  !# </galacticusTask>
  logical function Galacticus_Task_Start()
    use Input_Parameters
    use Galacticus_Output_Open
    implicit none
    logical, save :: doneStart=.false.
    integer       :: verbosityLevel

    if (.not.doneStart) then
       ! Get the verbosity level parameter.
       !@ <inputParameter>
       !@   <name>verbosityLevel</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The level of verbosity for \glc\ (higher values give more verbosity).
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('verbosityLevel',verbosityLevel,1)
       call Galacticus_Verbosity_Level_Set(verbosityLevel)
       
       ! Open the Galacticus output file.
       call Galacticus_Output_Open_File

       call Galacticus_Display_Indent('Starting task set')
       doneStart=.true.
    end if
    Galacticus_Task_Start=.false.
    return
  end function Galacticus_Task_Start

  !# <galacticusTask>
  !#  <unitName>Galacticus_Task_End</unitName>
  !#  <after>Galacticus_Task_Start</after>
  !# </galacticusTask>
  logical function Galacticus_Task_End()
    use Galacticus_Output_Open
    use Galacticus_Output_Merger_Tree
    implicit none
    logical, save :: doneEnd=.false.

    if (.not.doneEnd) then

       ! Finalize outputs.
       call Galacticus_Merger_Tree_Output_Finalize()

       ! Close the Galacticus output file.
       call Galacticus_Output_Close_File

       call Galacticus_Display_Unindent('Finished task set')
       doneEnd=.true.
    end if
    Galacticus_Task_End=.false.
    return
  end function Galacticus_Task_End

end module Galacticus_Tasks_Basic
