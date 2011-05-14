!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module of halo formation time methods.

module Tree_Node_Methods_Formation_Times_Cole2000
  !% Implement tracking of halo formation times.
  use Tree_Nodes
  use Components
  private
  public :: Tree_Node_Methods_Formation_Times_Cole2000_Initialize, Tree_Node_Methods_Formation_Times_Cole2000_Dump
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=0, dataCount=2, historyCount=0
  integer, parameter :: formationTimeIndex=1
  integer, parameter :: formationMassIndex=2

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Formation_Time</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical          :: methodSelected=.false.

  ! Factor by which mass must increase to trigger a new formation event.
  double precision :: haloReformationMassFactor

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Formation_Times_Cole2000_Initialize</unitName>
  !#  <optionName>treeNodeMethodFormationTimes</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Formation_Times_Cole2000_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node formation time tracking module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'Cole2000') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex    =componentTypeCount

       ! Display message.
       message='Cole2000 halo formation time method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Formation_Time              => Tree_Node_Formation_Time_Formation_Times_Cole2000
       Tree_Node_Formation_Time_Set          => null()
       Tree_Node_Formation_Time_Rate_Adjust  => null()
       Tree_Node_Formation_Time_Rate_Compute => Tree_Node_Formation_Time_Rate_Compute_Cole2000

       !@ <inputParameter>
       !@   <name>haloReformationMassFactor</name>
       !@   <defaultValue>2.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Factor by which halo mass must have increased to trigger a new formation event.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('haloReformationMassFactor',haloReformationMassFactor,defaultValue=2.0d0)
    end if
    return
  end subroutine Tree_Node_Methods_Formation_Times_Cole2000_Initialize

  double precision function Tree_Node_Formation_Time_Formation_Times_Cole2000(thisNode)
    !% Return the time of the last major merger.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
    
    if (.not.thisNode%componentExists(componentIndex)) call Tree_Node_Formation_Time_Create_Component(thisNode)
    thisIndex=Tree_Node_Formation_Times_Cole2000_Index(thisNode)
    Tree_Node_Formation_Time_Formation_Times_Cole2000=thisNode%components(thisIndex)%data(formationTimeIndex)
    return
  end function Tree_Node_Formation_Time_Formation_Times_Cole2000

  subroutine Tree_Node_Formation_Time_Rate_Compute_Cole2000(thisNode,interrupt,interruptProcedure)
    !% Compute the exponential disk node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    integer                                :: thisIndex
   
    ! If no formation time component exists, stop and make one.
    if (.not.thisNode%componentExists(componentIndex)) then
       interrupt=.true.
       interruptProcedure => Tree_Node_Formation_Time_Create_Component
       return
    end if

    ! Check if the halo has grown sufficiently in mass to trigger a new formation event.
    thisIndex=Tree_Node_Formation_Times_Cole2000_Index(thisNode)
    if (Tree_Node_Mass(thisNode) > haloReformationMassFactor*thisNode%components(thisIndex)%data(formationMassIndex)) then
       interrupt=.true.
       interruptProcedure => Tree_Node_Formation_Time_Create_Component
       return
    end if

    return
  end subroutine Tree_Node_Formation_Time_Rate_Compute_Cole2000

  subroutine Tree_Node_Formation_Time_Create_Component(thisNode)
    !% Creates a halo formation time component for {\tt thisNode}. This function is also used to ``reform'' the halo, since it
    !% simply resets the formation time and mass to the current values.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    use Events_Halo_Formation
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message
    integer                                     :: thisIndex

    ! Display a message.
    if (.not.thisNode%componentExists(componentIndex)) then
       message='Creating halo formation component for node '
       message=message//thisNode%index()
       call Galacticus_Display_Message(message,verbosityInfo)
    end if
    ! Get the index of the component (which will also ensure that the component is created).
    thisIndex=Tree_Node_Formation_Times_Cole2000_Index(thisNode)
    ! Set initial formation time and formation mass.
    thisNode%components(thisIndex)%data(formationTimeIndex)=Tree_Node_Time(thisNode)
    thisNode%components(thisIndex)%data(formationMassIndex)=Tree_Node_Mass(thisNode)
    ! Trigger a halo formation event.
    call Event_Halo_Formation(thisNode)
    return
  end subroutine Tree_Node_Formation_Time_Create_Component

  integer function Tree_Node_Formation_Times_Cole2000_Index(thisNode)
    !% Ensure the merging statistics component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Formation_Times_Cole2000_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Formation_Times_Cole2000_Index

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Formation_Times_Cole2000_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Formation_Times_Cole2000_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
 
    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'halo formation time component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'formation time:',Tree_Node_Formation_Time_Formation_Times_Cole2000(thisNode)
       else
          write (0,'(1x,a)'           ) 'halo formation time component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Formation_Times_Cole2000_Dump

end module Tree_Node_Methods_Formation_Times_Cole2000
