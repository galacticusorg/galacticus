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


!% Contains a module which handles computation and output of star formation histories for galaxies.

module Galacticus_Output_Star_Formation_Histories
  !% Handles computation and output of star formation histories for galaxies.
  use ISO_Varying_String
  use Tree_Nodes
  use Abundances_Structure
  use Histories
  use Kind_Numbers
  implicit none
  private
  public :: Star_Formation_History_Create, Star_Formation_History_Scales, Star_Formation_History_Record,&
       & Star_Formation_History_Output

  ! Flag indicating whether this module has been initialized.
  logical              :: starFormationHistoriesInitialized=.false.

  ! Name of the method to use.
  type(varying_string) :: starFormationHistoriesMethod

  ! Pointer to the subroutine that creates any history required for star formation histories.
  procedure(Star_Formation_History_Create_Template), pointer :: Star_Formation_History_Create_Do => null()
  abstract interface
     subroutine Star_Formation_History_Create_Template(thisNode,thisHistory)
       import treeNode, history
       type(treeNode), intent(inout), pointer :: thisNode
       type(history),  intent(inout)          :: thisHistory
     end subroutine Star_Formation_History_Create_Template
  end interface

  ! Pointer to the subroutine that sets scale factors for error control of star formation histories
  procedure(Star_Formation_History_Scales_Template), pointer :: Star_Formation_History_Scales_Do => null()
  abstract interface
     subroutine Star_Formation_History_Scales_Template(thisHistory,stellarMass,stellarAbundances)
       import abundancesStructure, history
       double precision,          intent(in)                  :: stellarMass
       type(abundancesStructure), intent(in)                  :: stellarAbundances
       type(history),             intent(inout)               :: thisHistory
     end subroutine Star_Formation_History_Scales_Template
  end interface

  ! Pointer to the subroutine that records the star formation history.
  procedure(Star_Formation_History_Record_Template), pointer :: Star_Formation_History_Record_Do => null()
  abstract interface
     subroutine Star_Formation_History_Record_Template(thisNode,thisHistory,fuelAbundances,starFormationRate)
       import treeNode, history, abundancesStructure
       type(treeNode),            intent(inout), pointer :: thisNode
       type(history),             intent(inout)          :: thisHistory
       type(abundancesStructure), intent(in)             :: fuelAbundances
       double precision,          intent(in)             :: starFormationRate
     end subroutine Star_Formation_History_Record_Template
  end interface

  ! Pointer to the subroutine that outputs the star formation history.
  procedure(Star_Formation_History_Output_Template), pointer :: Star_Formation_History_Output_Do => null()
  abstract interface
     subroutine Star_Formation_History_Output_Template(thisNode,nodePassesFilter,thisHistory,iOutput,treeIndex,componentLabel)
       import treeNode, history, abundancesStructure, kind_int8
       type(treeNode),          intent(inout), pointer :: thisNode
       logical,                 intent(in)             :: nodePassesFilter
       type(history),           intent(in)             :: thisHistory
       integer,                 intent(in)             :: iOutput
       integer(kind=kind_int8), intent(in)             :: treeIndex
       character(len=*),        intent(in)             :: componentLabel
     end subroutine Star_Formation_History_Output_Template
  end interface

contains

  subroutine Galacticus_Output_Star_Formation_Histories_Initialize
    !% Initialize the star formation histories module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationHistoriesMethod" type="moduleUse">
    include 'galacticus.output.merger_tree.star_formation.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Galacticus_Output_Star_Formation_Histories_Initialization) 
    ! Initialize if necessary.
    if (.not.starFormationHistoriesInitialized) then
       ! Get the star formation history method parameter.
       !@ <inputParameter>
       !@   <name>starFormationHistoriesMethod</name>
       !@   <defaultValue>null</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The method to use for computing and outputting star formation histories.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoriesMethod',starFormationHistoriesMethod,defaultValue='null')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="starFormationHistoriesMethod" type="code" action="subroutine">
       !#  <subroutineArgs>starFormationHistoriesMethod,Star_Formation_History_Create_Do,Star_Formation_History_Scales_Do,Star_Formation_History_Record_Do,Star_Formation_History_Output_Do</subroutineArgs>
       include 'galacticus.output.merger_tree.star_formation.inc'
       !# </include>
       if (.not.(associated(Star_Formation_History_Create_Do).and.associated(Star_Formation_History_Scales_Do).and.associated(Star_Formation_History_Record_Do).and.associated(Star_Formation_History_Output_Do))) &
            & call Galacticus_Error_Report('Galacticus_Output_Star_Formation_Histories_Initialize'&
            &,'method '//char(starFormationHistoriesMethod)//' is unrecognized')
       starFormationHistoriesInitialized=.true.
    end if
    !$omp end critical(Galacticus_Output_Star_Formation_Histories_Initialization) 

    return
  end subroutine Galacticus_Output_Star_Formation_Histories_Initialize

  subroutine Star_Formation_History_Create(thisNode,thisHistory)
    !% Create any history required for storing the star formation history.
    use Histories
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(history),  intent(inout)          :: thisHistory
  
    ! Ensure module is initialized.
    call Galacticus_Output_Star_Formation_Histories_Initialize

    ! Simply call the function which does the actual work.
    call Star_Formation_History_Create_Do(thisNode,thisHistory)

    return
  end subroutine Star_Formation_History_Create

  subroutine Star_Formation_History_Record(thisNode,thisHistory,fuelAbundances,starFormationRate)
    !% Record the star formation history for {\tt thisNode}.
    use Histories
    use Tree_Nodes
    implicit none
    type(treeNode),            intent(inout), pointer :: thisNode
    type(history),             intent(inout)          :: thisHistory
    type(abundancesStructure), intent(in)             :: fuelAbundances
    double precision,          intent(in)             :: starFormationRate
  
    ! Ensure module is initialized.
    call Galacticus_Output_Star_Formation_Histories_Initialize

    ! Simply call the function which does the actual work.
    call Star_Formation_History_Record_Do(thisNode,thisHistory,fuelAbundances,starFormationRate)

    return
  end subroutine Star_Formation_History_Record

  subroutine Star_Formation_History_Output(thisNode,nodePassesFilter,thisHistory,iOutput,treeIndex,componentLabel)
    !% Output the star formation history for {\tt thisNode}.
    use Histories
    use Tree_Nodes
    implicit none
    type(treeNode),          intent(inout), pointer :: thisNode
    logical,                 intent(in)             :: nodePassesFilter
    type(history),           intent(inout)          :: thisHistory
    integer,                 intent(in)             :: iOutput
    integer(kind=kind_int8), intent(in)             :: treeIndex
    character(len=*),        intent(in)             :: componentLabel
  
    ! Ensure module is initialized.
    call Galacticus_Output_Star_Formation_Histories_Initialize

    ! Simply call the function which does the actual work.
    call Star_Formation_History_Output_Do(thisNode,nodePassesFilter,thisHistory,iOutput,treeIndex,componentLabel)

    return
  end subroutine Star_Formation_History_Output

  subroutine Star_Formation_History_Scales(thisHistory,stellarMass,stellarAbundances)
    !% Set the scaling factors for error control on the absolute value of stellar population properties.
    implicit none
    double precision,          intent(in)    :: stellarMass
    type(abundancesStructure), intent(in)    :: stellarAbundances
    type(history),             intent(inout) :: thisHistory
    
    ! Ensure module is initialized.
    call Galacticus_Output_Star_Formation_Histories_Initialize

    ! Simply call the subroutine which does the actual work.
    call Star_Formation_History_Scales_Do(thisHistory,stellarMass,stellarAbundances)
    return
  end subroutine Star_Formation_History_Scales

end module Galacticus_Output_Star_Formation_Histories
