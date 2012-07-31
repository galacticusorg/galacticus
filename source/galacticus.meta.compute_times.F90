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


!% Contains a module which implements calculations of the time taken to process merger trees.

module Galacticus_Meta_Compute_Times
  !% Implements calculations of the time taken to process merger trees.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Galacticus_Time_Per_Tree
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: metaComputeTimesInitialized=.false.

  ! Name of tree timing method used.
  type(varying_string) :: timePerTreeMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Time_Per_Tree_Template), pointer :: Galacticus_Time_Per_Tree_Get => null()
  abstract interface
     double precision function Time_Per_Tree_Template(treeRootMass)
       double precision, intent(in) :: treeRootMass
     end function Time_Per_Tree_Template
  end interface

contains

  subroutine Galacticus_Time_Per_Tree_Initialize
    !% Initialize the time per tree module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="timePerTreeMethod" type="moduleUse">
    include 'galacticus.meta.compute_times.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Galacticus_Time_Per_Tree_Initialization) 
    ! Initialize if necessary.
    if (.not.metaComputeTimesInitialized) then
       ! Get the time per tree method parameter.
       !@ <inputParameter>
       !@   <name>timePerTreeMethod</name>
       !@   <defaultValue>file</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing the time per tree.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('timePerTreeMethod',timePerTreeMethod,defaultValue='file')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="timePerTreeMethod" type="code" action="subroutine">
       !#  <subroutineArgs>timePerTreeMethod,Galacticus_Time_Per_Tree_Get</subroutineArgs>
       include 'galacticus.meta.compute_times.inc'
       !# </include>
       if (.not.associated(Galacticus_Time_Per_Tree_Get)) call Galacticus_Error_Report('Galacticus_Time_Per_Tree_Initialize'&
            &,'method '//char(timePerTreeMethod)//' is unrecognized')
       metaComputeTimesInitialized=.true.
    end if
    !$omp end critical(Galacticus_Time_Per_Tree_Initialization) 
    
    return
  end subroutine Galacticus_Time_Per_Tree_Initialize

  double precision function Galacticus_Time_Per_Tree(treeRootMass)
    !% Returns the time (in seconds) to compute a tree of mass {\tt treeRootMass}.
    implicit none
    double precision, intent(in) :: treeRootMass

    ! Initialize the module.
    call Galacticus_Time_Per_Tree_Initialize

    ! Call the function to do the actual calculation.
    Galacticus_Time_Per_Tree=Galacticus_Time_Per_Tree_Get(treeRootMass)
    return
  end function Galacticus_Time_Per_Tree
  
end module Galacticus_Meta_Compute_Times
