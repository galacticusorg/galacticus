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


!% Contains a module which implements calculations of sizes of galactic components (or more general components).

module Galactic_Structure_Radii
  !% Implements calculations of sizes of galactic components (or more general components).
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Galactic_Structure_Radii_Solve

  ! Flag to indicate if this module has been initialized.  
  logical              :: galacticStructureRadiusSolverInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: galacticStructureRadiusSolverMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Galactic_Structure_Radii_Solve_Template), pointer :: Galactic_Structure_Radii_Solve_Do => null()
  abstract interface
     subroutine Galactic_Structure_Radii_Solve_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end subroutine Galactic_Structure_Radii_Solve_Template
  end interface
  
contains

  !# <preDerivativeTask>
  !# <unitName>Galactic_Structure_Radii_Solve</unitName>
  !# </preDerivativeTask>
  subroutine Galactic_Structure_Radii_Solve(thisNode)
    !% Solve for the radii of galactic components in {\tt thisNode}.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="galacticStructureRadiusSolverMethod" type="moduleUse">
    include 'galactic_structure.radius_solver.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    !$omp critical(Galactic_Structure_Radii_Initialization) 
    ! Initialize if necessary.
    if (.not.galacticStructureRadiusSolverInitialized) then
       ! Get the galactic structure radii solver method parameter.
       !@ <inputParameter>
       !@   <name>galacticStructureRadiusSolverMethod</name>
       !@   <defaultValue>adiabatic</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Selects the method to be used for solving for galactic structure.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('galacticStructureRadiusSolverMethod',galacticStructureRadiusSolverMethod,defaultValue='adiabatic')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="galacticStructureRadiusSolverMethod" type="code" action="subroutine">
       !#  <subroutineArgs>galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do</subroutineArgs>
       include 'galactic_structure.radius_solver.inc'
       !# </include>
       if (.not.associated(Galactic_Structure_Radii_Solve_Do)) &
            & call Galacticus_Error_Report('Galactic_Structure_Radii','method '//char(galacticStructureRadiusSolverMethod)//' is unrecognized')
       galacticStructureRadiusSolverInitialized=.true.
    end if
    !$omp end critical(Galactic_Structure_Radii_Initialization) 

    ! Call the routine to solve for the radii.
    !! <gfortran 4.6> We shouldn't need the OpenMP Critical section here, but it seems that the current version of gfortran 4.6
    !! has problems with the procedure pointers used in radius solving (Radius_Set, Radius_Get, Velocity_Set, Velocity_Get) which can
    !! cause them to become set by another thread, leading to strange results.
    !$omp critical(Galactic_Structure_Radii_Do)
    call Galactic_Structure_Radii_Solve_Do(thisNode)
    !$omp end critical(Galactic_Structure_Radii_Do)

    return
  end subroutine Galactic_Structure_Radii_Solve

end module Galactic_Structure_Radii
