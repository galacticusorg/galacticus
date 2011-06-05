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


!% Contains a module which implements calculations of dark matter halo mass accretion histories.

module Dark_Matter_Halo_Mass_Accretion_Histories
  !% Implements calculations of dark matter halo mass accretion histories.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Dark_Matter_Halo_Mass_Accretion_Time

  ! Flag to indicate if this module has been initialized.  
  logical              :: darkMatterAccretionHistoryInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: darkMatterAccretionHistoryMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Dark_Matter_Accretion_Template), pointer :: Dark_Matter_Halo_Mass_Accretion_Time_Get => null()
  abstract interface
     double precision function Dark_Matter_Accretion_Template(baseNode,nodeMass)
       import treeNode
       type(treeNode),   intent(inout), pointer :: baseNode
       double precision, intent(in)             :: nodeMass
     end function Dark_Matter_Accretion_Template
  end interface

contains

  subroutine Dark_Matter_Mass_Accretion_Initialize
    !% Initialize the dark matter mass accretion history module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="darkMatterAccretionHistoryMethod" type="moduleUse">
    include 'dark_matter_halos.mass_accretion_history.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Dark_Matter_Mass_Accretion_Initialization) 
    ! Initialize if necessary.
    if (.not.darkMatterAccretionHistoryInitialized) then
       ! Get the mass accretion history method parameter.
       !@ <inputParameter>
       !@   <name>darkMatterAccretionHistoryMethod</name>
       !@   <defaultValue>Wechsler 2002</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of dark matter halo mass accretion histories.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('darkMatterAccretionHistoryMethod',darkMatterAccretionHistoryMethod,defaultValue='Wechsler 2002')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="darkMatterAccretionHistoryMethod" type="code" action="subroutine">
       !#  <subroutineArgs>darkMatterAccretionHistoryMethod,Dark_Matter_Halo_Mass_Accretion_Time_Get</subroutineArgs>
       include 'dark_matter_halos.mass_accretion_history.inc'
       !# </include>
       if (.not.associated(Dark_Matter_Halo_Mass_Accretion_Time_Get)) &
            & call Galacticus_Error_Report('Dark_Matter_Mass_Accretion_Initialize','method ' //char(darkMatterAccretionHistoryMethod)//' is unrecognized')
       darkMatterAccretionHistoryInitialized=.true.
    end if
    !$omp end critical(Dark_Matter_Mass_Accretion_Initialization) 

    return
  end subroutine Dark_Matter_Mass_Accretion_Initialize

  double precision function Dark_Matter_Halo_Mass_Accretion_Time(baseNode,nodeMass)
    !% Returns the time for {\tt thisNode} in {\tt thisTree} according to the mass accretion history.
    implicit none
    type(treeNode),   pointer, intent(inout) :: baseNode
    double precision,          intent(in)    :: nodeMass

    ! Initialize the module.
    call Dark_Matter_Mass_Accretion_Initialize

    ! Get the time for the node.
    Dark_Matter_Halo_Mass_Accretion_Time=Dark_Matter_Halo_Mass_Accretion_Time_Get(baseNode,nodeMass)

    return
  end function Dark_Matter_Halo_Mass_Accretion_Time
  
end module Dark_Matter_Halo_Mass_Accretion_Histories
