!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations of dark matter halo bias.

module Dark_Matter_Halo_Biases
  !% Implements calculations of dark matter halo bias.
  use Tree_Nodes
  use ISO_Varying_String
  implicit none
  private
  public :: Dark_Matter_Halo_Bias

  ! Flag to indicate if this module has been initialized.  
  logical                                      :: haloBiasInitialized=.false.

  ! Name of halo bias method used.
  type(varying_string)                         :: darkMatterHaloBiasMethod

  ! Pointer to the function that returns halo bias.
  procedure(Halo_Bias_Template), pointer :: Dark_Matter_Halo_Bias_Get => null()
  abstract interface
     double precision function Halo_Bias_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Halo_Bias_Template
  end interface

contains

  subroutine Dark_Matter_Halo_Bias_Initialize
    !% Initalize the dark matter halo bias module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="darkMatterHaloBiasMethod" type="moduleUse">
    include 'structure_formation.CDM.halo_bias.modules.inc'
    !# </include>
    implicit none

    !$omp critical(haloBiasInitialize)
    if (.not.haloBiasInitialized) then
       ! Get the halo bias method parameter.
       !@ <inputParameter>
       !@   <name>darkMatterHaloBiasMethod</name>
       !@   <defaultValue>Tinker2010</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Selects which dark matter halo bias method to use.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('darkMatterHaloBiasMethod',darkMatterHaloBiasMethod,defaultValue='Tinker2010')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="darkMatterHaloBiasMethod" type="code" action="subroutine">
       !#  <subroutineArgs>darkMatterHaloBiasMethod,Dark_Matter_Halo_Bias_Get</subroutineArgs>
       include 'structure_formation.CDM.halo_bias.inc'
       !# </include>
       if (.not.associated(Dark_Matter_Halo_Bias_Get)) call&
            & Galacticus_Error_Report('Dark_Matter_Halo_Bias_Initialize','method '//char(darkMatterHaloBiasMethod)//' is unrecognized')
       ! Flag that the module is now initialized.
       haloBiasInitialized=.true.
    end if
    !$omp end critical(haloBiasInitialize)

    return
  end subroutine Dark_Matter_Halo_Bias_Initialize
  
  double precision function Dark_Matter_Halo_Bias(thisNode)
    !% Computes the bias for a dark matter halo.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Ensure the module is initalized.
    call Dark_Matter_Halo_Bias_Initialize

    ! Get the dark matter halo bias.
    Dark_Matter_Halo_Bias=Dark_Matter_Halo_Bias_Get(thisNode)

    return
  end function Dark_Matter_Halo_Bias

end module Dark_Matter_Halo_Biases
