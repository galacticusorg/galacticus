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


!% Contains a module which implements a cut off in the cooling rate at given redshift and virial velocity.

module Cooling_Rates_Modifier_Cut_Off
  !% Implements a cut off in the cooling rate at given redshift and virial velocity.
  implicit none
  private
  public :: Cooling_Rate_Modifier_Cut_Off

  ! Record of whether module has been initialized.
  logical          :: moduleInitialized=.false.

  ! Parameters controlling the time and velocity scale at which cooling is cut off.
  double precision :: coolingCutOffVelocity,coolingCutOffRedshift,coolingCutOffTime
  logical          :: coolingCutOffFormationNode

contains

  !# <coolingRateModifierMethod>
  !#  <unitName>Cooling_Rate_Modifier_Cut_Off</unitName>
  !# </coolingRateModifierMethod>
  subroutine Cooling_Rate_Modifier_Cut_Off(thisNode,coolingRate)
    !% Modify cooling rates by truncating them to zero below a given redshift and virial velocity.
    use Input_Parameters
    use Cosmology_Functions
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode)  , intent(inout), pointer :: thisNode
    double precision, intent(inout)          :: coolingRate
    double precision                         :: virialVelocity
    
    !$omp critical (Cooling_Rate_Modifier_Cut_Off_Initialize)
    if (.not.moduleInitialized) then
       !@ <inputParameter>
       !@   <name>coolingCutOffFormationNode</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether to use the virial velocity of the formation node or current node in the cooling rate ``cut-off'' modifier.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("coolingCutOffFormationNode",coolingCutOffFormationNode,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>coolingCutOffVelocity</name>
       !@   <defaultValue>0.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The velocity below which cooling is suppressed in the ``cut-off'' cooling rate modifier method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("coolingCutOffVelocity",coolingCutOffVelocity,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>coolingCutOffRedshift</name>
       !@   <defaultValue>0.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The redshift below which cooling is suppressed in the ``cut-off'' cooling rate modifier method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("coolingCutOffRedshift",coolingCutOffRedshift,defaultValue=0.0d0)
       coolingCutOffTime=Cosmology_Age(Expansion_Factor_from_Redshift(coolingCutOffRedshift))
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Cooling_Rate_Modifier_Cut_Off_Initialize)

    ! Test for halos where cooling should be cut off.
    select case (coolingCutOffFormationNode)
    case (.false.)
       virialVelocity=Dark_Matter_Halo_Virial_Velocity(thisNode              )
    case (.true. )
       virialVelocity=Dark_Matter_Halo_Virial_Velocity(thisNode%formationNode)
    end select
    if     (                                                  &
         &  Tree_Node_Time(thisNode) >= coolingCutOffTime     &
         &   .and.                                            &
         &  virialVelocity           <= coolingCutOffVelocity &
         & ) coolingRate=0.0d0
    return
  end subroutine Cooling_Rate_Modifier_Cut_Off
  
end module Cooling_Rates_Modifier_Cut_Off
