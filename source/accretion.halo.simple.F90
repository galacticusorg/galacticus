!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations of baryonic accretion onto halos using a simple truncation to mimic
!% reionization.

module Accretion_Halos_Simple
  !% Implements calculations of baryonic accretion onto halos using a simple truncation to mimic reionization.
  private
  public :: Accretion_Halos_Simple_Initialize

  ! Parameters controlling when accretion is suppressed.
  double precision :: reionizationSuppressionRedshift,reionizationSuppressionTime,reionizationSuppressionVelocity

contains

  !# <accretionHalosMethod>
  !#  <unitName>Accretion_Halos_Simple_Initialize</unitName>
  !# </accretionHalosMethod>
  subroutine Accretion_Halos_Simple_Initialize(accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get &
       &,Halo_Baryonic_Accreted_Mass_Get,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    use Cosmology_Functions
    implicit none
    type(varying_string),          intent(in)    :: accretionHalosMethod
    procedure(),          pointer, intent(inout) :: Halo_Baryonic_Accretion_Rate_Get,Halo_Baryonic_Accreted_Mass_Get&
         &,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get
    
    if (accretionHalosMethod == 'simple') then
       Halo_Baryonic_Accretion_Rate_Get        => Halo_Baryonic_Accretion_Rate_Simple_Get
       Halo_Baryonic_Accreted_Mass_Get         => Halo_Baryonic_Accreted_Mass_Simple_Get
       Halo_Baryonic_Failed_Accretion_Rate_Get => Halo_Baryonic_Failed_Accretion_Rate_Simple_Get
       Halo_Baryonic_Failed_Accreted_Mass_Get  => Halo_Baryonic_Failed_Accreted_Mass_Simple_Get
       !@ <inputParameter>
       !@   <name>reionizationSuppressionRedshift</name>
       !@   <defaultValue>9.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The redshift below which baryonic accretion is suppressed.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("reionizationSuppressionRedshift",reionizationSuppressionRedshift,defaultValue= 9.0d0)
       reionizationSuppressionTime=Cosmology_Age(Expansion_Factor_from_Redshift(reionizationSuppressionRedshift))
       !@ <inputParameter>
       !@   <name>reionizationSuppressionVelocity</name>
       !@   <defaultValue>30.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The velocity scale below which baryonic accretion is suppressed.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("reionizationSuppressionVelocity",reionizationSuppressionVelocity,defaultValue=30.0d0)
    end if
    return
  end subroutine Accretion_Halos_Simple_Initialize

  double precision function Halo_Baryonic_Accretion_Rate_Simple_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: growthRate,unaccretedMass

    if (thisNode%isSatellite()) then
       Halo_Baryonic_Accretion_Rate_Simple_Get=0.0d0
    else
       if (Tree_Node_Time(thisNode) > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
            & reionizationSuppressionVelocity) then
          Halo_Baryonic_Accretion_Rate_Simple_Get=0.0d0
       else
          Halo_Baryonic_Accretion_Rate_Simple_Get=(Omega_b()/Omega_0())*Tree_Node_Mass_Accretion_Rate(thisNode)
          unaccretedMass=Tree_Node_Hot_Halo_Unaccreted_Mass(thisNode)
          if (unaccretedMass > 0.0d0) then
             growthRate=Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)
             Halo_Baryonic_Accretion_Rate_Simple_Get=Halo_Baryonic_Accretion_Rate_Simple_Get+unaccretedMass*growthRate
          end if
       end if
    end if
    return
  end function Halo_Baryonic_Accretion_Rate_Simple_Get

  double precision function Halo_Baryonic_Accreted_Mass_Simple_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (Tree_Node_Time(thisNode) > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
         & reionizationSuppressionVelocity) then
       Halo_Baryonic_Accreted_Mass_Simple_Get=0.0d0
    else
       Halo_Baryonic_Accreted_Mass_Simple_Get=(Omega_b()/Omega_0())*Tree_Node_Mass(thisNode)
    end if
    return
  end function Halo_Baryonic_Accreted_Mass_Simple_Get
  
  double precision function Halo_Baryonic_Failed_Accretion_Rate_Simple_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: growthRate,unaccretedMass

    if (thisNode%isSatellite()) then
       Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=0.0d0
    else
       if (Tree_Node_Time(thisNode) > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
            & reionizationSuppressionVelocity) then
          Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=(Omega_b()/Omega_0())*Tree_Node_Mass_Accretion_Rate(thisNode)
       else
          unaccretedMass=Tree_Node_Hot_Halo_Unaccreted_Mass(thisNode)
          if (unaccretedMass > 0.0d0) then
             growthRate=Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)
             Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=-unaccretedMass*growthRate
          else
             Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=0.0d0
          end if
       end if
    end if
    return
  end function Halo_Baryonic_Failed_Accretion_Rate_Simple_Get

  double precision function Halo_Baryonic_Failed_Accreted_Mass_Simple_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (Tree_Node_Time(thisNode) > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
         & reionizationSuppressionVelocity) then
       Halo_Baryonic_Failed_Accreted_Mass_Simple_Get=(Omega_b()/Omega_0())*Tree_Node_Mass(thisNode)
    else
       Halo_Baryonic_Failed_Accreted_Mass_Simple_Get=0.0d0
    end if
    return
  end function Halo_Baryonic_Failed_Accreted_Mass_Simple_Get
  
end module Accretion_Halos_Simple
