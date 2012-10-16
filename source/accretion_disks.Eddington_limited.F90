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


!% Contains a module which implements calculations of properties of accretion disks which ignore the details of accretion physics
!% and assume that black hole jets have a power that is a fixed fraction of the Eddington luminosity.

module Accretion_Disks_Eddington
  !% Implements calculations of properties of accretion disks which ignore the details of accretion physics
  !% and assume that black hole jets have a power that is a fixed fraction of the Eddington luminosity.
  implicit none
  private
  public :: Accretion_Disks_Eddington_Initialize

  ! Efficiency parameters for the accretion disk model.
  double precision :: accretionDiskJetPowerEddington,accretionDiskRadiativeEfficiencyEddington

contains

  !# <accretionDisksMethod>
  !#  <unitName>Accretion_Disks_Eddington_Initialize</unitName>
  !# </accretionDisksMethod>
  subroutine Accretion_Disks_Eddington_Initialize(accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get&
       &,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: accretionDisksMethod
    procedure(double precision), pointer, intent(inout) :: Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get&
         &,Accretion_Disk_Jet_Power_Get
    character(len=30)                            :: accretionRateThin
    
    if (accretionDisksMethod == 'eddingtonLimited') then
       Accretion_Disk_Radiative_Efficiency_Get => Accretion_Disk_Radiative_Efficiency_Eddington
       Black_Hole_Spin_Up_Rate_Get             => Black_Hole_Spin_Up_Rate_Eddington
       Accretion_Disk_Jet_Power_Get            => Accretion_Disk_Jet_Power_Eddington
       !@ <inputParameter>
       !@   <name>accretionDiskJetPowerEddington</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The jet power produced by an Eddington limited accretion disk in units of the Eddington luminosity.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionDiskJetPowerEddington",accretionDiskJetPowerEddington,defaultValue=0.1d0)
       !@ <inputParameter>
       !@   <name>accretionDiskRadiativeEfficiencyEddington</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The radiative efficiency of an Eddington limited accretion disk.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionDiskRadiativeEfficiencyEddington",accretionDiskRadiativeEfficiencyEddington,defaultValue=0.1d0)
    end if
    return
  end subroutine Accretion_Disks_Eddington_Initialize

  double precision function Accretion_Disk_Radiative_Efficiency_Eddington(thisNode,massAccretionRate)
    !% Computes the radiative efficiency for an Eddington-limited accretion disk.
    use Black_Hole_Fundamentals
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    Accretion_Disk_Radiative_Efficiency_Eddington=accretionDiskRadiativeEfficiencyEddington
    return
  end function Accretion_Disk_Radiative_Efficiency_Eddington

  double precision function Accretion_Disk_Jet_Power_Eddington(thisNode,massAccretionRate)
    !% Computes the jet power for an Eddington-limited accretion disk.
    use Tree_Nodes
    use Black_Hole_Fundamentals
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    Accretion_Disk_Jet_Power_Eddington=accretionDiskJetPowerEddington*Black_Hole_Eddington_Accretion_Rate(thisNode)*(speedLight/kilo)**2
    return
  end function Accretion_Disk_Jet_Power_Eddington

  double precision function Black_Hole_Spin_Up_Rate_Eddington(thisNode,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisNode} due to accretion from an Eddington-limited accretion disk.
    !% disk. This is always zero, as no physical model is specified for this accretion disk method.
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    Black_Hole_Spin_Up_Rate_Eddington=0.0d0
    return
  end function Black_Hole_Spin_Up_Rate_Eddington

end module Accretion_Disks_Eddington
