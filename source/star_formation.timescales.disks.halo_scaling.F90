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


!% Contains a module which implements a star formation timescale for galactic disks which scales with halo virial velocity and
!% redshift.

module Star_Formation_Timescale_Disks_Halo_Scaling
  !% Implements a star formation timescale for galactic disks which scales with halo virial velocity and
  !% redshift.
  implicit none
  private
  public :: Star_Formation_Timescale_Disks_Halo_Scaling_Initialize

  ! Parameters of the timescale model.
  double precision :: starFormationTimescaleDisksHaloScalingTimescale&
       &,starFormationTimescaleDisksHaloScalingVirialVelocityExponent,starFormationTimescaleDisksHaloScalingRedshiftExponent
  
contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Halo_Scaling_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Halo_Scaling_Initialize(starFormationTimescaleDisksMethod,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``halo scaling'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Tree_Nodes
    implicit none
    type(varying_string),                 intent(in)    :: starFormationTimescaleDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Timescale_Disk_Get
    
    if (starFormationTimescaleDisksMethod == 'haloScaling') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Halo_Scaling
       ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksHaloScalingTimescale</name>
       !@   <defaultValue>1 Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The timescale for star formation in the halo scaling timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksHaloScalingTimescale',starFormationTimescaleDisksHaloScalingTimescale,defaultValue=1.0d0)
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksHaloScalingVirialVelocityExponent</name>
       !@   <defaultValue>$0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of virial velocity in the timescale for star formation in the halo scaling timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksHaloScalingVirialVelocityExponent',starFormationTimescaleDisksHaloScalingVirialVelocityExponent,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksHaloScalingRedshiftExponent</name>
       !@   <defaultValue>$0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of redshift in the timescale for star formation in the halo scaling timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksHaloScalingRedshiftExponent',starFormationTimescaleDisksHaloScalingRedshiftExponent,defaultValue=0.0d0)
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Halo_Scaling_Initialize

  double precision function Star_Formation_Timescale_Disk_Halo_Scaling(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\tt thisNode} in the halo scaling timescale model.
    use Tree_Nodes
    use Cosmology_Functions
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode)  , intent(inout), pointer :: thisNode
    double precision, parameter              :: virialVelocityNormalization=200.0d0
    double precision                         :: expansionFactor,virialVelocity

    ! Get virial velocity and expansion factor.
    virialVelocity =Dark_Matter_Halo_Virial_Velocity               (thisNode)
    expansionFactor=Expansion_Factor                (Tree_Node_Time(thisNode))

    ! Return the timescale.
    Star_Formation_Timescale_Disk_Halo_Scaling=                                                                        &
         &  starFormationTimescaleDisksHaloScalingTimescale                                                            &
         & *(virialVelocity/virialVelocityNormalization)**starFormationTimescaleDisksHaloScalingVirialVelocityExponent &
         & /expansionFactor**starFormationTimescaleDisksHaloScalingRedshiftExponent
    return
  end function Star_Formation_Timescale_Disk_Halo_Scaling
  
end module Star_Formation_Timescale_Disks_Halo_Scaling
