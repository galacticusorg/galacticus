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


!% Contains a module which implements calculations of satellite merging times using the \cite{lacey_merger_1993} method.

module Dynamical_Friction_Lacey_Cole
  !% Implements calculations of satellite merging times using the \cite{lacey_merger_1993} method.
  private
  public :: Satellite_Time_Until_Merging_Lacey_Cole_Initialize

  double precision :: mergingTimescaleMultiplier

contains

  !# <satelliteMergingMethod>
  !#  <unitName>Satellite_Time_Until_Merging_Lacey_Cole_Initialize</unitName>
  !# </satelliteMergingMethod>
  subroutine Satellite_Time_Until_Merging_Lacey_Cole_Initialize(satelliteMergingMethod,Satellite_Time_Until_Merging)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string), intent(in)    :: satelliteMergingMethod
    procedure(double precision), pointer, intent(inout) :: Satellite_Time_Until_Merging

    if (satelliteMergingMethod == 'Lacey-Cole') then
       Satellite_Time_Until_Merging => Satellite_Time_Until_Merging_Lacey_Cole
       !@ <inputParameter>
       !@   <name>mergingTimescaleMultiplier</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     A multiplier for the merging timescale in the ``Lacey-Cole'' dynamical friction timescale method.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergingTimescaleMultiplier',mergingTimescaleMultiplier,defaultValue=1.0d0)
    end if
    return
  end subroutine Satellite_Time_Until_Merging_Lacey_Cole_Initialize

  double precision function Satellite_Time_Until_Merging_Lacey_Cole(thisNode)
    !% Return the timescale for merging satellites using the \cite{lacey_merger_1993} method.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Virial_Orbits
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    type(treeNode),   pointer                :: hostNode
    logical,          parameter              :: acceptUnboundOrbits=.false.
    double precision, parameter              :: inverseTwoB1=1.169335453d0 ! 1/2/B(1).
    double precision                         :: angularMomentum,orbitalEnergy,equivalentCircularOrbitRadius,orbitalCircularity &
         &,velocityScale,radialScale,massRatio

    ! Find the host node.
    hostNode => thisNode%parentNode
    ! Get orbital parameters for this satellite.
    call Virial_Orbital_Parameters(thisNode,acceptUnboundOrbits,angularMomentum=angularMomentum,orbitalEnergy=orbitalEnergy)
    ! Get velocity scale.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    radialScale=Dark_Matter_Halo_Virial_Radius(hostNode)
    ! Compute radius of orbit with same energy.
    equivalentCircularOrbitRadius=exp(orbitalEnergy/velocityScale**2+0.5d0)
    ! Compute orbital circularity.
    orbitalCircularity=angularMomentum/velocityScale/radialScale/equivalentCircularOrbitRadius
    ! Compute mass ratio.
    massRatio=Tree_Node_Mass(hostNode)/Tree_Node_Mass(thisNode)
    ! Compute dynamical friction timescale.
    Satellite_Time_Until_Merging_Lacey_Cole=(orbitalCircularity**0.78d0)*(equivalentCircularOrbitRadius**2)&
         &*Dark_Matter_Halo_Dynamical_Timescale(hostNode)*mergingTimescaleMultiplier*inverseTwoB1*massRatio/dlog(massRatio)
    return
  end function Satellite_Time_Until_Merging_Lacey_Cole

end module Dynamical_Friction_Lacey_Cole
