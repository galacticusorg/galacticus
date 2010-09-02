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


!% Contains a module which implements calculations of bar instability based on the \cite{efstathiou_stability_1982} criterion.

module Galactic_Dynamics_Bar_Instabilities_ELN
  !% Implements calculations of bar instability based on the \cite{efstathiou_stability_1982} criterion.
  use Tree_Nodes
  private
  public :: Galactic_Dynamics_Bar_Instabilities_ELN_Initialize

  ! Stability parameters for stellar and gaseous disks.
  double precision :: stabilityThresholdStellar,stabilityThresholdGaseous

contains

  !# <barInstabilityMethod>
  !#  <unitName>Galactic_Dynamics_Bar_Instabilities_ELN_Initialize</unitName>
  !# </barInstabilityMethod>
  subroutine Galactic_Dynamics_Bar_Instabilities_ELN_Initialize(barInstabilityMethod,Bar_Instability_Timescale_Get)
    !% Initializes the ``ELN'' bar instability module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: barInstabilityMethod
    procedure(),          pointer, intent(inout) :: Bar_Instability_Timescale_Get
    
    if (barInstabilityMethod == 'ELN') then
       Bar_Instability_Timescale_Get => Bar_Instability_Timescale_ELN
       ! Read in stability threshold parameters.
       !@ <inputParameter>
       !@   <name>stabilityThresholdStellar</name>
       !@   <defaultValue>1.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The stability threshold in the \cite{efstathiou_stability_1982} algorithm for purely stellar disks.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('stabilityThresholdStellar',stabilityThresholdStellar,defaultValue=1.1d0)
       !@ <inputParameter>
       !@   <name>stabilityThresholdGaseous</name>
       !@   <defaultValue>0.9</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The stability threshold in the \cite{efstathiou_stability_1982} algorithm for purely gaseous disks.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('stabilityThresholdGaseous',stabilityThresholdGaseous,defaultValue=0.9d0)
    end if
  
    return
  end subroutine Galactic_Dynamics_Bar_Instabilities_ELN_Initialize

  double precision function Bar_Instability_Timescale_ELN(thisNode)
    !% Computes a timescale for depletion of a disk to a pseudo-bulge via bar instability based on the criterion of
    !% \cite{efstathiou_stability_1982}.
    use Tree_Nodes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, parameter              :: stabilityIsolatedDisk=0.6221297315d0
    ! Factor by which to boost velocity (evaluated at scale radius) to convert to maximum velocity (assuming an isolated disk) as
    ! appears in stability criterion.
    double precision, parameter              :: velocityBoostFactor          =1.180023758d0
    ! Maximum timescale (in dynamical times) allowed.
    double precision, parameter              :: timescaleDimensionlessMaximum=1.0d10
    double precision                         :: stabilityEstimator,stabilityThreshold,dynamicalTime,gasFraction,diskMass&
         &,timescaleDimensionless,stabilityIsolatedRelative,stabilityEstimatorRelative

    ! Assume infinite timescale (i.e. no instability) initially.
    Bar_Instability_Timescale_ELN=-1.0d0

    ! Compute the disk mass.
    diskMass=Tree_Node_Disk_Gas_Mass(thisNode)+Tree_Node_Disk_Stellar_Mass(thisNode)
    ! Return if there is no disk.
    if (diskMass <= 0.0d0) return

    ! Return if disk has unphysical angular momentum.
    if (Tree_Node_Disk_Angular_Momentum(thisNode) <= 0.0d0) return

    ! Return if disk has unphysical velocity or radius.
    if (Tree_Node_Disk_Velocity(thisNode) <= 0.0d0 .or. Tree_Node_Disk_Radius(thisNode) <= 0.0d0) return

    ! Compute the gas fraction in the disk.
    gasFraction=Tree_Node_Disk_Gas_Mass(thisNode)/diskMass

    ! Compute the stability threshold.
    stabilityThreshold=stabilityThresholdStellar*(1.0d0-gasFraction)+stabilityThresholdGaseous*gasFraction

    ! Compute the stability estimator for this node.
    stabilityEstimator=velocityBoostFactor*Tree_Node_Disk_Velocity(thisNode)/dsqrt(gravitationalConstantGalacticus*diskMass&
         &/Tree_Node_Disk_Radius(thisNode))
    
    ! Check if the disk is bar unstable.
    if (stabilityEstimator < stabilityThreshold) then
       ! Disk is unstable, compute a timescale for depletion.
       
       ! Begin by finding the disk dynamical time.
       dynamicalTime=(megaParsec/kilo/gigaYear)*Tree_Node_Disk_Radius(thisNode)/Tree_Node_Disk_Velocity(thisNode)

       ! Simple scaling which gives infinite timescale at the threshold, decreasing to dynamical time for a maximally unstable
       ! disk.
       stabilityIsolatedRelative =stabilityThreshold-stabilityIsolatedDisk
       stabilityEstimatorRelative=stabilityThreshold-stabilityEstimator
       if (stabilityIsolatedRelative > timescaleDimensionlessMaximum*stabilityEstimatorRelative) then
          timescaleDimensionless=timescaleDimensionlessMaximum
       else
          timescaleDimensionless=stabilityIsolatedRelative/stabilityEstimatorRelative
       end if
       Bar_Instability_Timescale_ELN=dynamicalTime*timescaleDimensionless

    end if

    return
  end function Bar_Instability_Timescale_ELN
  
end module Galactic_Dynamics_Bar_Instabilities_ELN
