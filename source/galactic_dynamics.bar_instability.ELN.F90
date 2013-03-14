!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of bar instability based on the \cite{efstathiou_stability_1982} criterion.

module Galactic_Dynamics_Bar_Instabilities_ELN
  !% Implements calculations of bar instability based on the \cite{efstathiou_stability_1982} criterion.
  use Galacticus_Nodes
  implicit none
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
    type(varying_string),                 intent(in)    :: barInstabilityMethod
    procedure(double precision), pointer, intent(inout) :: Bar_Instability_Timescale_Get
    
    if (barInstabilityMethod == 'ELN') then
       Bar_Instability_Timescale_Get => Bar_Instability_Timescale_ELN
       ! Read in stability threshold parameters.
       !@ <inputParameter>
       !@   <name>stabilityThresholdStellar</name>
       !@   <defaultValue>0.7</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The stability threshold in the \cite{efstathiou_stability_1982} algorithm for purely stellar disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stabilityThresholdStellar',stabilityThresholdStellar,defaultValue=0.7d0)
       !@ <inputParameter>
       !@   <name>stabilityThresholdGaseous</name>
       !@   <defaultValue>0.7</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The stability threshold in the \cite{efstathiou_stability_1982} algorithm for purely gaseous disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stabilityThresholdGaseous',stabilityThresholdGaseous,defaultValue=0.7d0)
    end if
  
    return
  end subroutine Galactic_Dynamics_Bar_Instabilities_ELN_Initialize

  double precision function Bar_Instability_Timescale_ELN(thisNode)
    !% Computes a timescale for depletion of a disk to a pseudo-bulge via bar instability based on the criterion of
    !% \cite{efstathiou_stability_1982}.
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    implicit none
    type (treeNode         ), intent(inout), pointer :: thisNode
    class(nodeComponentDisk),                pointer :: thisDiskComponent
    double precision        , parameter              :: stabilityIsolatedDisk=0.6221297315d0
    ! Factor by which to boost velocity (evaluated at scale radius) to convert to maximum velocity (assuming an isolated disk) as
    ! appears in stability criterion.
    double precision        , parameter              :: velocityBoostFactor          =1.180023758d0
    ! Maximum timescale (in dynamical times) allowed.
    double precision        , parameter              :: timescaleDimensionlessMaximum=1.0d10
    double precision                                 :: stabilityEstimator,stabilityThreshold,dynamicalTime,gasFraction,diskMass &
         &,timescaleDimensionless,stabilityIsolatedRelative,stabilityEstimatorRelative

    ! Assume infinite timescale (i.e. no instability) initially.
    Bar_Instability_Timescale_ELN=-1.0d0

    ! Get the disk.
    thisDiskComponent => thisNode%disk()

    ! Compute the disk mass.
    diskMass=thisDiskComponent%massGas()+thisDiskComponent%massStellar()
    ! Return if there is no disk.
    if (diskMass <= 0.0d0) return

    ! Return if disk has unphysical angular momentum.
    if (thisDiskComponent%angularMomentum() <= 0.0d0) return

    ! Return if disk has unphysical velocity or radius.
    if (thisDiskComponent%velocity() <= 0.0d0 .or. thisDiskComponent%radius() <= 0.0d0) return

    ! Compute the gas fraction in the disk.
    gasFraction=thisDiskComponent%massGas()/diskMass

    ! Compute the stability threshold.
    stabilityThreshold=stabilityThresholdStellar*(1.0d0-gasFraction)+stabilityThresholdGaseous*gasFraction

    ! Compute the stability estimator for this node.
    stabilityEstimator=max(&
         &                  stabilityIsolatedDisk               , &
         &                  velocityBoostFactor                   &
         &                 *thisDiskComponent      %velocity()    &
         &                 /sqrt(                                 &
         &                        gravitationalConstantGalacticus &
         &                       *diskMass                        &
         &                       /thisDiskComponent%radius  ()    &
         &                      )                                 &
         &                )
    
    ! Check if the disk is bar unstable.
    if (stabilityEstimator < stabilityThreshold) then
       ! Disk is unstable, compute a timescale for depletion.
       
       ! Begin by finding the disk dynamical time.
       dynamicalTime=(megaParsec/kilo/gigaYear)*thisDiskComponent%radius()/thisDiskComponent%velocity()

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
