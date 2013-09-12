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

!% Contains a module which implements calculations of bar instability based on the \cite{efstathiou_stability_1982} criterion, but
!% including the effects of tidal forces.

module Galactic_Dynamics_Bar_Instabilities_ELN_Tidal
  !% Implements calculations of bar instability based on the \cite{efstathiou_stability_1982} criterion, but including the effects
  !% of tidal forces.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galactic_Dynamics_Bar_Instabilities_ELN_Tidal_Initialize

  ! Stability parameters for stellar and gaseous disks.
  double precision :: stabilityThresholdGaseous, stabilityThresholdStellar

  ! Harrassment mass threshold.
  double precision :: harrassmentMassThreshold

contains

  !# <barInstabilityMethod>
  !#  <unitName>Galactic_Dynamics_Bar_Instabilities_ELN_Tidal_Initialize</unitName>
  !# </barInstabilityMethod>
  subroutine Galactic_Dynamics_Bar_Instabilities_ELN_Tidal_Initialize(barInstabilityMethod,Bar_Instability_Timescale_Get)
    !% Initializes the ``ELN+tidal'' bar instability module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                                   ), intent(in   )          :: barInstabilityMethod
    procedure(Bar_Instability_Timescale_ELN_Tidal              ), intent(inout), pointer :: Bar_Instability_Timescale_Get

    if (barInstabilityMethod == 'ELN+tidal') then
       Bar_Instability_Timescale_Get => Bar_Instability_Timescale_ELN_Tidal
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
       !@ <inputParameter>
       !@   <name>harrassmentMassThreshold</name>
       !@   <defaultValue>0.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The host halo mass threshold for harrassment to take effect.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('harrassmentMassThreshold',harrassmentMassThreshold,defaultValue=0.0d0)
    end if

    return
  end subroutine Galactic_Dynamics_Bar_Instabilities_ELN_Tidal_Initialize

  subroutine Bar_Instability_Timescale_ELN_Tidal(thisNode,barInstabilityTimeScale,barInstabilityExternalDrivingSpecificTorque)
    !% Computes a timescale for depletion of a disk to a pseudo-bulge via bar instability based on the criterion of
    !% \cite{efstathiou_stability_1982}, but including an additional term due to external tidal forces.
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Satellites_Tidal_Fields
    implicit none
    type            (treeNode             )           , intent(inout), pointer :: thisNode
    double precision                                  , intent(  out)          :: barInstabilityExternalDrivingSpecificTorque               , barInstabilityTimeScale
    class           (nodeComponentDisk    )                          , pointer :: thisDisk
    class           (nodeComponentSpheroid)                          , pointer :: thisSpheroid
    class           (nodeComponentBasic   )                          , pointer :: parentBasic
    double precision                       , parameter                         :: stabilityIsolatedDisk                      =0.6221297315d0
    ! Factor by which to boost velocity (evaluated at scale radius) to convert to maximum velocity (assuming an isolated disk) as
    ! appears in stability criterion.
    double precision                       , parameter                         :: velocityBoostFactor                        =1.180023758d0
    ! Maximum timescale (in dynamical times) allowed.
    double precision                       , parameter                         :: timescaleDimensionlessMaximum              =1.0d10
    double precision                                                           :: diskMass                                                  , dynamicalTime            , &
         &                                                                        gasFraction                                               , stabilityEstimator       , &
         &                                                                        stabilityEstimatorRelative                                , stabilityIsolatedRelative, &
         &                                                                        stabilityThreshold                                        , tidalField               , &
         &                                                                        timescaleDimensionless

    ! Assume infinite timescale (i.e. no instability) initially.
    barInstabilityTimeScale                    =-1.0d0
    barInstabilityExternalDrivingSpecificTorque= 0.0d0
    ! Get the disk and spheroid.
    thisDisk     => thisNode       %disk    ()
    thisSpheroid => thisNode       %spheroid()
    parentBasic  => thisNode%parent%basic   ()
    ! Compute the disk mass.
    diskMass=thisDisk%massGas()+thisDisk%massStellar()
    ! Return if there is no disk.
    if (diskMass                            <= 0.0d0                       ) return
    ! Return if disk has unphysical angular momentum.
    if (thisDisk%angularMomentum() <= 0.0d0                                ) return
    ! Return if disk has unphysical velocity or radius.
    if (thisDisk%velocity       () <= 0.0d0 .or. thisDisk%radius() <= 0.0d0) return
    ! Find the tidal field.
    if (thisNode%isSatellite() .and. parentBasic%mass() > harrassmentMassThreshold) then
       tidalField=Satellite_Tidal_Field(thisNode)
    else
       tidalField=0.0d0
    end if
    ! Compute the gas fraction in the disk.
    gasFraction=thisDisk%massGas()/diskMass
    ! Compute the stability threshold.
    stabilityThreshold=stabilityThresholdStellar*(1.0d0-gasFraction)+stabilityThresholdGaseous*gasFraction
    ! Compute the stability estimator for this node.
    stabilityEstimator=max(                                        &
         &                  stabilityIsolatedDisk                , &
         &                  velocityBoostFactor                    &
         &                 *thisDisk%velocity()                    &
         &                 /sqrt(                                  &
         &                        gravitationalConstantGalacticus  &
         &                       *diskMass                         &
         &                       /thisDisk%radius()                &
         &                       +max(                             &
         &                             tidalField                  &
         &                            *thisDisk%radius()**2      , &
         &                             0.0d0                       &
         &                           )                             &
         &                      )                                  &
         &                )
    ! Check if the disk is bar unstable.
    if (stabilityEstimator < stabilityThreshold) then
       ! Disk is unstable, compute a timescale for depletion.
       ! Begin by finding the disk dynamical time.
       dynamicalTime=(megaParsec/kilo/gigaYear)*thisDisk%radius()/thisDisk%velocity()
       ! Simple scaling which gives infinite timescale at the threshold, decreasing to dynamical time for a maximally unstable
       ! disk.
       stabilityIsolatedRelative =stabilityThreshold-stabilityIsolatedDisk
       stabilityEstimatorRelative=stabilityThreshold-stabilityEstimator
       if (stabilityIsolatedRelative > timescaleDimensionlessMaximum*stabilityEstimatorRelative) then
          timescaleDimensionless=timescaleDimensionlessMaximum
       else
          timescaleDimensionless=stabilityIsolatedRelative/stabilityEstimatorRelative
       end if
       barInstabilityTimeScale=dynamicalTime*timescaleDimensionless
       ! Compute the external torque.
       barInstabilityExternalDrivingSpecificTorque=tidalField*thisSpheroid%radius()**2
    end if
    return
  end subroutine Bar_Instability_Timescale_ELN_Tidal

end module Galactic_Dynamics_Bar_Instabilities_ELN_Tidal
