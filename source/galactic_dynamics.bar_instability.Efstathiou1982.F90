!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !!{
  Implementation of the \cite{efstathiou_stability_1982} model for galactic disk bar instability.
  !!}

  !![
  <galacticDynamicsBarInstability name="galacticDynamicsBarInstabilityEfstathiou1982">
   <description>
    A galactic dynamics bar instability class that uses the stability criterion of \cite{efstathiou_stability_1982} to estimate
    when disks are unstable to bar formation:
    \begin{equation}
     \epsilon \left( \equiv {V_\mathrm{peak} \over \sqrt{\G M_\mathrm{disk}/r_\mathrm{disk}}} \right) &lt; \epsilon_\mathrm{c},
    \end{equation}
    for stability, where $V_\mathrm{peak}$ is the peak velocity in the rotation curve\footnote{In practice, the velocity is
    evaluated at the disk scale radius and multiplied by the factor $1.1800237580$ which relates this velocity to the peak
    rotation velocity for an isolated, thin exponential disk.}, $M_\mathrm{disk}$ is the mass of the disk and $r_\mathrm{disk}$ is
    its scale length (assuming an exponential disk). The value of $\epsilon_\mathrm{c}$ is linearly interpolated in the disk gas
    fraction between values for purely gaseous and stellar disks as specified by {\normalfont \ttfamily
    [stabilityThresholdStellar]} and {\normalfont \ttfamily [stabilityThresholdGaseous]} respectively. For disks which are judged
    to be unstable, the timescale for bar formation is estimated to be
    \begin{equation}
     t_\mathrm{bar} = t_\mathrm{disk} \left( {\epsilon_\mathrm{c} - \epsilon_\mathrm{iso} \over \epsilon_\mathrm{c} - \epsilon}
     \right)^2,
    \end{equation}
    where $\epsilon_\mathrm{iso}$ is the value of $\epsilon$ for an isolated disk and $t_\mathrm{disk}$ is the disk dynamical
    time, defined as $r/V$, at one scale length. This form gives an infinite timescale at the stability threshold, reducing to
    a dynamical time for highly unstable disks, while also ensuring that the slope of $t_\mathrm{bar}$ is continuous at the
    instability threshold. This method returns zero external driving torque.
   </description>
  </galacticDynamicsBarInstability>
  !!]
  type, extends(galacticDynamicsBarInstabilityClass) :: galacticDynamicsBarInstabilityEfstathiou1982
     !!{
     Implementation of the \cite{efstathiou_stability_1982} model for galactic disk bar instability.
     !!}
     private
     double precision :: stabilityThresholdStellar          , stabilityThresholdGaseous              , &
          &              fractionAngularMomentumRetainedDisk, fractionAngularMomentumRetainedSpheroid, &
          &              timescaleMinimum
   contains
     !![
     <methods>
       <method description="Compute the stability estimator for the \cite{efstathiou_stability_1982} model for galactic disk bar instability." method="estimator" />
     </methods>
     !!]
     procedure :: timescale => efstathiou1982Timescale
     procedure :: estimator => efstathiou1982Estimator
  end type galacticDynamicsBarInstabilityEfstathiou1982

  interface galacticDynamicsBarInstabilityEfstathiou1982
     !!{
     Constructors for the \refClass{galacticDynamicsBarInstabilityEfstathiou1982} model for galactic disk bar instability class.
     !!}
     module procedure efstathiou1982ConstructorParameters
     module procedure efstathiou1982ConstructorInternal
  end interface galacticDynamicsBarInstabilityEfstathiou1982

  double precision, parameter :: stabilityDiskIsolated=0.6221297315d0

contains

  function efstathiou1982ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticDynamicsBarInstabilityEfstathiou1982} model for galactic disk bar instability class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticDynamicsBarInstabilityEfstathiou1982)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: stabilityThresholdStellar          , stabilityThresholdGaseous              , &
         &                                                                           fractionAngularMomentumRetainedDisk, fractionAngularMomentumRetainedSpheroid, &
         &                                                                           timescaleMinimum

    !![
    <inputParameter>
      <name>stabilityThresholdStellar</name>
      <defaultValue>1.1d0</defaultValue>
      <description>The stability threshold in the \cite{efstathiou_stability_1982} algorithm for purely stellar disks.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>stabilityThresholdGaseous</name>
      <defaultValue>0.7d0</defaultValue>
      <description>The stability threshold in the \cite{efstathiou_stability_1982} algorithm for purely gaseous disks.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timescaleMinimum</name>
      <defaultValue>1.0d-9</defaultValue>
      <description>The minimum absolute dynamical timescale (in Gyr) to use in the \cite{efstathiou_stability_1982} algorithm.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fractionAngularMomentumRetainedDisk</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The fraction of angular momentum of material depleted from the disk by bar instability which is retained in the disk.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fractionAngularMomentumRetainedSpheroid</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The fraction of angular momentum of material depleted from the disk by bar instability which is retained in the spheroid.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=galacticDynamicsBarInstabilityEfstathiou1982(stabilityThresholdStellar,stabilityThresholdGaseous,timescaleMinimum,fractionAngularMomentumRetainedDisk,fractionAngularMomentumRetainedSpheroid)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function efstathiou1982ConstructorParameters

  function efstathiou1982ConstructorInternal(stabilityThresholdStellar,stabilityThresholdGaseous,timescaleMinimum,fractionAngularMomentumRetainedDisk,fractionAngularMomentumRetainedSpheroid) result(self)
    !!{
    Internal constructor for the \refClass{galacticDynamicsBarInstabilityEfstathiou1982} model for galactic disk bar instability class.
    !!}
    implicit none
    type            (galacticDynamicsBarInstabilityEfstathiou1982)                :: self
    double precision                                              , intent(in   ) :: stabilityThresholdStellar          , stabilityThresholdGaseous              , &
         &                                                                           fractionAngularMomentumRetainedDisk, fractionAngularMomentumRetainedSpheroid, &
         &                                                                           timescaleMinimum
    !![
    <constructorAssign variables="stabilityThresholdStellar, stabilityThresholdGaseous, timescaleMinimum, fractionAngularMomentumRetainedDisk, fractionAngularMomentumRetainedSpheroid"/>
    !!]

    return
  end function efstathiou1982ConstructorInternal

  subroutine efstathiou1982Timescale(self,node,timescale,externalDrivingSpecificTorque,fractionAngularMomentumRetainedDisk,fractionAngularMomentumRetainedSpheroid)
    !!{
    Computes a timescale for depletion of a disk to a pseudo-bulge via bar instability based on the criterion of
    \cite{efstathiou_stability_1982}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDisk, treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear         , megaParsec
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (galacticDynamicsBarInstabilityEfstathiou1982), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(  out) :: externalDrivingSpecificTorque                      , timescale                              , &
         &                                                                           fractionAngularMomentumRetainedDisk                , fractionAngularMomentumRetainedSpheroid
    class           (nodeComponentDisk                           ), pointer       :: disk
    ! Maximum timescale (in dynamical times) allowed.
    double precision                                              , parameter     :: timescaleDimensionlessMaximum      =1.0000000000d10
    double precision                                                              :: massDisk                                           , timeDynamical                          , &
         &                                                                           fractionGas                                        , stabilityEstimator                     , &
         &                                                                           stabilityEstimatorRelative                         , stabilityIsolatedRelative              , &
         &                                                                           stabilityThreshold                                 , timescaleDimensionless

    ! Assume infinite timescale (i.e. no instability) initially.
    timescale                    =-1.0d0
    externalDrivingSpecificTorque= 0.0d0
    ! Set the fraction of angular momentum retained in the disk and spheroid.
    fractionAngularMomentumRetainedDisk    =self%fractionAngularMomentumRetainedDisk
    fractionAngularMomentumRetainedSpheroid=self%fractionAngularMomentumRetainedSpheroid
    ! Get the disk.
    disk => node%disk()
    ! Compute the disk mass.
    massDisk=disk%massGas()+disk%massStellar()
    ! Return if disk has unphysical angular momentum.
    if (disk%angularMomentum() <= 0.0d0                            ) return
    ! Return if disk has unphysical velocity or radius.
    if (disk%velocity       () <= 0.0d0 .or. disk%radius() <= 0.0d0) return
    ! Compute the gas fraction in the disk.
    fractionGas=disk%massGas()/massDisk
    ! Compute the stability threshold.
    stabilityThreshold=+self%stabilityThresholdStellar*(1.0d0-fractionGas) &
         &             +self%stabilityThresholdGaseous*       fractionGas
    ! Compute the stability estimator for this node.
    stabilityEstimator=self%estimator(node)
    ! Check if the disk is bar unstable.
    if (stabilityEstimator < stabilityThreshold) then
       ! Disk is unstable, compute a timescale for depletion.
       ! Begin by finding the disk dynamical time.
       timeDynamical=(megaParsec/kilo/gigaYear)*disk%radius()/min(disk%velocity(),speedLight/kilo)
       ! Simple scaling which gives infinite timescale at the threshold, decreasing to dynamical time for a maximally unstable
       ! disk.
       stabilityIsolatedRelative =stabilityThreshold-stabilityDiskIsolated
       stabilityEstimatorRelative=stabilityThreshold-stabilityEstimator
       if (stabilityIsolatedRelative > timescaleDimensionlessMaximum*stabilityEstimatorRelative) then
          timescaleDimensionless=timescaleDimensionlessMaximum
       else
          timescaleDimensionless=(stabilityIsolatedRelative/stabilityEstimatorRelative)**2
       end if
       timescale=max(timeDynamical,self%timescaleMinimum)*timescaleDimensionless
    end if
    return
  end subroutine efstathiou1982Timescale

  double precision function efstathiou1982Estimator(self,node)
    !!{
    Compute the stability estimator for the \cite{efstathiou_stability_1982} model for galactic disk bar instability.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDisk             , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (galacticDynamicsBarInstabilityEfstathiou1982), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    ! Factor by which to boost velocity (evaluated at the scale radius) to convert to maximum velocity (assuming an isolated
    ! exponential disk) as appears in stability criterion.
    double precision                                              , parameter     :: velocityBoostFactor=1.1800237580d0
    class           (nodeComponentDisk                           ), pointer       :: disk
    double precision                                                              :: massDisk                          , velocitySelf
    !$GLC attributes unused :: self

    ! Assume an extremely stable disk by default.
    efstathiou1982Estimator=huge(0.0d0)
    ! Get the disk.
    disk => node%disk()
    ! Check for a physically-plausible disk.
    if (disk%radius() <= 0.0d0) return
    ! Compute the disk mass.
    massDisk=disk%massGas()+disk%massStellar()
    if (massDisk < 0.0d0) return
    ! Compute the velocity due to the disk's self-gravity.
    velocitySelf=+sqrt(                                &
         &             +gravitationalConstant_internal &
         &             *massDisk                       &
         &             /disk%radius  ()                &
         &            )  
    if     (                                                                                           &
         &                                                          velocitySelf  <=            0.0d0  &
         &  .or.                                                                                       &
         &   exponent(velocityBoostFactor*disk%velocity())-exponent(velocitySelf) > maxExponent(0.0d0) &
         & ) return
    ! Compute the stability estimator for this node.
    efstathiou1982Estimator=max(                            &
         &                      +stabilityDiskIsolated    , &
         &                      +     velocityBoostFactor   &
         &                      *disk%velocity          ()  &
         &                      /     velocitySelf          &
         &                     )
    return
  end function efstathiou1982Estimator

