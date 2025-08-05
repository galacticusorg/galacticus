!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implementation of the \cite{efstathiou_stability_1982} model for galactic disk bar instability, but including the effects of
  tidal forces.
  !!}

  use :: Satellites_Tidal_Fields, only : satelliteTidalFieldClass

  !![
  <galacticDynamicsBarInstability name="galacticDynamicsBarInstabilityEfstathiou1982Tidal">
   <description>
    A galactic dynamics bar instability class that uses the stability criterion of \cite{efstathiou_stability_1982} to estimate
    when disks are unstable to bar formation, but includes an additional term due to external tidal forces:
    \begin{equation}
     \epsilon \left( \equiv {V_\mathrm{peak} \over \sqrt{\G M_\mathrm{disk}/r_\mathrm{disk}} + \hbox{max}(\mathcal{T},0) }
     \right) &lt; \epsilon_\mathrm{c},
    \end{equation}
    for stability, where $V_\mathrm{peak}$ is the peak velocity in the rotation curve (computed here assuming an isolated
    exponential disk), $M_\mathrm{disk}$ is the mass of the disk, $r_\mathrm{disk}$ is its scale length (assuming an
    exponential disk), $\mathcal{T} = \mathcal{F} r_\mathrm{disk}^2$ is the external driving specific torque, and $\mathcal{F}$
    is the external tidal field (evaluated using the selected method; see \refPhysics{satelliteTidalField}). The value of
    $\epsilon_\mathrm{c}$ is linearly interpolated in the disk gas fraction between values for purely gaseous and stellar disks
    as specified by {\normalfont \ttfamily stabilityThresholdStellar} and {\normalfont \ttfamily stabilityThresholdGaseous}
    respectively. For disks which are judged to be unstable, the timescale for bar formation is estimated to be
    \begin{equation}
     t_\mathrm{bar} = t_\mathrm{disk} {\epsilon_\mathrm{c} - \epsilon_\mathrm{iso} \over \epsilon_\mathrm{c} - \epsilon},
    \end{equation}
    where $\epsilon_\mathrm{iso}$ is the value of $\epsilon$ for an isolated disk and $t_\mathrm{disk}$ is the disk dynamical
    time, defined as $r/V$, at one scale length. This form gives an infinite timescale at the stability threshold, reducing to
    a dynamical time for highly unstable disks.
   </description>
  </galacticDynamicsBarInstability>
  !!]
  type, extends(galacticDynamicsBarInstabilityEfstathiou1982) :: galacticDynamicsBarInstabilityEfstathiou1982Tidal
     !!{
     Implementation of the \cite{efstathiou_stability_1982} model for galactic disk bar instability, but include the effects of tidal forces.
     !!}
     private
     class           (satelliteTidalFieldClass), pointer :: satelliteTidalField_     => null()
     double precision                                    :: massThresholdHarrassment
   contains
     !![
     <methods>
       <method description="Compute the radial term of the tidal tensor." method="tidalTensorRadial" />
     </methods>
     !!]
     final     ::                      efstathiou1982Destructor
     procedure :: tidalTensorRadial => efstathiou1982TidalTidalTensorRadial
     procedure :: timescale         => efstathiou1982TidalTimescale
     procedure :: estimator         => efstathiou1982TidalEstimator
  end type galacticDynamicsBarInstabilityEfstathiou1982Tidal

  interface galacticDynamicsBarInstabilityEfstathiou1982Tidal
     !!{
     Constructors for the \refClass{galacticDynamicsBarInstabilityEfstathiou1982Tidal} model for galactic disk bar instability class.
     !!}
     module procedure efstathiou1982TidalConstructorParameters
     module procedure efstathiou1982TidalConstructorInternal
  end interface galacticDynamicsBarInstabilityEfstathiou1982Tidal

contains

  function efstathiou1982TidalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticDynamicsBarInstabilityEfstathiou1982Tidal} model for galactic disk bar instability class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(galacticDynamicsBarInstabilityEfstathiou1982Tidal)                :: self
    type(inputParameters                                  ), intent(inout) :: parameters

    !![
    <inputParameter>
      <name>massThresholdHarrassment</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The host halo mass threshold for harassment to take effect.</description>
      <source>parameters</source>
    <variable>self%massThresholdHarrassment</variable>
    </inputParameter>
    <objectBuilder class="satelliteTidalField" name="self%satelliteTidalField_" source="parameters"/>
    !!]
    self%galacticDynamicsBarInstabilityEfstathiou1982=galacticDynamicsBarInstabilityEfstathiou1982(parameters)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function efstathiou1982TidalConstructorParameters

  function efstathiou1982TidalConstructorInternal(stabilityThresholdStellar,stabilityThresholdGaseous,timescaleMinimum,fractionAngularMomentumRetainedDisk_,fractionAngularMomentumRetainedSpheroid_,massThresholdHarrassment,satelliteTidalField_) result(self)
    !!{
    Internal constructor for the \refClass{galacticDynamicsBarInstabilityEfstathiou1982Tidal} model for galactic disk bar instability class.
    !!}
    implicit none
    type            (galacticDynamicsBarInstabilityEfstathiou1982Tidal)                        :: self
    double precision                                                   , intent(in   )         :: stabilityThresholdStellar, stabilityThresholdGaseous               , &
         &                                                                                        timescaleMinimum         , fractionAngularMomentumRetainedDisk_    , &
         &                                                                                        massThresholdHarrassment , fractionAngularMomentumRetainedSpheroid_
    class           (satelliteTidalFieldClass                         ), intent(in   ), target :: satelliteTidalField_
    !![
    <constructorAssign variables="massThresholdHarrassment, *satelliteTidalField_"/>
    !!]

    self%galacticDynamicsBarInstabilityEfstathiou1982=galacticDynamicsBarInstabilityEfstathiou1982(stabilityThresholdStellar,stabilityThresholdGaseous,timescaleMinimum,fractionAngularMomentumRetainedDisk_,fractionAngularMomentumRetainedSpheroid_)
    return
  end function efstathiou1982TidalConstructorInternal

  subroutine efstathiou1982Destructor(self)
    !!{
    Destructor for the \refClass{galacticDynamicsBarInstabilityEfstathiou1982Tidal} model for galactic disk bar instability class.
    !!}
    implicit none
    type(galacticDynamicsBarInstabilityEfstathiou1982Tidal), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalField_"/>
    !!]
    return
  end subroutine efstathiou1982Destructor

  subroutine efstathiou1982TidalTimescale(self,node,timescale,externalDrivingSpecificTorque,fractionAngularMomentumRetainedDisk,fractionAngularMomentumRetainedSpheroid)
    !!{
    Computes a timescale for depletion of a disk to a pseudo-bulge via bar instability based on the criterion of
    \cite{efstathiou_stability_1982}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    class           (galacticDynamicsBarInstabilityEfstathiou1982Tidal), intent(inout) :: self
    type            (treeNode                                         ), intent(inout) :: node
    double precision                                                   , intent(  out) :: externalDrivingSpecificTorque      , timescale                              , &
         &                                                                                fractionAngularMomentumRetainedDisk, fractionAngularMomentumRetainedSpheroid
    class           (nodeComponentSpheroid                            ), pointer       :: spheroid

    call self%galacticDynamicsBarInstabilityEfstathiou1982%timescale(node,timescale,externalDrivingSpecificTorque,fractionAngularMomentumRetainedDisk,fractionAngularMomentumRetainedSpheroid)
    ! Compute the external torque.
    if (timescale > 0.0d0) then
       spheroid                      =>  node    %spheroid         (    )
       externalDrivingSpecificTorque =  +self    %tidalTensorRadial(node)    &
            &                           *spheroid%radius           (    )**2
    end if
    return
  end subroutine efstathiou1982TidalTimescale

  double precision function efstathiou1982TidalEstimator(self,node)
    !!{
    Compute the stability estimator for the \cite{efstathiou_stability_1982} model for galactic disk bar instability.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDisk             , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (galacticDynamicsBarInstabilityEfstathiou1982Tidal), intent(inout) :: self
    type            (treeNode                                         ), intent(inout) :: node
    ! Factor by which to boost velocity (evaluated at scale radius) to convert to maximum velocity (assuming an isolated disk) as
    ! appears in stability criterion.
    double precision                                                   , parameter     :: velocityBoostFactor=1.1800237580d0
    class           (nodeComponentDisk                                ), pointer       :: disk
    double precision                                                                   :: massDisk

    ! Get the disk.
    disk => node%disk()
    ! Compute the disk mass.
    massDisk=disk%massGas()+disk%massStellar()
    ! Return perfect stability if there is no disk.
    if (massDisk <= 0.0d0) then
       efstathiou1982TidalEstimator=huge(0.0d0)
    else
       ! Compute the stability estimator for this node.
       efstathiou1982TidalEstimator=max(                                             &
            &                           +stabilityDiskIsolated                     , &
            &                           +velocityBoostFactor                         &
            &                           *           disk%velocity         (    )     &
            &                           /sqrt(                                       &
            &                                 +gravitationalConstant_internal        &
            &                                 *massDisk                              &
            &                                 /     disk%radius           (    )     &
            &                                 +max(                                  &
            &                                      +self%tidalTensorRadial(node)     &
            &                                      *disk%radius           (    )**2, &
            &                                      +0.0d0                            &
            &                                     )                                  &
            &                                )                                       &
            &                          )
    end if
    return
  end function efstathiou1982TidalEstimator

  double precision function efstathiou1982TidalTidalTensorRadial(self,node)
    !!{
    Compute the radial part of the tidal tensor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(galacticDynamicsBarInstabilityEfstathiou1982Tidal), intent(inout) :: self
    type (treeNode                                         ), intent(inout) :: node
    class(nodeComponentBasic                               ), pointer       :: basicParent

    basicParent => node%parent%basic()
    if (node%isSatellite() .and. basicParent%mass() > self%massThresholdHarrassment) then
       efstathiou1982TidalTidalTensorRadial=self%satelliteTidalField_%tidalTensorRadial(node)
    else
       efstathiou1982TidalTidalTensorRadial=0.0d0
    end if
    return
  end function efstathiou1982TidalTidalTensorRadial
