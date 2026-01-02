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
  Implementation of a simple tidal stripping class.
  !!}

  use :: Satellites_Tidal_Fields, only : satelliteTidalFieldClass

  !![
  <tidalStripping name="tidalStrippingSimple">
    <description>
      A simple model of tidal stripping.  Specifically, the mass loss rate is
      \begin{equation}
      \dot{M} = -\alpha M/\tau,
      \end{equation}
      where
      \begin{equation}
      \alpha = \beta F_\mathrm{tidal}/F_\mathrm{gravity},
      \end{equation}
      $F_\mathrm{tidal}=\mathcal{F}_\mathrm{tidal} r_{1/2}$, $\mathcal{F}_\mathrm{tidal}$ is the tidal field from the host halo
      (see \refPhysics{satelliteTidalField}),      
      \begin{equation}
      F_\mathrm{gravity} = V_{1/2}^2(r_{1/2})/r_{1/2}
      \end{equation}      
      is the gravitational restoring force at the half-mass radius, $r_\mathrm{1/2}$, and $\tau =
      r_\mathrm{s}/v_\mathrm{c}(r_\mathrm{s})$ is the dynamical time of the galactic component with $r_\mathrm{s}$ being the scale
      radius of the component, and $v_\mathrm{c}(r)$ the circular velocity of the component at radius $r$.
    </description>
  </tidalStripping>
  !!]
  type, extends(tidalStrippingClass) :: tidalStrippingSimple
     !!{
     Implementation of a simple model of tidal stripping.
     !!}
     private
     class           (satelliteTidalFieldClass), pointer :: satelliteTidalField_  => null()
     double precision                                    :: rateFractionalMaximum          , beta
   contains
     final     ::                 simpleDestructor
     procedure :: rateMassLoss => simpleRateMassLoss
  end type tidalStrippingSimple

  interface tidalStrippingSimple
     !!{
     Constructors for the \refClass{tidalStrippingSimple} model of tidal stripping class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface tidalStrippingSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{tidalStrippingSimple} tidal stripping class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (tidalStrippingSimple    )                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (satelliteTidalFieldClass), pointer       :: satelliteTidalField_
    double precision                                          :: rateFractionalMaximum, beta

    !![
    <inputParameter>
      <name>rateFractionalMaximum</name>
      <defaultValue>10.0d0</defaultValue>
      <description>The maximum fractional mass loss rate per dynamical time in the simple model of mass loss due to tidal stripping.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scaling factor which multiplies the tidal mass loss rate.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="satelliteTidalField" name="satelliteTidalField_" source="parameters"/>
    !!]
    self=tidalStrippingSimple(rateFractionalMaximum,beta,satelliteTidalField_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalField_"/>
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(rateFractionalMaximum,beta,satelliteTidalField_) result(self)
    !!{
    Internal constructor for the \refClass{tidalStrippingSimple} model of tidal stripping class.
    !!}
    implicit none
    type            (tidalStrippingSimple    )                        :: self
    double precision                          , intent(in   )         :: rateFractionalMaximum, beta
    class           (satelliteTidalFieldClass), intent(in   ), target :: satelliteTidalField_
    !![
    <constructorAssign variables="rateFractionalMaximum, beta, *satelliteTidalField_"/>
    !!]

    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{
    Destructor for the \refClass{tidalStrippingSimple} model of tidal stripping class.
    !!}
    implicit none
    type(tidalStrippingSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalField_"/>
    !!]
    return
  end subroutine simpleDestructor

  double precision function simpleRateMassLoss(self,component)
    !!{
    Computes the mass loss rate due to tidal stripping assuming a simple model.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDisk    , nodeComponentSpheroid, treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear             , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Mass_Distributions              , only : massDistributionClass
    implicit none
    class           (tidalStrippingSimple ), intent(inout) :: self
    class           (nodeComponent        ), intent(inout) :: component
    type            (treeNode             ), pointer       :: node
    class           (massDistributionClass), pointer       :: massDistribution_
    double precision                                       :: forceGravitational    , forceTidal    , &
         &                                                    rateMassLossFractional, radiusHalfMass, &
         &                                                    tidalTensorRadial     , timeDynamical , &
         &                                                    velocityRotation      , velocity      , &
         &                                                    massGas               , massStellar   , &
         &                                                    radius

    ! Assume no mass loss rate due to tidal by default.
    simpleRateMassLoss=0.0d0
    ! Return immediately if this is not a satellite.
    node => component%hostNode
    if (.not.node%isSatellite()) return
    ! Get component properties.
    select type (component)
    class is (nodeComponentDisk    )
       radius        =component%radius        ()
       radiusHalfMass=component%halfMassRadius()
       velocity      =component%velocity      ()
       massGas       =component%massGas       ()
       massStellar   =component%massStellar   ()
    class is (nodeComponentSpheroid)
       radius        =component%radius        ()
       radiusHalfMass=component%halfMassRadius()
       velocity      =component%velocity      ()
       massGas       =component%massGas       ()
       massStellar   =component%massStellar   ()
    class default
       radius        =0.0d0
       radiusHalfMass=0.0d0
       velocity      =0.0d0
       massGas       =0.0d0
       massStellar   =0.0d0
       call Error_Report('unsupported component'//{introspection:location})
    end select
    ! Get the tidal field due to the host halo.
    tidalTensorRadial =  self%satelliteTidalField_%tidalTensorRadial(node)
    ! Get the tidal force exerted at the half-mass radius.
    forceTidal        =  +tidalTensorRadial &
         &               *radiusHalfMass
    ! Return if the tidal field is compressive.
    if (forceTidal <= 0.0d0) return
    ! Compute the rotation curve.
    massDistribution_ => node             %massDistribution(              )
    velocityRotation  =  massDistribution_%rotationCurve   (radiusHalfMass)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Compute the gravitational restoring force at the half-mass radius.
    forceGravitational=+velocityRotation**2 &
         &             /radiusHalfMass
    ! Return zero rate if the gravitational force is zero.
    if (forceGravitational <= 0.0d0) return
    ! Compute the mass loss fraction per dynamical time.
    if (self%beta*forceTidal < self%rateFractionalMaximum*forceGravitational) then
       rateMassLossFractional=+self%beta          &
            &                 *forceTidal         &
            &                 /forceGravitational
    else
       rateMassLossFractional=+self%rateFractionalMaximum
    end if
    ! Compute the dynamical time.
    timeDynamical     =+megaParsec             &
         &             /kilo                   &
         &             /gigaYear               &
         &             *radius                 &
         &             /velocity
    ! Compute the mass loss rate.
    simpleRateMassLoss=+rateMassLossFractional &
         &             *(                      &
         &               +massGas              &
         &               +massStellar          &
         &              )                      &
         &             /timeDynamical
    return
  end function simpleRateMassLoss
