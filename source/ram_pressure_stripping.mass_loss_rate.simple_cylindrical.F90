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
  Implementation of a simple ram pressure stripping class for cylindrically symmetric systems.
  !!}

  use :: Hot_Halo_Ram_Pressure_Forces, only : hotHaloRamPressureForceClass

  !![
  <ramPressureStripping name="ramPressureStrippingSimpleCylindrical">
   <description>
    A ram pressure stripping class which applies to systems with cylindrical symmetry (e.g. disks), and computes the mass loss
    rate to be:
    \begin{equation}
     \dot{M}_\mathrm{gas, disk} = \hbox{min}\left({\beta \mathcal{F}_\mathrm{hot, host} \over 2 \pi \mathrm{G}
     \Sigma_\mathrm{gas}(r_\mathrm{half}) \Sigma_\mathrm{total}(r_\mathrm{half})}, R_\mathrm{maximum}\right) {M_\mathrm{gas,
     disk} \over \tau_\mathrm{dyn, disk}},
    \end{equation}
    where $\mathcal{F}_\mathrm{hot, host}$ is the ram pressure force due to the hot halo of the node's host (computed using the
    selected hot halo ram pressure force method; see \refPhysics{hotHaloRamPressureForce}), $\Sigma_\mathrm{gas}(r)$ is the gas
    surface density in the disk, $\Sigma_\mathrm{total}(r)$ is the total surface density in the disk, $r_\mathrm{half}$ is the
    disk half-mass radius, $M_\mathrm{gas, disk}$ is the total gas mass in the disk, $\tau_\mathrm{dyn, disk} =
    r_\mathrm{disk}/v_\mathrm{disk}$ is the dynamical time in the disk, $\beta=${\normalfont \ttfamily [beta]} scales the rate of
    mass loss, and $R_\mathrm{maximum}=${\normalfont \ttfamily [rateFractionalMaximum]} controls the maximum allowed rate of mass
    loss.
   </description>
  </ramPressureStripping>
  !!]
  type, extends(ramPressureStrippingClass) :: ramPressureStrippingSimpleCylindrical
     !!{
     Implementation of a simple model of ram pressure stripping of cylindrically symmetric systems.
     !!}
     private
     class           (hotHaloRamPressureForceClass), pointer :: hotHaloRamPressureForce_ => null()
     double precision                                        :: rateFractionalMaximum             , beta
   contains
     final     ::                 simpleCylindricalDestructor
     procedure :: rateMassLoss => simpleCylindricalRateMassLoss
  end type ramPressureStrippingSimpleCylindrical

  interface ramPressureStrippingSimpleCylindrical
     !!{
     Constructors for the \refClass{ramPressureStrippingSimpleCylindrical} model of ram pressure stripping of .
     !!}
     module procedure simpleCylindricalConstructorParameters
     module procedure simpleCylindricalConstructorInternal
  end interface ramPressureStrippingSimpleCylindrical

contains

  function simpleCylindricalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{ramPressureStrippingSimpleCylindrical} timescale for star formation feedback in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (ramPressureStrippingSimpleCylindrical)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (hotHaloRamPressureForceClass         ), pointer       :: hotHaloRamPressureForce_
    double precision                                                       :: rateFractionalMaximum   , beta

    !![
    <inputParameter>
      <name>rateFractionalMaximum</name>
      <defaultValue>10.0d0</defaultValue>
      <description>The maximum fractional mass loss rate per dynamical time in the simple model of mass loss in cylindrically symmetric systems due to ram pressure stripping.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scaling factor which multiplies the ram pressure mass loss rate.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="hotHaloRamPressureForce" name="hotHaloRamPressureForce_" source="parameters"/>
    !!]
    self=ramPressureStrippingSimpleCylindrical(rateFractionalMaximum,beta,hotHaloRamPressureForce_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloRamPressureForce_"/>
    !!]
    return
  end function simpleCylindricalConstructorParameters

  function simpleCylindricalConstructorInternal(rateFractionalMaximum,beta,hotHaloRamPressureForce_) result(self)
    !!{
    Internal constructor for the \refClass{ramPressureStrippingSimpleCylindrical} model of ram pressure stripping class.
    !!}
    implicit none
    type            (ramPressureStrippingSimpleCylindrical)                        :: self
    double precision                                       , intent(in   )         :: rateFractionalMaximum   , beta
    class           (hotHaloRamPressureForceClass         ), intent(in   ), target :: hotHaloRamPressureForce_
    !![
    <constructorAssign variables="rateFractionalMaximum, beta, *hotHaloRamPressureForce_"/>
    !!]

    return
  end function simpleCylindricalConstructorInternal

  subroutine simpleCylindricalDestructor(self)
    !!{
    Destructor for the \refClass{ramPressureStrippingSimpleCylindrical} model of ram pressure stripping class.
    !!}
    implicit none
    type(ramPressureStrippingSimpleCylindrical), intent(inout) :: self

    !![
    <objectDestructor name="self%hotHaloRamPressureForce_"/>
    !!]
    return
  end subroutine simpleCylindricalDestructor

  double precision function simpleCylindricalRateMassLoss(self,component)
    !!{
    Computes the mass loss rate from cylindrically-symmetric systems due to ram pressure stripping assuming a simple model. Specifically, the mass loss
    rate is
    \begin{equation}
    \dot{M}_\mathrm{gas} = -\alpha M_\mathrm{gas}/\tau,
    \end{equation}
    where
    \begin{equation}
    \alpha = \beta F_\mathrm{ram}/F_\mathrm{gravity},
    \end{equation}
    $F_\mathrm{ram}$ is the ram pressure force from the hot halo (see \refPhysics{hotHaloRamPressureForce}), and
    \begin{equation}
    F_\mathrm{gravity} = 2 \pi \mathrm{G} \Sigma_\mathrm{gas}(r_{1/2}) \Sigma_\mathrm{total}(r_{1/2})
    \end{equation}
    is the gravitational restoring force at the half-mass radius, $r_\mathrm{1/2}$.
    !!}
    use :: Coordinates                     , only : coordinateCylindrical, assignment(=)
    use :: Display                         , only : displayGreen         , displayBlue                   , displayMagenta, displayReset
    use :: Galactic_Structure_Options      , only : componentTypeDisk    , enumerationComponentTypeType  , massTypeAll   , massTypeGaseous
    use :: Galacticus_Nodes                , only : nodeComponentDisk    , treeNode
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear             , gravitationalConstant_internal, megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (ramPressureStrippingSimpleCylindrical), intent(inout) :: self
    class           (nodeComponent                        ), intent(inout) :: component
    type            (treeNode                             ), pointer       :: node
    class           (massDistributionClass                ), pointer       :: massDistributionGaseous, massDistributionTotal
    type            (enumerationComponentTypeType         )                :: componentType
    type            (coordinateCylindrical                )                :: coordinates
    double precision                                                       :: forceGravitational     , forceRamPressure   , &
         &                                                                    rateMassLossFractional , radiusHalfMass     , &
         &                                                                    surfaceDensityGas      , surfaceDensityTotal, &
         &                                                                    timeDynamical          , radius             , &
         &                                                                    radiusHalfMass         , velocity           , &
         &                                                                    massGas

    ! Assume no mass loss rate due to ram pressure by default.
    simpleCylindricalRateMassLoss=0.0d0
    ! Return immediately if this is not a satellite.
    node => component%hostNode
    if (.not.node%isSatellite()) return
    ! Get the ram pressure force due to the hot halo.
    forceRamPressure=self%hotHaloRamPressureForce_%force(node)
    ! Get required properties.
    select type (component)
    class is (nodeComponentDisk)
       componentType =componentTypeDisk
       radius        =component%radius        ()
       radiusHalfMass=component%halfMassRadius()
       velocity      =component%velocity      ()
       massGas       =component%massGas       ()
    class default
       componentType =componentTypeDisk
       radius        =0.0d0
       radiusHalfMass=0.0d0
       velocity      =0.
       massGas       =0.0d0
       call Error_Report(                                                                                                                                                                                                                                        &
            &            'only "'//displayBlue()//'disk'//displayReset()//'" components are supported by the "'//displayGreen()//'simpleCylindrical'//displayReset()//'" '//displayBlue()//'ramPressureStripping'//displayReset()//' class'//char(10)//          &
            &            displayGreen()//'HELP:'//displayReset()//' see '//displayMagenta()//'https://github.com/galacticusorg/galacticus/wiki/Troubleshooting:-Component-not-supported-by-ramPressureStripping-class'//displayReset()//{introspection:location} &
            &           )
    end select
    ! Compute the surface densities at the half mass radius.
    coordinates             =  [radiusHalfMass,0.0d0,0.0d0]
    massDistributionGaseous => node                   %massDistribution(massType=massTypeGaseous,componentType=componentType)
    massDistributionTotal   => node                   %massDistribution(massType=massTypeAll    ,componentType=componentType)
    surfaceDensityGas       =  massDistributionGaseous%surfaceDensity  (         coordinates                                )
    surfaceDensityTotal     =  massDistributionTotal  %surfaceDensity  (         coordinates                                )
    !![
    <objectDestructor name="massDistributionGaseous"/>
    <objectDestructor name="massDistributionTotal"  />
    !!]
    ! Compute the gravitational restoring force in the midplane.
    forceGravitational  =  +2.0d0                          &
         &                 *Pi                             &
         &                 *gravitationalConstant_internal &
         &                 *surfaceDensityGas              &
         &                 *surfaceDensityTotal
    ! Return zero rate if the gravitational force is zero.
    if (forceGravitational <= 0.0d0) return
    ! Compute the mass loss fraction per dynamical time.
    if (self%beta*forceRamPressure < self%rateFractionalMaximum*forceGravitational) then
       rateMassLossFractional=+self%beta          &
            &                 *forceRamPressure   &
            &                 /forceGravitational
    else
       rateMassLossFractional=self%rateFractionalMaximum
    end if
    ! Compute the dynamical time.
    timeDynamical                =+megaParsec             &
         &                        /kilo                   &
         &                        /gigaYear               &
         &                        *radius                 &
         &                        /velocity
    ! Compute the mass loss rate.
    simpleCylindricalRateMassLoss=+rateMassLossFractional &
         &                        *massGas                &
         &                        /timeDynamical
    return
  end function simpleCylindricalRateMassLoss
