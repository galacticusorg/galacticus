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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implementation of a tidal stripping class for galactic components which follows the model of :cite:t:`zentner_physics_2005`.
  !!}

  use :: Dark_Matter_Halo_Scales        , only : darkMatterHaloScaleClass
  use :: Satellite_Tidal_Stripping_Radii, only : satelliteTidalStrippingRadiusClass

  !![
  <tidalStripping name="tidalStrippingZentner2005" docformat="rst">
    <description>
    A model of tidal stripping of galactic components which applies the formalism of :cite:t:`zentner_physics_2005` (see
    :galacticus-class:`satelliteTidalStrippingZentner2005`) to the component, rather than to the dark matter halo. The mass loss
    rate is

    .. math::

       \dot{M} = -\alpha \frac{M(&gt;r_\mathrm{tidal})}{T_\mathrm{loss}},

    where :math:`\alpha=`\ ``[efficiency]``, :math:`M(&gt;r_\mathrm{tidal})` is the mass (gas plus stars) of the component outside of
    the tidal radius of the satellite, :math:`r_\mathrm{tidal}` (see :galacticus-class:`satelliteTidalStrippingRadius`), and
    :math:`T_\mathrm{loss}` is the time scale of mass loss. For disks, the mass outside of the tidal radius is computed in
    cylindrical radius (i.e. treating the disk as thin); for spheroids and nuclear star clusters it is computed in spherical
    radius. As such, no mass is lost until the tidal radius has shrunk to the point where it encloses less than the entire
    component.

    By default, :math:`T_\mathrm{loss}` is taken to be the orbital time scale

    .. math::

       T_\mathrm{orb} = {2 \pi \over \hbox{max}(\omega,v_\mathrm{r}/r)},

    where :math:`\omega` is the angular velocity of the satellite, :math:`v_\mathrm{r}` is the radial velocity, :math:`r` is the
    orbital radius. If ``[useDynamicalTimeScale]`` is set to true, :math:`T_\mathrm{loss}` is instead taken to be the dynamical
    time scale computed at the tidal radius

    .. math::

       T_\mathrm{dyn} = 2 \pi \sqrt{\frac{r_\mathrm{tidal}^3}{16 G M_\mathrm{sat}(r_\mathrm{tidal})}},

    where :math:`M_\mathrm{sat}(r_\mathrm{tidal})` is the total mass of the satellite enclosed within the tidal radius.

    Stripped material carries the mean specific angular momentum of the material outside of the tidal radius. Consistent with the
    assumption made by the galactic components in solving for their sizes, the rotation curve is assumed to be flat, such that the
    specific angular momentum is proportional to radius, :math:`j(r) \propto r`. This is assumed to hold out to a truncation
    radius, :math:`r_\mathrm{t}=`\ ``[radiusTruncateOverRadiusScale]``\ :math:`\times r_\mathrm{s}`, where :math:`r_\mathrm{s}` is
    the scale radius of the component, beyond which the specific angular momentum is taken to be constant, :math:`j(r) \propto
    r_\mathrm{t}`. The specific angular momentum of the stripped material, in units of the mean specific angular momentum of the
    component, is then

    .. math::

       \frac{j_\mathrm{stripped}}{\langle j \rangle} = \frac{\int_{r_\mathrm{tidal}}^\infty \hbox{min}(r,r_\mathrm{t}) \mathrm{d}M / M(&gt;r_\mathrm{tidal})}{\int_0^\infty \hbox{min}(r,r_\mathrm{t}) \mathrm{d}M / M},

    where radii are cylindrical for disks, and spherical otherwise. The truncation is required since the mean radius of some
    profiles (e.g. that of :cite:t:`hernquist_analytical_1990`) is formally infinite. Without it, material at very large radii
    would dominate the angular momentum removed, leaving the remnant with almost no angular momentum. With it, the ratio is
    bounded, :math:`j_\mathrm{stripped}/\langle j \rangle \le r_\mathrm{t}/\langle \hbox{min}(r,r_\mathrm{t}) \rangle`. Removing
    angular momentum in this way, rather than at the mean specific angular momentum of the component, results in the remnant
    being more compact.
    </description>
  </tidalStripping>
  !!]
  type, extends(tidalStrippingClass) :: tidalStrippingZentner2005
     !!{RST
     Implementation of a tidal stripping class for galactic components which follows the model of :cite:t:`zentner_physics_2005`.
     !!}
     private
     class           (satelliteTidalStrippingRadiusClass), pointer :: satelliteTidalStrippingRadius_ => null()
     class           (darkMatterHaloScaleClass          ), pointer :: darkMatterHaloScale_           => null()
     double precision                                              :: efficiency                              , radiusTruncateOverRadiusScale
     logical                                                       :: useDynamicalTimeScale
   contains
     final     ::                 zentner2005Destructor
     procedure :: rateMassLoss => zentner2005RateMassLoss
  end type tidalStrippingZentner2005

  interface tidalStrippingZentner2005
     !!{RST
     Constructors for the :galacticus-class:`tidalStrippingZentner2005` tidal stripping class.
     !!}
     module procedure zentner2005ConstructorParameters
     module procedure zentner2005ConstructorInternal
  end interface tidalStrippingZentner2005

contains

  function zentner2005ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`tidalStrippingZentner2005` tidal stripping class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (tidalStrippingZentner2005         )                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (satelliteTidalStrippingRadiusClass), pointer       :: satelliteTidalStrippingRadius_
    class           (darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_
    double precision                                                    :: efficiency                    , radiusTruncateOverRadiusScale
    logical                                                             :: useDynamicalTimeScale

    !![
    <inputParameter docformat="rst">
      <name>efficiency</name>
      <defaultValue>2.5d0</defaultValue>
      <description>
      The dimensionless rate coefficient appearing in the :cite:t:`zentner_physics_2005` expression for the tidal mass loss rate from galactic components.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>useDynamicalTimeScale</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, the mass outside the tidal radius is assumed to be lost on the dynamical time scale computed at the tidal radius. Otherwise, mass loss occurs on the orbital timescale of the satellite.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusTruncateOverRadiusScale</name>
      <defaultValue>5.0d0</defaultValue>
      <description>
      The radius, in units of the scale radius of the component, beyond which the specific angular momentum of material in the component is assumed to be constant when computing the specific angular momentum of tidally-stripped material.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="satelliteTidalStrippingRadius" name="satelliteTidalStrippingRadius_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"           name="darkMatterHaloScale_"           source="parameters"/>
    !!]
    self=tidalStrippingZentner2005(efficiency,radiusTruncateOverRadiusScale,useDynamicalTimeScale,satelliteTidalStrippingRadius_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStrippingRadius_"/>
    <objectDestructor name="darkMatterHaloScale_"          />
    !!]
    return
  end function zentner2005ConstructorParameters

  function zentner2005ConstructorInternal(efficiency,radiusTruncateOverRadiusScale,useDynamicalTimeScale,satelliteTidalStrippingRadius_,darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`tidalStrippingZentner2005` tidal stripping class.
    !!}
    implicit none
    type            (tidalStrippingZentner2005         )                        :: self
    class           (satelliteTidalStrippingRadiusClass), intent(in   ), target :: satelliteTidalStrippingRadius_
    class           (darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                    , intent(in   )         :: efficiency                    , radiusTruncateOverRadiusScale
    logical                                             , intent(in   )         :: useDynamicalTimeScale
    !![
    <constructorAssign variables="efficiency, radiusTruncateOverRadiusScale, useDynamicalTimeScale, *satelliteTidalStrippingRadius_, *darkMatterHaloScale_"/>
    !!]

    ! Ensure that required methods are supported.
    !![
    <componentPropertyAssert class="satellite" properties="position velocity" require="gettable"/>
    !!]
    return
  end function zentner2005ConstructorInternal

  subroutine zentner2005Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`tidalStrippingZentner2005` tidal stripping class.
    !!}
    implicit none
    type(tidalStrippingZentner2005), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalStrippingRadius_"/>
    <objectDestructor name="self%darkMatterHaloScale_"          />
    !!]
    return
  end subroutine zentner2005Destructor

  double precision function zentner2005RateMassLoss(self,component,ratioSpecificAngularMomentum)
    !!{RST
    Computes the mass loss rate due to tidal stripping of a galactic component using the formulation of :cite:t:`zentner_physics_2005`.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galactic_Structure_Options      , only : componentTypeDisk           , componentTypeSpheroid         , componentTypeNuclearStarCluster, massTypeBaryonic, &
         &                                          enumerationComponentTypeType
    use :: Galacticus_Nodes                , only : nodeComponentDisk           , nodeComponentSpheroid         , nodeComponentNSC               , treeNode
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear                    , gravitationalConstant_internal, megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Satellite_Orbits                , only : Satellite_Orbit_Period
    implicit none
    class           (tidalStrippingZentner2005   ), intent(inout)           :: self
    class           (nodeComponent               ), intent(inout)           :: component
    double precision                              , intent(  out), optional :: ratioSpecificAngularMomentum
    type            (treeNode                    ), pointer                 :: node
    class           (massDistributionClass       ), pointer                 :: massDistribution_           , massDistributionComponent_
    type            (enumerationComponentTypeType)                          :: componentType
    logical                                                                 :: isCylindrical
    double precision                                                        :: massComponent               , massEnclosedTidalRadius   , &
         &                                                                     massOuter                   , radiusTidal               , &
         &                                                                     timeScaleMassLoss           , massSatelliteTidalRadius  , &
         &                                                                     radiusTruncate              , momentOuter               , &
         &                                                                     momentTotal                 , massTruncate              , &
         &                                                                     radiusScale

    ! Assume no mass loss, with stripped material carrying the mean specific angular momentum of the component, by default.
    zentner2005RateMassLoss=0.0d0
    if (present(ratioSpecificAngularMomentum)) ratioSpecificAngularMomentum=1.0d0
    ! Return immediately if this is not a satellite.
    node => component%hostNode
    if (.not.node%isSatellite()) return
    ! Get component properties.
    select type (component)
    class is (nodeComponentDisk    )
       componentType=componentTypeDisk
       massComponent=component%massGas()+component%massStellar()
       radiusScale  =component%radius ()
       isCylindrical=.true.
    class is (nodeComponentSpheroid)
       componentType=componentTypeSpheroid
       massComponent=component%massGas()+component%massStellar()
       radiusScale  =component%radius ()
       isCylindrical=.false.
    class is (nodeComponentNSC     )
       componentType=componentTypeNuclearStarCluster
       massComponent=component%massGas()+component%massStellar()
       radiusScale  =component%radius ()
       isCylindrical=.false.
    class default
       componentType=componentTypeDisk
       massComponent=0.0d0
       radiusScale  =0.0d0
       isCylindrical=.false.
       call Error_Report('unsupported component'//{introspection:location})
    end select
    if (massComponent <= 0.0d0) return
    ! Find the mass of the component outside of the tidal radius.
    radiusTidal                =  self%satelliteTidalStrippingRadius_%radius(node)
    massDistributionComponent_ => node%massDistribution(componentType=componentType,massType=massTypeBaryonic)
    if (radiusTidal > 0.0d0) then
       massEnclosedTidalRadius   =max(0.0d0,massEnclosed(radiusTidal))
    else
       massEnclosedTidalRadius   =0.0d0
    end if
    massOuter=max(0.0d0,massComponent-massEnclosedTidalRadius)
    if (massOuter > 0.0d0) then
       ! Find the time scale for mass loss.
       timeScaleMassLoss=Satellite_Orbit_Period(node,self%darkMatterHaloScale_)
       if (self%useDynamicalTimeScale .and. radiusTidal > 0.0d0) then
          massDistribution_        => node             %massDistribution    (           )
          massSatelliteTidalRadius =  massDistribution_%massEnclosedBySphere(radiusTidal)
          !![
          <objectDestructor name="massDistribution_"/>
          !!]
          if (massSatelliteTidalRadius > 0.0d0)                          &
               & timeScaleMassLoss=+2.0d0                                &
               &                   *Pi                                   &
               &                   *sqrt(                                &
               &                         +radiusTidal**3                 &
               &                         /16.0d0                         &
               &                         /gravitationalConstant_internal &
               &                         /massSatelliteTidalRadius       &
               &                        )                                &
               &                   *megaParsec                           &
               &                   /kilo                                 &
               &                   /gigaYear
       end if
       ! Compute the rate of mass loss.
       zentner2005RateMassLoss=+self%efficiency   &
            &                  *massOuter         &
            &                  /timeScaleMassLoss
       ! Compute the specific angular momentum of the stripped material in units of the mean specific angular momentum of the
       ! component. For a flat rotation curve the specific angular momentum is proportional to radius. This is assumed to hold out
       ! to some multiple of the scale radius of the component, beyond which the specific angular momentum is held constant.
       if (present(ratioSpecificAngularMomentum) .and. radiusScale > 0.0d0) then
          radiusTruncate=+self%radiusTruncateOverRadiusScale &
               &         *radiusScale
          massTruncate  =max(0.0d0,massComponent-massEnclosed(radiusTruncate))
          momentTotal   =+momentRadial(0.0d0,radiusTruncate) &
               &         +massTruncate                       &
               &         *radiusTruncate
          if (radiusTidal < radiusTruncate) then
             momentOuter=+momentRadial(radiusTidal,radiusTruncate) &
                  &      +massTruncate                             &
                  &      *radiusTruncate
          else
             momentOuter=+massOuter                                &
                  &      *radiusTruncate
          end if
          if (momentTotal > 0.0d0) ratioSpecificAngularMomentum=+momentOuter   &
               &                                                /massOuter     &
               &                                                /momentTotal   &
               &                                                *massComponent
       end if
    end if
    !![
    <objectDestructor name="massDistributionComponent_"/>
    !!]
    return

  contains

    double precision function massEnclosed(radius)
      !!{RST
      Return the mass of the component enclosed within the given radius. Cylindrical radius is used for disks, spherical radius
      otherwise.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      if (isCylindrical) then
         massEnclosed=massDistributionComponent_%massEnclosedByCylinder(radius)
      else
         massEnclosed=massDistributionComponent_%massEnclosedBySphere  (radius)
      end if
      return
    end function massEnclosed

    double precision function momentRadial(radiusMinimum,radiusMaximum)
      !!{RST
      Return the mass-weighted radial moment, :math:`\int_{r_\mathrm{min}}^{r_\mathrm{max}} r \mathrm{d}M`, of the component.
      Cylindrical radius is used for disks, spherical radius otherwise.
      !!}
      implicit none
      double precision, intent(in   ) :: radiusMinimum, radiusMaximum

      if (isCylindrical) then
         momentRadial=+2.0d0                                                                                                                &
              &       *Pi                                                                                                                   &
              &       *massDistributionComponent_%surfaceDensityRadialMoment(2.0d0,radiusMinimum=radiusMinimum,radiusMaximum=radiusMaximum)
      else
         momentRadial=+4.0d0                                                                                                                &
              &       *Pi                                                                                                                   &
              &       *massDistributionComponent_%densityRadialMoment       (3.0d0,radiusMinimum=radiusMinimum,radiusMaximum=radiusMaximum)
      end if
      return
    end function momentRadial

  end function zentner2005RateMassLoss
