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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a dust attenuation class describing a clumpy torus of dust around an AGN.
  !!}

  use :: Black_Hole_Torus_Covering_Fractions, only : blackHoleTorusCoveringFractionClass
  use :: Dust_Extinction_Curves             , only : dustExtinctionCurveClass
  use :: Galactic_Inclinations              , only : galacticInclinationClass

  !![
  <dustAttenuation name="dustAttenuationAGNTorus" docformat="rst">
   <description>
   A clumpy torus of dust around the accretion disk of an AGN, following :cite:t:`nenkova_agn_2008`. Clouds of dust,
   each of :math:`V`-band optical depth ``[depthOpticalCloudV]``, are distributed about the equatorial plane of the
   torus such that the mean number met along a ray leaving the accretion disk at an angle :math:`\beta` from that plane
   is

   .. math::

      N(\beta) = N_0 \exp\left( -\frac{\beta^2}{\sigma^2} \right),

   where :math:`N_0` is ``[countCloudsEquatorial]``. For clouds placed at random (Poisson distributed) along the ray, the
   mean transmission is that of a screen of effective optical depth
   :math:`N(\beta)[1-\exp(-\tau_\lambda)]` :cite:p:`natta_extinction_1984`,

   .. math::

      T(\lambda,\beta) = \exp\left\{ -N(\beta) \left[ 1 - \exp\left( -\tau_\mathrm{V} \frac{k(\lambda)}{k_\mathrm{V}} \right) \right] \right\},

   with the wavelength dependence :math:`k(\lambda)/k_\mathrm{V}` taken from the ``[dustExtinctionCurve]`` object. One
   minus this transmission is counted as absorbed, which counts light scattered by the clouds as absorbed too.

   **Covering fraction.** The torus covers a fraction

   .. math::

      f_\mathrm{cov} = \int_0^{\pi/2} \left\{ 1 - \exp\left[ -N(\beta) \right] \right\} \cos\beta\,\mathrm{d}\beta

   of the sky as seen from the accretion disk---the probability that a ray meets at least one cloud. The angular width
   :math:`\sigma` is chosen for each galaxy so that this equals the covering fraction given by the
   ``[blackHoleTorusCoveringFraction]`` object, which may depend on the luminosity of the AGN (as in the receding torus
   model). The largest covering fraction a torus of :math:`N_0` clouds can have is :math:`1-\exp(-N_0)`, reached as
   :math:`\sigma\rightarrow\infty`; a larger one is reduced to it. Because the clouds thin out smoothly away from the
   equatorial plane, rather than ending at a sharp edge, the transmission is a smooth function of orientation, so that
   averages over orientation by quadrature---as made by :galacticus-class:`dustAttenuationInclinationAveraged` and when
   finding the luminosity absorbed by dust---converge rapidly.

   **What is attenuated.** Only the continuum emission of the accretion disk of the central black hole
   (``emissionSourceAccretionDisk`` emitted by ``componentTypeBlackHole``) is attenuated. The narrow-line region lies
   beyond the torus, so emission lines from the AGN (``emissionSourceNebular``) are not attenuated, nor is any emission
   from other components. To see AGN light through the dust of the host galaxy as well, combine this attenuator with
   others through :galacticus-class:`dustAttenuationSequence`, listing it first, since the torus is where the light is
   emitted.

   **Orientation.** The axis of the torus is taken to be that of its host galaxy, so the torus is seen at the
   inclination given by the ``[galacticInclination]`` object (or by the ``inclination`` passed in when averaging over
   orientation), with :math:`\beta=\pi/2-i`.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationAGNTorus
     !!{RST
     A dust attenuation class describing a clumpy torus of dust around an AGN.
     !!}
     private
     class           (dustExtinctionCurveClass           ), pointer                   :: dustExtinctionCurve_            => null()
     class           (galacticInclinationClass           ), pointer                   :: galacticInclination_            => null()
     class           (blackHoleTorusCoveringFractionClass), pointer                   :: blackHoleTorusCoveringFraction_ => null()
     double precision                                                                 :: countCloudsEquatorial                    , depthOpticalCloudV
     logical                                                                          :: inclinationAvailable
     ! Tabulation of the covering fraction as a function of the logarithm of the angular width of the torus, used to
     ! find the width giving a required covering fraction.
     double precision                                     , allocatable, dimension(:) :: widthAngularLogarithmic                  , coveringFractionTabulated
   contains
     !![
     <methods docformat="rst">
       <method method="coveringFractionWidth" description="Return the covering fraction of a torus of the given angular width, :math:`\sigma`, in radians."                   />
       <method method="widthAngular"          description="Return the angular width, :math:`\sigma` in radians, of the torus around the central black hole of the given node."/>
     </methods>
     !!]
     final     ::                           agnTorusDestructor
     procedure :: transmission           => agnTorusTransmission
     procedure :: isOrientationDependent => agnTorusIsOrientationDependent
     procedure :: coveringFractionWidth  => agnTorusCoveringFractionWidth
     procedure :: widthAngular           => agnTorusWidthAngular
  end type dustAttenuationAGNTorus

  interface dustAttenuationAGNTorus
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationAGNTorus` dust attenuation class.
     !!}
     module procedure agnTorusConstructorParameters
     module procedure agnTorusConstructorInternal
  end interface dustAttenuationAGNTorus

  ! Range and number of points of the tabulation of covering fraction against the angular width of the torus, in radians.
  double precision, parameter :: agnTorusWidthAngularMinimum=1.0d-3, agnTorusWidthAngularMaximum=1.0d+2
  integer         , parameter :: agnTorusCountTabulation    =1000

contains

  function agnTorusConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationAGNTorus` dust attenuation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustAttenuationAGNTorus            )                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (dustExtinctionCurveClass           ), pointer       :: dustExtinctionCurve_
    class           (galacticInclinationClass           ), pointer       :: galacticInclination_
    class           (blackHoleTorusCoveringFractionClass), pointer       :: blackHoleTorusCoveringFraction_
    double precision                                                     :: countCloudsEquatorial          , depthOpticalCloudV

    !![
    <inputParameter docformat="rst">
      <name>countCloudsEquatorial</name>
      <defaultValue>5.0d0</defaultValue>
      <defaultSource>:cite:t:`nenkova_agn_2008`, who find 5--15</defaultSource>
      <description>
      The mean number of clouds, :math:`N_0`, met along a ray in the equatorial plane of the torus.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>depthOpticalCloudV</name>
      <defaultValue>30.0d0</defaultValue>
      <defaultSource>:cite:t:`nenkova_agn_2008`, who find 30--100</defaultSource>
      <description>
      The :math:`V`-band optical depth, :math:`\tau_\mathrm{V}`, of each cloud.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustExtinctionCurve"            name="dustExtinctionCurve_"            source="parameters"/>
    <objectBuilder class="galacticInclination"            name="galacticInclination_"            source="parameters"/>
    <objectBuilder class="blackHoleTorusCoveringFraction" name="blackHoleTorusCoveringFraction_" source="parameters"/>
    !!]
    self=dustAttenuationAGNTorus(countCloudsEquatorial,depthOpticalCloudV,dustExtinctionCurve_,galacticInclination_,blackHoleTorusCoveringFraction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustExtinctionCurve_"           />
    <objectDestructor name="galacticInclination_"           />
    <objectDestructor name="blackHoleTorusCoveringFraction_"/>
    !!]
    return
  end function agnTorusConstructorParameters

  function agnTorusConstructorInternal(countCloudsEquatorial,depthOpticalCloudV,dustExtinctionCurve_,galacticInclination_,blackHoleTorusCoveringFraction_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationAGNTorus` dust attenuation class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (dustAttenuationAGNTorus            )                        :: self
    double precision                                     , intent(in   )         :: countCloudsEquatorial          , depthOpticalCloudV
    class           (dustExtinctionCurveClass           ), intent(in   ), target :: dustExtinctionCurve_
    class           (galacticInclinationClass           ), intent(in   ), target :: galacticInclination_
    class           (blackHoleTorusCoveringFractionClass), intent(in   ), target :: blackHoleTorusCoveringFraction_
    integer                                                                      :: i
    !![
    <constructorAssign variables="countCloudsEquatorial, depthOpticalCloudV, *dustExtinctionCurve_, *galacticInclination_, *blackHoleTorusCoveringFraction_"/>
    !!]

    if (countCloudsEquatorial <= 0.0d0) call Error_Report('`countCloudsEquatorial` must be positive'   //{introspection:location})
    if (depthOpticalCloudV    <  0.0d0) call Error_Report('`depthOpticalCloudV` must be non-negative'//{introspection:location})
    self%inclinationAvailable=self%galacticInclination_%isAvailable()
    ! Tabulate the covering fraction against the logarithm of the angular width. It increases monotonically with width.
    allocate(self%widthAngularLogarithmic  (agnTorusCountTabulation))
    allocate(self%coveringFractionTabulated(agnTorusCountTabulation))
    do i=1,agnTorusCountTabulation
       self%widthAngularLogarithmic  (i)=+log(agnTorusWidthAngularMinimum                            ) &
            &                            +log(agnTorusWidthAngularMaximum/agnTorusWidthAngularMinimum) &
            &                            *dble(i                      -1)                              &
            &                            /dble(agnTorusCountTabulation-1)
       self%coveringFractionTabulated(i)=self%coveringFractionWidth(exp(self%widthAngularLogarithmic(i)))
    end do
    return
  end function agnTorusConstructorInternal

  subroutine agnTorusDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationAGNTorus` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationAGNTorus), intent(inout) :: self

    !![
    <objectDestructor name="self%dustExtinctionCurve_"           />
    <objectDestructor name="self%galacticInclination_"           />
    <objectDestructor name="self%blackHoleTorusCoveringFraction_"/>
    !!]
    return
  end subroutine agnTorusDestructor

  double precision function agnTorusCoveringFractionWidth(self,widthAngular) result(coveringFraction)
    !!{RST
    Return the covering fraction of a torus of angular width ``widthAngular`` (:math:`\sigma`, in radians),
    :math:`\int_0^{\pi/2} \{1-\exp[-N(\beta)]\} \cos\beta\,\mathrm{d}\beta`.

    The integrand is negligible once :math:`N(\beta)` is, so the integral is taken only out to the angle at which
    :math:`N(\beta)=10^{-16}`, by a 64-point Gauss-Legendre rule. That keeps the integrand resolved however narrow the
    torus.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (dustAttenuationAGNTorus), intent(inout)               :: self
    double precision                         , intent(in   )               :: widthAngular
    integer                                  , parameter                   :: orderQuadrature  =64
    ! The mean number of clouds below which a ray is considered unobscured.
    double precision                         , parameter                   :: countCloudsSmall=1.0d-16
    double precision                         , allocatable  , dimension(:) :: abscissae               , weights
    double precision                                                       :: angleMaximum
    integer                                                                :: i

    call gaussLegendreRule(orderQuadrature,abscissae,weights)
    angleMaximum=min(                                                                &
         &           Pi/2.0d0,                                                       &
         &           widthAngular*sqrt(max(log(self%countCloudsEquatorial/countCloudsSmall),0.0d0)) &
         &          )
    coveringFraction=0.0d0
    do i=1,orderQuadrature
       coveringFraction=+coveringFraction                                                                           &
            &           +weights(i)                                                                                 &
            &           *(1.0d0-exp(-self%countCloudsEquatorial*exp(-(abscissae(i)*angleMaximum/widthAngular)**2))) &
            &           *cos(abscissae(i)*angleMaximum)
    end do
    ! The rule averages over [0,1], so scale by the width of the interval integrated over.
    coveringFraction=coveringFraction*angleMaximum
    return
  end function agnTorusCoveringFractionWidth

  double precision function agnTorusWidthAngular(self,node) result(widthAngular)
    !!{RST
    Return the angular width, :math:`\sigma` in radians, of the torus around the central black hole of ``node``: that
    for which its covering fraction equals the one required. The tabulated covering fraction is interpolated linearly
    in the logarithm of the width. A covering fraction outside the range tabulated is given the width at the end of the
    range---which, at the upper end, is the largest covering fraction the torus can have.
    !!}
    implicit none
    class           (dustAttenuationAGNTorus), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    double precision                                         :: coveringFraction, factor
    integer                                                  :: i               , iLower, &
         &                                                      iUpper

    coveringFraction=self%blackHoleTorusCoveringFraction_%coveringFraction(node)
    if      (coveringFraction <= self%coveringFractionTabulated(                      1)) then
       widthAngular=agnTorusWidthAngularMinimum
    else if (coveringFraction >= self%coveringFractionTabulated(agnTorusCountTabulation)) then
       widthAngular=agnTorusWidthAngularMaximum
    else
       ! Bisect to find the interval containing the required covering fraction.
       iLower=1
       iUpper=agnTorusCountTabulation
       do while (iUpper-iLower > 1)
          i=(iLower+iUpper)/2
          if (self%coveringFractionTabulated(i) > coveringFraction) then
             iUpper=i
          else
             iLower=i
          end if
       end do
       factor      =+(     coveringFraction                 -self%coveringFractionTabulated(iLower)) &
            &       /(self%coveringFractionTabulated(iUpper)-self%coveringFractionTabulated(iLower))
       widthAngular=exp(                                                     &
            &           +self%widthAngularLogarithmic(iLower)*(1.0d0-factor) &
            &           +self%widthAngularLogarithmic(iUpper)*       factor  &
            &          )
    end if
    return
  end function agnTorusWidthAngular

  function agnTorusTransmission(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the transmission through a clumpy torus of dust. Only accretion disk emission from the central black hole is
    attenuated; everything else is transmitted in full. The angular width of the torus, and the number of clouds along
    the line of sight, are found once for all parcels.
    !!}
    use :: Dust_Attenuation_Descriptors, only : emissionSourceAccretionDisk
    use :: Error                       , only : Error_Report
    use :: Galactic_Structure_Options  , only : componentTypeBlackHole
    use :: Numerical_Constants_Math    , only : Pi
    implicit none
    class           (dustAttenuationAGNTorus), intent(inout)                               :: self
    type            (treeNode               ), intent(inout), target                       :: node
    type            (emissionDescriptor     ), intent(in   ), dimension(:                ) :: descriptors
    double precision                         , intent(in   ), optional                     :: inclination
    double precision                                        , dimension(size(descriptors)) :: transmission
    double precision                                                                       :: inclination_, countClouds, &
         &                                                                                    widthAngular
    logical                                                                                :: computed
    integer                                                                                :: i

    transmission=1.0d0
    computed    =.false.
    countClouds =0.0d0
    do i=1,size(descriptors)
       if     (                                                             &
            &   descriptors(i)%componentType /= componentTypeBlackHole      &
            &  .or.                                                         &
            &   descriptors(i)%sourceType    /= emissionSourceAccretionDisk &
            & ) cycle
       if (.not.computed) then
          if (present(inclination)) then
             inclination_=inclination
          else if (self%inclinationAvailable) then
             inclination_=self%galacticInclination_%inclination(node)
          else
             inclination_=0.0d0
             call Error_Report('this attenuator depends on orientation, but no inclination is available: either set `galacticInclination` to a class which supplies one, or wrap this attenuator in `inclinationAveraged`'//{introspection:location})
          end if
          widthAngular=self%widthAngular(node)
          countClouds =+self%countCloudsEquatorial                      &
               &       *exp(-((Pi/2.0d0-inclination_)/widthAngular)**2)
          computed    =.true.
       end if
       transmission(i)=exp(                                                                                 &
            &              -countClouds                                                                     &
            &              *(                                                                               &
            &                +1.0d0                                                                         &
            &                -exp(                                                                          &
            &                     -self%depthOpticalCloudV                                                  &
            &                     *self%dustExtinctionCurve_%attenuationRelative(descriptors(i)%wavelength) &
            &                    )                                                                          &
            &               )                                                                               &
            &             )
    end do
    return
  end function agnTorusTransmission

  logical function agnTorusIsOrientationDependent(self) result(isOrientationDependent)
    !!{RST
    Return true: the transmission through the torus depends on the inclination at which it is seen.
    !!}
    implicit none
    class(dustAttenuationAGNTorus), intent(inout) :: self
    !$GLC attributes unused :: self

    isOrientationDependent=.true.
    return
  end function agnTorusIsOrientationDependent
