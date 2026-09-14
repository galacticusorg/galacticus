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
Contains a module which provides a class implementing attenuation of galactic emission by dust.
!!}

module Dust_Attenuations
  !!{RST
  Provides a class implementing attenuation of galactic emission by dust.

  A ``dustAttenuation`` object answers one question: of the light emitted by a given parcel of a given galaxy, what
  fraction escapes? It is handed a ``node`` and a list of ``emissionDescriptor`` parcels---each one
  wavelength, one component, one source type, one age range---and returns the fraction of that light which reaches an
  observer.

  Two design points are worth stating explicitly.

  **Transmission, not magnitudes.** The method returns a multiplicative factor rather than an attenuation in
  magnitudes or an optical depth. This is unambiguous, composes by multiplication (so a sequence of attenuators is
  simply a product), and remains well defined for models---such as radiative transfer through a specified
  geometry---whose result is not proportional to any single optical depth.

  **Vectorized over parcels.** ``transmission`` receives every parcel at once so that an implementation can hoist its
  per-node work---the :math:`V`-band optical depth of each component, the inclination, interpolation factors in a
  tabulated atlas---out of the loop over wavelength, without caching state on the object and risking staleness.

  How finely a producer must split its luminosity is negotiated: ``request`` reports the age bin boundaries an
  attenuator distinguishes and which other axes it depends upon, so that an age-independent screen applied to a
  broad-band luminosity is handed one parcel per component, rather than one per wavelength, age and metallicity.

  **Absorption.** Light which does not reach the observer has either been absorbed by dust, or scattered into some
  other direction. Only the former heats the dust, and so only it is re-emitted in the infrared. ``absorbedFractions``
  therefore reports what each *phase* of dust absorbs---birth clouds and the diffuse interstellar medium, for
  example---averaged over orientation, which is the quantity a dust emission model needs, and which differs from one
  minus the directional ``transmission`` of an orientation-dependent attenuator.

  Note that any constant which must be shared between implementations of this class has to be declared here rather
  than in an implementation file, because each implementation is generated into its own submodule and so cannot see
  module-level declarations made by its siblings.
  !!}
  use :: Dust_Attenuation_Descriptors, only : decompositionRequest  , emissionDescriptor
  use :: Dust_Properties             , only : componentGasProperties, densitySurfaceGasDepthOpticalVUnitMilkyWay, dustPropertiesClass
  use :: Galactic_Structure_Options  , only : componentTypeAll      , enumerationComponentTypeType
  use :: Galacticus_Nodes            , only : treeNode
  use :: ISO_Varying_String          , only : varying_string
  private
  ! Made public so that it survives into the object file: it is called only from the submodules into which the
  ! implementations of this class are generated, never from this module itself, and a private procedure with no
  ! caller in its own module can be discarded before those submodules are linked against it.
  public :: absorbedFractionsPhases, gaussLegendreRule, radiusSpheroidRelative

  !![
  <functionClass docformat="rst">
   <name>dustAttenuation</name>
   <descriptiveName>Dust Attenuation</descriptiveName>
   <description>
   Class computing the fraction of the emission from a galaxy which is transmitted through its dust.
   </description>
   <default>zero</default>
   <method name="transmission" >
    <description>
    Return the fraction of the emission transmitted through dust, for each of the given parcels of emission from the
    given ``node``. A value of unity indicates no attenuation. The result has one element per element of
    ``descriptors``, in the same order.

    The value is non-negative, and is normally at most unity---but not necessarily. This is a *directional*
    transmission, the fraction reaching an observer in one particular direction, and where dust scatters light it can
    redirect more into that direction than it removes from it, leaving the galaxy brighter at that angle than it
    would be with no dust at all. Radiative transfer calculations show this at low optical depth and low inclination:
    the atlas of :cite:t:`ferrara_atlas_1999` exceeds unity by up to 3% there. Consumers must therefore not assume an
    upper bound of one. Energy is still conserved, since what is gained along one line of sight is lost along
    others---a constraint on the average over orientation, not on any single direction.

    ``inclination`` overrides the angle at which an orientation-dependent attenuator is evaluated, in radians. It is
    how :galacticus-class:`dustAttenuationInclinationAveraged` drives its quadrature: the wrapped attenuator is
    asked for the transmission at each angle in turn, without any shared state being mutated, so the arrangement is
    safe under threading. An implementation which does not depend on orientation ignores it; one which does falls
    back to its own ``galacticInclination`` object when it is absent.
    </description>
    <type>double precision, dimension(size(descriptors))</type>
    <pass>yes</pass>
    <argument>type(treeNode          ), intent(inout), target       :: node       </argument>
    <argument>type(emissionDescriptor), intent(in   ), dimension(:) :: descriptors</argument>
    <argument>double precision        , intent(in   ), optional     :: inclination</argument>
   </method>
   <method name="request" >
    <description>
    Return the resolution which a decomposition must have for this attenuator to be applied correctly: the boundaries
    of the age bins it distinguishes, and whether it depends upon component, metallicity, or radius. A producer of
    luminosities uses this to split its output no more finely than necessary.

    The default requests no age resolution and no metallicity or radius resolution, but does request splitting by
    component---which is what any attenuator depending on the properties of an individual component needs.
    </description>
    <type>type(decompositionRequest)</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     dustAttenuationRequest%resolveComponents =.true.
     dustAttenuationRequest%resolveMetallicity=.false.
     dustAttenuationRequest%resolveRadius     =.false.
    </code>
   </method>
   <method name="supportsComponent" >
    <description>
    Return true if this attenuator can be applied to emission from the given component.

    The default accepts any individual component but rejects ``componentTypeAll``: different components are attenuated
    differently, so a luminosity summed over components can not meaningfully be attenuated---it must be decomposed,
    attenuated, and only then summed. This is checked when the attenuating property extractor is constructed, so that
    a misconfiguration is reported immediately rather than part way through a run.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(enumerationComponentTypeType), intent(in   ) :: componentType</argument>
    <code>
     !$GLC attributes unused :: self
     dustAttenuationSupportsComponent=(componentType /= componentTypeAll)
    </code>
   </method>
   <method name="countPhases" >
    <description>
    Return the number of distinct phases of dust---birth clouds and the diffuse interstellar medium, for example---into
    which this attenuator divides the absorption of light, so that the energy absorbed by each can be re-emitted with a
    spectrum of its own. The default is a single phase.
    </description>
    <type>integer</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     dustAttenuationCountPhases=1
    </code>
   </method>
   <method name="labelPhase" >
    <description>
    Return a label for the phase of dust with index ``indexPhase``, used to name the luminosity it absorbs. The default
    is the short name of the class.
    </description>
    <type>type(varying_string)</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: indexPhase</argument>
    <code>
     !$GLC attributes unused :: indexPhase
     dustAttenuationLabelPhase=self%objectType(short=.true.)
    </code>
   </method>
   <method name="transmissionPhases" >
    <description>
    Return the transmission through each phase of dust separately, for each of the given parcels of emission. The first
    index of the result runs over parcels and the second over phases, which are ordered along the path of the light:
    the first phase is the one in which the light is emitted. The product over phases is the transmission returned by
    ``transmission``, and ``inclination`` has the same meaning as there. The default is a single phase, whose
    transmission is that returned by ``transmission``.
    </description>
    <type>double precision, allocatable, dimension(:,:)</type>
    <pass>yes</pass>
    <argument>type(treeNode          ), intent(inout), target       :: node       </argument>
    <argument>type(emissionDescriptor), intent(in   ), dimension(:) :: descriptors</argument>
    <argument>double precision        , intent(in   ), optional     :: inclination</argument>
    <code>
     allocate(dustAttenuationTransmissionPhases(size(descriptors),1))
     dustAttenuationTransmissionPhases(:,1)=self%transmission(node,descriptors,inclination)
    </code>
   </method>
   <method name="isOrientationDependent" >
    <description>
    Return true if the transmission depends on the orientation of the galaxy relative to the observer. The luminosity
    absorbed by dust must then be averaged over orientation, since light scattered out of one line of sight escapes
    along another rather than being absorbed. The default is false.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     dustAttenuationIsOrientationDependent=.false.
    </code>
   </method>
   <method name="absorbedFractions" >
    <description>
    Return the fraction of the emission of each of the given parcels which is absorbed by each phase of dust, averaged
    over orientation. The first index of the result runs over parcels and the second over phases.

    For a phase :math:`k` of transmission :math:`T_k`, reached by light which has already passed through phases
    :math:`1` to :math:`k-1`, the absorbed fraction is :math:`f_k = (1-T_k)\prod_{j&lt;k} T_j`, so that the fractions
    sum to one minus the total transmission. Counting one minus the transmission as absorbed is exact for an attenuator
    which conserves energy once averaged over orientation---a radiative transfer atlas, for example---and is the
    standard interpretation of the *effective* attenuation laws applied through screens, such as that of
    :cite:t:`charlot_simple_2000`, whose transmission already allows on average for light scattered back into the line
    of sight. It over-estimates absorption if a screen is given the extinction curve of the grains themselves, since
    that counts scattered light as absorbed.

    Where the attenuator depends on orientation the fractions are computed at each abscissa of the quadrature rule
    given by ``cosineInclination`` and ``weight``---whose weights must sum to unity, so that the rule averages over
    :math:`\cos i` between zero and one---and then averaged. The product over phases is formed at each orientation
    before averaging, since the average of a product is not the product of the averages. An orientation-independent
    attenuator is evaluated once. Fractions are clamped to be non-negative: a directional transmission may exceed
    unity, and although its average over orientation should not, a tabulation is only so accurate.
    </description>
    <type>double precision, allocatable, dimension(:,:)</type>
    <pass>yes</pass>
    <argument>type(treeNode          ), intent(inout), target       :: node             </argument>
    <argument>type(emissionDescriptor), intent(in   ), dimension(:) :: descriptors      </argument>
    <argument>double precision        , intent(in   ), dimension(:) :: cosineInclination</argument>
    <argument>double precision        , intent(in   ), dimension(:) :: weight           </argument>
    <code>
     integer :: i
     if (self%isOrientationDependent()) then
        allocate(dustAttenuationAbsorbedFractions(size(descriptors),self%countPhases()))
        dustAttenuationAbsorbedFractions=0.0d0
        do i=1,size(weight)
           dustAttenuationAbsorbedFractions=dustAttenuationAbsorbedFractions+weight(i)*absorbedFractionsPhases(self%transmissionPhases(node,descriptors,acos(cosineInclination(i))))
        end do
     else
        dustAttenuationAbsorbedFractions=absorbedFractionsPhases(self%transmissionPhases(node,descriptors))
     end if
     dustAttenuationAbsorbedFractions=max(dustAttenuationAbsorbedFractions,0.0d0)
    </code>
   </method>
  </functionClass>
  !!]

contains

  function absorbedFractionsPhases(transmission) result(fraction)
    !!{RST
    Return the fraction of the emission of each parcel absorbed by each phase of dust, given the transmission through
    each phase in the order the light passes through them: :math:`f_k = (1-T_k)\prod_{j<k} T_j`. The first index of
    both arrays runs over parcels, and the second over phases.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:,:                                       ) :: transmission
    double precision               , dimension(size(transmission,1),size(transmission,2)) :: fraction
    double precision               , dimension(size(transmission,1)                     ) :: transmitted
    integer                                                                               :: k

    transmitted=1.0d0
    do k=1,size(transmission,2)
       fraction   (:,k)=+transmitted*(1.0d0-transmission(:,k))
       transmitted     =+transmitted*       transmission(:,k)
    end do
    return
  end function absorbedFractionsPhases

  subroutine gaussLegendreRule(order,abscissae,weights)
    !!{RST
    Return the abscissae and weights of the Gauss-Legendre rule of the given ``order`` on the interval
    :math:`[0,1]`.

    The nodes are the roots of the Legendre polynomial of that order, found by Newton iteration from the standard
    Chebyshev-like starting guess, with the polynomial and its derivative evaluated by the usual recurrence. Both are
    then mapped from :math:`[-1,1]` onto :math:`[0,1]`, and the weights scaled by the half-width of the interval so
    that they sum to unity---which is what makes the result an average rather than an integral.

    Averages over orientation are needed both by :galacticus-class:`dustAttenuationInclinationAveraged` and by any
    consumer of ``absorbedFractions``, so the rule lives here rather than in an implementation.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    integer                       , intent(in   )               :: order
    double precision, allocatable , intent(inout), dimension(:) :: abscissae                     , weights
    double precision              , parameter                   :: toleranceRelative     =1.0d-15
    integer                       , parameter                   :: countIterationsMaximum=100
    double precision                                            :: root                          , rootPrevious , &
         &                                                         legendre                      , legendreLower, &
         &                                                         legendreLowerLower            , derivative
    integer                                                     :: i                             , j            , &
         &                                                         countIterations

    if (allocated(abscissae)) deallocate(abscissae)
    if (allocated(weights  )) deallocate(weights  )
    allocate(abscissae(order))
    allocate(weights  (order))
    do i=1,order
       ! Initial guess for the i'th root of the Legendre polynomial of this order.
       root=cos(Pi*(dble(i)-0.25d0)/(dble(order)+0.5d0))
       countIterations=0
       do
          ! Evaluate the Legendre polynomial and its derivative at the current estimate by recurrence.
          legendre     =1.0d0
          legendreLower=0.0d0
          do j=1,order
             legendreLowerLower=legendreLower
             legendreLower     =legendre
             legendre          =(dble(2*j-1)*root*legendreLower-dble(j-1)*legendreLowerLower)/dble(j)
          end do
          derivative  =dble(order)*(root*legendre-legendreLower)/(root**2-1.0d0)
          rootPrevious=root
          root        =rootPrevious-legendre/derivative
          countIterations=countIterations+1
          if (abs(root-rootPrevious) <= toleranceRelative*abs(root) .or. countIterations >= countIterationsMaximum) exit
       end do
       ! Map from [-1,1] onto [0,1]. The weights are halved along with the interval, so that they sum to unity and
       ! the quadrature returns a mean.
       abscissae(i)=0.5d0*(1.0d0-root)
       weights  (i)=1.0d0/((1.0d0-root**2)*derivative**2)
    end do
    return
  end subroutine gaussLegendreRule

  double precision function radiusSpheroidRelative(node) result(radiusSpheroid)
    !!{RST
    Return the half-mass radius of the spheroid, in units of the disk scale length.

    Each radiative transfer atlas tabulates the spheroid along an axis of its own, and none of those axes is
    directly a model galaxy's spheroid radius: the atlas simulated a particular density profile, and labeled the
    axis with a particular radius of it. What is returned here is the one measure that means the same thing whatever
    profile either side assumes---the half-mass radius, taken from the stellar mass distribution of the spheroid
    rather than from its scale radius, since for a Hernquist profile the latter is smaller by
    :math:`(1+\sqrt{2})`. Each atlas then converts this to its own axis by dividing by the half-mass radius that one
    unit of that axis corresponds to.

    Matching two differently shaped profiles on a single radius is itself an approximation, and not the best one
    available: :cite:t:`bianchi_monte_carlo_1996` matched an :math:`R^{1/4}` profile to a Jaffe profile by fitting
    their enclosed luminosity, and obtained a relation differing by :math:`\approx 14` percent from what matching
    half-light radii would have given. Attenuations for a spheroid whose profile is not the one an atlas simulated
    should therefore not be relied upon at better than the ten percent level, whatever the optical depth.

    The disk is measured by its scale radius, which is what the atlases normalize to, and which is what the disk
    component's radius already is for an exponential profile.

    A galaxy with no disk has no scale to measure the spheroid against, and no dust either, so the value is
    immaterial: zero is returned, which callers clamp into the tabulated range.

    This lives here rather than in either atlas: each implementation of this class is generated into its own
    submodule, and while a submodule can see its parent module it can not see its siblings.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure_Options, only : componentTypeSpheroid, massTypeStellar
    use :: Galacticus_Nodes          , only : nodeComponentDisk
    use :: Mass_Distributions        , only : massDistributionClass, massDistributionSpherical
    implicit none
    type            (treeNode             ), intent(inout), target  :: node
    class           (nodeComponentDisk    )               , pointer :: disk
    class           (massDistributionClass)               , pointer :: massDistributionSpheroid
    double precision                                                :: radiusDisk

    disk       => node%disk  ()
    radiusDisk =  disk%radius()
    if (radiusDisk <= 0.0d0) then
       radiusSpheroid=0.0d0
       return
    end if
    massDistributionSpheroid => node%massDistribution(componentTypeSpheroid,massTypeStellar)
    select type (massDistributionSpheroid)
    class is (massDistributionSpherical)
       radiusSpheroid=+massDistributionSpheroid%radiusHalfMass() &
            &         /                         radiusDisk
    class default
       radiusSpheroid=0.0d0
       call Error_Report('a half-mass radius is needed for the spheroid, which requires a spherical mass distribution'//{introspection:location})
    end select
    !![
    <objectDestructor name="massDistributionSpheroid"/>
    !!]
    return
  end function radiusSpheroidRelative

end module Dust_Attenuations
