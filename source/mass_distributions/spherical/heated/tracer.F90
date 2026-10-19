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
  Implements a spherical mass distribution of a tracer population heated along with a heated reference mass distribution.
  !!}

  !![
  <massDistribution name="massDistributionSphericalHeatedTracer" docformat="rst">
   <description>
   A mass distribution class describing a tracer population (e.g. the stars of a galactic component) embedded in a reference
   mass distribution (e.g. a dark matter halo) which has been heated (see :galacticus-class:`massDistributionSphericalHeated`
   and :galacticus-class:`massDistributionSphericalHeatedMonotonic`). The tracer is assumed to be heated along with the reference
   distribution, such that a shell of tracer material initially at radius :math:`r` moves to the same final radius,
   :math:`r^\prime`, as the shell of the reference distribution initially at that radius. Since heating assumes no shell
   crossing, :math:`r^\prime` is found from :math:`M^\prime_\mathrm{ref}(r^\prime) = M_\mathrm{ref}(r)`, where
   :math:`M_\mathrm{ref}(r)` and :math:`M^\prime_\mathrm{ref}(r)` are the masses enclosed within radius :math:`r` in the
   unheated and heated reference distributions respectively. The mass of the tracer enclosed within :math:`r^\prime` is then
   :math:`M^\prime(r^\prime) = M(r)`, and its density is

   .. math::

      \rho^\prime(r^\prime) = \rho(r) {\rho^\prime_\mathrm{ref}(r^\prime) \over \rho_\mathrm{ref}(r)},

   where :math:`\rho(r)` is the density of the unheated tracer.

   For first-order heating, in which the specific energy injected, :math:`\Delta\epsilon = Q r^2`, is the same for all
   species, this is exact if the binding energy of each shell is set by the reference distribution. This is the tracer limit:
   the self-gravity of the tracer population is accounted for only insofar as it is included in the reference distribution.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalDecorator) :: massDistributionSphericalHeatedTracer
     !!{RST
     Implementation of a spherical mass distribution of a tracer population heated along with a reference mass distribution.
     !!}
     private
     class           (massDistributionSpherical), pointer :: massDistributionReference_ => null(), massDistributionReferenceUnheated_ => null()
     logical                                              :: isUnheated
     double precision                                     :: radiusFinalPrevious                 , radiusInitialPrevious
   contains
     !![
     <methods docformat="rst">
       <method method="radiusInitial" description="Compute the initial radius of the shell of tracer material which is at the given final radius after heating."/>
       <method method="radiusFinal"   description="Compute the final radius after heating of the shell of tracer material initially at the given radius."      />
     </methods>
     !!]
     final     ::                         sphericalHeatedTracerDestructor
     procedure :: radiusInitial        => sphericalHeatedTracerRadiusInitial
     procedure :: radiusFinal          => sphericalHeatedTracerRadiusFinal
     procedure :: density              => sphericalHeatedTracerDensity
     procedure :: massTotal            => sphericalHeatedTracerMassTotal
     procedure :: massEnclosedBySphere => sphericalHeatedTracerMassEnclosedBySphere
     procedure :: radiusEnclosingMass  => sphericalHeatedTracerRadiusEnclosingMass
  end type massDistributionSphericalHeatedTracer

  interface massDistributionSphericalHeatedTracer
     !!{RST
     Constructors for the :galacticus-class:`massDistributionSphericalHeatedTracer` mass distribution class.
     !!}
     module procedure sphericalHeatedTracerConstructorParameters
     module procedure sphericalHeatedTracerConstructorInternal
  end interface massDistributionSphericalHeatedTracer

contains

  function sphericalHeatedTracerConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`massDistributionSphericalHeatedTracer` mass distribution class which builds the
    object from a parameter set.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (massDistributionSphericalHeatedTracer)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(massDistributionClass                ), pointer       :: massDistribution_, massDistributionReference_

    !![
    <objectBuilder class="massDistribution"                                              name="massDistribution_"          source="parameters"/>
    <objectBuilder class="massDistribution" parameterName="massDistributionReference"    name="massDistributionReference_" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       select type (massDistributionReference_)
       class is (massDistributionSpherical)
          self=massDistributionSphericalHeatedTracer(massDistribution_,massDistributionReference_)
       class default
          call Error_Report('a spherically-symmetric reference mass distribution is required'//{introspection:location})
       end select
    class default
       call    Error_Report('a spherically-symmetric mass distribution is required'          //{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"         />
    <objectDestructor name="massDistributionReference_"/>
    !!]
    return
  end function sphericalHeatedTracerConstructorParameters

  function sphericalHeatedTracerConstructorInternal(massDistribution_,massDistributionReference_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`massDistributionSphericalHeatedTracer` mass distribution class. The
    reference distribution, ``massDistributionReference_``, must be a heated distribution (either
    :galacticus-class:`massDistributionSphericalHeated` or :galacticus-class:`massDistributionSphericalHeatedMonotonic`), from
    which the unheated reference distribution is found.
    !!}
    use :: Display, only : displayGreen, displayReset
    use :: Error  , only : Error_Report
    implicit none
    type (massDistributionSphericalHeatedTracer)                        :: self
    class(massDistributionSpherical            ), intent(in   ), target :: massDistribution_, massDistributionReference_
    !![
    <constructorAssign variables="*massDistribution_, *massDistributionReference_"/>
    !!]

    ! Find the unheated reference distribution, and determine whether any heating has occurred.
    select type (massDistributionReference_)
    class is (massDistributionSphericalHeated         )
       self%massDistributionReferenceUnheated_ => massDistributionReference_%massDistribution_
       self%isUnheated                         =  massDistributionReference_%massDistributionHeating_%specificEnergyIsEverywhereZero()
    class is (massDistributionSphericalHeatedMonotonic)
       self%massDistributionReferenceUnheated_ => massDistributionReference_%massDistribution_
       self%isUnheated                         =  massDistributionReference_%massDistributionHeating_%specificEnergyIsEverywhereZero()
    class default
       call Error_Report(                                                                                                                                                                     &
            &            'the reference mass distribution must be a heated mass distribution'                                                                                    //char(10)// &
            &            displayGreen()//'  HELP:'//displayReset()//' heating of a tracer population requires that its reference (e.g. dark matter) mass distribution is heated,'//char(10)// &
            &            '        e.g. by using a `heated` or `heatedMonotonic` dark matter profile'                                                                                       // &
            &            {introspection:location}                                                                                                                                             &
            &           )
    end select
    !![
    <referenceCountIncrement object="self%massDistributionReferenceUnheated_"/>
    !!]
    self%nonAnalyticSolver    =  nonAnalyticSolversNumerical
    self%radiusFinalPrevious  = -huge(0.0d0)
    self%radiusInitialPrevious= -huge(0.0d0)
    self%componentType        =  self%massDistribution_%componentType
    self%massType             =  self%massDistribution_%massType
    self%dimensionless        = .false.
    return
  end function sphericalHeatedTracerConstructorInternal

  subroutine sphericalHeatedTracerDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`massDistributionSphericalHeatedTracer` mass distribution class.
    !!}
    implicit none
    type(massDistributionSphericalHeatedTracer), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"                 />
    <objectDestructor name="self%massDistributionReference_"        />
    <objectDestructor name="self%massDistributionReferenceUnheated_"/>
    !!]
    return
  end subroutine sphericalHeatedTracerDestructor

  double precision function sphericalHeatedTracerRadiusInitial(self,radiusFinal) result(radiusInitial)
    !!{RST
    Compute the initial radius of the shell of tracer material which is at radius ``radiusFinal`` after heating. Shells which
    heating has unbound from the reference distribution are moved to infinity. To allow for this, all radii at or beyond
    ``radiusLarge`` are treated as being at infinity, and the initial radius is returned as ``radiusLarge`` for them, such that
    all tracer material is enclosed within that radius.
    !!}
    use :: Galactic_Structure_Options, only : radiusLarge
    implicit none
    class           (massDistributionSphericalHeatedTracer), intent(inout) :: self
    double precision                                       , intent(in   ) :: radiusFinal

    if (self%isUnheated) then
       radiusInitial=radiusFinal
    else if (radiusFinal <= 0.0d0) then
       radiusInitial=0.0d0
    else if (radiusFinal >= radiusLarge) then
       radiusInitial=radiusLarge
    else
       if (radiusFinal /= self%radiusFinalPrevious) then
          self%radiusFinalPrevious=radiusFinal
          select type (massDistributionReference_ => self%massDistributionReference_)
          class is (massDistributionSphericalHeated)
             ! The heated distribution provides the initial radius directly.
             self%radiusInitialPrevious=massDistributionReference_%radiusInitial(radiusFinal)
          class default
             ! Find the radius in the unheated reference distribution which encloses the same mass.
             self%radiusInitialPrevious=self%massDistributionReferenceUnheated_%radiusEnclosingMass(mass=self%massDistributionReference_%massEnclosedBySphere(radiusFinal))
          end select
       end if
       radiusInitial=self%radiusInitialPrevious
    end if
    return
  end function sphericalHeatedTracerRadiusInitial

  double precision function sphericalHeatedTracerRadiusFinal(self,radiusInitial) result(radiusFinal)
    !!{RST
    Compute the final radius after heating of the shell of tracer material initially at radius ``radiusInitial``.
    !!}
    implicit none
    class           (massDistributionSphericalHeatedTracer), intent(inout) :: self
    double precision                                       , intent(in   ) :: radiusInitial

    if (self%isUnheated) then
       radiusFinal=radiusInitial
    else if (radiusInitial <= 0.0d0) then
       radiusFinal=0.0d0
    else
       radiusFinal=self%massDistributionReference_%radiusEnclosingMass(mass=self%massDistributionReferenceUnheated_%massEnclosedBySphere(radiusInitial))
    end if
    return
  end function sphericalHeatedTracerRadiusFinal

  double precision function sphericalHeatedTracerDensity(self,coordinates) result(density)
    !!{RST
    Return the density at the specified ``coordinates`` in a heated tracer mass distribution.
    !!}
    use :: Coordinates               , only : coordinateSpherical, assignment(=)
    use :: Galactic_Structure_Options, only : radiusLarge
    implicit none
    class           (massDistributionSphericalHeatedTracer), intent(inout) :: self
    class           (coordinate                           ), intent(in   ) :: coordinates
    type            (coordinateSpherical                  )                :: coordinatesInitial
    double precision                                                       :: radiusInitial     , densityReferenceInitial

    if (self%isUnheated) then
       density=self%massDistribution_%density(coordinates)
       return
    end if
    radiusInitial=self%radiusInitial(coordinates%rSpherical())
    if (radiusInitial <= 0.0d0 .or. radiusInitial >= radiusLarge) then
       ! At zero radius the density is unchanged. Beyond the extent of the reference distribution there is no tracer material.
       if (radiusInitial <= 0.0d0) then
          density=self%massDistribution_%density(coordinates)
       else
          density=0.0d0
       end if
       return
    end if
    coordinatesInitial     =[radiusInitial,0.0d0,0.0d0]
    densityReferenceInitial=self%massDistributionReferenceUnheated_%density(coordinatesInitial)
    if (densityReferenceInitial > 0.0d0) then
       density=+self%massDistribution_         %density(coordinatesInitial) &
            &  *self%massDistributionReference_%density(coordinates       ) &
            &  /                                densityReferenceInitial
    else
       density=+0.0d0
    end if
    return
  end function sphericalHeatedTracerDensity

  double precision function sphericalHeatedTracerMassTotal(self) result(mass)
    !!{RST
    Return the total mass in a heated tracer mass distribution. Heating conserves mass, so this is the total mass of the
    unheated tracer distribution.
    !!}
    implicit none
    class(massDistributionSphericalHeatedTracer), intent(inout) :: self

    mass=self%massDistribution_%massTotal()
    return
  end function sphericalHeatedTracerMassTotal

  double precision function sphericalHeatedTracerMassEnclosedBySphere(self,radius) result(mass)
    !!{RST
    Computes the mass enclosed within a sphere of given ``radius`` for a heated tracer mass distribution.
    !!}
    use :: Galactic_Structure_Options, only : radiusLarge
    implicit none
    class           (massDistributionSphericalHeatedTracer), intent(inout), target :: self
    double precision                                       , intent(in   )         :: radius
    double precision                                                               :: radiusInitial

    radiusInitial=self%radiusInitial(radius)
    if (radiusInitial >= radiusLarge) then
       mass=self%massDistribution_%massTotal           (             )
    else
       mass=self%massDistribution_%massEnclosedBySphere(radiusInitial)
    end if
    return
  end function sphericalHeatedTracerMassEnclosedBySphere

  double precision function sphericalHeatedTracerRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{RST
    Computes the radius enclosing a given mass or mass fraction for a heated tracer mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalHeatedTracer), intent(inout), target   :: self
    double precision                                       , intent(in   ), optional :: mass, massFractional

    radius=self%radiusFinal(self%massDistribution_%radiusEnclosingMass(mass,massFractional))
    return
  end function sphericalHeatedTracerRadiusEnclosingMass
