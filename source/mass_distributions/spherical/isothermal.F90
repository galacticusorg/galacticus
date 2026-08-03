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

  !!{RST
  Implementation of an isothermal mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionIsothermal" docformat="rst">
    <description>
    An isothermal mass distribution class in which the density profile is given by:

    .. math::

       \rho(r) \propto r^{-2}.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionIsothermal
     !!{RST
     The isothermal mass distribution.
     !!}
     private
     double precision :: densityNormalization, lengthReference, &
          &              velocityRotation
   contains
     procedure :: massTotal                         => isothermalMassTotal
     procedure :: density                           => isothermalDensity
     procedure :: densityGradientRadial             => isothermalDensityGradientRadial
     procedure :: densityRadialMoment               => isothermalDensityRadialMoment
     procedure :: massEnclosedBySphere              => isothermalMassEnclosedBySphere
     procedure :: rotationCurve                     => isothermalRotationCurve
     procedure :: rotationCurveGradient             => isothermalRotationCurveGradient
     procedure :: velocityRotationCurveMaximum      => isothermalVelocityRotationCurveMaximum
     procedure :: radiusRotationCurveMaximum        => isothermalRadiusRotationCurveMaximum
     procedure :: radiusEnclosingMass               => isothermalRadiusEnclosingMass
     procedure :: radiusEnclosingDensity            => isothermalRadiusEnclosingDensity
     procedure :: radiusFromSpecificAngularMomentum => isothermalRadiusFromSpecificAngularMomentum
     procedure :: fourierTransform                  => isothermalFourierTransform
     procedure :: radiusFreefall                    => isothermalRadiusFreefall
     procedure :: radiusFreefallIncreaseRate        => isothermalRadiusFreefallIncreaseRate 
     procedure :: energyPotential                   => isothermalEnergyPotential
     procedure :: energyKinetic                     => isothermalEnergyKinetic
     procedure :: potentialIsAnalytic               => isothermalPotentialIsAnalytic
     procedure :: potential                         => isothermalPotential
     procedure :: positionSample                    => isothermalPositionSample
     procedure :: descriptor                        => isothermalDescriptor
  end type massDistributionIsothermal

  interface massDistributionIsothermal
     !!{RST
     Constructors for the :galacticus-class:`massDistributionIsothermal` mass distribution class.
     !!}
     module procedure massDistributionIsothermalConstructorParameters
     module procedure massDistributionIsothermalConstructorInternal
  end interface massDistributionIsothermal

contains

  function massDistributionIsothermalConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`massDistributionIsothermal` mass distribution class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    type            (massDistributionIsothermal)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    double precision                                            :: mass                , lengthReference, &
         &                                                         densityNormalization
    logical                                                     :: dimensionless
    type            (varying_string            )                :: componentType
    type            (varying_string            )                :: massType

    !![
    <inputParameter docformat="rst">
      <name>densityNormalization</name>
      <defaultValue>0.25d0/Pi</defaultValue>
      <description>
      The density normalization of the isothermal profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>lengthReference</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The reference radius (in Mpc) of the isothermal profile at which ``densityNormalization`` is evaluated; sets the physical length scale of the :math:`\rho \propto r^{-2}` profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The total mass (in :math:`\mathrm{M}_\odot`) enclosed within the reference radius, used to normalize the :math:`\rho \propto r^{-2}` density profile when ``densityNormalization`` is not supplied.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true the isothermal profile is considered to be dimensionless.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>
      The component type that this mass distribution represents.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>
      The mass type that this mass distribution represents.
      </description>
      <source>parameters</source>
    </inputParameter>
    <conditionalCall>
     <call>self=massDistributionIsothermal(componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="densityNormalization" value="densityNormalization" parameterPresent="parameters"/>
     <argument name="mass"                 value="mass"                 parameterPresent="parameters"/>
     <argument name="lengthReference"      value="lengthReference"      parameterPresent="parameters"/>
     <argument name="dimensionless"        value="dimensionless"        parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massDistributionIsothermalConstructorParameters

  function massDistributionIsothermalConstructorInternal(densityNormalization,mass,lengthReference,dimensionless,componentType,massType) result(self)
    !!{RST
    Internal constructor for "isothermal" mass distribution class.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Comparison            , only : Values_Differ
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    type            (massDistributionIsothermal   )                          :: self
    double precision                               , intent(in   ), optional :: densityNormalization, mass, &
         &                                                                      lengthReference
    logical                                        , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType      ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="componentType, massType"/>
    !!]

    ! Determine if profile is dimensionless.
    self%dimensionless=.false.
    if (present(dimensionless)) self%dimensionless=dimensionless
    ! If dimensionless, then set scale length and mass to unity.
    if (self%dimensionless) then
       if (present(lengthReference     )) then
          if (Values_Differ(lengthReference     ,1.0d0    ,absTol=1.0d-6)) call Error_Report('scaleLength should be unity for a dimensionless profile (or simply do not specify a scale length)'               //{introspection:location})
       end if
       if (present(mass                )) then
          if (Values_Differ(mass                ,1.0d0    ,absTol=1.0d-6)) call Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'                              //{introspection:location})
       end if
       if (present(densityNormalization)) then
          if (Values_Differ(densityNormalization,0.25d0/Pi,absTol=1.0d-6)) call Error_Report('densityNormalization should be π/4 for a dimensionless profile (or simply do not specify a densityNormalization)'//{introspection:location})
       end if
       self%lengthReference     =1.00d0
       self%densityNormalization=0.25d0/Pi
    else
       if      (present(lengthReference     )) then
          self%lengthReference=lengthReference
       else
          call Error_Report('"lengthReference" must be specified'//{introspection:location})
       end if
       if      (present(densityNormalization)) then
          self%densityNormalization=densityNormalization
       else if (present(mass                )) then
          self%densityNormalization=mass                /4.0d0/Pi/lengthReference**3
       else
          call Error_Report('one of "densityNormalization" or "mass" must be specified'//{introspection:location})
       end if
    end if
    ! Compute the rotation velocity.
    if (self%isDimensionless()) then
       self%velocityRotation=+1.0d0
    else
       self%velocityRotation=+sqrt(                                   &
            &                      +4.0d0                             &
            &                      *Pi                                &
            &                      *gravitationalConstant_internal    &
            &                      *self%lengthReference          **2 &
            &                      *self%densityNormalization         &
            &                     )
    end if
    return
  end function massDistributionIsothermalConstructorInternal

  double precision function isothermalMassTotal(self)
    !!{RST
    Return the total mass in an isothermal mass distribution.
    !!}
    implicit none
    class(massDistributionIsothermal), intent(inout) :: self

    isothermalMassTotal=huge(0.0d0)
    return
  end function isothermalMassTotal

  double precision function isothermalDensity(self,coordinates)
    !!{RST
    Return the density at the specified ``coordinates`` in an isothermal mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class(massDistributionIsothermal), intent(inout) :: self
    class(coordinate                ), intent(in   ) :: coordinates

    isothermalDensity=+ self       %densityNormalization   &
         &            /(                                   &
         &             +coordinates%rSpherical          () &
         &             /self       %lengthReference        &
         &            )**2
    return
  end function isothermalDensity

  double precision function isothermalDensityGradientRadial(self,coordinates,logarithmic)
    !!{RST
    Return the density at the specified ``coordinates`` in an isothermal mass distribution.
    !!}
    implicit none
    class           (massDistributionIsothermal), intent(inout), target   :: self
    class           (coordinate                ), intent(in   )           :: coordinates
    logical                                     , intent(in   ), optional :: logarithmic
    double precision                                                      :: radius
    logical                                                               :: logarithmicActual

    ! Set default options.
    logarithmicActual=.false.
    if (present(logarithmic)) logarithmicActual=logarithmic
    ! Get position in spherical coordinate system.
    radius=coordinates%rSpherical()
    ! Compute density gradient.
    if (logarithmicActual) then
       isothermalDensityGradientRadial=-2.0d0
    else
       isothermalDensityGradientRadial=-2.0d0                        &
            &                          *self%densityNormalization    &
            &                          *self%lengthReference     **2 &
            &                          /radius                   **3
    end if
    return
  end function isothermalDensityGradientRadial

  double precision function isothermalMassEnclosedBySphere(self,radius)
    !!{RST
    Computes the mass enclosed within a sphere of given ``radius`` for isothermal mass distributions.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionIsothermal), intent(inout), target :: self
    double precision                            , intent(in   )         :: radius

    isothermalMassEnclosedBySphere=+4.0d0                        &
         &                         *Pi                           &
         &                         *self%densityNormalization    &
         &                         *self%lengthReference     **2 &
         &                         *radius
    return
  end function isothermalMassEnclosedBySphere

  double precision function isothermalRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{RST
    Computes the radius enclosing a given mass or mass fraction for isothermal mass distributions.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionIsothermal), intent(inout), target   :: self
    double precision                            , intent(in   ), optional :: mass , massFractional
    double precision                                                      :: mass_

    mass_=0.0d0
    if (present(mass)) then
       mass_=mass
    else if (present(massFractional)) then
       call Error_Report('mass is unbounded, so mass fraction is undefined'//{introspection:location})
    else
       call Error_Report('either mass or massFractional must be supplied'//{introspection:location})
    end if
    radius=+     mass_                   &
         & /     4.0d0                   &
         & /     Pi                      &
         & /self%densityNormalization    &
         & /self%lengthReference     **2
    return
  end function isothermalRadiusEnclosingMass
  
  double precision function isothermalRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{RST
    Computes the radius enclosing a given mean density for isothermal mass distributions.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionIsothermal), intent(inout), target   :: self
    double precision                            , intent(in   )           :: density
    double precision                            , intent(in   ), optional :: radiusGuess

    radius=+      self%lengthReference      &
         & *sqrt(                           &
         &       +3.0d0                     &
         &       *self%densityNormalization &
         &       /     density              &
         &      )
    return
  end function isothermalRadiusEnclosingDensity
  
  double precision function isothermalRadiusFromSpecificAngularMomentum(self,angularMomentumSpecific) result(radius)
    !!{RST
    Computes the radius corresponding to a given specific angular momentum for isothermal mass distributions.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionIsothermal), intent(inout), target :: self
    double precision                            , intent(in   )         :: angularMomentumSpecific

    radius=+angularMomentumSpecific         &
         & /sqrt(                           &
         &       +4.0d0                     &
         &       *Pi                        &
         &       *self%densityNormalization &
         &       *self%lengthReference**2   &
         &      )
    if (.not.self%isDimensionless()) radius=+radius                               &
         &                                  /sqrt(gravitationalConstant_internal)
    return
  end function isothermalRadiusFromSpecificAngularMomentum
  
  double precision function isothermalDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{RST
    Returns a radial density moment for the Isothermal mass distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionIsothermal), intent(inout)           :: self
    double precision                            , intent(in   )           :: moment
    double precision                            , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                     , intent(  out), optional :: isInfinite
    double precision                                                      :: momentMinimum, momentMaximum

    isothermalDensityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    if (present(radiusMinimum)) then
       if (moment == 1.0d0) then
          momentMinimum=+log(+radiusMinimum                )
       else
          momentMinimum=     +radiusMinimum**(moment-1.0d0) &
               &        /                    (moment-1.0d0)
       end if
    else if (moment <= 1.0d0) then
       momentMinimum=+0.0d0
       if (present(isInfinite)) then
          isInfinite=.true.
          return
       else
          call Error_Report('radial moment is infinite'//{introspection:location})
       end if
    else
       momentMinimum=0.0d0
    end if
    if (present(radiusMaximum)) then
       if (moment == 1.0d0) then
          momentMaximum=+log(+radiusMaximum                )
       else
          momentMaximum=     +radiusMaximum**(moment-1.0d0) &
               &        /                    (moment-1.0d0)
       end if
    else if (moment >= 1.0d0) then
       momentMaximum=+0.0d0
       if (present(isInfinite)) then
          isInfinite=.true.
          return
       else
          call Error_Report('radial moment is infinite'//{introspection:location})
       end if
    else
       momentMaximum=0.0d0
    end if
    isothermalDensityRadialMoment=+self%densityNormalization    &
         &                        *self%lengthReference     **2 &
         &                        *(                            &
         &                          +momentMaximum              &
         &                          -momentMinimum              &
         &                         )
    return
  end function isothermalDensityRadialMoment

  double precision function isothermalRotationCurve(self,radius) result(rotationCurve)
    !!{RST
    Return the rotation curve for an isothermal mass distribution.
    !!}
    implicit none
    class           (massDistributionIsothermal), intent(inout) :: self
    double precision                            , intent(in   ) :: radius

    rotationCurve=self%velocityRotation
    return
  end function isothermalRotationCurve

  double precision function isothermalRotationCurveGradient(self,radius) result(rotationCurveGradient)
    !!{RST
    Return the rotation curve gradient (specifically, :math:`\mathrm{d}V^2_\mathrm{c}/\mathrm{d}r`) for an isothermal mass distribution.
    !!}
    implicit none
    class           (massDistributionIsothermal), intent(inout) :: self
    double precision                            , intent(in   ) :: radius
    !$GLC attributes unused :: self, radius
    
    rotationCurveGradient=0.0d0
    return
  end function isothermalRotationCurveGradient

  double precision function isothermalVelocityRotationCurveMaximum(self) result(velocity)
    !!{RST
    Return the peak velocity in the rotation curve for an isothermal mass distribution.
    !!}
    implicit none
    class(massDistributionIsothermal), intent(inout) :: self

    velocity=self%velocityRotation
    return
  end function isothermalVelocityRotationCurveMaximum

  double precision function isothermalRadiusRotationCurveMaximum(self) result(radius)
    !!{RST
    Return the peak velocity in the rotation curve for an isothermal mass distribution.
    !!}
    implicit none
    class(massDistributionIsothermal), intent(inout), target :: self
    !$GLC attributes unused :: self

    radius=1.0d0
    return
  end function isothermalRadiusRotationCurveMaximum

  logical function isothermalPotentialIsAnalytic(self) result(isAnalytic)
    !!{RST
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionIsothermal), intent(inout) :: self

    isAnalytic=.true.
    return
  end function isothermalPotentialIsAnalytic

  double precision function isothermalPotential(self,coordinates,status)
    !!{RST
    Return the potential at the specified ``coordinates`` in an isothermal mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess     , structureErrorCodeInfinite
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    use :: Error                           , only : Error_Report
    implicit none
    class(massDistributionIsothermal       ), intent(inout), target   :: self
    class(coordinate                       ), intent(in   )           :: coordinates
    type (enumerationStructureErrorCodeType), intent(  out), optional :: status

    if (present(status)) status=structureErrorCodeSuccess
    ! Compute the potential at this position.
    if (coordinates%rSpherical() <= 0.0d0) then
       isothermalPotential=0.0d0
       if (present(status)) then
          status=structureErrorCodeInfinite
          return
       else
          call Error_Report('potential is divergent at zero radius'//{introspection:location})
       end if
    end if
    isothermalPotential=+     4.0d0                                 &
         &              *Pi                                         &
         &              *     self       %densityNormalization      &
         &              *     self       %lengthReference       **2 &
         &              *log(                                       &
         &                   +coordinates%rSpherical          ()    &
         &                   /self       %lengthReference           &
         &                  )
    if (.not.self%isDimensionless()) isothermalPotential=+gravitationalConstant_internal &
         &                                               *isothermalPotential
    return
  end function isothermalPotential

  double precision function isothermalFourierTransform(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{RST
    Compute the Fourier transform of the density profile at the given ``wavenumber`` in an isothermal mass distribution.
    !!}
    use :: Exponential_Integrals, only : Sine_Integral
    implicit none
    class           (massDistributionIsothermal), intent(inout) :: self
    double precision                            , intent(in   ) :: radiusOuter        , wavenumber
    double precision                                            :: wavenumberScaleFree

    waveNumberScaleFree=+waveNumber  &
         &              *radiusOuter
    fourierTransform   =+Sine_Integral(waveNumberScaleFree) &
         &              /              waveNumberScaleFree
    return
  end function isothermalFourierTransform
  
  double precision function isothermalRadiusFreefall(self,time) result(radius)
    !!{RST
    Compute the freefall radius at the given ``time`` in an isothermal mass distribution. For an isothermal potential, the freefall radius, :math:`r_\mathrm{ff}(t)`, is:

    .. math::

       r_\mathrm{ff}(t) = \sqrt{{2 \over \pi}} V_\mathrm{virial} t.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    class           (massDistributionIsothermal), intent(inout) :: self
    double precision                            , intent(in   ) :: time

    radius=+sqrt(                   &
         &       +2.0d0             &
         &       /Pi                &
         &      )                   &
         & *self%velocityRotation   &
         & *     time               &
         & /MpcPerKmPerSToGyr
    return
  end function isothermalRadiusFreefall
  
  double precision function isothermalRadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{RST
    Compute the rate of increase of the freefall radius at the given ``time`` in an isothermal mass distribution. For an isothermal potential, the rate of increase of the freefall radius, :math:`\dot{r}_\mathrm{ff}(t)`, is:

    .. math::

       \dot{r}_\mathrm{ff}(t) = \sqrt{{2 \over \pi}} V_\mathrm{virial}.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    class           (massDistributionIsothermal), intent(inout) :: self
    double precision                            , intent(in   ) :: time
    !$GLC attributes unused :: time

    radiusIncreaseRate=+sqrt(                   &
         &                   +2.0d0             &
         &                   /Pi                &
         &                  )                   &
         &             *self%velocityRotation   &
         &             /MpcPerKmPerSToGyr
    return
  end function isothermalRadiusFreefallIncreaseRate
  
  double precision function isothermalEnergyPotential(self,radiusOuter) result(energy)
    !!{RST
    Compute the potential energy within a given ``radius`` in an isothermal mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionIsothermal), intent(inout) :: self
    double precision                            , intent(in   ) :: radiusOuter

    energy=-16.0d0                                 &
         & *Pi                                 **2 &
         & *     gravitationalConstant_internal    &
         & *self%lengthReference               **4 &
         & *self%densityNormalization          **2 &
         & *     radiusOuter
    return
  end function isothermalEnergyPotential

  double precision function isothermalEnergyKinetic(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{RST
    Compute the kinetic energy within a given ``radius`` in an isothermal mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class           (massDistributionIsothermal), intent(inout) :: self
    double precision                            , intent(in   ) :: radiusOuter
    class           (massDistributionClass     ), intent(inout) :: massDistributionEmbedding
    logical                                                     :: analytic
    type            (coordinateSpherical       )                :: coordinates

    analytic=.false.
    select type (massDistributionEmbedding)
    class is (massDistributionIsothermal)
       select type (kinematicsDistribution_ => massDistributionEmbedding%kinematicsDistribution_)
       class is (kinematicsDistributionIsothermal)
          analytic   =.true.
          coordinates=[radiusOuter,0.0d0,0.0d0]
          energy     =+1.5d0                                                                                    &
               &      *self                   %massEnclosedBySphere(radiusOuter                               ) &
               &      *kinematicsDistribution_%velocityDispersion1D(coordinates,self,massDistributionEmbedding)
       end select
    end select
    if (.not.analytic) energy=self%energyKineticNumerical(radiusOuter,massDistributionEmbedding)
    return
  end function isothermalEnergyKinetic
  
  function isothermalPositionSample(self,randomNumberGenerator_) result(position)
    !!{RST
    Computes the half-mass radius of a spherically symmetric mass distribution using numerical root finding.
    !!}
    use :: Coordinates             , only : coordinateCartesian, assignment(=)
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (coordinateCartesian        )               :: position
    class           (massDistributionIsothermal), intent(inout) :: self
    class           (randomNumberGeneratorClass), intent(inout) :: randomNumberGenerator_

    position=[0.0d0,0.0d0,0.0d0]
    call Error_Report('can not sample positions, mass is unbounded'//{introspection:location})
    return
  end function isothermalPositionSample

  subroutine isothermalDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionIsothermal), intent(inout)           :: self
    type     (inputParameters           ), intent(inout)           :: descriptor
    logical                              , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                    )                          :: parameterLabel
    type     (inputParameters           )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','isothermal')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%lengthReference
    call parameters%addParameter('lengthReference'     ,trim(adjustl(parameterLabel)))
    return
  end subroutine isothermalDescriptor
