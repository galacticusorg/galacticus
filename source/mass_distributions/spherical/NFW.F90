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
  Implementation of an NFW :cite:p:`navarro_structure_1996` mass distribution class.
  !!}

  use :: Tabulations_Inverse, only : tabulationInverse

  public :: massDistributionNFWStateStore, massDistributionNFWStateRestore
  
  !![
  <massDistribution name="massDistributionNFW" docformat="rst">
    <description>
    An NFW :cite:p:`navarro_structure_1996` mass distribution class. The density profile is given by:

    .. math::

        \rho_\mathrm{dark matter}(r) \propto \left({r\over r_\mathrm{s}}\right)^{-1} \left[1 + \left({r\over r_\mathrm{s}}\right)
       \right]^{-2}.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionNFW
     !!{RST
     The NFW :cite:p:`navarro_structure_1996` mass distribution.
     !!}
     private
     double precision :: densityNormalization      , scaleLength
     double precision :: enclosedMassRadiusPrevious, enclosedMassPrevious
   contains
     !![
     <methods docformat="rst">
       <method method="initialize"           description="(Re)initialize the parameters of the NFW mass distribution."                              />
       <method method="timeFreefallTabulate" description="Tabulate the freefall time as a function of radius in a scale-free NFW mass distribution."/>
     </methods>
     !!]
     procedure :: initialize                        => nfwInitialize
     procedure :: massTotal                         => nfwMassTotal
     procedure :: density                           => nfwDensity
     procedure :: densityGradientRadial             => nfwDensityGradientRadial
     procedure :: densitySlopeLogarithmicCentral    => nfwDensitySlopeLogarithmicCentral
     procedure :: densityRadialMoment               => nfwDensityRadialMoment
     procedure :: massEnclosedBySphere              => nfwMassEnclosedBySphere
     procedure :: velocityRotationCurveMaximum      => nfwVelocityRotationCurveMaximum
     procedure :: radiusRotationCurveMaximum        => nfwRadiusRotationCurveMaximum
     procedure :: radiusEnclosingMass               => nfwRadiusEnclosingMass
     procedure :: radiusEnclosingDensity            => nfwRadiusEnclosingDensity
     procedure :: radiusFromSpecificAngularMomentum => nfwRadiusFromSpecificAngularMomentum
     procedure :: fourierTransform                  => nfwFourierTransform
     procedure :: radiusFreefall                    => nfwRadiusFreefall
     procedure :: radiusFreefallIncreaseRate        => nfwRadiusFreefallIncreaseRate
     procedure :: timeFreefallTabulate              => nfwTimeFreefallTabulate
     procedure :: potentialIsAnalytic               => nfwPotentialIsAnalytic
     procedure :: potential                         => nfwPotential
     procedure :: energyPotential                   => nfwEnergyPotential
     procedure :: energyKinetic                     => nfwEnergyKinetic
     procedure :: descriptor                        => nfwDescriptor
  end type massDistributionNFW
  
  interface massDistributionNFW
     !!{RST
     Constructors for the :galacticus-class:`massDistributionNFW` mass distribution class.
     !!}
     module procedure massDistributionNFWConstructorParameters
     module procedure massDistributionNFWConstructorInternal
  end interface massDistributionNFW

  ! Tabulated solutions. These are held per thread rather than per object: the scale-free NFW profile has no shape parameter,
  ! so each of these tabulations is universal and is shared by every NFW distribution a thread constructs. (The corresponding
  ! tabulations of the Einasto, Burkert and Zhao1996 profiles are held per object precisely because each of those profiles
  ! *does* have a shape parameter, which makes them specific to the object.)
  !
  ! The abscissae are points of an absolute octave lattice. The bounds these tabulations reach are found by halving and
  ! doubling from a seed of unity, so they are exact powers of two and are therefore already points of that lattice: pinning
  ! costs nothing in range here, where anchoring to whole decades would pad every bound outward. Thirty points per octave is
  ! 99.66 per decade, against the hundred per decade used before.
  integer                      , parameter    :: countRadiiPerOctave=30
  type            (tabulationInverse)         :: densityScaleFree_                , angularMomentumSpecificScaleFree_, &
       &                                         timeFreefallScaleFree_
  !$omp threadprivate(densityScaleFree_,angularMomentumSpecificScaleFree_,timeFreefallScaleFree_)

contains

  function massDistributionNFWConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`massDistributionNFW` mass distribution class which builds the object from a parameter set.
    !!}
     use :: Input_Parameters          , only : inputParameter                , inputParameters
     use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
     use :: Numerical_Constants_Math  , only : Pi
    implicit none
    type            (massDistributionNFW)                :: self
    type            (inputParameters    ), intent(inout) :: parameters
    double precision                                     :: mass                , scaleLength  , &
         &                                                  densityNormalization, concentration, &
         &                                                  virialRadius
    logical                                              :: dimensionless
    type            (varying_string     )                :: componentType
    type            (varying_string     )                :: massType

    !![
    <inputParameter docformat="rst">
      <name>densityNormalization</name>
      <defaultValue>1.0d0/2.0d0/Pi/(log(4.0d0)-1.0d0)</defaultValue>
      <description>
      The density normalization of the NFW profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>scaleLength</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The NFW scale radius (in Mpc) :math:`r_\mathrm{s}` at which the density profile transitions from the inner :math:`\rho \propto r^{-1}` to the outer :math:`\rho \propto r^{-3}` slope.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The total mass (in :math:`\mathrm{M}_\odot`) of the NFW profile, used to set the density normalization when the concentration and virial radius are provided.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>concentration</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The halo concentration parameter :math:`c = r_\mathrm{vir}/r_\mathrm{s}` of the NFW profile, controlling how centrally concentrated the dark matter density profile is.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>virialRadius</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The virial radius (in Mpc) :math:`r_\mathrm{vir}` of the NFW halo, which defines the outer boundary of the profile at which the mean enclosed density equals the virial overdensity threshold.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true the NFW profile is considered to be dimensionless.
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
     <call>self=massDistributionNFW(componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="densityNormalization" value="densityNormalization" parameterPresent="parameters"/>
     <argument name="mass"                 value="mass"                 parameterPresent="parameters"/>
     <argument name="scaleLength"          value="scaleLength"          parameterPresent="parameters"/>
     <argument name="virialRadius"         value="virialRadius"         parameterPresent="parameters"/>
     <argument name="concentration"        value="concentration"        parameterPresent="parameters"/>
     <argument name="dimensionless"        value="dimensionless"        parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massDistributionNFWConstructorParameters

  function massDistributionNFWConstructorInternal(scaleLength,concentration,densityNormalization,mass,virialRadius,dimensionless,componentType,massType) result(self)
    !!{RST
    Internal constructor for "nfw" mass distribution class.
    !!}
    implicit none
    type            (massDistributionNFW         )                          :: self
    double precision                              , intent(in   ), optional :: scaleLength         , concentration, &
         &                                                                     densityNormalization, mass         , &
         &                                                                     virialRadius
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="componentType, massType"/>
    !!]

    call self%initialize(scaleLength,concentration,densityNormalization,mass,virialRadius,dimensionless)
    return
  end function massDistributionNFWConstructorInternal

  subroutine nfwInitialize(self,scaleLength,concentration,densityNormalization,mass,virialRadius,dimensionless)
    !!{RST
    Initialize the parameters of an NFW mass distribution.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionNFW), intent(inout)           :: self
    double precision                     , intent(in   ), optional :: scaleLength         , concentration, &
         &                                                            densityNormalization, mass         , &
         &                                                            virialRadius
    logical                              , intent(in   ), optional :: dimensionless
    double precision                                               :: radiusScaleFree

    ! Determine scale length
    if      (                            &
         &   present(scaleLength  )      &
         &  ) then
       self%scaleLength=scaleLength
    else if (                            &
         &   present(concentration).and. &
         &   present(virialRadius )      &
         &  ) then
       self%scaleLength=virialRadius/concentration
    else
       self%scaleLength=0.0d0
       call Error_Report('no means to determine scale length'//{introspection:location})
    end if
    ! Determine density normalization.
    if      (                                   &
         &   present(densityNormalization)      &
         &  ) then
       self%densityNormalization=densityNormalization
    else if (                                   &
         &   present(mass                ).and. &
         &   present(virialRadius        )      &
         &  ) then
       radiusScaleFree          =+virialRadius/self%scaleLength
       self%densityNormalization=+mass/4.0d0/Pi/self%scaleLength**3/(log(1.0d0+radiusScaleFree)-radiusScaleFree/(1.0d0+radiusScaleFree))
    else
       self%densityNormalization=+0.00
       call Error_Report('either "densityNormalization", or "mass" and "virialRadius" must be specified'//{introspection:location})
    end if
    ! Determine if profile is dimensionless.
    if      (present(dimensionless     )) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! (Re)initialize memoized results.
    self%enclosedMassPrevious      =-huge(0.0d0)
    self%enclosedMassRadiusPrevious=-huge(0.0d0)
    ! Defensively clear memoized state inherited from the spherical base class, and
    ! from any attached kinematics distribution, so that a re-initialized (pooled)
    ! object never serves stale tabulations built for a previous node. For NFW these
    ! paths are analytic and normally unused, but this guards against use of the
    ! numerical fallbacks (e.g. when the profile is embedded in a composite system).
    call self%tabulationReset()
    if (associated(self%kinematicsDistribution_)) call self%kinematicsDistribution_%reset()
    return
  end subroutine nfwInitialize

  double precision function nfwMassTotal(self)
    !!{RST
    Return the total mass in an NFW mass distribution.
    !!}
    implicit none
    class(massDistributionNFW), intent(inout) :: self

    nfwMassTotal=huge(0.0d0)
    return
  end function nfwMassTotal

  double precision function nfwDensity(self,coordinates)
    !!{RST
    Return the density at the specified ``coordinates`` in an NFW mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW), intent(inout) :: self
    class           (coordinate         ), intent(in   ) :: coordinates
    double precision                                     :: radiusScaleFree

    ! Compute the density at this position.
    radiusScaleFree=+coordinates%rSpherical          () &
         &          /self       %scaleLength
    nfwDensity     =+self       %densityNormalization   &
         &          /       radiusScaleFree             &
         &          /(1.0d0+radiusScaleFree)**2
    return
  end function nfwDensity

  double precision function nfwDensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{RST
    Return the density at the specified ``coordinates`` in an NFW :cite:p:`navarro_structure_1996` mass distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionNFW), intent(inout), target   :: self
    class           (coordinate         ), intent(in   )           :: coordinates
    logical                              , intent(in   ), optional :: logarithmic
    double precision                                               :: radiusScaleFree
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    densityGradientRadial=0.0d0
    radiusScaleFree      =+coordinates%rSpherical()         &
         &                /self       %scaleLength
    if (radiusScaleFree <= 0.0d0) then
       if (logarithmic_) then
          densityGradientRadial=-1.0d0
       else
          call Error_Report('gradient is divergent at r=0'//{introspection:location})
       end if
    else if (logarithmic_) then
       ! Closed form for the logarithmic slope: dln ρ/dln r = -(1+3x)/(1+x) with x=r/r_s.
       densityGradientRadial=-(1.0d0+3.0d0*radiusScaleFree)   &
            &                /(1.0d0+      radiusScaleFree)
    else
       densityGradientRadial=-self       %densityNormalization &
            &                /self       %scaleLength          &
            &                /             radiusScaleFree **2 &
            &                *(1.0d0+3.0d0*radiusScaleFree)    &
            &                /(1.0d0+      radiusScaleFree)**3
    end if
    return
  end function nfwDensityGradientRadial

  double precision function nfwDensitySlopeLogarithmicCentral(self) result(slope)
    !!{RST
    Return the central logarithmic slope of the density profile in an NFW :cite:p:`navarro_structure_1996` mass
    distribution. The profile is cuspy, so the density is not evaluable at zero radius.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class(massDistributionNFW ), intent(inout) :: self
    type (coordinateSpherical)                 :: coordinates

    coordinates=[0.0d0,0.0d0,0.0d0]
    slope      =self%densityGradientRadial(coordinates,logarithmic=.true.)
    return
  end function nfwDensitySlopeLogarithmicCentral

  double precision function nfwDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{RST
    Computes radial moments of the density in an NFW :cite:p:`navarro_structure_1996` mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW), intent(inout)           :: self
    double precision                     , intent(in   )           :: moment
    double precision                     , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    logical                              , intent(  out), optional :: isInfinite
    double precision                                               :: radialMomentMinimum, radialMomentMaximum

    densityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    if (present(radiusMinimum)) then
       radialMomentMinimum=radialMomentScaleFree(radiusMinimum/self%scaleLength)
    else
       radialMomentMinimum=radialMomentScaleFree(                         0.0d0)
    end if
    if (present(radiusMaximum)) then
       radialMomentMaximum=radialMomentScaleFree(radiusMaximum/self%scaleLength)
    else
       radialMomentMaximum=0.0d0
       if (moment >= 3.0d0) then
          if (present(isInfinite)) then
             isInfinite=.true.
             return
          else
             call Error_Report('moment is infinite'//{introspection:location})
          end if
       end if
    end if
    densityRadialMoment=+self%densityNormalization                 &
         &              *self%scaleLength         **(moment+1.0d0) &
         &              *(                                         &
         &                +radialMomentMaximum                     &
         &                -radialMomentMinimum                     &
         &               )    
    return

  contains

    double precision function radialMomentScaleFree(radius)
      !!{RST
      Provides the scale-free part of the radial moment of the NFW density profile.
      !!}
      use :: Hypergeometric_Functions, only : Hypergeometric_2F1
      use :: Numerical_Comparison    , only : Values_Agree
      implicit none
      double precision, intent(in   ) :: radius

      if (Values_Agree(moment,0.0d0,absTol=1.0d-6)) then
         ! Take the real part of this improper integral. The imaginary parts must cancel when taking differences to compute a
         ! proper integral.
         radialMomentScaleFree=+1.0d0/                 (1.0d0+      radius        ) &
              &                -2.0d0*real(atanh(dcmplx(1.0d0+2.0d0*radius,0.0d0)))
      else if (Values_Agree(moment,1.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=-1.0d0/                 (1.0d0      +radius        )
      else if (Values_Agree(moment,2.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=+1.0d0/                 (1.0d0      +radius        ) &
              &                +      log              (1.0d0      +radius        )
      else if (Values_Agree(moment,3.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=+                                    radius          &
              &                -1.0d0/                 (1.0d0      +radius        ) &
              &                -2.0d0*log              (1.0d0      +radius        )
      else
         radialMomentScaleFree=+(1.0d0+radius)**(moment-1.0d0)                                                     &
              &                /moment                                                                             &
              &                /                (moment-1.0d0)                                                     &
              &                *(                                                                                  &
              &                  - moment                                                                          &
              &                  *  Hypergeometric_2F1([1.0d0-moment,-moment],[2.0d0-moment],1.0d0/(1.0d0+radius)) &
              &                  +(1.0d0+radius)                                                                   &
              &                  *(moment-1.0d0)                                                                   &
              &                  *(                                                                                &
              &                    +(radius/(1.0d0+radius))**moment                                                &
              &                    -Hypergeometric_2F1([     -moment,-moment],[1.0d0-moment],1.0d0/(1.0d0+radius)) &
              &                  )                                                                                 &
              &                 )
      end if
      return
    end function radialMomentScaleFree

  end function nfwDensityRadialMoment

  double precision function nfwMassEnclosedBySphere(self,radius) result(mass)
    !!{RST
    Computes the mass enclosed within a sphere of given ``radius`` for nfw mass distributions.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionNFW ), intent(inout), target :: self
    double precision                      , intent(in   )         :: radius
    double precision                                              :: radiusScaleFree
    
    if (radius /= self%enclosedMassRadiusPrevious) then
       self%enclosedMassRadiusPrevious=+     radius
       radiusScaleFree                =+     radius                                    &
            &                          /self%scaleLength
       self%enclosedMassPrevious      =+     massEnclosedScaleFree(radiusScaleFree)    &
            &                          *self%densityNormalization                      &
            &                          *self%scaleLength                           **3
    end if
    mass=self%enclosedMassPrevious
    return
  end function nfwMassEnclosedBySphere
  
  double precision function nfwRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{RST
    Computes the radius enclosing a given mass or mass fraction for nfw mass distributions.
    !!}    
    use :: Numerical_Constants_Math, only : Pi
    use :: Lambert_Ws              , only : Lambert_W0
    use :: Error                   , only : Error_Report
    implicit none
    class           (massDistributionNFW ), intent(inout), target   :: self
    double precision                      , intent(in   ), optional :: mass                     , massFractional
    double precision                      , parameter               :: massScaleFreeSmall=3.0d-4
    double precision                                                :: mass_                    , massScaleFree

    mass_=0.0d0
    if (present(mass)) then
       mass_=mass
    else if (present(massFractional)) then
       call Error_Report('mass is unbounded, so mass fraction is undefined'//{introspection:location})
    else
       call Error_Report('either mass or massFractional must be supplied'  //{introspection:location})
    end if
    massScaleFree=+     mass_                &
         &        /     4.0d0                &
         &        /     Pi                   &
         &        /self%densityNormalization &
         &        /self%scaleLength**3
    if      (massScaleFree <= 0.0d0             ) then
       radius=+0.0d0
    else if (massScaleFree <  massScaleFreeSmall) then
       ! Use a series solution for very small radii.
       radius=+                     sqrt(2.0d0)*massScaleFree**0.5d0 &
            & +    4.0d0/     3.0d0            *massScaleFree        &
            & +   13.0d0/     9.0d0/sqrt(2.0d0)*massScaleFree**1.5d0 &
            & +   92.0d0/   135.0d0            *massScaleFree**2     &
            & +  313.0d0/   540.0d0/sqrt(2.0d0)*massScaleFree**2.5d0 &
            & + 1928.0d0/  8505.0d0            *massScaleFree**3     &
            & +56201.0d0/340200.0d0/sqrt(2.0d0)*massScaleFree**3.5d0 &
            & +  358.0d0/  1701.0d0            *massScaleFree**4
    else
       radius=-1.0d0                          &
            & /Lambert_W0(                    &
            &             -exp(               &
            &                  -1.0d0         &
            &                  -massScaleFree &
            &                 )               &
            &            )                    &
            & -1.0d0
    end if
    radius=+radius           &
         & *self%scaleLength
    return
  end function nfwRadiusEnclosingMass
  
  double precision function nfwRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{RST
    Computes the radius enclosing a given mean density for nfw mass distributions.
    !!}
    implicit none
    class           (massDistributionNFW), intent(inout), target       :: self
    double precision                     , intent(in   )               :: density
    double precision                     , intent(in   ), optional     :: radiusGuess
    double precision                     , allocatable  , dimension(:) :: radii
    double precision                                                   :: densityScaleFree
    integer                                                            :: i

    densityScaleFree=+density                   &
         &           /self%densityNormalization
    ! Extend the tabulation - by one octave at a time, and preserving everything already computed - until it spans the
    ! requested density. The mean enclosed density decreases with radius.
    if (densityScaleFree_%pointsPerOctave == 0) call densityScaleFree_%reset(countRadiiPerOctave,increasing=.false.)
    do while (.not.densityScaleFree_%brackets(densityScaleFree))
       call densityScaleFree_%expand(densityScaleFree)
       radii=densityScaleFree_%abscissae()
       do i=1,size(radii)
          if (densityScaleFree_%isComputed(i)) cycle
          call densityScaleFree_%set(i,densityEnclosedScaleFree(radii(i)))
       end do
       call densityScaleFree_%build()
    end do
    radius=+densityScaleFree_%invert     (densityScaleFree) &
         & *self             %scaleLength
    return
  end function nfwRadiusEnclosingDensity

  elemental double precision function massEnclosedScaleFree(radius) result(mass)
    !!{RST
    Evaluate the mass enclosed by a given radius in a scale-free NFW mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    double precision, parameter     :: minimumRadiusForExactSolution   =1.0d-6
    double precision, parameter     :: nfwNormalizationFactorUnitRadius=log(2.0d0)-0.5d0 ! Precomputed NFW normalization factor for unit radius.
    
    if      (radius == 1.0d0                        ) then
       mass=nfwNormalizationFactorUnitRadius
    else if (radius >= minimumRadiusForExactSolution) then
       mass=log(1.0d0+radius)-radius/(1.0d0+radius)
    else
       mass=radius**2*(0.5d0+radius*(-2.0d0/3.0d0+radius*(0.75d0+radius*(-0.8d0))))
    end if
    mass   =+4.0d0 &
         &  *Pi    &
         &  *mass
    return
  end function massEnclosedScaleFree

  elemental double precision function densityEnclosedScaleFree(radius) result(density)
    !!{RST
    Evaluate the mean enclosed density at a given radius in a scale-free NFW mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    
    density=+3.0d0                            &
         &  /4.0d0                            &
         &  /Pi                               &
         &  *massEnclosedScaleFree(radius)    &
         &  /                      radius **3
    return
  end function densityEnclosedScaleFree  

  double precision function nfwRadiusFromSpecificAngularMomentum(self,angularMomentumSpecific) result(radius)
    !!{RST
    Computes the radius corresponding to a given specific angular momentum for nfw mass distributions.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionNFW), intent(inout), target       :: self
    double precision                     , intent(in   )               :: angularMomentumSpecific
    double precision                     , allocatable  , dimension(:) :: radii
    double precision                                                   :: angularMomentumSpecificScaleFree
    integer                                                            :: i

    if (angularMomentumSpecific > 0.0d0) then
       angularMomentumSpecificScaleFree=+angularMomentumSpecific                 &
            &                           /sqrt(                                   &
            &                                 +gravitationalConstant_internal    &
            &                                 *self%densityNormalization         &
            &                                )                                   &
            &                           /      self%scaleLength              **2
       ! The enclosed specific angular momentum increases with radius.
       if (angularMomentumSpecificScaleFree_%pointsPerOctave == 0) call angularMomentumSpecificScaleFree_%reset(countRadiiPerOctave,increasing=.true.)
       do while (.not.angularMomentumSpecificScaleFree_%brackets(angularMomentumSpecificScaleFree))
          call angularMomentumSpecificScaleFree_%expand(angularMomentumSpecificScaleFree)
          radii=angularMomentumSpecificScaleFree_%abscissae()
          do i=1,size(radii)
             if (angularMomentumSpecificScaleFree_%isComputed(i)) cycle
             call angularMomentumSpecificScaleFree_%set(i,angularMomentumSpecificEnclosedScaleFree(radii(i)))
          end do
          call angularMomentumSpecificScaleFree_%build()
       end do
       radius=+angularMomentumSpecificScaleFree_%invert     (angularMomentumSpecificScaleFree) &
            & *self                             %scaleLength
    else
       radius=+0.0d0
    end if
    return
  end function nfwRadiusFromSpecificAngularMomentum

  elemental double precision function angularMomentumSpecificEnclosedScaleFree(radius) result(angularMomentumSpecific)
    !!{RST
    Evaluate the specific angular momentum at a given radius in a scale-free NFW mass distribution.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    angularMomentumSpecific=+sqrt(                               &
         &                        +massEnclosedScaleFree(radius) &
         &                        *                      radius  &
         &                       )
    return
  end function angularMomentumSpecificEnclosedScaleFree

  double precision function nfwVelocityRotationCurveMaximum(self) result(velocity)
    !!{RST
    Return the peak velocity in the rotation curve for an nfw mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionNFW), intent(inout) :: self
    double precision                     , parameter     :: circularVelocityMaximumScaleFree=0.4649909628174221d0 ! The circular velocity (in scale-free units) at the peak of the NFW rotation curve.
    !                                                                                                               Numerical value found using Mathematica.

    velocity=+circularVelocityMaximumScaleFree            &
         &   *sqrt(                                       &
         &         +4.0d0                                 &
         &         *Pi                                    &
         &         *self%densityNormalization             &
         &        )                                       &
         &   *      self%scaleLength
    if (.not.self%isDimensionless())                      &
         & velocity=+velocity                             &
         &          *sqrt(gravitationalConstant_internal)
    return
  end function nfwVelocityRotationCurveMaximum

  double precision function nfwRadiusRotationCurveMaximum(self) result(radius)
    !!{RST
    Return the peak velocity in the rotation curve for an nfw mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW         ), intent(inout), target :: self
    ! The radius (in scale-free units) at the peak of the NFW rotation curve. Numerical value found using Mathematica.
    double precision                              , parameter             :: radiusCircularVelocityMaximumScaleFree=2.162581587064612d0
    
    radius=+radiusCircularVelocityMaximumScaleFree &
         & *self%scaleLength
    return
  end function nfwRadiusRotationCurveMaximum

  logical function nfwPotentialIsAnalytic(self) result(isAnalytic)
    !!{RST
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionNFW), intent(inout) :: self

    isAnalytic=.true.
    return
  end function nfwPotentialIsAnalytic

  double precision function nfwPotential(self,coordinates,status) result(potential)
    !!{RST
    Return the potential at the specified ``coordinates`` in an nfw mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess     , structureErrorCodeInfinite
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Error                           , only : Error_Report
    implicit none
    class           (massDistributionNFW              ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                                             :: radiusScaleFree
    
    if (present(status)) status=structureErrorCodeSuccess
    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %scaleLength
    potential=+potentialScaleFree       (radiusScaleFree)    &
         &    *self%densityNormalization                     &
         &    *self%scaleLength                          **2
    if (.not.self%isDimensionless()) potential=+gravitationalConstant_internal &
         &                                     *potential
    return
  end function nfwPotential

  elemental double precision function potentialScaleFree(radius) result(potential)
    !!{RST
    Compute the potential in a scale-free NFW mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    double precision, parameter     :: radiusSmall=1.0d-6
    double precision                :: radiusTerm
    
    if (radius < radiusSmall) then
       ! Use a series solution for very small radii.
       radiusTerm=+1.0d0-radius/2.0d0+radius**2/3.0d0
    else
       ! Use the full expression for larger radii.
       radiusTerm=log(1.0d0+radius)/radius
    end if
    potential=-4.0d0      &
         &    *Pi         &
         &    *radiusTerm
    return
  end function potentialScaleFree

 double precision function potentialDifferenceScaleFree(radius1,radius2) result(potential)
    !!{RST
    Compute the potential difference in a scale-free NFW mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Comparison    , only : Values_Agree
    implicit none
    double precision, intent(in   ) :: radius1                            , radius2
    double precision, parameter     :: radiusSmall                 =1.0d-6
    double precision, parameter     :: toleranceRelative           =1.0d-3
    double precision                :: potentialGradientLogarithmic       , radiusDifferenceLogarithmic
    
    if (Values_Agree(radius1,radius2,relTol=toleranceRelative)) then
       if (radius1 < radiusSmall) then
          potentialGradientLogarithmic=-      radius1   / 2.0d0 &
               &                       +5.0d0*radius1**2/12.0d0
       else
          potentialGradientLogarithmic=-    1.0d0          &
               &                       +          radius1  &
               &                       /   (1.0d0+radius1) &
               &                       /log(1.0d0+radius1)
       end if
       radiusDifferenceLogarithmic=+1.0d0                                 &
            &                      -radius2                               &
            &                      /radius1
       potential                  =+potentialScaleFree          (radius1) &
            &                      *potentialGradientLogarithmic          &
            &                      *radiusDifferenceLogarithmic       
    else
       potential=+potentialScaleFree(radius1) &
            &    -potentialScaleFree(radius2)
    end if
    return
  end function potentialDifferenceScaleFree
  
  double precision function nfwFourierTransform(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{RST
    Compute the Fourier transform of the density profile at the given ``wavenumber`` in an NFW mass distribution, using the expression given in :cite:author:`cooray_halo_2002` (:cite:year:`cooray_halo_2002`; eqn. 81).
    !!}
    use :: Exponential_Integrals, only : Cosine_Integral, Sine_Integral
    implicit none
    class           (massDistributionNFW ), intent(inout) :: self
    double precision                      , intent(in   ) :: radiusOuter        , wavenumber
    double precision                                      :: wavenumberScaleFree, radiusOuterScaleFree

    waveNumberScaleFree =+waveNumber *self%scaleLength
    radiusOuterScaleFree=+radiusOuter/self%scaleLength
    fourierTransform    =+(                                                                                                                                                         &
         &                 +sin(+                     waveNumberScaleFree)*(Sine_Integral  ((1.0d0+radiusOuterScaleFree)*waveNumberScaleFree)-Sine_Integral  (waveNumberScaleFree)) &
         &                 -sin(+radiusOuterScaleFree*waveNumberScaleFree)/                 (1.0d0+radiusOuterScaleFree)/waveNumberScaleFree                                        &
         &                 +cos(+                     waveNumberScaleFree)*(Cosine_Integral((1.0d0+radiusOuterScaleFree)*waveNumberScaleFree)-Cosine_Integral(waveNumberScaleFree)) &
         &                )                                                                                                                                                         &
         &               /(log(1.0d0+radiusOuterScaleFree)-radiusOuterScaleFree/(1.0d0+radiusOuterScaleFree))
    return
  end function nfwFourierTransform
  
  double precision function nfwRadiusFreefall(self,time) result(radius)
    !!{RST
    Compute the freefall radius at the given ``time`` in an NFW mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr, gravitationalConstant_internal
    implicit none
    class           (massDistributionNFW), intent(inout) :: self
    double precision                     , intent(in   ) :: time
    double precision                                     :: timeScaleFree, timeScale

    if (time > 0.0d0) then
       timeScale    =+1.0d0/sqrt(                                &
            &                    +gravitationalConstant_internal &
            &                    *self%densityNormalization      &
            &                   )                                &
            &        *MpcPerKmPerSToGyr
       timeScaleFree=+time                                       &
            &        /timeScale
       call self%timeFreefallTabulate(timeScaleFree)
       radius=+timeFreefallScaleFree_%invert     (timeScaleFree) &
            & *self                  %scaleLength
    else
       ! For non-positive freefall times, return a zero freefall radius.
       radius=+0.0d0
    end if
    return   
  end function nfwRadiusFreefall
  
  double precision function nfwRadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{RST
    Compute the rate of increase of the freefall radius at the given ``time`` in an nfw mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr, gravitationalConstant_internal
    implicit none
    class           (massDistributionNFW), intent(inout) :: self
    double precision                     , intent(in   ) :: time
    double precision                                     :: timeScaleFree, timeScale

    if (time > 0.0d0) then
       timeScale    =+1.0d0/sqrt(                                &
            &                    +gravitationalConstant_internal &
            &                    *self%densityNormalization      &
            &                   )                                &
            &        *MpcPerKmPerSToGyr
       timeScaleFree=+time                                       &
            &        /timeScale
       call self%timeFreefallTabulate(timeScaleFree)
       radiusIncreaseRate=+timeFreefallScaleFree_%derivative(timeScaleFree) &
            &             *self                  %scaleLength               &
            &             /                       timeScale
    else
       ! For non-positive freefall times, return the limiting value for small radii.
       radiusIncreaseRate=+0.0d0
    end if
    return
  end function nfwRadiusFreefallIncreaseRate
  
  subroutine nfwTimeFreefallTabulate(self,timeScaleFree)
    !!{RST
    Tabulate the freefall radius at the given ``time`` in an NFW mass distribution.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (massDistributionNFW), intent(inout)               :: self
    double precision                     , intent(in   )               :: timeScaleFree
    double precision                     , allocatable  , dimension(:) :: radii
    double precision                                                   :: radiusStart
    integer                                                            :: i
    type            (integrator         )                              :: integrator_

    ! The freefall time increases with radius. Note that each point is an independent quadrature from the centre out to its
    ! own radius, so a point carried over by an extension is precisely the value which would be computed afresh.
    if (timeFreefallScaleFree_%pointsPerOctave == 0) call timeFreefallScaleFree_%reset(countRadiiPerOctave,increasing=.true.)
    if (.not.timeFreefallScaleFree_%brackets(timeScaleFree)) then
       integrator_=integrator(timeFreeFallIntegrand,toleranceRelative=1.0d-3)
       do while (.not.timeFreefallScaleFree_%brackets(timeScaleFree))
          call timeFreefallScaleFree_%expand(timeScaleFree)
          radii=timeFreefallScaleFree_%abscissae()
          do i=1,size(radii)
             if (timeFreefallScaleFree_%isComputed(i)) cycle
             call timeFreefallScaleFree_%set(i,timeFreefallScaleFree(radii(i)))
          end do
          call timeFreefallScaleFree_%build()
       end do
    end if
    return

  contains
    
    double precision function timeFreefallScaleFree(radius)
      !!{RST
      Evaluate the freefall time from a given radius in a scale-free NFW mass distribution.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: radius
      double precision, parameter     :: radiusSmall=4.0d-6

      if (radius > radiusSmall) then
         radiusStart          =                            radius
         timeFreefallScaleFree=integrator_%integrate(0.0d0,radius)
      else
         ! Use an approximation here, found by taking series expansions of the logarithms in the integrand and keeping only the
         ! first order terms.
         timeFreefallScaleFree=2.0d0*sqrt(radius/4.0d0/Pi)
      end if
      return
    end function timeFreefallScaleFree
    
    double precision function timeFreeFallIntegrand(radius)
      !!{RST
      Integrand used to find the freefall time in a scale-free NFW mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius
      double precision                :: potentialDifference
      
      if (radius == 0.0d0) then
         timeFreeFallIntegrand=+0.0d0
      else
         potentialDifference=+potentialDifferenceScaleFree(radiusStart,radius)
         if (potentialDifference > 0.0d0) then
            timeFreeFallIntegrand=+1.0d0                     &
                 &                /sqrt(                     &
                 &                      +2.0d0               &
                 &                      *potentialDifference &
                 &                     )
         else
            timeFreeFallIntegrand=+0.0d0
         end if
      end if
      return
    end function timeFreeFallIntegrand
    
  end subroutine nfwTimeFreefallTabulate
  
  double precision function nfwEnergyPotential(self,radiusOuter) result(energy)
    !!{RST
    Compute the potential energy within a given ``radius`` in an NFW mass distribution. This is

    .. math::

       W &=& - \frac{\mathrm{G}}{2} \rho_0^2 r_\mathrm{s}^5 \int_0^{x_\mathrm{out}} \frac{m^2(x)}{x^2} \mathrm{d} x, \nonumber \\
       &-& - \frac{\mathrm{G}}{2} \rho_0^2 r_\mathrm{s}^5 \left[ \frac{x}{1+x} - \frac{\log^2(1+x)}{x} + \frac{\left\{\log(1+x)-x/(1+x)\right\}^2}{x} \right],

    where :math:`x=r/r_\mathrm{s}` and :math:`m(x)` is the scale-free mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionNFW), intent(inout) :: self
    double precision                     , intent(in   ) :: radiusOuter
    double precision                                     :: radiusOuterScaleFree
    
    radiusOuterScaleFree=+     radiusOuter                                                                                   &
         &               /self%scaleLength
    energy=-gravitationalConstant_internal                                                                                   &
         & *self%scaleLength                           **5                                                                   &
         & *self%densityNormalization                  **2                                                                   &
         & *8.0d0                                                                                                            &
         & *Pi**2                                                                                                            &
         & *(                                                                                                                &
         &                                       +radiusOuterScaleFree/(1.0d0+radiusOuterScaleFree)                          &
         &   - log(1.0d0+radiusOuterScaleFree)**2                                                      /radiusOuterScaleFree &
         &   +(log(1.0d0+radiusOuterScaleFree)   -radiusOuterScaleFree/(1.0d0+radiusOuterScaleFree))**2/radiusOuterScaleFree &
         &  )
    return
  end function nfwEnergyPotential

  double precision function nfwEnergyKinetic(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{RST
    Compute the kinetic energy within a given ``radius`` in an NFW mass distribution. This is

    .. math::

       T &=& 6 \pi \mathrm{G} \rho_0^2 r_\mathrm{s}^5 \int_0^{x_\mathrm{out}} \rho(x) \sigma^2(x) x^2 \mathrm{d} x, \nonumber \\
       &=& \pi \mathrm{G} \rho_0^2 r_\mathrm{s}^5 \left[ 6 x^3 \text{Li}_2(-x)+x^3 (-\log (x))+\log (x+1) \left(3 x^3 \log (x+1)+((x-6) x+3) x-2\right)+\left(x \left(\pi ^2 x-7\right)+5\right) x+\frac{3}{x+1} \right],

    where :math:`x=r/r_\mathrm{s}`, :math:`\rho(x)` is the scale-free density, and :math:`\sigma^2(x)` is the scale-free velocity dispersion.
    !!}
    use :: Dilogarithms                    , only : Dilogarithm
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionNFW  ), intent(inout) :: self
    double precision                       , intent(in   ) :: radiusOuter
    class           (massDistributionClass), intent(inout) :: massDistributionEmbedding
    logical                                                :: analytic
    double precision                                       :: radiusOuterScaleFree

    analytic=.false.
    select type (massDistributionEmbedding)
    class is (massDistributionNFW)
       select type (kinematicsDistribution_ => massDistributionEmbedding%kinematicsDistribution_)
       class is (kinematicsDistributionNFW)
          analytic   =.true.
          radiusOuterScaleFree=+     radiusOuter                                                                      &
               &               /self%scaleLength          
          energy              =+gravitationalConstant_internal                                                        &
               &               *self%scaleLength                **5                                                   &
               &               *self%densityNormalization       **2                                                   &
               &               *4.0d0                                                                                 &
               &               *Pi**2                                                                                 &
               &               *(                                                                                     &
               &                 +(                                                                                   &
               &                   +2.0d0                                                                             &
               &                   +            radiusOuterScaleFree                                                  &
               &                   *(                                                                                 &
               &                     -2.0d0                                                                           &
               &                     +          radiusOuterScaleFree                                                  &
               &                     *(                                                                               &
               &                       -7.0d0                                                                         &
               &                       +Pi**2                                                                         &
               &                       *(1.0d0+radiusOuterScaleFree)                                                  &
               &                      )                                                                               &
               &                    )                                                                                 &         
               &                  )                                                                                   &
               &                 *                                                              radiusOuterScaleFree  &
               &                 +        radiusOuterScaleFree**4     *log        (+1.0d0+1.0d0/radiusOuterScaleFree) &
               &                 -        radiusOuterScaleFree**3     *log        (+            radiusOuterScaleFree) &
               &                 +(                                                                                   &
               &                   -    2.0d0                                                                         &
               &                   +          radiusOuterScaleFree                                                    &
               &                   -    3.0d0*radiusOuterScaleFree**2                                                 &
               &                   -    5.0d0*radiusOuterScaleFree**3                                                 &
               &                   +    3.0d0*radiusOuterScaleFree**3                                                 &
               &                   *   (1.0d0+radiusOuterScaleFree   )                                                &
               &                   *log(1.0d0+radiusOuterScaleFree   )                                                &
               &                  )                                                                                   &
               &                 *                                     log        (1.0d0+radiusOuterScaleFree)        &
               &                 +6.0d0                                                                               &
               &                 *             radiusOuterScaleFree**3                                                &
               &                 *     (1.0d0+radiusOuterScaleFree    )                                               &
               &                 *                                     Dilogarithm(     -radiusOuterScaleFree)        &
               &                )                                                                                     &
               &               /(1.0d0+radiusOuterScaleFree)
       end select
    end select
    if (.not.analytic) energy=self%energyKineticNumerical(radiusOuter,massDistributionEmbedding)
    return
  end function nfwEnergyKinetic

  subroutine nfwDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionNFW), intent(inout)           :: self
    type     (inputParameters    ), intent(inout)           :: descriptor
    logical                       , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18             )                          :: parameterLabel
    type     (inputParameters    )                          :: parameters
    !$GLC attributes unused :: includeFileModificationTimes
    
    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','NFW')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%scaleLength
    call parameters%addParameter('scaleLength'         ,trim(adjustl(parameterLabel)))
    return
  end subroutine nfwDescriptor

  !![
  <stateStoreTask function="massDistributionNFWStateStore"/>
  !!]
  subroutine massDistributionNFWStateStore(stateFile,gslStateFile,stateOperationID)
    !!{RST
    Write the tabulation state to file.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: massDistributionNFW',verbosity=verbosityLevelInfo)
    ! Store the tabulations in full - the lattice on which each was built, and the values computed on it. Recording only the
    ! range, as was done previously, left a restored object which reported that it spanned a value while holding no
    ! interpolator with which to invert it.
    call densityScaleFree_                %stateStore(stateFile)
    call angularMomentumSpecificScaleFree_%stateStore(stateFile)
    call timeFreefallScaleFree_           %stateStore(stateFile)
    return
  end subroutine massDistributionNFWStateStore

  !![
  <stateRetrieveTask function="massDistributionNFWStateRestore"/>
  !!]
  subroutine massDistributionNFWStateRestore(stateFile,gslStateFile,stateOperationID)
    !!{RST
    Retrieve the tabulation state from the file.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: massDistributionNFW',verbosity=verbosityLevelInfo)
    call densityScaleFree_                %stateRestore(stateFile)
    call angularMomentumSpecificScaleFree_%stateRestore(stateFile)
    call timeFreefallScaleFree_           %stateRestore(stateFile)
    return
  end subroutine massDistributionNFWStateRestore
  
