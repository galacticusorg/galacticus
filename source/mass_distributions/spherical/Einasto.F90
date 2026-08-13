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
  Implementation of an Einasto (e.g. :cite:author:`cardone_spherical_2005` :cite:year:`cardone_spherical_2005`) mass distribution class.
  !!}

  use :: Tabulations_Inverse, only : tabulationInverse
  
  !![
  <massDistribution name="massDistributionEinasto" docformat="rst">
    <description>
    An Einasto (e.g. :cite:author:`cardone_spherical_2005` :cite:year:`cardone_spherical_2005`) mass distribution class. The density profile is given by:

    .. math::

       \rho_\mathrm{dark matter}(r) = \rho_{-2} \exp \left( - {2 \over \alpha} \left[ \left( {r \over r_{-2}} \right)^\alpha - 1
       \right] \right).
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionEinasto
     !!{RST
     The Einasto (e.g. :cite:author:`cardone_spherical_2005` :cite:year:`cardone_spherical_2005`) mass distribution.
     !!}
     private
     double precision                    :: densityNormalization             , scaleLength           , &
          &                                 shapeParameter                   , massTotal_
     double precision                    :: enclosedMassRadiusPrevious       , enclosedMassPrevious
     ! Tabulations of the scale-free profile, built for inversion. These are held per object because the scale-free Einasto
     ! profile depends on the shape parameter, so each is specific to the object which owns it.
     type            (tabulationInverse) :: massScaleFree_                   , densityScaleFree_     , &
          &                                 angularMomentumSpecificScaleFree_, timeFreefallScaleFree_
   contains
     !![
     <methods docformat="rst">
       <method method="timeFreefallTabulate" description="Tabulate the freefall time as a function of radius in a scale-free Einasto mass distribution."/>
       <method method="timeFreefallMinimum"  description="Compute the minimum freefall time in a scale-free Einasto mass distribution."                 />
     </methods>
     !!]
     procedure :: massTotal                         => einastoMassTotal
     procedure :: density                           => einastoDensity
     procedure :: densityGradientRadial             => einastoDensityGradientRadial
     procedure :: densitySlopeLogarithmicCentral    => einastoDensitySlopeLogarithmicCentral
     procedure :: densityRadialMoment               => einastoDensityRadialMoment
     procedure :: massEnclosedBySphere              => einastoMassEnclosedBySphere
     procedure :: radiusEnclosingMass               => einastoRadiusEnclosingMass
     procedure :: radiusEnclosingDensity            => einastoRadiusEnclosingDensity
     procedure :: radiusFromSpecificAngularMomentum => einastoRadiusFromSpecificAngularMomentum
     procedure :: radiusFreefall                    => einastoRadiusFreefall
     procedure :: radiusFreefallIncreaseRate        => einastoRadiusFreefallIncreaseRate
     procedure :: timeFreefallTabulate              => einastoTimeFreefallTabulate
     procedure :: timeFreefallMinimum               => einastoTimeFreefallMinimum
     procedure :: potentialIsAnalytic               => einastoPotentialIsAnalytic
     procedure :: potential                         => einastoPotential
     procedure :: descriptor                        => einastoDescriptor
  end type massDistributionEinasto
  
  interface massDistributionEinasto
     !!{RST
     Constructors for the :galacticus-class:`massDistributionEinasto` mass distribution class.
     !!}
     module procedure massDistributionEinastoConstructorParameters
     module procedure massDistributionEinastoConstructorInternal
  end interface massDistributionEinasto

  ! Density of the tabulations. The bounds these reach are found by halving and doubling from a unit seed, so they are exact
  ! powers of two and are already points of an octave lattice. Thirty points per octave is 99.66 per decade, against the
  ! hundred per decade used before.
  integer, parameter :: countRadiiPerOctave=30

contains

  function massDistributionEinastoConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`massDistributionEinasto` mass distribution class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    use :: Numerical_Constants_Math  , only : Pi
    use :: Gamma_Functions           , only : Gamma_Function
    implicit none
    type            (massDistributionEinasto)                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    double precision                                         :: mass                , scaleLength  , &
         &                                                      densityNormalization, concentration, &
         &                                                      virialRadius        , shapeParameter
    logical                                                  :: dimensionless
    type            (varying_string         )                :: componentType
    type            (varying_string         )                :: massType

    !![
    <inputParameter docformat="rst">
      <name>shapeParameter</name>
      <description>
      The shape parameter, :math:`\alpha`, of the Einasto profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>densityNormalization</name>
      <defaultValue>shapeParameter/4.0d0/Pi*(2.0d0/shapeParameter)**(3.0d0/shapeParameter)*exp(-2.0d0/shapeParameter)/Gamma_Function(3.0d0/shapeParameter)</defaultValue>
      <description>
      The density normalization of the Einasto profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>scaleLength</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The scale radius of the Einasto profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The total mass (in :math:`\mathrm{M}_\odot`) of the Einasto profile, used to set the density normalization :math:`\rho_{-2}` when ``densityNormalization`` is not supplied.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>concentration</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The concentration of the Einasto profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>virialRadius</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The virial radius of the Einasto profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true the Einasto profile is considered to be dimensionless.
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
     <call>self=massDistributionEinasto(shapeParameter=shapeParameter,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
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
  end function massDistributionEinastoConstructorParameters

  function massDistributionEinastoConstructorInternal(shapeParameter,scaleLength,concentration,densityNormalization,mass,virialRadius,dimensionless,componentType,massType) result(self)
    !!{RST
    Internal constructor for "einasto" mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    use :: Gamma_Functions         , only : Gamma_Function
    implicit none
    type            (massDistributionEinasto     )                          :: self
    double precision                              , intent(in   )           :: shapeParameter
    double precision                              , intent(in   ), optional :: scaleLength         , concentration, &
         &                                                                     densityNormalization, mass         , &
         &                                                                     virialRadius
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree
    !![
    <constructorAssign variables="shapeParameter, componentType, massType"/>
    !!]

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
       self%densityNormalization=+densityNormalization
       self%massTotal_          =+4.0d0*Pi*densityNormalization*self%scaleLength**3/shapeParameter/(2.0d0/shapeParameter)**(3.0d0/shapeParameter)*Gamma_Function(3.0d0/shapeParameter)*exp(2.0d0/shapeParameter)
    else if (                                   &
         &   present(mass                ).and. &
         &   present(virialRadius        )      &
         &  ) then
       radiusScaleFree          =+virialRadius/self%scaleLength
       self%densityNormalization=+mass/4.0d0/Pi/self%scaleLength**3*shapeParameter*exp(-2.0d0/shapeParameter)*(2.0d0/shapeParameter)**(3.0d0/shapeParameter)/Gamma_Function(3.0d0/shapeParameter)
       self%massTotal_          =+mass
    else
       self%densityNormalization=+0.0d0
       self%massTotal_          =+0.0d0
       call Error_Report('either "densityNormalization", or "mass" and "virialRadius" must be specified'//{introspection:location})
    end if
    ! Determine if profile is dimensionless.
    if      (present(dimensionless     )) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! Initialize memoized results.
    self%enclosedMassPrevious      =-huge(0.0d0)
    self%enclosedMassRadiusPrevious=-huge(0.0d0)
    ! Initialize the tabulations. This is done here, where the shape parameter is set, so that an object can never serve a
    ! tabulation built for a different shape.
    call self%massScaleFree_                   %reset(countRadiiPerOctave,increasing=.true. )
    call self%densityScaleFree_                %reset(countRadiiPerOctave,increasing=.false.)
    call self%angularMomentumSpecificScaleFree_%reset(countRadiiPerOctave,increasing=.true. )
    call self%timeFreefallScaleFree_           %reset(countRadiiPerOctave,increasing=.true. )
    return
  end function massDistributionEinastoConstructorInternal

  double precision function einastoMassTotal(self)
    !!{RST
    Return the total mass in an Einasto mass distribution.
    !!}
    implicit none
    class(massDistributionEinasto), intent(inout) :: self
 
    einastoMassTotal=self%massTotal_   
    return
  end function einastoMassTotal

  double precision function einastoDensity(self,coordinates) result(density)
    !!{RST
    Return the density at the specified ``coordinates`` in an Einasto mass distribution.
    !!}
    implicit none
    class           (massDistributionEinasto     ), intent(inout) :: self
    class           (coordinate                  ), intent(in   ) :: coordinates
    double precision                                              :: radiusScaleFree

    ! Compute the density at this position.
    radiusScaleFree=+coordinates%rSpherical          () &
         &          /self       %scaleLength
    density        =+self       %densityNormalization            &
         &          *exp(                                        &
         &               -(2.0d0/self%shapeParameter)            &
         &               *(                                      &
         &                 +radiusScaleFree**self%shapeParameter &
         &                 -1.0d0                                &
         &                )                                      &
         &              )
    return
  end function einastoDensity

  double precision function einastoDensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{RST
    Return the density at the specified ``coordinates`` in an Einasto :cite:p:`navarro_structure_1996` mass distribution.
    !!}
    implicit none
    class           (massDistributionEinasto), intent(inout), target   :: self
    class           (coordinate             ), intent(in   )           :: coordinates
    logical                                  , intent(in   ), optional :: logarithmic
    double precision                                                   :: radiusScaleFree
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    radiusScaleFree      =+coordinates%rSpherical()             &
         &                /self       %scaleLength
    if (radiusScaleFree <= 0.0d0) then
       densityGradientRadial=+0.0d0
    else
       densityGradientRadial=-2.0d0                                &
            &                *radiusScaleFree**self%shapeParameter
       if (.not.logarithmic_) densityGradientRadial=+            densityGradientRadial              &
            &                                       *self       %density              (coordinates) &
            &                                       /coordinates%rSpherical           (           )
    end if
    return
  end function einastoDensityGradientRadial

  double precision function einastoDensitySlopeLogarithmicCentral(self) result(slope)
    !!{RST
    Return the central logarithmic slope of the density profile in an Einasto mass distribution. The profile has a finite
    central density, so the density is evaluable at zero radius.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class(massDistributionEinasto), intent(inout) :: self
    type (coordinateSpherical    )                :: coordinates

    coordinates=[0.0d0,0.0d0,0.0d0]
    slope      =self%densityGradientRadial(coordinates,logarithmic=.true.)
    return
  end function einastoDensitySlopeLogarithmicCentral

  double precision function einastoDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{RST
    Computes radial moments of the density in an Einasto :cite:p:`navarro_structure_1996` mass distribution.
    !!}
    implicit none
    class           (massDistributionEinasto), intent(inout)           :: self
    double precision                         , intent(in   )           :: moment
    double precision                         , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    logical                                  , intent(  out), optional :: isInfinite
    double precision                                                   :: radialMomentMinimum, radialMomentMaximum

    densityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    if ((.not.present(radiusMinimum) .or. radiusMinimum <= 0.0d0) .and. moment <= -1) then
       if (present(isInfinite)) then
          isInfinite=.true.
          return
       else
          call Error_Report('radial moment is infinite'//{introspection:location})
       end if
    end if
    if (present(radiusMinimum)) then
       radialMomentMinimum=radialMomentScaleFree(radiusMinimum/self%scaleLength)
    else
       radialMomentMinimum=radialMomentScaleFree(                         0.0d0)
    end if
    if (present(radiusMaximum)) then
       radialMomentMaximum=radialMomentScaleFree(radiusMaximum/self%scaleLength)
    else
       radialMomentMaximum=0.0d0
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
      Provides the scale-free part of the radial moment of the Einasto density profile.
      !!}
      use :: Gamma_Functions, only : Gamma_Function, Gamma_Function_Incomplete
      implicit none
      double precision, intent(in   ) :: radius

      radialMomentScaleFree=-exp                      ( 2.0d0        /self%shapeParameter                                                      )                                       &
           &                *Gamma_Function           ((1.0d0+moment)/self%shapeParameter                                                      )                                       &
           &                *Gamma_Function_Incomplete((1.0d0+moment)/self%shapeParameter,2.0d0*radius**self%shapeParameter/self%shapeParameter)                                       &
           &                /                                         self%shapeParameter                                                                                              &
           &                /                         ( 2.0d0        /self%shapeParameter                                                      )**((1.0d0+moment)/self%shapeParameter)
      return
    end function radialMomentScaleFree

  end function einastoDensityRadialMoment

  double precision function einastoMassEnclosedBySphere(self,radius) result(mass)
    !!{RST
    Computes the mass enclosed within a sphere of given ``radius`` for einasto mass distributions.
    !!}
    implicit none
    class           (massDistributionEinasto), intent(inout), target :: self
    double precision                         , intent(in   )         :: radius
    double precision                                                 :: radiusScaleFree
    
    if (radius /= self%enclosedMassRadiusPrevious) then
       self%enclosedMassRadiusPrevious=+     radius
       radiusScaleFree                =+     radius                                                        &
            &                          /self%scaleLength
       self%enclosedMassPrevious      =+     massEnclosedScaleFree(radiusScaleFree,self%shapeParameter)    &
            &                          *self%densityNormalization                                          &
            &                          *self%scaleLength                                               **3
    end if
    mass=self%enclosedMassPrevious
    return
  end function einastoMassEnclosedBySphere
  
  double precision function einastoRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{RST
    Computes the radius enclosing a given mass or mass fraction for einasto mass distributions.
    !!}    
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionEinasto), intent(inout), target       :: self
    double precision                         , intent(in   ), optional     :: mass         , massFractional
    double precision                         , allocatable  , dimension(:) :: radii
    double precision                                                       :: massScaleFree, mass_
    integer                                                                :: i

    mass_=0.0d0
    if (present(mass)) then
       mass_=                    mass
    else if (present(massFractional)) then
       mass_=massFractional*self%massTotal_
    else
       call Error_Report('either mass or massFractional must be supplied'//{introspection:location})
    end if
    massScaleFree=+     mass_                   &
         &        /self%densityNormalization    &
         &        /self%scaleLength         **3
    do while (.not.self%massScaleFree_%brackets(massScaleFree))
       call self%massScaleFree_%expand(massScaleFree)
       radii=self%massScaleFree_%abscissae()
       do i=1,size(radii)
          if (self%massScaleFree_%isComputed(i)) cycle
          call self%massScaleFree_%set(i,massEnclosedScaleFree(radii(i),self%shapeParameter))
       end do
       call self%massScaleFree_%build()
    end do
    radius=+self%massScaleFree_%invert(massScaleFree) &
         & *self%scaleLength
    return
  end function einastoRadiusEnclosingMass
  
  double precision function einastoRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{RST
    Computes the radius enclosing a given mean density for Einasto mass distributions.
    !!}
    implicit none
    class           (massDistributionEinasto), intent(inout), target       :: self
    double precision                         , intent(in   )               :: density
    double precision                         , intent(in   ), optional     :: radiusGuess
    double precision                         , allocatable  , dimension(:) :: radii
    double precision                                                       :: densityScaleFree
    integer                                                                :: i

    densityScaleFree=+density                   &
         &           /self%densityNormalization
    do while (.not.self%densityScaleFree_%brackets(densityScaleFree))
       call self%densityScaleFree_%expand(densityScaleFree)
       radii=self%densityScaleFree_%abscissae()
       do i=1,size(radii)
          if (self%densityScaleFree_%isComputed(i)) cycle
          call self%densityScaleFree_%set(i,densityEnclosedScaleFree(radii(i),self%shapeParameter))
       end do
       call self%densityScaleFree_%build()
    end do
    radius=+self%densityScaleFree_%invert( densityScaleFree) &
         & *self%scaleLength
    return    
  end function einastoRadiusEnclosingDensity

  double precision function massEnclosedScaleFree(radius,shapeParameter) result(mass)
    !!{RST
    Evaluate the mass enclosed by a given radius in a scale-free Einasto mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Gamma_Functions         , only : Gamma_Function, Gamma_Function_Incomplete_Complementary
    implicit none
    double precision, intent(in   ) :: radius, shapeParameter

    mass   =+4.0d0                                                                                                                             &
         &  *Pi                                                                                                                                &
         &  /                                              shapeParameter                                                                      &
         &  /                                       (2.0d0/shapeParameter                                            )**(3.0d0/shapeParameter) &
         &  *exp                                    (2.0d0/shapeParameter                                            )                         &
         &  *Gamma_Function                         (3.0d0/shapeParameter                                            )                         &
         &  *Gamma_Function_Incomplete_Complementary(3.0d0/shapeParameter,2.0d0*radius**shapeParameter/shapeParameter)
    return
  end function massEnclosedScaleFree

  double precision function densityEnclosedScaleFree(radius,shapeParameter) result(density)
    !!{RST
    Evaluate the mean enclosed density at a given radius in a scale-free Einasto mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius, shapeParameter
    
    density=+3.0d0                                           &
         &  /4.0d0                                           &
         &  /Pi                                              &
         &  *massEnclosedScaleFree(radius,shapeParameter)    &
         &  /                      radius                **3
    return
  end function densityEnclosedScaleFree
  
  double precision function einastoRadiusFromSpecificAngularMomentum(self,angularMomentumSpecific) result(radius)
    !!{RST
    Computes the radius corresponding to a given specific angular momentum for einasto mass distributions.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionEinasto), intent(inout), target       :: self
    double precision                         , intent(in   )               :: angularMomentumSpecific
    double precision                         , allocatable  , dimension(:) :: radii
    double precision                                                       :: angularMomentumSpecificScaleFree
    integer                                                                :: i

    if (angularMomentumSpecific > 0.0d0) then
       angularMomentumSpecificScaleFree=+angularMomentumSpecific                 &
            &                           /sqrt(                                   &
            &                                 +gravitationalConstant_internal    &
            &                                 *self%densityNormalization         &
            &                                )                                   &
            &                           /      self%scaleLength              **2
       do while (.not.self%angularMomentumSpecificScaleFree_%brackets(angularMomentumSpecificScaleFree))
          call self%angularMomentumSpecificScaleFree_%expand(angularMomentumSpecificScaleFree)
          radii=self%angularMomentumSpecificScaleFree_%abscissae()
          do i=1,size(radii)
             if (self%angularMomentumSpecificScaleFree_%isComputed(i)) cycle
             call self%angularMomentumSpecificScaleFree_%set(i,angularMomentumSpecificEnclosedScaleFree(radii(i),self%shapeParameter))
          end do
          call self%angularMomentumSpecificScaleFree_%build()
       end do
       radius=+self%angularMomentumSpecificScaleFree_%invert(angularMomentumSpecificScaleFree) &
            & *self%scaleLength
    else
       radius=+0.0d0
    end if
    return    
  end function einastoRadiusFromSpecificAngularMomentum

  double precision function angularMomentumSpecificEnclosedScaleFree(radius,shapeParameter) result(angularMomentumSpecific)
    !!{RST
    Evaluate the specific angular momentum at a given radius in a scale-free Einasto mass distribution.
    !!}
    implicit none
    double precision, intent(in   ) :: radius, shapeParameter
    
    angularMomentumSpecific=+sqrt(                                              &
         &                        +massEnclosedScaleFree(radius,shapeParameter) &
         &                        *                      radius                 &
         &                       )
    return
  end function angularMomentumSpecificEnclosedScaleFree

  logical function einastoPotentialIsAnalytic(self) result(isAnalytic)
    !!{RST
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionEinasto), intent(inout) :: self

    isAnalytic=.true.
    return
  end function einastoPotentialIsAnalytic

  double precision function einastoPotential(self,coordinates,status) result(potential)
    !!{RST
    Return the potential at the specified ``coordinates`` in an einasto mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess     , structureErrorCodeInfinite
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Error                           , only : Error_Report
    implicit none
    class           (massDistributionEinasto          ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                                             :: radiusScaleFree
    
    if (present(status)) status=structureErrorCodeSuccess
    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %scaleLength
    potential=+potentialScaleFree       (radiusScaleFree,self%shapeParameter)    &
         &    *self%densityNormalization                                         &
         &    *self%scaleLength                                              **2
    if (.not.self%isDimensionless()) potential=+gravitationalConstant_internal &
         &                                     *potential
    return
  end function einastoPotential

  double precision function potentialScaleFree(radius,shapeParameter) result(potential)
    !!{RST
    Compute the potential in a scale-free Einasto mass distribution. Uses the results from :cite:t:`retana-montenegro_analytical_2012`, their equations (19) and (20), but with different normalizations for the density and scale radius.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Gamma_Functions         , only : Gamma_Function, Gamma_Function_Incomplete, Gamma_Function_Incomplete_Complementary
    implicit none
    double precision, intent(in   ) :: radius, shapeParameter

    if (radius <= 0.0d0) then
       potential=-2.0d0                                  **(+1.0d0-2.0d0/shapeParameter                                            ) &
            &    *shapeParameter                         **(       2.0d0/shapeParameter                                            ) &
            &    *exp                                      (       2.0d0/shapeParameter                                            ) &
            &    *Pi                                                                                                                 &
            &    *Gamma_Function                           (+1.0d0+2.0d0/shapeParameter                                            )
    else
       potential=-4.0d0                                                                                                              &
            &    *Pi                                                                                                                 &
            &    *exp(2.0d0/shapeParameter)                                                                                          &
            &    /          shapeParameter                                                                                           &
            &    /                                                                            radius                                 &
            &    *(                                                                                                                  &
            &      +                                                                          radius                                 &
            &      *(2.0d0/shapeParameter)               **(      -2.0d0/shapeParameter                                            ) &
            &      *Gamma_Function_Incomplete              (      +2.0d0/shapeParameter,2.0d0*radius**shapeParameter/shapeParameter) &
            &      *Gamma_Function                         (      +2.0d0/shapeParameter                                            ) &
            &      +shapeParameter                       **(      +3.0d0/shapeParameter                                            ) &
            &      / 8.0d0                               **(      +1.0d0/shapeParameter                                            ) &
            &      *Gamma_Function_Incomplete_Complementary(      +3.0d0/shapeParameter,2.0d0*radius**shapeParameter/shapeParameter) &
            &      *Gamma_Function                         (      +3.0d0/shapeParameter                                            ) &
            &     )
    end if
    return
  end function potentialScaleFree

 double precision function potentialDifferenceScaleFree(radius1,radius2,shapeParameter) result(potential)
    !!{RST
    Compute the potential difference in a scale-free Einasto mass distribution.
    !!}
    implicit none
    double precision, intent(in   ) :: radius1       , radius2, &
         &                             shapeParameter
    
    potential=+potentialScaleFree(radius1,shapeParameter) &
         &    -potentialScaleFree(radius2,shapeParameter)
    return
  end function potentialDifferenceScaleFree
  
  double precision function einastoRadiusFreefall(self,time) result(radius)
    !!{RST
    Compute the freefall radius at the given ``time`` in an Einasto mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr, gravitationalConstant_internal
    implicit none
    class           (massDistributionEinasto), intent(inout) :: self
    double precision                         , intent(in   ) :: time
    double precision                                         :: timeScaleFree, timeScale
    
    timeScale    =+1.0d0                                &
         &        /sqrt(                                &
         &              +gravitationalConstant_internal &
         &              *self%densityNormalization      &
         &             )                                &
         &        *MpcPerKmPerSToGyr
    timeScaleFree=+time                                 &
         &        /timeScale
    if (timeScaleFree <= self%timeFreefallMinimum()) then
       radius=0.0d0
       return
    end if
    call self%timeFreefallTabulate(timeScaleFree)
    radius=+self%timeFreefallScaleFree_%invert(timeScaleFree) &
         & *self%scaleLength
    return   
  end function einastoRadiusFreefall
  
  double precision function einastoRadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{RST
    Compute the rate of increase of the freefall radius at the given ``time`` in an einasto mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr, gravitationalConstant_internal
    implicit none
    class           (massDistributionEinasto), intent(inout) :: self
    double precision                         , intent(in   ) :: time
    double precision                                         :: timeScaleFree, timeScale

    timeScale    =+1.0d0                                &
         &        /sqrt(                                &
         &              +gravitationalConstant_internal &
         &              *self%densityNormalization      &
         &             )                                &
         &        *MpcPerKmPerSToGyr
    timeScaleFree=+time                                 &
         &        /timeScale
    if (timeScaleFree <= self%timeFreefallMinimum()) then
       radiusIncreaseRate=0.0d0
       return
    end if
    call self%timeFreefallTabulate(timeScaleFree)
    radiusIncreaseRate=+self%timeFreefallScaleFree_%derivative(timeScaleFree) &
         &             *self%scaleLength                                      &
         &             /     timeScale
    return
  end function einastoRadiusFreefallIncreaseRate

  double precision function einastoTimeFreefallMinimum(self) result(timeScaleFreeMinimum)
    !!{RST
    Compute the minimum freefall time in a scale-free Einasto profile.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(massDistributionEinasto), intent(inout) :: self

    timeScaleFreeMinimum=+sqrt(                    &
         &                     + 3.0d0             &
         &                     /16.0d0             &
         &                     *Pi                 &
         &                   )                     &
         &               *exp(                     &
         &                    -1.0d0               &
         &                    /self%shapeParameter &
         &                   )
    return
  end function einastoTimeFreefallMinimum
  
  subroutine einastoTimeFreefallTabulate(self,timeScaleFree)
    !!{RST
    Tabulate the freefall radius at the given ``time`` in an Einasto mass distribution.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (massDistributionEinasto), intent(inout)               :: self
    double precision                         , intent(in   )               :: timeScaleFree
    double precision                         , allocatable  , dimension(:) :: radii
    double precision                                                       :: radiusStart
    integer                                                                :: i
    type            (integrator             )                              :: integrator_

    ! Each point is an independent quadrature from the centre out to its own radius, so a point carried over by an extension
    ! is precisely the value which would be computed afresh.
    if (.not.self%timeFreefallScaleFree_%brackets(timeScaleFree)) then
       integrator_=integrator(timeFreeFallIntegrand,toleranceRelative=1.0d-6)
       do while (.not.self%timeFreefallScaleFree_%brackets(timeScaleFree))
          call self%timeFreefallScaleFree_%expand(timeScaleFree)
          radii=self%timeFreefallScaleFree_%abscissae()
          do i=1,size(radii)
             if (self%timeFreefallScaleFree_%isComputed(i)) cycle
             call self%timeFreefallScaleFree_%set(i,timeFreefallScaleFree(radii(i)))
          end do
          call self%timeFreefallScaleFree_%build()
       end do
    end if
    return
    
  contains
    
    double precision function timeFreefallScaleFree(radius)
      !!{RST
      Evaluate the freefall time from a given radius in a scale-free Einasto mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      radiusStart          =                            radius
      timeFreefallScaleFree=integrator_%integrate(0.0d0,radius)
      return
    end function timeFreefallScaleFree
    
    double precision function timeFreeFallIntegrand(radius)
      !!{RST
      Integrand used to find the freefall time in a scale-free Einasto mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius
      double precision                :: potentialDifference
      
      if (radius == 0.0d0) then
         timeFreeFallIntegrand=+0.0d0
      else
         potentialDifference=+potentialDifferenceScaleFree(radiusStart,radius,self%shapeParameter)
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
    
  end subroutine einastoTimeFreefallTabulate

  subroutine einastoDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionEinasto), intent(inout)           :: self
    type     (inputParameters        ), intent(inout)           :: descriptor
    logical                           , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                 )                          :: parameterLabel
    type     (inputParameters        )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','Einasto')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%scaleLength
    call parameters%addParameter('scaleLength'         ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%shapeParameter
    call parameters%addParameter('shapeParameter'      ,trim(adjustl(parameterLabel)))
    return
  end subroutine einastoDescriptor

