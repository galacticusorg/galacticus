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
  Implementation of an Einasto (e.g. \citealt{cardone_spherical_2005}) mass distribution class.
  !!}

  use :: Numerical_Interpolation, only : interpolator
  
  !![
  <massDistribution name="massDistributionEinasto">
    <description>
      An Einasto (e.g. \citealt{cardone_spherical_2005}) mass distribution class. The density profile is given by:
      \begin{equation}
      \rho_\mathrm{dark matter}(r) = \rho_{-2} \exp \left( - {2 \over \alpha} \left[ \left( {r \over r_{-2}} \right)^\alpha - 1
      \right] \right).
      \end{equation}
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionEinasto
     !!{
     The Einasto (e.g. \citealt{cardone_spherical_2005}) mass distribution.
     !!}
     private
     double precision                            :: densityNormalization                         , scaleLength                                  , &
          &                                         shapeParameter                               , massTotal_
     double precision                            :: enclosedMassRadiusPrevious                   , enclosedMassPrevious
     double precision                            :: massScaleFreeRadiusMinimum                   , massScaleFreeRadiusMaximum
     double precision                            :: massScaleFreeMinimum                         , massScaleFreeMaximum
     type            (interpolator), allocatable :: massScaleFree_
     double precision                            :: densityScaleFreeRadiusMinimum                , densityScaleFreeRadiusMaximum
     double precision                            :: densityScaleFreeMinimum                      , densityScaleFreeMaximum
     type            (interpolator), allocatable :: densityScaleFree_
     double precision                            :: angularMomentumSpecificScaleFreeRadiusMinimum, angularMomentumSpecificScaleFreeRadiusMaximum
     double precision                            :: angularMomentumSpecificScaleFreeMinimum      , angularMomentumSpecificScaleFreeMaximum
     type            (interpolator), allocatable :: angularMomentumSpecificScaleFree_
     double precision                            :: timeFreefallScaleFreeRadiusMinimum           , timeFreefallScaleFreeRadiusMaximum
     double precision                            :: timeFreefallScaleFreeMinimum                 , timeFreefallScaleFreeMaximum
     type            (interpolator), allocatable :: timeFreefallScaleFree_
   contains
     !![
     <methods>
       <method method="timeFreefallTabulate" description="Tabulate the freefall time as a function of radius in a scale-free Einasto mass distribution."/>
       <method method="timeFreefallMinimum"  description="Compute the minimum freefall time in a scale-free Einasto mass distribution."                 />
     </methods>
     !!]
     procedure :: massTotal                         => einastoMassTotal
     procedure :: density                           => einastoDensity
     procedure :: densityGradientRadial             => einastoDensityGradientRadial
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
     !!{
     Constructors for the \refClass{massDistributionEinasto} mass distribution class.
     !!}
     module procedure massDistributionEinastoConstructorParameters
     module procedure massDistributionEinastoConstructorInternal
  end interface massDistributionEinasto

contains

  function massDistributionEinastoConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionEinasto} mass distribution class which builds the object from a parameter
    set.
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
    <inputParameter>
      <name>shapeParameter</name>
      <description>The shape parameter, $\alpha$, of the Einasto profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>shapeParameter/4.0d0/Pi*(2.0d0/shapeParameter)**(3.0d0/shapeParameter)*exp(-2.0d0/shapeParameter)/Gamma_Function(3.0d0/shapeParameter)</defaultValue>
      <description>The density normalization of the Einasto profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scaleLength</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scale radius of the Einasto profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the Einasto profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>concentration</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The concentration of the Einasto profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>virialRadius</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The virial radius of the Einasto profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the Einasto profile is considered to be dimensionless.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
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
    !!{
    Internal constructor for ``einasto'' mass distribution class.
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
       self%scaleLength         =scaleLength
    else if (                            &
         &   present(concentration).and. &
         &   present(virialRadius )      &
         &  ) then
       self%scaleLength=virialRadius/concentration
    else
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
       call Error_Report('either "densityNormalization", or "mass" and "virialRadius" must be specified'//{introspection:location})
    end if
    ! Determine if profile is dimensionless.
    if      (present(dimensionless     )) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if

    ! Initialize memoized results.
    self%enclosedMassPrevious                         =-huge(0.0d0)
    self%enclosedMassRadiusPrevious                   =-huge(0.0d0)
    self%densityScaleFreeMinimum                      =+huge(0.0d0)
    self%densityScaleFreeMaximum                      =-huge(0.0d0)
    self%densityScaleFreeRadiusMinimum                =+1.0d0
    self%densityScaleFreeRadiusMaximum                =+1.0d0
    self%angularMomentumSpecificScaleFreeMinimum      =+huge(0.0d0)
    self%angularMomentumSpecificScaleFreeMaximum      =-huge(0.0d0)
    self%angularMomentumSpecificScaleFreeRadiusMinimum=+1.0d0
    self%angularMomentumSpecificScaleFreeRadiusMaximum=+1.0d0
    self%timeFreefallScaleFreeMinimum                 =+huge(0.0d0)
    self%timeFreefallScaleFreeMaximum                 =-huge(0.0d0)
    self%timeFreefallScaleFreeRadiusMinimum           =+1.0d0
    self%timeFreefallScaleFreeRadiusMaximum           =+1.0d0
    self%massScaleFreeMinimum                         =+huge(0.0d0)
    self%massScaleFreeMaximum                         =-huge(0.0d0)
    self%massScaleFreeRadiusMinimum                   =+1.0d0
    self%massScaleFreeRadiusMaximum                   =+1.0d0
    return
  end function massDistributionEinastoConstructorInternal

  double precision function einastoMassTotal(self)
    !!{
    Return the total mass in an Einasto mass distribution.
    !!}
    implicit none
    class(massDistributionEinasto), intent(inout) :: self
 
    einastoMassTotal=self%massTotal_   
    return
  end function einastoMassTotal

  double precision function einastoDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an Einasto mass distribution.
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
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an Einasto \citep{navarro_structure_1996} mass distribution.
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

  double precision function einastoDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Computes radial moments of the density in an Einasto \citep{navarro_structure_1996} mass distribution.
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
      !!{
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
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for einasto mass distributions.
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
    !!{
    Computes the radius enclosing a given mass or mass fraction for einasto mass distributions.
    !!}    
    use :: Numerical_Ranges, only : Make_Range  , rangeTypeLogarithmic
    use :: Error           , only : Error_Report
    implicit none
    class           (massDistributionEinasto), intent(inout), target       :: self
    double precision                         , intent(in   ), optional     :: mass                       , massFractional
    double precision                                                       :: mass_
    double precision                         , allocatable  , dimension(:) :: radii                      , masses
    double precision                         , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                       :: massScaleFree              , mass_
    integer                                                                :: countRadii                 , i

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
    if     (                                            &
         &   massScaleFree <= self%massScaleFreeMinimum &
         &  .or.                                        &
         &   massScaleFree >  self%massScaleFreeMaximum &
         & ) then
       do while (massEnclosedScaleFree(self%massScaleFreeRadiusMinimum,self%shapeParameter) >= massScaleFree)
          self%massScaleFreeRadiusMinimum=0.5d0*self%massScaleFreeRadiusMinimum
       end do
       do while (massEnclosedScaleFree(self%massScaleFreeRadiusMaximum,self%shapeParameter) <  massScaleFree)
          self%massScaleFreeRadiusMaximum=2.0d0*self%massScaleFreeRadiusMaximum
       end do
       countRadii=int(log10(self%massScaleFreeRadiusMaximum/self%massScaleFreeRadiusMinimum)*countRadiiPerDecade)+1
       if (allocated(self%massScaleFree_)) deallocate(self%massScaleFree_)
       allocate(     radii         (countRadii))
       allocate(     masses        (countRadii))
       allocate(self%massScaleFree_            )
       radii                       =Make_Range(self%massScaleFreeRadiusMinimum,self%massScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
       do i=1,countRadii
          masses                (i)=massEnclosedScaleFree(           radii(i),self%shapeParameter)
       end do
       self%massScaleFreeMinimum   =masses               (         1                             )
       self%massScaleFreeMaximum   =masses               (countRadii                             )
       self%massScaleFree_         =interpolator         (masses    ,radii                       )
    end if
    radius=+self%massScaleFree_%interpolate(massScaleFree) &
         & *self%scaleLength
    return
  end function einastoRadiusEnclosingMass
  
  double precision function einastoRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mean density for Einasto mass distributions.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionEinasto), intent(inout), target       :: self
    double precision                         , intent(in   )               :: density
    double precision                         , intent(in   ), optional     :: radiusGuess
    double precision                         , allocatable  , dimension(:) :: radii                      , densities
    double precision                         , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                       :: densityScaleFree
    integer                                                                :: countRadii                 , i

    densityScaleFree=+density                   &
         &           /self%densityNormalization
    if     (                                                  &
         &   densityScaleFree <= self%densityScaleFreeMinimum &
         &  .or.                                              &
         &   densityScaleFree >  self%densityScaleFreeMaximum &
         & ) then
       do while (densityEnclosedScaleFree(self%densityScaleFreeRadiusMinimum,self%shapeParameter) <  densityScaleFree)
          self%densityScaleFreeRadiusMinimum=0.5d0*self%densityScaleFreeRadiusMinimum
       end do
       do while (densityEnclosedScaleFree(self%densityScaleFreeRadiusMaximum,self%shapeParameter) >= densityScaleFree)
          self%densityScaleFreeRadiusMaximum=2.0d0*self%densityScaleFreeRadiusMaximum
       end do
       countRadii=int(log10(self%densityScaleFreeRadiusMaximum/self%densityScaleFreeRadiusMinimum)*countRadiiPerDecade)+1
       if (allocated(self%densityScaleFree_)) deallocate(self%densityScaleFree_)
       allocate(     radii            (countRadii))
       allocate(     densities        (countRadii))
       allocate(self%densityScaleFree_            )
       radii                          = Make_Range(self%densityScaleFreeRadiusMinimum,self%densityScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
       do i=1,countRadii
          densities                (i)=-densityEnclosedScaleFree(radii     (i),self%shapeParameter)
       end do
       self%densityScaleFreeMinimum   =-densities               (countRadii                       )
       self%densityScaleFreeMaximum   =-densities               (         1                       )
       self%densityScaleFree_         = interpolator            (densities    ,radii              )
    end if
    radius=+self%densityScaleFree_%interpolate(-densityScaleFree) &
         & *self%scaleLength
    return    
  end function einastoRadiusEnclosingDensity

  double precision function massEnclosedScaleFree(radius,shapeParameter) result(mass)
    !!{
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
    !!{
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
    !!{
    Computes the radius corresponding to a given specific angular momentum for einasto mass distributions.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Ranges                , only : Make_Range                    , rangeTypeLogarithmic
    implicit none
    class           (massDistributionEinasto), intent(inout), target       :: self
    double precision                         , intent(in   )               :: angularMomentumSpecific
    double precision                         , allocatable  , dimension(:) :: radii                                   , angularMomentaSpecific
    double precision                         , parameter                   :: countRadiiPerDecade             =100.0d0
    double precision                                                       :: angularMomentumSpecificScaleFree
    integer                                                                :: countRadii                              , i

    if (angularMomentumSpecific > 0.0d0) then
       angularMomentumSpecificScaleFree=+angularMomentumSpecific                 &
            &                           /sqrt(                                   &
            &                                 +gravitationalConstant_internal    &
            &                                 *self%densityNormalization         &
            &                                )                                   &
            &                           /      self%scaleLength              **2
       if     (                                                                                  &
            &   angularMomentumSpecificScaleFree <= self%angularMomentumSpecificScaleFreeMinimum &
            &  .or.                                                                              &
            &   angularMomentumSpecificScaleFree >  self%angularMomentumSpecificScaleFreeMaximum &
            & ) then
          do while (angularMomentumSpecificEnclosedScaleFree(self%angularMomentumSpecificScaleFreeRadiusMinimum,self%shapeParameter) >= angularMomentumSpecificScaleFree)
             self%angularMomentumSpecificScaleFreeRadiusMinimum=0.5d0*self%angularMomentumSpecificScaleFreeRadiusMinimum
          end do
          do while (angularMomentumSpecificEnclosedScaleFree(self%angularMomentumSpecificScaleFreeRadiusMaximum,self%shapeParameter) <  angularMomentumSpecificScaleFree)
             self%angularMomentumSpecificScaleFreeRadiusMaximum=2.0d0*self%angularMomentumSpecificScaleFreeRadiusMaximum
          end do
          countRadii=int(log10(self%angularMomentumSpecificScaleFreeRadiusMaximum/self%angularMomentumSpecificScaleFreeRadiusMinimum)*countRadiiPerDecade)+1
          if (allocated(self%angularMomentumSpecificScaleFree_)) deallocate(self%angularMomentumSpecificScaleFree_)
          allocate(     radii                            (countRadii))
          allocate(     angularMomentaSpecific           (countRadii))
          allocate(self%angularMomentumSpecificScaleFree_            )
          radii                                          =Make_Range(self%angularMomentumSpecificScaleFreeRadiusMinimum,self%angularMomentumSpecificScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
          do i=1,countRadii
             angularMomentaSpecific                   (i)=angularMomentumSpecificEnclosedScaleFree(                       radii(i),self%shapeParameter)
          end do
          self%angularMomentumSpecificScaleFreeMinimum   =angularMomentaSpecific                  (                     1                             )
          self%angularMomentumSpecificScaleFreeMaximum   =angularMomentaSpecific                  (            countRadii                             )
          self%angularMomentumSpecificScaleFree_         =interpolator                            (angularMomentaSpecific,radii                       )
       end if
       radius=+self%angularMomentumSpecificScaleFree_%interpolate(angularMomentumSpecificScaleFree) &
            & *self%scaleLength
    else
       radius=+0.0d0
    end if
    return    
  end function einastoRadiusFromSpecificAngularMomentum

  double precision function angularMomentumSpecificEnclosedScaleFree(radius,shapeParameter) result(angularMomentumSpecific)
    !!{
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
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionEinasto), intent(inout) :: self

    isAnalytic=.true.
    return
  end function einastoPotentialIsAnalytic

  double precision function einastoPotential(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in an einasto mass distribution.
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
    !!{
    Compute the potential in a scale-free Einasto mass distribution. Uses the results from \cite{retana-montenegro_analytical_2012},
    their equations (19) and (20), but with different normalizations for the density and scale radius.
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
    !!{
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
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in an Einasto mass distribution.
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
    radius=+self%timeFreefallScaleFree_%interpolate(timeScaleFree) &
         & *self%scaleLength
    return   
  end function einastoRadiusFreefall
  
  double precision function einastoRadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in an einasto mass
    distribution.
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
    !!{
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
    !!{
    Tabulate the freefall radius at the given {\normalfont \ttfamily time} in an Einasto mass distribution.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges     , only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionEinasto), intent(inout)               :: self
    double precision                         , intent(in   )               :: timeScaleFree
    double precision                         , allocatable  , dimension(:) :: radii                      , timesFreefall
    double precision                         , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                       :: radiusStart
    integer                                                                :: countRadii                 , i
    type            (integrator             )                              :: integrator_

    if     (                                                    &
         &   timeScaleFree <= self%timeFreefallScaleFreeMinimum &
         &  .or.                                                &
         &   timeScaleFree >  self%timeFreefallScaleFreeMaximum &
         & ) then
       integrator_=integrator(timeFreeFallIntegrand,toleranceRelative=1.0d-6)
       do while (timeFreefallScaleFree(self%timeFreefallScaleFreeRadiusMinimum) >= timeScaleFree)
          self%timeFreefallScaleFreeRadiusMinimum=0.5d0*self%timeFreefallScaleFreeRadiusMinimum
       end do
       do while (timeFreefallScaleFree(self%timeFreefallScaleFreeRadiusMaximum) <  timeScaleFree)
          self%timeFreefallScaleFreeRadiusMaximum=2.0d0*self%timeFreefallScaleFreeRadiusMaximum
       end do
       countRadii=int(log10(self%timeFreefallScaleFreeRadiusMaximum/self%timeFreefallScaleFreeRadiusMinimum)*countRadiiPerDecade)+1
       if (allocated(self%timeFreefallScaleFree_)) deallocate(self%timeFreefallScaleFree_)
       allocate(     radii                 (countRadii))
       allocate(     timesFreefall         (countRadii))
       allocate(self%timeFreefallScaleFree_            )
       radii=Make_Range(self%timeFreefallScaleFreeRadiusMinimum,self%timeFreefallScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
       do i=1,countRadii
          timesFreefall(i)=timeFreefallScaleFree(radii(i))
       end do
       self%timeFreefallScaleFreeMinimum=timesFreefall(            1      )
       self%timeFreefallScaleFreeMaximum=timesFreefall(   countRadii      )
       self%timeFreefallScaleFree_      =interpolator (timesFreefall,radii)
    end if
    return
    
  contains
    
    double precision function timeFreefallScaleFree(radius)
      !!{
      Evaluate the freefall time from a given radius in a scale-free Einasto mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      radiusStart          =                            radius
      timeFreefallScaleFree=integrator_%integrate(0.0d0,radius)
      return
    end function timeFreefallScaleFree
    
    double precision function timeFreeFallIntegrand(radius)
      !!{
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
    !!{
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

