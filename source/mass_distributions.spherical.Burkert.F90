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
  Implementation of a Burkert \citep{burkert_structure_1995} mass distribution class.
  !!}

  use :: Numerical_Interpolation , only : interpolator
  use :: Numerical_Constants_Math, only : Pi
  
  !![
  <massDistribution name="massDistributionBurkert">
    <description>
    A mass distribution class which implements the \citep{burkert_structure_1995} density profile:
    \begin{equation}
      \rho_\mathrm{dark matter}(r) = \rho_0 \left(1+{r\over r_\mathrm{s}}\right)^{-1} \left(1+[{r\over
      r_\mathrm{s}}]^2\right)^{-1}.
    \end{equation}
    The mass enclosed within radius $r$ is given by
    \begin{equation}
    M(&lt;r) = \pi \rho_0 r_\mathrm{s}^3 \left[ 2 \log(1 + R) + \log(1 + R^2) -2 \tan^{-1}(R) \right]
    \end{equation}
    where $R=r/r_\mathrm{s}$. The associated gravitational potential is
    \begin{equation}
    \Phi(r) = - \frac{\mathrm{G} \pi \rho_0 r_\mathrm{s}^2}{R} \left[ (R-1) \log \left(R^2+1\right)-2 (R+1) \log (R+1)-2 (R+1) \cot^{-1}(R)+\pi \right]
    \end{equation}    
    The peak of the rotation curve occurs at $R=3.2446257246042642$ (found by numerical solution) at which point the rotation
    curve amplitude is 1.644297750532498, and the Fourier transform of the profile, $F(k) = \int_0^c 4 \pi r^2 \exp(-i k r)
    \rho(r) \mathrm{d} r / k r$ (needed in calculations of clustering using the halo model) is given by
    \begin{eqnarray}
    F(k) &amp;=&amp;  (1+i) \frac{\pi}{k m(c) } \left( \right.                                                              \nonumber \\
         &amp; &amp;   +      \exp( k) \left\{ -i \pi -\mathrm{E}_\mathrm{i}[-  k]+\mathrm{E}_\mathrm{i}[(-1+ic)k] \right\} \nonumber \\
         &amp; &amp;   +(1-i) \exp(-k) \left\{        +\mathrm{E}_\mathrm{i}[-i k]+\mathrm{E}_\mathrm{i}[(+i+ic)k] \right\} \nonumber \\
         &amp; &amp;   +   i  \exp(-k) \left\{        +\mathrm{E}_\mathrm{i}[+  k]+\mathrm{E}_\mathrm{i}[(+1+ic)k] \right\} \nonumber \\
         &amp; &amp;  \left. \right).
    \end{eqnarray}
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionBurkert
     !!{
     The \citep{burkert_structure_1995} mass distribution.
     !!}
     private
     double precision                            :: densityNormalization                         , scaleLength
     double precision                            :: densityScaleFreeRadiusMinimum                , densityScaleFreeRadiusMaximum
     double precision                            :: densityScaleFreeMinimum                      , densityScaleFreeMaximum
     type            (interpolator), allocatable :: densityScaleFree_
     double precision                            :: massScaleFreeRadiusMinimum                   , massScaleFreeRadiusMaximum
     double precision                            :: massScaleFreeMinimum                         , massScaleFreeMaximum
     type            (interpolator), allocatable :: massScaleFree_
     double precision                            :: angularMomentumSpecificScaleFreeRadiusMinimum, angularMomentumSpecificScaleFreeRadiusMaximum
     double precision                            :: angularMomentumSpecificScaleFreeMinimum      , angularMomentumSpecificScaleFreeMaximum
     type            (interpolator), allocatable :: angularMomentumSpecificScaleFree_
     double precision                            :: timeFreefallScaleFreeRadiusMinimum           , timeFreefallScaleFreeRadiusMaximum
     double precision                            :: timeFreefallScaleFreeMinimum                 , timeFreefallScaleFreeMaximum
     type            (interpolator), allocatable :: timeFreefallScaleFree_
   contains
     !![
     <methods>
       <method method="timeFreefallTabulate" description="Tabulate the freefall time as a function of radius in a scale-free Burkert mass distribution."/>
     </methods>
     !!]
     procedure :: massTotal                         => burkertMassTotal
     procedure :: density                           => burkertDensity
     procedure :: densityGradientRadial             => burkertDensityGradientRadial
     procedure :: densityRadialMoment               => burkertDensityRadialMoment
     procedure :: massEnclosedBySphere              => burkertMassEnclosedBySphere
     procedure :: velocityRotationCurveMaximum      => burkertVelocityRotationCurveMaximum
     procedure :: radiusRotationCurveMaximum        => burkertRadiusRotationCurveMaximum
     procedure :: radiusEnclosingMass               => burkertRadiusEnclosingMass
     procedure :: radiusEnclosingDensity            => burkertRadiusEnclosingDensity
     procedure :: radiusFromSpecificAngularMomentum => burkertRadiusFromSpecificAngularMomentum
     procedure :: fourierTransform                  => burkertFourierTransform
     procedure :: radiusFreefall                    => burkertRadiusFreefall
     procedure :: radiusFreefallIncreaseRate        => burkertRadiusFreefallIncreaseRate
     procedure :: timeFreefallTabulate              => burkertTimeFreefallTabulate
     procedure :: potentialIsAnalytic               => burkertPotentialIsAnalytic
     procedure :: potential                         => burkertPotential
     procedure :: energyPotential                   => burkertEnergyPotential
     procedure :: descriptor                        => burkertDescriptor
  end type massDistributionBurkert
  
  interface massDistributionBurkert
     !!{
     Constructors for the \refClass{massDistributionBurkert} mass distribution class.
     !!}
     module procedure massDistributionBurkertConstructorParameters
     module procedure massDistributionBurkertConstructorInternal
  end interface massDistributionBurkert

  ! The minimum (scale-free) freefall timescale in a Burkert profile.
  double precision , parameter :: timeFreefallScaleFreeMinimum=sqrt(3.0d0*Pi)/4.0d0

contains

  function massDistributionBurkertConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionBurkert} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    type            (massDistributionBurkert)                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    double precision                                         :: mass                , scaleLength  , &
         &                                                      densityNormalization, radiusOuter
    logical                                                  :: dimensionless
    type            (varying_string         )                :: componentType       , massType

    !![
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>1.0d0/Pi/(log(8.0d0)-Pi/2.0d0)</defaultValue>
      <description>The density normalization of the Burkert profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scaleLength</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scale radius of the Burkert profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the Burkert profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusOuter</name>
      <description>The outer radius of the Burkert profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the Burkert profile is considered to be dimensionless.</description>
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
     <call>self=massDistributionBurkert(scaleLength=scaleLength,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="densityNormalization" value="densityNormalization" parameterPresent="parameters"/>
     <argument name="mass"                 value="mass"                 parameterPresent="parameters"/>
     <argument name="radiusOuter"          value="radiusOuter"          parameterPresent="parameters"/>
     <argument name="dimensionless"        value="dimensionless"        parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massDistributionBurkertConstructorParameters

  function massDistributionBurkertConstructorInternal(scaleLength,densityNormalization,mass,radiusOuter,dimensionless,componentType,massType) result(self)
    !!{
    Internal constructor for ``burkert'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionBurkert     )                          :: self
    double precision                              , intent(in   )           :: scaleLength
    double precision                              , intent(in   ), optional :: densityNormalization, radiusOuter, &
         &                                                                     mass
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree
    !![
    <constructorAssign variables="scaleLength, componentType, massType"/>
    !!]

    ! Determine density normalization.
    if      (                                   &
         &   present(densityNormalization)      &
         &  ) then
       self%densityNormalization=densityNormalization
    else if (                                   &
         &   present(mass                ).and. &
         &   present(radiusOuter        )       &
         &  ) then
       radiusScaleFree          =+radiusOuter/self%scaleLength
       self%densityNormalization=+mass/Pi/self%scaleLength**3/(2.0d0*log(1.0d0+radiusScaleFree)+log(1.0d0+radiusScaleFree**2)-2.0d0*atan(radiusScaleFree))
    else
       self%densityNormalization=+0.0d0
       call Error_Report('either "densityNormalization", or "mass" and "radiusOuter" must be specified'//{introspection:location})
    end if
    ! Determine if profile is dimensionless.
    if      (present(dimensionless     )) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! Initialize memoized results.
    self%densityScaleFreeMinimum                      =+huge(0.0d0)
    self%densityScaleFreeMaximum                      =-huge(0.0d0)
    self%densityScaleFreeRadiusMinimum                =+1.0d0
    self%densityScaleFreeRadiusMaximum                =+1.0d0
    self%massScaleFreeMinimum                         =+huge(0.0d0)
    self%massScaleFreeMaximum                         =-huge(0.0d0)
    self%massScaleFreeRadiusMinimum                   =+1.0d0
    self%massScaleFreeRadiusMaximum                   =+1.0d0
    self%angularMomentumSpecificScaleFreeMinimum      =+huge(0.0d0)
    self%angularMomentumSpecificScaleFreeMaximum      =-huge(0.0d0)
    self%angularMomentumSpecificScaleFreeRadiusMinimum=+1.0d0
    self%angularMomentumSpecificScaleFreeRadiusMaximum=+1.0d0
    self%timeFreefallScaleFreeMinimum                 =+huge(0.0d0)
    self%timeFreefallScaleFreeMaximum                 =-huge(0.0d0)
    self%timeFreefallScaleFreeRadiusMinimum           =+1.0d0
    self%timeFreefallScaleFreeRadiusMaximum           =+1.0d0
    return
  end function massDistributionBurkertConstructorInternal

  double precision function burkertMassTotal(self)
    !!{
    Return the total mass in an Burkert mass distribution.
    !!}
    implicit none
    class(massDistributionBurkert), intent(inout) :: self
 
    burkertMassTotal=huge(0.0d0)
    return
  end function burkertMassTotal

  double precision function burkertDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an Burkert mass distribution.
    !!}
    implicit none
    class           (massDistributionBurkert), intent(inout) :: self
    class           (coordinate             ), intent(in   ) :: coordinates
    double precision                                         :: radiusScaleFree

    ! Compute the density at this position.
    radiusScaleFree=+coordinates%rSpherical          () &
         &          /self       %scaleLength
    burkertDensity =+self       %densityNormalization   &
         &          /(+1.0d0+radiusScaleFree   )        &
         &          /(+1.0d0+radiusScaleFree**2)
    return
  end function burkertDensity

  double precision function burkertDensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an Burkert \citep{burkert_structure_1995} mass distribution.
    !!}
    implicit none
    class           (massDistributionBurkert), intent(inout), target   :: self
    class           (coordinate             ), intent(in   )           :: coordinates
    logical                                  , intent(in   ), optional :: logarithmic
    double precision                                                   :: radiusScaleFree
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    radiusScaleFree      =+coordinates%rSpherical()         &
         &                /self       %scaleLength
    densityGradientRadial=-3.0d0                            &
         &                +1.0d0/(1.0d0+radiusScaleFree   ) &
         &                +2.0d0/(1.0d0+radiusScaleFree**2)
    if (.not.logarithmic_) densityGradientRadial=+            densityGradientRadial              &
         &                                       *self       %density              (coordinates) &
         &                                       /coordinates%rSpherical           (           )
    return
  end function burkertDensityGradientRadial

  double precision function burkertDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Computes radial moments of the density in an Burkert \citep{burkert_structure_1995} mass distribution.
    !!}
    implicit none
    class           (massDistributionBurkert), intent(inout)           :: self
    double precision                         , intent(in   )           :: moment
    double precision                         , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    logical                                  , intent(  out), optional :: isInfinite
    double precision                                                   :: radialMomentMinimum, radialMomentMaximum

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
      !!{
      Provides the scale-free part of the radial moment of the Burkert density profile.
      !!}
      use :: Hypergeometric_Functions, only : Hypergeometric_2F1
      use :: Numerical_Comparison    , only : Values_Agree
      implicit none
      double precision, intent(in   ) :: radius
      
      if (Values_Agree(moment,0.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=+0.50d0*atan(       radius   ) &
              &                +0.50d0*log (+1.0d0+radius   ) &
              &                -0.25d0*log (+1.0d0+radius**2)
      else if (Values_Agree(moment,1.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=+0.50d0*atan(       radius   ) &
              &                -0.50d0*log (+1.0d0+radius   ) &
              &                +0.25d0*log (+1.0d0+radius**2)
      else if (Values_Agree(moment,2.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=-0.50d0*atan(       radius   ) &
              &                +0.50d0*log (+1.0d0+radius   ) &
              &                +0.25d0*log (+1.0d0+radius**2)
      else if (Values_Agree(moment,3.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=+                   radius     &
              &                -0.50d0*atan(       radius   ) &
              &                -0.50d0*log (+1.0d0+radius   ) &
              &                -0.25d0*log (+1.0d0+radius**2)
      else
         radialMomentScaleFree=+radius**(1.0d0+moment)                                                                              &
              &                /   2.0d0                                                                                            &
              &                /  (1.0d0+moment)                                                                                    &
              &                /  (2.0d0+moment)                                                                                    &
              &                *(                                                                                                   &
              &                  +(2.0d0+moment)*Hypergeometric_2F1([1.0d0,0.5d0*(1.0d0+moment)],[0.5d0*(3.0d0+moment)],-radius**2) &
              &                  +(2.0d0+moment)*Hypergeometric_2F1([1.0d0,      (1.0d0+moment)],[      (2.0d0+moment)],-radius   ) &
              &                  +(1.0d0+moment)*Hypergeometric_2F1([1.0d0,0.5d0*(2.0d0+moment)],[0.5d0*(4.0d0+moment)],-radius**2) &
              &                 )
      end if
      return
    end function radialMomentScaleFree

  end function burkertDensityRadialMoment

  double precision function burkertMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for burkert mass distributions.
    !!}
    implicit none
    class           (massDistributionBurkert), intent(inout), target :: self
    double precision                         , intent(in   )         :: radius
    double precision                                                 :: radiusScaleFree
    
    radiusScaleFree=+      radius                              &
         &           /self%scaleLength
    mass           =+self%densityNormalization                 &
         &          *self%scaleLength                      **3 &
         &          *massEnclosedScaleFree(radiusScaleFree)
    return
  end function burkertMassEnclosedBySphere
  
  double precision function burkertRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for burkert mass distributions.
    !!}    
    use :: Numerical_Ranges, only : Make_Range  , rangeTypeLogarithmic
    use :: Error           , only : Error_Report
    implicit none
    class           (massDistributionBurkert), intent(inout), target       :: self
    double precision                         , intent(in   ), optional     :: mass                       , massFractional
    double precision                         , allocatable  , dimension(:) :: radii                      , masses
    double precision                         , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                       :: massScaleFree              , mass_
    integer                                                                :: countRadii

    mass_=0.0d0
    if (present(mass)) then
       mass_=mass
    else if (present(massFractional)) then
       call Error_Report('mass is unbounded, so mass fraction is undefined'//{introspection:location})
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
       do while (massEnclosedScaleFree(self%massScaleFreeRadiusMinimum) >= massScaleFree)
          self%massScaleFreeRadiusMinimum=0.5d0*self%massScaleFreeRadiusMinimum
       end do
       do while (massEnclosedScaleFree(self%massScaleFreeRadiusMaximum) <  massScaleFree)
          self%massScaleFreeRadiusMaximum=2.0d0*self%massScaleFreeRadiusMaximum
       end do
       countRadii=int(log10(self%massScaleFreeRadiusMaximum/self%massScaleFreeRadiusMinimum)*countRadiiPerDecade)+1
       if (allocated(self%massScaleFree_)) deallocate(self%massScaleFree_)
       allocate(     radii         (countRadii))
       allocate(     masses        (countRadii))
       allocate(self%massScaleFree_            )
       radii                    =Make_Range(self%massScaleFreeRadiusMinimum,self%massScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
       masses                   =massEnclosedScaleFree(           radii)
       self%massScaleFreeMinimum=masses               (         1      )
       self%massScaleFreeMaximum=masses               (countRadii      )
       self%massScaleFree_      =interpolator         (masses    ,radii)
    end if
    radius=+self%massScaleFree_%interpolate(massScaleFree) &
         & *self%scaleLength
    return
  end function burkertRadiusEnclosingMass
  
  double precision function burkertRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mean density for burkert mass distributions.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionBurkert), intent(inout), target       :: self
    double precision                         , intent(in   )               :: density
    double precision                         , intent(in   ), optional     :: radiusGuess
    double precision                         , allocatable  , dimension(:) :: radii                      , densities
    double precision                         , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                       :: densityScaleFree
    integer                                                                :: countRadii

    densityScaleFree=+density                   &
         &           /self%densityNormalization
    if     (                                                  &
         &   densityScaleFree <= self%densityScaleFreeMinimum &
         &  .or.                                              &
         &   densityScaleFree >  self%densityScaleFreeMaximum &
         & ) then
       do while (densityEnclosedScaleFree(self%densityScaleFreeRadiusMinimum) <  densityScaleFree)
          self%densityScaleFreeRadiusMinimum=0.5d0*self%densityScaleFreeRadiusMinimum
       end do
       do while (densityEnclosedScaleFree(self%densityScaleFreeRadiusMaximum) >= densityScaleFree)
          self%densityScaleFreeRadiusMaximum=2.0d0*self%densityScaleFreeRadiusMaximum
       end do
       countRadii=int(log10(self%densityScaleFreeRadiusMaximum/self%densityScaleFreeRadiusMinimum)*countRadiiPerDecade)+1
       if (allocated(self%densityScaleFree_)) deallocate(self%densityScaleFree_)
       allocate(     radii            (countRadii))
       allocate(     densities        (countRadii))
       allocate(self%densityScaleFree_            )
       radii                       = Make_Range(self%densityScaleFreeRadiusMinimum,self%densityScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
       densities                   =-densityEnclosedScaleFree(           radii)
       self%densityScaleFreeMinimum=-densities               (countRadii      )
       self%densityScaleFreeMaximum=-densities               (         1      )
       self%densityScaleFree_      = interpolator            (densities ,radii)
    end if
    radius=+self%densityScaleFree_%interpolate(-densityScaleFree) &
         & *self%scaleLength
    return    
  end function burkertRadiusEnclosingDensity

  elemental double precision function massEnclosedScaleFree(radius) result(mass)
    !!{
    Evaluate the mass enclosed by a given radius in a scale-free Burkert mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    double precision, parameter     :: minimumRadiusForExactSolution=1.0d-4

    if (radius < minimumRadiusForExactSolution) then
       ! Use a series solution for small radii.
       mass   =+4.0d0             &
            &  /3.0d0             &
            &  *Pi                &
            &  *        radius**3 &
            &  *(+1.0d0-radius)
    else
       ! Use the exact solution.
       mass   =+Pi                            &
            &  *(                             &
            &    +2.0d0*log (1.0d0+radius   ) &
            &    +      log (1.0d0+radius**2) &
            &    -2.0d0*atan(      radius   ) &
            &   )
    end if
    return
  end function massEnclosedScaleFree

  elemental double precision function densityEnclosedScaleFree(radius) result(density)
    !!{
    Evaluate the mean enclosed density at a given radius in a scale-free Burkert mass distribution.
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
  
  double precision function burkertRadiusFromSpecificAngularMomentum(self,angularMomentumSpecific) result(radius)
    !!{
    Computes the radius corresponding to a given specific angular momentum for burkert mass distributions.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Ranges                , only : Make_Range                    , rangeTypeLogarithmic
    implicit none
    class           (massDistributionBurkert), intent(inout), target       :: self
    double precision                         , intent(in   )               :: angularMomentumSpecific
    double precision                         , allocatable  , dimension(:) :: radii                                   , angularMomentaSpecific
    double precision                         , parameter                   :: countRadiiPerDecade             =100.0d0
    double precision                                                       :: angularMomentumSpecificScaleFree
    integer                                                                :: countRadii

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
          do while (angularMomentumSpecificEnclosedScaleFree(self%angularMomentumSpecificScaleFreeRadiusMinimum) >= angularMomentumSpecificScaleFree)
             self%angularMomentumSpecificScaleFreeRadiusMinimum=0.5d0*self%angularMomentumSpecificScaleFreeRadiusMinimum
          end do
          do while (angularMomentumSpecificEnclosedScaleFree(self%angularMomentumSpecificScaleFreeRadiusMaximum) <  angularMomentumSpecificScaleFree)
             self%angularMomentumSpecificScaleFreeRadiusMaximum=2.0d0*self%angularMomentumSpecificScaleFreeRadiusMaximum
          end do
          countRadii=int(log10(self%angularMomentumSpecificScaleFreeRadiusMaximum/self%angularMomentumSpecificScaleFreeRadiusMinimum)*countRadiiPerDecade)+1
          if (allocated(self%angularMomentumSpecificScaleFree_)) deallocate(self%angularMomentumSpecificScaleFree_)
          allocate(     radii                            (countRadii))
          allocate(     angularMomentaSpecific           (countRadii))
          allocate(self%angularMomentumSpecificScaleFree_            )
          radii                                       =Make_Range(self%angularMomentumSpecificScaleFreeRadiusMinimum,self%angularMomentumSpecificScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
          angularMomentaSpecific                      =angularMomentumSpecificEnclosedScaleFree(                       radii)
          self%angularMomentumSpecificScaleFreeMinimum=angularMomentaSpecific                  (                     1      )
          self%angularMomentumSpecificScaleFreeMaximum=angularMomentaSpecific                  (            countRadii      )
          self%angularMomentumSpecificScaleFree_      =interpolator                            (angularMomentaSpecific,radii)
       end if
       radius=+self%angularMomentumSpecificScaleFree_%interpolate(angularMomentumSpecificScaleFree) &
            & *self%scaleLength
    else
       radius=+0.0d0
    end if
    return    
  end function burkertRadiusFromSpecificAngularMomentum

  elemental double precision function angularMomentumSpecificEnclosedScaleFree(radius) result(angularMomentumSpecific)
    !!{
    Evaluate the specific angular momentum at a given radius in a scale-free Burkert mass distribution.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    angularMomentumSpecific=+sqrt(                               &
         &                        +massEnclosedScaleFree(radius) &
         &                        *                      radius  &
         &                       )
    return
  end function angularMomentumSpecificEnclosedScaleFree

  double precision function burkertVelocityRotationCurveMaximum(self) result(velocity)
    !!{
    Return the peak velocity in the rotation curve for an burkert mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionBurkert ), intent(inout) :: self
    double precision                          , parameter     :: circularVelocityMaximumScaleFree=1.644297750532498d0 ! The circular velocity (in scale-free units) at the peak of the Burkert rotation curve.
    !                                                                                                                   Numerical value found using Mathematica.

    velocity=+circularVelocityMaximumScaleFree            &
         &   *sqrt(                                       &
         &         +self%densityNormalization             &
         &        )                                       &
         &   *      self%scaleLength
    if (.not.self%isDimensionless())                      &
         & velocity=+velocity                             &
         &          *sqrt(gravitationalConstant_internal)
    return
  end function burkertVelocityRotationCurveMaximum

  double precision function burkertRadiusRotationCurveMaximum(self) result(radius)
    !!{
    Return the peak velocity in the rotation curve for an burkert mass distribution.
    !!}
    implicit none
    class           (massDistributionBurkert), intent(inout), target :: self
    ! The radius (in scale-free units) at the peak of the Burkert rotation curve. Numerical value found using Mathematica.
    double precision                         , parameter             :: radiusCircularVelocityMaximumScaleFree=3.244625724604264d0
    
    radius=+radiusCircularVelocityMaximumScaleFree &
         & *self%scaleLength
    return
  end function burkertRadiusRotationCurveMaximum

  logical function burkertPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionBurkert), intent(inout) :: self

    isAnalytic=.true.
    return
  end function burkertPotentialIsAnalytic

  double precision function burkertPotential(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in an burkert mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess     , structureErrorCodeInfinite
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Error                           , only : Error_Report
    implicit none
    class           (massDistributionBurkert          ), intent(inout), target   :: self
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
  end function burkertPotential

  elemental double precision function potentialScaleFree(radius) result(potential)
    !!{
    Compute the potential in a scale-free Burkert mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    double precision, parameter     :: radiusSmall=1.0d-3

    if (radius < radiusSmall) then
       ! Use a series solution for very small radii.
       potential=-             Pi**2           &
            &    +2.0d0/ 3.0d0*Pi   *radius**2 &
            &    -1.0d0/ 3.0d0*Pi   *radius**3 &
            &    +2.0d0/21.0d0*Pi   *radius**6
    else
       ! Use the full expression for larger radii.
       potential=-Pi                                           &
            &    /radius                                       & 
            &    *(                                            &
            &      +2.0d0*(1.0d0+radius)*log (1.0d0+radius   ) &
            &      +      (1.0d0-radius)*log (1.0d0+radius**2) &
            &      -2.0d0               *atan(      radius   ) &
            &      +2.0d0*       radius *atan(1.0d0/radius   ) &
            &     )
    end if
    return
  end function potentialScaleFree

 double precision function potentialDifferenceScaleFree(radius1,radius2) result(potential)
    !!{
    Compute the potential difference in a scale-free Burkert mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Comparison    , only : Values_Agree
    implicit none
    double precision, intent(in   ) :: radius1                            , radius2
    double precision, parameter     :: radiusSmall                 =1.0d-3
    double precision, parameter     :: toleranceRelative           =1.0d-3
    double precision                :: potentialGradientLogarithmic       , radiusDifferenceLogarithmic
    
    if (Values_Agree(radius1,radius2,relTol=toleranceRelative) .or. max(radius1,radius2) < radiusSmall) then
       if (radius1 < radiusSmall) then
          potentialGradientLogarithmic=-4.0d0/3.0d0/Pi   *radius1**2 &
               &                       +1.0d0      /Pi   *radius1**3 &
               &                       -8.0d0/9.0d0/Pi**2*radius1**4
       else
          potentialGradientLogarithmic=-(                                              &
               &                         +                      log (1.0d0+radius1**2) &
               &                         +2.0d0                *log (1.0d0+radius1   ) &
               &                         -2.0d0                *atan(      radius1   ) &
               &                        )                                              &
               &                       /(                                              &
               &                         +      (1.0d0-radius1)*log (1.0d0+radius1**2) &
               &                         +2.0d0*(1.0d0+radius1)*log (1.0d0+radius1   ) &
               &                         -2.0d0                *atan(      radius1   ) &
               &                         +2.0d0*       radius1 *atan(1.0d0/radius1   ) &
               &                        )
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
  
  double precision function burkertFourierTransform(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in an Burkert mass
    distribution.
    !!}
    use :: Exponential_Integrals   , only : Exponential_Integral
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionBurkert), intent(inout) :: self
    double precision                         , intent(in   ) :: radiusOuter        , wavenumber
    double precision                                         :: wavenumberScaleFree, radiusOuterScaleFree

    waveNumberScaleFree =+waveNumber *self%scaleLength
    radiusOuterScaleFree=+radiusOuter/self%scaleLength
    fourierTransform    =+dimag(                                                                                                                                                                                                                                                    &
         &                      +dcmplx(1.0d0,1.0d0)                                                                                                                                                                                                                                &
         &                      *Pi                                                                                                                                                                                                                                                 &
         &                      /wavenumberScaleFree                                                                                                                                                                                                                                &
         &                      *(                                                                                                                                                                                                                                                  &
         &                        +                     exp(+                    wavenumberScaleFree)*(dcmplx(0.0d0,-1.0d0)*Pi-Exponential_Integral(-                    wavenumberScaleFree)+Exponential_Integral(dcmplx(-1.0d0,      +radiusOuterScaleFree)*wavenumberScaleFree)) &
         &                        +dcmplx(1.0d0,-1.0d0)*exp(-dcmplx(0.0d0,1.0d0)*wavenumberScaleFree)*(                       +Exponential_Integral(+dcmplx(0.0d0,1.0d0)*wavenumberScaleFree)-Exponential_Integral(dcmplx(+0.0d0,+1.0d0+radiusOuterScaleFree)*wavenumberScaleFree)) &
         &                        +dcmplx(0.0d0,+1.0d0)*exp(-                    wavenumberScaleFree)*(                       +Exponential_Integral(+                    wavenumberScaleFree)-Exponential_Integral(dcmplx(+1.0d0,      +radiusOuterScaleFree)*wavenumberScaleFree)) &
         &                       )                                                                                                                                                                                                                                                  &
         &                     )                                                                                                                                                                                                                                                    &
         &               /massEnclosedScaleFree(radiusOuterScaleFree)
    return
  end function burkertFourierTransform
  
  double precision function burkertRadiusFreefall(self,time) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in an Burkert mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr, gravitationalConstant_internal
    implicit none
    class           (massDistributionBurkert), intent(inout) :: self
    double precision                         , intent(in   ) :: time
    double precision                                         :: timeScaleFree, timeScale
    
    timeScale    =+1.0d0/sqrt(                                &
         &                    +gravitationalConstant_internal &
         &                    *self%densityNormalization      &
         &                   )                                &
         &        *MpcPerKmPerSToGyr
    timeScaleFree=+time                                       &
         &        /timeScale
    if (timeScaleFree <= timeFreefallScaleFreeMinimum) then
       radius=0.0d0
       return
    end if
    call self%timeFreefallTabulate(timeScaleFree)
    radius=+self%timeFreefallScaleFree_%interpolate(timeScaleFree) &
         & *self%scaleLength
    return   
  end function burkertRadiusFreefall
  
  double precision function burkertRadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in an burkert mass
    distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr, gravitationalConstant_internal
    implicit none
    class           (massDistributionBurkert), intent(inout) :: self
    double precision                         , intent(in   ) :: time
    double precision                                         :: timeScaleFree, timeScale

    timeScale    =+1.0d0/sqrt(                                &
         &                    +gravitationalConstant_internal &
         &                    *self%densityNormalization      &
         &                   )                                &
         &        *MpcPerKmPerSToGyr
    timeScaleFree=+time                                       &
         &        /timeScale
    if (timeScaleFree <= timeFreefallScaleFreeMinimum) then
       radiusIncreaseRate=0.0d0
       return
    end if
    call self%timeFreefallTabulate(timeScaleFree)
    radiusIncreaseRate=+self%timeFreefallScaleFree_%derivative(timeScaleFree) &
         &             *self%scaleLength                                      &
         &             /     timeScale
    return
  end function burkertRadiusFreefallIncreaseRate
  
  subroutine burkertTimeFreefallTabulate(self,timeScaleFree)
    !!{
    Tabulate the freefall radius at the given {\normalfont \ttfamily time} in an Burkert mass distribution.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges     , only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionBurkert), intent(inout)               :: self
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
      Evaluate the freefall time from a given radius in a scale-free Burkert mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      radiusStart          =                            radius
      timeFreefallScaleFree=integrator_%integrate(0.0d0,radius)
      return
    end function timeFreefallScaleFree
    
    double precision function timeFreeFallIntegrand(radius)
      !!{
      Integrand used to find the freefall time in a scale-free Burkert mass distribution.
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
    
  end subroutine burkertTimeFreefallTabulate

  double precision function burkertEnergyPotential(self,radiusOuter) result(energy)
    !!{
    Compute the potential energy within a given {\normalfont \ttfamily radius} in a Burkert mass distribution. This is
    \begin{eqnarray}
      W &=& \frac{1}{24} \pi ^2 \left(48 \mathrm{G}G + 48 i \text{Li}_2\left(\left(\frac{1}{2}+\frac{i}{2}\right) (x+1)\right)-48 i \text{Li}_2\left(\left(\frac{1}{2}-\frac{i}{2}\right) (x+1)\right)+48 i \text{Li}_2\left(\frac{i+x}{-i+x}\right)-96 \text{Li}_2\left(\left(-\frac{1}{2}+\frac{i}{2}\right) (-i+x)\right)-96 \text{Li}_2\left(\left(-\frac{1}{2}-\frac{i}{2}\right) (i+x)\right)-48 i \text{Li}_2\left(i \exp(2 i \tan ^{-1}(x))\right)+12 \left(\log ^2\left(x^2+1\right)+4 \log (x) \log \left(x^2+1\right)-4 \log (x+1) \log \left(x^2+1\right)+(2+4 i) \pi  \log \left(x^2+1\right)+\tan ^{-1}(x) \left(4 \log \left(x^2+1\right)+8 \log \left(-\frac{2 i}{x-i}\right)+8 \log \left(1-i \exp(2 i \tan ^{-1}(x))\right)+2 i \pi \right)-4 \log ^2(x+1)-4 \log (x-i) \log ((1-i) (x+1))-4 \log (x+i) \log ((1+i) (x+1))+\log (64) \log (x-i)-4 \log (x) \log (x-i)+4 \log (x+1) \log (x-i)-3 i \pi  \log (x-i)+\log (64) \log (x+i)-4 \log (x) \log (x+i)+4 \log (x+1) \log (x+i)-5 i \pi  \log (x+i)+4 i \log (x+1) \log ((-1-i) (x+i))-4 \tan ^{-1}(x)^2+4 \pi  \log \left(1+\exp(-2 i \tan ^{-1}(x))\right)+2 \pi  \log \left(1-i \exp(2 i \tan ^{-1}(x))\right)-2 \pi  \log \left(\sin \left(\tan ^{-1}(x)+\frac{\pi }{4}\right)\right)-\log (2) (7 \pi +\log (4))\right)-48 i \log ((1+i)-(1-i) x) \log (x+1)-\pi ^2 (14-9 i)\right)
    \end{eqnarray}
    where $x=r/r_\mathrm{s}$ and $\mathrm{G}$ is Catalan's constant.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi                            , catalan
    use :: Dilogarithms                    , only : Dilogarithm
    implicit none
    class           (massDistributionBurkert), intent(inout) :: self
    double precision                         , intent(in   ) :: radiusOuter
    double precision                                         :: radiusOuterScaleFree
    
    radiusOuterScaleFree=+     radiusOuter                                                                                                                                                                                         &
         &               /self%scaleLength
    energy              =real(                                                                                                                                                                                                     &
         &                    -gravitationalConstant_internal                                                                                                                                                                      &
         &                    *self%scaleLength               **5                                                                                                                                                                  &
         &                    *self%densityNormalization      **2                                                                                                                                                                  &
         &                    *(                                                                                                                                                                                                   &
         &                      +Pi**2                                                                                                                                                                                             &
         &                      *(                                                                                                                                                                                                 &
         &                        +       48.0d0         *catalan                                                                                                                                                                  &
         &                        -dcmplx(14.0d0,- 9.0d0)*Pi**2                                                                                                                                                                    &
         &                        -dcmplx( 0.0d0,+48.0d0)*log(dcmplx(1.0d0,1.0d0)-dcmplx(1.0d0,-1.0d0)*radiusOuterScaleFree)*log(1.0d0+radiusOuterScaleFree)                                                                       &
         &                        +12.0d0                                                                                                                                                                                          &
         &                        *(                                                                                                                                                                                               &
         &                          -4.0d0*atan(radiusOuterScaleFree)**2                                                                                                                                                           &
         &                          - log(2.0d0)*(7.0d0*Pi + log(4.0d0))                                                                                                                                                           &
         &                          +4.0d0*Pi*log(1.0d0+                    exp(-2.0d0*dcmplx(0.0d0,1.0d0)*atan(radiusOuterScaleFree)))                                                                                            &
         &                          +2.0d0*Pi*log(1.0d0-dcmplx(0.0d0,1.0d0)*exp(+2.0d0*dcmplx(0.0d0,1.0d0)*atan(radiusOuterScaleFree)))                                                                                            &
         &                          - dcmplx(0.0d0,3.0d0)*Pi*log(dcmplx(0.0d0,-1.0d0) +radiusOuterScaleFree)                                                                                                                       &
         &                          +log   (64.0d0      )   *log(dcmplx(+0.0d0,-1.0d0)+radiusOuterScaleFree)-log(dcmplx(+0.0d0,-1.0d0)       +radiusOuterScaleFree) *       4.0d0          *log(      radiusOuterScaleFree   )     &
         &                          -dcmplx( 0.0d0,5.0d0)*Pi*log(dcmplx(+0.0d0,+1.0d0)+radiusOuterScaleFree)                                                                                                                       &
         &                          +log   (64.0d0      )   *log(dcmplx(+0.0d0,+1.0d0)+radiusOuterScaleFree)                                                                                                                       &
         &                          -              4.0d0     *log(                     +radiusOuterScaleFree)*log(dcmplx(+0.0d0,+1.0d0)       +radiusOuterScaleFree) +       4.0d0                                                 &
         &                                                                                                   *log(dcmplx(+0.0d0,-1.0d0)       +radiusOuterScaleFree)                        *log(1.0d0+radiusOuterScaleFree   )    &
         &                          +dcmplx( 0.0d0,4.0d0)    *log(dcmplx(-1.0d0,-1.0d0)                      *   (dcmplx(+0.0d0,+1.0d0)       +radiusOuterScaleFree))                       *log(1.0d0+radiusOuterScaleFree   )    &
         &                          +              4.0d0     *log(dcmplx(+0.0d0,+1.0d0)+radiusOuterScaleFree)*log(                       1.0d0+radiusOuterScaleFree )-       4.0d0          *log(1.0d0+radiusOuterScaleFree   )**2 &
         &                          -              4.0d0     *log(dcmplx(+0.0d0,-1.0d0)+radiusOuterScaleFree)*log(dcmplx(+1.0d0,-1.0d0)*(1.0d0+radiusOuterScaleFree))                                                              &
         &                          -              4.0d0     *log(dcmplx(+0.0d0,+1.0d0)+radiusOuterScaleFree)*log(dcmplx(+1.0d0,+1.0d0)*(1.0d0+radiusOuterScaleFree))+dcmplx(2.0d0,4.0d0)*Pi*log(1.0d0+radiusOuterScaleFree**2)    &
         &                          +              4.0d0     *log(                      radiusOuterScaleFree)                                                                               *log(1.0d0+radiusOuterScaleFree**2)    &
         &                          -              4.0d0     *log(       +1.0d0        +radiusOuterScaleFree)                                                                               *log(1.0d0+radiusOuterScaleFree**2)    &
         &                          +                                                                                                                                                        log(1.0d0+radiusOuterScaleFree**2)**2 &
         &                          +                                                                                          atan(radiusOuterScaleFree)                                                                          &
         &                          *(dcmplx(0.0d0,2.0d0)*Pi + 8.0d0*log(1 - dcmplx(0.0d0,1.0d0)*exp(2.0d0*dcmplx(0.0d0,1.0d0)*atan(radiusOuterScaleFree)))                                                                        &
         &                          +8.0d0*log(dcmplx(0.0d0,-2.0d0)/(dcmplx(0.0d0,-1.0d0)+radiusOuterScaleFree   ))                                                                                                                &
         &                          +4.0d0*log(       1.0d0                             +radiusOuterScaleFree**2))                                                                                                                 &
         &                          -2.0d0*Pi*log(sin(Pi/4.0d0+atan(radiusOuterScaleFree)))                                                                                                                                        &
         &                         )                                                                                                                                                                                               &
         &                        -dcmplx(0.0d0,48.0d0)*Dilogarithm(dcmplx(0.0d0,1.0d0)*exp(2.0d0*dcmplx(0.0d0,1.0d0)*atan(radiusOuterScaleFree)))                                                                                 &
         &                        -             96.0d0 *Dilogarithm( dcmplx(-0.5d0,+0.5d0)                      *(dcmplx(0.0d0,-1.0d0)+radiusOuterScaleFree))                                                                      &
         &                        -             96.0d0 *Dilogarithm( dcmplx(-0.5d0,-0.5d0)                      *(dcmplx(0.0d0,+1.0d0)+radiusOuterScaleFree))                                                                      &
         &                        +dcmplx(0.0d0,48.0d0)*Dilogarithm((dcmplx(+0.0d0,+1.0d0)+radiusOuterScaleFree)/(dcmplx(0.0d0,-1.0d0)+radiusOuterScaleFree))                                                                      &
         &                        -dcmplx(0.0d0,48.0d0)*Dilogarithm( dcmplx(+0.5d0,-0.5d0)                      *(             +1.0d0 +radiusOuterScaleFree))                                                                      &
         &                        +dcmplx(0.0d0,48.0d0)*Dilogarithm( dcmplx(+0.5d0,+0.5d0)                      *(             +1.0d0 +radiusOuterScaleFree))                                                                      &
         &                       )                                                                                                                                                                                                 &
         &                     )                                                                                                                                                                                                   &
         &                    /24.d0                                                                                                                                                                                               &
         &                   )
         return
  end function burkertEnergyPotential
  
  subroutine burkertDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionBurkert), intent(inout)           :: self
    type     (inputParameters    ), intent(inout)           :: descriptor
    logical                       , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18             )                          :: parameterLabel
    type     (inputParameters    )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','Burkert')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%scaleLength
    call parameters%addParameter('scaleLength'         ,trim(adjustl(parameterLabel)))
    return
  end subroutine burkertDescriptor

