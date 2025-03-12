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
  Implementation of the \cite{zhao_analytical_1996} mass distribution class.
  !!}

  use :: Numerical_Interpolation , only : interpolator

  !![
  <enumeration>
   <name>specialCase</name>
   <description>Special cases for {\normalfont \ttfamily zhao1996} dark matter halo profile class.</description>
   <entry label="general"    />
   <entry label="coredNFW"   />
   <entry label="gamma0_5NFW"/>
   <entry label="NFW"        />
   <entry label="gamma1_5NFW"/>
  </enumeration>
  !!]
  
  !![
  <massDistribution name="massDistributionZhao1996">
    <description>
    A mass distribution class which implements the \citep{zhao_analytical_1996} density profile:
    \begin{equation}
      \rho_\mathrm{dark matter}(r) = \rho_0 \left({r\over r_\mathrm{s}}\right)^{-\gamma} \left(1+[{r\over r_\mathrm{s}}]^\alpha\right)^{-(\beta-\gamma)/\alpha}.
    \end{equation}
    The mass enclosed within radius $r$ is given by
    \begin{equation}
    M(&lt;r) = \frac{4 \pi}{3-\gamma} \rho_0 r_\mathrm{s}^{3-\gamma} {}_2F_1\left[\left(\frac{3-\gamma}{\alpha}\right),\left(\frac{-\beta+\gamma}{\alpha},1+\frac{3-\gamma}{\alpha}\right),-r^\alpha\right]
    \end{equation}
    where $R=r/r_\mathrm{s}$. The associated gravitational potential is
    \begin{equation}
    \Phi(r) = - \frac{4 \pi \mathrm{G}}{-3+\gamma} \rho_0 r^{2-\gamma} \frac{\Gamma[(3+\alpha-\gamma)/\alpha]}{\Gamma[(3-\gamma)/\alpha]} \left( \Gamma\left[\frac{2-\gamma}{\alpha}\right] {}_p\tilde{F}F_q\left[\left\{\frac{2-\gamma}{\alpha},\frac{\beta-\gamma}{\alpha}\right\},\left\{\frac{2+\alpha-\gamma}{\alpha}\right\},-r\alpha\right] -  \Gamma\left[\frac{3-\gamma}{\alpha}\right] {}_p\tilde{F}F_q\left[\left\{\frac{3-\gamma}{\alpha},\frac{\beta-\gamma}{\alpha}\right\},\left\{\frac{3+\alpha-\gamma}{\alpha}\right\},-r\alpha\right] \right).
    \end{equation}
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionZhao1996
     !!{
     The \citep{zhao_analytical_1996} mass distribution.
     !!}
     private
     type            (enumerationSpecialCaseType)              :: specialCase
     double precision                                          :: densityNormalization                         , scaleLength                                  , &
          &                                                       alpha                                        , beta                                         , &
          &                                                       gamma
     double precision                                          :: densityScaleFreeRadiusMinimum                , densityScaleFreeRadiusMaximum
     double precision                                          :: densityScaleFreeMinimum                      , densityScaleFreeMaximum
     type            (interpolator              ), allocatable :: densityScaleFree_
     double precision                                          :: massScaleFreeRadiusMinimum                   , massScaleFreeRadiusMaximum
     double precision                                          :: massScaleFreeMinimum                         , massScaleFreeMaximum
     type            (interpolator              ), allocatable :: massScaleFree_
     double precision                                          :: angularMomentumSpecificScaleFreeRadiusMinimum, angularMomentumSpecificScaleFreeRadiusMaximum
     double precision                                          :: angularMomentumSpecificScaleFreeMinimum      , angularMomentumSpecificScaleFreeMaximum
     type            (interpolator              ), allocatable :: angularMomentumSpecificScaleFree_
     double precision                                          :: timeFreefallScaleFreeRadiusMinimum           , timeFreefallScaleFreeRadiusMaximum
     double precision                                          :: timeFreefallScaleFreeMinimum                 , timeFreefallScaleFreeMaximum
     type            (interpolator              ), allocatable :: timeFreefallScaleFree_
   contains
     !![
     <methods>
       <method method="timeFreefallTabulate" description="Tabulate the freefall time as a function of radius in a scale-free Zhao1996 mass distribution."/>
     </methods>
     !!]
     procedure :: massTotal                         => zhao1996MassTotal
     procedure :: density                           => zhao1996Density
     procedure :: densityGradientRadial             => zhao1996DensityGradientRadial
     procedure :: densityRadialMoment               => zhao1996DensityRadialMoment
     procedure :: massEnclosedBySphere              => zhao1996MassEnclosedBySphere
     procedure :: velocityRotationCurveMaximum      => zhao1996VelocityRotationCurveMaximum
     procedure :: radiusRotationCurveMaximum        => zhao1996RadiusRotationCurveMaximum
     procedure :: radiusEnclosingMass               => zhao1996RadiusEnclosingMass
     procedure :: radiusEnclosingDensity            => zhao1996RadiusEnclosingDensity
     procedure :: radiusFromSpecificAngularMomentum => zhao1996RadiusFromSpecificAngularMomentum
     procedure :: radiusFreefall                    => zhao1996RadiusFreefall
     procedure :: radiusFreefallIncreaseRate        => zhao1996RadiusFreefallIncreaseRate
     procedure :: timeFreefallTabulate              => zhao1996TimeFreefallTabulate
     procedure :: potentialIsAnalytic               => zhao1996PotentialIsAnalytic
     procedure :: potential                         => zhao1996Potential
     procedure :: fourierTransform                  => zhao1996FourierTransform
     procedure :: energyPotential                   => zhao1996EnergyPotential
     procedure :: energyKinetic                     => zhao1996EnergyKinetic
     procedure :: descriptor                        => zhao1996Descriptor
  end type massDistributionZhao1996
  
  interface massDistributionZhao1996
     !!{
     Constructors for the {\normalfont \ttfamily zhao1996} mass distribution class.
     !!}
     module procedure massDistributionZhao1996ConstructorParameters
     module procedure massDistributionZhao1996ConstructorInternal
  end interface massDistributionZhao1996

  class(massDistributionZhao1996), pointer :: self_
  !$omp threadprivate(self_)

  ! The minimum (scale-free) freefall timescale in a cored NFW profile.
  double precision , parameter :: timeFreefallScaleFreeMinimumCoredNFW=sqrt(3.0d0*Pi)/4.0d0
  
contains

  function massDistributionZhao1996ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily zhao1996} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    use :: Numerical_Constants_Math  , only : Pi
    use :: Hypergeometric_Functions  , only : Hypergeometric_2F1
    implicit none
    type            (massDistributionZhao1996)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    double precision                                          :: mass                , scaleLength, &
         &                                                       densityNormalization, radiusOuter, &
         &                                                       alpha               , beta       , &
         &                                                       gamma
    logical                                                   :: dimensionless
    type            (varying_string          )                :: componentType
    type            (varying_string          )                :: massType

    !![
    <inputParameter>
      <name>alpha</name>
      <description>The parameter $\alpha$ of the Zhao1996 profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <description>The parameter $\beta$ of the Zhao1996 profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <description>The parameter $\gamma$ of the Zhao1996 profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>(3.0d0-gamma)/4.0d0/Pi/Hypergeometric_2F1([(3.0d0-gamma)/alpha,(beta-gamma)/alpha],[1.0d0+(3.0d0-gamma)/alpha],-1.0d0)</defaultValue>
      <description>The density normalization of the Zhao1996 profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scaleLength</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scale radius of the Zhao1996 profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the Zhao1996 profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusOuter</name>
      <description>The outer radius of the Zhao1996 profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the Zhao1996 profile is considered to be dimensionless.</description>
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
     <call>self=massDistributionZhao1996(alpha=alpha,beta=beta,gamma=gamma,scaleLength=scaleLength,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="densityNormalization" value="densityNormalization" parameterPresent="parameters"/>
     <argument name="mass"                 value="mass"                 parameterPresent="parameters"/>
     <argument name="radiusOuter"          value="radiusOuter"          parameterPresent="parameters"/>
     <argument name="dimensionless"        value="dimensionless"        parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massDistributionZhao1996ConstructorParameters

  function massDistributionZhao1996ConstructorInternal(alpha,beta,gamma,scaleLength,densityNormalization,mass,radiusOuter,dimensionless,componentType,massType) result(self)
    !!{
    Internal constructor for ``zhao1996'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Comparison    , only : Values_Agree
    implicit none
    type            (massDistributionZhao1996     )                         :: self
    double precision                              , intent(in   )           :: alpha               , beta       , &
         &                                                                     gamma               , scaleLength
    double precision                              , intent(in   ), optional :: densityNormalization, radiusOuter, &
         &                                                                     mass
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree
    !![
    <constructorAssign variables="alpha, beta, gamma, scaleLength, componentType, massType"/>
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
       self%densityNormalization=+mass/self%scaleLength**3*(3.0d0-gamma)/4.0d0/Pi/radiusScaleFree**(3.0d0-gamma)/Hypergeometric_2F1([(3.0d0-gamma)/alpha,(beta-gamma)/alpha],[1.0d0+(3.0d0-gamma)/alpha],-radiusScaleFree**alpha)
    else
       call Error_Report('either "densityNormalization", or "mass" and "radiusOuter" must be specified'//{introspection:location})
    end if
    ! Determine if profile is dimensionless.
    if      (present(dimensionless     )) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! Validate parameters.
    if (gamma >= 3.0d0) call Error_Report('γ ≥ 3 gives divergent mass as r → 0'//{introspection:location})
    ! Detect special cases.
    if      (                                         &
         &    Values_Agree(alpha,1.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(beta ,3.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(gamma,1.0d0,absTol=1.0d-6) &
         &  ) then
       ! The "NFW" profile.
       self%specialCase=specialCaseNFW
    else if (                                         &
         &    Values_Agree(alpha,1.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(beta ,3.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(gamma,0.0d0,absTol=1.0d-6) &
         &  ) then
       ! The "cored NFW" profile.
       self%specialCase=specialCaseCoredNFW
    else if (                                         &
         &    Values_Agree(alpha,1.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(beta ,3.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(gamma,0.5d0,absTol=1.0d-6) &
         &  ) then
       ! The "γ=1/2 NFW" profile.
       self%specialCase=specialCaseGamma0_5NFW
    else if (                                         &
         &    Values_Agree(alpha,1.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(beta ,3.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(gamma,1.5d0,absTol=1.0d-6) &
         &  ) then
       ! The "γ=3/2 NFW" profile.
       self%specialCase=specialCaseGamma1_5NFW
    else
       ! Use general solutions.
       self%specialCase=specialCaseGeneral
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
  end function massDistributionZhao1996ConstructorInternal

  double precision function zhao1996MassTotal(self) result(massTotal)
    !!{
    Return the total mass in an Zhao1996 mass distribution.
    !!}
    use :: Gamma_Functions, only : Gamma_Function
    implicit none
    class(massDistributionZhao1996), intent(inout) :: self

    if (self%beta <= 3.0d0) then
       massTotal=+huge(0.0d0)
    else
       massTotal=+ 4.0d0                                                              &
            &    /(3.0d0-self%gamma)                                                  &
            &    *Gamma_Function((-3.0d0           +self%beta           )/self%alpha) &
            &    *Gamma_Function((+3.0d0+self%alpha          -self%gamma)/self%alpha) &
            &    /Gamma_Function((      -self%alpha+self%beta           )/self%alpha) &
            &    *self%densityNormalization                                           &
            &    *self%scaleLength         **3
    end if
    return
  end function zhao1996MassTotal

  double precision function zhao1996Density(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a Zhao1996 mass distribution.
    !!}
    implicit none
    class           (massDistributionZhao1996), intent(inout) :: self
    class           (coordinate              ), intent(in   ) :: coordinates
    double precision                                          :: radiusScaleFree

    ! Compute the density at this position.
    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %scaleLength
    select case (self%specialCase%ID)
    case (specialCaseGeneral%ID)
       density=+self%densityNormalization              &
            &  /  radiusScaleFree**self%gamma          &
            &  /(                                      &
            &    +1.0d0                                &
            &    +radiusScaleFree**self%alpha          &
            &   )**((self%beta-self%gamma)/self%alpha)
    case (specialCaseNFW%ID)
       density=+self%densityNormalization              &
            &   /  radiusScaleFree                     &
            &   /(                                     &
            &     +1.0d0                               &
            &     +radiusScaleFree                     &
            &    )**2
    case (specialCaseCoredNFW%ID)
       density=+self%densityNormalization              &
            &  /(                                      &
            &    +1.0d0                                &
            &    +radiusScaleFree                      &
            &   )**3
    case (specialCaseGamma0_5NFW%ID)
       density=+self%densityNormalization              &
            &  /sqrt(radiusScaleFree)                  &
            &  /(                                      &
            &    +1.0d0                                &
            &    +radiusScaleFree                      &
            &   )**2.5d0
    case (specialCaseGamma1_5NFW%ID)
       density=+self%densityNormalization              &
            &  /  radiusScaleFree**1.5d0               &
            &  /(                                      &
            &    +1.0d0                                &
            &    +radiusScaleFree                      &
            &   )**1.5d0
    case default
       density=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    return
  end function zhao1996Density

  double precision function zhao1996DensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an Zhao1996 \citep{zhao_analytical_1996} mass distribution.
    !!}
    implicit none
    class           (massDistributionZhao1996), intent(inout), target   :: self
    class           (coordinate              ), intent(in   )           :: coordinates
    logical                                   , intent(in   ), optional :: logarithmic
    double precision                                                    :: radiusScaleFree
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    radiusScaleFree      =+coordinates%rSpherical()                  &
         &                /self       %scaleLength
    densityGradientRadial=-(                                         &
         &                  +self%beta  *radiusScaleFree**self%alpha &
         &                  +self%gamma                              &
         &                 )                                         &
         &                /(                                         &
         &                  +1.0d0                                   &
         &                  +            radiusScaleFree**self%alpha &
         &                 )
    if (.not.logarithmic_) densityGradientRadial=+            densityGradientRadial              &
         &                                       *self       %density              (coordinates) &
         &                                       /coordinates%rSpherical           (           )
    return
  end function zhao1996DensityGradientRadial

  double precision function zhao1996DensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Computes radial moments of the density in an Zhao1996 \citep{zhao_analytical_1996} mass distribution.
    !!}
    implicit none
    class           (massDistributionZhao1996), intent(inout)           :: self
    double precision                          , intent(in   )           :: moment
    double precision                          , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    logical                                   , intent(  out), optional :: isInfinite
    double precision                                                    :: radialMomentMinimum, radialMomentMaximum, &
         &                                                                 radiusScaleFree

    densityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    if (present(radiusMinimum)) then
       radiusScaleFree     =+     radiusMinimum &
            &               /self%scaleLength
       radialMomentMinimum=+radialMomentIndefinite(radiusScaleFree)
    else
       radialMomentMinimum=+0.0d0
       if     (                            &
            &   0.0d0        >= self%alpha &
            &  .or.                        &
            &   1.0d0+moment <= self%gamma &
            & ) call Error_Report('radial moment is undefined'//{introspection:location})
    end if
    if (present(radiusMaximum)) then
       radiusScaleFree    =+     radiusMaximum &
            &              /self%scaleLength
       radialMomentMaximum=+radialMomentIndefinite(radiusScaleFree)
    else
       radialMomentMaximum=+0.0d0
       if     (                                          &
            &         self%alpha  > 0.0d0                &
            &  .and.                                     &
            &   1.0d0+     moment > self%alpha+self%beta &
            &  .and.                                     &
            &   1.0d0+     moment > self%gamma           &
            &  .and.                                     &
            &         self%beta   > self%gamma           &
            & ) call Error_Report('radial moment is undefined'//{introspection:location})
    end if
    densityRadialMoment=+(                                         &
         &                +radialMomentMaximum                     &
         &                -radialMomentMinimum                     &
         &               )                                         &
         &              *self%densityNormalization                 &
         &              *self%scaleLength         **(moment+1.0d0)
    return

  contains

    double precision function radialMomentIndefinite(radiusScaleFree)
      !!{
      Compute the indefinite radial moment.
      !!}
      use :: Hypergeometric_Functions, only : Hypergeometric_2F1
      implicit none
      double precision, intent(in   ) :: radiusScaleFree

      radialMomentIndefinite=+radiusScaleFree**   (1.0d0+moment-self%gamma)                                                                                                                         &
           &                 /                    (1.0d0+moment-self%gamma)                                                                                                                         &
           &                 *Hypergeometric_2F1([(1.0d0+moment-self%gamma)/self%alpha,(self%beta-self%gamma)/self%alpha],[1.0d0+(1.0d0+moment-self%gamma)/self%alpha],-radiusScaleFree**self%alpha)
      return
    end function radialMomentIndefinite

  end function zhao1996DensityRadialMoment

  double precision function zhao1996MassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for zhao1996 mass distributions.
    !!}
    implicit none
    class           (massDistributionZhao1996), intent(inout), target :: self
    double precision                          , intent(in   )         :: radius
    double precision                                                  :: radiusScaleFree
    
    self_           =>  self
    radiusScaleFree =  +      radius                              &
         &              /self%scaleLength
    mass            =  +self%densityNormalization                 &
         &             *self%scaleLength                      **3 &
         &             *massEnclosedScaleFree(radiusScaleFree)
    return
  end function zhao1996MassEnclosedBySphere

  double precision function zhao1996VelocityRotationCurveMaximum(self) result(velocity)
    !!{
    Return the peak velocity in the rotation curve for a Zhao1996 mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class(massDistributionZhao1996), intent(inout) :: self

    select case (self%specialCase%ID)
    case (specialCaseGeneral    %ID)
       velocity=+self%rotationCurve(self%radiusRotationCurveMaximum())
       return
    case (specialCaseNFW        %ID)
       velocity=+1.6483500453640064578d0
    case (specialCaseGamma0_5NFW%ID)
       velocity=+1.4026527358517898624d0
    case (specialCaseGamma1_5NFW%ID)
       velocity=+2.0932014912026087087d0
    case (specialCaseCoredNFW   %ID)
       velocity=+1.2414383571567440046d0
    case default
       velocity=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    velocity=+velocity                                    &
         &   *sqrt(                                       &
         &         +self%densityNormalization             &
         &        )                                       &
         &   *      self%scaleLength
    if (.not.self%isDimensionless())                      &
         & velocity=+velocity                             &
         &          *sqrt(gravitationalConstant_internal)
    return
  end function zhao1996VelocityRotationCurveMaximum

  double precision function zhao1996RadiusRotationCurveMaximum(self) result(radius)
    !!{
    Return the peak velocity in the rotation curve for a Zhao1996 mass distribution.
    !!}
    implicit none
    class(massDistributionZhao1996), intent(inout), target :: self
    
    select case (self%specialCase%ID)
    case (specialCaseGeneral    %ID)
       radius=+self%radiusRotationCurveMaximumNumerical()
       return
    case (specialCaseNFW        %ID)
       radius=+2.1625815870646098349d0
    case (specialCaseGamma0_5NFW%ID)
       radius=+3.2892765613841120232d0
    case (specialCaseGamma1_5NFW%ID)
       radius=+1.0549665718691230692d0
    case (specialCaseCoredNFW   %ID)
       radius=+4.4247006468722702269d0
    case default
       radius=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    radius=+radius           &
         & *self%scaleLength
    return
  end function zhao1996RadiusRotationCurveMaximum
  
  double precision function zhao1996RadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for zhao1996 mass distributions.
    !!}    
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    use :: Error           , only : Error_Report
    implicit none
    class           (massDistributionZhao1996), intent(inout), target       :: self
    double precision                          , intent(in   ), optional     :: mass                       , massFractional
    double precision                          , allocatable  , dimension(:) :: radii                      , masses
    double precision                          , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                        :: massScaleFree              , mass_
    integer                                                                 :: countRadii

    mass_=0.0d0
    if (present(mass)) then
       mass_=mass
    else if (present(massFractional)) then
       call Error_Report('mass is unbounded, so mass fraction is undefined'//{introspection:location})
    else
       call Error_Report('either mass or massFractional must be supplied'  //{introspection:location})
    end if
    massScaleFree=+     mass_                   &
         &        /self%densityNormalization    &
         &        /self%scaleLength         **3
    if     (                                            &
         &   massScaleFree <= self%massScaleFreeMinimum &
         &  .or.                                        &
         &   massScaleFree >  self%massScaleFreeMaximum &
         & ) then
       self_ => self
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
       radii                     =  Make_Range(self%massScaleFreeRadiusMinimum,self%massScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
       masses                    =  massEnclosedScaleFree(           radii)
       self%massScaleFreeMinimum =  masses               (         1      )
       self%massScaleFreeMaximum =  masses               (countRadii      )
       self%massScaleFree_       =  interpolator         (masses    ,radii)
    end if
    radius=+self%massScaleFree_%interpolate(massScaleFree) &
         & *self%scaleLength
    return
  end function zhao1996RadiusEnclosingMass
  
  double precision function zhao1996RadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mean density for zhao1996 mass distributions.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionZhao1996), intent(inout), target       :: self
    double precision                          , intent(in   )               :: density
    double precision                          , intent(in   ), optional     :: radiusGuess
    double precision                          , allocatable  , dimension(:) :: radii                      , densities
    double precision                          , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                        :: densityScaleFree
    integer                                                                 :: countRadii

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
       self_                        =>  self
       radii                        =   Make_Range(self%densityScaleFreeRadiusMinimum,self%densityScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
       densities                    =  -densityEnclosedScaleFree(           radii)
       self%densityScaleFreeMinimum =  -densities               (countRadii      )
       self%densityScaleFreeMaximum =  -densities               (         1      )
       self%densityScaleFree_       =   interpolator            (densities ,radii)
    end if
    radius=+self%densityScaleFree_%interpolate(-densityScaleFree) &
         & *self%scaleLength
    return    
  end function zhao1996RadiusEnclosingDensity

  impure elemental double precision function massEnclosedScaleFree(radius) result(mass)
    !!{
    Evaluate the mass enclosed by a given radius in a scale-free Zhao1996 mass distribution.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    implicit none
    double precision, intent(in   ) :: radius
    double precision, parameter     :: radiusTiny=1.0d-3
    
    select case (self_%specialCase%ID)
    case (specialCaseGeneral%ID)
       mass   =+4.0d0                                                                                                                                                   &
            &  *Pi                                                                                                                                                      &
            &  *radius          **  (3.0d0-self_%gamma)                                                                                                                 &
            &  *Hypergeometric_2F1([(3.0d0-self_%gamma)/self_%alpha,(self_%beta-self_%gamma)/self_%alpha],[1.0d0+(3.0d0-self_%gamma)/self_%alpha],-radius**self_%alpha) &
            &  /                    (3.0d0-self_%gamma)
    case (specialCaseNFW%ID)
       if (radius < radiusTiny) then
          ! Use series solution for small radii.
          mass   =+ 2.0d0      *Pi*radius**2 &
               &  - 8.0d0/3.0d0*Pi*radius**3 &
               &  + 3.0d0      *Pi*radius**4 &
               &  -16.0d0/5.0d0*Pi*radius**5 &
               &  +10.0d0/3.0d0*Pi*radius**6 &
               &  -24.0d0/7.0d0*Pi*radius**7
       else
          ! Use full solution.
          mass   =+4.0d0                &
               &  *Pi                   &
               &  *(                    &
               &    +log(+1.0d0+radius) &
               &    -           radius  &
               &    /   (+1.0d0+radius) &
               &   )
       end if
    case (specialCaseGamma0_5NFW%ID)
       if (radius <radiusTiny) then
          ! Use series solution for small radii.
          mass   =+  8.0d0/ 5.0d0*Pi*radius**2.5d0 &
               &  - 20.0d0/ 7.0d0*Pi*radius**3.5d0 &
               &  + 35.0d0/ 9.0d0*Pi*radius**4.5d0 &
               &  -105.0d0/22.0d0*Pi*radius**5.5d0
       else
          ! Use full solution.
          mass   =+(                             &
               &    -8.0d0                       &
               &    *Pi                          &
               &    *sqrt(radius)                &
               &    *(3.0d0+4.0d0*radius)        &
               &   )                             &
               &  /(                             &
               &    +3.0d0                       &
               &    *(1.0d0+      radius)**1.5d0 &
               &    )                            &
               &  +8.0d0                         &
               &  *Pi                            &
               &  *asinh(sqrt(radius))
       end if
    case (specialCaseGamma1_5NFW%ID)
       if (radius <radiusTiny) then
          ! Use series solution for small radii.
          mass   =+  8.0d0/  3.0d0*Pi*radius**1.5d0 &
               &  - 12.0d0/  5.0d0*Pi*radius**2.5d0 &
               &  + 15.0d0/  7.0d0*Pi*radius**3.5d0 &
               &  - 35.0d0/ 18.0d0*Pi*radius**4.5d0 &
               &  +315.0d0/176.0d0*Pi*radius**5.5d0
       else
          ! Use full solution.
          mass   =+8.0d0                               &
               &  *Pi                                  &
               &  *(                                   &
               &    -      sqrt(radius/(1.0d0+radius)) &
               &    +asinh(sqrt(radius              )) &
               &   )
       end if
    case (specialCaseCoredNFW%ID)
       if (radius <radiusTiny) then
          ! Use series solution for small radii.
          mass   =+ 4.0d0/3.0d0*Pi*radius**3 &
               &  - 3.0d0      *Pi*radius**4 &
               &  +24.0d0/5.0d0*Pi*radius**5 &
               &  -20.0d0/3.0d0*Pi*radius**6
       else
          ! Use full solution.
          mass   =+4.0d0                         &
               &  *Pi                            &
               &  *(                             &
               &    +log(+1.0d0+      radius)    &
               &    -                 radius     &
               &    *   (+2.0d0+3.0d0*radius)    &
               &    /     2.0d0                  &
               &    /   (+1.0d0+      radius)**2 &
               &   )
       end if
    case default
       mass   =+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    return
  end function massEnclosedScaleFree

  impure elemental double precision function densityEnclosedScaleFree(radius) result(density)
    !!{
    Evaluate the mean enclosed density at a given radius in a scale-free Zhao1996 mass distribution.
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
  
  double precision function zhao1996RadiusFromSpecificAngularMomentum(self,angularMomentumSpecific) result(radius)
    !!{
    Computes the radius corresponding to a given specific angular momentum for zhao1996 mass distributions.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Ranges                , only : Make_Range                    , rangeTypeLogarithmic
    implicit none
    class           (massDistributionZhao1996), intent(inout), target       :: self
    double precision                          , intent(in   )               :: angularMomentumSpecific
    double precision                          , allocatable  , dimension(:) :: radii                                   , angularMomentaSpecific
    double precision                          , parameter                   :: countRadiiPerDecade             =100.0d0
    double precision                                                        :: angularMomentumSpecificScaleFree
    integer                                                                 :: countRadii

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
          self_                                        => self
          radii                                        =  Make_Range(self%angularMomentumSpecificScaleFreeRadiusMinimum,self%angularMomentumSpecificScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
          angularMomentaSpecific                       =  angularMomentumSpecificEnclosedScaleFree(                       radii)
          self%angularMomentumSpecificScaleFreeMinimum =  angularMomentaSpecific                  (                     1      )
          self%angularMomentumSpecificScaleFreeMaximum =  angularMomentaSpecific                  (            countRadii      )
          self%angularMomentumSpecificScaleFree_       =  interpolator                            (angularMomentaSpecific,radii)
       end if
       radius=+self%angularMomentumSpecificScaleFree_%interpolate(angularMomentumSpecificScaleFree) &
            & *self%scaleLength
    else
       radius=+0.0d0
    end if
    return    
  end function zhao1996RadiusFromSpecificAngularMomentum

  impure elemental double precision function angularMomentumSpecificEnclosedScaleFree(radius) result(angularMomentumSpecific)
    !!{
    Evaluate the specific angular momentum at a given radius in a scale-free Zhao1996 mass distribution.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    angularMomentumSpecific=+sqrt(                               &
         &                        +massEnclosedScaleFree(radius) &
         &                        *                      radius  &
         &                       )
    return
  end function angularMomentumSpecificEnclosedScaleFree

  logical function zhao1996PotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionZhao1996), intent(inout) :: self

    isAnalytic=.true.
    return
  end function zhao1996PotentialIsAnalytic

  double precision function zhao1996Potential(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in an zhao1996 mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess     , structureErrorCodeInfinite
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Error                           , only : Error_Report
    implicit none
    class           (massDistributionZhao1996         ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                                             :: radiusScaleFree
    
    if (present(status)) status=structureErrorCodeSuccess
    self_           =>  self
    radiusScaleFree =  +coordinates%rSpherical () &
         &             /self       %scaleLength
    potential=+potentialScaleFree       (radiusScaleFree)    &
         &    *self%densityNormalization                     &
         &    *self%scaleLength                          **2
    if (.not.self%isDimensionless()) potential=+gravitationalConstant_internal &
         &                                     *potential
    return
  end function zhao1996Potential

  impure elemental double precision function potentialScaleFree(radius) result(potential)
    !!{
    Compute the potential in a scale-free Zhao1996 mass distribution.
    !!}
    use :: Gamma_Functions         , only : Gamma_Function
    use :: Numerical_Constants_Math, only : Pi
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1_Regularized
    implicit none
    double precision, intent(in   ) :: radius

    select case (self_%specialCase%ID)
    case (specialCaseGeneral    %ID)
       potential=+4.0d0                                                                                                                                                                                                                                   &
            &    *Pi                                                                                                                                                                                                                                      &
            &    *radius         **(2.0d0            -self_%gamma)                                                                                                                                                                                        &
            &    /                 (3.0d0            -self_%gamma)                                                                                                                                                                                        &
            &    *  Gamma_Function((3.0d0+self_%alpha-self_%gamma)/self_%alpha)                                                                                                                                                                           &
            &    /  Gamma_Function((3.0d0            -self_%gamma)/self_%alpha)                                                                                                                                                                           &
            &    *(                                                                                                                                                                                                                                       &
            &      +Gamma_Function((2.0d0            -self_%gamma)/self_%alpha)*Hypergeometric_2F1_Regularized([(2.0d0-self_%gamma)/self_%alpha,(self_%beta-self_%gamma)/self_%alpha],[(2.0d0+self_%alpha-self_%gamma)/self_%alpha],-radius**self_%alpha) &
            &      -Gamma_Function((3.0d0            -self_%gamma)/self_%alpha)*Hypergeometric_2F1_Regularized([(3.0d0-self_%gamma)/self_%alpha,(self_%beta-self_%gamma)/self_%alpha],[(3.0d0+self_%alpha-self_%gamma)/self_%alpha],-radius**self_%alpha) &
            &     )
    case (specialCaseNFW        %ID)
       potential=+4.0d0                &
            &    *Pi                   &
            &    *(                    &
            &      +     1.0d0         &
            &      -log(+1.0d0+radius) &
            &      /           radius  &
            &     )
    case (specialCaseGamma0_5NFW%ID)
       potential=+8.0d0                                &
            &    *Pi                                   &
            &    *(                                    &
            &      +                  (+3.0d0+radius)  &
            &      /                    3.0d0          &
            &      /      sqrt(radius*(+1.0d0+radius)) &
            &      -asinh(sqrt(radius               )) &
            &      /           radius                  &
            &     )
    case (specialCaseGamma1_5NFW%ID)
       potential=+8.0d0                                &
            &    *Pi                                   &
            &    *(                                    &
            &      +      sqrt(radius*(1.0d0+radius))  &
            &      -asinh(sqrt(radius               )) &
            &     )                                    &
            &    /radius
    case (specialCaseCoredNFW   %ID)
       potential=+2.0d0                &
            &    *Pi                   &
            &    *(                    &
            &      +   (+2.0d0+radius) &
            &      /   (+1.0d0+radius) &
            &      -     2.0d0         &
            &      *log(+1.0d0+radius) &
            &      /           radius  &
            &   )
    case default
       potential=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    return
  end function potentialScaleFree

 double precision function potentialDifferenceScaleFree(radius1,radius2) result(potential)
    !!{
    Compute the potential difference in a scale-free Zhao1996 mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Comparison    , only : Values_Agree
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    implicit none
    double precision, intent(in   ) :: radius1                            , radius2
    double precision, parameter     :: radiusSmall                 =1.0d-3
    double precision, parameter     :: toleranceRelative           =1.0d-3
    double precision                :: potentialGradientLogarithmic       , radiusDifferenceLogarithmic
    
    if (Values_Agree(radius1,radius2,relTol=toleranceRelative) .or. max(radius1,radius2) < radiusSmall) then
       if (radius1 < radiusSmall) then
          select case (self_%specialCase%ID)
          case (specialCaseGeneral    %ID)
             potentialGradientLogarithmic=-  1.0d0                                                                                                                                                          &
                  &                       /(                                                                                                                                                                &
                  &                         -1.0d0                                                                                                                                                          &
                  &                         +                    (3.0d0-self_%gamma)                                                                                                                        &
                  &                         /                    (2.0d0-self_%gamma)                                                                                                                        &
                  &                         *Hypergeometric_2F1([(2.0d0-self_%gamma)/self_%alpha,(self_%beta-self_%gamma)/self_%alpha],[(2.0d0+self_%alpha-self_%gamma)/self_%alpha],-radius1**self_%alpha) &
                  &                         /Hypergeometric_2F1([(3.0d0-self_%gamma)/self_%alpha,(self_%beta-self_%gamma)/self_%alpha],[(3.0d0+self_%alpha-self_%gamma)/self_%alpha],-radius1**self_%alpha) &
                  &                        )
          case (specialCaseNFW        %ID)
             potentialGradientLogarithmic=+1.0d0- 2.0d0*radius1/ 3.0d0+  5.0d0*radius1**2/  9.0d0-  67.0d0*radius1**3/ 135.0d0+    371.0d0*radius1**4/    810.0d0
          case (specialCaseGamma0_5NFW%ID)
             potentialGradientLogarithmic=+1.5d0-15.0d0*radius1/14.0d0+275.0d0*radius1**2/294.0d0-6525.0d0*radius1**3/7546.0d0+5067175.0d0*radius1**4/6180174.0d0
          case (specialCaseGamma1_5NFW%ID)
             potentialGradientLogarithmic=+0.5d0- 3.0d0*radius1/10.0d0+ 81.0d0*radius1**2/350.0d0- 341.0d0*radius1**3/1750.0d0+ 115477.0d0*radius1**4/ 673750.0d0
          case (specialCaseCoredNFW   %ID)
             potentialGradientLogarithmic=+2.0d0- 3.0d0*radius1/ 2.0d0+ 27.0d0*radius1**2/ 20.0d0-  51.0d0*radius1**3/  40.0d0+   3441.0d0*radius1**4/   2800.0d0
          case default
             potentialGradientLogarithmic=+0.0d0
             call Error_Report('unknown special case'//{introspection:location})
          end select
       else
          select case (self_%specialCase%ID)
          case (specialCaseGeneral    %ID)
             potentialGradientLogarithmic=+  1.0d0                                                                                                                                                          &
                  &                       /(                                                                                                                                                                &
                  &                         -1.0d0                                                                                                                                                          &
                  &                         +                    (3.0d0-self_%gamma)                                                                                                                        &
                  &                         /                    (2.0d0-self_%gamma)                                                                                                                        &
                  &                         *Hypergeometric_2F1([(2.0d0-self_%gamma)/self_%alpha,(self_%beta-self_%gamma)/self_%alpha],[(2.0d0+self_%alpha-self_%gamma)/self_%alpha],-radius1**self_%alpha) &
                  &                         /Hypergeometric_2F1([(3.0d0-self_%gamma)/self_%alpha,(self_%beta-self_%gamma)/self_%alpha],[(3.0d0+self_%alpha-self_%gamma)/self_%alpha],-radius1**self_%alpha) &
                  &                        )
          case (specialCaseNFW        %ID)
             potentialGradientLogarithmic=+  1.0d0                   &
                  &                       /(                         &
                  &                         -1.0d0                   &
                  &                         +            radius1 **2 &
                  &                         /(                       &
                  &                           -          radius1     &
                  &                           +   (1.0d0+radius1)    &
                  &                           *log(1.0d0+radius1)    &
                  &                          )                       &
                  &                        )
          case (specialCaseGamma0_5NFW%ID)
             potentialGradientLogarithmic=+(                                                                           &
                  &                         -        (3.0d0+4.0d0*radius1)   *      sqrt(radius1) *sqrt(1.0d0+radius1) &
                  &                         +  3.0d0*(1.0d0+      radius1)**2*asinh(sqrt(radius1))                     &
                  &                        )                                                                           &
                  &                       /(                                                                           &
                  &                         +        (1.0d0+      radius1)                                             &
                  &                         *(                                                                         &
                  &                           +      (3.0d0+      radius1)   *      sqrt(radius1) *sqrt(1.0d0+radius1) &
                  &                           -3.0d0*(1.0d0+      radius1)   *asinh(sqrt(radius1))                     &
                  &                          )                                                                         &
                  &                        )
          case (specialCaseGamma1_5NFW%ID)
             potentialGradientLogarithmic=-1.0d0                                 &
                  &                       /(                                     &
                  &                         +1.0d0                               &
                  &                         +1.0d0                               &
                  &                         /(                                   &
                  &                           +1.0d0                             &
                  &                           /                 radius1          &
                  &                           -      sqrt(1.0d0+radius1       )  &
                  &                           *asinh(sqrt(      radius1       )) &
                  &                           /                 radius1**1.5d0   &
                  &                          )                                   &
                  &                        )
          case (specialCaseCoredNFW   %ID)
             potentialGradientLogarithmic=+1.0d0                                   &
                  &                       /(                                       &
                  &                         -  1.0d0                               &
                  &                         +  radius1**3                          &
                  &                         /(                                     &
                  &                           -radius1*   (2.0d0+3.0d0*radius1)    &
                  &                           +2.0d0  *   (1.0d0+      radius1)**2 &
                  &                           *        log(1.0d0+      radius1)    &
                  &                          )                                     &
                  &                        )
          case default
             potentialGradientLogarithmic=+0.0d0
             call Error_Report('unknown special case'//{introspection:location})
          end select
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
    
  double precision function zhao1996RadiusFreefall(self,time) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in an Zhao1996 mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr, gravitationalConstant_internal
    implicit none
    class           (massDistributionZhao1996), intent(inout) :: self
    double precision                          , intent(in   ) :: time
    double precision                                          :: timeScaleFree, timeScale
    
    timeScale    =+1.0d0/sqrt(                                &
         &                    +gravitationalConstant_internal &
         &                    *self%densityNormalization      &
         &                   )                                &
         &        *MpcPerKmPerSToGyr
    timeScaleFree=+time                                       &
         &        /timeScale
    if (self%specialCase == specialCaseCoredNFW .and. timeScaleFree <= timeFreefallScaleFreeMinimumCoredNFW) then
       radius=0.0d0
       return
    end if
    call self%timeFreefallTabulate(timeScaleFree)
    radius=+self%timeFreefallScaleFree_%interpolate(timeScaleFree) &
         & *self%scaleLength
    return   
  end function zhao1996RadiusFreefall
  
  double precision function zhao1996RadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in an zhao1996 mass
    distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr, gravitationalConstant_internal
    implicit none
    class           (massDistributionZhao1996), intent(inout) :: self
    double precision                          , intent(in   ) :: time
    double precision                                          :: timeScaleFree, timeScale

    timeScale    =+1.0d0/sqrt(                                &
         &                    +gravitationalConstant_internal &
         &                    *self%densityNormalization      &
         &                   )                                &
         &        *MpcPerKmPerSToGyr
    timeScaleFree=+time                                       &
         &        /timeScale
    if (self%specialCase == specialCaseCoredNFW .and. timeScaleFree <= timeFreefallScaleFreeMinimumCoredNFW) then
       radiusIncreaseRate=0.0d0
       return
    end if
    call self%timeFreefallTabulate(timeScaleFree)
    radiusIncreaseRate=+self%timeFreefallScaleFree_%derivative(timeScaleFree) &
         &             *self%scaleLength                                      &
         &             /     timeScale
    return
  end function zhao1996RadiusFreefallIncreaseRate
  
  subroutine zhao1996TimeFreefallTabulate(self,timeScaleFree)
    !!{
    Tabulate the freefall radius at the given {\normalfont \ttfamily time} in an Zhao1996 mass distribution.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges     , only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionZhao1996), intent(inout), target       :: self
    double precision                          , intent(in   )               :: timeScaleFree
    double precision                          , allocatable  , dimension(:) :: radii                      , timesFreefall
    double precision                          , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                        :: radiusStart
    integer                                                                 :: countRadii                 , i
    type            (integrator              )                              :: integrator_

    if     (                                                    &
         &   timeScaleFree <= self%timeFreefallScaleFreeMinimum &
         &  .or.                                                &
         &   timeScaleFree >  self%timeFreefallScaleFreeMaximum &
         & ) then
       self_       => self
       integrator_ =  integrator(timeFreeFallIntegrand,toleranceRelative=1.0d-6)
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
      Evaluate the freefall time from a given radius in a scale-free Zhao1996 mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      radiusStart          =                            radius
      timeFreefallScaleFree=integrator_%integrate(0.0d0,radius)
      return
    end function timeFreefallScaleFree
    
    double precision function timeFreeFallIntegrand(radius)
      !!{
      Integrand used to find the freefall time in a scale-free Zhao1996 mass distribution.
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
    
  end subroutine zhao1996TimeFreefallTabulate

  double precision function zhao1996FourierTransform(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in an Zhao1996 mass
    distribution.
    !!}
    use :: Exponential_Integrals   , only : Exponential_Integral
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionZhao1996), intent(inout) :: self
    double precision                          , intent(in   ) :: radiusOuter        , wavenumber
    double precision                                          :: wavenumberScaleFree, radiusOuterScaleFree

    waveNumberScaleFree =+waveNumber *self%scaleLength
    radiusOuterScaleFree=+radiusOuter/self%scaleLength
    select case (self%specialCase%ID)
    case   (                           &
         &  specialCaseGeneral    %ID, &
         &  specialCaseGamma0_5NFW%ID, &
         &  specialCaseGamma1_5NFW%ID  &
         & )
       fourierTransform=+self%fourierTransformNumerical(radiusOuter,wavenumber)
       return
    case (specialCaseNFW        %ID)
       fourierTransform=+dimag(                                                                                                  &
            &                  +4.0d0                                                                                            &
            &                  *Pi                                                                                               &
            &                  *(                                                                                                &
            &                    -exp(dcmplx(0.0d0,1.0d0)*wavenumberScaleFree*(1.0d0+radiusOuterScaleFree))                      &
            &                    +                                            (1.0d0+radiusOuterScaleFree)                       &
            &                    *(                                                                                              &
            &                      +exp(dcmplx(0.0d0,1.0d0)*wavenumberScaleFree)                                                 &
            &                      -    dcmplx(0.0d0,1.0d0)*wavenumberScaleFree                                                  &
            &                      *(                                                                                            &
            &                        +Exponential_Integral(dcmplx(0.0d0,1.0d0)*wavenumberScaleFree                             ) &
            &                        -Exponential_Integral(dcmplx(0.0d0,1.0d0)*wavenumberScaleFree*(1.0d0+radiusOuterScaleFree)) &
            &                       )                                                                                            &
            &                     )                                                                                              &
            &                   )                                                                                                &
            &                  /exp(dcmplx(0.0d0,1.0d0)*wavenumberScaleFree)                                                     &
            &                  /wavenumberScaleFree                                                                              &
            &                  /(1.0d0+radiusOuterScaleFree)                                                                     &
            &                 )
    case (specialCaseCoredNFW   %ID)
       fourierTransform=+dimag(                                                                                                    &
            &                  +2.0d0                                                                                              &
            &                  *Pi                                                                                                 &
            &                  *(                                                                                                  &
            &                    +           1.0d0                                                                                 &
            &                    -    dcmplx(0.0d0,1.0d0)*                        wavenumberScaleFree                              &
            &                    +(                                                                                                &
            &                      +  dcmplx(0.0d0,1.0d0)*exp(dcmplx(0.0d0,1.0d0)*wavenumberScaleFree*radiusOuterScaleFree)        &
            &                      *(                                                                                              &
            &                        + dcmplx(0.0d0,1.0d0)+wavenumberScaleFree                                                     &
            &                        +(dcmplx(0.0d0,2.0d0)+wavenumberScaleFree)*radiusOuterScaleFree                               &
            &                       )                                                                                              &
            &                     )                                                                                                &
            &                    /(1.0d0+radiusOuterScaleFree)**2                                                                  &
            &                    -(                                                                                                &
            &                      +                     wavenumberScaleFree                                                       &
            &                      *(dcmplx(0.0d0,2.0d0)+wavenumberScaleFree)                                                      &
            &                      *(                                                                                              &
            &                        +Exponential_Integral(dcmplx(0.0d0,1.0d0)*wavenumberScaleFree)                                &
            &                        -Exponential_Integral(dcmplx(0.0d0,1.0d0)*wavenumberScaleFree*(1.0d0+radiusOuterScaleFree))   &
            &                       )                                                                                              &
            &                     )                                                                                                &
            &                    /exp(dcmplx(0.0d0,1.0d0)*wavenumberScaleFree)                                                     &
            &                   )                                                                                                  &
            &                  /wavenumberScaleFree                                                                                &
            &                 )
    case default
       fourierTransform=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    fourierTransform=+fourierTransform                            &
         &           /massEnclosedScaleFree(radiusOuterScaleFree)
    return
  end function zhao1996FourierTransform

  double precision function zhao1996EnergyPotential(self,radiusOuter) result(energy)
    !!{
    Compute the potential energy within a given {\normalfont \ttfamily radius} in a Zhao1996 mass distribution.
    \begin{eqnarray}
    \end{eqnarray}
    where $x=r/r_\mathrm{s}$ and $\mathrm{G}$ is Catalan's constant.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionZhao1996), intent(inout) :: self
    double precision                          , intent(in   ) :: radiusOuter
    double precision                                          :: radiusOuterScaleFree
    logical                                                   :: analytic
    
    analytic=.false.
    radiusOuterScaleFree=+     radiusOuter &
         &               /self%scaleLength
    select case (self_%specialCase%ID)
    case (specialCaseGeneral    %ID)
       analytic=.false.
    case (specialCaseNFW        %ID)
       analytic=.true.
       energy  =+8.0d0                                                                &
            &   *Pi**2                                                                &
            &   *(                                                                    &
            &     +             radiusOuterScaleFree *   (2.0d0+radiusOuterScaleFree) &
            &     -2.0d0*(1.0d0+radiusOuterScaleFree)*log(1.0d0+radiusOuterScaleFree) &
            &    )                                                                    &
            &   /(1.0d0+radiusOuterScaleFree)**2
    case (specialCaseGamma0_5NFW%ID)
       analytic=.true.
       energy  =-16.0d0                                                                        &
            &   *Pi**2                                                                         &
            &   *(                                                                             &
            &     -12.0d0                                                                      &
            &     *(1.0d0+radiusOuterScaleFree)                                                &
            &     *(                                                                           &
            &       -3.0d0*sqrt(radiusOuterScaleFree       )*sqrt(1.0d0+radiusOuterScaleFree)  &
            &       -4.0d0*     radiusOuterScaleFree**1.5d0 *sqrt(1.0d0+radiusOuterScaleFree)  &
            &       +3.0d0*sqrt(radiusOuterScaleFree        *    (1.0d0+radiusOuterScaleFree)) &
            &       +      sqrt(radiusOuterScaleFree**3     *    (1.0d0+radiusOuterScaleFree)) &
            &       +radiusOuterScaleFree                                                      &
            &       *(                                                                         &
            &         +3.0d0*sqrt(radiusOuterScaleFree   *(1.0d0+radiusOuterScaleFree))        &
            &         +      sqrt(radiusOuterScaleFree**3*(1.0d0+radiusOuterScaleFree))        &
            &        )                                                                         &
            &      )                                                                           &
            &     *asinh(sqrt(radiusOuterScaleFree))                                           &
            &     +radiusOuterScaleFree                                                        &
            &     *(                                                                           &
            &       +radiusOuterScaleFree                                                      &
            &       *(                                                                         &
            &         -6.0d0                                                                   &
            &         +radiusOuterScaleFree*(-3.0d0+5.0d0*radiusOuterScaleFree)                &
            &        )                                                                         &
            &       +6.0d0*(1.0d0+radiusOuterScaleFree)**3*log(1.0d0+radiusOuterScaleFree)     &
            &      )                                                                           &
            &    )                                                                             &
            &   /9.0d0                                                                         &
            &   /radiusOuterScaleFree                                                          &
            &   /(1.0d0+radiusOuterScaleFree)**3
    case (specialCaseGamma1_5NFW%ID)
       analytic=.true.
       energy  =+32.0d0                                                             &
            &   *Pi**2                                                              &
            &   *(                                                                  &
            &     -1.0d0                                                            &
            &     +1.0d0                                                            &
            &     /(1.0d0+radiusOuterScaleFree)                                     &
            &     +(                                                                &
            &       +2.0d0                                                          &
            &       *(                                                              &
            &         -radiusOuterScaleFree**1.5d0/sqrt(1.0d0+radiusOuterScaleFree) &
            &         +sqrt(radiusOuterScaleFree**3*(1.0d0+radiusOuterScaleFree))   &
            &        )                                                              &
            &       *asinh(sqrt(radiusOuterScaleFree))                              &
            &      )                                                                &
            &     /radiusOuterScaleFree**2                                          &
            &     +(                                                                &
            &       +(                                                              &
            &         +radiusOuterScaleFree**2                                      &
            &         +radiusOuterScaleFree**3                                      &
            &         -sqrt(radiusOuterScaleFree   *(1.0d0+radiusOuterScaleFree))   &
            &         *sqrt(radiusOuterScaleFree**3*(1.0d0+radiusOuterScaleFree))   &
            &        )                                                              &
            &       *asinh(sqrt(radiusOuterScaleFree))**2                           &
            &      )                                                                &
            &     /radiusOuterScaleFree**3                                          &
            &     /(1.0d0+radiusOuterScaleFree)                                     &
            &     -log(1.0d0+radiusOuterScaleFree)                                  &
            &    )
    case (specialCaseCoredNFW   %ID)
       analytic=.true.
       energy  =+2.0d0                                                                                                                  &
            &   *Pi**2                                                                                                                  &
            &   *(                                                                                                                      &
            &     +radiusOuterScaleFree*(12.0d0+radiusOuterScaleFree*(42.0d0+radiusOuterScaleFree*(40.0d0+7.0d0*radiusOuterScaleFree))) &
            &     -12.0d0*(1.0d0+radiusOuterScaleFree)**2*(1.0d0+2.0d0*radiusOuterScaleFree)*log(1.0d0+radiusOuterScaleFree)            &
            &   )                                                                                                                       &
            &  /(3.0d0*(1.0d0+radiusOuterScaleFree)**4)
    case default
       energy=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    if (analytic) then
       energy =-energy                            &
            &  *gravitationalConstant_internal    &
            &  *self%scaleLength              **5 &
            &  *self%densityNormalization     **2
    else
       energy =+self%energyPotentialNumerical(radiusOuter)
    end if
    return
  end function zhao1996EnergyPotential

  double precision function zhao1996EnergyKinetic(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the kinetic energy within a given {\normalfont \ttfamily radius} in a Zhao1996 mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    use :: Dilogarithms                    , only : Dilogarithm
    implicit none
    class           (massDistributionZhao1996), intent(inout) :: self
    double precision                          , intent(in   ) :: radiusOuter
    class           (massDistributionClass   ), intent(inout) :: massDistributionEmbedding
    logical                                                   :: analytic
    double precision                                          :: radiusOuterScaleFree
    
    analytic=.false.
    select type (massDistributionEmbedding)
    class is (massDistributionZhao1996)
       select type (kinematicsDistribution_ => massDistributionEmbedding%kinematicsDistribution_)
       class is (kinematicsDistributionZhao1996)
          radiusOuterScaleFree=+     radiusOuter &
               &               /self%scaleLength
          select case (self_%specialCase%ID)
          case (specialCaseGeneral    %ID)
             analytic=.false.
          case (specialCaseNFW        %ID)
             analytic=.true.
             energy  =+4.0d0                                                                                                                     &
                  &   *Pi**2                                                                                                                     &
                  &   *(                                                                                                                         &
                  &     +Pi**2                                                                                                                   &
                  &     *(+1.0d0+radiusOuterScaleFree   )                                                                                        &
                  &     *(+2.0d0+radiusOuterScaleFree**3)                                                                                        &
                  &     +radiusOuterScaleFree*(2.0d0-radiusOuterScaleFree*(2.0d0+7.0d0*radiusOuterScaleFree))                                    &
                  &     +2.0d0*radiusOuterScaleFree**4*atanh(1.0d0/(1.0d0+2.0d0*radiusOuterScaleFree))                                           &
                  &     -radiusOuterScaleFree**3*log(      radiusOuterScaleFree)                                                                 &
                  &     -2.0d0                  *log(1.0d0+radiusOuterScaleFree)                                                                 &
                  &     +(                                                                                                                       &
                  &       +              radiusOuterScaleFree                                                                                    &
                  &       - 3.0d0*       radiusOuterScaleFree**2                                                                                 &
                  &       - 5.0d0*       radiusOuterScaleFree**3                                                                                 &
                  &       +12.0d0*(1.0d0+radiusOuterScaleFree   )*log(radiusOuterScaleFree)                                                      &
                  &      )                                                                                                                       &
                  &     *log(1.0d0+radiusOuterScaleFree)                                                                                         &
                  &     + 3.0d0*(1.0d0+radiusOuterScaleFree)*(-2.0d0+radiusOuterScaleFree**3)*        log(       1.0d0+radiusOuterScaleFree )**2 &
                  &     + 6.0d0*(1.0d0+radiusOuterScaleFree)*(+2.0d0+radiusOuterScaleFree**3)*Dilogarithm(            -radiusOuterScaleFree )    &
                  &     -12.0d0*(1.0d0+radiusOuterScaleFree)                                 *Dilogarithm(1.0d0/(1.0d0+radiusOuterScaleFree))    &
                  &   )                                                                                                                          &
                  &   /(1.0d0+radiusOuterScaleFree)
          case (specialCaseGamma0_5NFW%ID)
             analytic=.true.
             energy  =+16.0d0                                                                                   &
                  &   /3.0d0                                                                                    &
                  &   *Pi**2                                                                                    &
                  &   *(                                                                                        &
                  &     -4.0d0*radiusOuterScaleFree**1.5d0                                                      &
                  &     *(-1.0d0+4.0d0*radiusOuterScaleFree+8.0d0*radiusOuterScaleFree**2)                      &
                  &     *asinh(sqrt(      radiusOuterScaleFree))                                                &
                  &     /      sqrt(1.0d0+radiusOuterScaleFree)                                                 &
                  &   +radiusOuterScaleFree                                                                     &
                  &   *(                                                                                        &
                  &     +  2.0d0+      radiusOuterScaleFree                                                     &
                  &     *(                                                                                      &
                  &       -5.0d0+2.0d0*radiusOuterScaleFree                                                     &
                  &       *(                                                                                    &
                  &         -11.0d0                                                                             &
                  &         +32.0d0*log(2.0d0)                                                                  &
                  &         + 8.0d0*radiusOuterScaleFree*(-1.0d0+radiusOuterScaleFree*log(16.0d0)+log(256.0d0)) &
                  &        )                                                                                    &
                  &      )                                                                                      &
                  &    )                                                                                        &
                  &   /  2.0d0                                                                                  &
                  &   /(+1.0d0+       radiusOuterScaleFree   )**2                                               &
                  &   +(-1.0d0+16.0d0*radiusOuterScaleFree**3)                                                  &
                  &   *log(1.0d0+radiusOuterScaleFree)                                                          &
                  &   )
          case (specialCaseGamma1_5NFW%ID)
             analytic=.true.
             energy  =-16.0d0                                                                                                   &
                  &   /5.0d0                                                                                                    &
                  &   *Pi**2                                                                                                    &
                  &   *(                                                                                                        &
                  &     -4.0d0                                                                                                  &
                  &     *(                                                                                                      &
                  &       +3.0d0*sqrt(radiusOuterScaleFree   *(1.0d0+radiusOuterScaleFree))                                     &
                  &       -4.0d0*sqrt(radiusOuterScaleFree**3*(1.0d0+radiusOuterScaleFree))                                     &
                  &       +8.0d0*sqrt(radiusOuterScaleFree**5*(1.0d0+radiusOuterScaleFree))                                     &
                  &      )                                                                                                      &
                  &     *asinh(sqrt(radiusOuterScaleFree))                                                                      &
                  &     +radiusOuterScaleFree*(7.0d0+4.0d0*radiusOuterScaleFree*(-3.0d0+8.0d0*radiusOuterScaleFree*log(2.0d0))) &
                  &     -4.0d0             *radiusOuterScaleFree**3 *log(radiusOuterScaleFree)                                  &
                  &     +5.0d0*(1.0d0+4.0d0*radiusOuterScaleFree**3)                                                            &
                  &     *log(1.0d0+radiusOuterScaleFree)                                                                        &
                  &   )
          case (specialCaseCoredNFW   %ID)
             analytic=.true.
             energy  =+Pi**2                                                                                          &
                  &   *(                                                                                              &
                  &     +radiusOuterScaleFree                                                                         &
                  &     *(                                                                                            &
                  &       +4.0d0+radiusOuterScaleFree                                                                 &
                  &       *(                                                                                          &
                  &         +10.0d0+radiusOuterScaleFree                                                              &
                  &         *(                                                                                        &
                  &           +35.0d0-4.0d0*Pi**2*(1.0d0+radiusOuterScaleFree)**3                                     &
                  &           + 6.0d0*radiusOuterScaleFree                                                            &
                  &           *(                                                                                      &
                  &             +9.0d0+4.0d0*radiusOuterScaleFree                                                     &
                  &            )                                                                                      &
                  &          )                                                                                        &
                  &        )                                                                                          &
                  &      )                                                                                            &
                  &     /(1.0d0+radiusOuterScaleFree)**3                                                              &
                  &     -4.0d0*log(1.0d0+radiusOuterScaleFree)                                                        &
                  &     *(                                                                                            &
                  &       +1.0d0+radiusOuterScaleFree                                                                 &
                  &       -3.0d0*radiusOuterScaleFree**2                                                              &
                  &       -6.0d0*radiusOuterScaleFree**3                                                              &
                  &       +3.0d0*radiusOuterScaleFree**3*(1.0d0+radiusOuterScaleFree)*log(1.0d0+radiusOuterScaleFree) &
                  &      )                                                                                            &
                  &     /(1.0d0+radiusOuterScaleFree)                                                                 &
                  &     -24.0d0*radiusOuterScaleFree**3*Dilogarithm(-radiusOuterScaleFree)                            &
                  &   )
          case default
             energy=+0.0d0
             call Error_Report('unknown special case'//{introspection:location})
          end select
       end select
    end select
    if (analytic) then
       energy=+energy                            &
            & *gravitationalConstant_internal    &
            & *self%scaleLength              **5 &
            & *self%densityNormalization     **2
    else
       energy=+self%energyKineticNumerical(radiusOuter,massDistributionEmbedding)
     end if
    return
  end function zhao1996EnergyKinetic
  
  subroutine zhao1996Descriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionZhao1996), intent(inout)           :: self
    type     (inputParameters         ), intent(inout)           :: descriptor
    logical                            , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                  )                          :: parameterLabel
    type     (inputParameters         )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','Zhao1996')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%scaleLength
    call parameters%addParameter('scaleLength'         ,trim(adjustl(parameterLabel)))
    return
  end subroutine zhao1996Descriptor
  
