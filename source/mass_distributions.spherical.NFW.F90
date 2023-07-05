!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implementation of an NFW \citep{navarro_structure_1996} mass distribution class.
  !!}

  use :: Numerical_Interpolation, only : interpolator
  
  !![
  <massDistribution name="massDistributionNFW">
   <description>An NFW \citep{navarro_structure_1996} mass distribution class.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionNFW
     !!{
     The NFW \citep{navarro_structure_1996} mass distribution.
     !!}
     private
     double precision                            :: densityNormalization                         , scaleLength
     double precision                            :: enclosedMassRadiusPrevious                   , enclosedMassPrevious
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
       <method method="timeFreefallTabulate" description="Tabulate the freefall time as a function of radius in a scale-free NFW mass distribution."/>
     </methods>
     !!]
     procedure :: massTotal                         => nfwMassTotal
     procedure :: density                           => nfwDensity
     procedure :: densityGradientRadial             => nfwDensityGradientRadial
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
     procedure :: potential                         => nfwPotential
     procedure :: energyPotential                   => nfwEnergyPotential
     procedure :: energyKinetic                     => nfwEnergyKinetic
     procedure :: descriptor                        => nfwDescriptor
  end type massDistributionNFW
  
  interface massDistributionNFW
     !!{
     Constructors for the {\normalfont \ttfamily nfw} mass distribution class.
     !!}
     module procedure massDistributionNFWConstructorParameters
     module procedure massDistributionNFWConstructorInternal
  end interface massDistributionNFW

contains

  function massDistributionNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily nfw} mass distribution class which builds the object from a parameter
    set.
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
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>0.5d0/Pi</defaultValue>
      <description>The density normalization of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scaleLength</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scale radius of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>concentration</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The concentration of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>virialRadius</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The virial radius of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the NFW profile is considered to be dimensionless.</description>
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
    !!{
    Internal constructor for ``nfw'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionNFW         )                          :: self
    double precision                              , intent(in   ), optional :: scaleLength         , concentration, &
         &                                                                     densityNormalization, mass         , &
         &                                                                     virialRadius
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree
    !![
    <constructorAssign variables="componentType, massType"/>
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
       self%densityNormalization=densityNormalization
    else if (                                   &
         &   present(mass                ).and. &
         &   present(virialRadius        )      &
         &  ) then
       radiusScaleFree          =+virialRadius/self%scaleLength
       self%densityNormalization=+mass/4.0d0/Pi/self%scaleLength**3/(log(1.0d0+radiusScaleFree)-radiusScaleFree/(1.0d0+radiusScaleFree))
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
    self%densityScaleFreeRadiusMinimum                =+2.0d0
    self%densityScaleFreeRadiusMaximum                =+0.5d0
    self%angularMomentumSpecificScaleFreeMinimum      =+huge(0.0d0)
    self%angularMomentumSpecificScaleFreeMaximum      =-huge(0.0d0)
    self%angularMomentumSpecificScaleFreeRadiusMinimum=+2.0d0
    self%angularMomentumSpecificScaleFreeRadiusMaximum=+0.5d0
    self%timeFreefallScaleFreeMinimum                 =+huge(0.0d0)
    self%timeFreefallScaleFreeMaximum                 =-huge(0.0d0)
    self%timeFreefallScaleFreeRadiusMinimum           =+2.0d0
    self%timeFreefallScaleFreeRadiusMaximum           =+0.5d0
    return
  end function massDistributionNFWConstructorInternal

  double precision function nfwMassTotal(self,componentType,massType)
    !!{
    Return the total mass in an NFW mass distribution.
    !!}
    implicit none
    class(massDistributionNFW         ), intent(inout)           :: self
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (self%matches(componentType,massType)) then
       nfwMassTotal=huge(0.0d0)
    else
       nfwMassTotal=0.0d0
    end if
    return
  end function nfwMassTotal

  double precision function nfwDensity(self,coordinates,componentType,massType)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an NFW mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW         ), intent(inout)           :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree

    if (.not.self%matches(componentType,massType)) then
       nfwDensity=0.0d0
       return
    end if    
    ! Compute the density at this position.
    radiusScaleFree=+coordinates%rSpherical          () &
         &          /self       %scaleLength
    nfwDensity     =+self       %densityNormalization   &
         &          /       radiusScaleFree             &
         &          /(1.0d0+radiusScaleFree)**2
    return
  end function nfwDensity

  double precision function nfwDensityGradientRadial(self,coordinates,logarithmic,componentType,massType) result(densityGradientRadial)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an NFW \citep{navarro_structure_1996} mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW         ), intent(inout), target   :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    logical                                       , intent(in   ), optional :: logarithmic
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    densityGradientRadial=0.0d0
    if (.not.self%matches(componentType,massType)) return
    radiusScaleFree      =+coordinates%rSpherical()         &
         &                /self       %scaleLength
    densityGradientRadial=-self       %densityNormalization &
         &                /self       %scaleLength          &
         &                /             radiusScaleFree **2 &
         &                *(1.0d0+3.0d0*radiusScaleFree)    &
         &                /(1.0d0+      radiusScaleFree)**3
    if (logarithmic_) densityGradientRadial=+            densityGradientRadial              &
         &                                  /self       %density              (coordinates) &
         &                                  *coordinates%rSpherical           (           )
    return
  end function nfwDensityGradientRadial

  double precision function nfwDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType) result(densityRadialMoment)
    !!{
    Computes radial moments of the density in an NFW \citep{navarro_structure_1996} mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW         ), intent(inout)           :: self
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    logical                                       , intent(  out), optional :: isInfinite
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radialMomentMinimum, radialMomentMaximum

    densityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    if (.not.self%matches(componentType,massType)) return
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

  double precision function nfwMassEnclosedBySphere(self,radius,componentType,massType) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for nfw mass distributions.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionNFW         ), intent(inout), target   :: self
    double precision                              , intent(in   )           :: radius
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree
    
    if (.not.self%matches(componentType,massType)) then
       mass=0.0d0
       return
    end if
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
  
  double precision function nfwRadiusEnclosingMass(self,mass,massFractional,componentType,massType) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for nfw mass distributions.
    !!}    
    use :: Numerical_Constants_Math, only : Pi
    use :: Lambert_Ws              , only : Lambert_W0
    use :: Error                   , only : Error_Report
    implicit none
    class           (massDistributionNFW         ), intent(inout), target   :: self
    double precision                              , intent(in   ), optional :: mass         , massFractional
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: mass_

    if (.not.self%matches(componentType,massType)) then
       radius=0.0d0
       return
    end if
    mass_=0.0d0
    if (present(mass)) then
       mass_=mass
    else if (present(massFractional)) then
       call Error_Report('mass is unbounded, so mass fraction is undefined'//{introspection:location})
    else
       call Error_Report('either mass or massFractional must be supplied'//{introspection:location})
    end if
    radius=-self%scaleLength                             &
         & *(                                            &
         &   +1.0d0                                      &
         &   /Lambert_W0(                                &
         &               -exp(                           &
         &                    -     1.0d0                &
         &                    -     mass_                &
         &                    /     4.0d0                &
         &                    /     Pi                   &
         &                    /self%densityNormalization &
         &                    /self%scaleLength**3 &
         &                   )                           &
         &              )                                &
         &   +1.0d0                                      &
         &  )
    return
  end function nfwRadiusEnclosingMass
  
  double precision function nfwRadiusEnclosingDensity(self,density,componentType,massType) result(radius)
    !!{
    Computes the radius enclosing a given mean density for nfw mass distributions.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionNFW         ), intent(inout), target       :: self
    double precision                              , intent(in   )               :: density
    type            (enumerationComponentTypeType), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional     :: massType
    double precision                              , allocatable  , dimension(:) :: radii                      , densities
    double precision                              , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                            :: densityScaleFree
    integer                                                                     :: countRadii

    densityScaleFree=+density                   &
         &           /self%densityNormalization
    if     (                                                 &
         &   densityScaleFree < self%densityScaleFreeMinimum &
         &  .or.                                             &
         &   densityScaleFree > self%densityScaleFreeMaximum &
         & ) then
       do while (densityEnclosedScaleFree(self%densityScaleFreeRadiusMinimum) < densityScaleFree)
          self%densityScaleFreeRadiusMinimum=0.5d0*self%densityScaleFreeRadiusMinimum
       end do
       do while (densityEnclosedScaleFree(self%densityScaleFreeRadiusMaximum) > densityScaleFree)
          self%densityScaleFreeRadiusMaximum=2.0d0*self%densityScaleFreeRadiusMaximum
       end do
       countRadii=int(log10(self%densityScaleFreeRadiusMaximum/self%densityScaleFreeRadiusMinimum)*countRadiiPerDecade)+1
       if (allocated(self%densityScaleFree_)) deallocate(self%densityScaleFree_)
       allocate(     radii            (countRadii))
       allocate(     densities        (countRadii))
       allocate(self%densityScaleFree_            )
       radii                       = Make_Range(self%densityScaleFreeRadiusMinimum,self%densityScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
       densities                   =-densityEnclosedScaleFree(           radii)
       self%densityScaleFreeMinimum= densities               (countRadii      )
       self%densityScaleFreeMaximum= densities               (         1      )
       self%densityScaleFree_      = interpolator            (densities ,radii)
    end if
    radius=+self%densityScaleFree_%interpolate(-densityScaleFree) &
         & *self%scaleLength
    return    
  end function nfwRadiusEnclosingDensity

  elemental double precision function massEnclosedScaleFree(radius) result(mass)
    !!{
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
    !!{
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
  
  double precision function nfwRadiusFromSpecificAngularMomentum(self,angularMomentumSpecific,componentType,massType) result(radius)
    !!{
    Computes the radius corresponding to a given specific angular momentum for nfw mass distributions.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Ranges                , only : Make_Range                     , rangeTypeLogarithmic
    implicit none
    class           (massDistributionNFW         ), intent(inout), target       :: self
    double precision                              , intent(in   )               :: angularMomentumSpecific
    type            (enumerationComponentTypeType), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional     :: massType
    double precision                              , allocatable  , dimension(:) :: radii                                   , angularMomentaSpecific
    double precision                              , parameter                   :: countRadiiPerDecade             =100.0d0
    double precision                                                            :: angularMomentumSpecificScaleFree
    integer                                                                     :: countRadii

    angularMomentumSpecificScaleFree=+angularMomentumSpecific                  &
         &                           /sqrt(                                    &
         &                                 +gravitationalConstantGalacticus    &
         &                                 *self%densityNormalization          &
         &                                )                                    &
         &                           /      self%scaleLength               **2
    if     (                                                                                 &
         &   angularMomentumSpecificScaleFree < self%angularMomentumSpecificScaleFreeMinimum &
         &  .or.                                                                             &
         &   angularMomentumSpecificScaleFree > self%angularMomentumSpecificScaleFreeMaximum &
         & ) then
       do while (angularMomentumSpecificEnclosedScaleFree(self%angularMomentumSpecificScaleFreeRadiusMinimum) > angularMomentumSpecificScaleFree)
          self%angularMomentumSpecificScaleFreeRadiusMinimum=0.5d0*self%angularMomentumSpecificScaleFreeRadiusMinimum
       end do
       do while (angularMomentumSpecificEnclosedScaleFree(self%angularMomentumSpecificScaleFreeRadiusMaximum) < angularMomentumSpecificScaleFree)
          self%angularMomentumSpecificScaleFreeRadiusMaximum=2.0d0*self%angularMomentumSpecificScaleFreeRadiusMaximum
       end do
       countRadii=int(log10(self%angularMomentumSpecificScaleFreeRadiusMaximum/self%angularMomentumSpecificScaleFreeRadiusMinimum)*countRadiiPerDecade)+1
       if (allocated(self%angularMomentumSpecificScaleFree_)) deallocate(self%angularMomentumSpecificScaleFree_)
       allocate(     radii                            (countRadii))
       allocate(     angularMomentaSpecific           (countRadii))
       allocate(self%angularMomentumSpecificScaleFree_            )
       radii                                       =Make_Range(self%angularMomentumSpecificScaleFreeRadiusMinimum,self%angularMomentumSpecificScaleFreeRadiusMaximum,countRadii,rangeTypeLogarithmic)
       angularMomentaSpecific                      =angularMomentumSpecificEnclosedScaleFree(                       radii)
       self%angularMomentumSpecificScaleFreeMinimum=angularMomentaSpecific                  (            countRadii      )
       self%angularMomentumSpecificScaleFreeMaximum=angularMomentaSpecific                  (                     1      )
       self%angularMomentumSpecificScaleFree_      =interpolator                            (angularMomentaSpecific,radii)
    end if
    radius=+self%angularMomentumSpecificScaleFree_%interpolate(angularMomentumSpecificScaleFree) &
         & *self%scaleLength
    return    
  end function nfwRadiusFromSpecificAngularMomentum

  elemental double precision function angularMomentumSpecificEnclosedScaleFree(radius) result(angularMomentumSpecific)
    !!{
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

  double precision function nfwVelocityRotationCurveMaximum(self,componentType,massType) result(velocity)
    !!{
    Return the peak velocity in the rotation curve for an nfw mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionNFW         ), intent(inout)           :: self
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType    
    double precision                              , parameter               :: circularVelocityMaximumScaleFree=0.4649909628174221d0 ! The circular velocity (in scale-free units) at the peak of the NFW rotation curve.
    !                                                                                                                                  Numerical value found using Mathematica.

    if (.not.self%matches(componentType,massType)) then
       velocity=0.0d0
       return
    end if
    velocity=+circularVelocityMaximumScaleFree             &
         &   *sqrt(                                        &
         &         +4.0d0                                  &
         &         *Pi                                     &
         &         *self%densityNormalization              &
         &        )                                        &
         &   *      self%scaleLength
    if (.not.self%isDimensionless())                       &
         & velocity=+velocity                              &
         &          *sqrt(gravitationalConstantGalacticus)
    return
  end function nfwVelocityRotationCurveMaximum

  double precision function nfwRadiusRotationCurveMaximum(self,componentType,massType) result(radius)
    !!{
    Return the peak velocity in the rotation curve for an nfw mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW         ), intent(inout), target   :: self
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    ! The radius (in scale-free units) at the peak of the NFW rotation curve. Numerical value found using Mathematica.
    double precision                              , parameter               :: radiusCircularVelocityMaximumScaleFree=2.162581587064612d0
    
    if (.not.self%matches(componentType,massType)) then
       radius=0.0d0
       return
    end if
    radius=+radiusCircularVelocityMaximumScaleFree &
         & *self%scaleLength
    return
  end function nfwRadiusRotationCurveMaximum

  double precision function nfwPotential(self,coordinates,componentType,massType,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in an nfw mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess      , structureErrorCodeInfinite
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Error                           , only : Error_Report
    implicit none
    class           (massDistributionNFW              ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                                             :: radiusScaleFree
    
    if (present(status)) status=structureErrorCodeSuccess
    if (.not.self%matches(componentType,massType)) then
       potential=0.0d0
       return
    end if
    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %scaleLength
    potential=+potentialScaleFree       (radiusScaleFree)    &
         &    *self%densityNormalization                     &
         &    *self%scaleLength                          **2
    if (.not.self%isDimensionless()) potential=+gravitationalConstantGalacticus &
         &                                     *potential
    return
  end function nfwPotential

  elemental double precision function potentialScaleFree(radius) result(potential)
    !!{
    Compute the potential in a scale-free NFW mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    double precision, parameter     :: radiusSmall=1.0d-10
    double precision                :: radiusTerm
    
    if (radius < radiusSmall) then
       ! Use a series solution for very small radii.
       radiusTerm=1.0d0-0.5d0*radius
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
    !!{
    Compute the potential difference in a scale-free NFW mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Comparison    , only : Values_Agree
    implicit none
    double precision, intent(in   ) :: radius1                             , radius2
    double precision, parameter     :: radiusSmall                 =1.0d-10
    double precision, parameter     :: toleranceRelative           =1.0d-03
    double precision                :: potentialGradientLogarithmic        , radiusDifferenceLogarithmic
    
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
  
  double precision function nfwFourierTransform(self,radiusOuter,wavenumber,componentType,massType) result(fourierTransform)
    !!{
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in an NFW mass
    distribution, using the expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    !!}
    use :: Exponential_Integrals, only : Cosine_Integral, Sine_Integral
    implicit none
    class           (massDistributionNFW         ), intent(inout)           :: self
    double precision                              , intent(in   )           :: radiusOuter        , wavenumber
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: wavenumberScaleFree, radiusOuterScaleFree

    fourierTransform=0.0d0
    if (.not.self%matches(componentType,massType)) return
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
  
  double precision function nfwRadiusFreefall(self,time,componentType,massType) result(radius)
    !!{
    Compute the freefall radius at the given {\normalfont \ttfamily time} in an NFW mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr, gravitationalConstantGalacticus
    implicit none
    class           (massDistributionNFW         ), intent(inout)               :: self
    double precision                              , intent(in   )               :: time
    type            (enumerationComponentTypeType), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional     :: massType
    double precision                                                            :: timeScaleFree, timeScale
    
    if (.not.self%matches(componentType,massType)) then
       radius=0.0d0
       return
    end if
    timeScale    =+1.0d0/sqrt(                                 &
         &              +gravitationalConstantGalacticus &
         &              *self%densityNormalization       &
         &             )                                 &
         &        *Mpc_per_km_per_s_To_Gyr
    timeScaleFree=+time                                  &
         &        /timeScale
    call self%timeFreefallTabulate(timeScaleFree)
    radius=+self%timeFreefallScaleFree_%interpolate(timeScaleFree) &
         & *self%scaleLength
    return   
  end function nfwRadiusFreefall
  
  double precision function nfwRadiusFreefallIncreaseRate(self,time,componentType,massType) result(radiusIncreaseRate)
    !!{
    Compute the rate of increase of the freefall radius at the given {\normalfont \ttfamily time} in an nfw mass
    distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr, gravitationalConstantGalacticus
    implicit none
    class           (massDistributionNFW         ), intent(inout)           :: self
    double precision                              , intent(in   )           :: time
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: timeScaleFree, timeScale

    if (.not.self%matches(componentType,massType)) then
       radiusIncreaseRate=0.0d0
       return
    end if
    timeScale    =+1.0d0/sqrt(                                 &
         &              +gravitationalConstantGalacticus &
         &              *self%densityNormalization       &
         &             )                                 &
         &        *Mpc_per_km_per_s_To_Gyr
    timeScaleFree=+time                                  &
         &        *timeScale
    call self%timeFreefallTabulate(timeScaleFree)
    radiusIncreaseRate=+self%timeFreefallScaleFree_%derivative(timeScaleFree) &
         &             *self%scaleLength                                      &
         &             /     timeScale
    return
  end function nfwRadiusFreefallIncreaseRate
  
  subroutine nfwTimeFreefallTabulate(self,timeScaleFree)
    !!{
    Tabulate the freefall radius at the given {\normalfont \ttfamily time} in an NFW mass distribution.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges     , only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionNFW), intent(inout)               :: self
    double precision                     , intent(in   )               :: timeScaleFree
    double precision                     , allocatable  , dimension(:) :: radii                      , timesFreefall
    double precision                     , parameter                   :: countRadiiPerDecade=100.0d0
    double precision                                                   :: radiusStart
    integer                                                            :: countRadii                 , i
    type            (integrator         )                              :: integrator_

    if     (                                                   &
         &   timeScaleFree < self%timeFreefallScaleFreeMinimum &
         &  .or.                                               &
         &   timeScaleFree > self%timeFreefallScaleFreeMaximum &
         & ) then
       integrator_=integrator(timeFreeFallIntegrand,toleranceRelative=1.0d-6)
       do while (timeFreefallScaleFree(self%timeFreefallScaleFreeRadiusMinimum) > timeScaleFree)
          self%timeFreefallScaleFreeRadiusMinimum=0.5d0*self%timeFreefallScaleFreeRadiusMinimum
       end do
       do while (timeFreefallScaleFree(self%timeFreefallScaleFreeRadiusMaximum) < timeScaleFree)
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
       self%timeFreefallScaleFreeMinimum=timesFreefall(   countRadii      )
       self%timeFreefallScaleFreeMaximum=timesFreefall(            1      )
       self%timeFreefallScaleFree_      =interpolator (timesFreefall,radii)
    end if
    return
    
  contains
    
    double precision function timeFreefallScaleFree(radius)
      !!{
      Evaluate the freefall time from a given radius in a scale-free NFW mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      radiusStart          =                            radius
      timeFreefallScaleFree=integrator_%integrate(0.0d0,radius)
      return
    end function timeFreefallScaleFree
    
    double precision function timeFreeFallIntegrand(radius)
      !!{
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
    !!{
    Compute the potential energy within a given {\normalfont \ttfamily radius} in an NFW mass distribution. This is
    \begin{eqnarray}
      W &=& - \frac{\mathrm{G}}{2} \rho_0^2 r_\mathrm{s}^5 \int_0^{x_\mathrm{out}} \frac{m^2(x)}{x^2} \mathrm{d} x, \nonumber \\
        &-& - \frac{\mathrm{G}}{2} \rho_0^2 r_\mathrm{s}^5 \left[ \frac{x}{1+x} - \frac{\log^2(1+x)}{x} + \frac{\left\{\log(1+x)-x/(1+x)\right\}^2}{x} \right],
    \end{eqnarray}
    where $x=r/r_\mathrm{s}$ and $m(x)$ is the scale-free mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionNFW), intent(inout) :: self
    double precision                     , intent(in   ) :: radiusOuter
    double precision                                     :: radiusOuterScaleFree
    
    radiusOuterScaleFree=+     radiusOuter                                                                                   &
         &               /self%scaleLength
    energy=-gravitationalConstantGalacticus                                                                                  &
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
    !!{
    Compute the kinetic energy within a given {\normalfont \ttfamily radius} in an NFW mass distribution. This is
    \begin{eqnarray}
      T &=& 6 \pi \mathrm{G} \rho_0^2 r_\mathrm{s}^5 \int_0^{x_\mathrm{out}} \rho(x) \sigma^2(x) x^2 \mathrm{d} x, \nonumber \\
        &=& \pi \mathrm{G} \rho_0^2 r_\mathrm{s}^5 \left[ 6 x^3 \text{Li}_2(-x)+x^3 (-\log (x))+\log (x+1) \left(3 x^3 \log (x+1)+((x-6) x+3) x-2\right)+\left(x \left(\pi ^2 x-7\right)+5\right) x+\frac{3}{x+1} \right],
    \end{eqnarray}
    where $x=r/r_\mathrm{s}$, $\rho(x)$ is the scale-free density, and $\sigma^2(x)$ is the scale-free velocity dispersion.
    !!}
    use :: Dilogarithms                    , only : Dilogarithm
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
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
          energy              =+gravitationalConstantGalacticus                                                       &
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

  subroutine nfwDescriptor(self,descriptor,includeClass)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionNFW), intent(inout)           :: self
    type     (inputParameters    ), intent(inout)           :: descriptor
    logical                       , intent(in   ), optional :: includeClass
    character(len=18             )                          :: parameterLabel
    type     (inputParameters    )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','NFW')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%scaleLength
    call parameters%addParameter('scaleLength'         ,trim(adjustl(parameterLabel)))
    return
  end subroutine nfwDescriptor

