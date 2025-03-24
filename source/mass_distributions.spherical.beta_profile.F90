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
  Implementation of a $\beta$-profile mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionBetaProfile">
   <description>An mass distribution class for $\beta$-profile distributions.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionBetaProfile
     !!{
     The $\beta$-profile: $\rho(r)=\rho_0/[1+(r/r_\mathrm{core})^2]^{3\beta/2}$
     !!}
     double precision :: beta                  , coreRadius           , densityNormalization  , &
          &              momentRadial2Previous , momentRadial3Previous, momentRadial2XPrevious, &
          &              momentRadial3XPrevious, outerRadius
     logical          :: betaIsTwoThirds       , truncateAtOuterRadius
   contains
     procedure :: density               => betaProfileDensity
     procedure :: densityGradientRadial => betaProfileDensityGradientRadial
     procedure :: densityRadialMoment   => betaProfileDensityRadialMoment
     procedure :: densitySquareIntegral => betaProfileDensitySquareIntegral
     procedure :: massEnclosedBySphere  => betaProfileMassEnclosedBySphere
     procedure :: potentialIsAnalytic   => betaProfilePotentialIsAnalytic
     procedure :: potential             => betaProfilePotential
     procedure :: descriptor            => betaProfileDescriptor
  end type massDistributionBetaProfile

  interface massDistributionBetaProfile
     !!{
     Constructors for the {\normalfont \ttfamily betaProfile} mass distribution class.
     !!}
     module procedure betaProfileConstructorParameters
     module procedure betaProfileConstructorInternal
  end interface massDistributionBetaProfile

contains

  function betaProfileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily betaProfile} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionBetaProfile)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: beta         , densityNormalization , &
         &                                                          mass         , outerRadius          , &
         &                                                          coreRadius
    logical                                                      :: dimensionless, truncateAtOuterRadius
    type            (varying_string             )                :: componentType
    type            (varying_string             )                :: massType

    !![
    <inputParameter>
      <name>beta</name>
      <defaultValue>2.0d0/3.0d0</defaultValue>
      <description>The value $\beta$ in a $\beta$-model mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The density normalization of a $\beta$-model mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The mass of a $\beta$-model mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outerRadius</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The outer radius of a $\beta$-model mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>coreRadius</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The core radius of a $\beta$-model mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.false.</defaultValue>
      <description>If true then the $\beta$-model mass distribution is considered to be in dimensionless units.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>truncateAtOuterRadius</name>
      <defaultValue>.false.</defaultValue>
      <description>If true then the $\beta$-model mass distribution is truncated beyond the outer radius.</description>
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
     <call>self=massDistributionBetaProfile(beta,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="densityNormalization"  value="densityNormalization"  parameterPresent="parameters"/>
     <argument name="mass"                  value="mass"                  parameterPresent="parameters"/>
     <argument name="outerRadius"           value="outerRadius"           parameterPresent="parameters"/>
     <argument name="coreRadius"            value="coreRadius"            parameterPresent="parameters"/>
     <argument name="dimensionless"         value="dimensionless"         parameterPresent="parameters"/>
     <argument name="truncateAtOuterRadius" value="truncateAtOuterRadius" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function betaProfileConstructorParameters

  function betaProfileConstructorInternal(beta,densityNormalization,mass,outerRadius,coreRadius,dimensionless,truncateAtOuterRadius,componentType,massType) result(self)
    !!{
    Constructor for {\normalfont \ttfamily betaProfile} convergence class.
    !!}
    use :: Display                 , only : displayIndent      , displayMessage, displayUnindent, displayVerbosity, &
          &                                 verbosityLevelDebug
    use :: Error                   , only : Error_Report
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Comparison    , only : Values_Agree       , Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionBetaProfile )                          :: self
    double precision                              , intent(in   )           :: beta
    double precision                              , intent(in   ), optional :: densityNormalization              , mass                      , &
         &                                                                     outerRadius                       , coreRadius
    logical                                       , intent(in   ), optional :: dimensionless                     , truncateAtOuterRadius
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                              , parameter               :: radiusTiny                  =1.0d-6
    double precision                                                        :: r
    character       (len=64                      )                          :: message
    double precision                              , save                    :: radiusCoreFractionalPrevious       , normalizationFactorStored
    !$omp threadprivate(radiusCoreFractionalPrevious,normalizationFactorStored)
    !![
    <constructorAssign variables="beta, densityNormalization, coreRadius, componentType, massType"/>
    !!]

    ! Check for special case of β=2/3.
    self%betaIsTwoThirds=Values_Agree(self%beta,2.0d0/3.0d0,relTol=1.0d-3)
    ! Determine if profile is dimensionless.
    self%dimensionless=.false.
    if (present(dimensionless)) self%dimensionless=dimensionless
    ! If dimensionless, then set scale length and mass to unity.
    if (self%dimensionless) then
       if (present(coreRadius          )) then
          if (Values_Differ(coreRadius          ,1.0d0,absTol=1.0d-6)) call Error_Report('coreRadius should be unity for a dimensionless profile (or simply do not specify a scale length)'                  //{introspection:location})
       end if
       if (present(densityNormalization)) then
          if (Values_Differ(densityNormalization,1.0d0,absTol=1.0d-6)) call Error_Report('densityNormalization should be unity for a dimensionless profile (or simply do not specify a densityNormalization)'//{introspection:location})
       end if
       if (present(mass                ))                              call Error_Report('mass cannot be specified for a dimensionless profile'                                                              //{introspection:location})
       if (present(outerRadius         ))                              call Error_Report('outer radius cannot be specified for a dimensionless profile'                                                      //{introspection:location})
       self%coreRadius          =1.0d0
       self%densityNormalization=1.0d0
    else
       ! Set core radius.
       if (.not.present(coreRadius)) call Error_Report('core radius must be specified for dimensionful profiles'//{introspection:location})
       self%coreRadius=coreRadius
       ! Determine density normalization.
       if      (                                   &
            &   present(densityNormalization)      &
            &  ) then
          self%densityNormalization=densityNormalization
       else if (                                   &
            &   present(mass                ).and. &
            &   present(outerRadius         )      &
            &  ) then
          r=outerRadius/coreRadius
          if (outerRadius > 0.0d0) then
             if (self%betaIsTwoThirds) then
                if (r /= radiusCoreFractionalPrevious) then
                   radiusCoreFractionalPrevious=r
                   if (r < radiusTiny) then
                      ! Use a series solution.
                      normalizationFactorStored=3.0d0/r**3+9.0d0/5.0d0/r-36.0d0*r/175.0d0
                   else
                      ! Use the full solution.
                      normalizationFactorStored=1.0d0/(r-atan(r))
                   end if
                end if
                self%densityNormalization=      mass/4.0d0/Pi/coreRadius **3*normalizationFactorStored
             else
                self%densityNormalization=3.0d0*mass/4.0d0/Pi/outerRadius**3/Hypergeometric_2F1([1.5d0,1.5d0*beta],[2.5d0],-r**2)
             end if
          else
             call Error_Report('unphysical outer radius'//{introspection:location})
          end if
          ! Assert that the mass within the outer radius equals that specified.
          if (displayVerbosity() >= verbosityLevelDebug) then
             if (.not.Values_Agree(self%massEnclosedBySphere(outerRadius),mass,relTol=1.0d-6,absTol=tiny(0.0d0))) then
                call displayIndent('beta-profile parameters:')
                write (message,'(a,e12.6)') '    coreRadius: ',coreRadius
                call displayMessage(message)
                write (message,'(a,e12.6)') '   outerRadius: ',outerRadius
                call displayMessage(message)
                write (message,'(a,e12.6)') '          mass: ',mass
                call displayMessage(message)
                write (message,'(a,e12.6)') '          beta: ',beta
                call displayMessage(message)
                write (message,'(a,e12.6)') 'mass(<r_outer): ',self%massEnclosedBySphere(outerRadius)
                call displayMessage(message)
                call displayUnindent('done')
                call Error_Report('profile normalization failed'//{introspection:location})
             end if
          end if
       end if
    end if
    ! Check for truncation.
    if (present(truncateAtOuterRadius)) then
       self%truncateAtOuterRadius=truncateAtOuterRadius
    else
       self%truncateAtOuterRadius=.false.
    end if
    if (self%truncateAtOuterRadius) then
       if (.not.present(outerRadius)) call Error_Report('can not truncate profile without an outer radius'//{introspection:location})
       self%outerRadius=outerRadius
    end if
    ! Initialize stored results.
    self%momentRadial2XPrevious=-1.0d0
    self%momentRadial3XPrevious=-1.0d0
    return
  end function betaProfileConstructorInternal

  double precision function betaProfileDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a $\beta$-profile mass distribution.
    !!}
    implicit none
    class           (massDistributionBetaProfile ), intent(inout) :: self
    class           (coordinate                  ), intent(in   ) :: coordinates
    double precision                                              :: radius

    ! Compute density.
    radius=coordinates%rSpherical()
    if (self%truncateAtOuterRadius .and. radius > self%outerRadius) then
       betaProfileDensity=0.0d0
    else
       betaProfileDensity=self%densityNormalization/(1.0d0+(radius/self%coreRadius)**2)**(1.5d0*self%beta)
    end if
    return
  end function betaProfileDensity

  double precision function betaProfileDensityGradientRadial(self,coordinates,logarithmic)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a $\beta$-profile mass distribution.
    !!}
    implicit none
    class           (massDistributionBetaProfile ), intent(inout), target   :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    logical                                       , intent(in   ), optional :: logarithmic
    double precision                                                        :: radius
    logical                                                                 :: logarithmicActual

    ! Set default options.
    logarithmicActual=.false.
    if (present(logarithmic)) logarithmicActual=logarithmic
    ! Get position in spherical coordinate system.
    radius=coordinates%rSpherical()
    ! Apply truncation.
    if (self%truncateAtOuterRadius .and. radius > self%outerRadius) then
       betaProfileDensityGradientRadial=0.0d0
       return
    end if
    ! Convert to dimensionless radius.
    radius=radius/self%coreRadius
    ! Compute density gradient.
    if (logarithmicActual) then
       betaProfileDensityGradientRadial=                  &
            & -3.0d0                                      &
            & *self%beta                                  &
            & * radius**2                                 &
            & /(radius**2+1.0d0)
    else
       betaProfileDensityGradientRadial=                  &
            & -3.0d0                                      &
            & *self%beta                                  &
            & *self%densityNormalization                  &
            & /self%coreRadius                            &
            & * radius                                    &
            & /(radius**2+1.0d0)**(1.5d0*self%beta+1.0d0)
    end if
    return
  end function betaProfileDensityGradientRadial

  double precision function betaProfileMassEnclosedBySphere(self,radius)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for $\beta$-profile mass distributions. Result computed
    using \href{http://www.wolframalpha.com/input/?i=integrate+4*pi*r^2*rho\%2F\%281\%2Br^2\%29^\%283*beta\%2F2\%29}{Wolfram Alpha}.
    !!}
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionBetaProfile ), intent(inout), target :: self
    double precision                              , intent(in   )         :: radius
    double precision                              , parameter             :: radiusTiny      =1.0d-6
    double precision                                                      :: fractionalRadius       , radius_

    if (self%truncateAtOuterRadius .and. radius > self%outerRadius) then
       radius_=self%outerRadius
    else
       radius_=radius
    end if
    fractionalRadius=radius_/self%coreRadius
    if (self%betaIsTwoThirds) then
       ! Solution for special case of β=2/3.
       if (fractionalRadius < radiusTiny) then
          ! Use a series solution.
          betaProfileMassEnclosedBySphere=                 &
               & +4.0d0                                    &
               & *Pi                                       &
               & *self%densityNormalization                &
               & *self%coreRadius                 **3      &
               & *                fractionalRadius**3      &
               & *(  +1.0d0/3.0d0+fractionalRadius**2      &
               & * ( -1.0d0/5.0d0+fractionalRadius**2      &
               & *  (+1.0d0/7.0d0                          &
               &    )                                      &
               &   )                                       &
               &  )
       else
          betaProfileMassEnclosedBySphere=                 &
               & +4.0d0                                    &
               & *Pi                                       &
               & *self%densityNormalization                &
               & *(                                        &
               &   +     fractionalRadius                  &
               &   -atan(fractionalRadius)                 &
               &  )                                        &
               & *self%coreRadius**3
       end if
    else
       ! General solution.
       betaProfileMassEnclosedBySphere=                    &
            & +4.0d0                                       &
            & /3.0d0                                       &
            & *Pi                                          &
            & *self%densityNormalization                   &
            & *radius_**3                                  &
            & *Hypergeometric_2F1(                         &
            &                     [1.5d0,1.5d0*self%beta], &
            &                     [2.5d0                ], &
            &                     -fractionalRadius**2     &
            &                    )
    end if
    return
  end function betaProfileMassEnclosedBySphere

  logical function betaProfilePotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionBetaProfile), intent(inout) :: self

    isAnalytic=.true.
    return
  end function betaProfilePotentialIsAnalytic

  double precision function betaProfilePotential(self,coordinates,status)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a $\beta$-profile mass distribution. Calculated using
    \href{http://www.wolframalpha.com/input/?i=integrate+4\%2F3+\%CF\%80+r+\%CF\%81+2F1\%283\%2F2\%2C+\%283+\%CE\%B2\%29\%2F2\%2C+5\%2F2\%2C+-r^2\%29}{Wolfram
    Alpha}.
    !!}
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Hypergeometric_Functions        , only : Hypergeometric_2F1
    use :: Numerical_Comparison            , only : Values_Agree
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionBetaProfile      ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                   , parameter               :: fractionalRadiusMinimum=1.0d-3
    double precision                                                             :: fractionalRadius              , radius

    if (present(status)) status=structureErrorCodeSuccess
    ! Compute the potential at this position.
    radius=coordinates%rSpherical()
    if (self%truncateAtOuterRadius .and. radius > self%outerRadius) then
       fractionalRadius=self%outerRadius/self%coreRadius
    else
       fractionalRadius=          radius/self%coreRadius
    end if
    if (Values_Agree(self%beta,2.0d0/3.0d0,absTol=1.0d-6)) then
       if (fractionalRadius < fractionalRadiusMinimum) then
          betaProfilePotential=                     &
               &  Pi                                &
               & *self%densityNormalization         &
               & *(                                 &
               &    4.0d0                           &
               &   +2.0d0                           &
               &   /3.0d0                           &
               &   *fractionalRadius**2             &
               &   -fractionalRadius**4             &
               &   /5.0d0                           &
               &  )
       else
          betaProfilePotential=                     &
               &  2.0d0                             &
               & *Pi                                &
               & *self%densityNormalization         &
               & *(                                 &
               &    log (                           &
               &          1.0d0                     &
               &         +fractionalRadius**2       &
               &        )                           &
               &   +2.0d0                           &
               &   *atan( fractionalRadius   )      &
               &   /      fractionalRadius          &
               &  )
       end if
    else
       if (fractionalRadius < fractionalRadiusMinimum) then
          betaProfilePotential=                     &
               &  Pi                                &
               & *self%densityNormalization         &
               & *fractionalRadius**3               &
               & *(                                 &
               &    4.0d0/3.0d0                     &
               &   -6.0d0                           &
               &   /5.0d0                           &
               &   *self%beta                       &
               &   *fractionalRadius**2             &
               &  )
       else
          betaProfilePotential=                                       &
               &  2.0d0                                               &
               & /3.0d0                                               &
               & *Pi                                                  &
               & *self%densityNormalization                           &
               & *     fractionalRadius**2                            &
               & *(                                                   &
               &    6.0d0                                             &
               &   *  (fractionalRadius**2+1.0d0)**(-1.5d0*self%beta) &
               &   *(                                                 &
               &      (fractionalRadius**2+1.0d0)**(+1.5d0*self%beta) &
               &     -(fractionalRadius**2+1.0d0)                     &
               &    )                                                 &
               &   /(3.0d0*self%beta-2.0d0)                           &
               &   /   fractionalRadius**2                            &
               &   -2.0d0                                             &
               &   *Hypergeometric_2F1(                               &
               &                       [1.5d0,1.5d0*self%beta],       &
               &                       [2.5d0                ],       &
               &                       -fractionalRadius**2           &
               &                      )                               &
               &  )
       end if
    end if
    if (.not.self%isDimensionless())         &
         & betaProfilePotential=             &
         &   betaProfilePotential            &
         &   *gravitationalConstant_internal
    return
  end function betaProfilePotential

  double precision function betaProfileDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Computes radial moments of the density in a $\beta$-profile mass distribution.
    !!}
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Comparison    , only : Values_Agree
    use :: Error                   , only : Error_Report
    use :: Gamma_Functions         , only : Gamma_Function
    implicit none
    class           (massDistributionBetaProfile ), intent(inout)           :: self
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum          , radiusMaximum
    logical                                       , intent(  out), optional :: isInfinite
    logical                                                                 :: haveRadiusMinimum      , haveRadiusMaximum
    double precision                                                        :: radiusMinimum_         , radiusMaximum_         , &
         &                                                                     fractionalRadiusMinimum, fractionalRadiusMaximum
    integer                                                                 :: specialCaseMoment

    ! Determine effective radii.
    haveRadiusMinimum=present(radiusMinimum)
    haveRadiusMaximum=present(radiusMaximum)
    radiusMinimum_   =                0.0d0
    radiusMaximum_   =huge   (        0.0d0)
    if (haveRadiusMinimum) then
       radiusMinimum_   =    radiusMinimum
    end if
    if (haveRadiusMaximum) then
       radiusMaximum_   =    radiusMaximum
    end if
    if (self%truncateAtOuterRadius) then
       radiusMinimum_   =min(radiusMinimum_,self%outerRadius)
       radiusMaximum_   =min(radiusMaximum_,self%outerRadius)
       haveRadiusMinimum=.true.
       haveRadiusMaximum=.true.
    end if
    ! Determine if special case solutions can be used.
    specialCaseMoment=-huge(0)
    if (self%betaIsTwoThirds) then
       if      (Values_Agree(moment,0.0d0,absTol=1.0d-6)) then
          specialCaseMoment=0
       else if (Values_Agree(moment,1.0d0,absTol=1.0d-6)) then
          specialCaseMoment=1
       else if (Values_Agree(moment,2.0d0,absTol=1.0d-6)) then
          specialCaseMoment=2
       else if (Values_Agree(moment,3.0d0,absTol=1.0d-6)) then
          specialCaseMoment=3
       end if
    end if
    if (present(isInfinite)) isInfinite=.false.
    if (haveRadiusMaximum) then
       fractionalRadiusMaximum=radiusMaximum_/self%coreRadius
       if (specialCaseMoment /= -huge(0)) then
          ! Special case for 0ᵗʰ, 1ˢᵗ, 2ⁿᵈ, and 3ʳᵈ moments of a β=2/3 distribution.
          betaProfileDensityRadialMoment=                    &
               & +radialMomentTwoThirds(specialCaseMoment,fractionalRadiusMaximum)
       else
          ! General solution.
          betaProfileDensityRadialMoment=                                         &
               & +fractionalRadiusMaximum**(moment+1.0d0)                         &
               & *Hypergeometric_2F1     (                                        &
               &                          [(moment+1.0d0)/2.0d0,1.5d0*self%beta], &
               &                          [(moment+3.0d0)/2.0d0                ], &
               &                          -fractionalRadiusMaximum**2             &
               &                         )                                        &
               & /                         (moment+1.0d0)
       end if
    else
       if (moment >= 3.0d0*self%beta-1.0d0) then
          betaProfileDensityRadialMoment=0.0d0
          if (present(isInfinite)) then
             isInfinite=.true.
             return
          else
             call Error_Report('radial moment is infinite at r=0 for m ≥ 3β-1'//{introspection:location})
          end if
       else
          betaProfileDensityRadialMoment=+Gamma_Function((+3.0d0+moment                )/2.0d0) &
               &                         *Gamma_Function((-1.0d0-moment+3.0d0*self%beta)/2.0d0) &
               &                         /Gamma_Function((             +3.0d0*self%beta)/2.0d0) &
               &                         /               (+1.0d0+moment                ) 
       end if
    end if
    if (haveRadiusMinimum) then
       fractionalRadiusMinimum=radiusMinimum_/self%coreRadius
       if (specialCaseMoment /= -huge(0)) then
          ! Special case for 0ᵗʰ, 1ˢᵗ, 2ⁿᵈ, and 3ʳᵈ moments of a β=2/3 distribution.
          betaProfileDensityRadialMoment=                                          &
               & +betaProfileDensityRadialMoment                                   &
               & -radialMomentTwoThirds(specialCaseMoment,fractionalRadiusMinimum)
       else
          betaProfileDensityRadialMoment=                                         &
               & +betaProfileDensityRadialMoment                                  &
               & -fractionalRadiusMinimum**(moment+1.0d0)                         &
               & *Hypergeometric_2F1     (                                        &
               &                          [(moment+1.0d0)/2.0d0,1.5d0*self%beta], &
               &                          [(moment+3.0d0)/2.0d0                ], &
               &                          -fractionalRadiusMinimum**2             &
               &                         )                                        &
               & /                         (moment+1.0d0)
       end if
    else
       if (moment <= -1.0d0) then
          betaProfileDensityRadialMoment=0.0d0
          if (present(isInfinite)) then
             isInfinite=.true.
             return
          else
             call Error_Report('radial moment is infinite at r=0 for m ≤ -1'//{introspection:location})
          end if
       end if
    end if
    ! Convert to dimensionful units.
    betaProfileDensityRadialMoment                         &
         & =betaProfileDensityRadialMoment                 &
         & *self%densityNormalization                      &
         & *self%coreRadius               **(1.0d0+moment)
    return

  contains

    double precision function radialMomentTwoThirds(moment,x)
      !!{
      Special case of radial moment for $\beta=2/3$ $\beta$-profile.
      !!}
      use :: Error, only : Error_Report
      implicit none
      integer         , intent(in   ) :: moment
      double precision, intent(in   ) :: x
      double precision, parameter     :: xSmall=1.0d-6

      if (x <= 0.0d0) then
         ! For zero radius the moment is always zero.
         radialMomentTwoThirds=0.0d0
      else
         ! For non-zero radius, compute the moment.
         if      (moment == 0) then
            ! 0ᵗʰ moment.
            radialMomentTwoThirds=atan(x)
         else if (moment == 1) then
            ! 1ˢᵗ moment.
            radialMomentTwoThirds=0.5d0*log(1.0d0+x**2)
         else if (moment == 2) then
            ! 2ⁿᵈ moment.
            if (x /= self%momentRadial2XPrevious) then
               self%momentRadial2XPrevious=x
               if (x < xSmall) then
                  ! Series solution for small x.
                  self%momentRadial2Previous=x**3*(1.0d0/3.0d0-x**2/5.0d0)
               else
                  ! General solution.
                  self%momentRadial2Previous=+     x &
                       &                     -atan(x)
               end if
            end if
            radialMomentTwoThirds=self%momentRadial2Previous
         else if (moment == 3) then
            ! 3ʳᵈ moment.
            if (x /= self%momentRadial3XPrevious) then
               self%momentRadial3XPrevious=x
               if (x < xSmall) then
                  ! Series solution for small x.
                  self%momentRadial3Previous=x**4*(1.0d0/4.0d0-x**2/6.0d0)
               else
                  ! General solution.
                  self%momentRadial3Previous=+0.5d0        &
                       &                     *(            &
                       &                       +x**2       &
                       &                       -log(       &
                       &                            +1.0d0 &
                       &                            +x**2  &
                       &                           )       &
                       &                      )
               end if
            end if
            radialMomentTwoThirds=self%momentRadial3Previous
         else
            radialMomentTwoThirds=0.0d0
            call Error_Report('unsupported moment'//{introspection:location})
         end if
      end if
      return
    end function radialMomentTwoThirds

  end function betaProfileDensityRadialMoment

  double precision function betaProfileDensitySquareIntegral(self,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Computes the integral of the square of the density in a $\beta$-profile mass distribution.
    !!}
    use :: Gamma_Functions         , only : Gamma_Function
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionBetaProfile ), intent(inout)           :: self
    double precision                              , intent(in   ), optional :: radiusMinimum          , radiusMaximum
    logical                                       , intent(  out), optional :: isInfinite
    logical                                                                 :: haveRadiusMinimum      , haveRadiusMaximum
    double precision                                                        :: radiusMinimum_         , radiusMaximum_         , &
         &                                                                     fractionalRadiusMinimum, fractionalRadiusMaximum

    if (present(isInfinite)) isInfinite=.false.
    betaProfileDensitySquareIntegral=0.0d0
    ! Determine effective radii.
    haveRadiusMinimum=present(radiusMinimum)
    haveRadiusMaximum=present(radiusMaximum)
    radiusMinimum_   =                0.0d0
    radiusMaximum_   =huge   (        0.0d0)
    if (haveRadiusMinimum) then
       radiusMinimum_   =    radiusMinimum
    end if
    if (haveRadiusMaximum) then
       radiusMaximum_   =    radiusMaximum
    end if
    if (self%truncateAtOuterRadius) then
       radiusMinimum_   =min(radiusMinimum_,self%outerRadius)
       radiusMaximum_   =min(radiusMaximum_,self%outerRadius)
       haveRadiusMinimum=.true.
       haveRadiusMaximum=.true.
    end if
    ! Determine if the special case solution can be used.
    if (self%betaIsTwoThirds) then
       ! Compute the integral for the case β=2/3.
       if (haveRadiusMinimum) then
          fractionalRadiusMinimum=+     radiusMinimum_ &
               &                  /self%coreRadius
          betaProfileDensitySquareIntegral=+betaProfileDensitySquareIntegral-4.0d0*Pi*(atan(fractionalRadiusMinimum)/2.0d0-fractionalRadiusMinimum/2.0d0/(1.0d0+fractionalRadiusMinimum**2))
       else
          betaProfileDensitySquareIntegral=+betaProfileDensitySquareIntegral+0.0d0
       end if
       if (haveRadiusMaximum) then
          fractionalRadiusMaximum=+     radiusMaximum_ &
               &                  /self%coreRadius
          betaProfileDensitySquareIntegral=betaProfileDensitySquareIntegral+4.0d0*Pi*(atan(fractionalRadiusMaximum)/2.0d0-fractionalRadiusMaximum/2.0d0/(1.0d0+fractionalRadiusMaximum**2))
       else
          betaProfileDensitySquareIntegral=betaProfileDensitySquareIntegral+Pi**2
       end if
    else
       ! Compute the integral for the general case.
       if (haveRadiusMinimum) then
          fractionalRadiusMinimum=+     radiusMinimum_ &
               &                  /self%coreRadius
          betaProfileDensitySquareIntegral=+betaProfileDensitySquareIntegral-4.0d0*Pi/3.0d0*fractionalRadiusMinimum**3*Hypergeometric_2F1([1.5d0,3.0d0*self%beta],[2.5d0],-fractionalRadiusMinimum**2)
       else
          betaProfileDensitySquareIntegral=+betaProfileDensitySquareIntegral+0.0d0
       end if
       if (haveRadiusMaximum) then
          fractionalRadiusMaximum=+     radiusMaximum_ &
               &                  /self%coreRadius
          betaProfileDensitySquareIntegral=+betaProfileDensitySquareIntegral+4.0d0*Pi/3.0d0*fractionalRadiusMaximum**3*Hypergeometric_2F1([1.5d0,3.0d0*self%beta],[2.5d0],-fractionalRadiusMaximum**2)
       else
          if (self%beta > 2.0d0/3.0d0) then
             betaProfileDensitySquareIntegral=+betaProfileDensitySquareIntegral+Pi**1.5d0*Gamma_Function(-1.5d0+3.0d0*self%beta)/Gamma_Function(3.0d0*self%beta)
          else
             if (present(isInfinite)) then
                isInfinite=.true.
                return
             else
                call Error_Report('integral diverges'//{introspection:location})
             end if
          end if
       end if
    end if
    ! Convert to dimensionful units.
    betaProfileDensitySquareIntegral=+betaProfileDensitySquareIntegral    &
         &                           *self%densityNormalization       **2 &
         &                           *self%coreRadius                 **3
    return
  end function betaProfileDensitySquareIntegral
  
  subroutine betaProfileDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionBetaProfile), intent(inout)           :: self
    type     (inputParameters            ), intent(inout)           :: descriptor
    logical                               , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                     )                          :: parameterLabel
    type     (inputParameters            )                          :: parameters
    !$GLC attributes unused :: includeFileModificationTimes
    
    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','betaProfile')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%beta
    call parameters%addParameter('beta'                ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%coreRadius
    call parameters%addParameter('coreRadius'          ,trim(adjustl(parameterLabel)))
    return
  end subroutine betaProfileDescriptor
