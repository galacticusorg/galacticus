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
  Implementation of mass distribution for the \cite{patej_simple_2015} model of the circumgalactic medium.
  !!}

  !![
  <massDistribution name="massDistributionPatejLoeb2015">
   <description>A mass distribution for the \cite{patej_simple_2015} model of the circumgalactic medium.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionPatejLoeb2015
     !!{
     The \cite{patej_simple_2015} model of the circumgalactic medium.
     !!}
     class           (massDistributionClass), pointer :: massDistribution_     => null()
     logical                                          :: truncateAtOuterRadius
     double precision                                 :: radiusShock                    , radiusOuter, &
          &                                              densityNormalization           , gamma      , &
          &                                              mass
   contains
     !![
     <methods>
       <method method="radiusDarkMatter"      description="Return the corresponding radius in the dark matter profile."     />
       <method method="coordinatesDarkMatter" description="Return the corresponding coordinates in the dark matter profile."/>
       </methods>
     !!]
     final     ::                          patejLoeb2015Destructor
     procedure :: density               => patejLoeb2015Density
     procedure :: densityGradientRadial => patejLoeb2015DensityGradientRadial
     procedure :: densityRadialMoment   => patejLoeb2015DensityRadialMoment
     procedure :: massEnclosedBySphere  => patejLoeb2015MassEnclosedBySphere
     procedure :: potentialIsAnalytic   => patejLoeb2015PotentialIsAnalytic
     procedure :: potential             => patejLoeb2015Potential
     procedure :: radiusDarkMatter      => patejLoeb2015RadiusDarkMatter
     procedure :: coordinatesDarkMatter => patejLoeb2015CoordinatesDarkMatter
  end type massDistributionPatejLoeb2015

  interface massDistributionPatejLoeb2015
     !!{
     Constructors for the {\normalfont \ttfamily patejLoeb2015} mass distribution class.
     !!}
     module procedure patejLoeb2015ConstructorParameters
     module procedure patejLoeb2015ConstructorInternal
  end interface massDistributionPatejLoeb2015

contains

  function patejLoeb2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily patejLoeb2015} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionPatejLoeb2015)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: radiusShock          , radiusOuter, &
         &                                                            mass                 , gamma      , &
         &                                                            densityNormalization
    class           (massDistributionClass        ), pointer       :: massDistribution_
    logical                                                        :: truncateAtOuterRadius
    type            (varying_string               )                :: componentType
    type            (varying_string               )                :: massType

    !![
    <inputParameter>
      <name>gamma</name>
      <defaultValue>1.15d0</defaultValue>
      <description>The parameter $\Gamma$ in the \cite{patej_simple_2015} mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The density normalization of the \cite{patej_simple_2015} mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The mass of the \cite{patej_simple_2015} mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusOuter</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The outer radius of the \cite{patej_simple_2015} mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>truncateAtOuterRadius</name>
      <defaultValue>.false.</defaultValue>
      <description>If true then the \cite{patej_simple_2015} mass distribution is truncated beyond the outer radius.</description>
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
    <objectBuilder class="massDistribution"  name="massDistribution_"  source="parameters"/>
    <conditionalCall>
     <call>self=massDistributionPatejLoeb2015(radiusShock=radiusShock,gamma=gamma,massDistribution_=massDistribution_,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="densityNormalization"  value="densityNormalization"  parameterPresent="parameters"/>
     <argument name="mass"                  value="mass"                  parameterPresent="parameters"/>
     <argument name="radiusOuter"           value="radiusOuter"           parameterPresent="parameters"/>
     <argument name="truncateAtOuterRadius" value="truncateAtOuterRadius" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function patejLoeb2015ConstructorParameters

  function patejLoeb2015ConstructorInternal(gamma,massDistribution_,densityNormalization,mass,radiusOuter,radiusShock,truncateAtOuterRadius,componentType,massType) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily patejLoeb2015} mass distribution class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (massDistributionPatejLoeb2015)                         :: self
    double precision                              , intent(in   )           :: gamma                , radiusShock
    double precision                              , intent(in   ), optional :: densityNormalization , mass       , &
         &                                                                     radiusOuter
    logical                                       , intent(in   ), optional :: truncateAtOuterRadius
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    class           (massDistributionClass       ), intent(in   ), target   :: massDistribution_
    !![
    <constructorAssign variables="gamma, densityNormalization, radiusShock, radiusOuter, truncateAtOuterRadius, componentType, massType, *massDistribution_"/>
    !!]

    ! This mass distribution is never dimensionless.
    self%dimensionless=.false.
    ! Determine density normalization.
    if      (                                   &
         &   present(densityNormalization)      &
         &  ) then
       self%densityNormalization=densityNormalization
    else if (                                   &
         &   present(mass                )      &
         &  ) then
       self%densityNormalization=mass/self%massDistribution_%massEnclosedBySphere(self%radiusDarkMatter(self%radiusShock))
    else
       call Error_Report('no way to determine density normalization'//{introspection:location})
    end if
    ! Check for truncation.
    if (present(truncateAtOuterRadius)) then
       self%truncateAtOuterRadius=truncateAtOuterRadius
    else
       self%truncateAtOuterRadius=.false.
    end if
    if (self%truncateAtOuterRadius) then
       if (.not.present(radiusOuter)) call Error_Report('can not truncate profile without an outer radius'//{introspection:location})
       self%radiusOuter=radiusOuter
    end if
    return
  end function patejLoeb2015ConstructorInternal

  subroutine patejLoeb2015Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily patejLoeb2015} mass distribution class.
    !!}
    type(massDistributionPatejLoeb2015), intent(inout) :: self
    implicit none

    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine patejLoeb2015Destructor
  
  double precision function patejLoeb2015RadiusDarkMatter(self,radius) result(radiusDarkMatter)
    !!{
    Return the corresponding radius in the dark matter distribution the specified {\normalfont \ttfamily radius} in a \cite{patej_simple_2015} mass distribution.
    !!}
    implicit none
    class           (massDistributionPatejLoeb2015), intent(inout) :: self
    double precision                               , intent(in   ) :: radius


    radiusDarkMatter=+  self%radiusShock &
         &           *(                  &
         &             +     radius      &
         &             /self%radiusShock &
         &           )**self%gamma
    return
  end function patejLoeb2015RadiusDarkMatter

  function patejLoeb2015CoordinatesDarkMatter(self,coordinates) result(coordinatesDarkMatter)
    !!{
    Return the corresponding coordinates in the dark matter distribution the specified {\normalfont \ttfamily radius} in a \cite{patej_simple_2015} mass distribution.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    type (coordinateSpherical          )                :: coordinatesDarkMatter
    class(massDistributionPatejLoeb2015), intent(inout) :: self
    class(coordinate                   ), intent(in   ) :: coordinates

    coordinatesDarkMatter=coordinates                                                      &
         &                *(coordinates%rSpherical()/self%radiusShock)**(self%gamma-1.0d0)
    return
  end function patejLoeb2015CoordinatesDarkMatter

  double precision function patejLoeb2015Density(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a \cite{patej_simple_2015} mass distribution.
    !!}
    implicit none
    class(massDistributionPatejLoeb2015), intent(inout) :: self
    class(coordinate                   ), intent(in   ) :: coordinates

    density=0.0d0
    if (self%truncateAtOuterRadius .and. coordinates%rSpherical() > self%radiusOuter) return
    ! Evaluate density using equation 12 of Patej & Loeb (2015).
    density=+self%gamma                                                              &
         &  *self%densityNormalization                                               &
         &  *(                                                                       &
         &    +coordinates%rSpherical ()                                             &
         &    /self       %radiusShock                                               &
         &   )**(3.0d0*self%gamma-3.0d0)                                             &
         &  *self%massDistribution_%density(self%coordinatesDarkMatter(coordinates))
    return
  end function patejLoeb2015Density

  double precision function patejLoeb2015DensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a \cite{patej_simple_2015} mass distribution.
    !!}
    implicit none
    class  (massDistributionPatejLoeb2015), intent(inout), target   :: self
    class  (coordinate                   ), intent(in   )           :: coordinates
    logical                               , intent(in   ), optional :: logarithmic
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    densityGradientRadial=0.0d0
    if (self%truncateAtOuterRadius .and. coordinates%rSpherical() > self%radiusOuter) return
    densityGradientRadial=+3.0d0                                                                                             &
         &                *(self%gamma-1.0d0)                                                                                &
         &                + self%gamma                                                                                       &
         &                * self%massDistribution_%densityGradientRadial(                                                    &
         &                                                                          self%coordinatesDarkMatter(coordinates), &
         &                                                              logarithmic=.true.                                   &
         &                                                            )
    if (.not.logarithmic_) densityGradientRadial=densityGradientRadial*self%density(coordinates)/coordinates%rSpherical()
    return
  end function patejLoeb2015DensityGradientRadial

  double precision function patejLoeb2015MassEnclosedBySphere(self,radius) result(massEnclosedBySphere)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a \cite{patej_simple_2015} mass distribution.
    !!}
    implicit none
    class           (massDistributionPatejLoeb2015), intent(inout), target :: self
    double precision                               , intent(in   )         :: radius
    double precision                                                       :: radius_

    ! Compute the enclosed mass (eqn. 4 of Patej & Loeb 2015).
    radius_=radius
    if (self%truncateAtOuterRadius) radius_=min(radius_,self%radiusOuter)
    massEnclosedBySphere=+self                  %densityNormalization                                 &
         &               *self%massDistribution_%massEnclosedBySphere(self%radiusDarkMatter(radius_))
    return
  end function patejLoeb2015MassEnclosedBySphere

  logical function patejLoeb2015PotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionPatejLoeb2015), intent(inout) :: self

    isAnalytic=.true.
    return
  end function patejLoeb2015PotentialIsAnalytic

  double precision function patejLoeb2015Potential(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a \cite{patej_simple_2015} mass distribution. The
    potential is given by
    \begin{equation}
      \phi(r) = - \int^\infty_r \mathrm{d}r \frac{\mathrm{G}M(r)}{r^2}.
    \end{equation}    
    Given that $M(r) = f M^\prime(r^\prime)$ where $M^\prime(r^\prime)$ is the mass profile of the dark matter distribution,
    $r^\prime = r_\mathrm{s} (r/r_\mathrm{s})^\Gamma$ (and, therefore, $r = r_\mathrm{s} (r^\prime/r_\mathrm{s})^{1/\Gamma}$), and
    $f$ is the normalization factor, we can write this as:
    \begin{eqnarray}
          \phi(r) &=& - \int^\infty_{r^\prime(r)} \mathrm{d}r^\prime  f \frac{\mathrm{d}r}{\mathrm{d}r^\prime} (r^\prime/r_\mathrm{s})^{-2/\Gamma} \frac{\mathrm{G}M^\prime(r^\prime)}{r_\mathrm{s}^2} \nonumber \\
          &=&   - \int^\infty_{r^\prime(r)} \mathrm{d}r^\prime  \frac{f}{\Gamma} (r^\prime/r_\mathrm{s})^{-1/\Gamma-1} \frac{\mathrm{G}M^\prime(r^\prime)}{r_\mathrm{s}^2}
    \end{eqnarray}
    which we can then integrate by parts to get:
    \begin{equation}
      \phi(r) = + f \left(\frac{r^\prime}{r_\mathrm{s}}\right)^{-1/\Gamma} \frac{\mathrm{G}M^\prime(r^\prime)}{r_\mathrm{s}} - 4 \pi f \frac{\mathrm{G} r_\mathrm{s}^{1/\Gamma} \mathcal{R}(r^\prime,\infty,2-1/\Gamma)}{r_\mathrm{s}},
    \end{equation}
    where $\mathcal{R}(a,b,m)$ is the $m^\mathrm{th}$ radial density moment between $r^\prime=a$ and $r^\prime=b$ of the dark
    matter distribution, and we have assumed that $(r^\prime/r_\mathrm{s})^{-1/\Gamma} M^\prime(r^\prime) \rightarrow 0$ as
    $r^\prime \rightarrow \infty$.    
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    implicit none
    class           (massDistributionPatejLoeb2015    ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    double precision                                                             :: radiusOuter          , radiusDarkMatterInner, &
         &                                                                          radiusDarkMatterOuter, potentialTermInner   , &
         &                                                                          potentialTermOuter   , potentialTermMass
    logical                                                                      :: computePotential
    
    potential=0.0d0
    if (present(status)) status=structureErrorCodeSuccess
    ! Note that we should check that the condition -1/Γ+3+α < 0, where α is the logarithmic slope of the dark matter density
    ! profile as r → ∞, holds to ensure that our assumption that the r=∞ term in the first term in the integration by parts is
    ! zero is a valid assumption. Currently we do not perform this check.
    if (self%truncateAtOuterRadius) then
       ! Add the potential due to a point mass outside of the outer radius.
       radiusOuter     =max(                           & 
            &               coordinates%rSpherical (), &
            &               self       %radiusOuter    &
            &              )
       potential       =-     gravitationalConstant_internal              &
            &           *self%massEnclosedBySphere          (radiusOuter) &
            &           /     radiusOuter
       computePotential=radiusOuter < self%radiusOuter
    else
       radiusOuter     =huge(0.0d0)
       computePotential=.true.
    end if
    if (.not.computePotential) return
    ! Add on the contribution from radii inside the mass distribution.
    radiusDarkMatterInner=self%radiusDarkMatter(coordinates%rSpherical())
    if (self%truncateAtOuterRadius) then
       radiusDarkMatterOuter=+self%radiusDarkMatter(radiusOuter)
       potentialTermMass    =-4.0d0                                                                                                          &
            &                *Pi                                                                                                             &
            &                *self                  %radiusShock**(1.0d0/self%gamma)                                                         &
            &                *self%massDistribution_%densityRadialMoment(2.0d0-1.0d0/self%gamma,radiusDarkMatterInner,radiusDarkMatterOuter)
       potentialTermOuter   =+(                                                                                                              &
            &                  +     radiusDarkMatterOuter                                                                                   &
            &                  /self%radiusShock                                                                                             &
            &                )**(-1.0d0/self%gamma)                                                                                          &
            &                *self%massDistribution_%massEnclosedBySphere(radiusDarkMatterOuter)
    else
       ! There is no outer truncation, so the outer radius is at infinity.
       potentialTermMass    =-4.0d0                                                                                                          &
            &                *Pi                                                                                                             &
            &                *self                  %radiusShock**(1.0d0/self%gamma)                                                         &
            &                *self%massDistribution_%densityRadialMoment(2.0d0-1.0d0/self%gamma,radiusDarkMatterInner                      )
       potentialTermOuter   =+0.0d0
    end if
    potentialTermInner      =+(                                                                                                               &
         &                     +     radiusDarkMatterInner                                                                                    &
         &                     /self%radiusShock                                                                                              &
         &                    )**(-1.0d0/self%gamma)                                                                                          &
         &                   *self%massDistribution_%massEnclosedBySphere(radiusDarkMatterInner)
    potential               =+     potential                                                                                                  &
         &                   +     gravitationalConstant_internal                                                                             &
         &                   *self%densityNormalization                                                                                       &
         &                   /self%radiusShock                                                                                                &
         &                   *(                                                                                                               &
         &                     +potentialTermOuter                                                                                            &
         &                     -potentialTermInner                                                                                            &
         &                     +potentialTermMass                                                                                             &
         &                    )
    return    
  end function patejLoeb2015Potential

  double precision function patejLoeb2015DensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{
    Computes radial moments of the density in a \cite{patej_simple_2015} mass distribution. For this profile we have:
    \begin{equation}
    \rho_\mathrm{g}(r) = f \Gamma (r/s)^{3 \Gamma - 3} \rho_\mathrm{DM}(s[r/s]^\Gamma).
    \end{equation}
    Defining $R=s[r/s]^\Gamma$, such that $r/s = (R/s)^{1/\Gamma}$, and $\mathrm{d}r = \Gamma^{-1} (R/s)^{1/\Gamma-1} \mathrm{d}R$, then
    \begin{equation}
    \int r^m \rho_\mathrm{g}(r) \mathrm{d}r = f s^{(m-2)(\Gamma-1)/\Gamma} \int R^{(2\Gamma-2+m)/\Gamma} \rho_\mathrm{DM}(R) \mathrm{d}R,
    \end{equation}
    or
    \begin{equation}
    \mathcal{R}_\mathrm{g}(r;m) = f s^{(m-2)(\Gamma-1)/\Gamma} \mathcal{R}_\mathrm{DM}(R;(2\Gamma-2+m)/\Gamma),
    \end{equation}
    where $\mathcal{R}(r;m)$ is the $m^\mathrm{th}$ radial moment of the density profile.
    !!}
    implicit none
    class           (massDistributionPatejLoeb2015), intent(inout)           :: self
    double precision                               , intent(in   )           :: moment
    double precision                               , intent(in   ), optional :: radiusMinimum  , radiusMaximum
    logical                                        , intent(  out), optional :: isInfinite
    double precision                                                         :: radiusMinimum_ , radiusMaximum_, &
         &                                                                      momentEffective

    if (present(radiusMinimum)) then
       radiusMinimum_=radiusMinimum
    else
       radiusMinimum_=0.0d0
    end if
    if (present(radiusMaximum)) then
       radiusMaximum_=radiusMaximum
    else
       radiusMaximum_=huge(0.0d0)
    end if
    if (self%truncateAtOuterRadius) then
       radiusMinimum_=min(radiusMinimum_,self%radiusOuter)
       radiusMaximum_=min(radiusMaximum_,self%radiusOuter)
    end if
    momentEffective=(2.0d0*self%gamma-2.0d0+moment)/self%gamma
    if (radiusMaximum_ < huge(0.0d0)) then
       densityRadialMoment=+self                  %densityNormalization                                                                                                                        &
            &              *self                  %radiusShock                                                                                                      **(moment-momentEffective) &
            &              *self%massDistribution_%densityRadialMoment (momentEffective,self%radiusDarkMatter(radiusMinimum_),self%radiusDarkMatter(radiusMaximum_))
    else
       densityRadialMoment=+self                  %densityNormalization                                                                                                                        &
            &              *self                  %radiusShock                                                                                                      **(moment-momentEffective) &
            &              *self%massDistribution_%densityRadialMoment (momentEffective,self%radiusDarkMatter(radiusMinimum_))
       end if
    return
  end function patejLoeb2015DensityRadialMoment
