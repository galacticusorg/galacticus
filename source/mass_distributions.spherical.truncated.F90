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

  !!{
  Implements a truncated spherical mass distribution.
  !!}

  !![
  <massDistribution name="massDistributionSphericalTruncated">
   <description>
     Implements a mass distribution in which the density is given by
     \begin{equation}
       \rho(r) = \rho^\prime(r) \left\{ \begin{array}{ll} 1 &amp; \hbox{ if } r &lt; r_\mathrm{min}, \\ 0 &amp; \hbox{ if } r &gt; r_\mathrm{max}, \\ 1-3 x^2 + 2x^3 &amp; \hbox{otherwise,} \end{array} \right.
     \end{equation}
     where
     \begin{equation}
       x=\frac{r-r_\mathrm{min}}{r_\mathrm{max}-r_\mathrm{min}},
     \end{equation}
     $\rho^\prime(r)$ is some other density profile, $r_\mathrm{min}=${\normalfont \ttfamily [radiusTruncateMinimum]}, and  $r_\mathrm{max}=${\normalfont \ttfamily [radiusTruncateMaximum]}.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalDecorator) :: massDistributionSphericalTruncated
     !!{
     Implementation of a truncated spherical mass distribution.
     !!}
     private
     double precision :: radiusTruncateMinimum, radiusTruncateMaximum, &
          &              massAtTruncation     , massTotal_
   contains
     !![
     <methods>
       <method method="truncationFunction" description="Compute the truncation fraction (and related quantities) from the radius."/>
     </methods>
     !!]
     final     ::                           sphericalTruncatedDestructor
     procedure :: density                => sphericalTruncatedDensity
     procedure :: densityGradientRadial  => sphericalTruncatedDensityGradientRadial
     procedure :: massTotal              => sphericalTruncatedMassTotal
     procedure :: massEnclosedBySphere   => sphericalTruncatedMassEnclosedBySphere
     procedure :: radiusEnclosingMass    => sphericalTruncatedRadiusEnclosingMass
     procedure :: truncationFunction     => sphericalTruncatedTruncationFunction
  end type massDistributionSphericalTruncated

  interface massDistributionSphericalTruncated
     !!{
     Constructors for the \refClass{massDistributionSphericalTruncated} mass distribution class.
     !!}
     module procedure sphericalTruncatedConstructorParameters
     module procedure sphericalTruncatedConstructorInternal
  end interface massDistributionSphericalTruncated

contains

  function sphericalTruncatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalTruncated} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalTruncated)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (massDistributionClass             ), pointer       :: massDistribution_
    type            (varying_string                    )                :: nonAnalyticSolver
    double precision                                                    :: radiusTruncateMinimum, radiusTruncateMaximum
    type            (varying_string                    )                :: componentType        , massType

    !![
    <inputParameter>
      <name>radiusTruncateMinimum</name>
      <source>parameters</source>
      <description>The minimum radius to begin truncating the density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusTruncateMaximum</name>
      <source>parameters</source>
      <description>The maximum radius to finish truncating the density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available.</description>
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
    <objectBuilder class="massDistribution" name="massDistribution_" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       self=massDistributionSphericalTruncated(radiusTruncateMinimum,radiusTruncateMaximum,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function sphericalTruncatedConstructorParameters
  
  function sphericalTruncatedConstructorInternal(radiusTruncateMinimum,radiusTruncateMaximum,nonAnalyticSolver,massDistribution_,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalTruncated} mass distribution class.
    !!}
    implicit none
    type            (massDistributionSphericalTruncated)                          :: self
    class           (massDistributionSpherical         ), intent(in   ), target   :: massDistribution_
    type            (enumerationNonAnalyticSolversType ), intent(in   )           :: nonAnalyticSolver
    double precision                                    , intent(in   )           :: radiusTruncateMinimum, radiusTruncateMaximum
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="radiusTruncateMinimum, radiusTruncateMaximum, nonAnalyticSolver, *massDistribution_, componentType, massType"/>
    !!]

    self%massAtTruncation=self%massDistribution_%massEnclosedBySphere(radiusTruncateMinimum)
    self%massTotal_      =self                  %massEnclosedBySphere(radiusTruncateMaximum)
    self%dimensionless   =self%massDistribution_%isDimensionless     (                     )
    return
  end function sphericalTruncatedConstructorInternal

  subroutine sphericalTruncatedDestructor(self)
    !!{
    Destructor for the abstract \refClass{massDistributionSphericalTruncated} class.
    !!}
    implicit none
    type(massDistributionSphericalTruncated), intent(inout) :: self
    
    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine sphericalTruncatedDestructor

  subroutine sphericalTruncatedTruncationFunction(self,radius,x,multiplier,multiplierGradient)
    !!{
    Return the scaled truncation radial coordinate, and the truncation multiplier.
    !!}
    implicit none
    class           (massDistributionSphericalTruncated), intent(inout)           :: self
    double precision                                    , intent(in   )           :: radius
    double precision                                    , intent(  out), optional :: x                 , multiplier, &
         &                                                                           multiplierGradient
    double precision                                                              :: x_

    if      (radius <= self%radiusTruncateMinimum) then
       if (present(x                 )) x                 =+0.0d0
       if (present(multiplier        )) multiplier        =+1.0d0
       if (present(multiplierGradient)) multiplierGradient=+0.0d0
    else if (radius >= self%radiusTruncateMaximum) then
       if (present(x                 )) x                 =+1.0d0
       if (present(multiplier        )) multiplier        =+0.0d0
       if (present(multiplierGradient)) multiplierGradient=+0.0d0
    else
       x_                                                 =+(+     radius               -self%radiusTruncateMinimum) &
            &                                              /(+self%radiusTruncateMaximum-self%radiusTruncateMinimum)
       if (present(x                 )) x                 =+x_
       if (present(multiplier        )) multiplier        =  +1.0d0                                                  &
            &                                                -3.0d0*x_**2                                            &
            &                                                +2.0d0*x_**3
       if (present(multiplierGradient)) multiplierGradient=+(                                                        &
            &                                                -6.0d0*x_                                               &
            &                                                +6.0d0*x_**2                                            &
            &                                               )                                                        &
            &                                              /(+self%radiusTruncateMaximum-self%radiusTruncateMinimum)
    end if
    return
  end subroutine sphericalTruncatedTruncationFunction

  double precision function sphericalTruncatedDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalTruncated), intent(inout) :: self
    class           (coordinate                        ), intent(in   ) :: coordinates
    double precision                                                    :: multiplier

    call self%truncationFunction(radius=coordinates%rSpherical(),multiplier=multiplier)
    density=+self%massDistribution_%density(coordinates) &
         &  *                       multiplier
    return
  end function sphericalTruncatedDensity

  double precision function sphericalTruncatedDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a truncated spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalTruncated), intent(inout), target   :: self
    class           (coordinate                        ), intent(in   )           :: coordinates
    logical                                             , intent(in   ), optional :: logarithmic
    double precision                                                              :: multiplier , multiplierGradient
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    call self%truncationFunction(radius=coordinates%rSpherical(),multiplier=multiplier,multiplierGradient=multiplierGradient)
    if (multiplier > 0.0d0) then
       densityGradient=+self%massDistribution_%densityGradientRadial(coordinates,logarithmic=.false.) &
            &          *                       multiplier                                             &
            &          +self%massDistribution_%density              (coordinates                    ) &
            &          *                       multiplierGradient
       if (logarithmic_) densityGradient=+            densityGradient              &
            &                            *coordinates%rSpherical     (           ) &
            &                            /self       %density        (coordinates)
    else
       densityGradient=+0.0d0
    end if
    return
  end function sphericalTruncatedDensityGradientRadial

  double precision function sphericalTruncatedMassTotal(self) result(mass)
    !!{
    Return the total mass in a truncated mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalTruncated), intent(inout) :: self

    mass=self%massTotal_
    return
  end function sphericalTruncatedMassTotal
  
  double precision function sphericalTruncatedMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for truncated mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalTruncated), intent(inout), target :: self
    double precision                                    , intent(in   )         :: radius

    if (radius <= self%radiusTruncateMinimum) then
       mass=self%massDistribution_%massEnclosedBySphere           (radius)
    else
       mass=self                  %massEnclosedBySphereNonAnalytic(radius)
    end if
    return
  end function sphericalTruncatedMassEnclosedBySphere
  
  double precision function sphericalTruncatedRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for truncated spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalTruncated), intent(inout), target   :: self
    double precision                                    , intent(in   ), optional :: mass , massFractional
    double precision                                                              :: mass_
    
    if (present(mass)) then
       mass_  =+     mass
    else if (present(massFractional)) then
       mass_  =+     massFractional   &
            &  *self%massTotal     ()
    else
       mass_  =+0.0d0
       call Error_Report('either `mass` or `massFractional` must be provided'//{introspection:location})
    end if
    if (mass_ <= self%massAtTruncation) then
       radius=self%massDistribution_%radiusEnclosingMass           (mass=mass_)
    else
       radius=self                  %radiusEnclosingMassNonAnalytic(mass=mass_)
    end if    
    return
  end function sphericalTruncatedRadiusEnclosingMass
