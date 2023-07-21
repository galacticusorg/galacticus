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

!+    Contributions to this file made by: Xiaolong Du, Andrew Benson.

  !!{
  Implements an exponentially truncated spherical mass distribution \cite{kazantzidis_2006}.
  !!}

  !![
  <massDistribution name="massDistributionSphericalTruncatedExponential">
   <description>
     Implements an exponentially truncated mass distribution \cite{kazantzidis_2006} in which the density is given by
     \begin{equation}
       \rho(r) = \rho^\prime(r_\mathrm{min}) \left\{ \begin{array}{ll} 1 &amp; \hbox{ if } r &lt; r_\mathrm{min}, \\ \rho^\prime(r_\mathrm{min} x^\kappa \exp\left(-\frac{x-1}{x_\mathrm{max}}\right) &amp; \hbox{otherwise,} \end{array} \right.
     \end{equation}
     where $x = r/r_\mathrm{min}$, $x_\mathrm{decay} = r_\mathrm{decay}/r_\mathrm{min}$, $\rho^\prime(r)$ is some other density
     profile, $r_\mathrm{min}=${\normalfont \ttfamily [radiusTruncateMinimum]}, $r_\mathrm{decay}=${\normalfont \ttfamily
     [radiusTruncateDecay]}, and
     \begin{equation}
     \kappa = \frac{r_\mathrm{min}}{r_\mathrm{decay}} + \frac{\mathrm{d}\log \rho^\prime}{\mathrm{d}\log r}(r_\mathrm{min})
     \end{equation}
     is chosen to ensure that the logarithmic gradient of the density profile is continuous across $r=r_\mathrm{min}$.
     </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalDecorator) :: massDistributionSphericalTruncatedExponential
     !!{
     Implementation of an exponentially-truncated spherical mass distribution.
     !!}
     private
     class           (massDistributionHeatingClass), pointer :: massDistributionHeating_ => null()
     double precision                                        :: radiusTruncateMinimum             , radiusTruncateDecay, &
          &                                                     massAtTruncation                  , massTotal_         , &
          &                                                     kappa
   contains
     final     ::                           sphericalTruncatedExponentialDestructor
     procedure :: density                => sphericalTruncatedExponentialDensity
     procedure :: densityGradientRadial  => sphericalTruncatedExponentialDensityGradientRadial
     procedure :: massTotal              => sphericalTruncatedExponentialMassTotal
     procedure :: massEnclosedBySphere   => sphericalTruncatedExponentialMassEnclosedBySphere
     procedure :: radiusEnclosingMass    => sphericalTruncatedExponentialRadiusEnclosingMass
  end type massDistributionSphericalTruncatedExponential

  interface massDistributionSphericalTruncatedExponential
     !!{
     Constructors for the {\normalfont \ttfamily sphericalTruncatedExponential} mass distribution class.
     !!}
     module procedure sphericalTruncatedExponentialConstructorParameters
     module procedure sphericalTruncatedExponentialConstructorInternal
  end interface massDistributionSphericalTruncatedExponential

contains

  function sphericalTruncatedExponentialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sphericalTruncatedExponential} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalTruncatedExponential)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (massDistributionClass                        ), pointer       :: massDistribution_
    type            (varying_string                               )                :: nonAnalyticSolver
    double precision                                                               :: radiusTruncateMinimum, radiusTruncateDecay         
    type            (varying_string                               )                :: componentType        , massType

    !![
    <inputParameter>
      <name>radiusTruncateMinimum</name>
      <source>parameters</source>
      <description>The minimum radius to begin truncating the density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusTruncateDecay</name>
      <source>parameters</source>
      <description>The exponential decay scale for truncating the density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
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
       self=massDistributionSphericalTruncatedExponential(radiusTruncateMinimum,radiusTruncateDecay,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function sphericalTruncatedExponentialConstructorParameters
  
  function sphericalTruncatedExponentialConstructorInternal(radiusTruncateMinimum,radiusTruncateDecay,nonAnalyticSolver,massDistribution_,componentType,massType) result(self)
    !!{
    Constructor for ``sphericalTruncatedExponential'' mass distribution class.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    type            (massDistributionSphericalTruncatedExponential)                          :: self
    class           (massDistributionSpherical                    ), intent(in   ), target   :: massDistribution_
    type            (enumerationNonAnalyticSolversType            ), intent(in   )           :: nonAnalyticSolver
    double precision                                               , intent(in   )           :: radiusTruncateMinimum            , radiusTruncateDecay
    type            (enumerationComponentTypeType                 ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType                      ), intent(in   ), optional :: massType
    double precision                                               , parameter               :: fractionRadialDecayMaximum=50.0d0
    type            (coordinateSpherical                          )                          :: coordinatesTruncateMinimum
    !![
    <constructorAssign variables="radiusTruncateMinimum, radiusTruncateDecay, nonAnalyticSolver, *massDistribution_, componentType, massType"/>
    !!]

    coordinatesTruncateMinimum=[self%radiusTruncateMinimum,0.0d0,0.0d0]
    self%massAtTruncation     =+self%massDistribution_%massEnclosedBySphere(                                               &
         &                                                                               +self%radiusTruncateMinimum       &
         &                                                                  )
    self%massTotal_           =+self                  %massEnclosedBySphere(                                               &
         &                                                                               +self%radiusTruncateMinimum       &
         &                                                                               +     fractionRadialDecayMaximum  &
         &                                                                               *self%radiusTruncateDecay         &
         &                                                                  )
    self%kappa                =+self%radiusTruncateMinimum                                                                 &
         &                     /self%radiusTruncateDecay                                                                   &
         &                     +self%massDistribution_%densityGradientRadial(                                              &
         &                                                                               coordinatesTruncateMinimum      , &
         &                                                                   logarithmic=.true.&
         &                                                                  )
    self%dimensionless        = self%massDistribution_%isDimensionless      (                                              &
         &                                                                  )
    return
  end function sphericalTruncatedExponentialConstructorInternal

  subroutine sphericalTruncatedExponentialDestructor(self)
    !!{
    Destructor for the abstract {\normalfont \ttfamily massDistributionSphericalTruncatedExponential} class.
    !!}
    implicit none
    type(massDistributionSphericalTruncatedExponential), intent(inout) :: self
    
    !![
    <objectDestcutor name="self%massDistribution_"/>
    !!]
    return
  end subroutine sphericalTruncatedExponentialDestructor

  double precision function sphericalTruncatedExponentialDensity(self,coordinates,componentType,massType) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an exponentially-truncated spherical mass distribution.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class(massDistributionSphericalTruncatedExponential), intent(inout)           :: self
    class(coordinate                                   ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType                 ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType                      ), intent(in   ), optional :: massType
    type (coordinateSpherical                          )                          :: coordinatesTruncateMinimum
    
    if (coordinates%rSpherical() <= self%radiusTruncateMinimum) then
       density                   =+self%massDistribution_%density(coordinates,componentType,massType)
    else
       coordinatesTruncateMinimum=[self%radiusTruncateMinimum,0.0d0,0.0d0]
       density                   =+self%massDistribution_%density(coordinatesTruncateMinimum,componentType,massType) &
            &                     *(                                                                                 &
            &                       +coordinates%rSpherical           ()                                             &
            &                       /self       %radiusTruncateMinimum                                               &
            &                      )**self%kappa                                                                     &
            &                     *exp(                                                                              &
            &                          -(                                                                            &
            &                            +coordinates%rSpherical           ()                                        &
            &                            -self       %radiusTruncateMinimum                                          &
            &                           )                                                                            &
            &                          /  self       %radiusTruncateDecay                                            &
            &                        )
    end if
    return
  end function sphericalTruncatedExponentialDensity

  double precision function sphericalTruncatedExponentialDensityGradientRadial(self,coordinates,logarithmic,componentType,massType) result(densityGradient)
    !!{
    Return the density gradient at the specified {\normalfont \ttfamily coordinates} in an exponentially-truncated spherical mass distribution.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class  (massDistributionSphericalTruncatedExponential), intent(inout), target   :: self
    class  (coordinate                                   ), intent(in   )           :: coordinates
    logical                                               , intent(in   ), optional :: logarithmic
    type   (enumerationComponentTypeType                 ), intent(in   ), optional :: componentType
    type   (enumerationMassTypeType                      ), intent(in   ), optional :: massType
    type   (coordinateSpherical                          )                          :: coordinatesTruncateMinimum
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

     if (coordinates%rSpherical() <= self%radiusTruncateMinimum) then
        densityGradient=+self%massDistribution_%densityGradientRadial(coordinates,logarithmic,componentType,massType)
     else
        coordinatesTruncateMinimum=[self%radiusTruncateMinimum,0.0d0,0.0d0]
        densityGradient=+self%massDistribution_%densityGradientRadial(coordinates               ,componentType=componentType,massType=massType,logarithmic=.false.) &
             &          +self%massDistribution_%density              (coordinatesTruncateMinimum,componentType=componentType,massType=massType                    ) &
             &          *(                                                                                                                                          &
             &             +coordinates%rSpherical           ()                                                                                                     &
             &             /self       %radiusTruncateMinimum                                                                                                       &
             &            )**self%kappa                                                                                                                             &
             &           *exp(                                                                                                                                      &
             &                -(                                                                                                                                    &
             &                  +coordinates%rSpherical           ()                                                                                                &
             &                  -self       %radiusTruncateMinimum                                                                                                  &
             &                 )                                                                                                                                    &
             &                /  self       %radiusTruncateDecay                                                                                                    &
             &               )                                                                                                                                      &
             &           *(                                                                                                                                         &
             &             +self       %kappa                                                                                                                       &
             &             -coordinates%rSpherical         ()                                                                                                       &
             &             /self       %radiusTruncateDecay                                                                                                         &
             &            )                                                                                                                                         &
             &           /  coordinates%rSpherical         ()
        if (logarithmic_) densityGradient=+            densityGradient                                     &
             &                            *coordinates%rSpherical     (                                  ) &
             &                            /self       %density        (coordinates,componentType,massType)
     end if
    return
  end function sphericalTruncatedExponentialDensityGradientRadial

  double precision function sphericalTruncatedExponentialMassTotal(self,componentType,massType) result(mass)
    !!{
    Return the total mass in a trunacted mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalTruncatedExponential), intent(inout)           :: self
    type (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType           ), intent(in   ), optional :: massType

    if (self%matches(componentType,massType)) then
       mass=self%massTotal_
    else
       mass=0.0d0
    end if
    return
  end function sphericalTruncatedExponentialMassTotal
  
  double precision function sphericalTruncatedExponentialMassEnclosedBySphere(self,radius,componentType,massType) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for truncatedExponential mass distributions.
    !!}
    use :: Coordinates    , only : coordinateSpherical                   , assignment(=)
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Unnormalized
    implicit none
    class           (massDistributionSphericalTruncatedExponential), intent(inout), target   :: self
    double precision                                               , intent(in   )           :: radius
    type            (enumerationComponentTypeType                 ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType                      ), intent(in   ), optional :: massType
    type            (coordinateSpherical                          )                          :: coordinatesTruncateMinimum
    
    if (radius <= self%radiusTruncateMinimum) then
       mass                      =+self%massDistribution_%massEnclosedBySphere           (     radius               ,componentType,massType)
    else
       coordinatesTruncateMinimum=[self%radiusTruncateMinimum,0.0d0,0.0d0]
       mass                      =+self%massDistribution_%massEnclosedBySphere           (self%radiusTruncateMinimum,componentType,massType)                  &
            &                     +4.0d0                                                                                                                      &
            &                     *Pi                                                                                                                         &
            &                     *self%massDistribution_%density(coordinatesTruncateMinimum)                                                                 &
            &                     *self%radiusTruncateDecay**3                                                                                                &
            &                     *                                     exp(                 self%radiusTruncateMinimum/self%radiusTruncateDecay)             &
            &                     /                                        (                 self%radiusTruncateMinimum/self%radiusTruncateDecay)**self%kappa &
            &                     *(                                                                                                                          &
            &                       +Gamma_Function_Incomplete_Unnormalized(3.0d0+self%kappa,self%radiusTruncateMinimum/self%radiusTruncateDecay)             &
            &                       +Gamma_Function_Incomplete_Unnormalized(3.0d0+self%kappa,     radius               /self%radiusTruncateDecay)             &
            &                     )
    end if
    return
  end function sphericalTruncatedExponentialMassEnclosedBySphere

  double precision function sphericalTruncatedExponentialRadiusEnclosingMass(self,mass,massFractional,componentType,massType) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for truncatedExponential spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalTruncatedExponential), intent(inout), target   :: self
    double precision                                    , intent(in   ), optional :: mass         , massFractional
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType
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
       radius=self%massDistribution_%radiusEnclosingMass           (mass=mass_,componentType=componentType,massType=massType)
    else
       radius=self                  %radiusEnclosingMassNonAnalytic(mass=mass_,componentType=componentType,massType=massType)
    end if    
    return
  end function sphericalTruncatedExponentialRadiusEnclosingMass
