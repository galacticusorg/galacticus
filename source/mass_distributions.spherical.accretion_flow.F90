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
  Implements a accretionFlow spherical mass distribution.
  !!}

  !![
  <massDistribution name="massDistributionSphericalAccretionFlow">
   <description>
      An implementation of a mass distribution which includes the accretion flow surrounding a halo. The density
      profile is modeled as
      \begin{equation}
      \rho(r) = f_\mathrm{trans}(r) \rho_\mathrm{halo}(r) + \rho_\mathrm{accretion}(r),
      \end{equation}
      where $\rho_\mathrm{halo}(r)$ is the halo mass distribution, $\rho_\mathrm{accretion}(r)$ is the accretion flow mass distribution,
      and
      \begin{equation}
      f_\mathrm{trans}(r) = \left( 1 + \left[\frac{r}{r_\mathrm{trans}}\right]^4 \right)^{-2}.
      \end{equation}
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalDecorator) :: massDistributionSphericalAccretionFlow
     !!{
     Implementation of an accretion flow spherical mass distribution.
     !!}
     private
     class           (massDistributionClass), pointer :: massDistributionAccretionFlow_ => null()
     double precision                                 :: radiusTransition
   contains
     !![
     <methods>
       <method method="transitionFunction" description="The transition function between halo and accretion flow."/>
     </methods>
     !!]
     final     ::                          sphericalAccretionFlowDestructor
     procedure :: density               => sphericalAccretionFlowDensity
     procedure :: densityGradientRadial => sphericalAccretionFlowDensityGradientRadial
     procedure :: energyPotential       => sphericalAccretionFlowEnergyPotential
     procedure :: energyKinetic         => sphericalAccretionFlowEnergyKinetic
     procedure :: transitionFunction    => sphericalAccretionFlowTransitionFunction
  end type massDistributionSphericalAccretionFlow

  interface massDistributionSphericalAccretionFlow
     !!{
     Constructors for the \refClass{massDistributionSphericalAccretionFlow} mass distribution class.
     !!}
     module procedure sphericalAccretionFlowConstructorParameters
     module procedure sphericalAccretionFlowConstructorInternal
  end interface massDistributionSphericalAccretionFlow

contains

  function sphericalAccretionFlowConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalAccretionFlow} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalAccretionFlow)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (massDistributionClass                 ), pointer       :: massDistribution_, massDistributionAccretionFlow_
    type            (varying_string                        )                :: nonAnalyticSolver
    double precision                                                        :: radiusTransition
    type            (varying_string                        )                :: componentType    , massType

    !![
    <inputParameter>
      <name>radiusTransition</name>
      <source>parameters</source>
      <description>The transition radius.</description>
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
    <objectBuilder class="massDistribution" name="massDistribution_"              source="parameters"                                              />
    <objectBuilder class="massDistribution" name="massDistributionAccretionFlow_" source="parameters" parameterName="massDistributionAccretionFlow"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       select type (massDistributionAccretionFlow_)
       class is (massDistributionSpherical)
          self=massDistributionSphericalAccretionFlow(radiusTransition,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,massDistributionAccretionFlow_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
       class default
          call Error_Report('a spherically-symmetric accretion flow mass distribution is required'//{introspection:location})
       end select
    class default
       call    Error_Report('a spherically-symmetric mass distribution is required'               //{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function sphericalAccretionFlowConstructorParameters
  
  function sphericalAccretionFlowConstructorInternal(radiusTransition,nonAnalyticSolver,massDistribution_,massDistributionAccretionFlow_,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalAccretionFlow} mass distribution class.
    !!}
    implicit none
    type            (massDistributionSphericalAccretionFlow)                          :: self
    class           (massDistributionSpherical             ), intent(in   ), target   :: massDistribution_, massDistributionAccretionFlow_
    double precision                                        , intent(in   )           :: radiusTransition
    type            (enumerationNonAnalyticSolversType     ), intent(in   )           :: nonAnalyticSolver
    type            (enumerationComponentTypeType          ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType               ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="radiusTransition, nonAnalyticSolver, *massDistribution_, *massDistributionAccretionFlow_, componentType, massType"/>
    !!]

    self%dimensionless=.false.
    return
  end function sphericalAccretionFlowConstructorInternal

  subroutine sphericalAccretionFlowDestructor(self)
    !!{
    Destructor for the abstract \refClass{massDistributionSphericalAccretionFlow} class.
    !!}
    implicit none
    type(massDistributionSphericalAccretionFlow), intent(inout) :: self
    
    !![
    <objectDestructor name="self%massDistribution_"             />
    <objectDestructor name="self%massDistributionAccretionFlow_"/>
    !!]
    return
  end subroutine sphericalAccretionFlowDestructor

  subroutine sphericalAccretionFlowTransitionFunction(self,radius,multiplier,multiplierGradient)
    !!{
    Return the scaled truncation radial coordinate, and the truncation multiplier.
    !!}
    implicit none
    class           (massDistributionSphericalAccretionFlow), intent(inout)           :: self
    double precision                                        , intent(in   )           :: radius
    double precision                                        , intent(  out), optional :: multiplier, multiplierGradient
    double precision                                                                  :: x
    
    x      =+     radius           &
         &  /self%radiusTransition
    if (present(multiplier        )) then
       multiplier        =+1.0d0   &
            &             /(       &
            &               +1.0d0 &
            &               +x**4  &
            &              )**2
    end if
    if (present(multiplierGradient)) then
       multiplierGradient=+8.0d0                 &
            &             *  x**3                &
            &             /(                     &
            &               +1.0d0               &
            &               +x**4                &
            &              )**3                  &
            &             /self%radiusTransition
    end if
    return
  end subroutine sphericalAccretionFlowTransitionFunction

  double precision function sphericalAccretionFlowDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalAccretionFlow), intent(inout) :: self
    class           (coordinate                            ), intent(in   ) :: coordinates
    double precision                                                        :: multiplier

    call self%transitionFunction(radius=coordinates%rSpherical(),multiplier=multiplier)
    density=+self%massDistribution_             %density(coordinates) &
         &  *       multiplier                                        &
         &  +self%massDistributionAccretionFlow_%density(coordinates) &
         &  *(1.0d0-multiplier)
    return
  end function sphericalAccretionFlowDensity

  double precision function sphericalAccretionFlowDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a accretionFlow spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalAccretionFlow), intent(inout), target   :: self
    class           (coordinate                            ), intent(in   )           :: coordinates
    logical                                                 , intent(in   ), optional :: logarithmic
    double precision                                                                  :: multiplier , multiplierGradient
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    call self%transitionFunction(radius=coordinates%rSpherical(),multiplier=multiplier,multiplierGradient=multiplierGradient)
    densityGradient=+  self%massDistribution_             %densityGradientRadial(coordinates,logarithmic=.false.) &
         &          *       multiplier                                                                            &
         &          +  self%massDistributionAccretionFlow_%densityGradientRadial(coordinates,logarithmic=.false.) &
         &          *(1.0d0-multiplier        )                                                                   &
         &          +(                                                                                            &
         &            +self%massDistribution_             %density              (coordinates                    ) &
         &            -self%massDistributionAccretionFlow_%density              (coordinates                    ) &
         &           )                                                                                            &
         &          +       multiplierGradient
    if (logarithmic_) densityGradient=+            densityGradient              &
         &                            *coordinates%rSpherical     (           ) &
         &                            /self       %density        (coordinates)
    return
  end function sphericalAccretionFlowDensityGradientRadial

  double precision function sphericalAccretionFlowEnergyPotential(self,radiusOuter) result(energy)
    !!{
    Compute the potential energy within a given {\normalfont \ttfamily radius} in a spherical accretion flow mass
    distribution. Note that this is defined to be the potential energy of the \emph{virialized} component of the mass
    distribution---the accretion flow itself is excluded.
    !!}
    implicit none
    class           (massDistributionSphericalAccretionFlow), intent(inout) :: self
    double precision                                        , intent(in   ) :: radiusOuter

    energy=self%massDistribution_%energyPotential(radiusOuter)
    return
  end function sphericalAccretionFlowEnergyPotential

  double precision function sphericalAccretionFlowEnergyKinetic(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the kinetic energy within a given {\normalfont \ttfamily radius} in a spherical accretion flow mass distribution. Note
    that this is defined to be the potential energy of the \emph{virialized} component of the mass distribution---the accretion
    flow itself is excluded.
    !!}
    implicit none
    class           (massDistributionSphericalAccretionFlow), intent(inout) :: self
    double precision                                        , intent(in   ) :: radiusOuter
    class           (massDistributionClass                 ), intent(inout) :: massDistributionEmbedding

    energy=self%massDistribution_%energyKinetic(radiusOuter,massDistributionEmbedding)
    return
  end function sphericalAccretionFlowEnergyKinetic
  
