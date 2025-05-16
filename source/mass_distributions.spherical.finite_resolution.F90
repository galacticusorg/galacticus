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
  Implements a finite resolution spherical mass distribution.
  !!}

  !![
  <massDistribution name="massDistributionSphericalFiniteResolution">
   <description>
     A mass distribution class which applies a finite resolution to some other mass distribution class, typically to mimic the
     effects of finite resolution in an N-body simulation. Specifically, the density profile is given by
     \begin{equation}
     \rho(r) = \rho^\prime(r) \left( 1 + \left[ \frac{\Delta x}{r} \right]^2 \right)^{-1/2},
     \end{equation}
     where $\Delta x$ is the larger of the resolution length, {\normalfont \ttfamily [lengthResolution]}, and the radius in the
     original profile enclosing the mass resolution, {\normalfont \ttfamily [massResolution]}.
     
     Note that this choice was constructed to give a constant density core in an NFW density profile. For a density profile,
     $\rho^\prime(r)$, which rises more steeply than $r^{-1}$ as $r \rightarrow 0$ we will still have a cuspy density profile
     under this model.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalDecorator) :: massDistributionSphericalFiniteResolution
     !!{
     Implementation of a finite resolution spherical mass distribution.
     !!}
     private
     double precision :: lengthResolution
   contains
     final     ::                          sphericalFiniteResolutionDestructor
     procedure :: density               => sphericalFiniteResolutionDensity
     procedure :: densityGradientRadial => sphericalFiniteResolutionDensityGradientRadial
  end type massDistributionSphericalFiniteResolution

  interface massDistributionSphericalFiniteResolution
     !!{
     Constructors for the \refClass{massDistributionSphericalFiniteResolution} mass distribution class.
     !!}
     module procedure sphericalFiniteResolutionConstructorParameters
     module procedure sphericalFiniteResolutionConstructorInternal
  end interface massDistributionSphericalFiniteResolution

contains

  function sphericalFiniteResolutionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalFiniteResolution} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalFiniteResolution)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (massDistributionClass                    ), pointer       :: massDistribution_
    type            (varying_string                           )                :: nonAnalyticSolver
    double precision                                                           :: lengthResolution
    type            (varying_string                           )                :: componentType    , massType

    !![
    <inputParameter>
      <name>lengthResolution</name>
      <source>parameters</source>
      <description>The resolution length scale.</description>
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
       self=massDistributionSphericalFiniteResolution(lengthResolution,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function sphericalFiniteResolutionConstructorParameters
  
  function sphericalFiniteResolutionConstructorInternal(lengthResolution,nonAnalyticSolver,massDistribution_,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalFiniteResolution} mass distribution class.
    !!}
    implicit none
    type            (massDistributionSphericalFiniteResolution)                          :: self
    class           (massDistributionSpherical                ), intent(in   ), target   :: massDistribution_
    type            (enumerationNonAnalyticSolversType        ), intent(in   )           :: nonAnalyticSolver
    double precision                                           , intent(in   )           :: lengthResolution
    type            (enumerationComponentTypeType             ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType                  ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="lengthResolution, nonAnalyticSolver, *massDistribution_, componentType, massType"/>
    !!]

    self%dimensionless=self%massDistribution_%isDimensionless()
    return
  end function sphericalFiniteResolutionConstructorInternal

  subroutine sphericalFiniteResolutionDestructor(self)
    !!{
    Destructor for the abstract \refClass{massDistributionSphericalFiniteResolution} class.
    !!}
    implicit none
    type(massDistributionSphericalFiniteResolution), intent(inout) :: self
    
    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine sphericalFiniteResolutionDestructor

  double precision function sphericalFiniteResolutionDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalFiniteResolution), intent(inout)           :: self
    class(coordinate                               ), intent(in   )           :: coordinates
 
    density=+self%massDistribution_%density(coordinates) &
         &  /sqrt(                                       &
         &        +1.0d0                                 &
         &        +(                                     &
         &          +self       %lengthResolution        &
         &          /coordinates%rSpherical      ()      &
         &         )**2                                  &
         &       )
      return
  end function sphericalFiniteResolutionDensity

  double precision function sphericalFiniteResolutionDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a finiteResolution spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSphericalFiniteResolution), intent(inout), target   :: self
    class  (coordinate                               ), intent(in   )           :: coordinates
    logical                                           , intent(in   ), optional :: logarithmic
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    densityGradient=+self%massDistribution_%densityGradientRadial(coordinates,logarithmic=.true.) &
         &          +(                                                                            &
         &            +  self       %lengthResolution                                             &
         &            /  coordinates%rSpherical      ()                                           &
         &           )  **2                                                                       &
         &          /(                                                                            &
         &            +1.0d0                                                                      &
         &            +(                                                                          &
         &              +self       %lengthResolution                                             &
         &              /coordinates%rSpherical      ()                                           &
         &             )**2                                                                       &
         &           )
    if (.not.logarithmic_) &
         densityGradient=+            densityGradient              &
         &               *self       %density        (coordinates) &
         &               /coordinates%rSpherical     (           )
   return
  end function sphericalFiniteResolutionDensityGradientRadial
