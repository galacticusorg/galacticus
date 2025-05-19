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
  Implementation of a kinematic distribution class for heated mass distributions.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionHeated">
   <description>A heated kinematic distribution class masses.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionHeated
     !!{
     A heated kinematic distribution.
     !!}
     logical                                    :: velocityDispersionApproximate
     type   (enumerationNonAnalyticSolversType) :: nonAnalyticSolver
   contains
     procedure :: isCollisional          => heatedIsCollisional
     procedure :: velocityDispersion1D   => heatedVelocityDispersion1D
     procedure :: jeansEquationIntegrand => heatedJeansEquationIntegrand
     procedure :: jeansEquationRadius    => heatedJeansEquationRadius
  end type kinematicsDistributionHeated

  interface kinematicsDistributionHeated
     !!{
     Constructors for the {\normalfont \ttfamily heated} kinematic distribution class.
     !!}
     module procedure heatedConstructorParameters
     module procedure heatedConstructorInternal
  end interface kinematicsDistributionHeated

  ! State used to indicate whether we are solving for a self-gravitating heated system or not.
  logical :: isSelfGravitating=.true.
  !$omp threadprivate(isSelfGravitating)

contains

  function heatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily isothermal} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionHeated)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    logical                                                       :: velocityDispersionApproximate
    type            (varying_string              )                :: nonAnalyticSolver
    double precision                                              :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum

    !![
    <inputParameter>
      <name>velocityDispersionApproximate</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>
	If {\normalfont \ttfamily true}, radial velocity dispersion is computed using an approximate method in which we assume
	that $\sigma_\mathrm{r}^2(r) \rightarrow \sigma_\mathrm{r}^2(r) - (2/3) \epsilon(r)$, where $\epsilon(r)$ is the specific
	heating energy. If {\normalfont \ttfamily false} then radial velocity dispersion is computed by numerically solving the
	Jeans equation.
      </description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersion</name>
      <defaultValue>1.0d-6</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersionMaximum</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>
	Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily
	fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then
	numerical solvers are used to find solutions.
      </description>
    </inputParameter>
    !!]
    self=kinematicsDistributionHeated(enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),velocityDispersionApproximate,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function heatedConstructorParameters
  
  function heatedConstructorInternal(nonAnalyticSolver,velocityDispersionApproximate,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum) result(self)
    !!{
    Constructor for {\normalfont \ttfamily heated} kinematics distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionHeated     )                :: self
    logical                                            , intent(in   ) :: velocityDispersionApproximate
    double precision                                   , intent(in   ) :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum
    type            (enumerationNonAnalyticSolversType), intent(in   ) :: nonAnalyticSolver
    !![
    <constructorAssign variables="nonAnalyticSolver, velocityDispersionApproximate, toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum"/>
    !!]

    return
  end function heatedConstructorInternal

  logical function heatedIsCollisional(self)
    !!{
    Return false indicating that the heated kinematic distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionHeated), intent(inout) :: self
    
    heatedIsCollisional=.false.
    return
  end function heatedIsCollisional

  double precision function heatedVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an heated kinematic distribution.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class           (kinematicsDistributionHeated), intent(inout)          :: self
    class           (coordinate                  ), intent(in   )          :: coordinates
    class           (massDistributionClass       ), intent(inout), target  :: massDistribution_       , massDistributionEmbedding
    class           (massDistributionClass       )               , pointer :: massDistribution__
    double precision                                                       :: radiusInitial           , energySpecific           , &
         &                                                                    velocityDispersionSquare, radius
    type            (coordinateSpherical         )                         :: coordinatesInitial
    
    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       ! For the case of a self-gravitating heated distribution we have an optimized numerical solution for the velocity dispersion.
       select type (massDistributionEmbedding)
       class is (massDistributionSphericalHeated)
          if (massDistributionEmbedding%massDistributionHeating_%specificEnergyIsEverywhereZero() .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
             ! Use the original, unheated profile velocity dispersion.
             velocityDispersion=massDistributionEmbedding%massDistribution_%kinematicsDistribution_%velocityDispersion1D(coordinates,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)
          else if (self%velocityDispersionApproximate) then
             ! Use the approximate solution for velocity dispersion.
             radius                        =+coordinates              %rSpherical                                                                    (                                                                                                          )
             radiusInitial                 =+massDistributionEmbedding%radiusInitial                                                                 (radius                                                                                                    )
             if (radius == radiusInitial .or. (radius > 0.0d0 .and. radiusInitial == 0.0d0)) then
                ! No change in radius (either there is no heating, or solution failed - as can happen in profiles that are very
                ! close to destruction), or no solution was found. In such cases return the velocity dispersion of the unheated
                ! profile.
                coordinatesInitial         = [radius,0.0d0,0.0d0]
                velocityDispersion         =+massDistributionEmbedding%massDistribution_      %kinematicsDistribution_%velocityDispersion1D          (coordinatesInitial,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)
                return
             end if
             coordinatesInitial            = [radiusInitial,0.0d0,0.0d0]
             energySpecific                =+massDistributionEmbedding%massDistributionHeating_                        %specificEnergy               (radiusInitial     ,massDistributionEmbedding%massDistribution_                                            )
             velocityDispersionSquare      =+massDistributionEmbedding%massDistribution_       %kinematicsDistribution_%velocityDispersion1D         (coordinatesInitial,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)**2 &
                  &                         -2.0d0/3.0d0*energySpecific
             velocityDispersion            =+sqrt(max(0.0d0,velocityDispersionSquare))
             if (velocityDispersion <= 0.0d0) &
                  & velocityDispersion     =+massDistributionEmbedding%massDistribution_       %kinematicsDistribution_%velocityDispersion1D         (coordinatesInitial,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)
          else
             ! Use a numerical solution.
             velocityDispersion            =+self                                                                      %velocityDispersion1DNumerical(coordinates       ,massDistributionEmbedding                  ,massDistributionEmbedding                  )
          end if
       class default
          velocityDispersion               =+0.0d0
          call Error_Report('mass distribution must be of the `massDistributionSphericalHeated` class but found `'//char(massDistributionEmbedding%objectType())//'`'//{introspection:location})
       end select
    else
       ! Our heated distribution is embedded in another distribution. We must compute the velocity dispersion numerically. We set
       ! state to indicate that we must solve for a non-self-gravitating system.
       isSelfGravitating =.false.
       velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
       isSelfGravitating =.true.
    end if
    return
  end function heatedVelocityDispersion1D
  
  double precision function heatedJeansEquationIntegrand(self,radius,massDistribution_,massDistributionEmbedding)
    !!{
    Integrand for the Jeans equation in a heated mass distribution. Here we do the integration with respect to the initial radius
    $r_i$.
    \begin{eqnarray}
     \sigma_r(r) &=& \frac{1}{\rho(r)}\int_r^{r^{\mathrm{max}}} \rho(r) \frac{\mathrm{G} M(r)}{r^2} \mathrm{d} r \nonumber \\
                 &=& \frac{1}{\rho(r)}\int_{r_i}^{r_{i}^{\mathrm{max}}} \rho_i(r_i) \frac{\mathrm{G} M(r_i)}{r_i^2}\left(\frac{r_i}{r}\right)^4 \mathrm{d} r_i.
    \end{eqnarray}
    Here $r$ can be written as a function of $r_i$
    \begin{equation}
     r=\frac{1}{1/r_i-2\epsilon(r_i)/(\mathrm{G}M(r_i))}.
    \end{equation}
    !!}
    use :: Error                           , only : Error_Report
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (kinematicsDistributionHeated), intent(inout) :: self
    double precision                              , intent(in   ) :: radius
    class           (massDistributionClass       ), intent(inout) :: massDistribution_, massDistributionEmbedding
    double precision                                              :: radiusFinal      , energySpecific           , &
         &                                                           massEnclosed
    type            (coordinateSpherical         )                :: coordinates

    if (isSelfGravitating) then
       select type (massDistributionEmbedding)
       class is (massDistributionSphericalHeated)
          massEnclosed  =+massDistributionEmbedding                         %massDistribution_%massEnclosedBySphere(radius                                            )
          energySpecific=+massDistributionEmbedding%massDistributionHeating_                  %specificEnergy      (radius,massDistributionEmbedding%massDistribution_)
          radiusFinal   =+1.0d0                                                              &
               &         /(                                                                  &
               &           +1.0d0/radius                                                     &
               &           -2.0d0*energySpecific/gravitationalConstant_internal/massEnclosed &
               &          )
          if (radiusFinal > 0.0d0) then
             coordinates                 = [radius,0.0d0,0.0d0]
             heatedJeansEquationIntegrand=+gravitationalConstant_internal                                   &
                  &                       *massEnclosed                                                     &
                  &                       *massDistributionEmbedding%massDistribution_%density(coordinates) &
                  &                       / radius             **2                                          &
                  &                       *(radius/radiusFinal)**4
          else
             heatedJeansEquationIntegrand=+0.0d0
          end if
          class default
          heatedJeansEquationIntegrand   =+0.0d0
          call Error_Report('mass distribution must be of the `massDistributionSphericalHeated` class'//{introspection:location})
       end select
    else
       heatedJeansEquationIntegrand=self%kinematicsDistributionClass%jeansEquationIntegrand(radius,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function heatedJeansEquationIntegrand

  double precision function heatedJeansEquationRadius(self,radius,massDistributionEmbedding)
    !!{
    Return the radius variable used in solving the Jeans equation that corresponds to a given physical radius.
    Here we do the integration with respect to the initial radius, so return the initial radius.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (kinematicsDistributionHeated), intent(inout) :: self
    double precision                              , intent(in   ) :: radius
    class           (massDistributionClass       ), intent(inout) :: massDistributionEmbedding

    if (isSelfGravitating) then
       select type (massDistributionEmbedding)
          class is (massDistributionSphericalHeated)
          heatedJeansEquationRadius=massDistributionEmbedding%radiusInitial(radius)
          class default
          heatedJeansEquationRadius=0.0d0
          call Error_Report('mass distribution must be of the `massDistributionSphericalHeated` class'//{introspection:location})
       end select
    else
       heatedJeansEquationRadius=self%kinematicsDistributionClass%jeansEquationRadius(radius,massDistributionEmbedding)
    end if
    return
  end function heatedJeansEquationRadius
