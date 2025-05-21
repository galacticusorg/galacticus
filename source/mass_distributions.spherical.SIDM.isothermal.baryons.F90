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
  Provides a mass distribution implementing the ``isothermal'' approximation to the effects of SIDM, including the baryonic
  potential, based on the model of \cite{jiang_semi-analytic_2023}.
  !!}

  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator
  use            :: Numerical_ODE_Solvers  , only : odeSolver

  public :: sphericalSIDMIsothermalBaryonsInitializor
  
  !![
  <massDistribution name="massDistributionSphericalSIDMIsothermalBaryons">
   <description>
      Mass distributions for self-interacting dark matter following the ``isothermal'' model of \cite{jiang_semi-analytic_2023}. This
      model assumes that the dark matter within the interaction radius, $r_1$, has thermalized and can therefore be described by a
      constant velocity dispersion, $\sigma_0$. Under this assumption the spherical Jeans equation has a solution of the form:
      \begin{equation}
      \rho(r) = \rho_0 \exp\left[-\frac{\phi(r)}{\sigma_0^2}\right],
      \end{equation}
      where $\rho(r)$ is the density $\rho_0$ is the density at $r=0$, and the gravitational potential satisfies \citep{jiang_semi-analytic_2023}:
      \begin{equation}
      \nabla^2 \phi(r) = 4 \pi \mathrm{G} \left[ \rho_0 \exp \left( - \frac{\phi(r)}{\sigma_0^2} \right) + \rho_\mathrm{b}(r) \right],
      \end{equation}
      where $\rho_\mathrm{b}(r)$ is the density of the baryonic component. This second-order differential equation is solved using the boundary conditions $\phi(r=0)=0$ and
      $\mathrm{d}\phi/\mathrm{d}r(r=0)=0$. The values of $\rho_0$ and $\sigma_0$ are then found by minimizing a function      
      \begin{equation}
      \delta^2(\rho_0,\sigma_0) = \left[ \frac{\rho(r_1)}{\rho^\prime(r_1)} - 1 \right]^2 + \left[ \frac{M(r_1)}{M^\prime(r_1)} - 1 \right]^2,
      \end{equation}
      where $M(r)$ is the mass contained within radius $r$, and primes indicate the profile prior to SIDM thermalization.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalSIDM) :: massDistributionSphericalSIDMIsothermalBaryons
     !!{
     A mass distribution implementing the ``isothermal'' approximation to the effects of SIDM, including the baryonic potential,
     based on the model of \cite{jiang_semi-analytic_2023}.
     !!}
     private
     double precision                                                                 :: velocityDispersionCentral
     class           (massDistributionClass                    ), pointer             :: massDistributionBaryonic  => null()
     type            (interpolator                             ), allocatable         :: densityProfile                     , massProfile
     ! Call-back function and arguments used for as-needed initialization of the baryonic component.
     logical                                                                          :: initialized
     procedure       (sphericalSIDMIsothermalBaryonsInitializor), pointer    , nopass :: initializationFunction
     class           (*                                        ), pointer             :: initializationSelf        => null(), initializationArgument => null()
   contains
     !![
     <methods>
       <method method="setBaryonicComponent" description="Set baryonic components in the mass distribution."         />
       <method method="computeSolution"      description="Compute a solution for the isothermal core of a SIDM halo."/>
     </methods>
     !!]
     final     ::                          sphericalSIDMIsothermalBaryonsDestructor
     procedure :: setBaryonicComponent  => sphericalSIDMIsothermalBaryonsSetBaryonicComponent
     procedure :: computeSolution       => sphericalSIDMIsothermalBaryonsComputeSolution
     procedure :: density               => sphericalSIDMIsothermalBaryonsDensity
     procedure :: densityGradientRadial => sphericalSIDMIsothermalBaryonsDensityGradientRadial
     procedure :: massEnclosedBySphere  => sphericalSIDMIsothermalBaryonsMassEnclosedBySphere
     procedure :: potentialIsAnalytic   => sphericalSIDMIsothermalBaryonsPotentialIsAnalytic
     procedure :: potential             => sphericalSIDMIsothermalBaryonsPotential
  end type massDistributionSphericalSIDMIsothermalBaryons

  interface massDistributionSphericalSIDMIsothermalBaryons
     !!{
     Constructors for the \refClass{massDistributionSphericalSIDMIsothermalBaryons} mass distribution class.
     !!}
     module procedure sphericalSIDMIsothermalBaryonsConstructorParameters
     module procedure sphericalSIDMIsothermalBaryonsConstructorInternal
  end interface massDistributionSphericalSIDMIsothermalBaryons

  abstract interface 
     subroutine sphericalSIDMIsothermalBaryonsInitializor(initializationSelf,initializationArgument,massDistributionBaryonic)
       !!{
       Interface for call-back functions for as-needed initialization of the baryonic component.
       !!}
       import massDistributionClass
       class(*                    ), intent(inout), target  :: initializationSelf      , initializationArgument
       class(massDistributionClass), intent(  out), pointer :: massDistributionBaryonic
     end subroutine sphericalSIDMIsothermalBaryonsInitializor
  end interface

contains

  function sphericalSIDMIsothermalBaryonsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalSIDMIsothermalBaryons} mass distribution class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalSIDMIsothermalBaryons)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (massDistributionClass                         ), pointer       :: massDistribution_     , massDistributionBaryonic
    class           (darkMatterParticleClass                       ), pointer       :: darkMatterParticle_
    procedure       (sphericalSIDMIsothermalBaryonsInitializor     ), pointer       :: initializationFunction
    class           (*                                             ), pointer       :: initializationSelf    , initializationArgument
    type            (varying_string                                )                :: componentType         , massType                , &
         &                                                                             nonAnalyticSolver
    double precision                                                                :: timeAge

    !![
    <inputParameter>
      <name>timeAge</name>
      <source>parameters</source>
      <description>The age of the halo (in Gyr).</description>
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
    <objectBuilder class="massDistribution"   name="massDistribution_"                                                 source="parameters"/>
    <objectBuilder class="massDistribution"   name="massDistributionBaryonic" parameterName="massDistributionBaryonic" source="parameters"/>
    <objectBuilder class="darkMatterParticle" name="darkMatterParticle_"                                               source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       initializationFunction => null()
       initializationSelf     => null()
       initializationArgument => null()
       self=massDistributionSphericalSIDMIsothermalBaryons(timeAge,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,massDistributionBaryonic,darkMatterParticle_,initializationFunction,initializationSelf,initializationArgument,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"  />
    <objectDestructor name="darkMatterParticle_"/>
    !!]
    return
  end function sphericalSIDMIsothermalBaryonsConstructorParameters

  function sphericalSIDMIsothermalBaryonsConstructorInternal(timeAge,nonAnalyticSolver,massDistribution_,massDistributionBaryonic,darkMatterParticle_,initializationFunction,initializationSelf,initializationArgument,componentType,massType) result(self)
    !!{
    Internal constructor for the \refClass{massDistributionSphericalSIDMIsothermalBaryons} mass distribution class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    implicit none
    type            (massDistributionSphericalSIDMIsothermalBaryons)                          :: self
    double precision                                                , intent(in   )           :: timeAge
    class           (massDistributionSpherical                     ), intent(in   ), target   :: massDistribution_
    class           (massDistributionClass                         ), intent(in   ), target   :: massDistributionBaryonic
    class           (darkMatterParticleClass                       ), intent(in   ), target   :: darkMatterParticle_
    type            (enumerationNonAnalyticSolversType             ), intent(in   )           :: nonAnalyticSolver
    procedure       (sphericalSIDMIsothermalBaryonsInitializor     ), intent(in   ), pointer  :: initializationFunction
    class           (*                                             ), intent(in   ), pointer  :: initializationSelf      , initializationArgument
    type            (enumerationComponentTypeType                  ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType                       ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="timeAge, nonAnalyticSolver, componentType, massType, *massDistribution_, *massDistributionBaryonic, *darkMatterParticle_, */initializationFunction, */initializationSelf, */initializationArgument"/>
    !!]

    ! Validate the dark matter particle type.
    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! This is as expected.
    class default
       call Error_Report('this class expects a self-interacting dark matter particle'//{introspection:location})
    end select
    ! Initialize state.
    self%dimensionless=.false.
    self%initialized  =.not.associated(initializationFunction)
    return
  end function sphericalSIDMIsothermalBaryonsConstructorInternal

  subroutine sphericalSIDMIsothermalBaryonsDestructor(self)
    !!{
    Destructor for the abstract \refClass{massDistributionSphericalSIDMIsothermalBaryons} class.
    !!}
    implicit none
    type(massDistributionSphericalSIDMIsothermalBaryons), intent(inout) :: self
    
    !![
    <objectDestructor name="self%massDistribution_"       />
    <objectDestructor name="self%massDistributionBaryonic"/>
    <objectDestructor name="self%darkMatterParticle_"     />
    !!]
    return
  end subroutine sphericalSIDMIsothermalBaryonsDestructor

  subroutine sphericalSIDMIsothermalBaryonsSetBaryonicComponent(self)
    !!{
    Set the baryonic component properties in an adiabatically-contracted spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalSIDMIsothermalBaryons), intent(inout) :: self
    class(massDistributionClass                         ), pointer       :: massDistributionBaryonic
  
    if (.not.self%initialized) then
       call self%initializationFunction(self%initializationSelf,self%initializationArgument,massDistributionBaryonic)
       self%massDistributionBaryonic => massDistributionBaryonic
       self%initialized              =  .true.
       call self%computeSolution()
    end if
    return
  end subroutine sphericalSIDMIsothermalBaryonsSetBaryonicComponent

  subroutine sphericalSIDMIsothermalBaryonsComputeSolution(self)
    !!{
    Compute a solution for the isothermal core of an SIDM halo.
    !!}
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_ODE_Solvers           , only : odeSolver
    use :: Numerical_Ranges                , only : Make_Range                    , rangeTypeLinear
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Multidimensional_Minimizer      , only : multiDMinimizer
    implicit none
    class           (massDistributionSphericalSIDMIsothermalBaryons), intent(inout)              :: self
    integer         (c_size_t                                      ), parameter                  :: propertyCount                =   2
    integer                                                         , parameter                  :: countTable                   =1000
    double precision                                                , parameter                  :: odeToleranceAbsolute         =1.0d-3, odeToleranceRelative     =1.0d-3
    double precision                                                , parameter                  :: fractionRadiusInitial        =1.0d-6
    double precision                                                , dimension(propertyCount+1) :: properties                          , propertyScales
    double precision                                                , dimension(countTable     ) :: radiusTable                         , densityTable                    , &
         &                                                                                          massTable
    double precision                                                , dimension(propertyCount  ) :: locationMinimum
    type            (odeSolver                                     )                             :: odeSolver_
    type            (multiDMinimizer                               )                             :: minimizer_
    integer                                                                                      :: i                                   , iteration
    logical                                                                                      :: converged
    double precision                                                                             :: densityCentral                      , velocityDispersionCentral       , &
         &                                                                                          densityInteraction                  , massInteraction                 , &
         &                                                                                          radiusInteraction                   , radius                          , &
         &                                                                                          velocityDispersionInteraction       , mass                            , &
         &                                                                                          density
    type            (coordinateSpherical                           )                             :: coordinatesInteraction

    ! Find the interaction radius.
    radiusInteraction            =self%radiusInteraction                     (                      )
    coordinatesInteraction       =[radiusInteraction,0.0d0,0.0d0]
    ! Properties of the original density profile at the interaction radius.
    densityInteraction           =self%massDistribution_%density             (coordinatesInteraction)
    massInteraction              =self%massDistribution_%massEnclosedBySphere(radiusInteraction     )
    ! Find the velocity dispersion scale.
    velocityDispersionInteraction=sqrt(gravitationalConstant_internal*massInteraction/radiusInteraction)
    ! Set ODE solver  scales.
    propertyScales               =[velocityDispersionInteraction**2,velocityDispersionInteraction**2/radiusInteraction,massInteraction]
    ! Construct an ODE solver.
    odeSolver_                   =odeSolver      (propertyCount+1,sidmIsothermalODEs     ,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative,scale=propertyScales)
    ! Construct a minimizer.
    minimizer_                   =multiDMinimizer(propertyCount  ,sidmIsothermalFitMetric                                                                                                   )
    ! Seek the solution.
    call minimizer_%set(x=[0.0d0,1.0d0],stepSize=[1.0d0,1.0d0])
    iteration=0
    converged=.false.
    do while (.not.converged .and. iteration < 100)
       call minimizer_%iterate()
       iteration=iteration+1
       converged=minimizer_%testSize(toleranceAbsolute=1.0d-3)
    end do
    locationMinimum          =minimizer_%x()
    densityCentral           =exp(locationMinimum(1))*densityInteraction
    velocityDispersionCentral=    locationMinimum(2) *velocityDispersionInteraction
    ! Tabulate solutions for density and mass.
    radiusTable    =Make_Range(rangeMinimum=0.0d0,rangeMaximum=radiusInteraction,rangeNumber=countTable,rangeType=rangeTypeLinear)
    densityTable(1)=densityCentral
    massTable   (1)=0.0d0
    do i=2,countTable
       radius    =fractionRadiusInitial*radiusInteraction
       properties=0.0d0
       call odeSolver_%solve(radius,radiusTable(i),properties)
       densityTable(i)=+densityCentral                    &
            &          *exp(                              &
            &               -properties(1)                &
            &               /velocityDispersionCentral**2 &
            &              )
       massTable   (i)=+     properties(3)
    end do
    allocate(self%densityProfile)
    allocate(self%   massProfile)
    self%           densityProfile=interpolator(radiusTable,             densityTable)
    self%              massProfile=interpolator(radiusTable,                massTable)
    self%velocityDispersionCentral=                         velocityDispersionCentral
    return
    
  contains
    
    double precision function sidmIsothermalFitMetric(propertiesCentral)
      !!{
      Evaluate the fit metric.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:)               :: propertiesCentral
      double precision               , dimension(propertyCount+1) :: properties
      double precision                                            :: radius
      
      ! Extract current parameters.
      densityCentral           =exp(propertiesCentral(1))*densityInteraction
      velocityDispersionCentral=    propertiesCentral(2) *velocityDispersionInteraction
      ! Solve the ODE to r₁.
      radius    =fractionRadiusInitial*radiusInteraction
      properties=0.0d0
      call odeSolver_%solve(radius,radiusInteraction,properties)
      ! Extract density and mass at r₁.
      density=+densityCentral                    &
           &  *exp(                              &
           &       -properties(1)                &
           &       /velocityDispersionCentral**2 &
           &      )
      mass     =+   properties(3)
      ! Evaluate the fit metric.
      sidmIsothermalFitMetric=+(density/densityInteraction-1.0d0)**2 &
           &                  +(   mass/   massInteraction-1.0d0)**2
      return
    end function sidmIsothermalFitMetric
  
    integer function sidmIsothermalODEs(radius,properties,propertiesRateOfChange)
      !!{
      Define the ODE system to solve for isothermal self-interacting dark matter cores.
      !!}
      use :: Interface_GSL                   , only : GSL_Success
      use :: Numerical_Constants_Math        , only : Pi
      use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
      implicit none
      double precision                     , intent(in   )               :: radius
      double precision                     , intent(in   ), dimension(:) :: properties
      double precision                     , intent(  out), dimension(:) :: propertiesRateOfChange
      double precision                                                   :: densityDarkMatter     , densityBaryons
      type            (coordinateSpherical)                              :: coordinates
      
      coordinates                     =[radius,0.0d0,0.0d0]
      densityDarkMatter               =+densityCentral                    &
           &                           *exp(                              &
           &                                -max(properties(1),0.0d0)     &
           &                                /velocityDispersionCentral**2 &
           &                               )
      densityBaryons                  =+self%massDistributionBaryonic%density(coordinates)
      propertiesRateOfChange       (1)=+properties(2)
      propertiesRateOfChange       (2)=+4.0d0                            &
           &                           *Pi                               &
           &                           *gravitationalConstant_internal   &
           &                           *(                                &
           &                             +densityDarkMatter              &
           &                             +densityBaryons                 &
           &                            )
      if (radius > 0.0d0)                                                &
           & propertiesRateOfChange(2)=+propertiesRateOfChange(2)        &
           &                           -2.0d0                            &
           &                           *properties            (2)        &
           &                           /radius
      propertiesRateOfChange       (3)=+4.0d0                            &
           &                           *Pi                               &
           &                           *radius**2                        &
           &                           *densityDarkMatter
      sidmIsothermalODEs              = GSL_Success
      return
    end function sidmIsothermalODEs
    
  end subroutine sphericalSIDMIsothermalBaryonsComputeSolution

  double precision function sphericalSIDMIsothermalBaryonsDensity(self,coordinates) result(density)
    !!{
    Compute the density at the specified {\normalfont \ttfamily coordinates} for the {\normalfont \ttfamily sphericalSIDMIsothermalBaryons}
    mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalSIDMIsothermalBaryons), intent(inout) :: self
    class(coordinate                                    ), intent(in   ) :: coordinates

    call self%setBaryonicComponent()
    if (coordinates%rSpherical() > self%radiusInteraction()) then
       density=self%massDistribution_%density    (coordinates             )
    else
       density=self%densityProfile   %interpolate(coordinates%rSpherical())
    end if
    return
  end function sphericalSIDMIsothermalBaryonsDensity

  double precision function sphericalSIDMIsothermalBaryonsDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a truncated spherical mass distribution.
    !!}
    implicit none
    class  (massDistributionSphericalSIDMIsothermalBaryons), intent(inout) , target   :: self
    class  (coordinate                                    ), intent(in   )            :: coordinates
    logical                                                , intent(in   ) , optional :: logarithmic
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    if (coordinates%rSpherical() > self%radiusInteraction()) then
       densityGradient=self%massDistribution_%densityGradientRadial(coordinates,logarithmic)
    else
       densityGradient=self%densityProfile%derivative(coordinates%rSpherical())
       if (logarithmic_) &
            & densityGradient=+            densityGradient                                      &
            &                 *coordinates%rSpherical                (                        ) &
            &                 /self       %densityProfile%interpolate(coordinates%rSpherical())
    end if
    return
  end function sphericalSIDMIsothermalBaryonsDensityGradientRadial
  
  double precision function sphericalSIDMIsothermalBaryonsMassEnclosedBySphere(self,radius) result(mass)
    !!{   
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for the {\normalfont \ttfamily sphericalSIDMIsothermalBaryons}
    mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalSIDMIsothermalBaryons), intent(inout) , target :: self
    double precision                                                , intent(in   )          :: radius

    call self%setBaryonicComponent()
    if (radius > self%radiusInteraction()) then
       mass=self%massDistribution_%massEnclosedBySphere(radius)
    else
       mass=self%massProfile      %interpolate (radius)
    end if
    return
  end function sphericalSIDMIsothermalBaryonsMassEnclosedBySphere

  logical function sphericalSIDMIsothermalBaryonsPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return if the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionSphericalSIDMIsothermalBaryons), intent(inout) :: self

    isAnalytic=self%massDistribution_%potentialIsAnalytic()
    return
  end function sphericalSIDMIsothermalBaryonsPotentialIsAnalytic

  double precision function sphericalSIDMIsothermalBaryonsPotential(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in an burkert mass distribution.
    !!}
    use :: Coordinates               , only : coordinateSpherical      , assignment(=)
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    implicit none
    class           (massDistributionSphericalSIDMIsothermalBaryons), intent(inout), target   :: self
    class           (coordinate                                    ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType             ), intent(  out), optional :: status
    type            (coordinateSpherical                           )                          :: coordinatesInteraction
    double precision                                                                          :: radiusInteraction
    
    call self%setBaryonicComponent()
    if (present(status)) status=structureErrorCodeSuccess
    if (coordinates%rSpherical() > self%radiusInteraction()) then
       potential             =+self%massDistribution_%potential(coordinates,status=status)
    else
       radiusInteraction     =+self%radiusInteraction()
       coordinatesInteraction=[radiusInteraction,0.0d0,0.0d0]
       potential             =+self%massDistribution_%potential(coordinatesInteraction)              &
            &                 -self%velocityDispersionCentral**2                                     &
            &                 *log(                                                                  &
            &                      +self%densityProfile%interpolate(coordinates%rSpherical       ()) &
            &                      /self%densityProfile%interpolate(            radiusInteraction  ) &
            &                 )
    end if  
    return
  end function sphericalSIDMIsothermalBaryonsPotential
