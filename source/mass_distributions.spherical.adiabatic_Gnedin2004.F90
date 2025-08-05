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
  Implements an adiabatically-contracted spherical mass distribution.
  !!}

  use :: Math_Exponentiation, only : fastExponentiator
  use :: Root_Finder        , only : rootFinder

  public :: sphericalAdiabaticGnedin2004Initializor
  
  ! Number of previous radius solutions to store.
  !![
  <scoping>
   <module variables="sphericalAdiabaticGnedin2004StoreCount"/>
  </scoping>
  !!]
  integer, parameter :: sphericalAdiabaticGnedin2004StoreCount=10
  
  !![
  <massDistribution name="massDistributionSphericalAdiabaticGnedin2004">
   <description>
    A dark matter profile class which applies adiabatic contraction to the halo as it responds to the presence of
    baryons. Adiabatic contraction follows the algorithm of \cite{gnedin_response_2004}. The parameters $A$ and $\omega$ of
    that model are specified via input parameters {\normalfont \ttfamily A} and {\normalfont \ttfamily omega} respectively.
  
    Given the final radius, $r_\mathrm{f}$, the corresponding initial radius, $r_\mathrm{i}$, is found by solving:
    \begin{equation}
    f_\mathrm{i} M_\mathrm{total,0}(\bar{r}_\mathrm{i}) r_\mathrm{i} = f_\mathrm{f} M_\mathrm{total,0}(\bar{r}_\mathrm{i})
    r_\mathrm{f} + V^2_\mathrm{b}(\bar{r}_\mathrm{f}) \bar{r}_\mathrm{f} r_\mathrm{f}/ \mathrm{G},
     \label{eq:adiabaticContractionGnedinSolution}
    \end{equation}
    where $M_\mathrm{total,0}(r)$ is the initial total matter profile, $V_\mathrm{b}(r)$ is the baryonic contribution to the
    rotation curve, $f_\mathrm{i}$, is the fraction of mass within the virial radius compared to the node mass\footnote{In
    \protect\glc\ the ``node mass'' refers to the total mass of the node, assuming it has the universal complement of
    baryons. Since some halos may contain less than the complete complement of baryons it is possible that $f_\mathrm{i}&lt;1$.},
    $f_\mathrm{f}=(\Omega_\mathrm{M}-\Omega_\mathrm{b})/\Omega_\mathrm{M}+M_\mathrm{satellite, baryonic}/M_\mathrm{total}$,
    $M_\mathrm{satellite, baryonic}$ is the baryonic mass in any satellite halos, $M_\mathrm{total}$ is the node mass, and
    \begin{equation}
    {\bar{r} \over r_0} = A \left({r \over r_0}\right)^\omega,
    \label{eq:adiabaticContractionGnedinPowerLaw}
    \end{equation}    
    where the pivot radius $r_0$ is set to $f_0 r_\mathrm{vir}$ where $f_0=${\normalfont \ttfamily [radiusFractionalPivot]}, and
    $r_\mathrm{vir}$ is the virial radius. The original \cite{gnedin_response_2004} assumed $f_0=1$, but the revised model of
    \cite{gnedin_halo_2011} found that $f_0=0.03$ lead to an improved model (less scatter in the best fit values of $(A,\omega)$
    when comparing to N-body simulations).

    Note that we explicitly assume that the initial, uncontracted total density profile has the same shape as the initial dark
    matter density profile, that contraction of the halo occurs with no shell crossing, and that satellite halos trace the dark
    matter profile of their host halo.  The derivative, $\mathrm{d} r_\mathrm{f}/\mathrm{d}d_\mathrm{i}\equiv r^\prime_\mathrm{i}$
    is found by taking the derivative of eqn.~(\ref{eq:adiabaticContractionGnedinSolution}) to give:    
    \begin{eqnarray}
     &amp; &amp; f_\mathrm{i} M_\mathrm{total,0}(\bar{r}_\mathrm{i}) r^\prime_\mathrm{i} + f_\mathrm{i} 4 \pi
     \bar{r}_\mathrm{i}^2 \rho_\mathrm{total,0}(\bar{r}_\mathrm{i}) {\mathrm{d} \bar{r}_\mathrm{i}\over\mathrm{d} r_\mathrm{i}}
     r_\mathrm{i} r^\prime_\mathrm{i} \nonumber \\
     &amp; = &amp; f_\mathrm{f} M_\mathrm{total,0}(\bar{r}_\mathrm{i}) + f_\mathrm{i} 4 \pi \bar{r}_\mathrm{i}^2
     \rho_\mathrm{total,0}(\bar{r}_\mathrm{i}) {\mathrm{d} \bar{r}_\mathrm{i}\over\mathrm{d} r_\mathrm{i}} r_\mathrm{f}
     r^\prime_\mathrm{i} \nonumber \\
     &amp; + &amp; V^2_\mathrm{b}(\bar{r}_\mathrm{f}) \bar{r}_\mathrm{f} / \mathrm{G} + V^2_\mathrm{b}(\bar{r}_\mathrm{f})
     {\mathrm{d}\bar{r}_\mathrm{f}\over \mathrm{d} r_\mathrm{f}} r_\mathrm{f}/ \mathrm{G} +
     {\mathrm{d}V^2_\mathrm{b}\over\mathrm{d} \bar{r}_\mathrm{f}}(\bar{r}_\mathrm{f}) {\mathrm{d}\bar{r}_\mathrm{f}\over
     \mathrm{d} r_\mathrm{f}} \bar{r}_\mathrm{f} r_\mathrm{f}/ \mathrm{G},
    \end{eqnarray}
    where
    \begin{equation}
     {\mathrm{d}\bar{r} \over \mathrm{d} r} = A \left({r \over r_0}\right)^{\omega-1},
    \end{equation}
    and which can then be solved numerically for $r^\prime_\mathrm{i}$.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalDecorator) :: massDistributionSphericalAdiabaticGnedin2004
     !!{
     An adiabatically-contracted spherical mass distribution.
     !!}
     private
     class           (massDistributionClass                  ), pointer                                                  :: massDistributionBaryonic      => null()
     type            (rootFinder                             )                                                           :: finder
     ! Parameters of the adiabatic contraction algorithm.
     double precision                                                                                                    :: A                                      , omega                               , &
          &                                                                                                                 radiusFractionalPivot
     ! Stored solutions for reuse.
     integer                                                                                                             :: radiusPreviousIndex                    , radiusPreviousIndexMaximum
     double precision                                         , dimension(sphericalAdiabaticGnedin2004StoreCount)        :: radiusPrevious                         , radiusInitialPrevious
     ! Quantities used in solving the initial radius root function.
     double precision                                                                                                    :: baryonicFinalTerm                      , baryonicFinalTermDerivative         , &
          &                                                                                                                 darkMatterDistributedFraction          , massFractionInitial                 , &
          &                                                                                                                 radiusFinal                            , radiusFinalMean            ,          &
          &                                                                                                                 darkMatterFraction                     , radiusVirial                        , &
          &                                                                                                                 toleranceRelative                      , massTotal_
     ! Call-back function and arguments used for as-needed initialization of the baryonic component.
     logical                                                                                                              :: initialized
     procedure       (sphericalAdiabaticGnedin2004Initializor), pointer                                          , nopass :: initializationFunction
     class           (*                                      ), pointer                                                   :: initializationSelf           => null(), initializationArgument      => null()
   contains
     !![
     <methods>
       <method description="Set baryonic components in the mass distribution."                        method="setBaryonicComponent"       />
       <method description="Compute factors needed for solving adiabatic contraction."                method="computeFactors"             />
       <method description="Compute the orbit-averaged radius for dark matter."                       method="radiusOrbitalMean"          />
       <method description="Compute the derivative of the orbit-averaged radius for dark matter."     method="radiusOrbitalMeanDerivative"/>
       <method description="Compute the initial radius in the dark matter profile."                   method="radiusInitial"              />
       <method description="Compute the derivative of the initial radius in the dark matter profile." method="radiusInitialDerivative"    />
     </methods>
     !!]
     final     ::                                sphericalAdiabaticGnedin2004Destructor
     procedure :: setBaryonicComponent        => sphericalAdiabaticGnedin2004SetBaryonicComponent
     procedure :: density                     => sphericalAdiabaticGnedin2004Density
     procedure :: massEnclosedBySphere        => sphericalAdiabaticGnedin2004MassEnclosedBySphere
     procedure :: radiusInitial               => sphericalAdiabaticGnedin2004RadiusInitial
     procedure :: radiusInitialDerivative     => sphericalAdiabaticGnedin2004RadiusInitialDerivative
     procedure :: computeFactors              => sphericalAdiabaticGnedin2004ComputeFactors
     procedure :: radiusOrbitalMean           => sphericalAdiabaticGnedin2004RadiusOrbitalMean
     procedure :: radiusOrbitalMeanDerivative => sphericalAdiabaticGnedin2004RadiusOrbitalMeanDerivative
  end type massDistributionSphericalAdiabaticGnedin2004

  interface massDistributionSphericalAdiabaticGnedin2004
     !!{
     Constructors for the \refClass{massDistributionSphericalAdiabaticGnedin2004} mass distribution class.
     !!}
     module procedure sphericalAdiabaticGnedin2004ConstructorParameters
     module procedure sphericalAdiabaticGnedin2004ConstructorInternal
  end interface massDistributionSphericalAdiabaticGnedin2004
    
  ! Module-scope quantities used in solving the initial radius root function.
  double precision                                              , parameter   :: toleranceAbsolute  =0.0d0
  class           (massDistributionSphericalAdiabaticGnedin2004), pointer     :: self_
  !$omp threadprivate(self_)

  ! Module-scope shared fast exponentiator.
  type            (fastExponentiator                           ), allocatable :: radiusExponentiator
  double precision                                                            :: omegaPrevious      =-huge(0.0d0)
  !$omp threadprivate(radiusExponentiator,omegaPrevious)

  abstract interface 
     subroutine sphericalAdiabaticGnedin2004Initializor(initializationSelf,initializationArgument,massDistributionBaryonic,darkMatterDistributedFraction,massFractionInitial)
       !!{
       Interface for call-back functions for as-needed initialization of the baryonic component.
       !!}
       import massDistributionClass
       class           (*                    ), intent(inout), target  :: initializationSelf           , initializationArgument
       class           (massDistributionClass), intent(  out), pointer :: massDistributionBaryonic
       double precision                       , intent(  out)          :: darkMatterDistributedFraction, massFractionInitial
     end subroutine sphericalAdiabaticGnedin2004Initializor
  end interface

contains

  function sphericalAdiabaticGnedin2004ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalAdiabaticGnedin2004} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalAdiabaticGnedin2004)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (massDistributionClass                       ), pointer       :: massDistribution_            , massDistributionBaryonic
    procedure       (sphericalAdiabaticGnedin2004Initializor     ), pointer       :: initializationFunction
    class           (*                                           ), pointer       :: initializationSelf           , initializationArgument
    double precision                                                              :: A                            , omega                   , &
         &                                                                           radiusFractionalPivot        , toleranceRelative       , &
         &                                                                           radiusVirial                 , darkMatterFraction      , &
         &                                                                           darkMatterDistributedFraction, massFractionInitial
    type            (varying_string                              )                :: componentType                , massType                , &
         &                                                                           nonAnalyticSolver

    !![
    <inputParameter>
      <name>A</name>
      <defaultSource>(\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultSource>
      <defaultValue>0.80d0</defaultValue>
      <description>The parameter $A$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>omega</name>
      <defaultSource>(\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultSource>
      <defaultValue>0.77d0</defaultValue>
      <description>The parameter $\omega$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusFractionalPivot</name>
      <defaultSource>\citep{gnedin_response_2004}</defaultSource>
      <defaultValue>1.0d0</defaultValue>
      <description>The pivot radius (in units of the virial radius), $r_0$, appearing in equation~(\ref{eq:adiabaticContractionGnedinPowerLaw}).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusVirial</name>
      <description>The virial radius, $r_\mathrm{v}$, appearing in equation~(\ref{eq:adiabaticContractionGnedinPowerLaw}).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>darkMatterFraction</name>
      <description>The universal dark matter fraction.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>darkMatterDistributedFraction</name>
      <description>The fraction of matter assumed to be distributed as the dark matter.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massFractionInitial</name>
      <description>The fraction of matter assumed to be initially distributed as the dark matter.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelative</name>
      <defaultValue>1.0d-2</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in solving for the initial radius in the adiabatically-contracted dark matter profile.</description>
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
    <objectBuilder class="massDistribution" name="massDistribution_"                                                 source="parameters"/>
    <objectBuilder class="massDistribution" name="massDistributionBaryonic" parameterName="massDistributionBaryonic" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       initializationFunction => null()
       initializationSelf     => null()
       initializationArgument => null()
       self=massDistributionSphericalAdiabaticGnedin2004(A,omega,radiusVirial,radiusFractionalPivot,darkMatterFraction,darkMatterDistributedFraction,massFractionInitial,toleranceRelative,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,massDistributionBaryonic,initializationFunction,initializationSelf,initializationArgument,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"       />
    <objectDestructor name="massDistributionBaryonic"/>
    !!]
    return
  end function sphericalAdiabaticGnedin2004ConstructorParameters
  
  function sphericalAdiabaticGnedin2004ConstructorInternal(A,omega,radiusVirial,radiusFractionalPivot,darkMatterFraction,darkMatterDistributedFraction,massFractionInitial,toleranceRelative,nonAnalyticSolver,massDistribution_,massDistributionBaryonic,initializationFunction,initializationSelf,initializationArgument,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalAdiabaticGnedin2004} mass distribution class.
    !!}
    implicit none
    type            (massDistributionSphericalAdiabaticGnedin2004)                          :: self
    double precision                                              , intent(in   )           :: A                       , omega                        , &
         &                                                                                     radiusVirial            , radiusFractionalPivot        , &
         &                                                                                     darkMatterFraction      , darkMatterDistributedFraction, &
         &                                                                                     massFractionInitial     , toleranceRelative
    class           (massDistributionSpherical                   ), intent(in   ), target   :: massDistribution_
    class           (massDistributionClass                       ), intent(in   ), target   :: massDistributionBaryonic
    type            (enumerationNonAnalyticSolversType           ), intent(in   )           :: nonAnalyticSolver
    procedure       (sphericalAdiabaticGnedin2004Initializor     ), intent(in   ), pointer  :: initializationFunction
    class           (*                                           ), intent(in   ), pointer  :: initializationSelf      , initializationArgument
    type            (enumerationComponentTypeType                ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType                     ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="A, omega, radiusVirial, radiusFractionalPivot, darkMatterFraction, darkMatterDistributedFraction, massFractionInitial, toleranceRelative, nonAnalyticSolver, *massDistribution_, *massDistributionBaryonic, */initializationFunction, */initializationSelf, */initializationArgument, componentType, massType"/>
    !!]

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    ! Evaluate the original total mass.
    self%massTotal_=self%massDistribution_%massEnclosedBySphere(radiusVirial)
    ! Construct a root finder.
    self%finder=rootFinder(                                                      &
         &                 rootFunction     =sphericalAdiabaticGnedin2004Solver, &
         &                 toleranceAbsolute=toleranceAbsolute                 , &
         &                 toleranceRelative=toleranceRelative                   &
         &                )    
    self%dimensionless=self%massDistribution_%isDimensionless()
    ! Initialize state.
    self%radiusPreviousIndex       = 0
    self%radiusPreviousIndexMaximum= 0
    self%radiusPrevious            =-1.0d0
    self%initialized               =.not.associated(initializationFunction)
    return
  end function sphericalAdiabaticGnedin2004ConstructorInternal

  subroutine sphericalAdiabaticGnedin2004Destructor(self)
    !!{
    Destructor for the abstract \refClass{massDistributionSphericalAdiabaticGnedin2004} class.
    !!}
    implicit none
    type(massDistributionSphericalAdiabaticGnedin2004), intent(inout) :: self
    
    !![
    <objectDestructor name="self%massDistribution_"       />
    <objectDestructor name="self%massDistributionBaryonic"/>
    !!]
    return
  end subroutine sphericalAdiabaticGnedin2004Destructor

  subroutine sphericalAdiabaticGnedin2004SetBaryonicComponent(self)
    !!{
    Set the baryonic component properties in an adiabatically-contracted spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalAdiabaticGnedin2004), intent(inout) :: self
    class           (massDistributionClass                       ), pointer       :: massDistributionBaryonic
    double precision                                                              :: darkMatterDistributedFraction, massFractionInitial

    if (.not.self%initialized) then
    call self%initializationFunction(self%initializationSelf,self%initializationArgument,massDistributionBaryonic,darkMatterDistributedFraction,massFractionInitial)
    self%massDistributionBaryonic      => massDistributionBaryonic
    self%darkMatterDistributedFraction =  darkMatterDistributedFraction
    self%massFractionInitial           =  massFractionInitial
    self%initialized                   =  .true.
    end if
    return
  end subroutine sphericalAdiabaticGnedin2004SetBaryonicComponent
  
  double precision function sphericalAdiabaticGnedin2004Density(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an adiabatically-contracted spherical mass distribution.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class           (massDistributionSphericalAdiabaticGnedin2004), intent(inout) :: self
    class           (coordinate                                  ), intent(in   ) :: coordinates
    type            (coordinateSpherical                         )                :: coordinatesInitial
    double precision                                                              :: radiusInitial     , radiusInitialDerivative, &
         &                                                                           densityInitial

    radiusInitial     =self                  %radiusInitial(coordinates       %rSpherical())
    coordinatesInitial=[radiusInitial,0.0d0,0.0d0]
    densityInitial    =self%massDistribution_%density      (coordinatesInitial             )
    if (coordinates%rSpherical() == radiusInitial) then
       density=+self%darkMatterFraction &
            &  *     densityInitial
    else
       radiusInitialDerivative=+self%radiusInitialDerivative(coordinates%rSpherical())
       density                =+self%darkMatterFraction       &
            &                  *densityInitial                &
            &                  *(                             &
            &                    +            radiusInitial   &
            &                    /coordinates%rSpherical   () &
            &                   )**2                          &
            &                  *radiusInitialDerivative
    end if
    return
  end function sphericalAdiabaticGnedin2004Density

  double precision function sphericalAdiabaticGnedin2004MassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for adiabatically-contracted mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalAdiabaticGnedin2004), intent(inout), target :: self
    double precision                                              , intent(in   )         :: radius
    
    mass   =+self                  %darkMatterFraction                               &
         &  *self%massDistribution_%massEnclosedBySphere(self%radiusInitial(radius))
    return
  end function sphericalAdiabaticGnedin2004MassEnclosedBySphere

  double precision function sphericalAdiabaticGnedin2004RadiusInitial(self,radius)
    !!{
    Compute the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    \cite{gnedin_response_2004}.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (massDistributionSphericalAdiabaticGnedin2004), intent(inout) :: self
    double precision                                              , intent(in   ) :: radius
    integer                                                                       :: i               , j           , &
         &                                                                           iMod
    double precision                                                              :: radiusUpperBound, massEnclosed

    call self%setBaryonicComponent()
    ! Handle non-positive radii.
    if (radius <= 0.0d0) then
       sphericalAdiabaticGnedin2004RadiusInitial=0.0d0
       return
    end if
    ! Check for a previously computed solution.
    if (self%radiusPreviousIndexMaximum > 0 .and. any(self%radiusPrevious(1:self%radiusPreviousIndexMaximum) == radius)) then
       sphericalAdiabaticGnedin2004RadiusInitial=0.0d0
       do i=1,self%radiusPreviousIndexMaximum
          if (self%radiusPrevious(i) == radius) then
             sphericalAdiabaticGnedin2004RadiusInitial=self%radiusInitialPrevious(i)
             exit
          end if
       end do
       return
    end if
    ! Return radius unchanged if larger than the virial radius.
    if (radius >= self%radiusVirial) then
       sphericalAdiabaticGnedin2004RadiusInitial=radius
       return
    end if
    ! Compute the various factors needed by this calculation.
    call self%computeFactors(radius,computeGradientFactors=.false.)
    !! Note that even if no baryons are present at this radius we can not assume that the initial radius is unchanged because it
    !! is possible that the initial fraction of baryons and the fraction of mass distributed as the dark matter are not equal, fᵢ
    !! ≠ fᵪ.
    ! Check that solution is within bounds.
    if (sphericalAdiabaticGnedin2004Solver(self%radiusVirial) < 0.0d0) then
       sphericalAdiabaticGnedin2004RadiusInitial=self%radiusVirial
       return
    end if
    j=-1
    if (self%radiusPreviousIndexMaximum > 0) then
       ! No exact match exists, look for approximate matches.
       do i=1,self%radiusPreviousIndexMaximum
          iMod=modulo(self%radiusPreviousIndex-i,sphericalAdiabaticGnedin2004StoreCount)+1
          if (abs(radius-self%radiusPrevious(iMod))/self%radiusPrevious(iMod) < self%toleranceRelative) then
             j=iMod
             exit
          end if
       end do
    end if
    ! Find the solution for initial radius.
    if (j == -1) then
       ! No previous solution to use as an initial guess. Instead, we make an estimate of the initial radius under the
       ! assumption that the mass of dark matter (in the initial profile) enclosed within the mean initial radius is the
       ! same as enclosed within the mean final radius. Since the initial and final radii are typically not too
       ! different, and since the mean radius is a weak (ѡ<1) function of the radius this is a useful
       ! approximation. Furthermore, since it will underestimate the actual mass within the initial mean radius it gives
       ! an overestimate of the initial radius. This means that we have a bracketing of the initial radius which we can
       ! use in the solver.
       massEnclosed=+self%massDistribution_%massEnclosedBySphere(self%radiusOrbitalMean(self%radiusFinal))
       if (massEnclosed > 0.0d0) then
          radiusUpperBound=+(                                    &
               &             +self%baryonicFinalTerm             &
               &             /     massEnclosed                  &
               &             +self%darkMatterDistributedFraction &
               &             *self%radiusFinal                   &
               &            )                                    &
               &           /  self%massFractionInitial
          if (radiusUpperBound < radius) radiusUpperBound=radius
       else
          radiusUpperBound=radius
       end if
       call self%finder%rangeExpand(                                                             &
            &                       rangeExpandUpward            =1.1d0                        , &
            &                       rangeExpandDownward          =0.9d0                        , &
            &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                       rangeExpandType              =rangeExpandMultiplicative      &
            &                      )
       sphericalAdiabaticGnedin2004RadiusInitial=self%finder%find(rootRange=[radius,radiusUpperBound])
    else
       ! Use previous solution as an initial guess.
       call self%finder%rangeExpand(                                                                        &
            &                       rangeExpandDownward          =1.0d0/sqrt(1.0d0+self%toleranceRelative), &
            &                       rangeExpandUpward            =1.0d0*sqrt(1.0d0+self%toleranceRelative), &
            &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative           , &
            &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive           , &
            &                       rangeExpandType              =rangeExpandMultiplicative                 &
            &                      )
       sphericalAdiabaticGnedin2004RadiusInitial=self%finder%find(                                                                             &
            &                                                     rootRange=[                                                                  &
            &                                                                self%radiusInitialPrevious(j)/sqrt(1.0d0+self%toleranceRelative), &
            &                                                                self%radiusInitialPrevious(j)*sqrt(1.0d0+self%toleranceRelative)  &
            &                                                               ]                                                                  &
            &                                                    )
    end if
    ! Store this solution.
    self%radiusPreviousIndex                                 =modulo(self%radiusPreviousIndex         ,sphericalAdiabaticGnedin2004StoreCount)+1
    self%radiusPreviousIndexMaximum                          =min   (self%radiusPreviousIndexMaximum+1,sphericalAdiabaticGnedin2004StoreCount)
    self%radiusPrevious            (self%radiusPreviousIndex)=radius
    self%radiusInitialPrevious     (self%radiusPreviousIndex)=sphericalAdiabaticGnedin2004RadiusInitial
    return
  end function sphericalAdiabaticGnedin2004RadiusInitial

  double precision function sphericalAdiabaticGnedin2004RadiusInitialDerivative(self,radius) result(radiusInitialDerivative)
    !!{
    Compute the derivative of the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    \cite{gnedin_response_2004}.
    !!}
    use :: Display                 , only : displayMessage     , displayIndent, displayUnindent
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    use :: Coordinates             , only : coordinateSpherical, assignment(=)
    implicit none
    class           (massDistributionSphericalAdiabaticGnedin2004), intent(inout) :: self
    double precision                                              , intent(in   ) :: radius
    type            (varying_string                              ), save          :: message
    !$omp threadprivate(message)
    character       (len=12                                      )                :: label
    double precision                                                              :: radiusInitial                  , radiusInitialMean            , &
         &                                                                           massDarkMatterInitial          , densityDarkMatterInitial     , &
         &                                                                           radiusInitialMeanSelfDerivative, radiusFinalMeanSelfDerivative, &
         &                                                                           numerator                      , denominator
    type            (coordinateSpherical                         )                :: coordinatesInitialMean
    
    call self%setBaryonicComponent()
    ! Compute the various factors needed by this calculation.
    call self%computeFactors(radius,computeGradientFactors=.true.)
    ! Return unit derivative if radius is larger than the virial radius.
    if (radius >= self%radiusVirial) then
       radiusInitialDerivative=1.0d0
       return
    end if
    ! Validate.
    if (radius <= 0.0d0) call Error_Report('non-positive radius')
    ! Compute initial radius, and derivatives of initial and final mean radii.
    radiusInitial                  =self                  %radiusInitial              (radius                )
    radiusInitialMeanSelfDerivative=self                  %radiusOrbitalMeanDerivative(radiusInitial         )
    radiusFinalMeanSelfDerivative  =self                  %radiusOrbitalMeanDerivative(radius                )
    ! Find the initial mean orbital radius.
    radiusInitialMean              =self                  %radiusOrbitalMean          (radiusInitial         )
    coordinatesInitialMean         =[radiusInitialMean,0.0d0,0.0d0]
    ! Get the mass of dark matter inside the initial radius.
    massDarkMatterInitial          =self%massDistribution_%massEnclosedBySphere       (radiusInitialMean     )
    ! Get the mass of dark matter inside the initial radius.
    densityDarkMatterInitial       =self%massDistribution_%density                    (coordinatesInitialMean)    
    ! Find the solution for the derivative of the initial radius.
    numerator                      =+(                                                      &
         &                            +massDarkMatterInitial                                &
         &                            *self%darkMatterDistributedFraction                   &
         &                            +self%baryonicFinalTerm                               &
         &                            *(                                                    &
         &                              +1.0d0                        /     radius          &
         &                              +radiusFinalMeanSelfDerivative/self%radiusFinalMean &
         &                             )                                                    &
         &                            +self%baryonicFinalTermDerivative                     &
         &                           )
    denominator                    =+(                                                      &
         &                            +massDarkMatterInitial*self%massFractionInitial       &
         &                            +(                                                    &
         &                              +self%massFractionInitial          *radiusInitial   &
         &                              -self%darkMatterDistributedFraction*radius          &
         &                             )                                                    &
         &                            *4.0d0                                                &
         &                            *Pi                                                   &
         &                            *radiusInitialMean**2                                 &
         &                            *densityDarkMatterInitial                             &
         &                            *radiusInitialMeanSelfDerivative                      &
         &                           )
    if (exponent(numerator)-exponent(denominator) > maxExponent(0.0d0)) then
       call displayIndent  ('Radius derivative calculation')
       write (label,'(e12.6)')      radius
       message='r_final                       = '//label//' Mpc'
       call displayMessage (message)
       write (label,'(e12.6)')      radiusInitial
       message='r_initial                     = '//label//' Mpc'
       call displayMessage (message)
       write (label,'(e12.6)') self%radiusFinalMean
       message='⟨r_final⟩                     = '//label//' Mpc'
       call displayMessage (message)
       write (label,'(e12.6)')      radiusInitialMean
       message='⟨r_initial⟩                   = '//label//' Mpc'
       call displayMessage (message)
       write (label,'(e12.6)')      radiusInitialMeanSelfDerivative
       message='d⟨r_initial⟩/dr_initial       = '//label
       call displayMessage (message)
       write (label,'(e12.6)')      radiusFinalMeanSelfDerivative
       message='d⟨r_final⟩/dr_final           = '//label
       call displayMessage (message)
       write (label,'(e12.6)')      massDarkMatterInitial
       message='M_dark,initial                = '//label//' M☉'
       call displayMessage (message)
       write (label,'(e12.6)') self%baryonicFinalTerm
       message='M_baryonic,final              = '//label//' M☉'
       call displayMessage (message)
       write (label,'(e12.6)') self%baryonicFinalTermDerivative
       message='dM_braryonic,final/dr_initial = '//label//' M☉/Mpc'
       call displayMessage (message)
       call displayUnindent(''     )
       call Error_Report('Overflow in initial radius derivative calculation'//{introspection:location})
    end if
    radiusInitialDerivative=+numerator   &
         &                  /denominator
    return
  end function sphericalAdiabaticGnedin2004RadiusInitialDerivative

  subroutine sphericalAdiabaticGnedin2004ComputeFactors(self,radius,computeGradientFactors)
    !!{
    Compute various factors needed when solving for the initial radius in the dark matter halo using the adiabatic contraction
    algorithm of \cite{gnedin_response_2004}.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionSphericalAdiabaticGnedin2004), intent(inout), target  :: self
    double precision                                              , intent(in   )          :: radius
    logical                                                       , intent(in   )          :: computeGradientFactors
    double precision                                                                       :: velocityCircularSquaredGradient, velocityCircularSquared

    ! Set module-scope pointer to self.
    self_ => self
    ! Store the final radius and its orbit-averaged mean.
    self%radiusFinal    =                       radius
    self%radiusFinalMean=self%radiusOrbitalMean(radius)
    ! Compute the baryonic contribution to the rotation curve.
    velocityCircularSquared=self%massDistributionBaryonic%rotationCurve(self%radiusFinalMean)**2
    self%baryonicFinalTerm=velocityCircularSquared*self%radiusFinalMean*self%radiusFinal/gravitationalConstant_internal
    ! Compute the baryonic contribution to the rotation curve.
    if (computeGradientFactors) then
       velocityCircularSquaredGradient =+self%massDistributionBaryonic%rotationCurveGradient(self%radiusFinalMean)
       self%baryonicFinalTermDerivative=+     velocityCircularSquaredGradient                                      &
            &                           *self%radiusOrbitalMeanDerivative(self%radiusFinal)                        &
            &                           *self%radiusFinalMean                                                      &
            &                           *self%radiusFinal                                                          &
            &                           /     gravitationalConstant_internal
    end if
    return
  end subroutine sphericalAdiabaticGnedin2004ComputeFactors

  double precision function sphericalAdiabaticGnedin2004RadiusOrbitalMean(self,radius)
    !!{
    Returns the orbit averaged radius for dark matter corresponding the given {\normalfont \ttfamily radius} using the model of
    \cite{gnedin_response_2004}.
    !!}
    implicit none
    class           (massDistributionSphericalAdiabaticGnedin2004), intent(inout) :: self
    double precision                                              , intent(in   ) :: radius

    if (self%omega /= omegaPrevious) then
       if (allocated(radiusExponentiator)) deallocate(radiusExponentiator)
       allocate(radiusExponentiator)
       radiusExponentiator=fastExponentiator(1.0d-3,1.0d0,self%omega,1.0d4,.false.)
       omegaPrevious      =                               self%omega
    end if
    sphericalAdiabaticGnedin2004RadiusOrbitalMean=+self               %A                                                 &
         &                                        *self               %radiusFractionalPivot                             &
         &                                        *self               %radiusVirial                                      &
         &                                        *radiusExponentiator%exponentiate         (                            &
         &                                                                                   +     radius                &
         &                                                                                   /self%radiusFractionalPivot &
         &                                                                                   /self%radiusVirial          &
         &                                                                                  )
    return
  end function sphericalAdiabaticGnedin2004RadiusOrbitalMean

  double precision function sphericalAdiabaticGnedin2004RadiusOrbitalMeanDerivative(self,radius)
    !!{
    Returns the derivative of the orbit averaged radius for dark matter corresponding the given {\normalfont \ttfamily radius} using the model of
    \cite{gnedin_response_2004}.
    !!}
    implicit none
    class           (massDistributionSphericalAdiabaticGnedin2004), intent(inout) :: self
    double precision                                              , intent(in   ) :: radius

    sphericalAdiabaticGnedin2004RadiusOrbitalMeanDerivative=+self%A                       &
         &                                                  *self%omega                   &
         &                                                  *(                            &
         &                                                    +     radius                &
         &                                                    /self%radiusFractionalPivot &
         &                                                    /self%radiusVirial          &
         &                                                   )**(self%omega-1.0d0)
    return
  end function sphericalAdiabaticGnedin2004RadiusOrbitalMeanDerivative

  double precision function sphericalAdiabaticGnedin2004Solver(radiusInitial)
    !!{
    Root function used in finding the initial radius in the dark matter halo when solving for adiabatic contraction.
    !!}
    implicit none
    double precision, intent(in   ) :: radiusInitial
    double precision                :: massDarkMatterInitial, radiusInitialMean

    ! Find the initial mean orbital radius.
    radiusInitialMean    =self_                      %radiusOrbitalMean(radiusInitial    )
    ! Get the mass of dark matter inside the initial radius.
    massDarkMatterInitial=self_%massDistribution_%massEnclosedBySphere(radiusInitialMean)
    ! Compute the root function.
    sphericalAdiabaticGnedin2004Solver=+massDarkMatterInitial                                      &
         &                             *(                                                          &
         &                               +self_%massFractionInitial*                 radiusInitial &
         &                               -self_%darkMatterDistributedFraction *self_%radiusFinal   &
         &                              )                                                          &
         &                             -self_%baryonicFinalTerm
    return
  end function sphericalAdiabaticGnedin2004Solver
