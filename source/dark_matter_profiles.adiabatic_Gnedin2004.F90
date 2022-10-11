!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  An implementation of adiabaticGnedin2004 dark matter halo profiles.
  !!}

  use :: Cosmology_Parameters        , only : cosmologyParameters              , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScale              , darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO    , only : darkMatterProfileDMO             , darkMatterProfileDMOClass
  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversType, enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough
  use :: Galactic_Structure_Options  , only : componentTypeAll                 , massTypeBaryonic                   , radiusLarge                         , weightByMass                 , &
          &                                   weightIndexNull                  , enumerationComponentTypeType       , enumerationMassTypeType             , enumerationWeightByType
  use :: Math_Exponentiation         , only : fastExponentiator
  use :: Root_Finder                 , only : rootFinder

  ! Number of previous radius solutions to store.
  !![
  <scoping>
   <module variables="adiabaticGnedin2004StoreCount"/>
  </scoping>
  !!]
  integer, parameter :: adiabaticGnedin2004StoreCount=10

  !![
  <darkMatterProfile name="darkMatterProfileAdiabaticGnedin2004" recursive="yes">
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
   <deepCopy>
     <ignore variables="recursiveSelf"/>
   </deepCopy>
   <stateStore>
     <stateStore variables="galacticStructure_" store="galacticStructureStateStore_" restore="galacticStructureStateRestore_" module="Functions_Global"/>
   </stateStore>
  </darkMatterProfile>
  !!]
  type, extends(darkMatterProfileClass) :: darkMatterProfileAdiabaticGnedin2004
     !!{
     A dark matter halo profile class implementing adiabaticGnedin2004 dark matter halos.
     !!}
     private
     logical                                                                                          :: isRecursive                              , parentDeferred
     class           (darkMatterProfileAdiabaticGnedin2004), pointer                                  :: recursiveSelf                   => null()
     class           (cosmologyParametersClass            ), pointer                                  :: cosmologyParameters_            => null()
     class           (darkMatterProfileDMOClass           ), pointer                                  :: darkMatterProfileDMO_           => null()
     class           (*                                   ), pointer                                  :: galacticStructure_              => null()
     type            (enumerationNonAnalyticSolversType   )                                           :: nonAnalyticSolver
     type            (rootFinder                          )                                           :: finder
     ! Parameters of the adiabatic contraction algorithm.
     double precision                                                                                 :: A                                        , omega                     , &
          &                                                                                              radiusFractionalPivot
     ! Stored solutions for reuse.
     integer         (kind=kind_int8                      )                                           :: lastUniqueID
     integer                                                                                          :: radiusPreviousIndex                      , radiusPreviousIndexMaximum
     double precision                                      , dimension(adiabaticGnedin2004StoreCount) :: radiusPrevious                           , radiusInitialPrevious
     type            (fastExponentiator                   )                                           :: radiusExponentiator
     ! Quantities used in solving the initial radius root function.
     double precision                                                                                 :: baryonicFinalTerm                        , baryonicFinalTermDerivative, &
          &                                                                                              darkMatterDistributedFraction            , initialMassFraction        , &
          &                                                                                              radiusFinal                              , radiusFinalMean            , &
          &                                                                                              darkMatterFraction                       , radiusVirial
     logical                                                                                          :: massesComputed
   contains
     !![
     <methods>
       <method description="Reset memoized calculations."                                             method="calculationReset"           />
       <method description="Compute factors needed for solving adiabatic contraction."                method="computeFactors"             />
       <method description="Compute the orbit-averaged radius for dark matter."                       method="radiusOrbitalMean"          />
       <method description="Compute the derivative of the orbit-averaged radius for dark matter."     method="radiusOrbitalMeanDerivative"/>
       <method description="Compute the initial radius in the dark matter profile."                   method="radiusInitial"              />
       <method description="Compute the derivative of the initial radius in the dark matter profile." method="radiusInitialDerivative"    />
     </methods>
     !!]
     final     ::                                      adiabaticGnedin2004Destructor
     procedure :: autoHook                          => adiabaticGnedin2004AutoHook
     procedure :: calculationReset                  => adiabaticGnedin2004CalculationReset
     procedure :: density                           => adiabaticGnedin2004Density
     procedure :: densityLogSlope                   => adiabaticGnedin2004DensityLogSlope
     procedure :: radiusEnclosingDensity            => adiabaticGnedin2004RadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => adiabaticGnedin2004RadiusEnclosingMass
     procedure :: radialMoment                      => adiabaticGnedin2004RadialMoment
     procedure :: enclosedMass                      => adiabaticGnedin2004EnclosedMass
     procedure :: potential                         => adiabaticGnedin2004Potential
     procedure :: circularVelocity                  => adiabaticGnedin2004CircularVelocity
     procedure :: circularVelocityMaximum           => adiabaticGnedin2004CircularVelocityMaximum
     procedure :: radialVelocityDispersion          => adiabaticGnedin2004RadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => adiabaticGnedin2004RadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => adiabaticGnedin2004RotationNormalization
     procedure :: energy                            => adiabaticGnedin2004Energy
     procedure :: kSpace                            => adiabaticGnedin2004KSpace
     procedure :: freefallRadius                    => adiabaticGnedin2004FreefallRadius
     procedure :: freefallRadiusIncreaseRate        => adiabaticGnedin2004FreefallRadiusIncreaseRate
     procedure :: radiusInitial                     => adiabaticGnedin2004RadiusInitial
     procedure :: radiusInitialDerivative           => adiabaticGnedin2004RadiusInitialDerivative
     procedure :: computeFactors                    => adiabaticGnedin2004ComputeFactors
     procedure :: radiusOrbitalMean                 => adiabaticGnedin2004RadiusOrbitalMean
     procedure :: radiusOrbitalMeanDerivative       => adiabaticGnedin2004RadiusOrbitalMeanDerivative
     procedure :: deepCopy                          => adiabaticGnedin2004DeepCopy
     procedure :: deepCopyReset                     => adiabaticGnedin2004DeepCopyReset
     procedure :: deepCopyFinalize                  => adiabaticGnedin2004DeepCopyFinalize
  end type darkMatterProfileAdiabaticGnedin2004

  interface darkMatterProfileAdiabaticGnedin2004
     !!{
     Constructors for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter halo profile class.
     !!}
     module procedure adiabaticGnedin2004ConstructorParameters
     module procedure adiabaticGnedin2004ConstructorInternal
  end interface darkMatterProfileAdiabaticGnedin2004
    
  ! Module-scope quantities used in solving the initial radius root function.
  type             (enumerationComponentTypeType       ), parameter :: adiabaticGnedin2004ComponentType=componentTypeAll
  type             (enumerationMassTypeType            ), parameter :: adiabaticGnedin2004MassType     =massTypeBaryonic
  type             (enumerationWeightByType            ), parameter :: adiabaticGnedin2004WeightBy     =weightByMass
  integer                                               , parameter :: adiabaticGnedin2004WeightIndex  =weightIndexNull
  double precision                                      , parameter :: toleranceAbsolute               =0.0d0           , toleranceRelative             =1.0d-2
  type            (treeNode                            ), pointer   :: adiabaticGnedin2004Node
  class           (darkMatterProfileAdiabaticGnedin2004), pointer   :: adiabaticGnedin2004Self
  !$omp threadprivate(adiabaticGnedin2004Self,adiabaticGnedin2004Node)
    
contains

  recursive function adiabaticGnedin2004ConstructorParameters(parameters,recursiveConstruct,recursiveSelf) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Functions_Global, only : galacticStructureConstruct_, galacticStructureDestruct_
    implicit none
    type            (darkMatterProfileAdiabaticGnedin2004)                          :: self
    type            (inputParameters                     ), intent(inout)           :: parameters
    logical                                               , intent(in   ), optional :: recursiveConstruct
    class           (darkMatterProfileClass              ), intent(in   ), optional :: recursiveSelf
    class           (cosmologyParametersClass            ), pointer                 :: cosmologyParameters_
    class           (darkMatterHaloScaleClass            ), pointer                 :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass           ), pointer                 :: darkMatterProfileDMO_
    class           (*                                   ), pointer                 :: galacticStructure_
    type            (varying_string                      )                          :: nonAnalyticSolver
    double precision                                                                :: A                    , omega, &
          &                                                                            radiusFractionalPivot
    !![
    <optionalArgument name="recursiveConstruct" defaultsTo=".false." />
    !!]
    
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
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring adiabatic contraction by baryons is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    if (recursiveConstruct_) then
       galacticStructure_ => null()
    else
       call galacticStructureConstruct_(parameters,galacticStructure_)
    end if
    self=darkMatterProfileAdiabaticGnedin2004(A,omega,radiusFractionalPivot,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_,galacticStructure_,recursiveConstruct,recursiveSelf)
    !![
    <inputParametersValidate source="parameters" extraAllowedNames="galacticStructure"/>
    <objectDestructor name="cosmologyParameters_" />
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    if (associated(galacticStructure_)) call galacticStructureDestruct_(galacticStructure_)
    return
  end function adiabaticGnedin2004ConstructorParameters

  function adiabaticGnedin2004ConstructorInternal(A,omega,radiusFractionalPivot,nonAnalyticSolver,cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_,galacticStructure_,recursiveConstruct,recursiveSelf) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter profile class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (darkMatterProfileAdiabaticGnedin2004)                                  :: self
    double precision                                      , intent(in   )                   :: A                    , omega, &
         &                                                                                     radiusFractionalPivot
    class           (cosmologyParametersClass            ), intent(in   ), target           :: cosmologyParameters_
    class           (darkMatterProfileDMOClass           ), intent(in   ), target           :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), intent(in   ), target           :: darkMatterHaloScale_
    class           (*                                   ), intent(in   ), target           :: galacticStructure_
    type            (enumerationNonAnalyticSolversType   ), intent(in   )                   :: nonAnalyticSolver
    logical                                               , intent(in   )        , optional :: recursiveConstruct
    class           (darkMatterProfileClass              ), intent(in   ), target, optional :: recursiveSelf
    !![
    <optionalArgument name="recursiveConstruct" defaultsTo=".false." />
    <constructorAssign variables="A, omega, radiusFractionalPivot, nonAnalyticSolver, *cosmologyParameters_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *galacticStructure_"/>
    !!]
    
    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    ! Construct the object.
    self%lastUniqueID       =-1_kind_int8
    self%genericLastUniqueID=-1_kind_int8
    self%massesComputed     =.false.
    self%radiusExponentiator=fastExponentiator(1.0d-3,1.0d0,omega,1.0d4,.false.)
    ! Evaluate the dark matter fraction.
    self%darkMatterFraction=+1.0d0                                   &
         &                  -self%cosmologyParameters_%OmegaBaryon() &
         &                  /self%cosmologyParameters_%OmegaMatter()
    ! Construct a root finder.
    self%finder=rootFinder(                                             &
         &                 rootFunction     =adiabaticGnedin2004Solver, &
         &                 toleranceAbsolute=toleranceAbsolute        , &
         &                 toleranceRelative=toleranceRelative          &
         &                )
    ! Handle recursive construction.
    self%isRecursive=recursiveConstruct_
    if (recursiveConstruct_) then
       if (.not.present(recursiveSelf)) call Error_Report('recursiveSelf not present'//{introspection:location})
       select type (recursiveSelf)
       class is (darkMatterProfileAdiabaticGnedin2004)
          self%recursiveSelf => recursiveSelf
       class default
          call Error_Report('recursiveSelf is of incorrect class'//{introspection:location})
       end select
    end if
    self%parentDeferred=.false.
    return
  end function adiabaticGnedin2004ConstructorInternal

  subroutine adiabaticGnedin2004AutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self

    call calculationResetEvent%attach(self,adiabaticGnedin2004CalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileAdiabaticGnedin2004')
    return
  end subroutine adiabaticGnedin2004AutoHook

  subroutine adiabaticGnedin2004Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter halo profile class.
    !!}
    use :: Events_Hooks    , only : calculationResetEvent
    use :: Functions_Global, only : galacticStructureDestruct_
    implicit none
    type(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    if (associated                      (self%galacticStructure_                 )) call galacticStructureDestruct_  (self%galacticStructure_                 )
    if (calculationResetEvent%isAttached(self,adiabaticGnedin2004CalculationReset)) call calculationResetEvent%detach(self,adiabaticGnedin2004CalculationReset)
    return
  end subroutine adiabaticGnedin2004Destructor

  subroutine adiabaticGnedin2004CalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    ! Reset calculations for this profile.
    self%lastUniqueID                                =node%uniqueID()
    self%genericLastUniqueID                         =node%uniqueID()
    self%radiusPreviousIndex                         = 0
    self%radiusPreviousIndexMaximum                  = 0
    self%radiusPrevious                              =-1.0d0
    self%massesComputed                              =.false.
    self%genericEnclosedMassRadiusMinimum            =+huge(0.0d0)
    self%genericEnclosedMassRadiusMaximum            =-huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMinimum=+huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMaximum=-huge(0.0d0)
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine adiabaticGnedin2004CalculationReset

  double precision function adiabaticGnedin2004EnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004EnclosedMass=self%recursiveSelf%enclosedMass(node,radius)
       return
    end if
    adiabaticGnedin2004EnclosedMass=+self%darkMatterFraction                                                       &
         &                          *self%darkMatterProfileDMO_%enclosedMass(node,self%radiusInitial(node,radius))
    return
  end function adiabaticGnedin2004EnclosedMass

  double precision function adiabaticGnedin2004Density(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    double precision                                                      :: radiusInitial , radiusInitialDerivative, &
         &                                                                   densityInitial

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004Density=self%recursiveSelf%density(node,radius)
       return
    end if
    radiusInitial =self                      %radiusInitial(node,radius       )
    densityInitial=self%darkMatterProfileDMO_%density      (node,radiusInitial)
    if (radius == radiusInitial) then
       adiabaticGnedin2004Density=+self%darkMatterFraction &
            &                     *densityInitial
    else
       radiusInitialDerivative   = self%radiusInitialDerivative(node,radius)
       adiabaticGnedin2004Density=+self%darkMatterFraction    &
            &                     *densityInitial             &
            &                     *(                          &
            &                       +radiusInitial            &
            &                       /radius                   &
            &                      )                      **2 &
            &                     *radiusInitialDerivative
    end if
    return
  end function adiabaticGnedin2004Density

  double precision function adiabaticGnedin2004DensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    logical                                                               :: fallThrough

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004DensityLogSlope=self%recursiveSelf%densityLogSlope(node,radius)
       return
    end if
    fallThrough=self%nonAnalyticSolver == nonAnalyticSolversFallThrough
    if (.not.fallThrough) fallThrough=radius == self%radiusInitial(node,radius)
    if (fallThrough) then
       adiabaticGnedin2004DensityLogSlope=self%darkMatterProfileDMO_%densityLogSlope         (node,radius)
    else
       adiabaticGnedin2004DensityLogSlope=self                      %densityLogSlopeNumerical(node,radius)
    end if
    return
  end function adiabaticGnedin2004DensityLogSlope

  double precision function adiabaticGnedin2004RadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: density

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004RadiusEnclosingDensity=self%recursiveSelf%radiusEnclosingDensity(node,density)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity         (node,density)
    else
       adiabaticGnedin2004RadiusEnclosingDensity=self                      %radiusEnclosingDensityNumerical(node,density)
    end if
    return
  end function adiabaticGnedin2004RadiusEnclosingDensity

  double precision function adiabaticGnedin2004RadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: mass

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004RadiusEnclosingMass=self%recursiveSelf%radiusEnclosingDensity(node,mass)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass         (node,mass)
    else
       adiabaticGnedin2004RadiusEnclosingMass=self                      %radiusEnclosingMassNumerical(node,mass)
    end if
    return
  end function adiabaticGnedin2004RadiusEnclosingMass

  double precision function adiabaticGnedin2004RadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the radial moment of the density in the dark matter profile of {\normalfont \ttfamily node} between the given
    {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    double precision                                      , intent(in   )           :: moment
    double precision                                      , intent(in   ), optional :: radiusMinimum, radiusMaximum

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004RadialMoment=self%recursiveSelf%radialMoment(node,moment,radiusMinimum,radiusMaximum)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadialMoment=self%darkMatterProfileDMO_%radialMoment         (node,moment,radiusMinimum,radiusMaximum)
    else
       adiabaticGnedin2004RadialMoment=self                      %radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    end if
    return
  end function adiabaticGnedin2004RadialMoment

  double precision function adiabaticGnedin2004Potential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target   :: self
    type            (treeNode                            ), intent(inout), target   :: node
    double precision                                      , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType   ), intent(  out), optional :: status
    logical                                                                         :: fallThrough

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004Potential=self%recursiveSelf%potential(node,radius,status)
       return
    end if
    fallThrough=self%nonAnalyticSolver == nonAnalyticSolversFallThrough
    if (.not.fallThrough) fallThrough=radius == self%radiusInitial(node,radius)
    if (fallThrough) then
       ! No adiabatic contraction - use the dark-matter-only result.
       adiabaticGnedin2004Potential=+self%darkMatterFraction                                           &
            &                       *self%darkMatterProfileDMO_%potential         (node,radius,status)
    else
       ! Adiabatic contraction is present - fall back to using a numerical calculation.
       adiabaticGnedin2004Potential=+self                      %potentialNumerical(node,radius,status)
    end if
    return
  end function adiabaticGnedin2004Potential

  double precision function adiabaticGnedin2004CircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004CircularVelocity=self%recursiveSelf%circularVelocity(node,radius)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004CircularVelocity=self%darkMatterProfileDMO_%circularVelocity         (node,radius)
    else
       adiabaticGnedin2004CircularVelocity=self                      %circularVelocityNumerical(node,radius)
    end if
    return
  end function adiabaticGnedin2004CircularVelocity

  double precision function adiabaticGnedin2004CircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004CircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum         (node)
    else
       adiabaticGnedin2004CircularVelocityMaximum=self                      %circularVelocityMaximumNumerical(node)
    end if
    return
  end function adiabaticGnedin2004CircularVelocityMaximum

  double precision function adiabaticGnedin2004RadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004RadialVelocityDispersion=self%recursiveSelf%radialVelocityDispersion(node,radius)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadialVelocityDispersion=self%darkMatterProfileDMO_%radialVelocityDispersion         (node,radius)
    else
       adiabaticGnedin2004RadialVelocityDispersion=self                      %radialVelocityDispersionNumerical(node,radius)
    end if
    return
  end function adiabaticGnedin2004RadialVelocityDispersion

  double precision function adiabaticGnedin2004RadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target :: self
    type            (treeNode                            ), intent(inout)         :: node
    double precision                                      , intent(in   )         :: specificAngularMomentum

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004RadiusFromSpecificAngularMomentum=self%recursiveSelf%radiusFromSpecificAngularMomentum(node,specificAngularMomentum)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum         (node,specificAngularMomentum)
    else
       adiabaticGnedin2004RadiusFromSpecificAngularMomentum=self                      %radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    end if
    return
  end function adiabaticGnedin2004RadiusFromSpecificAngularMomentum

  double precision function adiabaticGnedin2004RotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004RotationNormalization=self%recursiveSelf%rotationNormalization(node)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RotationNormalization=self%darkMatterProfileDMO_%rotationNormalization         (node)
    else
       adiabaticGnedin2004RotationNormalization=self                      %rotationNormalizationNumerical(node)
    end if
    return
  end function adiabaticGnedin2004RotationNormalization

  double precision function adiabaticGnedin2004Energy(self,node)
    !!{
    Return the energy of a adiabaticGnedin2004 halo density profile.
    !!}
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004Energy=self%recursiveSelf%energy(node)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004Energy=self%darkMatterProfileDMO_%energy         (node)
    else
       adiabaticGnedin2004Energy=self                      %energyNumerical(node)
    end if
    return
  end function adiabaticGnedin2004Energy

  double precision function adiabaticGnedin2004KSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the adiabaticGnedin2004 density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$), using the expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout)         :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: waveNumber

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004KSpace=self%recursiveSelf%kSpace(node,wavenumber)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004KSpace=self%darkMatterProfileDMO_%kSpace         (node,waveNumber)
    else
       adiabaticGnedin2004KSpace=self                      %kSpaceNumerical(node,waveNumber)
    end if
    return
  end function adiabaticGnedin2004KSpace

  double precision function adiabaticGnedin2004FreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the adiabaticGnedin2004 density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: time

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004FreefallRadius=self%recursiveSelf%freefallRadius(node,time)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004FreefallRadius=self%darkMatterProfileDMO_%freefallRadius         (node,time)
    else
       adiabaticGnedin2004FreefallRadius=self                      %freefallRadiusNumerical(node,time)
    end if
    return
  end function adiabaticGnedin2004FreefallRadius

  double precision function adiabaticGnedin2004FreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the adiabaticGnedin2004 density profile at the specified
    {\normalfont \ttfamily time} (given in Gyr).
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: time

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       adiabaticGnedin2004FreefallRadiusIncreaseRate=self%recursiveSelf%freefallRadiusIncreaseRate(node,time)
       return
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004FreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate         (node,time)
    else
       adiabaticGnedin2004FreefallRadiusIncreaseRate=self                      %freefallRadiusIncreaseRateNumerical(node,time)
    end if
    return
  end function adiabaticGnedin2004FreefallRadiusIncreaseRate

  double precision function adiabaticGnedin2004RadiusInitial(self,node,radius)
    !!{
    Compute the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    \cite{gnedin_response_2004}.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    integer                                                               :: i               , j           , &
         &                                                                   iMod
    double precision                                                      :: radiusUpperBound, massEnclosed

    ! Reset stored solutions if the node has changed.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check for a previously computed solution.
    if (self%radiusPreviousIndexMaximum > 0 .and. any(self%radiusPrevious(1:self%radiusPreviousIndexMaximum) == radius)) then
       adiabaticGnedin2004RadiusInitial=0.0d0
       do i=1,self%radiusPreviousIndexMaximum
          if (self%radiusPrevious(i) == radius) then
             adiabaticGnedin2004RadiusInitial=self%radiusInitialPrevious(i)
             exit
          end if
       end do
       return
    end if
    ! Get the virial radius of the node.
    self%radiusVirial=self%darkMatterHaloScale_%radiusVirial(node)
    ! Return radius unchanged if larger than the virial radius.
    if (radius >= self%radiusVirial) then
       adiabaticGnedin2004RadiusInitial=radius
       return
    end if
    ! Compute the various factors needed by this calculation.
    call self%computeFactors(node,radius,computeGradientFactors=.false.)
    !! Note that even if no baryons are present at this radius we can not assume that the initial radius is unchanged because it
    !! is possible that the initial fraction of baryons and the fraction of mass distributed as the dark matter are not equal, fᵢ
    !! ≠ fᵪ.
    ! Check that solution is within bounds.
    if (adiabaticGnedin2004Solver(self%radiusVirial) < 0.0d0) then
       adiabaticGnedin2004RadiusInitial=self%radiusVirial
       return
    end if
    j=-1
    if (self%radiusPreviousIndexMaximum > 0) then
       ! No exact match exists, look for approximate matches.
       do i=1,self%radiusPreviousIndexMaximum
          iMod=modulo(self%radiusPreviousIndex-i,adiabaticGnedin2004StoreCount)+1
          if (abs(radius-self%radiusPrevious(iMod))/self%radiusPrevious(iMod) < toleranceRelative) then
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
       massEnclosed=+self%darkMatterProfileDMO_%enclosedMass(node,self%radiusOrbitalMean(self%radiusFinal))
       if (massEnclosed > 0.0d0) then
          radiusUpperBound=+(                                    &
               &             +self%baryonicFinalTerm             &
               &             /     massEnclosed                  &
               &             +self%darkMatterDistributedFraction &
               &             *self%radiusFinal                   &
               &            )                                    &
               &           /  self%initialMassFraction
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
       adiabaticGnedin2004RadiusInitial=self%finder%find(rootRange=[radius,radiusUpperBound])
    else
       ! Use previous solution as an initial guess.
       call self%finder%rangeExpand(                                                                   &
            &                       rangeExpandDownward          =1.0d0/sqrt(1.0d0+toleranceRelative), &
            &                       rangeExpandUpward            =1.0d0*sqrt(1.0d0+toleranceRelative), &
            &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative      , &
            &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive      , &
            &                       rangeExpandType              =rangeExpandMultiplicative            &
            &                      )
       adiabaticGnedin2004RadiusInitial=self%finder%find(                                                                        &
            &                                            rootRange=[                                                             &
            &                                                       self%radiusInitialPrevious(j)/sqrt(1.0d0+toleranceRelative), &
            &                                                       self%radiusInitialPrevious(j)*sqrt(1.0d0+toleranceRelative)  &
            &                                                      ]                                                             &
            &                                           )
    end if
    ! Store this solution.
    self%radiusPreviousIndex                                 =modulo(self%radiusPreviousIndex         ,adiabaticGnedin2004StoreCount)+1
    self%radiusPreviousIndexMaximum                          =min   (self%radiusPreviousIndexMaximum+1,adiabaticGnedin2004StoreCount)
    self%radiusPrevious            (self%radiusPreviousIndex)=radius
    self%radiusInitialPrevious     (self%radiusPreviousIndex)=adiabaticGnedin2004RadiusInitial
    return
  end function adiabaticGnedin2004RadiusInitial

  double precision function adiabaticGnedin2004RadiusInitialDerivative(self,node,radius)
    !!{
    Compute the derivative of the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    \cite{gnedin_response_2004}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    double precision                                                      :: radiusInitial                  , radiusInitialMean            , &
         &                                                                   massDarkMatterInitial          , densityDarkMatterInitial     , &
         &                                                                   radiusInitialMeanSelfDerivative, radiusFinalMeanSelfDerivative

    ! Reset stored solutions if the node has changed.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Compute the various factors needed by this calculation.
    call self%computeFactors(node,radius,computeGradientFactors=.true.)
    ! Return unit derivative if radius is larger than the virial radius.
    if (radius >= self%radiusVirial) then
       adiabaticGnedin2004RadiusInitialDerivative=1.0d0
       return
    end if
    ! Compute initial radius, and derivatives of initial and final mean radii.
    radiusInitial                  =self                      %radiusInitial              (node,radius           )
    radiusInitialMeanSelfDerivative=self                      %radiusOrbitalMeanDerivative(     radiusInitial    )
    radiusFinalMeanSelfDerivative  =self                      %radiusOrbitalMeanDerivative(     radius           )
    ! Find the initial mean orbital radius.
    radiusInitialMean              =self                      %radiusOrbitalMean          (     radiusInitial    )
    ! Get the mass of dark matter inside the initial radius.
    massDarkMatterInitial          =self%darkMatterProfileDMO_%enclosedMass               (node,radiusInitialMean)
    ! Get the mass of dark matter inside the initial radius.
    densityDarkMatterInitial       =self%darkMatterProfileDMO_%density                    (node,radiusInitialMean)
    ! Find the solution for the derivtive of the initial radius.
    adiabaticGnedin2004RadiusInitialDerivative=+(                                                      &
         &                                       +massDarkMatterInitial                                &
         &                                       *self%darkMatterDistributedFraction                   &
         &                                       +self%baryonicFinalTerm                               &
         &                                       *(                                                    &
         &                                         +1.0d0                        /     radius          &
         &                                         +radiusFinalMeanSelfDerivative/self%radiusFinalMean &
         &                                        )                                                    &
         &                                       +self%baryonicFinalTermDerivative                     &
         &                                      )                                                      &
         &                                     /(                                                      &
         &                                       +massDarkMatterInitial*self%initialMassFraction       &
         &                                       +(                                                    &
         &                                         +self%initialMassFraction          *radiusInitial   &
         &                                         -self%darkMatterDistributedFraction*radius          &
         &                                        )                                                    &
         &                                       *4.0d0                                                &
         &                                       *Pi                                                   &
         &                                       *radiusInitialMean**2                                 &
         &                                       *densityDarkMatterInitial                             &
         &                                       *radiusInitialMeanSelfDerivative                      &
         &                                      )
    return
  end function adiabaticGnedin2004RadiusInitialDerivative

  subroutine adiabaticGnedin2004ComputeFactors(self,node,radius,computeGradientFactors)
    !!{
    Compute various factors needed when solving for the initial radius in the dark matter halo using the adiabatic contraction
    algorithm of \cite{gnedin_response_2004}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic             , optimizeForEnclosedMassSummation  , optimizeForRotationCurveGradientSummation , optimizeForRotationCurveSummation, &
          &                                         reductionSummation             , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Functions_Global                , only : galacticStructureMassEnclosed_ , galacticStructureVelocityRotation_, galacticStructureVelocityRotationGradient_
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target  :: self
    type            (treeNode                            ), intent(inout), target  :: node
    double precision                                      , intent(in   )          :: radius
    logical                                               , intent(in   )          :: computeGradientFactors
    type            (treeNode                            )               , pointer :: nodeCurrent                    , nodeHost
    class           (nodeComponentBasic                  )               , pointer :: basic
    double precision                                                               :: massBaryonicSelfTotal          , massBaryonicTotal      , &
         &                                                                            velocityCircularSquaredGradient, velocityCircularSquared, &
         &                                                                            

    ! Set module-scope pointers to node and self.
    adiabaticGnedin2004Node => node
    adiabaticGnedin2004Self => self
    ! Store the final radius and its orbit-averaged mean.
    self%radiusFinal    =                       radius
    self%radiusFinalMean=self%radiusOrbitalMean(radius)
    ! Compute the baryonic contribution to the rotation curve.
    velocityCircularSquared=galacticStructureVelocityRotation_(self%galacticStructure_,node,self%radiusFinalMean,adiabaticGnedin2004ComponentType,adiabaticGnedin2004MassType)**2
    self%baryonicFinalTerm=velocityCircularSquared*self%radiusFinalMean*self%radiusFinal/gravitationalConstantGalacticus
    ! Compute the baryonic contribution to the rotation curve.
    if (computeGradientFactors) then
       velocityCircularSquaredGradient =+galacticStructureVelocityRotationGradient_(self%galacticStructure_,node,self%radiusFinalMean,adiabaticGnedin2004ComponentType,adiabaticGnedin2004MassType) &
            &                           *2.0d0                                                                                                                                                      &
            &                           *sqrt(velocityCircularSquared)
       self%baryonicFinalTermDerivative=+     velocityCircularSquaredGradient                                                                                                                       &
            &                           *self%radiusOrbitalMeanDerivative(self%radiusFinal)                                                                                                         &
            &                           *self%radiusFinalMean                                                                                                                                       &
            &                           *self%radiusFinal                                                                                                                                           &
            &                           /     gravitationalConstantGalacticus
    end if
    ! Compute the initial baryonic contribution from this halo, and any satellites.
    if (.not.self%massesComputed) then
       massBaryonicTotal     =  0.0d0
       massBaryonicSelfTotal =  0.0d0
       nodeCurrent           => node
       nodeHost              => node
       do while (associated(nodeCurrent))
          massBaryonicTotal=+massBaryonicTotal                                                &
               &            +galacticStructureMassEnclosed_(                                  &
               &                                            self%galacticStructure_         , &
               &                                            nodeCurrent                     , &
               &                                            radiusLarge                     , &
               &                                            adiabaticGnedin2004ComponentType, &
               &                                            adiabaticGnedin2004MassType     , &
               &                                            adiabaticGnedin2004WeightBy     , &
               &                                            adiabaticGnedin2004WeightIndex    &
               &                                           )
          if (associated(nodeCurrent,nodeHost)) then
             massBaryonicSelfTotal=massBaryonicTotal
             do while (associated(nodeCurrent%firstSatellite))
                nodeCurrent => nodeCurrent%firstSatellite
             end do
             if (associated(nodeCurrent,nodeHost)) nodeCurrent => null()
          else
             if (associated(nodeCurrent%sibling)) then
                nodeCurrent => nodeCurrent%sibling
                do while (associated(nodeCurrent%firstSatellite))
                   nodeCurrent => nodeCurrent%firstSatellite
                end do
             else
                nodeCurrent => nodeCurrent%parent
                if (associated(nodeCurrent,node)) nodeCurrent => null()
             end if
          end if
       end do
       ! Limit masses to physical values.
       massBaryonicSelfTotal=max(massBaryonicSelfTotal,0.0d0)
       massBaryonicTotal    =max(massBaryonicTotal    ,0.0d0)
       ! Compute the fraction of matter assumed to be distributed like the dark matter.
       basic                              => node%basic()
       self%darkMatterDistributedFraction =min((self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon())/self%cosmologyParameters_%OmegaMatter()+(massBaryonicTotal-massBaryonicSelfTotal)/basic%mass(),1.0d0)
       ! Compute the initial mass fraction.
       self%initialMassFraction           =min((self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon())/self%cosmologyParameters_%OmegaMatter()+ massBaryonicTotal                       /basic%mass(),1.0d0)
       ! Record that masses (and mass fractions) have been computed.
       self%massesComputed=.true.
    end if
    return
  end subroutine adiabaticGnedin2004ComputeFactors

  double precision function adiabaticGnedin2004RadiusOrbitalMean(self,radius)
    !!{
    Returns the orbit averaged radius for dark matter corresponding the given {\normalfont \ttfamily radius} using the model of
    \cite{gnedin_response_2004}.
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius

    adiabaticGnedin2004RadiusOrbitalMean=+self%A                                                            &
         &                               *self%radiusFractionalPivot                                        &
         &                               *self%radiusVirial                                                 &
         &                               *self%radiusExponentiator%exponentiate(                            &
         &                                                                      +     radius                &
         &                                                                      /self%radiusFractionalPivot &
         &                                                                      /self%radiusVirial          &
         &                                                                     )
    return
  end function adiabaticGnedin2004RadiusOrbitalMean

  double precision function adiabaticGnedin2004RadiusOrbitalMeanDerivative(self,radius)
    !!{
    Returns the derivative of the orbit averaged radius for dark matter corresponding the given {\normalfont \ttfamily radius} using the model of
    \cite{gnedin_response_2004}.
    !!}
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius

    adiabaticGnedin2004RadiusOrbitalMeanDerivative=+self%A                       &
         &                                         *self%omega                   &
         &                                         *(                            &
         &                                           +     radius                &
         &                                           /self%radiusFractionalPivot &
         &                                           /self%radiusVirial          &
         &                                          )**(self%omega-1.0d0)
    return
  end function adiabaticGnedin2004RadiusOrbitalMeanDerivative

  double precision function adiabaticGnedin2004Solver(radiusInitial)
    !!{
    Root function used in finding the initial radius in the dark matter halo when solving for adiabatic contraction.
    !!}
    implicit none
    double precision, intent(in   ) :: radiusInitial
    double precision                :: massDarkMatterInitial, radiusInitialMean

    ! Find the initial mean orbital radius.
    radiusInitialMean    =adiabaticGnedin2004Self                      %radiusOrbitalMean(                        radiusInitial    )
    ! Get the mass of dark matter inside the initial radius.
    massDarkMatterInitial=adiabaticGnedin2004Self%darkMatterProfileDMO_%enclosedMass     (adiabaticGnedin2004Node,radiusInitialMean)
    ! Compute the root function.
    adiabaticGnedin2004Solver=+massDarkMatterInitial                                                                          &
         &                    *(                                                                                              &
         &                      +adiabaticGnedin2004Self%initialMassFraction*                                   radiusInitial &
         &                      -adiabaticGnedin2004Self%darkMatterDistributedFraction *adiabaticGnedin2004Self%radiusFinal   &
         &                     )                                                                                              &
         &                    -adiabaticGnedin2004Self%baryonicFinalTerm
    return
  end function adiabaticGnedin2004Solver

  subroutine adiabaticGnedin2004DeepCopyReset(self)
    !!{
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : galacticStructureDeepCopyReset_
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    
    self                           %   copiedSelf => null()
    if (.not.self%isRecursive) self%recursiveSelf => null()
    if (associated(self%cosmologyParameters_ )) call self%cosmologyParameters_ %deepCopyReset()
    if (associated(self%darkMatterHaloScale_ )) call self%darkMatterHaloScale_ %deepCopyReset()
    if (associated(self%darkMatterProfileDMO_)) call self%darkMatterProfileDMO_%deepCopyReset()
    if (associated(self%galacticStructure_   )) call galacticStructureDeepCopyReset_(self%galacticStructure_)
    return
  end subroutine adiabaticGnedin2004DeepCopyReset
  
  subroutine adiabaticGnedin2004DeepCopyFinalize(self)
    !!{
    Finalize a deep reset of the object.
    !!}
    use :: Functions_Global, only : galacticStructureDeepCopyFinalize_
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self

    if (self%isRecursive) call adiabaticGnedin2004FindParent(self)
    if (associated(self%cosmologyParameters_ )) call self%cosmologyParameters_ %deepCopyFinalize()
    if (associated(self%darkMatterHaloScale_ )) call self%darkMatterHaloScale_ %deepCopyFinalize()
    if (associated(self%darkMatterProfileDMO_)) call self%darkMatterProfileDMO_%deepCopyFinalize()
    if (associated(self%galacticStructure_   )) call galacticStructureDeepCopyFinalize_(self%galacticStructure_)
    return
  end subroutine adiabaticGnedin2004DeepCopyFinalize
  
  subroutine adiabaticGnedin2004DeepCopy(self,destination)
    !!{
    Perform a deep copy of the object.
    !!}
    use :: Error             , only : Error_Report
    use :: Functions_Global  , only : galacticStructureDeepCopy_
#ifdef OBJECTDEBUG
    use :: Display           , only : displayMessage            , verbosityLevelSilent
    use :: MPI_Utilities     , only : mpiSelf
    use :: Function_Classes  , only : debugReporting
    use :: ISO_Varying_String, only : operator(//)              , var_str
    use :: String_Handling   , only : operator(//)
#endif
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    class(darkMatterProfileClass              ), intent(inout) :: destination

    call self%darkMatterProfileClass%deepCopy(destination)
    select type (destination)
    type is (darkMatterProfileAdiabaticGnedin2004)
       destination%finder                         =self%finder                          
       destination%A                              =self%A                                       
       destination%omega                          =self%omega                           
       destination%radiusFractionalPivot          =self%radiusFractionalPivot                           
       destination%lastUniqueID                   =self%lastUniqueID                    
       destination%radiusPreviousIndex            =self%radiusPreviousIndex                     
       destination%radiusPreviousIndexMaximum     =self%radiusPreviousIndexMaximum      
       destination%radiusPrevious                 =self%radiusPrevious                          
       destination%radiusInitialPrevious          =self%radiusInitialPrevious           
       destination%radiusExponentiator            =self%radiusExponentiator             
       destination%baryonicFinalTerm              =self%baryonicFinalTerm                       
       destination%baryonicFinalTermDerivative    =self%baryonicFinalTermDerivative
       destination%darkMatterDistributedFraction  =self%darkMatterDistributedFraction           
       destination%initialMassFraction            =self%initialMassFraction             
       destination%radiusFinal                    =self%radiusFinal                             
       destination%radiusFinalMean                =self%radiusFinalMean                 
       destination%radiusVirial                   =self%radiusVirial                            
       destination%darkMatterFraction             =self%darkMatterFraction              
       destination%massesComputed                 =self%massesComputed                  
       destination%isRecursive                    =self%isRecursive
       destination%parentDeferred                 =.false.
       if (self%isRecursive) then
          if (associated(self%recursiveSelf%recursiveSelf)) then
             ! If the parent self's recursiveSelf pointer is set, it indicates that it was deep-copied, and the pointer points to
             ! that copy. In that case we set the parent self of our destination to that copy.
             destination%recursiveSelf  => self%recursiveSelf%recursiveSelf
          else
             ! The parent self does not appear to have been deep-copied yet. Retain the same parent self pointer in our copy, but
             ! indicate that we need to look for the new parent later.
             destination%recursiveSelf  => self%recursiveSelf
             destination%parentDeferred =  .true.
          end if
       else
          ! This is a parent of a recursively-constructed object. Record the location of our copy so that it can be used to set
          ! the parent in deep copies of the child object.
          call adiabaticGnedin2004DeepCopyAssign(self,destination)
          destination%recursiveSelf     => null()
       end if
       nullify(destination%cosmologyParameters_)
       if (associated(self%cosmologyParameters_)) then
          if (associated(self%cosmologyParameters_%copiedSelf)) then
             select type(s => self%cosmologyParameters_%copiedSelf)
                class is (cosmologyParametersClass)
                destination%cosmologyParameters_ => s
                class default
                call Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%cosmologyParameters_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%cosmologyParameters_,mold=self%cosmologyParameters_)
             call self%cosmologyParameters_%deepCopy(destination%cosmologyParameters_)
             self%cosmologyParameters_%copiedSelf => destination%cosmologyParameters_
             call destination%cosmologyParameters_%autoHook()
          end if
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): cosmologyparameters : [destination] : ')//loc(destination)//' : '//loc(destination%cosmologyParameters_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
       nullify(destination%darkMatterProfileDMO_)
       if (associated(self%darkMatterProfileDMO_)) then
          if (associated(self%darkMatterProfileDMO_%copiedSelf)) then
             select type(s => self%darkMatterProfileDMO_%copiedSelf)
                class is (darkMatterProfileDMOClass)
                destination%darkMatterProfileDMO_ => s
                class default
                call Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%darkMatterProfileDMO_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%darkMatterProfileDMO_,mold=self%darkMatterProfileDMO_)
             call self%darkMatterProfileDMO_%deepCopy(destination%darkMatterProfileDMO_)
             self%darkMatterProfileDMO_%copiedSelf => destination%darkMatterProfileDMO_
             call destination%darkMatterProfileDMO_%autoHook()
          end if
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): darkmatterprofiledmo_ : [destination] : ')//loc(destination)//' : '//loc(destination%darkMatterProfileDMO_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
       nullify(destination%galacticStructure_)
       if (associated(self%galacticStructure_)) then
          allocate(destination%galacticStructure_,mold=self%galacticStructure_)
          call galacticStructureDeepCopy_(self%galacticStructure_,destination%galacticStructure_)
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): galacticstructure : [destination] : ')//loc(destination)//' : '//loc(destination%galacticStructure_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
       nullify(destination%darkMatterHaloScale_)
       if (associated(self%darkMatterHaloScale_)) then
          if (associated(self%darkMatterHaloScale_%copiedSelf)) then
             select type(s => self%darkMatterHaloScale_%copiedSelf)
                class is (darkMatterHaloScaleClass)
                destination%darkMatterHaloScale_ => s
                class default
                call Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%darkMatterHaloScale_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%darkMatterHaloScale_,mold=self%darkMatterHaloScale_)
             call self%darkMatterHaloScale_%deepCopy(destination%darkMatterHaloScale_)
             self%darkMatterHaloScale_%copiedSelf => destination%darkMatterHaloScale_
             call destination%darkMatterHaloScale_%autoHook()
          end if
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): darkmatterhaloscale : [destination] : ')//loc(destination)//' : '//loc(destination%darkMatterHaloScale_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
       call destination%finder%deepCopyActions()
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine adiabaticGnedin2004DeepCopy

  subroutine adiabaticGnedin2004DeepCopyAssign(self,destination)
    !!{
    Perform pointer assignment during a deep copy of the object.
    !!}
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout)         :: self
    class(darkMatterProfileClass              ), intent(inout), target :: destination

    select type (destination)
    type is (darkMatterProfileAdiabaticGnedin2004)
       self%recursiveSelf => destination
    end select
    return
  end subroutine adiabaticGnedin2004DeepCopyAssign

  subroutine adiabaticGnedin2004FindParent(self)
    !!{
    Find the deep-copied parent of a recursive child.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self

    if (self%parentDeferred) then
       if (associated(self%recursiveSelf%recursiveSelf)) then
          self%recursiveSelf => self%recursiveSelf%recursiveSelf
       else
         call Error_Report("recursive child's parent was not copied"//{introspection:location})
       end if
       self%parentDeferred=.false.
    end if
    return
  end subroutine adiabaticGnedin2004FindParent
