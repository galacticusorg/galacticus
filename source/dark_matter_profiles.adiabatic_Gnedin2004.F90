!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  use :: Cosmology_Parameters      , only : cosmologyParameters              , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale              , darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMO             , darkMatterProfileDMOClass
  use :: Mass_Distributions        , only : enumerationNonAnalyticSolversType

  !![
  <darkMatterProfile name="darkMatterProfileAdiabaticGnedin2004">
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
  </darkMatterProfile>
  !!]
  type, extends(darkMatterProfileClass) :: darkMatterProfileAdiabaticGnedin2004
     !!{
     A dark matter halo profile class implementing adiabaticGnedin2004 dark matter halos.
     !!}
     private
     class           (cosmologyParametersClass            ), pointer  :: cosmologyParameters_   => null()
     class           (darkMatterProfileDMOClass           ), pointer  :: darkMatterProfileDMO_  => null()
     class           (darkMatterHaloScaleClass            ), pointer  :: darkMatterHaloScale_   => null()
     integer         (kind_int8                           )           :: lastUniqueID
     type            (enumerationNonAnalyticSolversType   )           :: nonAnalyticSolver
     double precision                                                 :: A                               , omega               , &
          &                                                              radiusFractionalPivot           , toleranceRelative   , &
          &                                                              darkMatterFraction              , massBaryonicSubhalos
   contains    
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset"/>
     </methods>
     !!]
     final     ::                     adiabaticGnedin2004Destructor
     procedure :: autoHook         => adiabaticGnedin2004AutoHook
     procedure :: calculationReset => adiabaticGnedin2004CalculationReset
     procedure :: get              => adiabaticGnedin2004Get
     procedure :: initialize       => adiabaticGnedin2004Initialize
  end type darkMatterProfileAdiabaticGnedin2004

  interface darkMatterProfileAdiabaticGnedin2004
     !!{
     Constructors for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter halo profile class.
     !!}
     module procedure adiabaticGnedin2004ConstructorParameters
     module procedure adiabaticGnedin2004ConstructorInternal
  end interface darkMatterProfileAdiabaticGnedin2004

contains

  function adiabaticGnedin2004ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter halo profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    use :: Input_Parameters  , only : inputParameters
    implicit none
    type            (darkMatterProfileAdiabaticGnedin2004)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (cosmologyParametersClass            ), pointer       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass            ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass           ), pointer       :: darkMatterProfileDMO_
    type            (varying_string                      )                :: nonAnalyticSolver
    double precision                                                      :: A                    , omega            , &
          &                                                                  radiusFractionalPivot, toleranceRelative
    
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
      <name>toleranceRelative</name>
      <defaultValue>1.0d-2</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in solving for the initial radius in the adiabatically-contracted dark matter profile.</description>
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
    self=darkMatterProfileAdiabaticGnedin2004(A,omega,radiusFractionalPivot,toleranceRelative,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_" />
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function adiabaticGnedin2004ConstructorParameters

  function adiabaticGnedin2004ConstructorInternal(A,omega,radiusFractionalPivot,toleranceRelative,nonAnalyticSolver,cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversIsValid
    use :: Error             , only : Error_Report
    implicit none
    type            (darkMatterProfileAdiabaticGnedin2004)                        :: self
    double precision                                      , intent(in   )         :: A                    , omega            , &
         &                                                                           radiusFractionalPivot, toleranceRelative
    class           (cosmologyParametersClass            ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterProfileDMOClass           ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), intent(in   ), target :: darkMatterHaloScale_
    type            (enumerationNonAnalyticSolversType   ), intent(in   )         :: nonAnalyticSolver
    !![
    <constructorAssign variables="A, omega, radiusFractionalPivot, toleranceRelative, nonAnalyticSolver, *cosmologyParameters_, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]
    
    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    ! Evaluate the dark matter fraction.
    self%darkMatterFraction=+1.0d0                                   &
         &                  -self%cosmologyParameters_%OmegaBaryon() &
         &                  /self%cosmologyParameters_%OmegaMatter()
    ! Initialize memoization state.
    self%lastUniqueID        =-1_kind_int8
    self%massBaryonicSubhalos=-huge(0.0d0)
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
    implicit none
    type(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine adiabaticGnedin2004Destructor

  subroutine adiabaticGnedin2004CalculationReset(self,node,uniqueID)
    !!{
    Reset the dark matter profile calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type   (treeNode                            ), intent(inout) :: node
    integer(kind_int8                           ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    if (uniqueID /= self%lastUniqueID) then
       self%lastUniqueID        =uniqueID
       self%massBaryonicSubhalos=-huge(0.0d0)
    end if
    return
  end subroutine adiabaticGnedin2004CalculationReset

  function adiabaticGnedin2004Get(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo                       , massTypeDark                       , massTypeBaryonic         , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalAdiabaticGnedin2004, kinematicsDistributionCollisionless, massDistributionSpherical
    implicit none
    class           (massDistributionClass               ), pointer                 :: massDistribution_
    type            (kinematicsDistributionCollisionless ), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout)           :: self
    type            (treeNode                            ), intent(inout), target   :: node
    type            (enumerationWeightByType             ), intent(in   ), optional :: weightBy
    integer                                               , intent(in   ), optional :: weightIndex
    class           (massDistributionClass               ), pointer                 :: massDistributionDecorated    , massDistributionBaryonic
    double precision                                                                :: massBaryonicSelfTotal        , massBaryonicTotal       , &
         &                                                                             darkMatterDistributedFraction, initialMassFraction
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Set the baryonic component to zero - we will compute this later during initialization.
    massDistributionBaryonic      => null()
    massBaryonicTotal             =  0.0d0
    massBaryonicSelfTotal         =  0.0d0
    darkMatterDistributedFraction =  0.0d0
    initialMassFraction           =  0.0d0
    ! Create the mass distribution.
    allocate(massDistributionSphericalAdiabaticGnedin2004 :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalAdiabaticGnedin2004)
       massDistributionDecorated => self%darkMatterProfileDMO_%get(node,weightBy,weightIndex)
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalAdiabaticGnedin2004(                                                                                             &amp;
	      &amp;                                        A                            =self                     %A                                  , &amp;
	      &amp;                                        omega                        =self                     %omega                              , &amp;
	      &amp;                                        radiusVirial                 =self%darkMatterHaloScale_%radiusVirial                 (node), &amp;
	      &amp;                                        radiusFractionalPivot        =self                     %radiusFractionalPivot              , &amp;
	      &amp;                                        darkMatterFraction           =self                     %darkMatterFraction                 , &amp;
	      &amp;                                        darkMatterDistributedFraction=                          darkMatterDistributedFraction      , &amp;
	      &amp;                                        massFractionInitial          =                          initialMassFraction                , &amp;
	      &amp;                                        nonAnalyticSolver            =self                     %nonAnalyticSolver                  , &amp;
	      &amp;                                        toleranceRelative            =self                     %toleranceRelative                  , &amp;
	      &amp;                                        massDistribution_            =                          massDistributionDecorated          , &amp;
	      &amp;                                        massDistributionBaryonic     =                          massDistributionBaryonic           , &amp;
              &amp;                                        componentType                =                          componentTypeDarkHalo              , &amp;
              &amp;                                        massType                     =                          massTypeDark                         &amp;
              &amp;                                       )
	    </constructor>
          </referenceConstruct>
          <objectDestructor name="massDistributionDecorated"/>
          !!]
       class default
          call Error_Report('expected a spherical mass distribution'//{introspection:location})
       end select
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionCollisionless( &amp;
	 &amp;                             )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function adiabaticGnedin2004Get

  subroutine adiabaticGnedin2004Initialize(self,node,massDistribution_)
    !!{
    Initialize the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : massTypeBaryonic
    use :: Mass_Distributions        , only : massDistributionSphericalAdiabaticGnedin2004
    use :: Error                     , only : Error_Report
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout)         :: self
    type            (treeNode                            ), intent(inout), target :: node
    class           (massDistributionClass               ), intent(inout)         :: massDistribution_
    class           (massDistributionClass               ), pointer               :: massDistributionBaryonic
    type            (treeNode                            ), pointer               :: nodeCurrent
    class           (nodeComponentBasic                  ), pointer               :: basic
    double precision                                                              :: massBaryonicSelf             , massBaryonicTotal  , &
         &                                                                           darkMatterDistributedFraction, initialMassFraction, &
         &                                                                           massBaryonicSubhalos
    
    select type (massDistribution_)
    type is (massDistributionSphericalAdiabaticGnedin2004)
       ! Compute the initial baryonic contribution from this halo, and any satellites.
       massDistributionBaryonic => node%massDistribution(massType=massTypeBaryonic)
       massBaryonicSelf         =  node%massBaryonic    (                         )
       ! Recompute baryonic mass in subhalos if necessary.
       if (self%massBaryonicSubhalos < 0.0d0 .or. node%uniqueID() /= self%lastUniqueID) then
          massBaryonicSubhalos =  0.0d0
          nodeCurrent          => node
          do while (associated(nodeCurrent%firstSatellite))
             nodeCurrent => nodeCurrent%firstSatellite
          end do
          if (associated(nodeCurrent,node)) nodeCurrent => null()
          do while (associated(nodeCurrent))
             massBaryonicSubhalos=+massBaryonicSubhalos                &
                  &               +nodeCurrent          %massBaryonic()
             if (associated(nodeCurrent%sibling)) then
                nodeCurrent => nodeCurrent%sibling
                do while (associated(nodeCurrent%firstSatellite))
                   nodeCurrent => nodeCurrent%firstSatellite
                end do
             else
                nodeCurrent => nodeCurrent%parent
                if (associated(nodeCurrent,node)) nodeCurrent => null()
             end if
          end do
       end if
       ! Compute the total baryonic mass.
       massBaryonicTotal=+massBaryonicSubhalos &
            &            +massBaryonicSelf
              ! Limit masses to physical values.
       massBaryonicSelf =max(massBaryonicSelf ,0.0d0)
       massBaryonicTotal=max(massBaryonicTotal,0.0d0)
       ! Compute the fraction of matter assumed to be distributed like the dark matter.
       basic                         => node%basic()
       darkMatterDistributedFraction =  min(self%darkMatterFraction+(massBaryonicTotal-massBaryonicSelf)/basic%mass(),1.0d0)
       ! Compute the initial mass fraction.
       initialMassFraction           =  min(self%darkMatterFraction+ massBaryonicTotal                  /basic%mass(),1.0d0)
       ! Set the baryonic component in the mass distribution.
       call massDistribution_%setBaryonicComponent(massDistributionBaryonic,darkMatterDistributedFraction,initialMassFraction)
       !![
       <objectDestructor name="massDistributionBaryonic"/>
       !!]
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine adiabaticGnedin2004Initialize
