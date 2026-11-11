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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Implements a merger tree branching probability class which conditions a branch on reaching a specified progenitor.
!!}

  use :: Cosmological_Density_Field        , only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions               , only : cosmologyFunctionsClass
  use :: ISO_Varying_String                , only : varying_string
  use :: Linear_Growth                     , only : linearGrowthClass
  use :: Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolutionClass
  use :: Numerical_Integration             , only : integrator
  use :: Root_Finder                       , only : rootFinder

  !![
  <mergerTreeBranchingProbability name="mergerTreeBranchingProbabilityConstraintConditioned" docformat="rst">
   <description>
   A merger tree branching probability class for the branch of a tree which is required to pass through a progenitor of mass
   :math:`M_1` at critical overdensity :math:`\delta_1`. Branching rates are those of a wrapped
   :galacticus-class:`mergerTreeBranchingProbabilityGeneralizedPressSchechter` object, which must use the full range of its
   distribution function, conditioned on the branch reaching the constraint.

   Under the tree process a halo of mass :math:`M` draws a mass :math:`D` at rate :math:`n(D)/2` per unit :math:`\delta`, and
   splits into :math:`D` and :math:`M-D`, where :math:`n` is the rate given by the wrapped object. The expected number,
   :math:`H(M,\delta)`, of progenitors in a narrow window about :math:`M_1` at :math:`\delta_1` satisfies the backward equation

   .. math::

      -{\partial H \over \partial \delta} = -v {\partial H \over \partial M} + \int_{M_\mathrm{res}}^M {n(D) \over 2} \left[ H(D) + H(M-D) - H(M) \right] \mathrm{d}D,

   where :math:`v` is the mass loss rate from progenitors too small to resolve, and :math:`H(M-D)` is omitted if :math:`M-D` is
   below the resolution. This equation is solved numerically, using the rates of the wrapped object. The branch then moves to mass
   :math:`M^\prime`, with complement :math:`M_2=M-M^\prime`, at rate

   .. math::

      R(M^\prime) = {1 \over 2} \left[ n(M^\prime) + n(M_2) \right] {H(M^\prime,\delta) \over H(M,\delta)},

   where the :math:`n(M_2)` term is omitted if :math:`M_2` is below the resolution. Sub-resolution mass loss is unchanged by the
   conditioning, and so is taken directly from the wrapped object.

   Since :math:`H` is an expected number, trees are weighted by the number of progenitors that they have in the window. This is
   equivalent to conditioning on the tree having such a progenitor when that number is rarely more than one, that is, when
   :math:`H(M_0,\delta_0) \ll 1` for the root halo, or when no two progenitors can both lie in the window.

   The tabulation of :math:`H` is shared by all copies of the object, and is rebuilt only if a larger range is needed. A separate
   table is needed for each distinct mass resolution, so a resolution which varies from tree to tree, as with
   :galacticus-class:`mergerTreeMassResolutionScaled`, can make the cost of tabulation significant. When tabulating, the relation
   between critical overdensity and time is evaluated at :math:`M_1`, so a critical overdensity which depends on mass is treated
   only approximately.
   </description>
   <stateStorable>
    <exclude variables="table_"/>
   </stateStorable>
  </mergerTreeBranchingProbability>
  !!]

  type :: constraintConditionedTable
     !!{RST
     A tabulation of the expected number of progenitors in the constraint window, :math:`H(M,\delta)`. Tables are immutable once
     built, so they can be shared freely between copies of an object.
     !!}
     type            (varying_string            )                              :: hash
     double precision                                                          :: massResolution          , logMassStep   , &
          &                                                                       logMassMinimum          , logMassMaximum, &
          &                                                                       deltaMaximum            , deltaMinimum  , &
          &                                                                       deltaStep
     integer                                                                   :: countMass               , countDelta
     double precision                            , allocatable, dimension(:,:) :: H
     type            (constraintConditionedTable), pointer                     :: next           => null()
  end type constraintConditionedTable

  type, extends(mergerTreeBranchingProbabilityClass) :: mergerTreeBranchingProbabilityConstraintConditioned
     !!{RST
     A merger tree branching probability class which conditions a branch on reaching a specified progenitor.
     !!}
     private
     class           (mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
     class           (mergerTreeMassResolutionClass      ), pointer :: mergerTreeMassResolution_       => null()
     class           (cosmologyFunctionsClass            ), pointer :: cosmologyFunctions_             => null()
     class           (criticalOverdensityClass           ), pointer :: criticalOverdensity_            => null()
     class           (cosmologicalMassVarianceClass      ), pointer :: cosmologicalMassVariance_       => null()
     class           (linearGrowthClass                  ), pointer :: linearGrowth_                   => null()
     type            (constraintConditionedTable         ), pointer :: table_                          => null()
     type            (rootFinder                         )          :: finder
     type            (integrator                         )          :: integrator_
     type            (varying_string                     )          :: hash
     ! The constraint, and the numerical controls for the tabulation.
     double precision                                               :: criticalOverdensityConstrained           , &
          &                                                            varianceConstrained                      , &
          &                                                            massConstrained                          , &
          &                                                            redshiftConstrained                      , &
          &                                                            windowLogMass                            , &
          &                                                            stepLogMass                              , &
          &                                                            stepDelta                                , &
          &                                                            massComplementMinimum
     ! The parent halo for which integrals are currently being evaluated.
     double precision                                               :: massParent                               , &
          &                                                            deltaParent                              , &
          &                                                            timeParent                               , &
          &                                                            massResolutionParent                     , &
          &                                                            probabilitySeek
   contains
     !![
     <methods docformat="rst">
       <method method="progenitorExpectation" description="Return the expected number of progenitors in the window about the constraint, for a halo of given mass and critical overdensity."/>
     </methods>
     !!]
     final     ::                          constraintConditionedDestructor
     procedure :: rate                  => constraintConditionedRate
     procedure :: probability           => constraintConditionedProbability
     procedure :: probabilityBound      => constraintConditionedProbabilityBound
     procedure :: fractionSubresolution => constraintConditionedFractionSubresolution
     procedure :: massBranch            => constraintConditionedMassBranch
     procedure :: stepMaximum           => constraintConditionedStepMaximum
     procedure :: progenitorExpectation => constraintConditionedProgenitorExpectation
  end type mergerTreeBranchingProbabilityConstraintConditioned

  interface mergerTreeBranchingProbabilityConstraintConditioned
     !!{RST
     Constructors for the :galacticus-class:`mergerTreeBranchingProbabilityConstraintConditioned` merger tree branching probability class.
     !!}
     module procedure constraintConditionedConstructorParameters
     module procedure constraintConditionedConstructorInternal
  end interface mergerTreeBranchingProbabilityConstraintConditioned

  ! All tables built so far. A table is never removed, because other copies of the object may still point to it. Tables are few,
  ! as one is needed for each combination of object and mass resolution, plus any extensions to larger halo masses.
  type            (constraintConditionedTable                         ), pointer   :: tables           => null()

  ! Module-scope pointers used in integrands and root functions.
  class           (mergerTreeBranchingProbabilityConstraintConditioned), pointer   :: self_
  type            (treeNode                                           ), pointer   :: node_
  !$omp threadprivate(self_,node_)

  ! Factors by which a table is made larger than is currently needed, so that extensions are rare.
  double precision                                                     , parameter :: massMarginFactor =  4.0d0 , deltaMargin=0.2d0

contains

  function constraintConditionedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`mergerTreeBranchingProbabilityConstraintConditioned` merger tree branching probability class which reads parameters from a provided parameter list.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBranchingProbabilityConstraintConditioned)                :: self
    type            (inputParameters                                    ), intent(inout) :: parameters
    class           (mergerTreeBranchingProbabilityClass                ), pointer       :: mergerTreeBranchingProbability_
    class           (mergerTreeMassResolutionClass                      ), pointer       :: mergerTreeMassResolution_
    class           (cosmologyFunctionsClass                            ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass                           ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                      ), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass                                  ), pointer       :: linearGrowth_
    double precision                                                                     :: criticalOverdensityConstrained , &
         &                                                                                  varianceConstrained            , &
         &                                                                                  redshiftConstrained            , &
         &                                                                                  massConstrained                , &
         &                                                                                  timeConstrained                , &
         &                                                                                  timePresent                    , &
         &                                                                                  windowLogMass                  , &
         &                                                                                  stepLogMass                    , &
         &                                                                                  stepDelta

    !![
    <inputParameter docformat="rst">
      <name>windowLogMass</name>
      <defaultValue>0.05d0</defaultValue>
      <description>
      The half-width, in :math:`\ln M`, of the window about the constrained mass into which the branch is conditioned to arrive.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>stepLogMass</name>
      <defaultValue>0.1d0*windowLogMass</defaultValue>
      <description>
      The step in :math:`\ln M` of the tabulation of the expected number of progenitors at the constraint, which also sets the
      spacing of the quadrature used in the tabulation. Must resolve the window, so may be at most one fifth of
      ``[windowLogMass]``.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>stepDelta</name>
      <defaultValue>0.01d0</defaultValue>
      <description>
      The step in critical overdensity of the tabulation of the expected number of progenitors at the constraint.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    <objectBuilder class="mergerTreeMassResolution"       name="mergerTreeMassResolution_"       source="parameters"/>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"/>
    <objectBuilder class="criticalOverdensity"            name="criticalOverdensity_"            source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="cosmologicalMassVariance_"       source="parameters"/>
    <objectBuilder class="linearGrowth"                   name="linearGrowth_"                   source="parameters"/>
    !!]
    timePresent=cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    if      (parameters%isPresent('criticalOverdensityConstrained')) then
       if     (                                                                                                                                                                       &
            &  .not.parameters%isPresent('varianceConstrained')                                                                                                                       &
            & ) call Error_Report('both "criticalOverdensityConstrained" and "varianceConstrained" must be provided'                                      //{introspection:location})
       if     (                                                                                                                                                                       &
            &       parameters%isPresent('redshiftConstrained')                                                                                                                       &
            &  .or.                                                                                                                                                                   &
            &       parameters%isPresent('massConstrained'    )                                                                                                                       &
            & ) call Error_Report('can not mix "criticalOverdensityConstrained/varianceConstrained" and "redshiftConstrained/massConstrained" constraints'//{introspection:location})
       !![
       <inputParameter docformat="rst">
         <name>criticalOverdensityConstrained</name>
         <source>parameters</source>
         <description>
         The critical overdensity :math:`\delta_1` (extrapolated to the present epoch) at which the branch must reach the constraint.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>varianceConstrained</name>
         <source>parameters</source>
         <description>
         The mass variance :math:`\sigma^2(M_1)`, at the present epoch, of the mass which the branch must reach.
         </description>
       </inputParameter>
       !!]
    else if (parameters%isPresent('redshiftConstrained'           )) then
       if     (                                                                                                                                                                       &
            &  .not.parameters%isPresent('massConstrained'    )                                                                                                                       &
            & ) call Error_Report('both "redshiftConstrained" and "massConstrained" must be provided'                                                     //{introspection:location})
       if     (                                                                                                                                                                       &
            &       parameters%isPresent('varianceConstrained')                                                                                                                       &
            & ) call Error_Report('can not mix "criticalOverdensityConstrained/varianceConstrained" and "redshiftConstrained/massConstrained" constraints'//{introspection:location})
       !![
       <inputParameter docformat="rst">
         <name>redshiftConstrained</name>
         <source>parameters</source>
         <description>
         The redshift at which the branch must reach the constraint.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>massConstrained</name>
         <source>parameters</source>
         <description>
         The mass (:math:`\mathrm{M}_\odot`) which the branch must reach.
         </description>
       </inputParameter>
       !!]
       timeConstrained               =+cosmologyFunctions_      %cosmicTime  (cosmologyFunctions_%expansionFactorFromRedshift(redshiftConstrained))
       criticalOverdensityConstrained=+criticalOverdensity_     %value       (time=timeConstrained,mass=massConstrained)                            &
            &                         /linearGrowth_            %value       (time=timeConstrained                     )
       varianceConstrained           =+cosmologicalMassVariance_%rootVariance(time=timePresent    ,mass=massConstrained)**2
    else
       criticalOverdensityConstrained=0.0d0
       varianceConstrained           =0.0d0
       call Error_Report('must provide either [criticalOverdensityConstrained] and [varianceConstrained], or [redshiftConstrained] and [massConstrained]'//{introspection:location})
    end if
    self=mergerTreeBranchingProbabilityConstraintConditioned(criticalOverdensityConstrained,varianceConstrained,windowLogMass,stepLogMass,stepDelta,mergerTreeBranchingProbability_,mergerTreeMassResolution_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    <objectDestructor name="mergerTreeMassResolution_"      />
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="criticalOverdensity_"           />
    <objectDestructor name="cosmologicalMassVariance_"      />
    <objectDestructor name="linearGrowth_"                  />
    !!]
    return
  end function constraintConditionedConstructorParameters

  function constraintConditionedConstructorInternal(criticalOverdensityConstrained,varianceConstrained,windowLogMass,stepLogMass,stepDelta,mergerTreeBranchingProbability_,mergerTreeMassResolution_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`mergerTreeBranchingProbabilityConstraintConditioned` merger tree branching probability class.
    !!}
    use :: Error                , only : Error_Report
    use :: ISO_Varying_String   , only : assignment(=)
    use :: Numerical_Integration, only : GSL_Integ_Gauss15
    implicit none
    type            (mergerTreeBranchingProbabilityConstraintConditioned)                        :: self
    class           (mergerTreeBranchingProbabilityClass                ), intent(in   ), target :: mergerTreeBranchingProbability_
    class           (mergerTreeMassResolutionClass                      ), intent(in   ), target :: mergerTreeMassResolution_
    class           (cosmologyFunctionsClass                            ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass                           ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                      ), intent(in   ), target :: cosmologicalMassVariance_
    class           (linearGrowthClass                                  ), intent(in   ), target :: linearGrowth_
    double precision                                                     , intent(in   )         :: criticalOverdensityConstrained        , &
         &                                                                                          varianceConstrained                   , &
         &                                                                                          windowLogMass                         , &
         &                                                                                          stepLogMass                           , &
         &                                                                                          stepDelta
    double precision                                                     , parameter             :: toleranceRelativeRoot          =1.0d-09, &
         &                                                                                          toleranceIntegral              =1.0d-02, &
         &                                                                                          toleranceIntegralAbsolute      =1.0d-10, &
         &                                                                                          stepLogMassWindowMaximum       =0.2d+00
    !![
    <constructorAssign variables="criticalOverdensityConstrained, varianceConstrained, windowLogMass, stepLogMass, stepDelta, *mergerTreeBranchingProbability_, *mergerTreeMassResolution_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *linearGrowth_"/>
    !!]

    ! A coarser tabulation does not resolve the window, and the expectation is then not converged.
    if (stepLogMass > stepLogMassWindowMaximum*windowLogMass) call Error_Report('[stepLogMass] must be at most one fifth of [windowLogMass]'//{introspection:location})
    ! Only the full-range generalized Press-Schechter process has the two-progenitor structure assumed by the conditioning.
    select type (wrapped => self%mergerTreeBranchingProbability_)
    class is (mergerTreeBranchingProbabilityGeneralizedPressSchechter)
       if (     wrapped%distributionFunctionLowerHalfOnly) call Error_Report('the wrapped generalized Press-Schechter object must have [distributionFunctionLowerHalfOnly]=false'//{introspection:location})
       if (.not.wrapped%distributionFunctionNormalize    ) call Error_Report('the wrapped generalized Press-Schechter object must have [distributionFunctionNormalize]=true'     //{introspection:location})
       ! σ(M) is never needed for a complement far below the minimum mass used by the wrapped object for sub-resolution accretion.
       self%massComplementMinimum=1.0d-3*wrapped%massMinimum
    class default
       call Error_Report('the wrapped merger tree branching probability object must be of the "generalizedPressSchechter" class'//{introspection:location})
    end select
    ! Find the mass and redshift of the constraint, as is done by the constrained merger tree build controller.
    self%massConstrained    =self%cosmologicalMassVariance_%mass                       (time               =self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0),rootVariance=sqrt(self%varianceConstrained))
    self%redshiftConstrained=self%cosmologyFunctions_      %redshiftFromExpansionFactor(                                                                                                                           &
         &                   self%cosmologyFunctions_      %expansionFactor            (                                                                                                                           &
         &                   self%criticalOverdensity_     %timeOfCollapse             (criticalOverdensity=self%criticalOverdensityConstrained,mass=self%massConstrained)                                         &
         &                                                                             )                                                                                                                           &
         &                                                                             )
    self%hash               =''
    self%finder             =rootFinder(                                                       &
         &                              rootFunction     =constraintConditionedMassBranchRoot, &
         &                              toleranceAbsolute=0.0d0                              , &
         &                              toleranceRelative=toleranceRelativeRoot                &
         &                             )
    ! Far from the constraint, the expectation in the window can underflow, so that a purely relative tolerance is unattainable
    ! there. The integral is a probability per unit change in δ_crit, so an absolute tolerance of this size is negligible.
    self%integrator_        =integrator(                                                       &
         &                                                constraintConditionedIntegrand     , &
         &                              toleranceAbsolute=toleranceIntegralAbsolute          , &
         &                              toleranceRelative=toleranceIntegral                  , &
         &                              integrationRule  =GSL_Integ_Gauss15                    &
         &                             )
    return
  end function constraintConditionedConstructorInternal

  subroutine constraintConditionedDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`mergerTreeBranchingProbabilityConstraintConditioned` merger tree branching probability class.
    !!}
    implicit none
    type(mergerTreeBranchingProbabilityConstraintConditioned), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    <objectDestructor name="self%mergerTreeMassResolution_"      />
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%criticalOverdensity_"           />
    <objectDestructor name="self%cosmologicalMassVariance_"      />
    <objectDestructor name="self%linearGrowth_"                  />
    !!]
    return
  end subroutine constraintConditionedDestructor

  double precision function constraintConditionedRate(self,mass,deltaCritical,time,massBranch,node)
    !!{RST
    Return the rate per unit mass and per unit change in :math:`\delta_\mathrm{crit}` at which the branch moves to ``massBranch``.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout), target :: self
    double precision                                                     , intent(in   )         :: mass          , deltaCritical, &
         &                                                                                          time          , massBranch
    type            (treeNode                                           ), intent(inout), target :: node
    double precision                                                                             :: massResolution

    massResolution=self%mergerTreeMassResolution_%resolution(node%hostTree)
    call constraintConditionedTableEnsure(self,mass,deltaCritical,massResolution,node)
    constraintConditionedRate=constraintConditionedRateEvaluate(self,mass,deltaCritical,time,massBranch,massResolution,node)
    return
  end function constraintConditionedRate

  double precision function constraintConditionedRateEvaluate(self,mass,deltaCritical,time,massBranch,massResolution,node) result(rate)
    !!{RST
    Evaluate the conditioned rate at which a branch of the given mass moves to ``massBranch``.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout)         :: self
    double precision                                                     , intent(in   )         :: mass          , deltaCritical, &
         &                                                                                          time          , massBranch   , &
         &                                                                                          massResolution
    type            (treeNode                                           ), intent(inout), target :: node
    double precision                                                                             :: massComplement, ratio        , &
         &                                                                                          massDrawn

    rate          =0.0d0
    massComplement=mass-massBranch
    if (massComplement <= 0.0d0) return
    ratio=constraintConditionedLookup(self,log(mass      ),deltaCritical)
    if (ratio <= 0.0d0) return
    ratio=constraintConditionedLookup(self,log(massBranch),deltaCritical)/ratio
    if (ratio <= 0.0d0) return
    ! The wrapped object evaluates σ(M) for the complement, so keep that mass away from zero.
    massDrawn=min(massBranch,mass-self%massComplementMinimum)
    rate     =self%mergerTreeBranchingProbability_%rate(mass,deltaCritical,time,massDrawn,node)
    if (massComplement >= massResolution) rate=+rate+self%mergerTreeBranchingProbability_%rate(mass,deltaCritical,time,massComplement,node)
    rate     =+0.5d0*ratio*rate
    return
  end function constraintConditionedRateEvaluate

  double precision function constraintConditionedIntegrand(massBranch) result(integrand)
    !!{RST
    Integrand for the probability of branching.
    !!}
    implicit none
    double precision, intent(in   ) :: massBranch

    integrand=constraintConditionedRateEvaluate(self_,self_%massParent,self_%deltaParent,self_%timeParent,massBranch,self_%massResolutionParent,node_)
    return
  end function constraintConditionedIntegrand

  double precision function constraintConditionedIntegrate(massLower,massUpper) result(integral)
    !!{RST
    Integrate the conditioned branching rate between the given masses, splitting the range at the sharp features of the integrand.
    !!}
    implicit none
    double precision, intent(in   ) :: massLower, massUpper
    double precision, dimension(4)  :: breaks
    double precision                :: massStart, value_
    integer                         :: i        , j

    breaks=[                                                                    &
         &        self_%massConstrained*exp(-3.0d0*self_%windowLogMass       ), &
         &        self_%massConstrained*exp(+3.0d0*self_%windowLogMass       ), &
         &  0.5d0*self_%massParent                                            , &
         &        self_%massParent     -           self_%massResolutionParent   &
         & ]
    do i=2,size(breaks)
       value_=breaks(i)
       j     =i-1
       do while (j >= 1)
          if (breaks(j) <= value_) exit
          breaks(j+1)=breaks(j)
          j          =j-1
       end do
       breaks(j+1)=value_
    end do
    integral =0.0d0
    massStart=massLower
    do i=1,size(breaks)
       if (breaks(i) <= massStart .or. breaks(i) >= massUpper) cycle
       integral =+integral+self_%integrator_%integrate(massStart,breaks(i))
       massStart=breaks(i)
    end do
    if (massUpper > massStart) integral=+integral+self_%integrator_%integrate(massStart,massUpper)
    return
  end function constraintConditionedIntegrate

  double precision function constraintConditionedProbability(self,haloMass,deltaCritical,time,massResolution,node)
    !!{RST
    Return the probability per unit change in :math:`\delta_\mathrm{crit}` that the branch undergoes a resolved branching.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout), target :: self
    double precision                                                     , intent(in   )         :: deltaCritical , haloMass, &
         &                                                                                          massResolution, time
    type            (treeNode                                           ), intent(inout), target :: node

    constraintConditionedProbability=0.0d0
    if (haloMass <= massResolution) return
    call constraintConditionedTableEnsure(self,haloMass,deltaCritical     ,massResolution,node)
    call constraintConditionedParentSet  (self,haloMass,deltaCritical,time,massResolution,node)
    constraintConditionedProbability=constraintConditionedIntegrate(massResolution,haloMass)
    return
  end function constraintConditionedProbability

  double precision function constraintConditionedProbabilityBound(self,haloMass,deltaCritical,time,massResolution,bound,node)
    !!{RST
    Return a bound on the probability per unit change in :math:`\delta_\mathrm{crit}` that the branch undergoes a resolved branching. The probability is computed exactly, so this is used for both bounds.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout)         :: self
    double precision                                                     , intent(in   )         :: deltaCritical , haloMass, &
         &                                                                                          massResolution, time
    type            (enumerationMergerTreeBranchingBoundType            ), intent(in   )         :: bound
    type            (treeNode                                           ), intent(inout), target :: node
    !$GLC attributes unused :: bound

    constraintConditionedProbabilityBound=self%probability(haloMass,deltaCritical,time,massResolution,node)
    return
  end function constraintConditionedProbabilityBound

  double precision function constraintConditionedFractionSubresolution(self,haloMass,deltaCritical,time,massResolution,node)
    !!{RST
    Return the fraction of mass accreted in sub-resolution halos per unit change in :math:`\delta_\mathrm{crit}`. This is unaffected by the conditioning, as such halos are side branches, separate from the mass element which reaches the constraint.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout), target :: self
    double precision                                                     , intent(in   )         :: deltaCritical , haloMass, &
         &                                                                                          massResolution, time
    type            (treeNode                                           ), intent(inout), target :: node

    constraintConditionedFractionSubresolution=self%mergerTreeBranchingProbability_%fractionSubresolution(haloMass,deltaCritical,time,massResolution,node)
    return
  end function constraintConditionedFractionSubresolution

  double precision function constraintConditionedStepMaximum(self,haloMass,deltaCritical,time,massResolution)
    !!{RST
    Return the maximum allowed step in :math:`\delta_\mathrm{crit}`, as given by the wrapped object.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout) :: self
    double precision                                                     , intent(in   ) :: deltaCritical , haloMass, &
         &                                                                                  massResolution, time

    constraintConditionedStepMaximum=self%mergerTreeBranchingProbability_%stepMaximum(haloMass,deltaCritical,time,massResolution)
    return
  end function constraintConditionedStepMaximum

  double precision function constraintConditionedMassBranch(self,haloMass,deltaCritical,time,massResolution,probabilityFraction,randomNumberGenerator_,node)
    !!{RST
    Return the mass to which the branch moves, given the cumulative branching probability ``probabilityFraction``.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout), target :: self
    double precision                                                     , intent(in   )         :: deltaCritical                  , &
         &                                                                                          haloMass                       , &
         &                                                                                          massResolution                 , &
         &                                                                                          probabilityFraction            , &
         &                                                                                          time
    class           (randomNumberGeneratorClass                         ), intent(inout)         :: randomNumberGenerator_
    type            (treeNode                                           ), intent(inout), target :: node
    double precision                                                     , parameter             :: smallProbabilityFraction=1.0d-3
    double precision                                                                             :: rootFunctionLower              , &
         &                                                                                          rootFunctionUpper
    !$GLC attributes unused :: randomNumberGenerator_

    call constraintConditionedTableEnsure(self,haloMass,deltaCritical     ,massResolution,node)
    call constraintConditionedParentSet  (self,haloMass,deltaCritical,time,massResolution,node)
    self%probabilitySeek=probabilityFraction
    rootFunctionLower   =constraintConditionedMassBranchRoot(massResolution)
    rootFunctionUpper   =constraintConditionedMassBranchRoot(haloMass      )
    if (rootFunctionLower*rootFunctionUpper >= 0.0d0) then
       ! The probability sought lies at the top of the range. This is acceptable only if it does so to within numerical accuracy.
       if (rootFunctionUpper > probabilityFraction*smallProbabilityFraction) call Error_Report('numerical accuracy lost in root finding calculation'//{introspection:location})
       constraintConditionedMassBranch=haloMass
       return
    end if
    constraintConditionedMassBranch=self%finder%find(rootRange=[massResolution,haloMass],rootRangeValues=[rootFunctionLower,rootFunctionUpper])
    return
  end function constraintConditionedMassBranch

  double precision function constraintConditionedMassBranchRoot(massMaximum)
    !!{RST
    Root function used in solving for the branch mass.
    !!}
    implicit none
    double precision, intent(in   ) :: massMaximum

    constraintConditionedMassBranchRoot=+self_%probabilitySeek                                                  &
         &                              -constraintConditionedIntegrate(self_%massResolutionParent,massMaximum)
    return
  end function constraintConditionedMassBranchRoot

  double precision function constraintConditionedProgenitorExpectation(self,mass,deltaCritical,massResolution,node)
    !!{RST
    Return the expected number of progenitors in the window about the constraint, for a halo of the given mass and critical overdensity.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout)         :: self
    double precision                                                     , intent(in   )         :: mass          , deltaCritical, &
         &                                                                                          massResolution
    type            (treeNode                                           ), intent(inout), target :: node

    call constraintConditionedTableEnsure(self,mass,deltaCritical,massResolution,node)
    constraintConditionedProgenitorExpectation=constraintConditionedLookup(self,log(mass),deltaCritical)
    return
  end function constraintConditionedProgenitorExpectation

  subroutine constraintConditionedParentSet(self,mass,deltaCritical,time,massResolution,node)
    !!{RST
    Record the parent halo for which integrals are to be evaluated.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout), target :: self
    double precision                                                     , intent(in   )         :: mass, deltaCritical , &
         &                                                                                          time, massResolution
    type            (treeNode                                           ), intent(inout), target :: node

    self_                     => self
    node_                     => node
    self%massParent           =  mass
    self%deltaParent          =  deltaCritical
    self%timeParent           =  time
    self%massResolutionParent =  massResolution
    return
  end subroutine constraintConditionedParentSet

  double precision function constraintConditionedLookup(self,logMass,delta) result(H)
    !!{RST
    Interpolate in the table of the expected number of progenitors at the constraint. This is zero below the mass resolution.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(in   ) :: self
    double precision                                                     , intent(in   ) :: logMass, delta
    double precision                                                                     :: u      , v
    integer                                                                              :: i      , j

    H      =0.0d0
    u      =(logMass-self%table_%logMassMinimum)/self%table_%logMassStep
    if (u < 0.0d0) return
    i      =min(int(u),self%table_%countMass-2)
    u      =min(u-dble(i),1.0d0)
    v      =max(min((self%table_%deltaMaximum-delta)/self%table_%deltaStep,dble(self%table_%countDelta-1)),0.0d0)
    j      =min(int(v),self%table_%countDelta-2)
    v      =v-dble(j)
    H      =+(1.0d0-v)*((1.0d0-u)*self%table_%H(i+1,j+1)+u*self%table_%H(i+2,j+1)) &
         &  +       v *((1.0d0-u)*self%table_%H(i+1,j+2)+u*self%table_%H(i+2,j+2))
    return
  end function constraintConditionedLookup

  subroutine constraintConditionedTableEnsure(self,mass,deltaCritical,massResolution,node)
    !!{RST
    Ensure that this object points to a table which covers the given mass and critical overdensity at the given resolution. Thread-safe, as the first thread to need a table builds it for all others.
    !!}
    use :: ISO_Varying_String, only : assignment(=), operator(==)
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout)         :: self
    double precision                                                     , intent(in   )         :: mass          , deltaCritical, &
         &                                                                                          massResolution
    type            (treeNode                                           ), intent(inout), target :: node
    type            (constraintConditionedTable                         ), pointer               :: table
    double precision                                                                             :: massCovered   , deltaCovered

    if (associated(self%table_)) then
       if     (                                                          &
            &   self%table_%massResolution ==             massResolution &
            &  .and.                                                     &
            &   log(mass)                  <= self%table_%logMassMaximum &
            &  .and.                                                     &
            &   deltaCritical              >= self%table_%deltaMinimum   &
            & ) return
    end if
    if (self%hash == '') self%hash=self%hashedDescriptor()
    !$omp critical(mergerTreeBranchingConstraintConditioned)
    massCovered =mass
    ! Tables always extend from the constraint, and lookups above it are clamped to the terminal column.
    deltaCovered=min(deltaCritical,self%criticalOverdensityConstrained)
    nullify(self%table_)
    table => tables
    do while (associated(table))
       if (table%hash == self%hash .and. table%massResolution == massResolution) then
          if (log(mass) <= table%logMassMaximum .and. deltaCritical >= table%deltaMinimum) then
             self%table_ => table
             exit
          end if
          ! Any replacement must also cover everything that this table covers.
          massCovered =max(massCovered ,exp(table%logMassMaximum))
          deltaCovered=min(deltaCovered,    table%deltaMinimum   )
       end if
       table => table%next
    end do
    if (.not.associated(self%table_)) then
       allocate(table)
       call constraintConditionedTableBuild(self,table,massCovered*massMarginFactor,deltaCovered-deltaMargin,massResolution,node)
       table  %hash  =  self%hash
       table  %next  =>      tables
       tables        =>      table
       self  %table_ =>      table
    end if
    !$omp end critical(mergerTreeBranchingConstraintConditioned)
    return
  end subroutine constraintConditionedTableEnsure

  subroutine constraintConditionedTableBuild(self,table,massMaximum,deltaMinimum,massResolution,node)
    !!{RST
    Solve the backward equation for the expected number of progenitors in the window about the constraint, stepping from the
    constraint toward lower critical overdensity. Both grids are anchored (at the mass resolution, and at the constraint), so a
    larger table has the same values at the points that it shares with a smaller one, since the solution at each mass depends only
    on lower masses.
    !!}
    use :: Display                        , only : displayMessage              , verbosityLevelWorking
    use :: Numerical_Integration_Tabulated, only : integrateTrapezoidal
    use :: Tables                         , only : table1DLinearMonotoneCSpline
    implicit none
    class           (mergerTreeBranchingProbabilityConstraintConditioned), intent(inout)               :: self
    type            (constraintConditionedTable                         ), intent(inout)               :: table
    double precision                                                     , intent(in   )               :: massMaximum                 , &
         &                                                                                                deltaMinimum                , &
         &                                                                                                massResolution
    type            (treeNode                                           ), intent(inout), target       :: node
    integer                                                              , parameter                   :: countDrift           =200   , &
         &                                                                                                stepsRebuild         =25    , &
         &                                                                                                countSample          =20    , &
         &                                                                                                kindLower            =1     , &
         &                                                                                                kindUpper            =2     , &
         &                                                                                                kindSecondary        =3
    ! Draws of complements smaller than a fraction of the parent are treated as a continuous drift. That fraction must be small
    ! compared to the window, as otherwise the drift smears out the edges of the window.
    double precision                                                     , parameter                   :: fractionDriftMaximum =0.03d0, &
         &                                                                                                fractionDriftWindow  =0.10d0, &
         &                                                                                                stepLogMassCoarse    =0.02d0, &
         &                                                                                                stepLogSubresolution =0.20d0
    double precision                                                                                   :: fractionDrift               , &
         &                                                                                                logMassEdge
    double precision                                                     , allocatable  , dimension(:) :: logMass                     , &
         &                                                                                                H                           , &
         &                                                                                                gain                        , &
         &                                                                                                loss                        , &
         &                                                                                                drift                       , &
         &                                                                                                HShifted                    , &
         &                                                                                                logDrawn                    , &
         &                                                                                                logComplement               , &
         &                                                                                                weight                      , &
         &                                                                                                fractionSubresolution       , &
         &                                                                                                growthFactor
    integer                                                              , allocatable  , dimension(:) :: offset
    type            (table1DLinearMonotoneCSpline                       )                              :: spline
    integer                                                                                            :: i                           , &
         &                                                                                                k                           , &
         &                                                                                                n                           , &
         &                                                                                                count_                      , &
         &                                                                                                countNodes                  , &
         &                                                                                                indexActive                 , &
         &                                                                                                strideSubresolution
    double precision                                                                                   :: delta                       , &
         &                                                                                                timePresent                 , &
         &                                                                                                growthCurrent

    call displayMessage('tabulating progenitor expectation for constrained branch',verbosityLevelWorking)
    table%massResolution=massResolution
    table%logMassStep   =self%stepLogMass
    table%logMassMinimum=log(massResolution)
    table%countMass     =ceiling((log(massMaximum)-table%logMassMinimum)/table%logMassStep)+3
    table%logMassMaximum=table%logMassMinimum+dble(table%countMass-3)*table%logMassStep
    table%deltaStep     =self%stepDelta
    table%deltaMaximum  =self%criticalOverdensityConstrained
    table%countDelta    =ceiling((table%deltaMaximum-deltaMinimum)/table%deltaStep)+2
    table%deltaMinimum  =table%deltaMaximum-dble(table%countDelta-2)*table%deltaStep
    allocate(table%H(table%countMass,table%countDelta))
    allocate(logMass              (table%countMass  ))
    allocate(H                    (table%countMass  ))
    allocate(gain                 (table%countMass  ))
    allocate(loss                 (table%countMass  ))
    allocate(drift                (table%countMass  ))
    allocate(HShifted             (table%countMass  ))
    allocate(fractionSubresolution(table%countMass  ))
    allocate(growthFactor         (table%countMass  ))
    allocate(offset               (table%countMass+1))
    fractionDrift      =min(fractionDriftMaximum,fractionDriftWindow*self%windowLogMass)
    strideSubresolution=max(1,nint(stepLogSubresolution/table%logMassStep))
    ! Halos can not gain mass going back in time, so the expectation vanishes, at all critical overdensities, for masses below the
    ! window. Two cells of margin allow for the terminal condition and the interpolation. Masses below this edge need no quadrature.
    logMassEdge=log(self%massConstrained)-self%windowLogMass-2.0d0*table%logMassStep
    do i=1,table%countMass
       logMass(i)=table%logMassMinimum+dble(i-1)*table%logMassStep
    end do
    call spline%create(table%logMassMinimum,logMass(table%countMass),table%countMass)
    indexActive=table%countMass+1
    do i=1,table%countMass
       if (logMass(i) >= logMassEdge) then
          indexActive=i
          exit
       end if
    end do
    countNodes=0
    do i=indexActive,table%countMass
       call quadratureMass(exp(logMass(i)),countNodes)
    end do
    allocate(logDrawn     (countNodes))
    allocate(logComplement(countNodes))
    allocate(weight       (countNodes))
    ! The terminal condition is the fraction of each cell which lies inside the window.
    do i=1,table%countMass
       count_=0
       do k=-countSample,countSample
          if (abs(logMass(i)+0.5d0*dble(k)/dble(countSample)*table%logMassStep-log(self%massConstrained)) < self%windowLogMass) count_=count_+1
       end do
       H(i)=dble(count_)/dble(2*countSample+1)
    end do
    table%H(:,1)=H
    timePresent =self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    do n=2,table%countDelta
       delta=table%deltaMaximum-dble(n-1)*table%deltaStep
       if (mod(n-2,stepsRebuild) == 0) call rebuild(delta)
       call spline%populate(H)
       do i=1,table%countMass
          gain(i)=0.0d0
          do k=offset(i),offset(i+1)-1
             gain(i)=+gain(i)+weight(k)*(interpolate(logDrawn(k))+interpolate(logComplement(k)))
          end do
          HShifted(i)=interpolate(logMass(i)-drift(i)*table%deltaStep)
       end do
       H           =max(HShifted+table%deltaStep*(gain-loss*H),0.0d0)
       table%H(:,n)=H
    end do
    call spline%destroy()
    return

  contains

    subroutine rebuild(deltaNow)
      !!{RST
      Compute the quadrature weights, loss rates, and drift rates for all masses at the given critical overdensity. Draws of small complements are treated as a continuous drift.
      !!}
      double precision, intent(in   )            :: deltaNow
      double precision, dimension(countDrift)    :: logComplementDrift, integrandDrift
      double precision                           :: time              , massParent    , &
           &                                        massLower         , massUpper     , &
           &                                        step              , x             , &
           &                                        integral
      integer                                    :: iMass             , j             , &
           &                                        k_                , iLower

      time=self%criticalOverdensity_%timeOfCollapse(criticalOverdensity=deltaNow,mass=self%massConstrained,node=node)
      ! Branching rates are per unit critical overdensity at the current time, while the critical overdensity used here is scaled
      ! to the present day. The tree builder converts between these using the growth factor of the root variance, and so must we.
      do iMass=1,table%countMass
         growthFactor(iMass)=+self%cosmologicalMassVariance_%rootVariance(time=time       ,mass=exp(logMass(iMass))) &
              &              /self%cosmologicalMassVariance_%rootVariance(time=timePresent,mass=exp(logMass(iMass)))
      end do
      ! The sub-resolution fraction is expensive, but smooth in mass, so is evaluated only on a coarse subset of the masses.
      fractionSubresolution=0.0d0
      do iMass=indexActive,table%countMass
         if (mod(iMass-indexActive,strideSubresolution) == 0 .or. iMass == table%countMass)                                                                                              &
              & fractionSubresolution(iMass)=+self%mergerTreeBranchingProbability_%fractionSubresolution(exp(logMass(iMass)),deltaNow,time,min(massResolution,exp(logMass(iMass))),node) &
              &                              *growthFactor(iMass)
      end do
      do iMass=indexActive,table%countMass
         if (mod(iMass-indexActive,strideSubresolution) == 0 .or. iMass == table%countMass) cycle
         iLower=indexActive+strideSubresolution*((iMass-indexActive)/strideSubresolution)
         j     =min(iLower+strideSubresolution,table%countMass)
         fractionSubresolution(iMass)=+fractionSubresolution(iLower)                            &
              &                       +(fractionSubresolution(j)-fractionSubresolution(iLower)) &
              &                       *dble(iMass-iLower)                                        &
              &                       /dble(j    -iLower)
      end do
      k_=1
      do iMass=1,table%countMass
         massParent   =exp(logMass(iMass))
         growthCurrent=growthFactor(iMass)
         offset(iMass)=k_
         loss  (iMass)=0.0d0
         drift (iMass)=0.0d0
         if (iMass < indexActive) cycle
         call quadratureMass(massParent,k_,deltaNow,time,loss(iMass))
         ! Drift from small complements, integrated over the logarithm of their mass.
         integral =0.0d0
         massLower=log(max(1.0d-10*massParent,self%massComplementMinimum))
         massUpper=log(fractionDrift*massParent)
         if (massUpper > massLower) then
            step=(massUpper-massLower)/dble(countDrift-1)
            do j=1,countDrift
               logComplementDrift(j)=massLower+dble(j-1)*step
               x                    =exp(logComplementDrift(j))
               integrandDrift    (j)=x*drawRate(massParent,deltaNow,time,massParent-x)*x
            end do
            integral=integrateTrapezoidal(logComplementDrift,integrandDrift)
         end if
         drift(iMass)=integral/massParent+fractionSubresolution(iMass)
      end do
      offset(table%countMass+1)=k_
      return
    end subroutine rebuild

    subroutine quadratureMass(massParent,countNodes,deltaNow,time,loss)
      !!{RST
      Count, or, if ``deltaNow`` is present, also set, the quadrature nodes for a parent halo: draws in the lower half, draws in the
      upper half above the drift zone (in terms of the complement), and resolved secondaries inside the drift zone, which can host
      progenitors while the drawn mass itself is handled by the drift.
      !!}
      double precision, intent(in   )           :: massParent
      integer         , intent(inout)           :: countNodes
      double precision, intent(in   ), optional :: deltaNow  , time
      double precision, intent(inout), optional :: loss

      if (              massParent > 2.0d0*massResolution) call quadratureSegment(massParent,log(massResolution          ),log(0.5d0        *massParent),kindLower    ,countNodes,deltaNow,time,loss)
      call                                                      quadratureSegment(massParent,log(fractionDrift*massParent),log(0.5d0        *massParent),kindUpper    ,countNodes,deltaNow,time,loss)
      if (fractionDrift*massParent >       massResolution) call quadratureSegment(massParent,log(massResolution          ),log(fractionDrift*massParent),kindSecondary,countNodes,deltaNow,time,loss)
      return
    end subroutine quadratureMass

    subroutine quadratureSegment(massParent,logMassLower,logMassUpper,kind,countNodes,deltaNow,time,loss)
      !!{RST
      Count, or set, the trapezoidal quadrature nodes for one segment, in :math:`\ln x` where :math:`x` is the drawn mass (lower
      half) or the complement (otherwise). Above the window edge, where the expectation at :math:`x` can be non-zero, nodes are
      spaced no more widely than the mass grid, so that the window is resolved and the quadrature is refined along with the grid.
      Below it only the expectation at :math:`M-x` matters, which varies more slowly with :math:`\ln x` by a factor
      :math:`x/(M-x)`, so coarser nodes suffice, and secondaries there can host no progenitors.
      !!}
      double precision, intent(in   )           :: massParent  , logMassLower, &
           &                                       logMassUpper
      integer         , intent(in   )           :: kind
      integer         , intent(inout)           :: countNodes
      double precision, intent(in   ), optional :: deltaNow    , time
      double precision, intent(inout), optional :: loss
      double precision                          :: logMassSplit, a           , &
           &                                       b           , spacing     , &
           &                                       step        , x
      integer                                   :: part        , count       , &
           &                                       j

      logMassSplit=max(min(logMassEdge,logMassUpper),logMassLower)
      do part=1,2
         if (part == 1) then
            if (kind == kindSecondary) cycle
            a      =logMassLower
            b      =logMassSplit
            spacing=min(stepLogMassCoarse,table%logMassStep*(massParent-exp(b))/exp(b))
         else
            a      =logMassSplit
            b      =logMassUpper
            spacing=table%logMassStep
         end if
         if (b <= a) cycle
         count=max(ceiling((b-a)/spacing),1)+1
         if (.not.present(deltaNow)) then
            countNodes=countNodes+count
            cycle
         end if
         step=(b-a)/dble(count-1)
         do j=1,count
            x=exp(a+dble(j-1)*step)
            select case (kind)
            case (kindLower    )
               logDrawn     (countNodes)=log(x)
               logComplement(countNodes)=log(massParent-x)
               weight       (countNodes)=drawRate(massParent,deltaNow,time,x)*x*step
            case (kindUpper    )
               logDrawn     (countNodes)=log(massParent-x)
               logComplement(countNodes)=log(x)
               weight       (countNodes)=drawRate(massParent,deltaNow,time,massParent-x)*x*step
            case (kindSecondary)
               logDrawn     (countNodes)=-huge(0.0d0)
               logComplement(countNodes)=log(x)
               weight       (countNodes)=drawRate(massParent,deltaNow,time,massParent-x)*x*step
            end select
            if (j == 1 .or. j == count) weight(countNodes)=0.5d0*weight(countNodes)
            if (kind /= kindSecondary) loss=loss+weight(countNodes)
            countNodes=countNodes+1
         end do
      end do
      return
    end subroutine quadratureSegment

    double precision function drawRate(massParent,deltaNow,time,massDrawn)
      !!{RST
      The rate per unit mass and per unit present-day-scaled critical overdensity at which the given mass is drawn.
      !!}
      double precision, intent(in   ) :: massParent, deltaNow , &
           &                             time      , massDrawn

      drawRate=0.5d0*self%mergerTreeBranchingProbability_%rate(massParent,deltaNow,time,massDrawn,node)*growthCurrent
      return
    end function drawRate

    double precision function interpolate(x)
      !!{RST
      Monotone cubic interpolation of the current solution. Zero below the mass resolution.
      !!}
      double precision, intent(in   ) :: x

      interpolate=0.0d0
      if (x < table%logMassMinimum) return
      interpolate=spline%interpolate(x)
      return
    end function interpolate

  end subroutine constraintConditionedTableBuild
