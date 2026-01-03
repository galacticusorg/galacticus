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

!+    Contributions to this file made by: Andrew Benson, Ethan Nadler.

!!{
Implements an excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint
method to perform the integrations \citep{du_substructure_2017}, and with a Brownian bridge constraint.
!!}
  
  use :: Cosmological_Density_Field, only : criticalOverdensityClass
  use :: Linear_Growth             , only : linearGrowthClass 

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingFarahiMidpointBrownianBridge">
    <description>
      An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method
      to perform the integrations \citep{du_substructure_2017}, and with a \href{https://en.wikipedia.org/wiki/Brownian_bridge}{Brownian bridge} constraint.

      Specifically, the trajectories are constrained to pass through a point $(S_2,\delta_2)$ (specified by the parameter
      {\normalfont \ttfamily [varianceConstrained]} and {\normalfont \ttfamily [criticalOverdensityConstrained]}, or equivalently
      by the parameters {\normalfont \ttfamily [massConstrained]} and {\normalfont \ttfamily [timeConstrained]}---note that
      $(S_2,\delta_2)$ here follow the convention in excursion set literature that $S_2$ is the variance evaluated at the present
      day, while $\delta_2$ the the critical overdensity for collapse divided by the linear growth factor), and, of course, always
      pass through the initial point $(S_1,\delta_1)$ corresponding to the current halo.

      For a Brownian bridge the distribution of $\delta$ at some $S$ (where $S_1 \le S \le S_2$), $P_0(\delta|S)$, is given by a normal distribution with mean
      \begin{equation}
      \mu(S) = \delta_1 + \frac{\delta_2-\delta_1}{S_2-S_1}(S - S_1)
      \end{equation}
      and variance
      \begin{equation}
      \mathrm{Var}(S) =  \frac{(S_2-S)(S-S_1)}{S_2-S_1},
      \end{equation}
      and the covariance between two points $S_\mathrm{a}$ and $S_\mathrm{b} > S_\mathrm{a}$ is given by,
      \begin{equation}
      \mathrm{Cov}(S_\mathrm{a},S_\mathrm{b}) = \frac{(S_2-S_\mathrm{b})(S_\mathrm{a}-S_1)}{S_2-S_1}.
      \end{equation}      
      Therefore, the same approach to solving for the first crossing distribution as was utilized by \cite{benson_dark_2012} and
      improved by \citep{du_substructure_2017} can be used (see \refClass{excursionSetFirstCrossingFarahi} and
      \refClass{excursionSetFirstCrossingFarahiMidpoint} for details), with just the appropriate change in the effective offset,
      $\Delta \delta$, and residual variance, $\Delta S$.

      Considering two points $(S,\delta)$ and $(\tilde{S},\tilde{\delta})$ the effective offset is just the difference in their offsets relative to their local means:
      \begin{equation}
      \Delta \delta = [ \delta - \mu(S) ] - [ \tilde{\delta} - \mu(\tilde{S}) ] = \delta - \tilde{\delta} - \frac{S-\tilde{S}}{S_2-S_1}(\delta_2-\delta_1),
      \end{equation}
      while the residual variance is, as always, just the variance at $(S,\delta)$ minus the covariance between the two points:
      \begin{equation}
      \Delta S = \mathrm{Var}(S) - \mathrm{Cov}(B({\tilde{S}}),\delta) = \frac{(S_2-S)(S-S_1)}{S_2-S_1} - \frac{(S_2-S)(\tilde{S}-S_1)}{S_2-S_1}.
      \end{equation}
      which simplifies to
      \begin{equation}
      \Delta S = \frac{(S_2-S)(S-\tilde{S})}{S_2-S_1}.
      \end{equation}

      Note that, in solving for the first crossing distribution we must also evaluate terms of the form
      \begin{equation}
      \int_{-\infty}^{\delta} P_{0}(\delta^\prime,S) \mathrm{d}\delta^\prime = \mathrm{erf}\left( \frac{\delta^\prime - \mu(S)}{\sqrt{2 \mathrm{Var}(S)}}\right).
      \end{equation}
      In these cases we still use the residual variance since $\Delta S \rightarrow \mathrm{Var}(S)$ as $\tilde{S} \rightarrow S_1$.

      When computing the distribution, $p(\delta,s)$, of trajectories at variance $S$, given that the first crossed the barrier,
      $B(\tilde{S})$, at some smaller variance, $\tilde{S}$, (equation A2 of \citealt{benson_dark_2012}) we must condition
      \emph{both} the residual variance and drift term on the intermediate point, $\tilde{S},B(\tilde{S})$. Fortunately, given any
      Brownian random walk (including Brownian bridges) for which we know two points, the distribution of trajectories between
      those points is simply another Brownian bridge. Therefore, we can write:      
      \begin{equation}
      \Delta \delta = \delta - \tilde{\delta} - \frac{S-\tilde{S}}{S_2-\tilde{S}}(\delta_2-\tilde{\delta}),
      \end{equation}
      and:
      \begin{equation}
      \Delta S = \frac{(S_2-S)(S-\tilde{S})}{S_2-\tilde{S}}.
      \end{equation}
      
      This class provides functions implementing these modified effective offset and residual variance.
    </description>
  </excursionSetFirstCrossing>
  !!]
  type, extends(excursionSetFirstCrossingFarahiMidpoint) :: excursionSetFirstCrossingFarahiMidpointBrownianBridge
     !!{
     An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method
     to perform the integrations \citep{du_substructure_2017}, and with a Brownian bridge constraint.
     !!}
     private
     class           (excursionSetFirstCrossingClass), pointer :: excursionSetFirstCrossing_     => null()
     class           (linearGrowthClass             ), pointer :: linearGrowth_                  => null()
     class           (criticalOverdensityClass      ), pointer :: criticalOverdensity_           => null()
     double precision                                          :: criticalOverdensityConstrained          , varianceConstrained, &
          &                                                       timeConstrained                         , massConstrained    , &
          &                                                       redshiftConstrained
   contains
     final     ::                     farahiMidpointBrownianBridgeDestructor
     procedure :: rate             => farahiMidpointBrownianBridgeRate
     procedure :: rateNonCrossing  => farahiMidpointBrownianBridgeRateNonCrossing
     procedure :: varianceLimit    => farahiMidpointBrownianBridgeVarianceLimit
     procedure :: varianceResidual => farahiMidpointBrownianBridgeVarianceResidual
     procedure :: offsetEffective  => farahiMidpointBrownianBridgeOffsetEffective
     procedure :: fileWrite        => farahiMidpointBrownianBridgeFileWrite
  end type excursionSetFirstCrossingFarahiMidpointBrownianBridge

  interface excursionSetFirstCrossingFarahiMidpointBrownianBridge
     !!{
     Constructors for the Farahi-midpoint Brownian bridge excursion set barrier class.
     !!}
     module procedure farahiMidpointBrownianBridgeConstructorParameters
     module procedure farahiMidpointBrownianBridgeConstructorInternal
  end interface excursionSetFirstCrossingFarahiMidpointBrownianBridge

contains

  function farahiMidpointBrownianBridgeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the Farahi-midpoint excursion set class first crossing class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (excursionSetFirstCrossingFarahiMidpointBrownianBridge)                :: self
    type            (inputParameters                                      ), intent(inout) :: parameters
    double precision                                                                       :: criticalOverdensityConstrained, varianceConstrained, &
         &                                                                                    timeConstrained               , massConstrained    , &
         &                                                                                    timePresent                   , redshiftConstrained, &
         &                                                                                    expansionFactor

    self%excursionSetFirstCrossingFarahiMidpoint=excursionSetFirstCrossingFarahiMidpoint(parameters)
    !![
    <objectBuilder class="linearGrowth"              name="self%linearGrowth_"              source="parameters"/>
    <objectBuilder class="criticalOverdensity"       name="self%criticalOverdensity_"       source="parameters"/>
    <objectBuilder class="excursionSetFirstCrossing" name="self%excursionSetFirstCrossing_" source="parameters"/>
    !!]
    timePresent=self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    if      (parameters%isPresent('criticalOverdensityConstrained')) then
       if     (                                                                                                                                                                       &
            &  .not.parameters%isPresent('varianceConstrained')                                                                                                                       &
            & ) call Error_Report('both "criticalOverdensityConstrained" and "varianceConstrained" must be provided'                                      //{introspection:location})
       if     (                                                                                                                                                                       &
            &       parameters%isPresent('redshiftConstrained'           )                                                                                                            &
            &  .or.                                                                                                                                                                   &
            &       parameters%isPresent('massConstrained'               )                                                                                                            &
            & ) call Error_Report('can not mix "criticalOverdensityConstrained/varianceConstrained" and "redshiftConstrained/massConstrained" constraints'//{introspection:location})
       !![
       <inputParameter>
         <name>criticalOverdensityConstrained</name>
         <source>parameters</source>
         <description>The critical overdensity at the end of the Brownian bridge.</description>
       </inputParameter>
       <inputParameter>
         <name>varianceConstrained</name>
         <source>parameters</source>
         <description>The variance at the end of the Brownian bridge.</description>
       </inputParameter>
       !!]
       massConstrained=self%cosmologicalMassVariance_%mass          (time               =timePresent                   ,rootVariance=sqrt(varianceConstrained))
       timeConstrained=self%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=criticalOverdensityConstrained,mass        =     massConstrained     )
    else if (parameters%isPresent('redshiftConstrained           ')) then
       if     (                                                                                                                                                                       &
            &  .not.parameters%isPresent('massConstrained'    )                                                                                                                       &
            & ) call Error_Report('both "redshiftConstrained" and "massConstrained" must be provided'                                                     //{introspection:location})
       if     (                                                                                                                                                                       &
            &       parameters%isPresent('criticalOverdensityConstrained')                                                                                                            &
            &  .or.                                                                                                                                                                   &
            &       parameters%isPresent('varianceConstrained'           )                                                                                                            &
            & ) call Error_Report('can not mix "criticalOverdensityConstrained/varianceConstrained" and "redshiftConstrained/massConstrained" constraints'//{introspection:location})
       !![
       <inputParameter>
         <name>redshiftConstrained</name>
         <source>parameters</source>
         <description>The redshift at the end of the Brownian bridge.</description>
       </inputParameter>
       <inputParameter>
         <name>massConstrained</name>
         <source>parameters</source>
         <description>The halo mass at the end of the Brownian bridge.</description>
       </inputParameter>
       !!]
       expansionFactor               =+self%cosmologyFunctions_      %expansionFactorFromRedshift(redshift       =redshiftConstrained                 )
       timeConstrained               =+self%cosmologyFunctions_      %cosmicTime                 (expansionFactor=expansionFactor                     )
       criticalOverdensityConstrained=+self%criticalOverdensity_     %value                      (time           =timeConstrained,mass=massConstrained)    &
            &                         /self%linearGrowth_            %value                      (time           =timeConstrained                     )
       varianceConstrained           =+self%cosmologicalMassVariance_%rootVariance               (time           =timePresent    ,mass=massConstrained)**2
    else
       criticalOverdensityConstrained=0.0d0
       varianceConstrained           =0.0d0
       timeConstrained               =0.0d0
       massConstrained               =0.0d0
       call Error_Report('must provide either [criticalOverdensityConstrained] and [varianceConstrained], or [timeConstrained] and [massConstrained]')
    end if
    self%criticalOverdensityConstrained=criticalOverdensityConstrained
    self%varianceConstrained           =varianceConstrained
    self%timeConstrained               =timeConstrained
    self%massConstrained               =massConstrained
    self%redshiftConstrained           =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeConstrained))    
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function farahiMidpointBrownianBridgeConstructorParameters

  function farahiMidpointBrownianBridgeConstructorInternal(varianceConstrained,criticalOverdensityConstrained,fractionalTimeStep,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,varianceNumberPerDecadeNonCrossing,timeNumberPerDecade,varianceIsUnlimited,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_,criticalOverdensity_,linearGrowth_,excursionSetFirstCrossing_) result(self)
    !!{
    Internal constructor for the Farahi-midpoint excursion set class first crossing class.
    !!}
    implicit none
    type            (excursionSetFirstCrossingFarahiMidpointBrownianBridge)                        :: self
    double precision                                                       , intent(in   )         :: varianceConstrained               , criticalOverdensityConstrained, &
         &                                                                                            fractionalTimeStep
    integer                                                                , intent(in   )         :: varianceNumberPerUnitProbability  , varianceNumberPerUnit         , &
         &                                                                                            timeNumberPerDecade               , varianceNumberPerDecade       , &
         &                                                                                            varianceNumberPerDecadeNonCrossing
    logical                                                                , intent(in   )         :: varianceIsUnlimited
    type            (varying_string                                       ), intent(in   )         :: fileName
    class           (cosmologyFunctionsClass                              ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetBarrierClass                             ), intent(in   ), target :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass                        ), intent(in   ), target :: cosmologicalMassVariance_
    class           (linearGrowthClass                                    ), intent(in   ), target :: linearGrowth_
    class           (criticalOverdensityClass                             ), intent(in   ), target :: criticalOverdensity_
    class           (excursionSetFirstCrossingClass                       ), intent(in   ), target :: excursionSetFirstCrossing_
    double precision                                                                               :: timePresent                     , expansionFactor
    !![
    <constructorAssign variables="varianceConstrained, criticalOverdensityConstrained, *criticalOverdensity_, *linearGrowth_, *excursionSetFirstCrossing_"/>
    !!]
    
    self%excursionSetFirstCrossingFarahiMidpoint=excursionSetFirstCrossingFarahiMidpoint(fractionalTimeStep,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,varianceNumberPerDecadeNonCrossing,timeNumberPerDecade,varianceIsUnlimited,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_)
    ! Find mass and time corresponding to the constraint point.
    timePresent             =self%cosmologyFunctions_      %cosmicTime                 (expansionFactor    =1.0d0                                                                          )
    self%massConstrained    =self%cosmologicalMassVariance_%mass                       (time               =timePresent                        ,rootVariance=sqrt(self%varianceConstrained))
    self%timeConstrained    =self%criticalOverdensity_     %timeOfCollapse             (criticalOverdensity=self%criticalOverdensityConstrained,mass        =     self%massConstrained     )
    expansionFactor         =self%cosmologyFunctions_      %expansionFactor            (time               =self%timeConstrained                                                           )
    self%redshiftConstrained=self%cosmologyFunctions_      %redshiftFromExpansionFactor(expansionFactor    =expansionFactor                                                                )
    return
  end function farahiMidpointBrownianBridgeConstructorInternal

  subroutine farahiMidpointBrownianBridgeDestructor(self)
    !!{
    Destructor for the piecewise Farahi excursion Brownian bridge set first crossing class.
    !!}
    implicit none
    type(excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self

    !![
    <objectDestructor name="self%linearGrowth_"             />
    <objectDestructor name="self%criticalOverdensity_"      />
    <objectDestructor name="self%excursionSetFirstCrossing_"/>
    !!]
    return
  end subroutine farahiMidpointBrownianBridgeDestructor

  double precision function farahiMidpointBrownianBridgeRate(self,variance,varianceProgenitor,time,node)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    double precision                                                       , intent(in   ) :: variance                      , varianceProgenitor , &
         &                                                                                    time
    type            (treeNode                                             ), intent(inout) :: node
    double precision                                                                       :: rootVarianceConstrained       , varianceConstrained, &
         &                                                                                    criticalOverdensityConstrained

    ! For progenitor variances less than or equal to the original variance, return zero.
    if (varianceProgenitor <= variance) then
       farahiMidpointBrownianBridgeRate=0.0d0
       return
    end if
    ! Determine effective constraint point at this epoch.
    rootVarianceConstrained       =+self%cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained           =+rootVarianceConstrained**2
    criticalOverdensityConstrained=+self%criticalOverdensityConstrained             &
         &                         *sqrt(                                           &
         &                               +     varianceConstrained                  &
         &                               /self%varianceConstrained                  &
         &                              )
    ! Determine whether to use the conditioned or unconditioned solutions.
    if (self%excursionSetBarrier_%barrier(varianceProgenitor,time,node,rateCompute=.true.) > criticalOverdensityConstrained) then
       ! The time corresponds to a barrier above the constrained point. Therefore we want the unconstrained solution.
       farahiMidpointBrownianBridgeRate=self%excursionSetFirstCrossing_%rate           (variance,varianceProgenitor,time,node)
    else if (varianceProgenitor >= varianceConstrained) then
       ! For progenitor variances in excess of the constrained variance the first crossing rate must be zero.
       farahiMidpointBrownianBridgeRate=0.0d0
    else
       ! Use the constrained solution.
       farahiMidpointBrownianBridgeRate=self                           %rateInterpolate(variance,varianceProgenitor,time,node)
    end if
    return
  end function farahiMidpointBrownianBridgeRate

  double precision function farahiMidpointBrownianBridgeRateNonCrossing(self,variance,massMinimum,time,node)
    !!{
    Return the rate for excursion set non-crossing.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    double precision                                                       , intent(in   ) :: time                          , variance           , &
         &                                                                                    massMinimum
    type            (treeNode                                             ), intent(inout) :: node
    double precision                                                                       :: rootVarianceConstrained       , varianceConstrained, &
         &                                                                                    criticalOverdensityConstrained

    ! Determine effective constraint point at this epoch.
    rootVarianceConstrained       =+self%cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained           =+rootVarianceConstrained**2
    criticalOverdensityConstrained=+self%criticalOverdensityConstrained             &
         &                         *sqrt(                                           &
         &                               +     varianceConstrained                  &
         &                               /self%varianceConstrained                  &
         &                              )
    ! Determine whether to use the conditioned or unconditioned solutions.
    if (self%excursionSetBarrier_%barrier(variance,time,node,rateCompute=.true.) > criticalOverdensityConstrained) then
       ! The time corresponds to a barrier above the constrained point. Therefore we want the unconstrained solution.
       farahiMidpointBrownianBridgeRateNonCrossing=self%excursionSetFirstCrossing_%rateNonCrossing           (variance,     massMinimum    ,time,node)
    else
       ! Use the constrained solution. By definition, all trajectories cross the barrier before the constrained point. Therefore,
       ! the non-crossing rate is always zero.
       farahiMidpointBrownianBridgeRateNonCrossing=0.0d0
    end if
    return
  end function farahiMidpointBrownianBridgeRateNonCrossing
  
  double precision function farahiMidpointBrownianBridgeVarianceLimit(self,varianceProgenitor)
    !!{
    Return the maximum variance to which to tabulate. For the case of a Brownian bridge the variance must not be allowed to exceed
    the variance $S_2$ at the end of the bridge---all trajectories must have crossed the barrier by this variance by construction.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    double precision                                                       , intent(in   ) :: varianceProgenitor

    farahiMidpointBrownianBridgeVarianceLimit=min(                                     &
         &                                        max(                                 &
         &                                                   self%varianceMaximumRate, &
         &                                             2.0d0*     varianceProgenitor   &
         &                                            )                              , &
         &                                                   self%varianceConstrained  &
         &                                       )
    return
  end function farahiMidpointBrownianBridgeVarianceLimit

  function farahiMidpointBrownianBridgeVarianceResidual(self,time,varianceCurrent,varianceProgenitor,varianceIntermediate,cosmologicalMassVariance_) result(varianceResidual)
    !!{
    Return the residual variance between two points for a Brownian bridge.
    !!}
    use :: Kind_Numbers, only : kind_quad
    implicit none
    real (kind_quad                                            )                :: varianceResidual
    class(excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    real (kind_quad                                            ), intent(in   ) :: varianceCurrent          , varianceIntermediate, &
         &                                                                         varianceProgenitor
    double precision                                            , intent(in   ) :: time
    class           (cosmologicalMassVarianceClass             ), intent(inout) :: cosmologicalMassVariance_
    real (kind_quad                                            )                :: rootVarianceConstrained  , varianceConstrained
    
    ! Note that this solver follows the convention used through Galacticus that σ(M) grows following linear theory. That is:
    !
    !  • the root-variance of the density field smoothed on a mass scale M is a function of time, σ(M,t) = σ(M,t₀) D(t)/D(t₀),
    !    where D(t) is the linear growth factor (which may also be scale-dependent);
    !  • the critical overdensity for collapse does not include a factor of the linear growth factor, i.e. δ_c ≅ 1.686 at all
    !    epochs (varying only due to the weak dependence on the epoch-dependent cosmological parameters).
    !
    ! This differs from standard treatments of the excursion set problem in which typically the root-variance, σ(M), is evaluated
    ! at z=0, and the critical overdensity for collapse is replaced with δ_c(t)/D(t). Mathematically these two approaches are
    ! equivalent, but it can be important to keep these distinctions in mind.

    ! In this function the following translations between internal variable names and math symbols are used:
    !
    !   S₂ = varianceConstrained
    !   S₁ = varianceCurrent
    !   S̃ = varianceIntermediate+varianceCurrent
    !   S  = varianceProgenitor  +varianceCurrent
    !    
    ! Note that the variables "varianceIntermediate" and "varianceProgenitor" are defined to be the variances in excess of S₁ - which is why they
    ! appear with "varianceCurrent" added to them in the above.
    !
    ! This function is used in the calculation of the distribution of δ at some S for trajectories originating from (S₁,δ₁) and
    ! which did not cross the barrier at any intermediate variance. As such suffixes in variable names have the following
    ! meanings:
    !
    !   "Current"      - refers to the current halo being considered for branching, i.e. the halo existing at point (S₁,δ₁);
    !   "Progenitor"   - refers to the potential progenitor halo being considered, i.e. the halo corresponding to some variance S > S₁;
    !   "Intermediate" - refers to the intermediate variance, S̃ (with S₁ < S̃ < S).

    ! Find the variance corresponding to the end of the bridge at the current time. Note that this is necessary because the
    ! variances used throughout the excursion set solver code grow in time according to linear perturbation theory. The fixed
    ! quantity in our Brownian bridge is the mass and time of the halo at the end of the bridge. Therefore, we must compute the
    ! variance corresponding to this mass at the requested epoch.
    rootVarianceConstrained=+cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained    =+rootVarianceConstrained**2
    ! Compute the residual variance.
    varianceResidual       =+(varianceConstrained-varianceProgenitor  -varianceCurrent) &
         &                  *(varianceProgenitor -varianceIntermediate                ) &
         &                  /(varianceConstrained-varianceIntermediate-varianceCurrent)
    return
  end function farahiMidpointBrownianBridgeVarianceResidual

  function farahiMidpointBrownianBridgeOffsetEffective(self,time,varianceCurrent,varianceProgenitor,varianceIntermediate,deltaCurrent,deltaProgenitor,deltaIntermediate,cosmologicalMassVariance_) result(offsetEffective)
    !!{
    Return the residual variance between two points for a Brownian bridge.
    !!}
    use :: Kind_Numbers, only : kind_quad
    implicit none
    real            (kind_quad                                            )                :: offsetEffective
    class           (excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    real            (kind_quad                                            ), intent(in   ) :: deltaCurrent                  , deltaIntermediate  , &
         &                                                                                    deltaProgenitor               , varianceCurrent    , &
         &                                                                                    varianceIntermediate          , varianceProgenitor
    double precision                                                       , intent(in   ) :: time
    class           (cosmologicalMassVarianceClass                        ), intent(inout) :: cosmologicalMassVariance_
    real            (kind_quad                                            )                :: rootVarianceConstrained       , varianceConstrained, &
         &                                                                                    criticalOverdensityConstrained

    ! Note that this solver follows the convention used through Galacticus that σ(M) grows following linear theory. That is:
    !
    !  • the root-variance of the density field smoothed on a mass scale M is a function of time, σ(M,t) = σ(M,t₀) D(t)/D(t₀),
    !    where D(t) is the linear growth factor (which may also be scale-dependent);
    !  • the critical overdensity for collapse does not include a factor of the linear growth factor, i.e. δ_c ≅ 1.686 at all
    !    epochs (varying only due to the weak dependence on the epoch-dependent cosmological parameters).
    !
    ! This differs from standard treatments of the excursion set problem in which typically the root-variance, σ(M), is evaluated
    ! at z=0, and the critical overdensity for collapse is replaced with δ_c(t)/D(t). Mathematically these two approaches are
    ! equivalent, but it can be important to keep these distinctions in mind.

    ! In this function the following translations between internal variable names and math symbols are used:
    !
    !   S₂ = varianceConstrained
    !   S₁ = varianceCurrent
    !   S̃ = varianceIntermediate+varianceCurrent
    !   S  = varianceProgenitor  +varianceCurrent
    !   δ₂ = criticalOverdensityConstrained
    !   δ₁ = deltaCurrent
    !   δ̃  = deltaIntermediate   +deltaCurrent
    !   δ  = deltaProgenitor     +deltaCurrent
    !    
    ! Note that the variables "varianceIntermediate" and "varianceProgenitor" are defined to be the variances in excess of S₁ - which is why they
    ! appear with "varianceCurrent" added to them in the above. Similarly, "deltaIntermediate" and "deltaProgenitor" are defined relative to "deltaCurrent".
    !
    ! This function is used in the calculation of the distribution of δ at some S for trajectories originating from (S₁,δ₁) and
    ! which did not cross the barrier at any intermediate variance. As such suffixes in variable names have the following
    ! meanings:
    !
    !   "Current"      - refers to the current halo being considered for branching, i.e. the halo existing at point (S₁,δ₁);
    !   "Progenitor"   - refers to the potential progenitor halo being considered, i.e. the halo corresponding to some variance S > S₁;
    !   "Intermediate" - refers to the intermediate variance, S̃ (with S₁ < S̃ < S).

    ! Find the variance corresponding to the end of the bridge at the current time. Note that this is necessary because the
    ! variances used throughout the excursion set solver code grow in time according to linear perturbation theory. The fixed
    ! quantity in our Brownian bridge is the mass and time of the halo at the end of the bridge. Therefore, we must compute the
    ! variance corresponding to this mass at the requested epoch.
    !
    ! Also note that the drift term is always computed from the difference between the constraint point and the intermediate
    ! point. This is because, for a trajectory having arrived at the intermediate point, the set of all possible trajectories that
    ! take it from there to the constrained point is precisely a Brownian bridge starting from the intermediate point and ending
    ! at the constraint point.
    rootVarianceConstrained       =+cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained           =+rootVarianceConstrained**2
    criticalOverdensityConstrained=+self%criticalOverdensityConstrained                    &
         &                         *sqrt(                                                  &
         &                               +     varianceConstrained                         &
         &                               /self%varianceConstrained                         &
         &                              )
    offsetEffective               =+(+deltaProgenitor               -deltaIntermediate                   ) &
         &                         -(+criticalOverdensityConstrained-deltaIntermediate   -deltaCurrent   ) &
         &                         *(+varianceProgenitor            -varianceIntermediate                ) &
         &                         /(+varianceConstrained           -varianceIntermediate-varianceCurrent)
    return
  end function farahiMidpointBrownianBridgeOffsetEffective
 
  subroutine farahiMidpointBrownianBridgeFileWrite(self)
    !!{
    Write additional data on excursion set first crossing probabilities to file for the case of the Brownian bridge. Specifically,
    linear growth factors are written to the file as a convenience useful for interpreting the results.
    !!}
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout)               :: self
    type            (hdf5Object                                           )                              :: dataFile          , dataGroup
    double precision                                                       , allocatable  , dimension(:) :: linearGrowthFactor
    integer                                                                                              :: i
    
    ! Write the primary data.
    call self%excursionSetFirstCrossingFarahiMidpoint%fileWrite()
    ! Don't write anything if neither table is initialized.
    if (.not.self%tableInitializedRate) return
    ! Open the data file.
    !$ call hdf5Access%set()
    call dataFile%openFile(char(self%fileName),overWrite=.false.)
    ! Check if the rate table is populated.
    if (self%tableInitializedRate) then
       allocate(linearGrowthFactor(size(self%timeRate)))
       do i=1,size(self%timeRate)
          linearGrowthFactor(i)=self%linearGrowth_%value(time=self%timeRate(i))
       end do
       dataGroup=dataFile%openGroup("rate")
       call dataGroup%writeDataset(linearGrowthFactor,'linearGrowthFactor','The linear growth factors at the times at which results are tabulated.')
       call dataGroup%close()
    end if
    ! Close the data file.
    call dataFile%close()
    !$ call hdf5Access%unset()
    return
  end subroutine farahiMidpointBrownianBridgeFileWrite
