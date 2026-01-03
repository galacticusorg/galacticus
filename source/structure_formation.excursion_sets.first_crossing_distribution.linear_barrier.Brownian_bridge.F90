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

!!{
Implements a excursion set first crossing statistics class for linear barriers with constrained branching
described by a Brownian bridge solution.
!!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Linear_Growth             , only : linearGrowthClass 
  use :: Excursion_Sets_Barriers   , only : excursionSetBarrierClass

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingLinearBarrierBrownianBridge">
   <description>
     An excursion set first crossing statistics class for linear barriers and where the trajectories are constrained to follow a
     Brownian bridge.

     If we consider the Brownian bridge to originate from $(0,0)$ (i.e. we apply the usual shift of coordinates to move our
     starting point to the origin), and to end at $(\delta_2,S_2)$ then we can transform this Brownian bridge into the standard
     bridge with non-zero drift through the transformations:    
     \begin{eqnarray}
     \tau &amp;=&amp; \frac{S}{S_2}, \\
     b    &amp;=&amp; \frac{\delta_2}{\sqrt{S_2}}.
     \end{eqnarray}
     To find the first crossing time distribution we then follow the general approach outlined by \cite{kiwiakos_answer_2014}, but
     with an important difference that we will detail below.

     The standard Brownian bridge (with no drift), $Y_0$, can be written in terms of a standard Weiner process, $W$, through a
     change of variables     
     \begin{equation}
     Y_0(t) = (1-t) W\left(\frac{t}{1-t}\right).
     \end{equation}
     The first crossing time distribution for our Brownian bridge can therefore be expressed as:
     \begin{equation}
     \tau_{Y}(B) = \mathrm{inf}\left\{ t : Y(t) = B(t)\right\} = \mathrm{inf}\left\{ \mu(t) + (1-t) W\left(\frac{t}{1-t}\right) \right\} = B(t) = \mathrm{inf}\left\{ t : W\left(\frac{t}{1-t}\right) = B(t) - \mu(t) \right\},
     \end{equation}
     where $B(t)$ is our barrier, and $\mu(t)$ is the drift term in the Brownian bridge.

     As can be seen from the above, for the case of a linear barrier, a Brownian bridge with non-zero drift effectively results in
     a new linear barrier equal to the original one minus the drift term, i.e.:
     \begin{equation}
     B^\prime(t) \rightarrow B(t) - \mu(t),
     \end{equation}
     where $B(S)$ is the barrier and $\mu(S)$ is the Brownian bridge term. This means that the first-crossing time of the
     Brownian bridge is just the hitting time of this time changed Weiner process. That is, is
     \begin{equation}
     \tau_W(B) = \mathrm{inf}\left\{ s : W(s) = B\right\},
     \end{equation}
     then
     \begin{equation}
     \frac{\tau_Y(B)}{1 - \tau_Y(B)} = \tau_W(B) \implies \tau_Y(B) = \frac{\tau_W(B)}{1+\tau_W(B)}.
     \end{equation}
     Here is where the solution presented by \cite{kiwiakos_answer_2014} is slightly wrong. We must use the first crossing time
     solution for the standard Weiner process, but with a \emph{linear} barrier (because, even if the actual barrier is constant,
     the effective barrier is linear due to the Brownian bridge drift term). Therefore (e.g. \citealt{zhang_random_2006}):
     \begin{equation}
     f_{W}(\tau_{W}) = B(0) \exp\left( - \frac{B(\tau_{W})^2}{2\ tau_{W}} \right) / \sqrt{2 pi \tau_{W}^3}.
     \end{equation}
     We then have that
     \begin{equation}
     f_{Y}(\tau_{Y}) \mathrm{d}\tau_{Y} = f_{W}(\tau_{W}) \mathrm{d}\tau_{W},
     \end{equation}
     such that
     \begin{equation}
     f_{Y}(\tau_{Y}) = \frac{B(0)}{\sqrt{2 \pi \tau_{Y}^3 (1 - \tau_{Y})}} \exp\left( \frac{B^{\prime 2}(\tau_{Y})}{2 \tau_{Y} (1-\tau_{Y})} \right),
     \end{equation}
     or, expressed in our usual variables
     \begin{equation}
     f(S) = \frac{S(0)}{\sqrt{2 \pi S^3 (1 - S/S_2)}} \exp\left( \frac{[B(S)-\mu(S)]^2}{2 S (1-S/S_2)} \right).
     \end{equation}
   </description>
  </excursionSetFirstCrossing>
  !!]
  type, extends(excursionSetFirstCrossingClass) :: excursionSetFirstCrossingLinearBarrierBrownianBridge
     !!{
     A linearBarrierBrownianBridge excursion set barrier class.
     !!}
     private
     class           (cosmologyFunctionsClass       ), pointer :: cosmologyFunctions_            => null()
     class           (excursionSetBarrierClass      ), pointer :: excursionSetBarrier_           => null()
     class           (excursionSetFirstCrossingClass), pointer :: excursionSetFirstCrossing_     => null()
     class           (cosmologicalMassVarianceClass ), pointer :: cosmologicalMassVariance_      => null()
     class           (linearGrowthClass             ), pointer :: linearGrowth_                  => null()
     class           (criticalOverdensityClass      ), pointer :: criticalOverdensity_           => null()
     double precision                                          :: criticalOverdensityConstrained          , varianceConstrained, &
          &                                                       timeConstrained                         , massConstrained    , &
          &                                                       redshiftConstrained
     ! The fractional step in time used to compute barrier crossing rates.
     double precision                                          :: fractionalTimeStep
   contains
     final     ::                    linearBarrierBrownianBridgeDestructor
     procedure :: probability     => linearBarrierBrownianBridgeProbability
     procedure :: rate            => linearBarrierBrownianBridgeRate
     procedure :: rateNonCrossing => linearBarrierBrownianBridgeRateNonCrossing
  end type excursionSetFirstCrossingLinearBarrierBrownianBridge

  interface excursionSetFirstCrossingLinearBarrierBrownianBridge
     !!{
     Constructors for the linearBarrierBrownianBridge excursion set barrier class.
     !!}
     module procedure linearBarrierBrownianBridgeConstructorParameters
     module procedure linearBarrierBrownianBridgeConstructorInternal
  end interface excursionSetFirstCrossingLinearBarrierBrownianBridge

contains

  function linearBarrierBrownianBridgeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (excursionSetFirstCrossingLinearBarrierBrownianBridge)                :: self
    type            (inputParameters                                     ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                             ), pointer       :: cosmologyFunctions_
    class           (excursionSetFirstCrossingClass                      ), pointer       :: excursionSetFirstCrossing_
    class           (excursionSetBarrierClass                            ), pointer       :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass                       ), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass                                   ), pointer       :: linearGrowth_
    class           (criticalOverdensityClass                            ), pointer       :: criticalOverdensity_
    double precision                                                                      :: criticalOverdensityConstrained, varianceConstrained, &
         &                                                                                   timeConstrained               , massConstrained    , &
         &                                                                                   timePresent                   , redshiftConstrained, &
         &                                                                                   expansionFactor               , fractionalTimeStep
    
    !![
    <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="linearGrowth"              name="linearGrowth_"              source="parameters"/>
    <objectBuilder class="criticalOverdensity"       name="criticalOverdensity_"       source="parameters"/>
    <objectBuilder class="excursionSetBarrier"       name="excursionSetBarrier_"       source="parameters"/>
    <objectBuilder class="excursionSetFirstCrossing" name="excursionSetFirstCrossing_" source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"  name="cosmologicalMassVariance_"  source="parameters"/>
    <inputParameter>
      <name>fractionalTimeStep</name>
      <defaultValue>0.01d0</defaultValue>
      <source>parameters</source>
      <description>The fractional time step used when computing barrier crossing rates (i.e. the step used in finite difference calculations).</description>
    </inputParameter>
    !!]
    timePresent=cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
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
       massConstrained=cosmologicalMassVariance_%mass          (time               =timePresent                   ,rootVariance=sqrt(varianceConstrained))
       timeConstrained=criticalOverdensity_     %timeOfCollapse(criticalOverdensity=criticalOverdensityConstrained,mass        =     massConstrained     )
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
       expansionFactor               =+cosmologyFunctions_      %expansionFactorFromRedshift(redshift       =redshiftConstrained                 )
       timeConstrained               =+cosmologyFunctions_      %cosmicTime                 (expansionFactor=expansionFactor                     )
       criticalOverdensityConstrained=+criticalOverdensity_     %value                      (time           =timeConstrained,mass=massConstrained)    &
            &                         /linearGrowth_            %value                      (time           =timeConstrained                     )
       varianceConstrained           =+cosmologicalMassVariance_%rootVariance               (time           =timePresent    ,mass=massConstrained)**2
    else
       criticalOverdensityConstrained=0.0d0
       varianceConstrained           =0.0d0
       timeConstrained               =0.0d0
       massConstrained               =0.0d0
       call Error_Report('must provide either [criticalOverdensityConstrained] and [varianceConstrained], or [timeConstrained] and [massConstrained]')
    end if
    self=excursionSetFirstCrossingLinearBarrierBrownianBridge(varianceConstrained,criticalOverdensityConstrained,fractionalTimeStep,cosmologyFunctions_,excursionSetFirstCrossing_,excursionSetBarrier_,cosmologicalMassVariance_,criticalOverdensity_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="excursionSetBarrier_"      />
    <objectDestructor name="cosmologicalMassVariance_" />
    <objectDestructor name="cosmologyFunctions_"       />
    <objectDestructor name="excursionSetFirstCrossing_"/>
    <objectDestructor name="criticalOverdensity_"      />
    <objectDestructor name="linearGrowth_"             />
    !!]
    return
  end function linearBarrierBrownianBridgeConstructorParameters

  function linearBarrierBrownianBridgeConstructorInternal(varianceConstrained,criticalOverdensityConstrained,fractionalTimeStep,cosmologyFunctions_,excursionSetFirstCrossing_,excursionSetBarrier_,cosmologicalMassVariance_,criticalOverdensity_,linearGrowth_) result(self)
    !!{
    Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    !!}
    implicit none
    type            (excursionSetFirstCrossingLinearBarrierBrownianBridge)                        :: self
    double precision                                                      , intent(in   )         :: varianceConstrained       , criticalOverdensityConstrained, &
         &                                                                                           fractionalTimeStep
    class           (cosmologyFunctionsClass                             ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetFirstCrossingClass                      ), intent(in   ), target :: excursionSetFirstCrossing_
    class           (excursionSetBarrierClass                            ), intent(in   ), target :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass                       ), intent(in   ), target :: cosmologicalMassVariance_
    class           (linearGrowthClass                                   ), intent(in   ), target :: linearGrowth_
    class           (criticalOverdensityClass                            ), intent(in   ), target :: criticalOverdensity_
    double precision                                                                              :: timePresent               , expansionFactor
     !![
    <constructorAssign variables="varianceConstrained, criticalOverdensityConstrained, fractionalTimeStep, *cosmologyFunctions_, *excursionSetFirstCrossing_, *excursionSetBarrier_, *cosmologicalMassVariance_, *linearGrowth_, *criticalOverdensity_"/>
    !!]

    ! Find mass and time corresponding to the constraint point.
    timePresent             =self%cosmologyFunctions_      %cosmicTime                 (expansionFactor    =1.0d0                                                                          )
    self%massConstrained    =self%cosmologicalMassVariance_%mass                       (time               =timePresent                        ,rootVariance=sqrt(self%varianceConstrained))
    self%timeConstrained    =self%criticalOverdensity_     %timeOfCollapse             (criticalOverdensity=self%criticalOverdensityConstrained,mass        =     self%massConstrained     )
    expansionFactor         =self%cosmologyFunctions_      %expansionFactor            (time               =self%timeConstrained                                                           )
    self%redshiftConstrained=self%cosmologyFunctions_      %redshiftFromExpansionFactor(expansionFactor    =expansionFactor                                                                )
    return
  end function linearBarrierBrownianBridgeConstructorInternal

  subroutine linearBarrierBrownianBridgeDestructor(self)
    !!{
    Destructor for the critical overdensity excursion set barrier class.
    !!}
    implicit none
    type(excursionSetFirstCrossingLinearBarrierBrownianBridge), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"       />
    <objectDestructor name="self%excursionSetFirstCrossing_"/>
    <objectDestructor name="self%excursionSetBarrier_"      />
    <objectDestructor name="self%cosmologicalMassVariance_" />
    <objectDestructor name="self%linearGrowth_"             />
    <objectDestructor name="self%criticalOverdensity_"      />
    !!]
    return
  end subroutine linearBarrierBrownianBridgeDestructor

  double precision function linearBarrierBrownianBridgeProbability(self,variance,time,node) result(probability)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (excursionSetFirstCrossingLinearBarrierBrownianBridge), intent(inout) :: self
    double precision                                                      , intent(in   ) :: variance           , time
    type            (treeNode                                            ), intent(inout) :: node
    double precision                                                                      :: varianceConstrained, criticalOverdensityConstrained, &
         &                                                                                   barrierConstrained , rootVarianceConstrained
    
    ! Determine effective constraint point at this epoch.
    rootVarianceConstrained       =+self%cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained           =+rootVarianceConstrained**2
    criticalOverdensityConstrained=+self%criticalOverdensityConstrained &
         &                         *sqrt(                               &
         &                               +     varianceConstrained      &
         &                               /self%varianceConstrained      &
         &                              )
    ! Determine whether to use the conditioned or unconditioned solutions.
    if      (self%excursionSetBarrier_%barrier(variance,time,node,rateCompute=.true.) > criticalOverdensityConstrained) then
       ! The time corresponds to a barrier above the constrained point. Therefore we want the unconstrained solution.
       probability       =+self%excursionSetFirstCrossing_%probability(variance,time,node)
    else if (                                  variance                               >= varianceConstrained          ) then
       ! For variances in excess of the constrained variance the first crossing probability must be zero.
       probability       =+0.0d0
    else
       barrierConstrained=+(                                                                                                    &
            &               +  self%excursionSetBarrier_%barrier                       (variance,time,node,rateCompute=.false.) &
            &               -(                                                                                                  &
            &                 +                          criticalOverdensityConstrained                                         &
            &                 -self%excursionSetBarrier_%barrier                       (variance,time,node,rateCompute=.false.) &
            &                )                                                                                                  &
            &               *variance                                                                                           &
            &               /varianceConstrained                                                                                &
            &              )                                                                                                    &
            &             /(                                                                                                    &
            &               +1.0d0                                                                                              &
            &               -variance                                                                                           &
            &               /varianceConstrained                                                                                &
            &             )
       probability       =+     self%excursionSetBarrier_%barrier(   0.0d0,time,node,rateCompute=.false.)                          &
            &             *exp(                                                                                                    &
            &                  -0.5d0                                                                                              &
            &                  *(                                                                                                  &
            &                    +self%excursionSetBarrier_%barrier                       (variance,time,node,rateCompute=.false.) &
            &                    -                          criticalOverdensityConstrained                                         &
            &                    *variance                                                                                         &
            &                    /varianceConstrained                                                                              &
            &                   )**2                                                                                               &
            &                  *(                                                                                                  &
            &                    +1.0d0                                                                                            &
            &                    -variance                                                                                         &
            &                    /varianceConstrained                                                                              &
            &                   )                                                                                                  &
            &                   / variance                                                                                         &
            &                  )                                                                                                   &
            &                 /sqrt(                                                                                               &
            &                       +2.0d0                                                                                         &
            &                       *Pi                                                                                            &
            &                       *  variance           **3                                                                      &
            &                       *(                                                                                             &
            &                         +1.0d0                                                                                       &
            &                         -variance                                                                                    &
            &                         /varianceConstrained                                                                         &
            &                        )                                                                                             &
            &                      )
    end if
    return
  end function linearBarrierBrownianBridgeProbability

  double precision function linearBarrierBrownianBridgeRate(self,variance,varianceProgenitor,time,node) result(rate)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (excursionSetFirstCrossingLinearBarrierBrownianBridge), intent(inout) :: self
    double precision                                                      , intent(in   ) :: variance                      , varianceProgenitor     , &
         &                                                                                   time
    type            (treeNode                                            ), intent(inout) :: node
    double precision                                                                      :: timeProgenitor                , massProgenitor         , &
         &                                                                                   growthFactorEffective         , varianceConstrained    , &
         &                                                                                   criticalOverdensityConstrained, varianceDifference     , &
         &                                                                                   varianceConstrainedDifference , rootVarianceConstrained
    
    ! Determine effective constraint point at this epoch.
    rootVarianceConstrained       =+self%cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained           =+rootVarianceConstrained**2
    criticalOverdensityConstrained=+self%criticalOverdensityConstrained &
         &                         *sqrt(                               &
         &                               +     varianceConstrained      &
         &                               /self%varianceConstrained      &
         &                              )
      ! Determine whether to use the conditioned or unconditioned solutions.
    if (self%excursionSetBarrier_%barrier(varianceProgenitor,time,node,rateCompute=.true.) > criticalOverdensityConstrained) then
       ! The time corresponds to a barrier above the constrained point. Therefore we want the unconstrained solution.
       rate=self%excursionSetFirstCrossing_%rate           (variance,varianceProgenitor,time,node)
    else if (                             varianceProgenitor                               >= varianceConstrained          ) then
       ! For progenitor variances in excess of the constrained variance the first crossing rate must be zero.
       rate=0.0d0
    else
       ! * To estimate the rate we use a finite difference method - we compute the effective barrier for a slightly earlier time,
       !   compute the fraction of trajectories which will have crossed that effective barrier, and divide by the time difference.
       !
       ! * In Galacticus, the time evolution due to linear growth is included in the root-variance of the density field, *not* in
       !   the barrier height as is often done in the literature. As such, when computing the barrier at some earlier time we must
       !   account for the fact that, at a fixed mass, the root variance will be smaller at that earlier time. Since the solution
       !   to the excursion set problem must always be a function of δc(M,t)/√S(M,t) then we can simply scale δc by the ratio of
       !   root-variances for the progenitor at the current and earlier times.
       timeProgenitor               =+time*(1.0d0-self%fractionalTimeStep)
       massProgenitor               =+self%cosmologicalMassVariance_%mass        (sqrt(varianceProgenitor),time          )
       growthFactorEffective        =+self%cosmologicalMassVariance_%rootVariance(         massProgenitor ,time          ) &
            &                        /self%cosmologicalMassVariance_%rootVariance(         massProgenitor ,timeProgenitor)
       varianceDifference           =+varianceProgenitor -variance
       varianceConstrainedDifference=+varianceConstrained-variance
       rate                         =+     barrierEffective           (variance,time,variance          ,timeProgenitor)    &
           &                         *exp(                                                                                 &
           &                              -0.5d0                                                                           &
           &                              *barrierEffectiveConstrained(variance,time,varianceProgenitor,timeProgenitor)**2 &
           &                              *(                                                                               &
           &                                +1.0d0                                                                         &
           &                                -varianceDifference                                                            &
           &                                /varianceConstrainedDifference                                                 &
           &                               )                                                                               &
           &                              /  varianceDifference                                                            &
           &                             )                                                                                 &
           &                         /sqrt(                                                                                &
           &                               +2.0d0                                                                          &
           &                               *Pi                                                                             &
           &                               *  varianceDifference           **3                                             &
           &                               *(                                                                              &
           &                                 +1.0d0                                                                        &
           &                                 -varianceDifference                                                           &
           &                                 /varianceConstrainedDifference                                                &
           &                                )                                                                              &
           &                              )                                                                                &
           &                         /time                                                                                 &
           &                         /self%fractionalTimeStep
      rate                          =max(                                                                                  &
           &                             rate                                                                            , &
           &                             0.0d0                                                                             &
           &                            )
   end if
   return

  contains

    double precision function barrierEffective(variance0,time0,variance1,time1)
      !!{
      The effective barrier for conditional excursion sets.
      !!}
      implicit none
      double precision, intent(in   ) :: time1    , time0    , &
           &                             variance1, variance0
      
      barrierEffective=+self%excursionSetBarrier_%barrier(variance1,time1,node,rateCompute=.false.)*growthFactorEffective &
           &           -self%excursionSetBarrier_%barrier(variance0,time0,node,rateCompute=.false.)
      return
    end function barrierEffective

    double precision function barrierEffectiveConstrained(variance0,time0,variance1,time1)
      !!{
      The effective barrier for conditional excursion sets.
      !!}
      implicit none
      double precision, intent(in   ) :: time1    , time0    , &
           &                             variance1, variance0
      
      barrierEffectiveConstrained=+(                                                                                                                            &
           &                        +(                                                                                                                          &
           &                          +self%excursionSetBarrier_%barrier                       (variance1,time1,node,rateCompute=.false.)*growthFactorEffective &
           &                          -self%excursionSetBarrier_%barrier                       (variance0,time0,node,rateCompute=.false.)                       &
           &                         )                                                                                                                          &
           &                        -(                                                                                                                          &
           &                          +                          criticalOverdensityConstrained                                                                 &
           &                          -self%excursionSetBarrier_%barrier                       (variance0,time0,node,rateCompute=.false.)                       &
           &                         )                                                                                                                          &
           &                        *varianceDifference                                                                                                         &
           &                        /varianceConstrainedDifference                                                                                              &
           &                       )                                                                                                                            &
           &                      /(                                                                                                                            &
           &                        +1.0d0                                                                                                                      &
           &                        -varianceDifference                                                                                                         &
           &                        /varianceConstrainedDifference                                                                                              &
           &                       )
      return
    end function barrierEffectiveConstrained

  end function linearBarrierBrownianBridgeRate

  double precision function linearBarrierBrownianBridgeRateNonCrossing(self,variance,massMinimum,time,node) result(rateNonCrossing)
    !!{
    Return the rate for excursion set non-crossing assuming a linear barrier.
    !!}
    implicit none
    class           (excursionSetFirstCrossingLinearBarrierBrownianBridge), intent(inout) :: self
    double precision                                                      , intent(in   ) :: time                          , variance           , &
         &                                                                                   massMinimum
    type            (treeNode                                            ), intent(inout) :: node
    double precision                                                                      :: rootVarianceConstrained       , varianceConstrained, &
         &                                                                                   criticalOverdensityConstrained

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
       rateNonCrossing=self%excursionSetFirstCrossing_%rateNonCrossing(variance,massMinimum,time,node)
    else
       ! Use the constrained solution. By definition, all trajectories cross the barrier before the constrained point. Therefore,
       ! the non-crossing rate is always zero.
       rateNonCrossing=0.0d0
    end if
    return
  end function linearBarrierBrownianBridgeRateNonCrossing
