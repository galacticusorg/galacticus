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
Implements a excursion set first crossing statistics class for linear barriers.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Excursion_Sets_Barriers   , only : excursionSetBarrierClass

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingLinearBarrier">
   <description>
    An excursion set first crossing statistics class for linear barriers. Specifically, the first crossing distribution is
    \begin{equation}
     f(S,t) = B(0,t) \exp(- B(S,t)^2/2S)/S/\sqrt{2 pi S},
    \end{equation}
    where $B(S,t)$ is the (assumed-to-be-linear-in-$S$) barrier at time $t$ and variance $S$. The first crossing rate is
    computed using a finite difference approximation between two closely-spaced times. The non-crossing rate is zero.
   </description>
  </excursionSetFirstCrossing>
  !!]
  type, extends(excursionSetFirstCrossingClass) :: excursionSetFirstCrossingLinearBarrier
     !!{
     A linearBarrier excursion set barrier class.
     !!}
     private
     class           (excursionSetBarrierClass     ), pointer :: excursionSetBarrier_      => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     ! The fractional step in time used to compute barrier crossing rates.
     double precision                                         :: fractionalTimeStep
   contains
     final     ::                    linearBarrierDestructor
     procedure :: probability     => linearBarrierProbability
     procedure :: rate            => linearBarrierRate
     procedure :: rateNonCrossing => linearBarrierRateNonCrossing
  end type excursionSetFirstCrossingLinearBarrier

  interface excursionSetFirstCrossingLinearBarrier
     !!{
     Constructors for the linearBarrier excursion set barrier class.
     !!}
     module procedure linearBarrierConstructorParameters
     module procedure linearBarrierConstructorInternal
  end interface excursionSetFirstCrossingLinearBarrier

contains

  function linearBarrierConstructorParameters(parameters) result(self)
    !!{
    Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (excursionSetFirstCrossingLinearBarrier)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (excursionSetBarrierClass              ), pointer       :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass         ), pointer       :: cosmologicalMassVariance_
    double precision                                                        :: fractionalTimeStep

    !![
    <inputParameter>
      <name>fractionalTimeStep</name>
      <defaultValue>0.01d0</defaultValue>
      <source>parameters</source>
      <description>The fractional time step used when computing barrier crossing rates (i.e. the step used in finite difference calculations).</description>
    </inputParameter>
    <objectBuilder class="excursionSetBarrier"      name="excursionSetBarrier_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=excursionSetFirstCrossingLinearBarrier(fractionalTimeStep,excursionSetBarrier_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="excursionSetBarrier_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function linearBarrierConstructorParameters

  function linearBarrierConstructorInternal(fractionalTimeStep,excursionSetBarrier_,cosmologicalMassVariance_) result(self)
    !!{
    Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    !!}
    implicit none
    type            (excursionSetFirstCrossingLinearBarrier)                        :: self
    double precision                                        , intent(in   )         :: fractionalTimeStep
    class           (excursionSetBarrierClass              ), intent(in   ), target :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass         ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="fractionalTimeStep, *excursionSetBarrier_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function linearBarrierConstructorInternal

  subroutine linearBarrierDestructor(self)
    !!{
    Destructor for the critical overdensity excursion set barrier class.
    !!}
    implicit none
    type(excursionSetFirstCrossingLinearBarrier), intent(inout) :: self

    !![
    <objectDestructor name="self%excursionSetBarrier_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine linearBarrierDestructor

  double precision function linearBarrierProbability(self,variance,time,node)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (excursionSetFirstCrossingLinearBarrier), intent(inout) :: self
    double precision                                        , intent(in   ) :: variance, time
    type            (treeNode                              ), intent(inout) :: node

    linearBarrierProbability=+     self%excursionSetBarrier_%barrier(   0.0d0,time,node,rateCompute=.false.)    &
         &                   *exp(                                                                              &
         &                        -0.5d0                                                                        &
         &                        *self%excursionSetBarrier_%barrier(variance,time,node,rateCompute=.false.)**2 &
         &                        /variance                                                                     &
         &                       )                                                                              &
         &                   /variance                                                                          &
         &                   /sqrt(                                                                             &
         &                         +2.0d0                                                                       &
         &                         *Pi                                                                          &
         &                         *variance                                                                    &
         &                        )
    return
  end function linearBarrierProbability

  double precision function linearBarrierRate(self,variance,varianceProgenitor,time,node)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (excursionSetFirstCrossingLinearBarrier), intent(inout) :: self
    double precision                                        , intent(in   ) :: variance                   , varianceProgenitor, &
         &                                                                     time
    type            (treeNode                              ), intent(inout) :: node
    double precision                                                        :: timeProgenitor             , massProgenitor    , &
         &                                                                     growthFactorEffective

    if (variance >= varianceProgenitor) then
       linearBarrierRate=0.0d0
    else
       ! * To estimate the rate we use a finite difference method - we compute the effective barrier for a slightly earlier time,
       !   compute the fraction of trajectories which will have crossed that effective barrier, and divide by the time difference.
       !
       ! * In Galacticus, the time evolution due to linear growth is included in the root-variance of the density field, *not* in
       !   the barrier height as is often done in the literature. As such, when computing the barrier at some earlier time we must
       !   account for the fact that, at a fixed mass, the root variance will be smaller at that earlier time. Since the solution
       !   to the excursion set problem must always be a function of δc(M,t)/√S(M,t) then we can simply scale δc by the ratio of
       !   root-variances for the progenitor at the current and earlier times.
       timeProgenitor       =+time*(1.0d0-self%fractionalTimeStep)
       massProgenitor       =+self%cosmologicalMassVariance_%mass        (sqrt(varianceProgenitor),time          )
       growthFactorEffective=+self%cosmologicalMassVariance_%rootVariance(         massProgenitor ,time          ) &
            &                /self%cosmologicalMassVariance_%rootVariance(         massProgenitor ,timeProgenitor)
       linearBarrierRate    =+     barrierEffective(variance,time,variance          ,timeProgenitor)    &
            &                *exp(                                                                      &
            &                     -0.5d0                                                                &
            &                     *barrierEffective(variance,time,varianceProgenitor,timeProgenitor)**2 &
            &                     / (+varianceProgenitor-variance)                                      &
            &                     )                                                                     &
            &                /      (+varianceProgenitor-variance)                                      &
            &                /sqrt(                                                                     &
            &                      +2.0d0                                                               &
            &                      *Pi                                                                  &
            &                      *(+varianceProgenitor-variance)                                      &
            &                     )                                                                     &
            &                /time                                                                      &
            &                /self%fractionalTimeStep
       linearBarrierRate    =max(                   &
            &                    linearBarrierRate, &
            &                    0.0d0              &
            &                   )
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

  end function linearBarrierRate

  double precision function linearBarrierRateNonCrossing(self,variance,massMinimum,time,node)
    !!{
    Return the rate for excursion set non-crossing assuming a linear barrier.
    !!}
    use :: Error_Functions, only : Error_Function
    implicit none
    class           (excursionSetFirstCrossingLinearBarrier), intent(inout) :: self
    double precision                                        , intent(in   ) :: time                           , variance                    , &
         &                                                                     massMinimum
    type            (treeNode                              ), intent(inout) :: node
    double precision                                                        :: varianceMaximum                , timeProgenitor              , &
         &                                                                     growthFactorEffective          , barrierEffectiveZeroVariance, &
         &                                                                     barrierEffectiveMaximumVariance, barrierEffectiveGradient    , &
         &                                                                     varianceDifference

    varianceMaximum=self%cosmologicalMassVariance_%rootVariance(massMinimum,time)**2
    if (variance < varianceMaximum) then
       timeProgenitor                 =+time*(1.0d0-self%fractionalTimeStep)
       growthFactorEffective          =+self%cosmologicalMassVariance_%rootVariance(         massMinimum ,time          ) &
            &                          /self%cosmologicalMassVariance_%rootVariance(         massMinimum ,timeProgenitor)
       barrierEffectiveZeroVariance   =+self%excursionSetBarrier_%barrier        (variance       ,timeProgenitor,node,rateCompute=.false.)*growthFactorEffective &
            &                          -self%excursionSetBarrier_%barrier        (variance       ,time          ,node,rateCompute=.false.)
       barrierEffectiveMaximumVariance=+self%excursionSetBarrier_%barrier        (varianceMaximum,timeProgenitor,node,rateCompute=.false.)*growthFactorEffective &
            &                          -self%excursionSetBarrier_%barrier        (variance       ,time          ,node,rateCompute=.false.)
       barrierEffectiveGradient       =+self%excursionSetBarrier_%barrierGradient(varianceMaximum,timeProgenitor,node,rateCompute=.false.)*growthFactorEffective
       varianceDifference             =+varianceMaximum-variance
       linearBarrierRateNonCrossing   =+0.5d0                                                                                                                 &
            &                          *(                                                                                                                     &
            &                            +1.0d0                                                                                                               &
            &                            -exp(-2.0d0*barrierEffectiveZeroVariance*barrierEffectiveGradient)                                                   &
            &                            +Error_Function(                                    barrierEffectiveMaximumVariance /sqrt(2.0d0*varianceDifference)) &
            &                            +exp(-2.0d0*barrierEffectiveZeroVariance*barrierEffectiveGradient)                                                   &
            &                            *Error_Function((2.0d0*barrierEffectiveZeroVariance-barrierEffectiveMaximumVariance)/sqrt(2.0d0*varianceDifference)) &
            &                           )                                                                                                                     &
            &                          /time                                                                                                                  &
            &                          /self%fractionalTimeStep
    else
       linearBarrierRateNonCrossing   =+1.0d0                   &
            &                          /time                    &
            &                          /self%fractionalTimeStep
    end if
    return
  end function linearBarrierRateNonCrossing
