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
  Implements a excursion set first crossing statistics class utilizing the algorithm of \cite{zhang_random_2006}.
  !!}

  use :: Excursion_Sets_Barriers, only : excursionSetBarrierClass
  use :: Numerical_Interpolation, only : interpolator

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingZhangHui">
   <description>
    An excursion set first crossing statistics class utilizing the algorithm of \cite{zhang_random_2006}. First crossing (and
    non-crossing) rates are not supported by this method.
   </description>
  </excursionSetFirstCrossing>
  !!]
  type, extends(excursionSetFirstCrossingClass) :: excursionSetFirstCrossingZhangHui
     !!{
     An excursion set first crossing statistics class utilizing the algorithm of \cite{zhang_random_2006}.
     !!}
     private
     class           (excursionSetBarrierClass), pointer                     :: excursionSetBarrier_         => null()
     double precision                                                        :: timeMaximum                  =  0.0d0  , timeMinimum             =0.0d0 , &
          &                                                                     varianceMaximum              =  0.0d0
     integer                                                                 :: timeTableCount                         , varianceTableCount
     double precision                          , allocatable, dimension(:,:) :: firstCrossingProbabilityTable
     double precision                          , allocatable, dimension(:  ) :: timeTable                              , varianceTable
     double precision                                                        :: varianceTableStep
     logical                                                                 :: tableInitialized             =  .false.
     type            (interpolator            ), allocatable                 :: interpolatorTime                       , interpolatorVariance
     ! Stored values.
     double precision                                                        :: variancePrevious                       , timePrevious                   , &
          &                                                                     barrierStored                          , barrierGradientStored
     ! Stored values for Delta function.
     integer                                                                 :: iDeltaPrevious                         , jDeltaPrevious
     double precision                                                        :: timeDeltaPrevious                      , deltaStored
     ! Variables used in integrations.
     double precision                                                        :: barrierIntegrand                       , barrierGradientIntegrand       , &
          &                                                                     timeIntegrand                          , varianceIntegrand
   contains
     !![
     <methods>
       <method description="Returns the function $g_1(S)$ \citep{zhang_random_2006}." method="g1" />
       <method description="Returns the function $g_2(S,S^\prime)$ \citep{zhang_random_2006}." method="g2" />
       <method description="Returns the function $g_2(S,S^\prime)$ integrated over a range $\Delta S$ \citep{zhang_random_2006}." method="g2Integrated" />
       <method description="Returns the function $g_2(S,S^\prime)$ integrated over a range $\Delta S$ \citep{zhang_random_2006}." method="delta" />
     </methods>
     !!]
     final     ::                    zhangHuiDestructor
     procedure :: probability     => zhangHuiProbability
     procedure :: rate            => zhangHuiRate
     procedure :: rateNonCrossing => zhangHuiRateNonCrossing
     procedure :: g1              => zhangHuiG1
     procedure :: g2              => zhangHuiG2
     procedure :: g2Integrated    => zhangHuiG2Integrated
     procedure :: delta           => zhangHuiDelta
  end type excursionSetFirstCrossingZhangHui

  interface excursionSetFirstCrossingZhangHui
     !!{
     Constructors for the \cite{zhang_random_2006} excursion set barrier class.
     !!}
     module procedure zhangHuiConstructorParameters
     module procedure zhangHuiConstructorInternal
  end interface excursionSetFirstCrossingZhangHui

  ! Table granularity.
  integer, parameter :: timeTableNumberPerDecade=40, varianceTableNumberPerUnit=400

contains

  function zhangHuiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (excursionSetFirstCrossingZhangHui)                :: self
    type (inputParameters                  ), intent(inout) :: parameters
    class(excursionSetBarrierClass         ), pointer       :: excursionSetBarrier_

    !![
    <objectBuilder class="excursionSetBarrier" name="excursionSetBarrier_" source="parameters"/>
    !!]
    self=excursionSetFirstCrossingZhangHui(excursionSetBarrier_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="excursionSetBarrier_"/>
    !!]
    return
  end function zhangHuiConstructorParameters

  function zhangHuiConstructorInternal(excursionSetBarrier_) result(self)
    !!{
    Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    !!}
    implicit none
    type (excursionSetFirstCrossingZhangHui)                        :: self
    class(excursionSetBarrierClass         ), intent(in   ), target :: excursionSetBarrier_
    !![
    <constructorAssign variables="*excursionSetBarrier_"/>
    !!]

    self% tableInitialized=.false.
    self%  varianceMaximum=-huge(0.0d0)
    self%      timeMinimum=+huge(0.0d0)
    self%      timeMaximum=-huge(0.0d0)
    self% variancePrevious=-huge(0.0d0)
    self%     timePrevious=-huge(0.0d0)
    self%timeDeltaPrevious=-huge(0.0d0)
    self%   iDeltaPrevious=-     1
    self%   jDeltaPrevious=-     1
    return
  end function zhangHuiConstructorInternal

  subroutine zhangHuiDestructor(self)
    !!{
    Destructor for the critical overdensity excursion set barrier class.
    !!}
    implicit none
    type(excursionSetFirstCrossingZhangHui), intent(inout) :: self

    !![
    <objectDestructor name="self%excursionSetBarrier_" />
    !!]
    return
  end subroutine zhangHuiDestructor

  double precision function zhangHuiProbability(self,variance,time,node)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    use            :: Display          , only : displayCounter       , displayCounterClear, displayIndent       , displayUnindent, &
          &                                     verbosityLevelWorking
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    use            :: Numerical_Ranges , only : Make_Range           , rangeTypeLinear    , rangeTypeLogarithmic
    implicit none
    class           (excursionSetFirstCrossingZhangHui), intent(inout)  :: self
    double precision                                   , intent(in   )  :: variance         , time
    type            (treeNode                         ), intent(inout)  :: node
    double precision                                   , dimension(0:1) :: hTime            , hVariance
    logical                                                             :: makeTable
    integer         (c_size_t                         )                 :: iTime            , iVariance
    integer                                                             :: i                , j        , &
         &                                                                 jTime            , jVariance
    double precision                                                    :: summedProbability

    ! Determine if we need to make the table.
    makeTable=           .not. self%tableInitialized  &
         &    .or.                                    &
         &     (variance >     self%varianceMaximum ) &
         &    .or.                                    &
         &     (time     <     self%timeMinimum     ) &
         &    .or.                                    &
         &     (time     >     self%timeMaximum     )
    if (makeTable) then
       ! Construct the table of variance on which we will solve for the first crossing distribution.
       if (allocated(self%varianceTable                )) deallocate(self%varianceTable                )
       if (allocated(self%timeTable                    )) deallocate(self%timeTable                    )
       if (allocated(self%firstCrossingProbabilityTable)) deallocate(self%firstCrossingProbabilityTable)
       self%varianceMaximum   =max(self%varianceMaximum,variance)
       self%varianceTableCount=int(self%varianceMaximum*dble(varianceTableNumberPerUnit))
       if (self%tableInitialized) then
          self%timeMinimum=min(self%timeMinimum,time)
          self%timeMaximum=max(self%timeMaximum,time)
       else
          self%timeMinimum=0.5d0*time
          self%timeMaximum=2.0d0*time
       end if
       self%timeTableCount=int(log10(self%timeMaximum/self%timeMinimum)*dble(timeTableNumberPerDecade))+1
       allocate(self%varianceTable                (0:self%varianceTableCount                    ))
       allocate(self%timeTable                    (                          self%timeTableCount))
       allocate(self%firstCrossingProbabilityTable(0:self%varianceTableCount,self%timeTableCount))
       self%timeTable        =Make_Range(self%timeMinimum,self%timeMaximum    ,self%timeTableCount      ,rangeType=rangeTypeLogarithmic)
       self%varianceTable    =Make_Range(0.0d0           ,self%varianceMaximum,self%varianceTableCount+1,rangeType=rangeTypeLinear     )
       self%varianceTableStep=+self%varianceTable(1) &
            &                 -self%varianceTable(0)
       ! Loop through the table and solve for the first crossing distribution.
       call displayIndent("solving for excursion set barrier crossing probabilities",verbosityLevelWorking)
       do iTime=1,self%timeTableCount
          do i=0,self%varianceTableCount
             call displayCounter(int(100.0d0*dble(i+(iTime-1)*self%varianceTableCount)/dble(self%varianceTableCount*self%timeTableCount)),i==0 .and. iTime==1,verbosityLevelWorking)
             if      (i  > 2) then
                summedProbability=0.0d0
                do j=1,i-1
                   summedProbability                       =+summedProbability                                                                                                                               &
                        &                                   +    self%firstCrossingProbabilityTable(    j                                                                       ,               iTime)       &
                        &                                   *(                                                                                                                                               &
                        &                                     +  self%delta                        (i  ,j  ,self%varianceTable(i),self%varianceTable(j  ),self%varianceTableStep,self%timeTable(iTime),node)&
                        &                                     +  self%delta                        (i  ,j+1,self%varianceTable(i),self%varianceTable(j+1),self%varianceTableStep,self%timeTable(iTime),node) &
                        &                                    )
                end do
                self%firstCrossingProbabilityTable(i,iTime)=+(                                                                                                                                               &
                     &                                        +  self%g1                           (        self%varianceTable(i)                                               ,self%timeTable(iTime),node) &
                     &                                        +summedProbability                                                                                                                             &
                     &                                       )                                                                                                                                               &
                     &                                      /(                                                                                                                                               &
                     &                                        +1.0d0                                                                                                                                         &
                     &                                        -  self%delta                        (i  ,i  ,self%varianceTable(i),self%varianceTable(i  ),self%varianceTableStep,self%timeTable(iTime),node) &
                     &                                       )
             else if (i == 2) then
                self%firstCrossingProbabilityTable(i,iTime)=+(                                                                                                                                               &
                     &                                        +  self%g1                           (        self%varianceTable(i)                                               ,self%timeTable(iTime),node) &
                     &                                        +  self%firstCrossingProbabilityTable(i-1                                                                         ,               iTime)       &
                     &                                        *(                                                                                                                                             &
                     &                                          +self%delta                        (i  ,1  ,self%varianceTable(i),self%varianceTable(1  ),self%varianceTableStep,self%timeTable(iTime),node) &
                     &                                          +self%delta                        (i  ,2  ,self%varianceTable(i),self%varianceTable(2  ),self%varianceTableStep,self%timeTable(iTime),node) &
                     &                                         )                                                                                                                                             &
                     &                                       )                                                                                                                                               &
                     &                                      /(                                                                                                                                               &
                     &                                        +1.0d0                                                                                                                                         &
                     &                                        -  self%delta                        (i  ,i  ,self%varianceTable(i),self%varianceTable(i  ),self%varianceTableStep,self%timeTable(iTime),node) &
                     &                                       )
             else if (i == 1) then
                self%firstCrossingProbabilityTable(i,iTime)=+    self%g1                           (        self%varianceTable(i)                                               ,self%timeTable(iTime),node) &
                     &                                      /(                                                                                                                                               &
                     &                                        +1.0d0                                                                                                                                         &
                     &                                        -  self%delta                        (1  ,  1,self%varianceTable(1),self%varianceTable(1  ),self%varianceTableStep,self%timeTable(iTime),node) &
                     &                                       )
             else if (i == 0) then
                self%firstCrossingProbabilityTable(i,iTime)=+0.0d0
             end if
          end do
       end do
       call displayCounterClear(verbosityLevelWorking)
       call displayUnindent("done",verbosityLevelWorking)
       ! Build the interpolators.
       if (allocated(self%interpolatorVariance)) deallocate(self%interpolatorVariance)
       if (allocated(self%interpolatorTime    )) deallocate(self%interpolatorTime    )
       allocate(self%interpolatorVariance)
       allocate(self%interpolatorTime    )
       self%interpolatorVariance=interpolator(self%varianceTable)
       self%interpolatorTime    =interpolator(self%timeTable    )
       ! Record that the table is now built.
       self%tableInitialized          =.true.
    end if
    ! Get interpolating factors.
    call self%interpolatorTime    %linearFactors(time    ,iTime    ,hTime    )
    call self%interpolatorVariance%linearFactors(variance,iVariance,hVariance)
    ! Compute first crossing probability by interpolating.
    zhangHuiProbability=0.0d0
    do jTime=0,1
       do jVariance=0,1
          zhangHuiProbability=zhangHuiProbability+hTime(jTime)*hVariance(jVariance)*self%firstCrossingProbabilityTable(iVariance-1+jVariance,iTime+jTime)
       end do
    end do
    return
  end function zhangHuiProbability

  double precision function zhangHuiRate(self,variance,varianceProgenitor,time,node)
    !!{
    Return the excursion set barrier at the given variance and time. This method is not implemented.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (excursionSetFirstCrossingZhangHui), intent(inout) :: self
    double precision                                   , intent(in   ) :: variance, varianceProgenitor, &
         &                                                                time
    type            (treeNode                         ), intent(inout) :: node
    !$GLC attributes unused :: self, time, variance, varianceProgenitor, node

    zhangHuiRate=0.0d0
    call Error_Report('barrier crossing rates are not implemented for this method [too slow]'//{introspection:location})
    return
  end function zhangHuiRate

  double precision function zhangHuiRateNonCrossing(self,variance,massMinimum,time,node)
    !!{
    Return the rate for excursion set non-crossing. This method is not implemented.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (excursionSetFirstCrossingZhangHui), intent(inout) :: self
    double precision                                   , intent(in   ) :: time           , variance, &
         &                                                                massMinimum
    type            (treeNode                         ), intent(inout) :: node
    !$GLC attributes unused :: self, time, variance, massMinimum, node

    zhangHuiRateNonCrossing=0.0d0
    call Error_Report('barrier non-crossing rates are not implemented for this method [too slow]'//{introspection:location})
    return
  end function zhangHuiRateNonCrossing

  double precision function zhangHuiG1(self,variance,time,node)
    !!{
    Returns the function $g_1(S)$ in the \cite{zhang_random_2006} algorithm for excursion set barrier crossing probabilities.
    !!}
    use :: Math_Distributions_Gaussian, only : Gaussian_Distribution
    implicit none
    class           (excursionSetFirstCrossingZhangHui), intent(inout) :: self
    double precision                                   , intent(in   ) :: time                , variance
    double precision                                                   :: barrier
    type            (treeNode                         ), intent(inout) :: node

    barrier   =+  self%excursionSetBarrier_%barrier        (variance,time,node,rateCompute=.false.)
    zhangHuiG1=+(                                                                                   &
         &       +barrier                                                                           &
         &       /variance                                                                          &
         &       -2.0d0                                                                             &
         &       *self%excursionSetBarrier_%barrierGradient(variance,time,node,rateCompute=.false.) &
         &      )                                                                                   &
         &     *Gaussian_Distribution(barrier,sqrt(variance))
    return
  end function zhangHuiG1

  double precision function zhangHuiG2(self,variance,variancePrimed,time,node)
    !!{
    Returns the function $g_2(S,S^\prime)$ in the \cite{zhang_random_2006} algorithm for excursion set barrier crossing probabilities.
    !!}
    use :: Math_Distributions_Gaussian, only : Gaussian_Distribution
    implicit none
    class           (excursionSetFirstCrossingZhangHui), intent(inout) :: self
    double precision                                   , intent(in   ) :: time          , variance, &
         &                                                                variancePrimed
    type            (treeNode                         ), intent(inout) :: node
    double precision                                                   :: barrierPrimed

    ! Compute the barriers.
    if (variance /= self%variancePrevious .or. time /= self%timePrevious) then
       self%variancePrevious     =variance
       self%    timePrevious     =time
       self%barrierStored        =self%excursionSetBarrier_%barrier        (variance,time,node,rateCompute=.false.)
       self%barrierGradientStored=self%excursionSetBarrier_%barrierGradient(variance,time,node,rateCompute=.false.)
    end if
    barrierPrimed=self%excursionSetBarrier_%barrier(variancePrimed,time,node,rateCompute=.false.)
    ! Compute the function.
    zhangHuiG2=+(                                                                                               &
         &       +2.0d0                                                                                         &
         &       *                     self%barrierGradientStored                                               &
         &       -                    (self%barrierStored        -barrierPrimed)/    (variance-variancePrimed)  &
         &      )                                                                                               &
         &     *Gaussian_Distribution( self%barrierStored        -barrierPrimed ,sqrt(variance-variancePrimed))
    return
  end function zhangHuiG2

  double precision function zhangHuiDelta(self,i,j,iVariance,jVariance,deltaVariance,time,node)
    !!{
    Returns the factor $\Delta{i,j}$ in the \cite{zhang_random_2006} algorithm for excursion set barrier crossing probabilities.
    !!}
    implicit none
    class           (excursionSetFirstCrossingZhangHui), intent(inout) :: self
    integer                                            , intent(in   ) :: i               , j
    double precision                                   , intent(in   ) :: deltaVariance   , iVariance          , jVariance, &
         &                                                                time
    type            (treeNode                         ), intent(inout) :: node

    if (.not.(i == self%iDeltaPrevious .and. j == self%jDeltaPrevious .and. time == self%timeDeltaPrevious)) then
       self%   iDeltaPrevious=i
       self%   jDeltaPrevious=j
       self%timeDeltaPrevious=time
       if (i == j) then
          ! In this case integrate over the range to get an average value.
          self%deltaStored=0.5d0              *self%g2Integrated(iVariance,         +deltaVariance      ,time,node)
       else
          ! Compute the appropriate factor.
          self%deltaStored=0.5d0*deltaVariance*self%g2          (iVariance,jVariance-deltaVariance/2.0d0,time,node)
       end if
    end if
    zhangHuiDelta=self%deltaStored
    return
  end function zhangHuiDelta

  double precision function zhangHuiG2Integrated(self,variance,deltaVariance,time,node)
    !!{
    Integrated function $g_2(S,S^\prime)$ in the \cite{zhang_random_2006} algorithm for excursion set barrier crossing probabilities.
    !!}
    use :: Numerical_Comparison , only : Values_Differ
    use :: Numerical_Integration, only : GSL_Integ_Gauss15, integrator
    implicit none
    class           (excursionSetFirstCrossingZhangHui), intent(inout) :: self
    double precision                                   , intent(in   ) :: deltaVariance                 , time           , &
         &                                                                variance
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , parameter     :: gradientChangeTolerance=1.0d-3
    double precision                                                   :: smallStep                     , barrierGradient, &
         &                                                                barrier
    type            (integrator                       )                :: integrator_

    ! Store variables needed in the integrand.
    barrier        =self%excursionSetBarrier_%barrier        (variance,time,node,rateCompute=.false.)
    barrierGradient=self%excursionSetBarrier_%barrierGradient(variance,time,node,rateCompute=.false.)
    ! Find a suitably small step in variance that allows us to compute the divergent part of the integral with an analytic
    ! approximation. The approximation used assumes that the barrier gradient, dB/dS, is constant, so find a step over which
    ! the gradient is constant to within a specified tolerance.
    smallStep=deltaVariance
    do while (Values_Differ(self%excursionSetBarrier_%barrierGradient(variance-smallStep,time,node,rateCompute=.false.),barrierGradient,relTol=gradientChangeTolerance))
       smallStep=0.5d0*smallStep
    end do
    ! Compute the non-divergent part of the integral numerically.
    integrator_         =integrator           (                                          &
         &                                                        zhangHuiG2Integrand  , &
         &                                     toleranceAbsolute=1.0d-50               , &
         &                                     toleranceRelative=1.0d-06               , &
         &                                     hasSingularities =.true.                , &
         &                                     integrationRule  =GSL_Integ_Gauss15       &
         &                                    )
    zhangHuiG2Integrated=integrator_%integrate(                                          &
         &                                                       variance-deltaVariance, &
         &                                                       variance-smallStep      &
         &                                    )
    ! Compute the divergent part of the integral with an analytic approximation.
    zhangHuiG2Integrated=+zhangHuiG2Integrated               &
         &               +erf(                               &
         &                    +barrierGradient &
         &                    *sqrt(                         &
         &                          +0.5d0                   &
         &                          *smallStep               &
         &                         )                         &
         &                   )
    return

  contains

    double precision function zhangHuiG2Integrand(variancePrimed)
      !!{
      Integrand function used in computing $\Delta_{i,i}$ in the \cite{zhang_random_2006} algorithm for excursion set barrier
      crossing probabilities.
      !!}
      use :: Math_Distributions_Gaussian, only : Gaussian_Distribution
      implicit none
      double precision, intent(in   ) :: variancePrimed
      double precision                :: barrierPrimed

      if (variancePrimed >= variance) then
         zhangHuiG2Integrand=0.0d0
      else
         barrierPrimed      =+self%excursionSetBarrier_%barrier(variancePrimed,time,node,rateCompute=.false.)
         zhangHuiG2Integrand=+(                                                                             &
              &                +2.0d0                                                                       &
              &                *barrierGradient                                                             &
              &                -                    (+barrier-barrierPrimed)/    (variance-variancePrimed)  &
              &               )                                                                             &
              &              *Gaussian_Distribution( +barrier-barrierPrimed ,sqrt(variance-variancePrimed))
      end if
      return
    end function zhangHuiG2Integrand

  end function zhangHuiG2Integrated
