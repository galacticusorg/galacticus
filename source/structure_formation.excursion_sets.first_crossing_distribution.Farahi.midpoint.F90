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

!+    Contributions to this file made by: Andrew Benson, Christoph Behrens, Xiaolong Du.

!!{
Implements an excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint
method to perform the integrations \citep{du_substructure_2017}.
!!}

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingFarahiMidpoint">
    <description>
      An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method
      to perform the integrations \citep{du_substructure_2017}.

      Specifically, in the method used by \cite{benson_dark_2012} the integral equation for $f(S)$ (see
      equation~\ref{eq:OldExcursionMethod}) is solved using the trapezoidal rule (see \refClass{excursionSetFirstCrossingFarahi}
      for complete details):
      \begin{equation}
       \int_0^{S_j} f(\tilde{S})K(S_j,\tilde{S}) \mathrm{d}\tilde{S} = \sum_{j=0}^{i-1} \frac{f(S_j)K(S_i,S_j)+f(S_{j+1})K(S_i,S_{j+1})}{2} \Delta S_j,
     \end{equation}
     where
     \begin{equation}
       K(S_i,\tilde{S}) = \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(\tilde{S}), S_i, \tilde{S}]}{\sqrt{2 \Delta S[S_i,\tilde{S}]}}\right\},
     \end{equation}
     is the kernel of the integral equation for $f(S)$ which is being solved.

     The method used in this class increases the precision of the solver by replacing the trapezoidal integration rule in the
     above with a mid-point integration rule:
     \begin{equation}
       \int_0^{S_j} f(\tilde{S})K(S_j,\tilde{S}) \mathrm{d}\tilde{S} = \sum_{j=0}^{i-1} f(S_{j+1/2})K(S_i,S_{j+1/2}) \Delta S_j,
      \end{equation}
      such that the first crossing distribution at $S_{i-1/2}$ is given by
      \begin{equation}
       f(S_{i-1/2}) = \frac{1}{K(S_i,S_{i-1/2})\Delta S_{i-1/2}} \left( \hbox{erfc}\left\{ \frac{\Delta \delta [B(S_i),B(S_1),S_i,S_1]}{\sqrt{2 \Delta S[S_i,S_1]}} \right\} - \sum_{j=0}^{i-2} f(S_{j+1/2} K(S_i,S_{j+1/2}) \Delta S_j) \right).
      \end{equation}
    </description>
  </excursionSetFirstCrossing>
  !!]
  type, extends(excursionSetFirstCrossingFarahi) :: excursionSetFirstCrossingFarahiMidpoint
     !!{
     An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method
     to perform the integrations \citep{du_substructure_2017}.
     !!}
     private
   contains
     procedure :: probability  => farahiMidpointProbability
     procedure :: rateTabulate => farahiMidpointRateTabulate
  end type excursionSetFirstCrossingFarahiMidpoint

  interface excursionSetFirstCrossingFarahiMidpoint
     !!{
     Constructors for the Farahi-midpoint excursion set barrier class.
     !!}
     module procedure farahiMidpointConstructorParameters
     module procedure farahiMidpointConstructorInternal
  end interface excursionSetFirstCrossingFarahiMidpoint

contains

  function farahiMidpointConstructorParameters(parameters,inputParametersValidate) result(self)
    !!{
    Constructor for the Farahi-midpoint excursion set class first crossing class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (excursionSetFirstCrossingFarahiMidpoint)                          :: self
    type   (inputParameters                        ), intent(inout)           :: parameters
    logical                                         , intent(in   ), optional :: inputParametersValidate
    !![
    <optionalArgument name="inputParametersValidate" defaultsTo=".true."/>
    !!]
    
    self%excursionSetFirstCrossingFarahi=excursionSetFirstCrossingFarahi(parameters,inputParametersValidate)
    if (inputParametersValidate_) then
       !![
       <inputParametersValidate source="parameters"/>
       !!]
    end if
    return
  end function farahiMidpointConstructorParameters

  function farahiMidpointConstructorInternal(fractionalTimeStep,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,varianceNumberPerDecadeNonCrossing,timeNumberPerDecade,varianceIsUnlimited,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the Farahi-midpoint excursion set class first crossing class.
    !!}
    implicit none
    type            (excursionSetFirstCrossingFarahiMidpoint)                        :: self
    double precision                                         , intent(in   )         :: fractionalTimeStep
    integer                                                  , intent(in   )         :: varianceNumberPerUnitProbability  , varianceNumberPerUnit  , &
         &                                                                              timeNumberPerDecade               , varianceNumberPerDecade, &
         &                                                                              varianceNumberPerDecadeNonCrossing
    logical                                                  , intent(in   )         :: varianceIsUnlimited
    type            (varying_string                         ), intent(in   )         :: fileName
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetBarrierClass               ), intent(in   ), target :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass          ), intent(in   ), target :: cosmologicalMassVariance_

    self%excursionSetFirstCrossingFarahi=excursionSetFirstCrossingFarahi(fractionalTimeStep,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,varianceNumberPerDecadeNonCrossing,timeNumberPerDecade,varianceIsUnlimited,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_)
    return
  end function farahiMidpointConstructorInternal

  double precision function farahiMidpointProbability(self,variance,time,node)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    use :: Display         , only : displayCounter              , displayCounterClear  , displayIndent       , displayMessage, &
          &                         displayUnindent             , verbosityLevelWorking
    use :: Error_Functions , only : Error_Function_Complementary
    use :: File_Utilities  , only : File_Lock                   , File_Unlock          , lockDescriptor
    use :: Kind_Numbers    , only : kind_dble                   , kind_quad
    use :: MPI_Utilities   , only : mpiBarrier                  , mpiSelf
    use :: Numerical_Ranges, only : Make_Range                  , rangeTypeLinear      , rangeTypeLogarithmic
    use :: Table_Labels    , only : extrapolationTypeFix
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpoint), intent(inout)                 :: self
    double precision                                         , intent(in   )                 :: variance                       , time
    type            (treeNode                               ), intent(inout)                 :: node
    double precision                                                        , dimension(0:1) :: hTime                          , hVariance
    double precision                                         , parameter                     :: varianceTolerance       =1.0d-6
    double precision                                         , allocatable  , dimension( : ) :: varianceMidpoint               , barrier         , &
         &                                                                                      barrierMidpoint
    double precision                                                                         :: barrierTest
    class           (excursionSetBarrierClass               ), pointer                       :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass          ), pointer                       :: cosmologicalMassVariance_
    logical                                                                                  :: makeTable
    integer         (c_size_t                               )                                :: iTime                          , iVariance       , &
         &                                                                                      loopCount                      , loopCountTotal  , &
         &                                                                                      i                              , j               , &
         &                                                                                      jTime                          , jVariance
    double precision                                                                         :: probabilityCrossingPrior
    double precision                                                                         :: integralKernel
    real            (kind_quad                              )                                :: offsetEffective                , varianceResidual
    character       (len =9                                 )                                :: label
    type            (varying_string                         )                                :: message
    type            (lockDescriptor                         )                                :: fileLock
    
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
    !  Sᵢ                 = self%variance                (i      )
    !  B(Sᵢ)              =      barrier                 (i      )
    !  f(Sᵢ,t)            = self%firstCrossingProbability(i,iTime)
    !  Δδ[t,S₁,S₂,δ₁,δ₂] = self%offsetEffective         (self%time(iTime),0,S1,S2,0,barrier1,barrier2)
    !  ΔS[t,S₁,S₂]       = self%varianceResidual        (self%time(iTime),0,S1,S2                    )
    
    ! The suffix "Midpoint" refers to a quantity computed at the midpoint between two tabulated points.

    ! Read tables from file if possible.
    do i=1,2
       makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+varianceTolerance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
       if (i == 1 .and. self%useFile .and. makeTable) then
          call self%fileNameInitialize()
          call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       else
          exit
       end if
    end do
#ifdef USEMPI
    if (self%coordinatedMPI_) call mpiBarrier()
#endif
    if (makeTable) then
       !$omp critical(farahiMidpointProbabilityTabulate)
       ! Attempt to read the file again now that we are within the critical section. If another thread made the file while we were waiting we may be able to skip building the table.
       if (self%useFile) then
          call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       end if
       makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+varianceTolerance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
       if (makeTable) then
          ! Construct the table of variance on which we will solve for the first crossing distribution.
          if (allocated(self%variance                )) deallocate(self%variance                )
          if (allocated(self%time                    )) deallocate(self%time                    )
          if (allocated(self%firstCrossingProbability)) deallocate(self%firstCrossingProbability)
          self%varianceMaximum=max(self%varianceMaximum,variance)
          self%countVariance  =int(self%varianceMaximum*dble(self%varianceNumberPerUnitProbability))
          if (self%tableInitialized) then
             self%timeMinimum=min(self%timeMinimum,time/10.0d0**(2.0d0/dble(self%timeNumberPerDecade)))
             self%timeMaximum=max(self%timeMaximum,time*10.0d0**(2.0d0/dble(self%timeNumberPerDecade)))
          else
             self%timeMinimum=                     time/10.0d0**(2.0d0/dble(self%timeNumberPerDecade))
             self%timeMaximum=                     time*10.0d0**(2.0d0/dble(self%timeNumberPerDecade))
          end if
          self%countTime=max(2,int(log10(self%timeMaximum/self%timeMinimum)*dble(self%timeNumberPerDecade))+1)
          allocate(self%variance                (0:self%countVariance               ))
          allocate(self%time                    (                     self%countTime))
          allocate(self%firstCrossingProbability(0:self%countVariance,self%countTime))
          allocate(     varianceMidpoint        (0:self%countVariance               ))
          self%time        =Make_Range(self%timeMinimum,self%timeMaximum    ,self%countTime      ,rangeType=rangeTypeLogarithmic)
          self%variance    =Make_Range(0.0d0           ,self%varianceMaximum,self%countVariance+1,rangeType=rangeTypeLinear     )
          self%varianceStep=+self%variance(1) &
               &            -self%variance(0)
          ! Compute the variance at the mid-points.
          varianceMidpoint(0)=0.0d0
          forall(i=1:self%countVariance)
             varianceMidpoint(i)=(self%variance(i-1)+self%variance(i))/2.0d0
          end forall
          ! Loop through the table and solve for the first crossing distribution.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayIndent("solving for excursion set barrier crossing probabilities",verbosityLevelWorking)
             message="    time: "
             write (label,'(f6.3)') self%timeMinimum
             message=message//label//" to "
             write (label,'(f6.3)') self%timeMaximum
             message=message//label
             call displayMessage(message,verbosityLevelWorking)
             message="variance: "
             write (label,'(f9.3)') self%varianceMaximum
             message=message//label
             call displayMessage(message,verbosityLevelWorking)
#ifdef USEMPI
          end if
#endif
#ifdef USEMPI
          if (mpiSelf%isMaster() .and. self%coordinatedMPI_) then
             loopCountTotal=(int(self%countTime,kind=c_size_t)/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t)*(int(self%countVariance-1,kind=c_size_t)*int(self%countVariance,kind=c_size_t))/2_c_size_t
          else
#endif
             loopCountTotal= int(self%countTime,kind=c_size_t)                                               *(int(self%countVariance-1,kind=c_size_t)*int(self%countVariance,kind=c_size_t))/2_c_size_t
#ifdef USEMPI
          end if
#endif
          loopCount=0
#ifdef USEMPI
          if (self%coordinatedMPI_) self%firstCrossingProbability=0.0d0
#endif
          ! Make a call to the barrier function at maximum variance for the minimum and maximum times so that the barrier function
          ! is initialized and covers the whole range in which we are interested.
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMinimum,node,rateCompute=.false.)
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMaximum,node,rateCompute=.false.)
          ! Enter an OpenMP parallel region. Each parallel thread will solve for the first crossing distribution at a different epoch.
          !$omp parallel private(iTime,i,j,probabilityCrossingPrior,integralKernel,excursionSetBarrier_,cosmologicalMassVariance_,barrier,barrierMidpoint,offsetEffective,varianceResidual) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
          allocate(excursionSetBarrier_     ,mold=self%excursionSetBarrier_     )
          allocate(cosmologicalMassVariance_,mold=self%cosmologicalMassVariance_)
          !$omp critical(excursionSetsSolverFarahiMidpointDeepCopy)
          !![
          <deepCopyReset variables="self%excursionSetBarrier_ self%cosmologicalMassVariance_"/>
          <deepCopy source="self%excursionSetBarrier_"      destination="excursionSetBarrier_"     />
          <deepCopy source="self%cosmologicalMassVariance_" destination="cosmologicalMassVariance_"/>
          <deepCopyFinalize variables="excursionSetBarrier_ cosmologicalMassVariance_"/>
	  !!]
          !$omp end critical(excursionSetsSolverFarahiMidpointDeepCopy)
          allocate(barrier        (0:self%countVariance))
          allocate(barrierMidpoint(0:self%countVariance))
          !$omp do schedule(dynamic)
          do iTime=1,self%countTime
#ifdef USEMPI
             if (self%coordinatedMPI_ .and. mod(iTime-1,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
             ! Construct the barrier table.
             do i=0,self%countVariance
                barrier        (i)=excursionSetBarrier_%barrier(self%variance        (i),self%time(iTime),node,rateCompute=.false.)
                barrierMidpoint(i)=excursionSetBarrier_%barrier(     varianceMidpoint(i),self%time(iTime),node,rateCompute=.false.)
             end do
             ! Find the first crossing distribution at the first grid point.
             offsetEffective                       =+self%offsetEffective (self%time(iTime),0.0_kind_quad,real(self%variance(1),kind_quad),real(varianceMidpoint(1),kind_quad),0.0_kind_quad,real(barrier(1),kind_quad),real(barrierMidpoint(1),kind_quad),cosmologicalMassVariance_)
             varianceResidual                      =+self%varianceResidual(self%time(iTime),0.0_kind_quad,real(self%variance(1),kind_quad),real(varianceMidpoint(1),kind_quad)                                                                            ,cosmologicalMassVariance_)
             self%firstCrossingProbability(0,iTime)=+0.0d0
             integralKernel                        =+real(                                                                    &
                  &                                       Error_Function_Complementary(                                       &
                  &                                                                    +offsetEffective                       &
                  &                                                                    /sqrt(2.0_kind_quad*varianceResidual)  &
                  &                                                                   )                                     , &
                  &                                        kind=kind_dble                                                     &
                  &                                       )
             offsetEffective                       =+self%offsetEffective (self%time(iTime),0.0_kind_quad,real(self%variance(1),kind_quad),0.0_kind_quad                      ,0.0_kind_quad,real(barrier(1),kind_quad),0.0_kind_quad                      ,cosmologicalMassVariance_)
             varianceResidual                      =+self%varianceResidual(self%time(iTime),0.0_kind_quad,real(self%variance(1),kind_quad),0.0_kind_quad                                                                                                   ,cosmologicalMassVariance_)
             self%firstCrossingProbability(1,iTime)=+Error_Function_Complementary(                                   &
                  &                                                               +barrier(1)                        &
                  &                                                               /sqrt(2.0d0*self%variance(1)) &
                  &                                                              )                                   &
                  &                                 /self%varianceStep                                               &
                  &                                 /integralKernel
             ! Iterate over remaining variance grid points and solve for the first crossing rate.
             do i=2,self%countVariance
#ifdef USEMPI
                if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
                   call displayCounter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityLevelWorking)
#ifdef USEMPI
                end if
#endif
                !$omp atomic
                loopCount               =loopCount+(i-1)
                ! Iterate over smaller variances, accumulating their contribution to the first crossing rate at the current variance.
                probabilityCrossingPrior=0.0d0
                do j=1,i-1
                   offsetEffective         =+self%offsetEffective (self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),real(varianceMidpoint(j),kind_quad),0.0_kind_quad,real(barrier(i),kind_quad),real(barrierMidpoint(j),kind_quad),cosmologicalMassVariance_)
                   varianceResidual        =+self%varianceResidual(self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),real(varianceMidpoint(j),kind_quad)                                                                            ,cosmologicalMassVariance_)
                   probabilityCrossingPrior=+probabilityCrossingPrior                                                &
                        &                   +self%firstCrossingProbability(j,iTime)                                  &
                        &                   *real(                                                                   &
                        &                        Error_Function_Complementary(                                       &
                        &                                                     +offsetEffective                       &
                        &                                                     /sqrt(2.0_kind_quad*varianceResidual)  &
                        &                                                    )                                     , &
                        &                         kind=kind_dble                                                     &
                        &                        )
                end do
                ! Compute the resulting first crossing distribution.
                offsetEffective =+self%offsetEffective (self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),real(varianceMidpoint(i),kind_quad),0.0_kind_quad,real(barrier(i),kind_quad),real(barrierMidpoint(i),kind_quad),cosmologicalMassVariance_)
                varianceResidual=+self%varianceResidual(self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),real(varianceMidpoint(i),kind_quad)                                                                            ,cosmologicalMassVariance_)
                integralKernel  =+real(                                                                   &
                        &             Error_Function_Complementary(                                       &
                        &                                          +offsetEffective                       &
                        &                                          /sqrt(2.0_kind_quad*varianceResidual)  &
                        &                                         )                                     , &
                        &              kind=kind_dble                                                     &
                        &             )
                if (integralKernel==0.0d0) then
                   self%firstCrossingProbability(i,iTime)=0.0d0
                else
                   offsetEffective                       =+self%offsetEffective (self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),0.0_kind_quad,0.0_kind_quad,real(barrier(i),kind_quad),0.0_kind_quad,cosmologicalMassVariance_)
                   varianceResidual                      =+self%varianceResidual(self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),0.0_kind_quad                                                        ,cosmologicalMassVariance_)
                   self%firstCrossingProbability(i,iTime)=+max(                                                                            &
                        &                                      +0.0d0,                                                                     &
                        &                                      +(                                                                          &
                        &                                        +real(                                                                    &
                        &                                              Error_Function_Complementary(                                       &
                        &                                                                           +offsetEffective                       &
                        &                                                                           /sqrt(2.0_kind_quad*varianceResidual)  &
                        &                                                                          )                                     , &
                        &                                              kind=kind_dble                                                      &
                        &                                             )                                                                    &
                        &                                        /self%varianceStep                                                        &
                        &                                        -probabilityCrossingPrior                                                 &
                        &                                       )                                                                          &
                        &                                      /integralKernel                                                             &
                        &                                     )
                end if
             end do
             ! Force the probability at maximum variance to zero.
             self%firstCrossingProbability(self%countVariance,iTime)=0.0d0
          end do
          !$omp end do
          !![
          <objectDestructor name="excursionSetBarrier_"     />
          <objectDestructor name="cosmologicalMassVariance_"/>
	  !!]
          deallocate(barrier        )
          deallocate(barrierMidPoint)
          !$omp end parallel
          ! Update the variance table to reflect the variances at the midpoints. Note that the first crossing probability is computed
          ! at the mid-points. The last element of the variance table is unchanged to ensure that its value equals
          ! varianceMaximum. This will not affect the result because the probability at maximum variance is set to zero anyway.
          self%variance(1:self%countVariance-1)=varianceMidpoint(1:self%countVariance-1)
          deallocate(varianceMidpoint)
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayCounterClear(verbosityLevelWorking)
             call displayUnindent("done",verbosityLevelWorking)
#ifdef USEMPI
          end if
          if (self%coordinatedMPI_) then
             call mpiBarrier()
             self%firstCrossingProbability=mpiSelf%sum(self%firstCrossingProbability)
          end if
#endif
          ! Build the interpolators.
          if (allocated(self%interpolatorVariance)) deallocate(self%interpolatorVariance)
          if (allocated(self%interpolatorTime    )) deallocate(self%interpolatorTime    )
          allocate(self%interpolatorVariance)
          allocate(self%interpolatorTime    )
          self%interpolatorVariance=interpolator(self%variance,extrapolationType=extrapolationTypeFix)
          self%interpolatorTime    =interpolator(self%time    ,extrapolationType=extrapolationTypeFix)
          ! Record that the table is now built.
          self%tableInitialized=.true.
          ! Write the table to file if possible.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             if (self%useFile) then
                call File_Lock(char(self%fileName),fileLock,lockIsShared=.false.)
                call self%fileWrite()
                call File_Unlock(fileLock)
             end if
#ifdef USEMPI
          end if
#endif
       end if
       !$omp end critical(farahiMidpointProbabilityTabulate)
    end if
    ! Get interpolating factors.
    call self%interpolatorTime%linearFactors    (time    ,iTime    ,hTime    )
    call self%interpolatorVariance%linearFactors(variance,iVariance,hVariance)
    ! Compute first crossing probability by interpolating.
    farahiMidpointProbability=0.0d0
    do jTime=0,1
       do jVariance=0,1
          farahiMidpointProbability=+farahiMidpointProbability                                        &
               &                    +hTime                        (                            jTime) &
               &                    *hVariance                    (            jVariance            ) &
               &                    *self%firstCrossingProbability(iVariance-1+jVariance,iTime+jTime)
       end do
    end do
    return
  end function farahiMidpointProbability

  subroutine farahiMidpointRateTabulate(self,varianceProgenitor,time,node)
    !!{
    Tabulate the excursion set crossing rate.
    !!}
    use :: Display         , only : displayCounter              , displayCounterClear  , displayIndent       , displayMagenta  , &
          &                         displayMessage              , displayReset         , displayUnindent     , displayVerbosity, &
          &                         verbosityLevelWarn          , verbosityLevelWorking
    use :: Error_Functions , only : Error_Function_Complementary
    use :: File_Utilities  , only : File_Lock                   , File_Unlock          , lockDescriptor
    use :: Kind_Numbers    , only : kind_dble                   , kind_quad
    use :: MPI_Utilities   , only : mpiBarrier                  , mpiSelf
    use :: Numerical_Ranges, only : Make_Range                  , rangeTypeLinear      , rangeTypeLogarithmic
    use :: Table_Labels    , only : extrapolationTypeFix
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpoint), intent(inout)                   :: self
    double precision                                         , intent(in   )                   :: time                                         , varianceProgenitor
    type            (treeNode                               ), intent(inout)                   :: node
    double precision                                         , parameter                       :: varianceMinimumDefault    =1.0d-2
    double precision                                         , parameter                       :: varianceTolerance         =1.0d-6
    double precision                                         , parameter                       :: massLarge                 =1.0d16
    real            (kind=kind_quad                         ), allocatable  , dimension(:    ) :: firstCrossingRateQuad                        , varianceCurrentRateQuad , &
         &                                                                                        varianceProgenitorRateQuad                   , varianceMidpointRateQuad, &
         &                                                                                        barrierRateQuad                              , barrierMidpointRateQuad
    double precision                                         , allocatable  , dimension(:,:  ) :: nonCrossingRate
    double precision                                         , allocatable  , dimension(:,:,:) :: firstCrossingRate
    double precision                                                                           :: barrierRateTest
    class           (excursionSetBarrierClass               ), pointer                         :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass          ), pointer                         :: cosmologicalMassVariance_
    real            (kind=kind_quad                         ), parameter                       :: nonCrossingFractionTiny   =1.0e-002_kind_quad
    real            (kind=kind_quad                         ), parameter                       :: firstCrossingRateHuge     =1.0e+100_kind_quad
#ifdef USEMPI
    integer                                                                                    :: taskCount
#endif
    logical                                                                                    :: makeTable
    integer         (c_size_t                               )                                  :: loopCount                                    , loopCountTotal
    integer                                                                                    :: i                                            , iTime                   , &
         &                                                                                        iVariance                                    , j                       , &
         &                                                                                        countNewLower                                , countNewUpper           , &
         &                                                                                        countTimeNew                                 , iCompute                , &
         &                                                                                        countVarianceCurrentRate
    double precision                                                                           :: timeProgenitor                               , varianceMinimumRate     , &
         &                                                                                        massProgenitor                               , timeMinimumRate         , &
         &                                                                                        timeMaximumRate                              , varianceMaximumRateLimit
    character       (len=64                                 )                                  :: label
    type            (varying_string                         )                                  :: message                                      , reasonRemake
    type            (lockDescriptor                         )                                  :: fileLock
    real            (kind=kind_quad                         )                                  :: crossingFraction                             , effectiveBarrierInitial , &
         &                                                                                        probabilityCrossingPrior                     , varianceStepRate        , &
         &                                                                                        barrier                                      , integralKernelRate_     , &
         &                                                                                        growthFactorEffective                        , erfcArgumentNumerator   , &
         &                                                                                        erfcArgumentDenominator                      , erfcValue               , &
         &                                                                                        crossingFractionNew                          , varianceResidual        , &
         &                                                                                        offsetEffective
    logical                                                                                    :: varianceMaximumChanged
    
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
    !   S₁                = varianceCurrent
    !   S̃                = varianceProgenitor  +varianceCurrent
    !   S                 = varianceIntermediate+varianceCurrent
    !   B(Sᵢ)              = barrier(i)
    !   f(Sᵢ,t)            = self%firstCrossingProbability(i,iTime)
    !   Δδ[t,S₁,S₂,δ₁,δ₂] = self%offsetEffective         (self%time(iTime),0,S1,S2,0,barrier1,barrier2)
    !   ΔS[t,S₁,S₂]       = self%varianceResidual        (self%time(iTime),0,S1,S2                    )

    ! Note that the variables "varianceIntermediate" and "varianceProgenitor" are defined to be the variances in excess of S₁ - which is why they
    ! appear with "varianceCurrent" added to them in the above.
    !
    ! This function is used in the calculation of the distribution of δ at some S for trajectories originating from (S₁,δ₁) and
    ! which did not cross the barrier at any intermediate variance. As such suffixes in variable names have the following
    ! meanings:
    !
    !   "Current"      - refers to the current halo being considered for branching, i.e. the halo existing at point (S₁,δ₁);
    !   "Progenitor"   - refers to the potential progenitor halo being considered, i.e. the halo corresponding to some variance S > S₁;
    !   "Intermediate" - refers to the intermediate variance, S̃ (with S₁ < S̃ < S);
    !   "Quad"         - refers to a quantity computed in quad-precision;
    !   "Midpoint"     - refers to a quantity computed at the midpoint between two tabulated points.

    ! Determine if we need to make the table.
    !
    !! Read tables from file if possible. We make two passes through the logic that determines if the table needs to be remade. On
    !! the first pass if the table does need to be remade we attempt to read it from file. If the file is read, then we re-check if
    !! the table needs to be remade.
    do i=1,2
       makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
       if (i == 1 .and. self%useFile .and. makeTable) then
          call self%fileNameInitialize()
          call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       else
          exit
       end if
    end do
#ifdef USEMPI
    if (self%coordinatedMPI_) call mpiBarrier()
#endif
    if (makeTable.or.self%retabulateRateNonCrossing) then
       !$omp critical(farahiMidpointRateTabulate)
       ! Attempt to read the file again now that we are within the critical section. If another thread made the file while we were waiting we may be able to skip building the table.
       if (self%useFile) then
          call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       end if
       makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
       if (makeTable.or.self%retabulateRateNonCrossing) then
          ! Determine reason(s) for the remake.
          reasonRemake=''
          if (.not.self%tableInitializedRate) then
             reasonRemake=reasonRemake//' table does not exist;'
          else
             if (varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)) reasonRemake=reasonRemake//' progenitor variance exceeds previous maximum; '
             if (time               < self%timeMinimumRate                              ) reasonRemake=reasonRemake//' time exceeds previous minimum; '
             if (time               > self%timeMaximumRate                              ) reasonRemake=reasonRemake//' time exceeds previous minimum; '
             if (                     self%retabulateRateNonCrossing                    ) reasonRemake=reasonRemake//' non-crossing rates need to be retabulated;'
          end if
          ! Construct or expand the range of times to tabulate.
          countNewLower=0
          countNewUpper=0
          if (makeTable) then
             if (self%tableInitializedRate) then
                varianceMaximumChanged=varianceProgenitor > self%varianceMaximumRate
                timeMinimumRate=min(time/10.0d0**(2.0d0/dble(self%timeNumberPerDecade)),self%timeMinimumRate)
                timeMaximumRate=max(time*10.0d0**(2.0d0/dble(self%timeNumberPerDecade)),self%timeMaximumRate)
                ! Determine how many points the table must be extended by in each direction to span the new required range.
                if (self%timeMinimumRate > timeMinimumRate) countNewLower=int(+log10(self%timeMinimumRate/timeMinimumRate)*dble(self%timeNumberPerDecade)+1.0d0)
                if (self%timeMaximumRate < timeMaximumRate) countNewUpper=int(-log10(self%timeMaximumRate/timeMaximumRate)*dble(self%timeNumberPerDecade)+1.0d0)
                self%countTimeRate=self%countTimeRate+countNewLower+countNewUpper
                ! Adjust the limits of the table by an integer number of steps.
                self%timeMinimumRate=self%timeMinimumRate/10.0d0**(dble(countNewLower)/dble(self%timeNumberPerDecade))
                self%timeMaximumRate=self%timeMaximumRate*10.0d0**(dble(countNewUpper)/dble(self%timeNumberPerDecade))
             else
                varianceMaximumChanged=.true.
                self%timeMinimumRate  =time/10.0d0**(2.0d0/dble(self%timeNumberPerDecade))
                self%timeMaximumRate  =time*10.0d0**(2.0d0/dble(self%timeNumberPerDecade))
                self%countTimeRate    =max(int(log10(self%timeMaximumRate/self%timeMinimumRate)*dble(self%timeNumberPerDecade))+2,2)
                ! Ensure the maximum of the table is precisely an integer number of steps above the minimum.
                self%timeMaximumRate   =self%timeMinimumRate*10.0d0**(dble(self%countTimeRate-1)/dble(self%timeNumberPerDecade))
             end if
             ! Set the default minimum variance.
             varianceMinimumRate=varianceMinimumDefault
             ! Next reduce the variance if necessary such that the typical amplitude of fluctuations is less (by a factor of 10) than
             ! the effective barrier height at zero variance for the minimum and maximum times that we must consider. We use some
             ! suitably large mass to estimate the growth of fluctuations on large scales (since we can't assume infinitely large
             ! scales).
             allocate(excursionSetBarrier_     ,mold=self%excursionSetBarrier_     )
             allocate(cosmologicalMassVariance_,mold=self%cosmologicalMassVariance_)
             !$omp critical(excursionSetsSolverFarahiMidpointDeepCopy)
             !![
             <deepCopyReset variables="self%excursionSetBarrier_ self%cosmologicalMassVariance_"/>
             <deepCopy source="self%excursionSetBarrier_"      destination="excursionSetBarrier_"     />
             <deepCopy source="self%cosmologicalMassVariance_" destination="cosmologicalMassVariance_"/>
             <deepCopyFinalize variables="excursionSetBarrier_ cosmologicalMassVariance_"/>
             !!]
             !$omp end critical(excursionSetsSolverFarahiMidpointDeepCopy)
             growthFactorEffective=+cosmologicalMassVariance_%rootVariance(massLarge,self%timeMaximumRate                                ) &
                  &                /cosmologicalMassVariance_%rootVariance(massLarge,self%timeMaximumRate*(1.0d0-self%fractionalTimeStep))
             varianceMinimumRate  =min(                                                                                                                      &
                  &                    +varianceMinimumRate                                                                                                , &
                  &                    +1.0d-2                                                                                                               &
                  &                    *(                                                                                                                    &
                  &                      +excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate*(1.0d0-self%fractionalTimeStep),node,rateCompute=.true.)  &
                  &                      *dble(growthFactorEffective)                                                                                        &
                  &                      -excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate                                ,node,rateCompute=.true.)  &
                  &                     )**2                                                                                                                 &
                  &                   )
             !![
             <objectDestructor name="excursionSetBarrier_"     />
             <objectDestructor name="cosmologicalMassVariance_"/>
             !!]
             self%varianceMaximumRate=self%varianceLimit(varianceProgenitor)
             self%countVarianceProgenitorRate        =int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecade           ))+1
             self%countVarianceCurrentRate           =int(self%varianceMaximumRate*dble(self%varianceNumberPerUnit))
             self%countVarianceCurrentRateNonCrossing=int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecadeNonCrossing))+1
          else
             varianceMaximumChanged=.false.
             ! The progenitor variance table stores the values at the mid-points, thus a factor of two is added.
             varianceMinimumRate   =2.0d0*self%varianceProgenitorRate(1)
          end if
          ! Store copies of the current tables if these will be used later.
          if (.not.varianceMaximumChanged.and.     makeTable                     ) then
             call move_alloc(self%firstCrossingRate,firstCrossingRate)
          else
             allocate(firstCrossingRate(0,0,0))
          end if
          if (.not.varianceMaximumChanged.and..not.self%retabulateRateNonCrossing) then
             call move_alloc(self%nonCrossingRate  ,nonCrossingRate  )
          else
             allocate(  nonCrossingRate(  0,0))
          end if
          if (makeTable) then
             if (allocated(self%firstCrossingRate)) deallocate(self%firstCrossingRate)
             allocate(self%firstCrossingRate(0:self%countVarianceProgenitorRate,0:self%countVarianceCurrentRate,self%countTimeRate))
          end if
          if (allocated(self%varianceProgenitorRate        )) deallocate(self%varianceProgenitorRate        )
          if (allocated(self%varianceCurrentRate           )) deallocate(self%varianceCurrentRate           )
          if (allocated(self%varianceCurrentRateNonCrossing)) deallocate(self%varianceCurrentRateNonCrossing)
          if (allocated(self%timeRate                      )) deallocate(self%timeRate                      )
          if (allocated(self%nonCrossingRate               )) deallocate(self%nonCrossingRate               )
          allocate(self%varianceProgenitorRate        (0:self%countVarianceProgenitorRate                                                              ))
          allocate(self%varianceCurrentRate           (                                   0:self%countVarianceCurrentRate                              ))
          allocate(self%varianceCurrentRateNonCrossing(                                   0:self%countVarianceCurrentRateNonCrossing                   ))
          allocate(self%timeRate                      (                                                                              self%countTimeRate))
          allocate(self%nonCrossingRate               (                                   0:self%countVarianceCurrentRateNonCrossing,self%countTimeRate))
          ! For the variance table, the zeroth point is always zero, higher points are distributed uniformly in variance.
          self%varianceProgenitorRate        (0                                         )=0.0d0
          self%varianceProgenitorRate        (1:self%countVarianceProgenitorRate        )=self%varianceRange(varianceMinimumRate ,self%varianceMaximumRate,self%countVarianceProgenitorRate          ,exponent =1.0d0               )
          self%varianceCurrentRate           (0:self%countVarianceCurrentRate           )=Make_Range        (0.0d0               ,self%varianceMaximumRate,self%countVarianceCurrentRate           +1,rangeType=rangeTypeLinear     )
          self%varianceCurrentRateNonCrossing(0                                         )=0.0d0
          self%varianceCurrentRateNonCrossing(1:self%countVarianceCurrentRateNonCrossing)=self%varianceRange(varianceMinimumRate ,self%varianceMaximumRate,self%countVarianceCurrentRateNonCrossing  ,exponent =1.0d0               )
          ! The time table is logarithmically distributed in time.
          self%timeRate                                                                  =Make_Range        (self%timeMinimumRate,self%timeMaximumRate    ,self%countTimeRate                        ,rangeType=rangeTypeLogarithmic)
          ! Allocate temporary arrays used in quad-precision solver for barrier crossing rates.
          allocate(varianceProgenitorRateQuad(0:self%countVarianceProgenitorRate))
          varianceProgenitorRateQuad=self%varianceProgenitorRate
          ! Compute the variance at the mid-points.
          allocate(varianceMidpointRateQuad(0:self%countVarianceProgenitorRate))
          varianceMidpointRateQuad(0)=0.0_kind_quad
          forall(i=1:self%countVarianceProgenitorRate)
             varianceMidpointRateQuad(i)=(varianceProgenitorRateQuad(i-1)+varianceProgenitorRateQuad(i))/2.0_kind_quad
          end forall
          ! Loop through the table and solve for the first crossing distribution.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayIndent("solving for excursion set barrier crossing rates",verbosityLevelWorking)
             message="reason(s):"//reasonRemake
             call displayMessage(message,verbosityLevelWorking)
             message="     time: "
             write (label,'(f6.3)') self%timeMinimumRate
             message=message//trim(label)//" to "
             write (label,'(f6.3)') self%timeMaximumRate
             message=message//trim(label)
             call displayMessage(message,verbosityLevelWorking)
             message=" variance: "
             write (label,'(f9.3)') self%varianceMaximumRate
             message=message//trim(label)
             call displayMessage(message,verbosityLevelWorking)
#ifdef USEMPI
          end if
#endif
          countTimeNew=self%countTimeRate
          if (.not.varianceMaximumChanged.and..not.self%retabulateRateNonCrossing) countTimeNew=countNewLower+countNewUpper
          loopCountTotal   = int(countTimeNew,kind=c_size_t)*int(self%countVarianceCurrentRateNonCrossing+1,kind=c_size_t)
          if (makeTable) then
             countTimeNew=self%countTimeRate
             if (.not.varianceMaximumChanged) countTimeNew=countNewLower+countNewUpper
             loopCountTotal=+loopCountTotal                                                                                &
                  &         +int(countTimeNew,kind=c_size_t)*int(self%countVarianceCurrentRate           +1,kind=c_size_t)
          end if
#ifdef USEMPI
          if (mpiSelf%isMaster() .and. self%coordinatedMPI_) then
             loopCountTotal=loopCountTotal/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t
          end if
#endif
          loopCount=0
#ifdef USEMPI
          if (self%coordinatedMPI_) then
             if (makeTable) then
                self%firstCrossingRate=0.0d0
             end if
             self%nonCrossingRate=0.0d0
          end if
          taskCount=-1
#endif
          ! Make a call to the barrier function at maximum variance for the minimum and maximum times so that the barrier function
          ! is initialized and covers the whole range in which we are interested.
          barrierRateTest=self%excursionSetBarrier_%barrier(self%varianceMaximumRate,self%timeMinimumRate*(1.0d0-self%fractionalTimeStep),node,rateCompute=.true.)
          barrierRateTest=self%excursionSetBarrier_%barrier(self%varianceMaximumRate,self%timeMaximumRate                                ,node,rateCompute=.true.)
          ! Enter an OpenMP parallel region. Each parallel thread will solve for the first crossing rate at a different epoch.
          !$omp parallel private(iTime,timeProgenitor,iVariance,varianceStepRate,i,j,iCompute,probabilityCrossingPrior,integralKernelRate_,crossingFraction,crossingFractionNew,barrier,effectiveBarrierInitial,firstCrossingRateQuad,excursionSetBarrier_,cosmologicalMassVariance_,barrierRateQuad,barrierMidpointRateQuad,varianceCurrentRateQuad,massProgenitor,growthFactorEffective,offsetEffective,varianceResidual,erfcArgumentNumerator,erfcArgumentDenominator,erfcValue,message,label,varianceMaximumRateLimit) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
          allocate(excursionSetBarrier_     ,mold=self%excursionSetBarrier_     )
          allocate(cosmologicalMassVariance_,mold=self%cosmologicalMassVariance_)
          !$omp critical(excursionSetsSolverFarahiMidpointDeepCopy)
          !![
          <deepCopyReset variables="self%excursionSetBarrier_ self%cosmologicalMassVariance_"/>
          <deepCopy source="self%excursionSetBarrier_"      destination="excursionSetBarrier_"     />
          <deepCopy source="self%cosmologicalMassVariance_" destination="cosmologicalMassVariance_"/>
          <deepCopyFinalize variables="excursionSetBarrier_ cosmologicalMassVariance_"/>
          !!]
          !$omp end critical(excursionSetsSolverFarahiMidpointDeepCopy)
          allocate(barrierRateQuad        (self%countVarianceProgenitorRate))
          allocate(barrierMidpointRateQuad(self%countVarianceProgenitorRate))
          ! In the first run, first crossing rates are computed. In the second run, non-crossing rates are computed at different
          ! grid points.
          do iCompute=1,2
             if (iCompute == 1) then
                if (.not.makeTable) cycle
                countVarianceCurrentRate=self%countVarianceCurrentRate
                allocate(varianceCurrentRateQuad(0:self%countVarianceCurrentRate           ))
                varianceCurrentRateQuad =self%varianceCurrentRate
             else
                countVarianceCurrentRate=self%countVarianceCurrentRateNonCrossing
                if (allocated(varianceCurrentRateQuad)) deallocate(varianceCurrentRateQuad)
                allocate(varianceCurrentRateQuad(0:self%countVarianceCurrentRateNonCrossing))
             end if
             do iTime=1,self%countTimeRate
                ! Skip if this time was already computed.
                if (iCompute == 1) then
                   if (.not.varianceMaximumChanged                                        .and.(iTime > countNewLower .and. self%countTimeRate+1-iTime > countNewUpper)) cycle
                   varianceMaximumRateLimit=self%varianceMaximumRate
                else
                   if (.not.varianceMaximumChanged.and..not.self%retabulateRateNonCrossing.and.(iTime > countNewLower .and. self%countTimeRate+1-iTime > countNewUpper)) cycle
                   if (self%massMinimumRateNonCrossing > 0.0d0) then
                      varianceMaximumRateLimit=cosmologicalMassVariance_%rootVariance(self%massMinimumRateNonCrossing,self%timeRate(iTime))**2
                   else
                      varianceMaximumRateLimit=self%varianceMaximumRate
                   end if
                   ! For computing non-crossing rates, the results are tabulated with respect to $S_{\rm max}-S$ so that interpolation
                   ! is more accurate when $S$ approaches $S_{\rm max}$.
                   do iVariance=0,countVarianceCurrentRate
                      varianceCurrentRateQuad(iVariance)=max(varianceMaximumRateLimit-self%varianceCurrentRateNonCrossing(iVariance),0.0d0)
                   end do
                end if
                ! Allocate workspace table.
                if (.not.allocated(firstCrossingRateQuad)) allocate(firstCrossingRateQuad(0:self%countVarianceProgenitorRate))
                ! Compute a suitable progenitor time.
                timeProgenitor=self%timeRate(iTime)*(1.0d0-self%fractionalTimeStep)
                ! Iterate over variances of the current halo.
                !$omp do schedule(dynamic)
                do iVariance=0,countVarianceCurrentRate
#ifdef USEMPI
                   taskCount=taskCount+1
                   if (self%coordinatedMPI_ .and. mod(taskCount,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
#ifdef USEMPI
                   if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
                      call displayCounter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityLevelWorking)
#ifdef USEMPI
                   end if
#endif
                   !$omp atomic
                   loopCount=loopCount+1_c_size_t
                   ! Construct the barrier table as a function of the progenitor variance.
                   do i=1,self%countVarianceProgenitorRate
                      massProgenitor            =+cosmologicalMassVariance_%mass        (real(sqrt(+varianceProgenitorRateQuad(i)+varianceCurrentRateQuad(iVariance)),kind=8),self%timeRate      (iTime)                        )
                      growthFactorEffective     =+cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,self%timeRate      (iTime)                        ) &
                           &                     /cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,     timeProgenitor                               )
                      barrierRateQuad        (i)=+excursionSetBarrier_     %barrier     (real(     +varianceProgenitorRateQuad(i)+varianceCurrentRateQuad(iVariance) ,kind=8),     timeProgenitor       ,node,rateCompute=.true.) &
                           &                     *growthFactorEffective
                      massProgenitor            =+cosmologicalMassVariance_%mass        (real(sqrt(+varianceMidpointRateQuad  (i)+varianceCurrentRateQuad(iVariance)),kind=8),self%timeRate      (iTime)                        )
                      growthFactorEffective     =+cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,self%timeRate      (iTime)                        ) &
                           &                     /cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,     timeProgenitor                               )
                      barrierMidpointRateQuad(i)=+excursionSetBarrier_     %barrier     (real(     +varianceMidpointRateQuad  (i)+varianceCurrentRateQuad(iVariance),kind=8 ),     timeProgenitor       ,node,rateCompute=.true.) &
                           &                     *growthFactorEffective
                   end do
                   ! For zero variance, the rate is initialized to zero.
                   firstCrossingRateQuad(0)=0.0_kind_quad
                   ! Compute the step in variance across this first grid point.
                   varianceStepRate=varianceProgenitorRateQuad(1)-varianceProgenitorRateQuad(0)
                   ! Compute the barrier for the descendant.
                   barrier=real(excursionSetBarrier_%barrier(real(varianceCurrentRateQuad(iVariance),kind=8),self%timeRate(iTime),node,rateCompute=.true.),kind=kind_quad)
                   ! Compute the first crossing distribution at the first grid point.
                   if (varianceProgenitorRateQuad(1)+varianceCurrentRateQuad(iVariance) >= varianceMaximumRateLimit) then
                      firstCrossingRateQuad(1)= 0.0_kind_quad
                   else
                      offsetEffective         =+self%offsetEffective (self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(1),varianceMidpointRateQuad(1),0.0_kind_quad,barrierRateQuad(1),barrierMidpointRateQuad(1),cosmologicalMassVariance_)
                      varianceResidual        =+self%varianceResidual(self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(1),varianceMidpointRateQuad(1)                                                            ,cosmologicalMassVariance_)
                      integralKernelRate_     =+integralKernelRate(varianceResidual,offsetEffective)
                      ! If the integral kernel is zero (to machine precision) then simply assume no crossing rate.
                      if (integralKernelRate_ <= 0.0d0) then
                         firstCrossingRateQuad=0.0d0
                         cycle
                      end if
                      ! Compute the first crossing rate at the first grid point.
                      offsetEffective         =+self%offsetEffective (self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(1),0.0_kind_quad,0.0_kind_quad,barrierRateQuad(1),barrier,cosmologicalMassVariance_)
                      varianceResidual        =+self%varianceResidual(self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(1),0.0_kind_quad                                         ,cosmologicalMassVariance_)
                      firstCrossingRateQuad(1)=+integralKernelRate(varianceResidual,offsetEffective) &
                           &                   /varianceStepRate                                     &
                           &                   /integralKernelRate_
                   end if
                   varianceStepRate=+varianceProgenitorRateQuad(1) &
                        &           -varianceProgenitorRateQuad(0)
                   crossingFraction=+firstCrossingRateQuad     (1) &
                        &           *varianceStepRate
                   ! Iterate over remaining progenitor variances.
                   do i=2,self%countVarianceProgenitorRate
                      if (varianceProgenitorRateQuad(i)+varianceCurrentRateQuad(iVariance) >= varianceMaximumRateLimit) then
                         firstCrossingRateQuad(i)=0.0_kind_quad
                      else
                         effectiveBarrierInitial=+barrierRateQuad(i) &
                              &                  -barrier
                         offsetEffective        =+self%offsetEffective (self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),varianceMidpointRateQuad(i),barrier,effectiveBarrierInitial,barrierMidpointRateQuad(i)-barrier,cosmologicalMassVariance_)
                         varianceResidual       =+self%varianceResidual(self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),varianceMidpointRateQuad(i)                                                                   ,cosmologicalMassVariance_)
                         integralKernelRate_=integralKernelRate(varianceResidual,offsetEffective)
                         if (integralKernelRate_ == 0.0_kind_quad) then
                            firstCrossingRateQuad(i)=0.0_kind_quad
                         else
                            ! Iterate over all smaller variances, computing the contribution from trajectories that crossed the barrier at those variances.
                            probabilityCrossingPrior=+0.0_kind_quad
                            do j=1,i-1
                               offsetEffective        =self%offsetEffective (self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),varianceMidpointRateQuad(j),barrier,effectiveBarrierInitial,barrierMidpointRateQuad(j)-barrier,cosmologicalMassVariance_)
                               varianceResidual       =self%varianceResidual(self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),varianceMidpointRateQuad(j)                                                                   ,cosmologicalMassVariance_)
                               varianceStepRate        =+varianceProgenitorRateQuad(j  )                      &
                                    &                   -varianceProgenitorRateQuad(j-1)
                               probabilityCrossingPrior=+probabilityCrossingPrior                             &
                                    &                   +firstCrossingRateQuad     (j  )                      &
                                    &                   *varianceStepRate                                     &
                                    &                   *integralKernelRate(varianceResidual,offsetEffective)
                            end do
                            varianceStepRate        =+varianceProgenitorRateQuad(i  ) &
                                 &                   -varianceProgenitorRateQuad(i-1)
                            offsetEffective         =+self%offsetEffective (self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),0.0_kind_quad,barrier,effectiveBarrierInitial,0.0_kind_quad,cosmologicalMassVariance_)
                            varianceResidual        =+self%varianceResidual(self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),0.0_kind_quad                                              ,cosmologicalMassVariance_)
                            firstCrossingRateQuad(i)=+integralKernelRate(varianceResidual,offsetEffective) &
                                 &                   -probabilityCrossingPrior
                            if     (                                                                                                                          &
                                 &   firstCrossingRateQuad(i)                                                                   > 0.0d0                       &
                                 &  .and.                                                                                                                     &
                                 &   integralKernelRate_                                                                        > 0.0d0                       &
                                 &  .and.                                                                                                                     &
                                 &   exponent(firstCrossingRateQuad(i))-exponent(varianceStepRate)-exponent(integralKernelRate_) < maxExponent(0.0_kind_quad) &
                                 & ) then
                               firstCrossingRateQuad(i)=+firstCrossingRateQuad(i) &
                                    &                   /varianceStepRate         &
                                    &                   /integralKernelRate_
                            else
                               firstCrossingRateQuad(i)=0.0d0
                            end if
                            ! Accumulate the crossing fraction for use in the following check.
                            varianceStepRate   =+varianceProgenitorRateQuad(i  ) &
                                 &              -varianceProgenitorRateQuad(i-1)
                            crossingFractionNew=+crossingFraction                &
                                 &              +firstCrossingRateQuad     (i  ) &
                                 &              *varianceStepRate
                            ! Remove unphysical values. Force the crossing rate at points close to maximum variance, or where most
                            ! trajectories have already crossed the barrier to zero if its value is one order of magnitude larger
                            ! than the value at previous point.
                            if     (                                                                                                                              &
                                 &  (                                                                                                                             &
                                 &    varianceProgenitorRateQuad(i)+varianceCurrentRateQuad(iVariance) > varianceMaximumRateLimit-10.0_kind_quad*varianceStepRate &
                                 &   .or.                                                                                                                         &
                                 &    crossingFractionNew                                              > (1.0_kind_quad-nonCrossingFractionTiny)                  &
                                 &   )                                                                                                                            &
                                 &  .and.                                                                                                                         &
                                 &  (                                                                                                                             &
                                 &   firstCrossingRateQuad      (i)                                    > 10.0_kind_quad*firstCrossingRateQuad(i-1)                &
                                 &   .or.                                                                                                                         &
                                 &   firstCrossingRateQuad      (i)                                    >                firstCrossingRateHuge                     &
                                 &   .or.                                                                                                                         &
                                 &   firstCrossingRateQuad      (i)                                    <  0.0_kind_quad                                           &
                                 &   )                                                                                                                            &
                                 & )                                                                                                                              &
                                 & firstCrossingRateQuad(i)=0.0_kind_quad
                            if (abs(firstCrossingRateQuad(i)) > firstCrossingRateHuge .and. displayVerbosity() >= verbosityLevelWarn) then
                               message=         displayMagenta()//"WARNING:"//displayReset()//" unphysical solution for crossing rate:"//char(10)
                               write (label,'(e15.6," (",e15.6,")")') crossingFraction,crossingFractionNew
                               message=message//"    crossing fraction = "//trim(label)//char(10)
                               write (label,'(e15.6)') firstCrossingRateQuad(i)
                               message=message//"  first crossing rate = "//trim(label)//char(10)
                               write (label,'(i8)'   ) iTime
                               message=message//"           index time = "//trim(label)//char(10)
                               write (label,'(i8)'   ) iVariance
                               message=message//"       index variance = "//trim(label)//char(10)
                               write (label,'(i8)'   ) i
                               message=message//"                index = "//trim(label)
                               call displayMessage(message,verbosityLevelWarn)
                            end if
                            crossingFraction     =+crossingFractionNew
                         end if
                      end if
                   end do
                   if (iCompute == 1) then
                      ! Store the compute crossing rate in our table.
                      self%firstCrossingRate(:,iVariance,iTime)=real(firstCrossingRateQuad,kind=kind_dble)
                      ! Divide through by the time step to get the rate of barrier crossing.
                      self%firstCrossingRate(:,iVariance,iTime)=+self%firstCrossingRate (:,iVariance,iTime) &
                           &                                    /self%timeRate          (            iTime) &
                           &                                    /self%fractionalTimeStep
                   else
                      ! Compute the fraction of trajectories which never cross the barrier.
                      crossingFraction=0.0_kind_quad
                      do j=1,self%countVarianceProgenitorRate
                         if (varianceCurrentRateQuad(iVariance)+varianceProgenitorRateQuad(j) <= varianceMaximumRateLimit) then
                            varianceStepRate=+varianceProgenitorRateQuad(j  ) &
                                 &           -varianceProgenitorRateQuad(j-1)
                            crossingFraction=+crossingFraction                &
                                 &           +firstCrossingRateQuad     (j  ) &
                                 &           *varianceStepRate
                         end if
                      end do
                      ! Compute the rate for trajectories which never cross the barrier.
                      self%nonCrossingRate(iVariance,iTime)=real(                                     &
                           &                                     +max(                                &
                           &                                          1.0_kind_quad-crossingFraction, &
                           &                                          0.0_kind_quad                   &
                           &                                         )                                &
                           &                                     /self%timeRate(iTime)                &
                           &                                     /self%fractionalTimeStep           , &
                           &                                     kind=kind_dble                       &
                           &                                    )
                   end if
                end do
                !$omp end do
             end do
          end do
          !![
          <objectDestructor name="excursionSetBarrier_"     />
          <objectDestructor name="cosmologicalMassVariance_"/>
          !!]
          ! Deallocate work arrays.
          deallocate(barrierRateQuad        )
          deallocate(barrierMidpointRateQuad)
          deallocate(varianceCurrentRateQuad)
          if (allocated(firstCrossingRateQuad)) deallocate(firstCrossingRateQuad)
          !$omp end parallel
          ! Update the variance table to reflect the variances at the midpoints. Note that the first crossing probability is computed
          ! at the mid-points. The last element of the variance table is unchanged to ensure that its value equals
          ! varianceMaximum. This will not affect the result because the crossing rate at maximum variance is set to zero anyway.
          self%varianceProgenitorRate(1:self%countVarianceProgenitorRate-1)=real(varianceMidpointRateQuad(1:self%countVarianceProgenitorRate-1),kind=kind_dble)
          deallocate(varianceMidpointRateQuad  )
          deallocate(varianceProgenitorRateQuad)
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayCounterClear(       verbosityLevelWorking)
             call displayUnindent    ("done",verbosityLevelWorking)
#ifdef USEMPI
          end if
          if (self%coordinatedMPI_) then
             call mpiBarrier()
             if (makeTable) then
                self%firstCrossingRate=mpiSelf%sum(self%firstCrossingRate)
             end if
             self%nonCrossingRate=mpiSelf%sum(self%nonCrossingRate)
          end if
#endif
          ! If only times have changed then copy results previously computed to the tables.
          if (.not.varianceMaximumChanged) then
             if (makeTable) then
                self%firstCrossingRate(:,:,countNewLower+1:countNewLower+size(firstCrossingRate,dim=3))=firstCrossingRate
                deallocate(firstCrossingRate)
             end if
             if (.not.self%retabulateRateNonCrossing) then
                self%nonCrossingRate  (  :,countNewLower+1:countNewLower+size(nonCrossingRate  ,dim=2))=nonCrossingRate
                deallocate(nonCrossingRate  )
             end if
          end if
          ! Build the interpolators.
          if (allocated(self%interpolatorVarianceRate                  )) deallocate(self%interpolatorVarianceRate                  )
          if (allocated(self%interpolatorVarianceCurrentRate           )) deallocate(self%interpolatorVarianceCurrentRate           )
          if (allocated(self%interpolatorVarianceCurrentRateNonCrossing)) deallocate(self%interpolatorVarianceCurrentRateNonCrossing)
          if (allocated(self%interpolatorTimeRate                      )) deallocate(self%interpolatorTimeRate                      )
          allocate(self%interpolatorVarianceRate                  )
          allocate(self%interpolatorVarianceCurrentRate           )
          allocate(self%interpolatorVarianceCurrentRateNonCrossing)
          allocate(self%interpolatorTimeRate                      )
          self%interpolatorVarianceRate                  =interpolator(self%varianceProgenitorRate        ,extrapolationType=extrapolationTypeFix)
          self%interpolatorVarianceCurrentRate           =interpolator(self%varianceCurrentRate           ,extrapolationType=extrapolationTypeFix)
          self%interpolatorVarianceCurrentRateNonCrossing=interpolator(self%varianceCurrentRateNonCrossing,extrapolationType=extrapolationTypeFix)
          self%interpolatorTimeRate                      =interpolator(self%timeRate                      ,extrapolationType=extrapolationTypeFix)
          ! Set previous variance and time to unphysical values to force recompute of interpolation factors on next call.
          self%variancePreviousRate=-1.0d0
          self%timePreviousRate    =-1.0d0
          ! Record that the table is now built.
          self%tableInitializedRate     =.true.
          self%retabulateRateNonCrossing=.false.
          ! Write the table to file if possible.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             if (self%useFile) then
                call File_Lock(char(self%fileName),fileLock,lockIsShared=.false.)
                call self%fileWrite()
                call File_Unlock(fileLock)
             end if
#ifdef USEMPI
          end if
#endif
       end if
       !$omp end critical(farahiMidpointRateTabulate)
    end if
    return
  end subroutine farahiMidpointRateTabulate

  function integralKernelRate(varianceResidual,offsetEffective)
    !!{
    Compute the integral kernel rate.
    !!}
    use :: Error_Functions, only : Error_Function_Complementary
    implicit none
    real(kind_quad)                :: integralKernelRate
    real(kind_quad), intent(in   ) :: varianceResidual  , offsetEffective
    real(kind_quad)                :: denominator

    if      (varianceResidual <  0.0_kind_quad) then
       integralKernelRate=0.0_kind_quad
    else if (varianceResidual == 0.0_kind_quad) then
       ! Zero residual variance - the first crossing rate is either 0 or 1, depending on the sign of the offset.
       if (offsetEffective > 0.0_kind_quad) then
          integralKernelRate=0.0_kind_quad
       else
          integralKernelRate=2.0_kind_quad
       end if
    else
       denominator=sqrt(2.0_kind_quad*varianceResidual)
       if (offsetEffective == 0.0_kind_quad .or. exponent(offsetEffective)-exponent(denominator) > maxExponent(0.0_kind_quad)) then
          integralKernelRate=1.0_kind_quad
       else
          integralKernelRate=Error_Function_Complementary(                 &
               &                                          +offsetEffective &
               &                                          /denominator     &
               &                                         )
       end if
    end if
    return
  end function integralKernelRate
