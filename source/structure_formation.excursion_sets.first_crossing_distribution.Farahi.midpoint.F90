!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method to perform the integrations \citep{du_substructure_2017}.
!!}

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingFarahiMidpoint">
   <description>An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method to perform the integrations \citep{du_substructure_2017}.</description>
  </excursionSetFirstCrossing>
  !!]
  type, extends(excursionSetFirstCrossingFarahi) :: excursionSetFirstCrossingFarahiMidpoint
     !!{
     An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method to perform the integrations \citep{du_substructure_2017}.
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

  function farahiMidpointConstructorParameters(parameters) result(self)
    !!{
    Constructor for the Farahi-midpoint excursion set class first crossing class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(excursionSetFirstCrossingFarahiMidpoint)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self%excursionSetFirstCrossingFarahi=excursionSetFirstCrossingFarahi(parameters)
    return
  end function farahiMidpointConstructorParameters

  function farahiMidpointConstructorInternal(timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the Farahi-midpoint excursion set class first crossing class.
    !!}
    implicit none
    type            (excursionSetFirstCrossingFarahiMidpoint)                        :: self
    double precision                                         , intent(in   )         :: timeStepFractional
    integer                                                  , intent(in   )         :: varianceNumberPerUnitProbability, varianceNumberPerUnit  , &
         &                                                                              timeNumberPerDecade             , varianceNumberPerDecade
    type            (varying_string                         ), intent(in   )         :: fileName
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetBarrierClass               ), intent(in   ), target :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass          ), intent(in   ), target :: cosmologicalMassVariance_

    self%excursionSetFirstCrossingFarahi=excursionSetFirstCrossingFarahi(timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_)
    return
  end function farahiMidpointConstructorInternal

  double precision function farahiMidpointProbability(self,variance,time,node)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    use :: Display          , only : displayCounter , displayCounterClear  , displayIndent       , displayMessage, &
          &                          displayUnindent, verbosityLevelWorking
    use :: Error_Functions  , only : erfApproximate
    use :: File_Utilities   , only : File_Lock      , File_Unlock          , lockDescriptor
    use :: Kind_Numbers     , only : kind_dble      , kind_quad
    use :: MPI_Utilities    , only : mpiBarrier     , mpiSelf
    use :: Memory_Management, only : allocateArray  , deallocateArray
    use :: Numerical_Ranges , only : Make_Range     , rangeTypeLinear      , rangeTypeLogarithmic
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpoint), intent(inout)                 :: self
    double precision                                         , intent(in   )                 :: variance                     , time
    type            (treeNode                               ), intent(inout)                 :: node
    double precision                                                        , dimension(0:1) :: hTime                        , hVariance
    double precision                                         , parameter                     :: varianceTableTolerance=1.0d-6
    double precision                                         , allocatable  , dimension( : ) :: varianceMidTable             , barrierTable  , &
         &                                                                                      barrierMidTable
    double precision                                                                         :: barrierTest
    class           (excursionSetBarrierClass               ), pointer                       :: excursionSetBarrier_
    logical                                                                                  :: makeTable
    integer         (c_size_t                               )                                :: iTime                        , iVariance     , &
         &                                                                                      loopCount                    , loopCountTotal, &
         &                                                                                      i                            , j             , &
         &                                                                                      jTime                        , jVariance
    double precision                                                                         :: sigma1f
    real            (kind=kind_quad                         )                                :: integralKernel
    character       (len =9                                 )                                :: label
    type            (varying_string                         )                                :: message
    type            (lockDescriptor                         )                                :: fileLock

    ! Read tables from file if possible.
    do i=1,2
       makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+varianceTableTolerance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
       if (i == 1 .and. self%useFile .and. makeTable) then
          call self%fileNameInitialize()
          call File_Lock(char(self%fileName),fileLock)
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
       makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+varianceTableTolerance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
       if (makeTable) then
          ! Construct the table of variance on which we will solve for the first crossing distribution.
          if (allocated(self%varianceTable                )) call deallocateArray(self%varianceTable                )
          if (allocated(self%timeTable                    )) call deallocateArray(self%timeTable                    )
          if (allocated(self%firstCrossingProbabilityTable)) call deallocateArray(self%firstCrossingProbabilityTable)
          self%varianceMaximum   =max(self%varianceMaximum,variance)
          self%varianceTableCount=int(self%varianceMaximum*dble(self%varianceNumberPerUnitProbability))
          if (self%tableInitialized) then
             self%timeMinimum=min(self%timeMinimum,time/10.0d0**(2.0d0/dble(self%timeNumberPerDecade)))
             self%timeMaximum=max(self%timeMaximum,time*10.0d0**(2.0d0/dble(self%timeNumberPerDecade)))
          else
             self%timeMinimum=                     time/10.0d0**(2.0d0/dble(self%timeNumberPerDecade))
             self%timeMaximum=                     time*10.0d0**(2.0d0/dble(self%timeNumberPerDecade))
          end if
          self%timeTableCount=max(2,int(log10(self%timeMaximum/self%timeMinimum)*dble(self%timeNumberPerDecade))+1)
          call allocateArray(self%varianceTable                ,[1+self%varianceTableCount                    ],lowerBounds=[0  ])
          call allocateArray(self%timeTable                    ,[                          self%timeTableCount]                  )
          call allocateArray(self%firstCrossingProbabilityTable,[1+self%varianceTableCount,self%timeTableCount],lowerBounds=[0,1])
          call allocateArray(     varianceMidTable             ,[1+self%varianceTableCount                    ],lowerBounds=[0  ])
          self%timeTable        =Make_Range(self%timeMinimum,self%timeMaximum    ,self%timeTableCount      ,rangeType=rangeTypeLogarithmic)
          self%varianceTable    =Make_Range(0.0d0           ,self%varianceMaximum,self%varianceTableCount+1,rangeType=rangeTypeLinear     )
          self%varianceTableStep=self%varianceTable(1)-self%varianceTable(0)
          ! Compute the variance at the mid-points.
          varianceMidTable(0)=0.0d0
          forall(i=1:self%varianceTableCount)
             varianceMidTable(i)=(self%varianceTable(i-1)+self%varianceTable(i))/2.0d0
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
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
             loopCountTotal=(int(self%timeTableCount,kind=c_size_t)/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t)*(int(self%varianceTableCount-1,kind=c_size_t)*int(self%varianceTableCount,kind=c_size_t))/2_c_size_t
          else
#endif
             loopCountTotal= int(self%timeTableCount,kind=c_size_t)                                               *(int(self%varianceTableCount-1,kind=c_size_t)*int(self%varianceTableCount,kind=c_size_t))/2_c_size_t
#ifdef USEMPI
          end if
#endif
          loopCount=0
#ifdef USEMPI
          if (self%coordinatedMPI_) self%firstCrossingProbabilityTable=0.0d0
#endif
          ! Make a call to the barrier function at maximum variance for the minimum and maximum times so that the barrier function
          ! is initialized and covers the whole range we are intereseted in.
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMinimum,node,rateCompute=.false.)
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMaximum,node,rateCompute=.false.)
          !$omp parallel private(iTime,i,j,sigma1f,integralKernel,excursionSetBarrier_,barrierTable,barrierMidTable) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
          allocate(excursionSetBarrier_,mold=self%excursionSetBarrier_)
          !$omp critical(excursionSetsSolverFarahiMidpointDeepCopy)
          !![
          <deepCopyReset variables="self%excursionSetBarrier_"/>
	  <deepCopy source="self%excursionSetBarrier_" destination="excursionSetBarrier_"/>
	  <deepCopyFinalize variables="excursionSetBarrier_"/>
	  !!]
          !$omp end critical(excursionSetsSolverFarahiMidpointDeepCopy)
          call allocateArray(barrierTable   ,[1+self%varianceTableCount],lowerBounds=[0])
          call allocateArray(barrierMidTable,[1+self%varianceTableCount],lowerBounds=[0])
          !$omp do schedule(dynamic)
          do iTime=1,self%timeTableCount
#ifdef USEMPI
             if (self%coordinatedMPI_ .and. mod(iTime-1,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
             ! Construct the barrier table.
             do i=0,self%varianceTableCount
                barrierTable   (i)=excursionSetBarrier_%barrier(self%varianceTable   (i),self%timeTable(iTime),node,rateCompute=.false.)
                barrierMidTable(i)=excursionSetBarrier_%barrier(     varianceMidTable(i),self%timeTable(iTime),node,rateCompute=.false.)
             end do
             self%firstCrossingProbabilityTable(0,iTime)=0.0d0
             integralKernel                             =+1.0_kind_quad                                                                   &
                  &                                      -erfApproximate(                                                                 &
                  &                                                      +(                                                               &
                  &                                                        +barrierTable   (1)                                            &
                  &                                                        -barrierMidTable(1)                                            &
                  &                                                       )                                                               &
                  &                                                      /sqrt(2.0_kind_quad*(self%varianceTable(1)-varianceMidTable(1))) &
                  &                                                     )
             self%firstCrossingProbabilityTable(1,iTime)=real(                                                              &
                  &                                           +(                                                            &
                  &                                             +1.0_kind_quad                                              &
                  &                                             -erfApproximate(                                            &
                  &                                                             +barrierTable(1)                            &
                  &                                                             /sqrt(2.0_kind_quad*self%varianceTable(1))  &
                  &                                                            )                                            &
                  &                                            )                                                            &
                  &                                           /self%varianceTableStep                                       &
                  &                                           /integralKernel                                             , &
                  &                                           kind=kind_dble                                                &
                  &                                          )
             do i=2,self%varianceTableCount
#ifdef USEMPI
                if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
                   call displayCounter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityLevelWorking)
#ifdef USEMPI
                end if
#endif
                !$omp atomic
                loopCount=loopCount+(i-1)
                sigma1f  =0.0d0
                do j=1,i-1
                   sigma1f=+sigma1f                                                                                  &
                        &  +self%firstCrossingProbabilityTable(j,iTime)                                              &
                        &  *real(                                                                                    &
                        &         1.0_kind_quad                                                                      &
                        &        -erfApproximate(                                                                    &
                        &                          (                                                                 &
                        &                           +barrierTable   (i)                                              &
                        &                           -barrierMidTable(j)                                              &
                        &                          )                                                                 &
                        &                          /sqrt(2.0_kind_quad*(self%varianceTable(i)-varianceMidTable(j)))  &
                        &                         )                                                                , &
                        &        kind=kind_dble                                                                      &
                        &       )
                end do
                integralKernel=+1.0_kind_quad                                                                   &
                     &         -erfApproximate(                                                                 &
                     &                         +(                                                               &
                     &                           +barrierTable   (i)                                            &
                     &                           -barrierMidTable(i)                                            &
                     &                          )                                                               &
                     &                         /sqrt(2.0_kind_quad*(self%varianceTable(i)-varianceMidTable(i))) &
                     &                        )
                if (integralKernel==0.0_kind_quad) then
                   self%firstCrossingProbabilityTable(i,iTime)=0.0d0
                else
                   self%firstCrossingProbabilityTable(i,iTime)=real(                                                                    &
                        &                                           max(                                                                &
                        &                                               +0.0_kind_quad,                                                 &
                        &                                               +(                                                              &
                        &                                                 +(                                                            &
                        &                                                   +1.0_kind_quad                                              &
                        &                                                   -erfApproximate(                                            &
                        &                                                                   +barrierTable(i)                            &
                        &                                                                   /sqrt(2.0_kind_quad*self%varianceTable(i))  &
                        &                                                                  )                                            &
                        &                                                  )                                                            &
                        &                                                 /self%varianceTableStep                                       &
                        &                                                 -1.0_kind_quad                                                &
                        &                                                 *sigma1f                                                      &
                        &                                                )                                                              &
                        &                                               /integralKernel                                                 &
                        &                                              )                                                              , &
                        &                                           kind=kind_dble                                                      &
                        &                                          )
                end if
             end do
             ! Force the probability at maximum variance to zero.
             self%firstCrossingProbabilityTable(self%varianceTableCount,iTime)=0.0d0
          end do
          !$omp end do
          !![
          <objectDestructor name="excursionSetBarrier_"/>
	  !!]
          call deallocateArray(barrierTable   )
          call deallocateArray(barrierMidTable)
          !$omp end parallel
          ! Update the variance table to reflect the variances at the midpoints. Note that the first crossing probability is computed
          ! at the mid-points. The last element of the variance table is unchanged to ensure that its value equals
          ! varianceMaximum. This will not affect the result becasue the probability at maximum variance is set to zero anyway.
          self%varianceTable(1:self%varianceTableCount-1)=varianceMidTable(1:self%varianceTableCount-1)
          call deallocateArray(varianceMidTable)
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayCounterClear(verbosityLevelWorking)
             call displayUnindent("done",verbosityLevelWorking)
#ifdef USEMPI
          end if
          if (self%coordinatedMPI_) then
             call mpiBarrier()
             self%firstCrossingProbabilityTable=mpiSelf%sum(self%firstCrossingProbabilityTable)
          end if
#endif
          ! Build the interpolators.
          if (allocated(self%interpolatorVariance)) deallocate(self%interpolatorVariance)
          if (allocated(self%interpolatorTime    )) deallocate(self%interpolatorTime    )
          allocate(self%interpolatorVariance)
          allocate(self%interpolatorTime    )
          self%interpolatorVariance=interpolator(self%varianceTable)
          self%interpolatorTime    =interpolator(self%timeTable    )
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
          farahiMidpointProbability=+farahiMidpointProbability                                             &
               &                    +hTime                             (                            jTime) &
               &                    *hVariance                         (            jVariance            ) &
               &                    *self%firstCrossingProbabilityTable(iVariance-1+jVariance,iTime+jTime)
       end do
    end do
    return
  end function farahiMidpointProbability

  subroutine farahiMidpointRateTabulate(self,varianceProgenitor,time,node)
    !!{
    Tabulate the excursion set crossing rate.
    !!}
    use :: Display          , only : displayCounter , displayCounterClear  , displayIndent       , displayMessage, &
          &                          displayUnindent, verbosityLevelWorking
    use :: Error_Functions  , only : erfApproximate
    use :: File_Utilities   , only : File_Lock      , File_Unlock          , lockDescriptor
    use :: Kind_Numbers     , only : kind_dble      , kind_quad
    use :: MPI_Utilities    , only : mpiBarrier     , mpiSelf
    use :: Memory_Management, only : allocateArray  , deallocateArray
    use :: Numerical_Ranges , only : Make_Range     , rangeTypeLinear      , rangeTypeLogarithmic
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpoint), intent(inout)                   :: self
    double precision                                         , intent(in   )                   :: time                             , varianceProgenitor
    type            (treeNode                               ), intent(inout)                   :: node
    double precision                                         , parameter                       :: varianceMinimumDefault    =1.0d-2
    double precision                                         , parameter                       :: varianceTolerance         =1.0d-6
    double precision                                         , parameter                       :: massLarge                 =1.0d16
    real            (kind=kind_quad                         ), allocatable  , dimension(:    ) :: firstCrossingTableRateQuad       , varianceTableRateBaseQuad, &
         &                                                                                        varianceTableRateQuad            , varianceMidTableRateQuad , &
         &                                                                                        barrierTableRateQuad             , barrierMidTableRateQuad
    double precision                                         , allocatable  , dimension(:,:  ) :: nonCrossingTableRate
    double precision                                         , allocatable  , dimension(:,:,:) :: firstCrossingTableRate
    double precision                                                                           :: barrierRateTest
    class           (excursionSetBarrierClass               ), pointer                         :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass          ), pointer                         :: cosmologicalMassVariance_
#ifdef USEMPI
    integer                                                                                    :: taskCount
#endif
    logical                                                                                    :: makeTable
    integer         (c_size_t                               )                                  :: loopCount                        , loopCountTotal
    integer                                                                                    :: i                                , iTime                    , &
         &                                                                                        iVariance                        , j                        , &
         &                                                                                        countNewLower                    , countNewUpper            , &
         &                                                                                        timeTableCountNew
    double precision                                                                           :: timeProgenitor                   , varianceMinimumRate      , &
         &                                                                                        massProgenitor                   , timeMinimumRate          , &
         &                                                                                        timeMaximumRate
    character       (len=9                                  )                                  :: label
    type            (varying_string                         )                                  :: message
    type            (lockDescriptor                         )                                  :: fileLock
    real            (kind=kind_quad                         )                                  :: crossingFraction                 , effectiveBarrierInitial  , &
         &                                                                                        sigma1f                          , varianceTableStepRate    , &
         &                                                                                        barrier                          , integralKernelRate       , &
         &                                                                                        growthFactorEffective            , erfArgumentNumerator     , &
         &                                                                                        erfArgumentDenominator           , erfValue
    logical                                                                                    :: varianceMaximumChanged

    ! Determine if we need to make the table.
    !
    !! Read tables from file if possible. We make two passes through the logic that determines if the table needs to be remade. On
    !! the first pass if the table does need to be remade we attempt to read it from file. If the file is read, then we re-check if
    !! the table needs to be remade.
    do i=1,2
       makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
       if (i == 1 .and. self%useFile .and. makeTable) then
          call self%fileNameInitialize()
          call File_Lock(char(self%fileName),fileLock)
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
       !$omp critical(farahiMidpointRateTabulate)
       ! Attempt to read the file again now that we are within the critical section. If another thread made the file while we were waiting we may be able to skip building the table.
       if (self%useFile) then
          call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       end if
       makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
       if (makeTable) then
          ! Construct or expand the range of times to tabulate.
          countNewLower=0
          countNewUpper=0
          if (self%tableInitializedRate) then
             varianceMaximumChanged=varianceProgenitor > self%varianceMaximumRate
             timeMinimumRate=min(time/10.0d0**(2.0d0/dble(self%timeNumberPerDecade)),self%timeMinimumRate)
             timeMaximumRate=max(time*10.0d0**(2.0d0/dble(self%timeNumberPerDecade)),self%timeMaximumRate)
             ! Determine how many points the table must be extended by in each direction to span the new required range.
             if (self%timeMinimumRate > timeMinimumRate) countNewLower=int(+log10(self%timeMinimumRate/timeMinimumRate)*dble(self%timeNumberPerDecade)+1.0d0)
             if (self%timeMaximumRate < timeMaximumRate) countNewUpper=int(-log10(self%timeMaximumRate/timeMaximumRate)*dble(self%timeNumberPerDecade)+1.0d0)
             self%timeTableCountRate=self%timeTableCountRate+countNewLower+countNewUpper
             ! Adjust the limits of the table by an integer number of steps.
             self%timeMinimumRate=self%timeMinimumRate/10.0d0**(dble(countNewLower)/dble(self%timeNumberPerDecade))
             self%timeMaximumRate=self%timeMaximumRate*10.0d0**(dble(countNewUpper)/dble(self%timeNumberPerDecade))
          else
             varianceMaximumChanged =.true.
             self%timeMinimumRate   =time/10.0d0**(2.0d0/dble(self%timeNumberPerDecade))
             self%timeMaximumRate   =time*10.0d0**(2.0d0/dble(self%timeNumberPerDecade))
             self%timeTableCountRate=max(int(log10(self%timeMaximumRate/self%timeMinimumRate)*dble(self%timeNumberPerDecade))+2,2)
             ! Ensure the maximum of the table is precisely an integer number of steps above the minimum.
             self%timeMaximumRate   =self%timeMinimumRate*10.0d0**(dble(self%timeTableCountRate-1)/dble(self%timeNumberPerDecade))
          end if
          ! Set the default minimum variance.
          varianceMinimumRate       =varianceMinimumDefault
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
          growthFactorEffective          =+cosmologicalMassVariance_%rootVariance(massLarge,self%timeMaximumRate                                ) &
               &                          /cosmologicalMassVariance_%rootVariance(massLarge,self%timeMaximumRate*(1.0d0-self%timeStepFractional))
          varianceMinimumRate            =min(                                                                                                                      &
               &                              +varianceMinimumRate                                                                                                , &
               &                              +1.0d-2                                                                                                               &
               &                              *(                                                                                                                    &
               &                                +excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate*(1.0d0-self%timeStepFractional),node,rateCompute=.true.)  &
               &                                *dble(growthFactorEffective)                                                                                        &
               &                                -excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate                                ,node,rateCompute=.true.)  &
               &                               )**2                                                                                                                 &
               &                             )
          !![
          <objectDestructor name="excursionSetBarrier_"     />
          <objectDestructor name="cosmologicalMassVariance_"/>
          !!]
          self%varianceMaximumRate       =max(self%varianceMaximumRate,varianceProgenitor)
          self%varianceTableCountRate    =int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecade))+1
          self%varianceTableCountRateBase=int(self%varianceMaximumRate*dble(self%varianceNumberPerUnit))
          ! Store copies of the current tables if these will be used later.
          if (.not.varianceMaximumChanged) then
             call move_alloc(self%firstCrossingTableRate,firstCrossingTableRate)
             call move_alloc(self%nonCrossingTableRate  ,nonCrossingTableRate  )
          end if
          if (allocated(self%varianceTableRate     )) call deallocateArray(self%varianceTableRate     )
          if (allocated(self%varianceTableRateBase )) call deallocateArray(self%varianceTableRateBase )
          if (allocated(self%timeTableRate         )) call deallocateArray(self%timeTableRate         )
          if (allocated(self%firstCrossingTableRate)) call deallocateArray(self%firstCrossingTableRate)
          if (allocated(self%nonCrossingTableRate  )) call deallocateArray(self%nonCrossingTableRate  )
          call allocateArray(self%varianceTableRate     ,[1+self%varianceTableCountRate                                                          ],lowerBounds=[0    ])
          call allocateArray(self%varianceTableRateBase ,[                              1+self%varianceTableCountRateBase                        ],lowerBounds=[0    ])
          call allocateArray(self%timeTableRate         ,[                                                                self%timeTableCountRate]                    )
          call allocateArray(self%firstCrossingTableRate,[1+self%varianceTableCountRate,1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[0,0,1])
          call allocateArray(self%nonCrossingTableRate  ,[                              1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[  0,1])
          ! If only times have changed then pre-populate the tables with results previously computed.
          if (.not.varianceMaximumChanged) then
             self%firstCrossingTableRate(:,:,countNewLower+1:countNewLower+size(firstCrossingTableRate,dim=3))=firstCrossingTableRate
             self%  nonCrossingTableRate(  :,countNewLower+1:countNewLower+size(  nonCrossingTableRate,dim=2))=  nonCrossingTableRate
             deallocate(firstCrossingTableRate)
             deallocate(  nonCrossingTableRate)
          end if
          ! For the variance table, the zeroth point is always zero, higher points are distributed uniformly in variance.
          self%varianceTableRate    (0                                )=0.0d0
          self%varianceTableRate    (1:self%varianceTableCountRate    )=self%varianceRange(varianceMinimumRate,self%varianceMaximumRate,self%varianceTableCountRate      ,exponent =1.0d0          )
          self%varianceTableRateBase(0:self%varianceTableCountRateBase)=Make_Range        (0.0d0              ,self%varianceMaximumRate,self%varianceTableCountRateBase+1,rangeType=rangeTypeLinear)
          ! Allocate temporary arrays used in quad-precision solver for barrier crossing rates.
          allocate(varianceTableRateQuad     (0:self%varianceTableCountRate    ))
          varianceTableRateQuad    =self%varianceTableRate
          allocate(varianceTableRateBaseQuad (0:self%varianceTableCountRateBase))
          varianceTableRateBaseQuad=self%varianceTableRateBase
          ! The time table is logarithmically distributed in time.
          self%timeTableRate=Make_Range(self%timeMinimumRate,self%timeMaximumRate,self%timeTableCountRate,rangeType=rangeTypeLogarithmic)
          ! Compute the variance at the mid-points.
          call allocateArray(varianceMidTableRateQuad,[1+self%varianceTableCountRate],lowerBounds=[0])
          varianceMidTableRateQuad(0)=0.0_kind_quad
          forall(i=1:self%varianceTableCountRate)
             varianceMidTableRateQuad(i)=(varianceTableRateQuad(i-1)+varianceTableRateQuad(i))/2.0_kind_quad
          end forall
          ! Loop through the table and solve for the first crossing distribution.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayIndent("solving for excursion set barrier crossing rates",verbosityLevelWorking)
             message="    time: "
             write (label,'(f6.3)') self%timeMinimumRate
             message=message//label//" to "
             write (label,'(f6.3)') self%timeMaximumRate
             message=message//label
             call displayMessage(message,verbosityLevelWorking)
             message="variance: "
             write (label,'(f9.3)') self%varianceMaximumRate
             message=message//label
             call displayMessage(message,verbosityLevelWorking)
#ifdef USEMPI
          end if
#endif
          timeTableCountNew=self%timeTableCountRate
          if (.not.varianceMaximumChanged) timeTableCountNew=countNewLower+countNewUpper
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
             loopCountTotal=int(timeTableCountNew,kind=c_size_t)*int(self%varianceTableCountRateBase+1,kind=c_size_t)/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t
          else
#endif
             loopCountTotal=int(timeTableCountNew,kind=c_size_t)*int(self%varianceTableCountRateBase+1,kind=c_size_t)
#ifdef USEMPI
          end if
#endif
          loopCount=0
#ifdef USEMPI
          if (self%coordinatedMPI_) then
             self%firstCrossingTableRate=0.0d0
             self%nonCrossingTableRate  =0.0d0
          end if
          taskCount=-1
#endif
          ! Make a call to the barrier function at maximum variance for the minimum and maximum times so that the barrier function
          ! is initialized and covers the whole range we are intereseted in.
          barrierRateTest=self%excursionSetBarrier_%barrier(self%varianceMaximumRate,self%timeMinimumRate*(1.0d0-self%timeStepFractional),node,rateCompute=.true.)
          barrierRateTest=self%excursionSetBarrier_%barrier(self%varianceMaximumRate,self%timeMaximumRate                                ,node,rateCompute=.true.)
          !$omp parallel private(iTime,timeProgenitor,iVariance,varianceTableStepRate,i,j,sigma1f,integralKernelRate,crossingFraction,barrier,effectiveBarrierInitial,firstCrossingTableRateQuad,excursionSetBarrier_,cosmologicalMassVariance_,barrierTableRateQuad,barrierMidTableRateQuad,massProgenitor,growthFactorEffective) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
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
          call allocateArray(barrierTableRateQuad   ,[self%varianceTableCountRate])
          call allocateArray(barrierMidTableRateQuad,[self%varianceTableCountRate])
          do iTime=1,self%timeTableCountRate
             ! Skip if this time was already computed.
             if (.not.varianceMaximumChanged.and.(iTime > countNewLower .and. self%timeTableCountRate+1-iTime > countNewUpper)) cycle
             ! Allocate workspace table.
             if (.not.allocated(firstCrossingTableRateQuad)) allocate(firstCrossingTableRateQuad(0:self%varianceTableCountRate))
             ! Compute a suitable progenitor time.
             timeProgenitor=self%timeTableRate(iTime)*(1.0d0-self%timeStepFractional)
             ! Loop through the starting variances.
             !$omp do schedule(dynamic)
             do iVariance=0,self%varianceTableCountRateBase
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
                ! Construct the barrier table.
                do i=1,self%varianceTableCountRate
                   massProgenitor            =+cosmologicalMassVariance_%mass        (real(sqrt(+varianceTableRateQuad   (i)+varianceTableRateBaseQuad(iVariance)),kind=8),self%timeTableRate(iTime)                        )
                   growthFactorEffective     =+cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,self%timeTableRate(iTime)                        ) &
                        &                     /cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,     timeProgenitor                              )
                   barrierTableRateQuad   (i)=+excursionSetBarrier_     %barrier     (real(     +varianceTableRateQuad   (i)+varianceTableRateBaseQuad(iVariance) ,kind=8),     timeProgenitor      ,node,rateCompute=.true.) &
                        &                     *growthFactorEffective
                   massProgenitor            =+cosmologicalMassVariance_%mass        (real(sqrt(+varianceMidTableRateQuad(i)+varianceTableRateBaseQuad(iVariance)),kind=8),self%timeTableRate(iTime)                        )
                   growthFactorEffective     =+cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,self%timeTableRate(iTime)                        ) &
                        &                     /cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,     timeProgenitor                              )
                   barrierMidTableRateQuad(i)=+excursionSetBarrier_     %barrier     (real(     +varianceMidTableRateQuad(i)+varianceTableRateBaseQuad(iVariance),kind=8 ),     timeProgenitor      ,node,rateCompute=.true.) &
                        &                     *growthFactorEffective
                end do
                ! For zero variance, the rate is initialized to zero.
                firstCrossingTableRateQuad(0)=0.0_kind_quad
                ! Compute the step in variance across this first grid cell.
                varianceTableStepRate=varianceTableRateQuad(1)-varianceTableRateQuad(0)
                ! Compute the barrier for the descendent.
                barrier=real(excursionSetBarrier_%barrier(real(varianceTableRateBaseQuad(iVariance),kind=8),self%timeTableRate(iTime),node,rateCompute=.true.),kind=kind_quad)
                ! Compute the first crossing distribution at the first grid cell.
                if (varianceTableRateQuad(1)+varianceTableRateBaseQuad(iVariance) > self%varianceMaximumRate) then
                   firstCrossingTableRateQuad(1)= 0.0_kind_quad
                else
                   integralKernelRate           =+1.0_kind_quad                                                                              &
                        &                        -erfApproximate(                                                                            &
                        &                                        +(                                                                          &
                        &                                          +(barrierTableRateQuad   (1)-barrier)                                     &
                        &                                          -(barrierMidTableRateQuad(1)-barrier)                                     &
                        &                                         )                                                                          &
                        &                                        /sqrt(2.0_kind_quad*(varianceTableRateQuad(1)-varianceMidTableRateQuad(1))) &
                        &                                       )
                   ! If the integral kernel is zero (to machine precision) then simply assume no crossing rate.
                   if (integralKernelRate <= 0.0d0) then
                      firstCrossingTableRateQuad=0.0d0
                      cycle
                   end if
                   firstCrossingTableRateQuad(1)=+(                                                              &
                        &                          +1.0_kind_quad                                                &
                        &                          -erfApproximate(                                              &
                        &                                          +(barrierTableRateQuad(1)-barrier)            &
                        &                                          /sqrt(2.0_kind_quad*varianceTableRateQuad(1)) &
                        &                                         )                                              &
                        &                         )                                                              &
                        &                         /varianceTableStepRate                                         &
                        &                         /integralKernelRate
                end if
                do i=2,self%varianceTableCountRate
                   if (varianceTableRateQuad(i)+varianceTableRateBaseQuad(iVariance) > self%varianceMaximumRate) then
                      firstCrossingTableRateQuad(i)=0.0_kind_quad
                   else
                      effectiveBarrierInitial=+barrierTableRateQuad(i)-barrier
                      integralKernelRate=+1.0_kind_quad                                                                              &
                           &             -erfApproximate(                                                                            &
                           &                             +(                                                                          &
                           &                               +effectiveBarrierInitial                                                  &
                           &                               -barrierMidTableRateQuad(i)                                               &
                           &                               +barrier                                                                  &
                           &                              )                                                                          &
                           &                             /sqrt(2.0_kind_quad*(varianceTableRateQuad(i)-varianceMidTableRateQuad(i))) &
                           &                            )
                      if (integralKernelRate==0.0_kind_quad) then
                         firstCrossingTableRateQuad(i)=0.0_kind_quad
                      else
                         sigma1f=+0.0_kind_quad
                         do j=1,i-1
                            erfArgumentNumerator=+effectiveBarrierInitial    &
                                 &               -barrierMidTableRateQuad(j) &
                                 &               +barrier
                            if (erfArgumentNumerator == 0.0d0) then
                               erfValue=0.0_kind_quad
                            else
                               erfArgumentDenominator=sqrt(2.0_kind_quad*(varianceTableRateQuad(i)-varianceMidTableRateQuad(j)))
                               erfValue=erfApproximate(erfArgumentNumerator/erfArgumentDenominator)
                            end if
                            if (erfValue < 1.0_kind_quad) then
                               varianceTableStepRate=varianceTableRateQuad(j)-varianceTableRateQuad(j-1)
                               sigma1f=+sigma1f                       &
                                    &  +firstCrossingTableRateQuad(j) &
                                    &  *varianceTableStepRate         &
                                    &  *(                             &
                                    &    +1.0_kind_quad               &
                                    &    -erfValue                    &
                                    &   )
                            end if
                         end do
                         varianceTableStepRate=varianceTableRateQuad(i)-varianceTableRateQuad(i-1)
                         firstCrossingTableRateQuad(i)=max(                                                                &
                              &                            +0.0_kind_quad,                                                 &
                              &                            +(                                                              &
                              &                              +1.0_kind_quad                                                &
                              &                              -erfApproximate(                                              &
                              &                                              +effectiveBarrierInitial                      &
                              &                                              /sqrt(2.0_kind_quad*varianceTableRateQuad(i)) &
                              &                                             )                                              &
                              &                              -sigma1f                                                      &
                              &                             )                                                              &
                              &                            /varianceTableStepRate                                          &
                              &                            /integralKernelRate                                             &
                              &                           )
                         ! Remove unphysical values. Force the crossing rate at points close to maximum variance to zero if its value is one order of magnitude
                         ! larger than the value at previous point.
                         if     (                                                                                                                               &
                              &   varianceTableRateQuad(i)+varianceTableRateBaseQuad(iVariance) > self%varianceMaximumRate-10.0_kind_quad*varianceTableStepRate &
                              &  .and.                                                                                                                          &
                              &   firstCrossingTableRateQuad(i) > 10.0_kind_quad*firstCrossingTableRateQuad(i-1)                                                &
                              & )                                                                                                                               &
                              & firstCrossingTableRateQuad(i)=0.0_kind_quad
                      end if
                   end if
                end do
                ! Compute the fraction of trajectories which never cross the barrier.
                crossingFraction=0.0_kind_quad
                do j=1,self%varianceTableCountRate
                   varianceTableStepRate=varianceTableRateQuad(j)-varianceTableRateQuad(j-1)
                   crossingFraction     =+crossingFraction              &
                        &                +firstCrossingTableRateQuad(j) &
                        &                *varianceTableStepRate
                end do
                ! Compute the rate for trajectories which never cross the barrier.
                self%nonCrossingTableRate(iVariance,iTime)=real(                                   &
                     &                                          +(1.0_kind_quad-crossingFraction)  &
                     &                                          /self%timeTableRate(iTime)         &
                     &                                          /self%timeStepFractional         , &
                     &                                          kind=kind_dble                     &
                     &                                         )
                ! Store the compute crossing rate in our table.
                self%firstCrossingTableRate(:,iVariance,iTime)=real(firstCrossingTableRateQuad,kind=kind_dble)
             end do
             !$omp end do
             ! Divide through by the time step to get the rate of barrier crossing.
             !$omp single
             self%firstCrossingTableRate(:,:,iTime)=+self%firstCrossingTableRate(:,:,iTime) &
                  &                                 /self%timeTableRate         (    iTime) &
                  &                                 /self%timeStepFractional
             !$omp end single
          end do
          !![
          <objectDestructor name="excursionSetBarrier_"     />
          <objectDestructor name="cosmologicalMassVariance_"/>
          !!]
          call deallocateArray(barrierTableRateQuad   )
          call deallocateArray(barrierMidTableRateQuad)
          !$omp end parallel
          ! Update the variance table to reflect the variances at the midpoints. Note that the first crossing probability is computed
          ! at the mid-points. The last element of the variance table is unchanged to ensure that its value equals
          ! varianceMaximum. This will not affect the result becasue the crossing rate at maximum variance is set to zero anyway.
          self%varianceTableRate(1:self%varianceTableCountRate-1)=real(varianceMidTableRateQuad(1:self%varianceTableCountRate-1),kind=kind_dble)
          call deallocateArray(varianceMidTableRateQuad)
          ! Deallocate work arrays.
          deallocate(varianceTableRateBaseQuad )
          deallocate(varianceTableRateQuad     )
          if (allocated(firstCrossingTableRateQuad)) deallocate(firstCrossingTableRateQuad)
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayCounterClear(       verbosityLevelWorking)
             call displayUnindent     ("done",verbosityLevelWorking)
#ifdef USEMPI
          end if
          if (self%coordinatedMPI_) then
             call mpiBarrier()
             self%firstCrossingTableRate=mpiSelf%sum(self%firstCrossingTableRate)
             self%  nonCrossingTableRate=mpiSelf%sum(self%  nonCrossingTableRate)
          end if
#endif
          ! Build the interpolators.
          if (allocated(self%interpolatorVarianceRate    )) deallocate(self%interpolatorVarianceRate    )
          if (allocated(self%interpolatorVarianceRateBase)) deallocate(self%interpolatorVarianceRateBase)
          if (allocated(self%interpolatorTimeRate        )) deallocate(self%interpolatorTimeRate        )
          allocate(self%interpolatorVarianceRate    )
          allocate(self%interpolatorVarianceRateBase)
          allocate(self%interpolatorTimeRate        )
          self%interpolatorVarianceRate    =interpolator(self%varianceTableRate    )
          self%interpolatorVarianceRateBase=interpolator(self%varianceTableRateBase)
          self%interpolatorTimeRate        =interpolator(self%timeTableRate        )
          ! Set previous variance and time to unphysical values to force recompute of interpolation factors on next call.
          self%varianceRatePrevious=-1.0d0
          self%timeRatePrevious    =-1.0d0
          ! Record that the table is now built.
          self%tableInitializedRate=.true.
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
