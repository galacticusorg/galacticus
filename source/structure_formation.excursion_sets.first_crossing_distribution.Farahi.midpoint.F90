!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements a excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method to perform the integrations \citep{du_substructure_2017}.

  !# <excursionSetFirstCrossing name="excursionSetFirstCrossingFarahiMidpoint">
  !#  <description>An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method to perform the integrations \citep{du_substructure_2017}.</description>
  !# </excursionSetFirstCrossing>
  type, extends(excursionSetFirstCrossingFarahi) :: excursionSetFirstCrossingFarahiMidpoint
     !% An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method to perform the integrations \citep{du_substructure_2017}.
     private
   contains
  end type excursionSetFirstCrossingFarahiMidpoint
  
  interface excursionSetFirstCrossingFarahiMidpoint
     !% Constructors for the Farahi-midpoint excursion set barrier class.
     module procedure farahiMidpointConstructorParameters
     module procedure farahiMidpointConstructorInternal
  end interface excursionSetFirstCrossingFarahiMidpoint

contains

  function farahiMidpointConstructorParameters(parameters) result(self)
    !% Constructor for the Farahi-midpoint excursion set class first crossing class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(excursionSetFirstCrossingFarahiMidpoint)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self%excursionSetFirstCrossingFarahi=excursionSetFirstCrossingFarahi(parameters)
    return
  end function farahiMidpointConstructorParameters

  function farahiMidpointConstructorInternal(timeStepFractional,fileName,cosmologyFunctions_,excursionSetBarrier_) result(self)
    !% Internal constructor for the Farahi-midpoint excursion set class first crossing class.
    use Input_Parameters
    implicit none
    type            (excursionSetFirstCrossingFarahiMidpoint)                        :: self
    double precision                                         , intent(in   )         :: timeStepFractional
    type            (varying_string                         ), intent(in   )         :: fileName
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetBarrierClass               ), intent(in   ), target :: excursionSetBarrier_

    self%excursionSetFirstCrossingFarahi=excursionSetFirstCrossingFarahi(timeStepFractional,fileName,cosmologyFunctions_,excursionSetBarrier_)
    return
  end function farahiMidpointConstructorInternal

  double precision function farahiMidpointProbability(self,variance,time,node)
    !% Return the excursion set barrier at the given variance and time.
    use Numerical_Ranges
    use Numerical_Interpolation
    use Memory_Management
    use Galacticus_Display
    use Kind_Numbers
    use Error_Functions
    use MPI_Utilities
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpoint), intent(inout)                 :: self
    double precision                                         , intent(in   )                 :: variance                     , time
    type            (treeNode                               ), intent(inout)                 :: node
    double precision                                                        , dimension(0:1) :: hTime                        , hVariance
    double precision                                         , parameter                     :: varianceTableTolerance=1.0d-6
    double precision                                         , allocatable  , dimension( : ) :: varianceMidTable
    class           (excursionSetBarrierClass               ), allocatable                   :: excursionSetBarrier_
    logical                                                                                  :: makeTable
    integer         (c_size_t                               )                                :: iTime                        , iVariance
    integer                                                                                  :: i                            , j             , &
         &                                                                                      jTime                        , jVariance     , &
         &                                                                                      loopCount                    , loopCountTotal
    double precision                                                                         :: sigma1f
    real            (kind=kind_quad                         )                                :: integralKernel                        
    character       (len =6                                 )                                :: label
    type            (varying_string                         )                                :: message
    logical                                                                                  :: locked

    ! Read tables from file if possible.
    locked=.false.
    if (self%useFile.and..not.self%tableInitialized) then
       call File_Lock(char(self%fileName),farahiFileLock)
       locked=.true.
       call self%fileRead()
    end if
    ! Construct the table if necessary.
    makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+varianceTableTolerance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
#ifdef USEMPI
    if (self%coordinatedMPI_) then
       if (locked) then
          call File_Unlock(farahiFileLock)
          locked=.false.
       end if
       call mpiBarrier()
    end if
#endif
    if (makeTable) then
#ifdef USEMPI
       ! If coordinating under MPI then only the rank-0 process locks the file.
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          ! Construct the table of variance on which we will solve for the first crossing distribution.
          if (self%useFile.and..not.locked) then
             call File_Lock(char(self%fileName),farahiFileLock)
             locked=.true.
          end if
#ifdef USEMPI
       end if
#endif
       if (allocated(self%varianceTable                )) call deallocateArray(self%varianceTable                )
       if (allocated(self%timeTable                    )) call deallocateArray(self%timeTable                    )
       if (allocated(self%firstCrossingProbabilityTable)) call deallocateArray(self%firstCrossingProbabilityTable)
       self%varianceMaximum   =max(self%varianceMaximum,variance)
       self%varianceTableCount=int(self%varianceMaximum*dble(farahiVarianceNumberPerUnityProbability))
       if (self%tableInitialized) then
          self%timeMinimum=min(self%timeMinimum,0.5d0*time)
          self%timeMaximum=max(self%timeMaximum,2.0d0*time)
       else
          self%timeMinimum=                     0.5d0*time
          self%timeMaximum=                     2.0d0*time
       end if
       self%timeTableCount=max(2,int(log10(self%timeMaximum/self%timeMinimum)*dble(farahiTimeNumberPerDecade))+1)
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
          call Galacticus_Display_Indent("solving for excursion set barrier crossing probabilities",verbosityWorking)
          message="    time: "
          write (label,'(f6.3)') self%timeMinimum
          message=message//label//" to "
          write (label,'(f6.3)') self%timeMaximum
          message=message//label
          call Galacticus_Display_Message(message,verbosityWorking)
          message="variance: "
          write (label,'(f6.3)') self%varianceMaximum
          message=message//label
          call Galacticus_Display_Message(message,verbosityWorking)
#ifdef USEMPI
       end if
#endif
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
          loopCountTotal=(self%timeTableCount/mpiSelf%count()+1)*((self%varianceTableCount-1)*self%varianceTableCount)/2
       else
#endif
          loopCountTotal= self%timeTableCount                   *((self%varianceTableCount-1)*self%varianceTableCount)/2
#ifdef USEMPI
       end if
#endif
       loopCount     =0
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) self%firstCrossingProbabilityTable=0.0d0
#endif
       !$omp parallel private(iTime,i,j,sigma1f,integralKernel,excursionSetBarrier_) if (.not.mpiSelf%isActive())
       allocate(excursionSetBarrier_,mold=self%excursionSetBarrier_)
       !# <deepCopy source="self%excursionSetBarrier_" destination="excursionSetBarrier_"/>
       !$omp do schedule(dynamic)
       do iTime=1,self%timeTableCount
#ifdef USEMPI
          if (self%coordinatedMPI_ .and. mod(iTime-1,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
          self%firstCrossingProbabilityTable(0,iTime)=0.0d0
          integralKernel                             =+1.0_kind_quad                                                                                                           &
               &                                      -erfApproximate(                                                                                                         &
               &                                                      +(                                                                                                       &
               &                                                        +excursionSetBarrier_%barrier(self%varianceTable   (1),self%timeTable(iTime),node,rateCompute=.false.) &
               &                                                        -excursionSetBarrier_%barrier(     varianceMidTable(1),self%timeTable(iTime),node,rateCompute=.false.) &
               &                                                       )                                                                                                       &
               &                                                      /sqrt(2.0_kind_quad*(self%varianceTable(1)-varianceMidTable(1)))                                         &
               &                                                     )
          self%firstCrossingProbabilityTable(1,iTime)=real(                                                                                                                       &
               &                                           +(                                                                                                                     &
               &                                             +1.0_kind_quad                                                                                                       &
               &                                             -erfApproximate(                                                                                                     &
               &                                                             +excursionSetBarrier_%barrier(self%varianceTable(1),self%timeTable(iTime),node,rateCompute=.false.)  &
               &                                                             /sqrt(2.0_kind_quad*self%varianceTable(1))                                                           &
               &                                                            )                                                                                                     &
               &                                            )                                                                                                                     &
               &                                           /self%varianceTableStep                                                                                                &
               &                                           /integralKernel                                                                                                      , &
               &                                           kind=kind_dble                                                                                                         &
               &                                          )
          do i=2,self%varianceTableCount
#ifdef USEMPI
             if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
                call Galacticus_Display_Counter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityWorking)
#ifdef USEMPI
             end if
#endif
             !$omp atomic
             loopCount=loopCount+(i-1)
             sigma1f  =0.0d0
             do j=1,i-1
                sigma1f=+sigma1f                                                                                                                         &
                     &  +self%firstCrossingProbabilityTable(j,iTime)                                                                                     &
                     &  *real(                                                                                                                           &
                     &         1.0_kind_quad                                                                                                             &
                     &        -erfApproximate(                                                                                                           &
                     &                          (                                                                                                        &
                     &                            excursionSetBarrier_%barrier(self%varianceTable   (i),self%timeTable(iTime),node,rateCompute=.false.)  &
                     &                           -excursionSetBarrier_%barrier(     varianceMidTable(j),self%timeTable(iTime),node,rateCompute=.false.)  &
                     &                          )                                                                                                        &
                     &                          /sqrt(2.0_kind_quad*(self%varianceTable(i)-varianceMidTable(j)))                                         &
                     &                         )                                                                                                       , &
                     &        kind=kind_dble                                                                                                             &
                     &       )
             end do
             integralKernel=+1.0_kind_quad                                                                                                           &
                  &         -erfApproximate(                                                                                                         &
                  &                         +(                                                                                                       &
                  &                           +excursionSetBarrier_%barrier(self%varianceTable   (i),self%timeTable(iTime),node,rateCompute=.false.) &
                  &                           -excursionSetBarrier_%barrier(     varianceMidTable(i),self%timeTable(iTime),node,rateCompute=.false.) &
                  &                          )                                                                                                       &
                  &                         /sqrt(2.0_kind_quad*(self%varianceTable(i)-varianceMidTable(i)))                                         &
                  &                        )
             if (integralKernel==0.0_kind_quad) then
               self%firstCrossingProbabilityTable(i,iTime)=0.0d0
             else
               self%firstCrossingProbabilityTable(i,iTime)=real(                                                                                                                             &
                    &                                           max(                                                                                                                         &
                    &                                               +0.0_kind_quad,                                                                                                          &
                    &                                               +(                                                                                                                       &
                    &                                                 +(                                                                                                                     &
                    &                                                   +1.0_kind_quad                                                                                                       &
                    &                                                   -erfApproximate(                                                                                                     &
                    &                                                                   +excursionSetBarrier_%barrier(self%varianceTable(i),self%timeTable(iTime),node,rateCompute=.false.)  &
                    &                                                                   /sqrt(2.0_kind_quad*self%varianceTable(i))                                                           &
                    &                                                                  )                                                                                                     &
                    &                                                  )                                                                                                                     &
                    &                                                 /self%varianceTableStep                                                                                                &
                    &                                                 -1.0_kind_quad                                                                                                         &
                    &                                                 *sigma1f                                                                                                               &
                    &                                                )                                                                                                                       &
                    &                                               /integralKernel                                                                                                          &
                    &                                              )                                                                                                                       , &
                    &                                           kind=kind_dble                                                                                                               &
                    &                                          )
             end if
          end do
          ! Force the probability at maximum variance to zero.
          self%firstCrossingProbabilityTable(self%varianceTableCount,iTime)=0.0d0
       end do
       !$omp end do
       !$omp end parallel       
       ! Update the variance table to reflect the variances at the midpoints. Note that the first crossing probability is computed
       ! at the mid-points. The last element of the variance table is unchanged to ensure that its value equals
       ! varianceMaximum. This will not affect the result becasue the probability at maximum variance is set to zero anyway.
       self%varianceTable(1:self%varianceTableCount-1)=varianceMidTable(1:self%varianceTableCount-1)
       call deallocateArray(varianceMidTable)
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          call Galacticus_Display_Counter_Clear(verbosityWorking)
          call Galacticus_Display_Unindent("done",verbosityWorking)
#ifdef USEMPI
       end if
       if (self%coordinatedMPI_) then
          call mpiBarrier()
          self%firstCrossingProbabilityTable=mpiSelf%sum(self%firstCrossingProbabilityTable)
       end if
#endif
       ! Reset the interpolators.
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorVariance,reset=self%interpolationResetVariance)
       call Interpolate_Done(interpolationAccelerator=self%interpolationAcceleratorTime    ,reset=self%interpolationResetTime    )
       self%interpolationResetVariance=.true.
       self%interpolationResetTime    =.true.
       ! Record that the table is now built.
       self%tableInitialized=.true.
       ! Write the table to file if possible.
#ifdef USEMPI
       if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
          if (self%useFile) call self%fileWrite()       
#ifdef USEMPI
       end if
#endif
    end if
    if (locked) call File_Unlock(farahiFileLock)
    ! Get interpolation in time.
    iTime    =Interpolate_Locate                 (self%timeTable    ,self%interpolationAcceleratorTime    ,time    ,reset=self%interpolationResetTime    )
    hTime    =Interpolate_Linear_Generate_Factors(self%timeTable    ,iTime    ,time    )
    ! Get interpolation in variance.
    iVariance=Interpolate_Locate                 (self%varianceTable,self%interpolationAcceleratorVariance,variance,reset=self%interpolationResetVariance)
    hVariance=Interpolate_Linear_Generate_Factors(self%varianceTable,iVariance,variance)
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
