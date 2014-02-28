!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines methods for a Gaussian regression likelihood function.

function likelihoodGaussianRegressionConstructor(simulatorLikelihoodDefinition,emulatorRebuildCount,polynomialOrder,sigmaBuffer,logLikelihoodBuffer,logLikelihoodErrorTolerance,reportCount,emulateOutliers,configFileName)
  !% Constructor for Gaussian regression likelihood class.
  use FoX_DOM
  implicit none
  type            (likelihoodGaussianRegression)                          :: likelihoodGaussianRegressionConstructor
  type            (node                        ), intent(in   ), pointer  :: simulatorLikelihoodDefinition
  integer                                       , intent(in   )           :: emulatorRebuildCount                   , polynomialOrder    , &
       &                                                                     reportCount
  double precision                              , intent(in   )           :: sigmaBuffer                            , logLikelihoodBuffer, &
       &                                                                     logLikelihoodErrorTolerance
  logical                                       , intent(in   )           :: emulateOutliers
  type            (varying_string              ), intent(in   ), optional :: configFileName
  
  ! Set emulator properties.
  likelihoodGaussianRegressionConstructor%simulatorLikelihoodFunction => likelihoodNew(simulatorLikelihoodDefinition,configFileName)
  likelihoodGaussianRegressionConstructor%emulatorRebuildCount        =  emulatorRebuildCount
  likelihoodGaussianRegressionConstructor%polynomialOrder             =  polynomialOrder
  likelihoodGaussianRegressionConstructor%sigmaBuffer                 =  sigmaBuffer
  likelihoodGaussianRegressionConstructor%logLikelihoodBuffer         =  logLikelihoodBuffer
  likelihoodGaussianRegressionConstructor%logLikelihoodErrorTolerance =  logLikelihoodErrorTolerance
  likelihoodGaussianRegressionConstructor%reportCount                 =  reportCount
  likelihoodGaussianRegressionConstructor%emulateOutliers             =  emulateOutliers
  ! Initialize state and counters.
  likelihoodGaussianRegressionConstructor%initialized                 =  .false.
  likelihoodGaussianRegressionConstructor%accumulatedStateCount       =  0
  likelihoodGaussianRegressionConstructor%simulationCount             =  0
  likelihoodGaussianRegressionConstructor%evaluationCount             =  0
  likelihoodGaussianRegressionConstructor%emulatorCheckCount          =  0
  likelihoodGaussianRegressionConstructor%emulatorFailCount           =  0
  return
end function likelihoodGaussianRegressionConstructor

subroutine likelihoodGaussianRegressionDestructor(self)
  !% Destructor for Gaussian regression likelihood class.
  implicit none
  type(likelihoodGaussianRegression), intent(inout) :: self
    
  deallocate(self%simulatorLikelihoodFunction)
  call FGSL_Matrix_Free     (self%regressionMatrixLU)
  call FGSL_Permutation_Free(self%permutations      )
  return
end subroutine likelihoodGaussianRegressionDestructor

double precision function likelihoodGaussianRegressionEvaluate(self,simulationState,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  !% Return the log-likelihood for a Gaussian regression likelihood function.
  use, intrinsic :: ISO_C_Binding
  use Constraints_State
  use Constraints_Convergence
  use Memory_Management
  use MPI_Utilities
  use Galacticus_Display
  use Error_Functions
  implicit none
  class           (likelihoodGaussianRegression  ), intent(inout)                   :: self
  class           (state                         ), intent(in   )                   :: simulationState
  class           (convergence                   ), intent(inout)                   :: simulationConvergence
  double precision                                , intent(in   )                   :: temperature                      , logLikelihoodCurrent, &
       &                                                                               logPriorCurrent                  ,&
       & logPriorProposed
  double precision                                , allocatable  , dimension(:,:,:) :: states
  double precision                                , allocatable  , dimension(  :,:) :: likelihoods
  double precision                                , allocatable  , dimension(    :) :: likelihoodsCombined              , workspace            , &
       &                                                                               likelihoodsFitted                , stateCurrent         , &
       &                                                                               separations                      , semiVariances
  integer                                         , allocatable  , dimension(    :) :: stateCount
  double precision                                , parameter                       :: emulatorFailureSignificance=3.0d0
  type            (polynomialIterator            )                                  :: iterator1                        , iterator2
  type            (vector                        )                                  :: likelihoodSums                   , coefficients
  type            (matrix                        )                                  :: stateSums
  type            (fgsl_vector                   )                                  :: stateOffsetVector                , weightVector
  integer         (kind=FGSL_Int                 )                                  :: decompositionSign                , status
  integer                                                                           :: accumulatedStateCount            , i                    , &
       &                                                                               j                                , k                    , &
       &                                                                               evaluationsTotal                 , simulationsTotal     , &
       &                                                                               stateCountAccept                 , stateCountAcceptChain, &
       &                                                                               determinantSign                  , checksTotal          , &
       &                                                                               failuresTotal
  double precision                                                                  :: separation                       , likelihoodError      , & 
       &                                                                               likelihoodFit                    , likelihoodEmulated
  logical                                                                           :: likelihoodIsSimulated
  character       (len=8                         )                                  :: label
  type            (varying_string                )                                  :: message

  ! Report on emulation efficiency.
  if (mod(self%evaluationCount,self%reportCount) == 0 .and. self%evaluationCount > 0 .and. Galacticus_Verbosity_Level() >= verbosityInfo) then
     evaluationsTotal=mpiSelf%sum(self%evaluationCount)
     simulationsTotal=mpiSelf%sum(self%simulationCount)
     if (mpiSelf%isMaster()) then
        write (label,'(i8)') simulationsTotal
        message='Simulated/Emulated '//trim(adjustl(label))//'/'
        write (label,'(i8)') evaluationsTotal-simulationsTotal
        message=message//trim(adjustl(label))//' of '
        write (label,'(i8)') evaluationsTotal
        message=message//trim(adjustl(label))//' states ('
        write (label,'(f6.2)') 100.0d0*dble(simulationsTotal)/dble(evaluationsTotal)
        message=message//trim(adjustl(label))//'/'
        write (label,'(f6.2)') 100.0d0*dble(evaluationsTotal-simulationsTotal)/dble(evaluationsTotal)
        message=message//trim(adjustl(label))//'%)'
        call Galacticus_Display_Message(message)
     end if
  end if
  ! Rebuild emulator if necessary.
  accumulatedStateCount=mpiSelf%sum(self%accumulatedStateCount)
  if (accumulatedStateCount >= self%emulatorRebuildCount) then
     ! Initialize Gaussian regression workspace.
     if (.not.self%initialized) then
        ! Determine the size of the regression matrix.
        self%regressionMatrixSize=self%emulatorRebuildCount+1
        ! Count number of terms in the N-dimensional polynomial
        self%polynomialCoefficientCount=0
        do i=0,self%polynomialOrder
           self%polynomialCoefficientCount=                                  &
                &  self%polynomialCoefficientCount                           &
                & +polynomialCoefficientCount(i,simulationState%dimension())
        end do
        call Alloc_Array(self%polynomialCoefficient,    [self%polynomialCoefficientCount                                ] )
        call Alloc_Array(self%likelihoodSums       ,    [self%polynomialCoefficientCount                                ] )
        call Alloc_Array(self%stateSums            ,    [self%polynomialCoefficientCount,self%polynomialCoefficientCount] )
        call Alloc_Array(self%coefficients         ,    [self%polynomialCoefficientCount                                ] )
        call Alloc_Array(self%regressionMatrix     ,int([self%regressionMatrixSize      ,self%regressionMatrixSize      ]))
        call Alloc_Array(self%stateOffset          ,int([self%regressionMatrixSize                                      ]))
        call Alloc_Array(self%weight               ,int([self%regressionMatrixSize                                      ]))
        call Alloc_Array(self%likelihoodResiduals  ,    [self%emulatorRebuildCount                                      ] )
        call Alloc_Array(self%statesCombined       ,    [simulationState%dimension()    ,self%emulatorRebuildCount      ] )
        call Alloc_Array(self%stateScales          ,    [simulationState%dimension()                                    ] )
        self%initialized=.true.
     end if
     ! Gather the simulator state and likelihood from other chains.
     allocate(states     (simulationState%dimension(),self%emulatorRebuildCount,0:mpiSelf%count()-1))
     allocate(likelihoods(                            self%emulatorRebuildCount,0:mpiSelf%count()-1))
     allocate(stateCount (                                                      0:mpiSelf%count()-1))
     states                    =mpiSelf%gather(self%simulationState    )
     likelihoods               =mpiSelf%gather(self%simulatorLikelihood)
     stateCount                =0
     stateCount(mpiSelf%rank())=self   %accumulatedStateCount
     stateCount                =mpiSelf%sum   (stateCount              )
     ! Allocate workspace arrays.
     allocate(likelihoodsCombined( self%emulatorRebuildCount                                 ))
     allocate(likelihoodsFitted  ( self%emulatorRebuildCount                                 ))
     allocate(workspace          ( self%emulatorRebuildCount                                 ))
     allocate(separations        ((self%emulatorRebuildCount*(self%emulatorRebuildCount-1))/2))
     allocate(semiVariances      ((self%emulatorRebuildCount*(self%emulatorRebuildCount-1))/2))
     ! Extract likelihoods and states. We may have more states stored than we actually need. The following algorithm combines
     ! states from chains such that we balance the states used across chains as much as possible. That is, we find the minimum
     ! number of states available from all chains and take that number of states from each chain. If more states are needed, we
     ! find the minimum number of remaining states available from all chains and take that number of states from each chain. The
     ! process is repeated until we have accumulated enough states.
     j=0
     do while (j < self%emulatorRebuildCount)
        stateCountAccept=minval(stateCount,mask=stateCount > 0)
        do i=0,mpiSelf%count()-1
           stateCountAcceptChain=min(stateCountAccept,self%emulatorRebuildCount-j,stateCount(i))
           likelihoodsCombined(  j+1:j+stateCountAcceptChain)=likelihoods(  stateCount(i)+1-stateCountAcceptChain:stateCount(i),i)
           self%statesCombined(:,j+1:j+stateCountAcceptChain)=states     (:,stateCount(i)+1-stateCountAcceptChain:stateCount(i),i)
           stateCount(i)=stateCount(i)-stateCountAcceptChain
           j            =j            +stateCountAcceptChain
        end do
     end do
     ! Deallocate workspace.
     deallocate(likelihoods)
     deallocate(stateCount ) 
     deallocate(states     )
     ! Evaluate sums needed in fitting polynomial trend model.
     iterator1=polynomialIterator(self%polynomialOrder,simulationState%dimension())
     iterator2=polynomialIterator(self%polynomialOrder,simulationState%dimension())
     do while (iterator1%iterate())
        workspace=likelihoodsCombined
        do j=1,iterator1%currentOrder()
           workspace=workspace*self%statesCombined(iterator1%index(j),:)
        end do
        self%likelihoodSums(iterator1%counter())=sum(workspace)
        call iterator2%reset()
        do while (iterator2%iterate())
           workspace=1.0d0
           do j=1,iterator1%currentOrder()
              workspace=workspace*self%statesCombined(iterator1%index(j),:)
           end do
           do j=1,iterator2%currentOrder()
              workspace=workspace*self%statesCombined(iterator2%index(j),:)
           end do
           self%stateSums(iterator1%counter(),iterator2%counter())=sum(workspace)
        end do
     end do
     ! Assign vector and matrix in our linear system.
     likelihoodSums   =self%likelihoodSums
     stateSums        =self%stateSums
     ! Solve for regression coefficients.
     coefficients     =stateSums%linearSystemSolve(likelihoodSums)
     self%coefficients=coefficients
     ! Compute fitted likelihoods.
     likelihoodsFitted=0.0d0
     call iterator1%reset()
     do while (iterator1%iterate())    
        workspace=self%coefficients(iterator1%counter())
        do j=1,iterator1%currentOrder()
           workspace=workspace*self%statesCombined(iterator1%index(j),:)
        end do
        likelihoodsFitted=likelihoodsFitted+workspace
     end do
     self%likelihoodResiduals=likelihoodsCombined-likelihoodsFitted
     ! Determine suitable scales for each dimension.
     self%stateScales=maxval(self%statesCombined,dim=2)-minval(self%statesCombined,dim=2)
     ! Compute the variogram.
     k=0
     do i=1,self%emulatorRebuildCount-1
        do j=i+1,self%emulatorRebuildCount
           k=k+1
           separations  (k)=likelihoodGaussianRegressionSeparation(self,self%statesCombined(:,i),self%statesCombined(:,j))
           semiVariances(k)= 0.5d0                         &
                &           *(                             &
                &             +self%likelihoodResiduals(i) &
                &             -self%likelihoodResiduals(j) &
                &            )**2
        end do
     end do
     ! Fit the variogram model.
     call likelihoodGaussianRegressionFitVariogram(self,separations,semiVariances)
     ! Compute regression matrix.
     do i=1,self%emulatorRebuildCount
        do j=1,self%emulatorRebuildCount
           separation                =likelihoodGaussianRegressionSeparation(self,self%statesCombined(:,i),self%statesCombined(:,j))
           self%regressionMatrix(i,j)=likelihoodGaussianRegressionCorrelation(self,separation)
       end do
     end do
     ! Find the LU decomposition of the regression matrix for later use.
     self%regressionMatrix(1:self%emulatorRebuildCount  ,  self%emulatorRebuildCount+1)=1.0d0
     self%regressionMatrix(  self%emulatorRebuildCount+1,1:self%emulatorRebuildCount  )=1.0d0
     self%regressionMatrix(  self%emulatorRebuildCount+1,  self%emulatorRebuildCount+1)=0.0d0
     if (FGSL_Well_Defined(self%regressionMatrixLU)) call FGSL_Matrix_Free     (self%regressionMatrixLU)
     if (FGSL_Well_Defined(self%permutations      )) call FGSL_Permutation_Free(self%permutations      )
     self%regressionMatrixLU=FGSL_Matrix_Init      (type=1.0_FGSL_Double)
     self%permutations      =FGSL_Permutation_Alloc(self%regressionMatrixSize)
     status                 =FGSL_Matrix_Align     (self%regressionMatrix,self%regressionMatrixSize,self%regressionMatrixSize,self%regressionMatrixSize,self%regressionMatrixLU)
     status                 =FGSL_LinAlg_LU_Decomp (self%regressionMatrixLU,self%permutations,decompositionSign)     
     determinantSign        =FGSL_Linalg_LU_SgnDet (self%regressionMatrixLU                  ,decompositionSign)
     self%regressionMatrixIsSingular=(determinantSign == 0)
     ! Retain every second of the currently stored states. This avoids losing all information and only replacing it with data from
     ! regions where the new emulator is bad.
     do i=1,self%accumulatedStateCount/2
        self%simulatorLikelihood(  i)=self%simulatorLikelihood(  2*i)
        self%simulationState    (:,i)=self%simulationState    (:,2*i)
     end do
     self%accumulatedStateCount=self%accumulatedStateCount/2
     ! Deallocate workspace.
     deallocate(likelihoodsCombined)
     deallocate(likelihoodsFitted  )
     deallocate(workspace          )
     deallocate(semiVariances      )
     deallocate(separations        )
  end if
  ! Count evaluations.
  self%evaluationCount=self%evaluationCount+1
  ! If prior probability is impossible, then don't even try to emulate/simulate.
  if (logPriorProposed <= logImpossible) then
     likelihoodGaussianRegressionEvaluate=logImpossible
     return
  end if
  ! Ensure arrays are allocated.
  if (.not.allocated(self%simulatorLikelihood)) then
     call Alloc_Array(self%simulatorLikelihood,[                            self%emulatorRebuildCount])
     call Alloc_Array(self%simulationState    ,[simulationState%dimension(),self%emulatorRebuildCount])
  end if
  likelihoodIsSimulated=.false.
  if (self%initialized.and..not.self%regressionMatrixIsSingular) then     
     ! Compute vector D.
     stateCurrent=simulationState%get()
     do i=1,self%emulatorRebuildCount
        separation         =likelihoodGaussianRegressionSeparation(self,stateCurrent,self%statesCombined(:,i))
        self%stateOffset(i)=likelihoodGaussianRegressionCorrelation(self,separation)
     end do
     self%stateOffset(self%emulatorRebuildCount+1)=1.0d0
     ! Solve the linear system.
     stateOffsetVector=FGSL_Vector_Init(type=1.0_FGSL_Double)
     weightVector     =FGSL_Vector_Init(type=1.0_FGSL_Double)
     status           =FGSL_Vector_Align(self%stateOffset,self%regressionMatrixSize,stateOffsetVector,self%regressionMatrixSize,0_FGSL_Size_T,1_FGSL_Size_T)
     status           =FGSL_Vector_Align(self%weight     ,self%regressionMatrixSize,weightVector     ,self%regressionMatrixSize,0_FGSL_Size_T,1_FGSL_Size_T)
     status           =FGSL_LinAlg_LU_Solve(self%regressionMatrixLU,self%permutations,stateOffsetVector,weightVector)
     call FGSL_Vector_Free(weightVector     )
     call FGSL_Vector_Free(stateOffsetVector)
     ! Compute the likelihood and variance.
     likelihoodGaussianRegressionEvaluate=sum(self%likelihoodResiduals*self%weight(1:self%emulatorRebuildCount))
     likelihoodError                     =sqrt(likelihoodGaussianRegressionCorrelation(self,0.0d0)-sum(self%weight*self%stateOffset))
     iterator1=polynomialIterator(self%polynomialOrder,simulationState%dimension())
     do while (iterator1%iterate())
        likelihoodFit=self%coefficients(iterator1%counter())
        do j=1,iterator1%currentOrder()
           likelihoodFit=likelihoodFit*stateCurrent(iterator1%index(j))
        end do
        likelihoodGaussianRegressionEvaluate=likelihoodGaussianRegressionEvaluate+likelihoodFit        
     end do
     ! Test likelihood emulation.
     if (mod(self%evaluationCount,self%reportCount) == 0 .and. self%evaluationCount > 0 .and. Galacticus_Verbosity_Level() >= verbosityInfo) then
        ! Every so many steps we evaluate the simulated likelihood and check that our emulator is reliable.
        likelihoodEmulated                  =likelihoodGaussianRegressionEvaluate
        likelihoodGaussianRegressionEvaluate=self%simulatorLikelihoodFunction%evaluate(simulationState,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
        likelihoodIsSimulated               =.true.
        ! Count number of emulator checks and the failure rate.
        self%emulatorCheckCount=self%emulatorCheckCount+1
        if (abs(likelihoodEmulated-likelihoodGaussianRegressionEvaluate) > emulatorFailureSignificance*likelihoodError) self%emulatorFailCount=self%emulatorFailCount+1
             checksTotal  =mpiSelf%sum(self%emulatorCheckCount)
        failuresTotal=mpiSelf%sum(self%emulatorFailCount )
        if (mpiSelf%isMaster()) then
           write (label,'(i8)') failuresTotal
           message='Emulator failed '//trim(adjustl(label))//' times out of '
           write (label,'(i8)') checksTotal
           message=message//trim(adjustl(label))//' checks ('
           write (label,'(f6.2)') 100.0d0*dble(failuresTotal)/dble(checksTotal)
           message=message//trim(adjustl(label))//'% - expect '
           write (label,'(f6.2)') 100.0d0*(1.0d0-Error_Function(emulatorFailureSignificance/sqrt(2.0d0)))
           message=message//trim(adjustl(label))//'% for perfect emulator)'
           call Galacticus_Display_Message(message)
        end if
     else
        ! If simulation is converged, this is an outlier chain, and we're told to emulate such outliers, then return.
        if (simulationConvergence%isConverged().and.simulationConvergence%chainIsOutlier(mpiSelf%rank()).and.self%emulateOutliers) return
        ! If likelihood is well below current likelihood (and this should be sufficient that changes in priors won't affect the
        ! conclusion), then use the emulated likelihood. We scale the likelihood buffer by the current temperature to reflect the
        ! fact that transitions between states are easier when the temperature is high.
        if (likelihoodGaussianRegressionEvaluate+logPriorProposed+self%sigmaBuffer*likelihoodError < logLikelihoodCurrent+logPriorCurrent-self%logLikelihoodBuffer        *temperature) return
        ! Return if the error is below the tolerance. We increase the tolerance value in proportion to temperature since the
        ! likelihoods will be divided by this amount when evaluating transition probabilties.
        if (                                                                       likelihoodError <                                      self%logLikelihoodErrorTolerance*temperature) return
     end if
  end if
  ! Evaluate the likelihood using the simulator, unless we have already done so.
  self%simulationCount=self%simulationCount+1
  if (.not.likelihoodIsSimulated) likelihoodGaussianRegressionEvaluate=self%simulatorLikelihoodFunction%evaluate(simulationState,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  ! Store the likelihood and state.
  self%accumulatedStateCount                            =self%accumulatedStateCount+1
  self%simulatorLikelihood(  self%accumulatedStateCount)=likelihoodGaussianRegressionEvaluate
  self%simulationState    (:,self%accumulatedStateCount)=simulationState%get()
 return
end function likelihoodGaussianRegressionEvaluate

elemental double precision function likelihoodGaussianRegressionEvaluateVariogram(self,h)
  !% Compute the variogram at separation {\tt h}.
  implicit none
  class           (likelihoodGaussianRegression), intent(in   ) :: self
  double precision                              , intent(in   ) :: h
  
  if (h == 0.0d0) then
     likelihoodGaussianRegressionEvaluateVariogram=0.0d0
  else if (h < self%CR) then
     likelihoodGaussianRegressionEvaluateVariogram=self%C0+self%C1*(1.5d0*h/self%CR-0.5d0*(h/self%CR)**3)
  else
     likelihoodGaussianRegressionEvaluateVariogram=self%C0+self%C1
  end if
  return
end function likelihoodGaussianRegressionEvaluateVariogram

elemental double precision function likelihoodGaussianRegressionCorrelation(self,h)
  !% Compute correlation of the variogram model.
  implicit none
  class           (likelihoodGaussianRegression), intent(in   ) :: self
  double precision                              , intent(in   ) :: h
  
  likelihoodGaussianRegressionCorrelation=self%C0+self%C1-likelihoodGaussianRegressionEvaluateVariogram(self,h)
  return
end function likelihoodGaussianRegressionCorrelation

double precision function likelihoodGaussianRegressionSeparation(self,state1,state2)
  !% Determine the separation between two state vectors.
  implicit none
  class           (likelihoodGaussianRegression)              , intent(in   ) :: self
  double precision                              , dimension(:), intent(in   ) :: state1, state2

  likelihoodGaussianRegressionSeparation=sqrt(sum(((state1-state2)/self%stateScales)**2))
  return
end function likelihoodGaussianRegressionSeparation

subroutine likelihoodGaussianRegressionFitVariogram(self,separations,semiVariances)
  !% Compute best fit coefficients for the variogram model.
  use, intrinsic :: ISO_C_Binding
  use Galacticus_Error
  use FGSL
  implicit none
  class           (likelihoodGaussianRegression), intent(inout)               :: self
  double precision                              , intent(in   ), dimension(:) :: separations                , semiVariances
  double precision                              , target       , dimension(3) :: C
  real            (fgsl_double                 ), pointer      , dimension(:) :: cPtr
  integer                                       , parameter                   :: iterationsMaximum   =100
  double precision                              , parameter                   :: gradientTolerance   =1.0d-2
  type            (fgsl_multimin_function_fdf  )                              :: minimizationFunction
  type            (fgsl_multimin_fdfminimizer  )                              :: minimizer
  type            (fgsl_vector                 )                              :: cVector
  type            (c_ptr                       )                              :: parameters
  integer                                                                     :: status                     , iteration
  double precision                                                            :: currentMinimum

  ! Allocate workspace.  
  allocate(separationsNormalized  (size(separations)))
  allocate(separationsLimited     (size(separations)))
  allocate(semiVariancesNormalized(size(separations)))
  ! Compute normalized separations and variances.
  separationNormalization  =sum(separations  )/dble(size(separations))
  semiVarianceNormalization=sum(semiVariances)/dble(size(separations))
  separationsNormalized    =separations       /separationNormalization
  semiVariancesNormalized  =semiVariances     /semiVarianceNormalization
  ! Initialize the minimizer.
  minimizationFunction=FGSL_MultiMin_Function_FDF_Init(                                              &
       &                                               likelihoodGaussianRegressionVariogramModelF , &
       &                                               likelihoodGaussianRegressionVariogramModelD , &
       &                                               likelihoodGaussianRegressionVariogramModelFD, &
       &                                               3_FGSL_Size_T                               , &
       &                                               parameters                                    &
       &                                              )
  minimizer=FGSL_MultiMin_FDFMinimizer_Alloc(FGSL_MultiMin_FDFMinimizer_Conjugate_PR,3_FGSL_Size_T)
  cVector  =FGSL_Vector_Init(type=1.0_FGSL_Double)
  C        =[0.5d0,1.0d0,0.5d0] ! Initial guess for the parameters.
  status   =FGSL_Vector_Align(C,3_FGSL_Size_T,cVector,3_FGSL_Size_T,0_FGSL_Size_T,1_FGSL_Size_T)
  if (status /= FGSL_Success) call Galacticus_Error_Report('likelihoodGaussianRegressionFitVariogram','failed to initialize parameter vector')
  status   =FGSL_MultiMin_FDFMinimizer_Set(minimizer,minimizationFunction,cVector,0.01d0,0.1d0)
  if (status /= FGSL_Success) call Galacticus_Error_Report('likelihoodGaussianRegressionFitVariogram','failed to set minimizer')
  ! Iterate the minimizer until a sufficiently good solution is found.
  currentMinimum=0.0d0
  iteration     =0
  do while (                                                                                                                                   &
       &     FGSL_MultiMin_Test_Gradient(FGSL_MultiMin_FDFMinimizer_Gradient(minimizer),gradientTolerance*currentMinimum) == FGSL_CONTINUE     &
       &    .and.                                                                                                                              &
       &     iteration                                                                                                    <  iterationsMaximum &
       &   )
     iteration     =iteration+1
     status        =FGSL_MultiMin_FDFMinimizer_Iterate(minimizer)
     if (status == FGSL_ENoProg) exit
     if (status /= FGSL_Success) call Galacticus_Error_Report('likelihoodGaussianRegressionFitVariogram','failed to iterate minimizer')
     currentMinimum=FGSL_MultiMin_FDFMinimizer_Minimum(minimizer)
  end do
  ! Extract the best fit parameters.
  call FGSL_Vector_Free(cVector)
  cVector=FGSL_Vector_Init(type=1.0_FGSL_Double)
  cVector=FGSL_Multimin_FDFMinimizer_x(minimizer   )
  status =FGSL_Vector_Align           (cPtr,cVector)
  self%C0=cPtr(1)*semiVarianceNormalization
  self%C1=cPtr(2)*semiVarianceNormalization
  self%CR=cPtr(3)*separationNormalization
  ! Clean up.
  call FGSL_MultiMin_FDFMinimizer_Free(minimizer           )
  call FGSL_MultiMin_Function_FDF_Free(minimizationFunction)
  deallocate(separationsNormalized  )
  deallocate(separationsLimited     )
  deallocate(semiVariancesNormalized)
  return
end subroutine likelihoodGaussianRegressionFitVariogram

function likelihoodGaussianRegressionVariogramModelF(x, parameters) bind(c)
  use, intrinsic :: ISO_C_Binding
  implicit none
  type(c_ptr), value :: x, parameters
  real(FGSL_double), pointer, dimension(:) :: xx
  real(c_double) :: likelihoodGaussianRegressionVariogramModelF
  type(FGSL_vector) :: vec
  integer(FGSL_int) :: status
  call FGSL_obj_c_ptr(vec, x)
  status = FGSL_vector_align(xx, vec)
  
  if (any(xx < 0.0d0)) then
     likelihoodGaussianRegressionVariogramModelF=1.0d30
  else
     where (separationsNormalized > xx(3))
        separationsLimited=xx(3)
     elsewhere
        separationsLimited=separationsNormalized
     end where
     likelihoodGaussianRegressionVariogramModelF=sum((xx(1)+xx(2)*(1.5d0*separationsLimited/xx(3)-0.5d0*(separationsLimited/xx(3))**3)-semiVariancesNormalized)**2)
  end if
end function likelihoodGaussianRegressionVariogramModelF

subroutine likelihoodGaussianRegressionVariogramModelD(x, parameters, df) bind(c)
  use, intrinsic :: ISO_C_Binding
  implicit none
  type(c_ptr), value :: x, parameters, df
  type(FGSL_vector) :: vec, grad
  real(c_double), pointer, dimension(:) :: xx,ddf
  integer(FGSL_int) :: status
  
  
   call FGSL_obj_c_ptr(vec, x)
   call FGSL_obj_c_ptr(grad, df)
   status = FGSL_vector_align(xx, vec)
   status = FGSL_vector_align(ddf, grad)
   
   where (separationsNormalized > xx(3))
      separationsLimited=xx(3)
   elsewhere
      separationsLimited=separationsNormalized
   end where
   ddf(1)=2.0d0*sum((xx(1)+xx(2)*(1.5d0*separationsLimited/xx(3)-0.5d0*(separationsLimited/xx(3))**3)-semiVariancesNormalized))
   ddf(2)=2.0d0*sum((xx(1)+xx(2)*(1.5d0*separationsLimited/xx(3)-0.5d0*(separationsLimited/xx(3))**3)-semiVariancesNormalized)*(1.5d0*separationsLimited/xx(3)-0.5d0*(separationsLimited/xx(3))**3))
   ddf(3)=2.0d0*sum((xx(1)+xx(2)*(1.5d0*separationsLimited/xx(3)-0.5d0*(separationsLimited/xx(3))**3)-semiVariancesNormalized)*xx(2)*(-1.5d0*separationsLimited/xx(3)+1.5d0*(separationsLimited/xx(3))**3)/xx(3))
   where(abs(ddf) < 1.0d-30)
      ddf=1.0d-30
   end where
 end subroutine  likelihoodGaussianRegressionVariogramModelD
 
 subroutine likelihoodGaussianRegressionVariogramModelFD(x, parameters, f, df) bind(c)
   use, intrinsic :: ISO_C_Binding
   implicit none
   real(c_double) :: f
   type(c_ptr), value :: x, parameters, df
  
      f = likelihoodGaussianRegressionVariogramModelF(x, parameters)
    call likelihoodGaussianRegressionVariogramModelD(x, parameters, df)
 
end subroutine likelihoodGaussianRegressionVariogramModelFD

integer function polynomialCoefficientCount(n,d)
  !% Return the number of coefficients at {\tt n}$^{\rm th}$ order in polynomial of dimension {\tt d}.
  use Factorials
  implicit none
  integer, intent(in   ) :: n,d

  polynomialCoefficientCount=nint(Factorial(d+n-1)/Factorial(n)/Factorial(d-1))
  return
end function polynomialCoefficientCount

function polynomialIteratorConstructor(order,rank)
  !% Create a polynomial iterator for a poynomial of specified {\tt order} and {\tt rank}.
  use Memory_Management
  implicit none
  type   (polynomialIterator)                :: polynomialIteratorConstructor
  integer                    , intent(in   ) :: order                        , rank

  call Alloc_Array(polynomialIteratorConstructor%indices,[order])
  polynomialIteratorConstructor     %order=order
  polynomialIteratorConstructor     %rank =rank
  call polynomialIteratorConstructor%reset()
  return
end function polynomialIteratorConstructor

subroutine polynomialIteratorReset(self)
  !% Reset a polynomial iterator.
  implicit none
  class(polynomialIterator), intent(inout) :: self

  self%orderCurrent=-1
  self%stateCurrent=-1
  self%count       = 0
  self%indices     = 0
  return
end subroutine polynomialIteratorReset

logical function polynomialIteratorIterate(self)
  !% Iterate over polynomial coefficients.
  implicit none
  class  (polynomialIterator), intent(inout) :: self
  integer                                    :: j

  polynomialIteratorIterate=.true.
  self%count=self%count+1
  if     (                                                                              &
       &   self%stateCurrent <  0                                                       &
       &  .or.                                                                          &
       &   self%stateCurrent >= polynomialCoefficientCount(self%orderCurrent,self%rank) &
       & ) then
     ! Move to next polynomial order.
     self%orderCurrent=self%orderCurrent+1
     if (self%orderCurrent > self%order) then
        polynomialIteratorIterate=.false.
     else
        self%stateCurrent                     =1
        self%indices     (1:self%orderCurrent)=1
     end if
  else
     self%stateCurrent=self%stateCurrent+1
     j=self%orderCurrent
     do while (j > 0)
        self%indices(j)=self%indices(j)+1
        if (self%indices(j) > self%rank) then
           self%indices(j)=self%indices(j-1)+1
           j=j-1
        else
           j=0
        end if
     end do
  end if
  return
end function polynomialIteratorIterate

integer function polynomialIteratorIndex(self,i)
  !% Return the requested index of a polynomial iterator.
  implicit none
  class  (polynomialIterator), intent(inout) :: self
  integer                    , intent(in   ) :: i

  polynomialIteratorIndex=self%indices(i)
  return
end function polynomialIteratorIndex

integer function polynomialIteratorCurrentOrder(self)
  !% Return the current order of a polynomial iterator.
  implicit none
  class  (polynomialIterator), intent(inout) :: self

  polynomialIteratorCurrentOrder=self%orderCurrent
  return
end function polynomialIteratorCurrentOrder

integer function polynomialIteratorCounter(self)
  !% Return the current count of a polynomial iterator.
  implicit none
  class  (polynomialIterator), intent(inout) :: self

  polynomialIteratorCounter=self%count
  return
end function polynomialIteratorCounter
