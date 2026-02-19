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

  !!{
  Implementation of a posterior sampling convergence class which implements the Gelman-Rubin statistic.
  !!}

  use :: ISO_Varying_String, only : varying_string

  !![
  <posteriorSampleConvergence name="posteriorSampleConvergenceGelmanRubin">
   <description>
    This class adopts the convergence criterion proposed by
    \citeauthor{gelman_a._inference_1992}~(\citeyear{gelman_a._inference_1992}; see also \citealt{brooks_general_1998}), which
    compares the variance in parameter values within chains to that between chains. Outlier detection is applied to the chains using a
    standard Grubb's outlier test. The behavior of this criterion is controlled by the following subparameters:
    \begin{description}
    \item [{\normalfont \ttfamily Rhat}] The correlation coefficient, $\hat{R}$, value at which to declare convergence.
    \item [{\normalfont \ttfamily burnCount}] Set number of steps to burn before applying the convergence test.
    \item [{\normalfont \ttfamily testCount}] Set the number of steps between successive applications of the convergence test.
    \item [{\normalfont \ttfamily outlierSignificance}] The significance level required in outlier detection.
    \item [{\normalfont \ttfamily outlierLogLikelihoodOffset}] The offset in log-likelihood from the current maximum likelihood chain
      required for a chain to be declared to be an outlier.
    \item [{\normalfont \ttfamily outlierCountMaximum}] The maximum number of outlier chains allowed.
    \end{description}
   </description>
  </posteriorSampleConvergence>
  !!]
  type, extends(posteriorSampleConvergenceClass) :: posteriorSampleConvergenceGelmanRubin
     !!{
     Implementation of a posterior sampling convergence class which implements the Gelman-Rubin statistic.
     !!}
     private
     double precision                                            :: thresholdHatR             , outlierSignificance         , &
          &                                                         outlierLogLikelihoodOffset
     integer                                                     :: burnCount                 , testCount                   , &
          &                                                         stepCount                 , outlierCountMaximum         , &
          &                                                         reportCount               , estimateCount               , &
          &                                                         logFileUnit               , convergedAtStepCount
     logical                                                     :: converged                 , logFileIsOpen       =.false.
     type            (varying_string)                            :: logFileName
     double precision                , allocatable, dimension(:) :: correctedHatR
     logical                         , allocatable, dimension(:) :: chainMask
   contains
     !![
     <methods>
       <method description="Return the current convergence measure, $\hat{R}$." method="convergenceMeasure" />
       <method description="Return the target convergence measure, $\hat{R}$." method="convergenceMeasureTarget" />
     </methods>
     !!]
     final     ::                             gelmanRubinDestructor
     procedure :: isConverged              => gelmanRubinIsConverged
     procedure :: convergedAtStep          => gelmanRubinConvergedAtStep
     procedure :: reset                    => gelmanRubinReset
     procedure :: logReport                => gelmanRubinLogReport
     procedure :: stateIsOutlier           => gelmanRubinStateIsOutlier
     procedure :: convergenceMeasure       => gelmanRubinConvergenceMeasure
     procedure :: convergenceMeasureTarget => gelmanRubinConvergenceMeasureTarget
  end type posteriorSampleConvergenceGelmanRubin

  interface posteriorSampleConvergenceGelmanRubin
     !!{
     Constructors for the \refClass{posteriorSampleConvergenceGelmanRubin} posterior sampling convergence class.
     !!}
     module procedure gelmanRubinConstructorParameters
     module procedure gelmanRubinConstructorInternal
  end interface posteriorSampleConvergenceGelmanRubin

contains

  function gelmanRubinConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleConvergenceGelmanRubin} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: ISO_Varying_String, only : varying_string
    use :: Input_Parameters  , only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleConvergenceGelmanRubin)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: thresholdHatR             , outlierSignificance, &
         &                                                                    outlierLogLikelihoodOffset
    integer                                                                :: burnCount                 , testCount          , &
         &                                                                    outlierCountMaximum       , reportCount
    type            (varying_string                       )                :: logFileName

    !![
    <inputParameter>
      <name>thresholdHatR</name>
      <defaultValue>1.2d0</defaultValue>
      <description>The $\hat{R}$ value at which convergence is declared.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>burnCount</name>
      <defaultValue>0</defaultValue>
      <description>The number of steps to burn before computing convergence.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>testCount</name>
      <defaultValue>10</defaultValue>
      <description>The interval in number of steps at which to check convergence.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outlierCountMaximum</name>
      <defaultValue>0</defaultValue>
      <description>The maximum number of outlier states allowed.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outlierSignificance</name>
      <defaultValue>0.05d0</defaultValue>
      <description>The significance at which to declare a state an outlier.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outlierLogLikelihoodOffset</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The log-likelihood offset at which to declare a state an outlier.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reportCount</name>
      <defaultValue>10</defaultValue>
      <description>The interval in number of steps at which to report on convergence status.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>logFileName</name>
      <description>The name of the file to which convergence state should be logged.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleConvergenceGelmanRubin(thresholdHatR,burnCount,testCount,outlierCountMaximum,outlierSignificance,outlierLogLikelihoodOffset,reportCount,logFileName)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function gelmanRubinConstructorParameters

  function gelmanRubinConstructorInternal(thresholdHatR,burnCount,testCount,outlierCountMaximum,outlierSignificance,outlierLogLikelihoodOffset,reportCount,logFileName) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleConvergenceGelmanRubin} convergence class.
    !!}
    use :: Error            , only : Error_Report
    use :: MPI_Utilities    , only : mpiSelf
    type            (posteriorSampleConvergenceGelmanRubin)                :: self
    double precision                                       , intent(in   ) :: thresholdHatR             , outlierSignificance, &
         &                                                                    outlierLogLikelihoodOffset
    integer                                                , intent(in   ) :: burnCount                 , testCount          , &
         &                                                                    outlierCountMaximum       , reportCount
    type            (varying_string                       ), intent(in   ) :: logFileName
    !![
    <constructorAssign variables="thresholdHatR,burnCount,testCount,outlierCountMaximum,outlierSignificance,outlierLogLikelihoodOffset,reportCount,logFileName"/>
    !!]

    allocate(self%chainMask(mpiSelf%count()))
    self%estimateCount       = 0
    self%converged           =.false.
    self%convergedAtStepCount=-1
    self%chainMask           =.true.
    self%stepCount           =testCount
    ! Validate.
    if (mpiSelf%count()-self%outlierCountMaximum < 3) call Error_Report('maximum number of outliers is too large'//{introspection:location})
    ! Open log file.
    if (mpiSelf%isMaster()) then
       open(newUnit=self%logFileUnit,file=char(logFileName),form='formatted',status='unknown')
       self%logFileIsOpen=.true.
    else
       self%logFileIsOpen=.false.
    end if
    return
  end function gelmanRubinConstructorInternal

  subroutine gelmanRubinDestructor(self)
    !!{
    Destroy a Gelman-Rubin convergence object.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    type(posteriorSampleConvergenceGelmanRubin), intent(inout) :: self

    ! Close the log file.
    if (mpiSelf%isMaster().and.self%logFileIsOpen) close(self%logFileUnit)
    return
  end subroutine gelmanRubinDestructor

  logical function gelmanRubinIsConverged(self,simulationState,logLikelihood)
    !!{
    Return whether the simulation is converged.
    !!}
    use :: Display                 , only : displayMessage
    use :: ISO_Varying_String      , only : varying_string
    use :: MPI_Utilities           , only : mpiBarrier                    , mpiSelf
    use :: Posterior_Sampling_State, only : posteriorSampleStateClass
    use :: Statistics_Distributions, only : distributionFunction1DStudentT
    use :: String_Handling         , only : operator(//)
    implicit none
    class           (posteriorSampleConvergenceGelmanRubin), intent(inout)               :: self
    class           (posteriorSampleStateClass            ), intent(inout), optional     :: simulationState
    double precision                                       , intent(in   ), optional     :: logLikelihood
    double precision                                       , allocatable  , dimension(:) :: chainMean                    , chainVariance          , &
         &                                                                                  interchainMean               , interchainMeanVariance , &
         &                                                                                  B                            , W                      , &
         &                                                                                  varianceW                    , covarianceVarianceMean , &
         &                                                                                  covarianceVarianceMeanSquared, interchainMeanSquared  , &
         &                                                                                  hatV                         , varianceHatV           , &
         &                                                                                  degreesOfFreedom             , hatR                   , &
         &                                                                                  currentStateVariance         , chainDeviation         , &
         &                                                                                  chainDeviationMaximum        , currentStateMean       , &
         &                                                                                  currentStateMeanSquared      , logLikelihoods
    integer                                                , allocatable  , dimension(:) :: chainDeviationIndex
    logical                                                , allocatable  , dimension(:) :: hasOutlier
    logical                                                                              :: newOutliersFound
    integer                                                                              :: outlierCount                 , activeChainCount       , &
         &                                                                                  i                            , deviationMaximumChain
    double precision                                                                     :: grubbsCriticalValue          , tStatisticCriticalValue, &
         &                                                                                  logLikelihoodMaximum         , deviationMaximum
    type            (distributionFunction1DStudentT       )                              :: studentT
    type            (varying_string                       )                              :: message
    character       (len=16                               )                              :: label

    ! If no arguments were provided, return current convergence status without updating.
    if (.not.(present(simulationState).and.present(logLikelihood))) then
       gelmanRubinIsConverged=self%converged
       return
    end if
    ! Check if we have enough steps to check convergence.
    if (simulationState%count() < self%burnCount) then
       gelmanRubinIsConverged=.false.
       return
    end if
    ! Decide if we should check convergence.
    self%stepCount=self%stepCount+1
    if (self%stepCount < self%testCount) then
       gelmanRubinIsConverged=self%converged
       return
    end if
    self%stepCount=0
    ! Allocate Rhat array.
    if (.not.allocated(self%correctedHatR)) allocate(self%correctedHatR(simulationState%dimension()))
    ! Find outlier chains using a Grubb's outlier test.
    ! Initialize chain mask to all chains accepted.
    self%chainMask=.true.
    ! Find the log likelihoods of all states.
    logLikelihoods=mpiSelf%gather(logLikelihood)
    ! Find the maximum likelihood over current states.
    logLikelihoodMaximum=mpiSelf%maxval(logLikelihood)
    ! Begin iterative removal of outliers.
    outlierCount    =0
    activeChainCount=mpiSelf%count()
    newOutliersFound=.true.
    allocate(currentStateMean       (simulationState%dimension()))
    allocate(currentStateMeanSquared(simulationState%dimension()))
    do while (outlierCount < self%outlierCountMaximum .and. newOutliersFound)
       ! Find interchain mean and variance of current states.
       currentStateMean       =mpiSelf%average(simulationState%get()   ,self%chainMask)
       currentStateMeanSquared=mpiSelf%average(simulationState%get()**2,self%chainMask)
       currentStateVariance   = (currentStateMeanSquared-currentStateMean**2) &
            &                  *dble(activeChainCount  )                      &
            &                  /dble(activeChainCount-1)
       allocate(chainDeviation(size(currentStateVariance)))
       where (currentStateVariance > 0.0d0)
          chainDeviation         = abs(simulationState%get()-currentStateMean)   &
               &                  /sqrt(currentStateVariance)
       elsewhere
          chainDeviation         = 0.0d0
       end where
       ! Find the maximum value of the deviation and which chain has this deviation.
       chainDeviationMaximum=mpiSelf%maxval(chainDeviation,self%chainMask)
       chainDeviationIndex  =mpiSelf%maxloc(chainDeviation,self%chainMask)+1
       deallocate(chainDeviation)
       ! Evaluate the critical value for outlier rejection.
       studentT               =distributionFunction1DStudentT( dble(activeChainCount-2))
       tStatisticCriticalValue=studentT%inverseUpper         (                            &
            &                                                 +(                          &
            &                                                   +1.0d0                    &
            &                                                   -self%outlierSignificance &
            &                                                  )                          &
            &                                                 /2.0d0                      &
            &                                                 /dble(activeChainCount  )   &
            &                                                )
       grubbsCriticalValue= dble(                              &
            &                         activeChainCount-1       &
            &                   )                              &
            &              /sqrt(                              &
            &                    dble(activeChainCount)        &
            &                   )                              &
            &              *sqrt(                              &
            &                       tStatisticCriticalValue**2 &
            &                    /(                            &
            &                      +  activeChainCount         &
            &                      -2.0d0                      &
            &                      +tStatisticCriticalValue**2 &
            &                     )                            &
            &                   )
       ! Test for outliers using the Grubb's statistic.
       hasOutlier=                                       &
            &        chainDeviationMaximum               &
            &      >                                     &
            &        grubbsCriticalValue                 &
            &     .and.                                  &
            &        logLikelihoods(chainDeviationIndex) &
            &      <                                     &
            &        logLikelihoodMaximum                &
            &       -                                    &
            &        self%outlierLogLikelihoodOffset
       ! Mask chains with outliers.
       newOutliersFound=any(hasOutlier)
       if (newOutliersFound) then
          deviationMaximum     =0.0d0
          deviationMaximumChain=-1
          do i=1,size(hasOutlier)
             if (hasOutlier(i).and.self%chainMask(chainDeviationIndex(i)).and.chainDeviationMaximum(i) > deviationMaximum) then
                deviationMaximum     =chainDeviationMaximum(i)
                deviationMaximumChain=chainDeviationIndex  (i)
             end if
          end do
          self%chainMask(deviationMaximumChain)=.false.
          outlierCount                         =outlierCount    +1
          activeChainCount                     =activeChainCount-1
       end if
       ! Proceed to next iteration.
       call mpiBarrier()
    end do
    ! Get the mean and variance of each parameter in our chain.
    chainMean    =simulationState%mean    ()
    chainVariance=simulationState%variance()
    ! Get the interchain mean.
    interchainMean       =mpiSelf%average(chainMean   ,self%chainMask)
    ! Get the interchain mean squared.
    interchainMeanSquared=mpiSelf%average(chainMean**2,self%chainMask)
    ! Get the interchain variance.
    interchainMeanVariance=(                          &
         &                  +interchainMeanSquared    &
         &                  -interchainMean       **2 &
         &                 )                          &
         &                 *dble(activeChainCount  )  &
         &                 /dble(activeChainCount-1)
    ! Compute B from Brooks & Gelman (section 1.2).
    B                     =dble(simulationState%count())*interchainMeanVariance
    ! Compute W from Brooks & Gelman (section 1.2).
    W                     = mpiSelf%average(chainVariance   ,self%chainMask)
    ! Check for zero chain variance.
    if (any(W <= 0.0d0)) then
       self%estimateCount=self%estimateCount+1
       if (mpiSelf%isMaster() .and. mod(self%estimateCount,self%reportCount) == 0) then
          write (label,'(i16)') simulationState%count()
          message="Gelman-Rubin statistic at "//trim(adjustl(label))
          message=message//" steps cannot be computed (zero variances)"
          call displayMessage(message)
       end if
       gelmanRubinIsConverged=.false.
       return
    end if
    ! Compute variance of chain variances.
    varianceW             =(                                                  &
         &                  +mpiSelf%average(chainVariance**2,self%chainMask) &
         &                  -                            W**2                 &
         &                 )                                                  &
         &                 *dble(activeChainCount  )                          &
         &                 /dble(activeChainCount-1)
    ! Find the covariance of chain variances and means.
    covarianceVarianceMean       = (mpiSelf%average(chainVariance*chainMean   ,self%chainMask)-W*interchainMean       ) &
         &                        *dble(activeChainCount  )                                                             &
         &                        /dble(activeChainCount-2)
    covarianceVarianceMeanSquared= (mpiSelf%average(chainVariance*chainMean**2,self%chainMask)-W*interchainMeanSquared) &
         &                        *dble(activeChainCount  )                                                             &
         &                        /dble(activeChainCount-2)
    ! Estimate Vhat and variance in Vhat from Gelman & Rubin (eqn. 4).
    hatV=  +W*dble(simulationState%count()-1)/dble(simulationState%count()) &
         & +B/dble(simulationState%count()  )                               &
         & +B/dble(simulationState%count()  )/dble(activeChainCount       )
    varianceHatV=                                      &
         &        (                                    &
         &          dble(simulationState%count()-1)    &
         &         /dble(simulationState%count()  )    &
         &        )**2                                 &
         &       *varianceW                            &
         &       /  dble(activeChainCount         )    &
         &       +(                                    &
         &          dble(activeChainCount       +1)    &
         &         /dble(activeChainCount         )    &
         &         /dble(simulationState%count()  )    &
         &        )**2                                 &
         &       *(                                    &
         &         2.0d0                               &
         &         /dble(activeChainCount       -1)    &
         &        )                                    &
         &       *B**2                                 &
         &       +2.0d0                                &
         &       *  dble(activeChainCount       +1)    &
         &       *  dble(simulationState%count()-1)    &
         &       /  dble(activeChainCount         )    &
         &       /  dble(simulationState%count()  )**2 &
         &       *(                                    &
         &          dble(simulationState%count()  )    &
         &         /dble(activeChainCount         )    &
         &        )                                    &
         &       *(                                    &
         &         +covarianceVarianceMeanSquared      &
         &         -2.0d0                              &
         &         *interchainMean                     &
         &         *covarianceVarianceMean             &
         &        )
    ! Compute degrees of freedom.
    degreesOfFreedom=2.0d0*hatV**2/varianceHatV
    ! Compute R-hat statistic.
    hatR=                                               &
         & sqrt(                                        &
         &      max(                                    &
         &           0.0d0                            , &
         &           (                                  &
         &             dble(activeChainCount        +1) &
         &            /dble(activeChainCount          ) &
         &           )                                  &
         &          *(                                  &
         &            (                                 &
         &              dble(simulationState%count()-1) &
         &             /dble(simulationState%count()  ) &
         &            )                                 &
         &            *W                                &
         &            +B                                &
         &            / dble(simulationState%count()  ) &
         &           )                                  &
         &           /W                                 &
         &           -  dble(simulationState%count()-1) &
         &           /  dble(simulationState%count()  ) &
         &           /  dble(activeChainCount         ) &
         &         )                                    &
         &     )
    ! Compute corrected R-hat statistic.
    self%correctedHatR=(degreesOfFreedom+3.0d0)*hatR/(degreesOfFreedom+1.0d0)
    ! Check for convergence.
    if (.not.self%converged.and.all(self%correctedHatR < self%thresholdHatR)) then
       self%converged           =.true.
       self%convergedAtStepCount=simulationState%count()
    end if
    gelmanRubinIsConverged=self%converged
    ! Report.
    self%estimateCount=self%estimateCount+1
    if (mpiSelf%isMaster() .and. mod(self%estimateCount,self%reportCount) == 0) then
       if (.not.self%converged) then
          write (label,'(i16)') simulationState%count()
          message="Gelman-Rubin ̂R at "//trim(adjustl(label))
          write (label,'(f6.2)') minval(self%correctedHatR)
          message=message//" steps min/max="//trim(adjustl(label))//"/"
          write (label,'(f6.2)') maxval(self%correctedHatR)
          message=message//trim(adjustl(label))//")"
          call displayMessage(message)
          if (activeChainCount < mpiSelf%count()) then
             message='outlier chains:'
             label=''
             do i=0,mpiSelf%count()-1
                if (.not.self%chainMask(i+1)) then
                   message=message//trim(label)//' '//i
                   label=","
                end if
             end do
             call displayMessage(message)
          else
             call displayMessage('no outlier chains')
          end if
       end if
       write (self%logFileUnit,*) "outliers    ",simulationState%count(),self%chainMask
       write (self%logFileUnit,*) "convergence ",simulationState%count(),minval(self%correctedHatR),maxval(self%correctedHatR),self%correctedHatR
       call flush(self%logFileUnit)
    end if
    return
  end function gelmanRubinIsConverged

  integer function gelmanRubinConvergedAtStep(self)
    !!{
    Return the step at which the simulation converged.
    !!}
    implicit none
    class(posteriorSampleConvergenceGelmanRubin), intent(inout) :: self

    gelmanRubinConvergedAtStep=self%convergedAtStepCount
    return
  end function gelmanRubinConvergedAtStep

  subroutine gelmanRubinReset(self)
    !!{
    Reset the convergence object.
    !!}
    implicit none
    class(posteriorSampleConvergenceGelmanRubin), intent(inout) :: self

    self%converged           =.false.
    self%convergedAtStepCount=-1
    return
  end subroutine gelmanRubinReset

  subroutine gelmanRubinLogReport(self,fileUnit)
    !!{
    Write a convergence report to the given {\normalfont \ttfamily fileUnit}.
    !!}
    implicit none
    class    (posteriorSampleConvergenceGelmanRubin), intent(inout) :: self
    integer                                         , intent(in   ) :: fileUnit
    character(len=25                               )                :: label

    if (size(self%correctedHatR) > 1) then
       write (label   ,'(a,i4.4,a)') '(a,f5.2,',size(self%correctedHatR)-1,'(", ",f5.2))'
    else
       write (label   ,'(a       )') '(a,f5.2)'
    end if
    write    (fileUnit,label       ) 'Gelman-Rubin covergence ̂R: ',self%correctedHatR
    write    (label   ,'(a,i4.4,a)') '(a,l1,'  ,size(self%chainMask    )-1,'(", ",l1))'
    write    (fileUnit,label       ) 'Gelman-Rubin chain mask  : ',self%chainMask
    return
  end subroutine gelmanRubinLogReport

  logical function gelmanRubinStateIsOutlier(self,stateIndex)
    !!{
    Return true if the specified chain is deemed to be an outlier.
    !!}
    implicit none
    class  (posteriorSampleConvergenceGelmanRubin), intent(inout) :: self
    integer                                       , intent(in   ) :: stateIndex

    gelmanRubinStateIsOutlier=.not.self%chainMask(stateIndex+1)
    return
  end function gelmanRubinStateIsOutlier

  double precision function gelmanRubinConvergenceMeasure(self)
    !!{
    Return the current maximum $\hat{R}$ convergence measure.
    !!}
    implicit none
    class           (posteriorSampleConvergenceGelmanRubin), intent(inout) :: self
    double precision                                       , parameter     :: convergenceMeasureLarge=100.0d0

    if (allocated(self%correctedHatR)) then
       gelmanRubinConvergenceMeasure=maxval(self%correctedHatR)
    else
       ! Convergence has not yet been computed - return a suitably large value.
       gelmanRubinConvergenceMeasure=convergenceMeasureLarge
    end if
    return
  end function gelmanRubinConvergenceMeasure

  double precision function gelmanRubinConvergenceMeasureTarget(self)
    !!{
    Return the target $\hat{R}$ convergence measure.
    !!}
    implicit none
    class(posteriorSampleConvergenceGelmanRubin), intent(inout) :: self

    gelmanRubinConvergenceMeasureTarget=self%thresholdHatR
    return
  end function gelmanRubinConvergenceMeasureTarget
