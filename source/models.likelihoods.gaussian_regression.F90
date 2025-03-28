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
  Implementation of a posterior sampling likelihood class which implements a likelihood using Gaussian regression to emulate
  another likelihood.
  !!}

  use, intrinsic :: ISO_C_Binding        , only : c_size_t
  use            :: Linear_Algebra       , only : matrixLU
  use            :: Statistics_Variograms, only : variogramClass

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodGaussianRegression">
   <description>
    The likelihood is computed either using another likelihood function (the ``simulator'') or via Gaussian regression emulation of
    that simulator. The details of the emulation algorithm are specified by the following sub-parameters:
    \begin{description}
    \item[{\normalfont \ttfamily emulatorRebuildCount}] The number of simulator evaluations from which the emulator is built;
    \item[{\normalfont \ttfamily polynomialOrder}] The order of the polynomial fitted to the simulator likelihoods prior to Gaussian regression;
    \item[{\normalfont \ttfamily sigmaBuffer}] See below;
    \item[{\normalfont \ttfamily logLikelihoodBuffer}] See below;
    \item[{\normalfont \ttfamily logLikelihoodErrorTolerance}] See below;
    \item[{\normalfont \ttfamily reportCount}] The number of likelihood evaluations between successive reports on the status of the emulator;
    \item[{\normalfont \ttfamily emulateOutliers}] If true, then outlier chains are always emulated post-convergence (this is safe if such chains are not
     used in constructing proposals for non-outlier chains);
    \item[{\normalfont \ttfamily simulatorLikelihood}] Contains another likelihood function definition which will be used to construct the simulator.
    \end{description}
    
    In detail, this likelihood function first collects {\normalfont \ttfamily emulatorRebuildCount} likelihood evaluations from the
    simulator. It then fits a polynomial of order {\normalfont \ttfamily polynomialOrder} and of dimension equal to the dimension of
    the state vector to the simulated likelihoods. Gaussian regression is performed on the residuals of the simulated likelihoods
    after this polynomial fit is removed. Once the emulator has been built in this way every second simulated state is discarded, and
    accumulation of new simulated states continues. Once {\normalfont \ttfamily emulatorRebuildCount} simulated states have once again
    been accumulated a new simulator is built. This ensures that the emulator does not lose all information used in building the
    previous emulator\footnote{This would be unfortunate as the second emulator to be built would then contain information on only
      those regions of the state space that were poorly emulated before.}, instead information from older emulators decays
    exponentially.
    
    Once an emulator has been built, on each successive likelihood evaluation the emulated log-likelihood $\log\mathcal{L}_\mathrm{e}$
    and its error estimate $\sigma_{\log\mathcal{L}_\mathrm{e}}$ are computed. The emulated likelihood is then returned if:
    \begin{equation}
    \log\mathcal{P}^\prime + \log\mathcal{L}_\mathrm{e} + N \sigma_{\log\mathcal{L}_\mathrm{e}} &lt; \log\mathcal{P} + \log\mathcal{L} - T \Delta\log\mathcal{L},
    \end{equation}
    where $N=${\normalfont \ttfamily sigmaBuffer}, $\Delta\log\mathcal{L}=${\normalfont \ttfamily logLikelihoodBuffer}, $T$ is the
    temperature, $\log\mathcal{L}$ is the current log-likelihood, $\log\mathcal{P}$ is the current log-prior probability, and
    $\log\mathcal{P}^\prime$ is the proposed log-prior probability, or if
    \begin{equation}
    \sigma_{\log\mathcal{L}_\mathrm{e}} &lt; T \sigma_{\log\mathcal{L}},
    \end{equation}
    where $\sigma_{\log\mathcal{L}}=${\normalfont \ttfamily logLikelihoodErrorTolerance}, otherwise the simulator is used to compute
    the exact likelihood. In this way, the emulated likelihood is used if it is sufficiently below the current likelihood that, even
    accounting for the emulation error, transition to the new state is highly unlikely, or if the error on the likelihood emulation is
    sufficiently small that it will not have a significant effect on the transition probability to the proposed state.
    
    If verbosity is set to {\normalfont \ttfamily info} or greater than a report will be issued every {\normalfont \ttfamily
    reportCount} evaluations. The report will give the proportions of simulated vs. emulated evaluations. Additionally, during the
    evaluation where the report is issued, both the emulated and simulated log-likelihoods are evaluated and are tested to see if
    they lie within $3 \sigma_{\log\mathcal{L}_\mathrm{e}}$ of each other. The rate of failures (i.e. where the two differ by more
    than this amount) is then reported.
   </description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodGaussianRegression
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood using Gaussian regression to emulate
     another likelihood.
     !!}
     private
     class           (posteriorSampleLikelihoodClass), pointer                     :: posteriorSampleLikelihood_ => null()
     class           (variogramClass                ), pointer                     :: variogram_                 => null()
     integer                                                                       :: accumulatedStateCount               , emulatorRebuildCount       , &
          &                                                                           polynomialOrder                     , polynomialCoefficientCount , &
          &                                                                           reportCount                         , simulationCount            , &
          &                                                                           evaluationCount                     , emulatorCheckCount         , &
          &                                                                           emulatorFailCount                   , dumpEmulatorCount
     integer         (c_size_t                      )                              :: regressionMatrixSize
     double precision                                , allocatable, dimension(:  ) :: simulatorLikelihood                 , polynomialCoefficient      , &
          &                                                                           likelihoodSums                      , coefficients               , &
          &                                                                           stateOffset                         , weight                     , &
          &                                                                           likelihoodResiduals                 , stateScales                , &
          &                                                                           stateMeans
     double precision                                , allocatable, dimension(:,:) :: simulationState                     , stateSums                  , &
          &                                                                           regressionMatrix                    , statesCombined
     logical                                                                       :: initialized                         , regressionMatrixIsSingular , &
          &                                                                           isGood                              , assumeZeroVarianceAtZeroLag
     type            (matrixLU                      ), allocatable                 :: regressionMatrixLU
     double precision                                                              :: C0                                  , C1                         , &
          &                                                                           CR                                  , sigmaBuffer                , &
          &                                                                           logLikelihoodBuffer                 , logLikelihoodErrorTolerance
     logical                                                                       :: emulateOutliers                     , dumpEmulator               , &
          &                                                                           dummyEmulator
     type            (varying_string                 )                             :: dumpEmulatorFileRoot
     ! Workspaces used when fitting the semi-variogram.
     double precision                                 , allocatable, dimension(:) :: separationsNormalized                , semiVariancesNormalized    , &
          &                                                                          separationsBinned                    , semiVariancesBinned        , &
          &                                                                          separationsLimited
     double precision                                                             :: separationNormalization              , semiVarianceNormalization
     integer                                                                      :: binCount
   contains
     !![
     <methods>
       <method method="separation" description="Determine the separation between two state vectors."/>
       <method method="emulate"    description="Evaluate the model emulator."                       />
     </methods>
     !!]
     final     ::                    gaussianRegressionDestructor
     procedure :: evaluate        => gaussianRegressionEvaluate
     procedure :: functionChanged => gaussianRegressionFunctionChanged
     procedure :: willEvaluate    => gaussianRegressionWillEvaluate
     procedure :: restore         => gaussianRegressionRestore
     procedure :: separation      => gaussianRegressionSeparation
     procedure :: emulate         => gaussianRegressionEmulate
  end type posteriorSampleLikelihoodGaussianRegression

  type(posteriorSampleLikelihoodGaussianRegression), pointer :: self_
  !$omp threadprivate(self_)

  interface posteriorSampleLikelihoodGaussianRegression
     !!{
     Constructors for the {\normalfont \ttfamily gaussianRegression} posterior sampling likelihood class.
     !!}
     module procedure gaussianRegressionConstructorParameters
     module procedure gaussianRegressionConstructorInternal
  end interface posteriorSampleLikelihoodGaussianRegression

  type :: polynomialIterator
     !!{
     An object used for iterating over coefficients of polynomials.
     !!}
     private
     integer                            :: order       , rank        , &
          &                                orderCurrent, stateCurrent, &
          &                                count
     integer, allocatable, dimension(:) :: indices
   contains
     !![
     <methods>
       <method description="Reset the iterator object to the start of its sequence."                                                                                       method="reset"       />
       <method description="Move to the next iteration of polynomial coefficient indices. Returns true if successful. If no more iterations are available, returns false." method="iterate"     />
       <method description="Return the $i^\mathrm{th}$ index of the polynomial coefficient."                                                                               method="index"       />
       <method description="Return the current order of the polynomial coefficient."                                                                                       method="currentOrder"/>
       <method description="Return an incremental counter (i.e. begins at $0$ and increases by $1$ on each iteration)."                                                    method="counter"     />
     </methods>
     !!]
     procedure :: index        => polynomialIteratorIndex
     procedure :: currentOrder => polynomialIteratorCurrentOrder
     procedure :: counter      => polynomialIteratorCounter
     procedure :: iterate      => polynomialIteratorIterate
     procedure :: reset        => polynomialIteratorReset
  end type polynomialIterator

  interface polynomialIterator
     module procedure polynomialIteratorConstructor
  end interface polynomialIterator

  ! Parameter controlling tolerance used in judging if the emulator has comparable variance to the simulator.
  double precision, parameter :: sigmaSimulatorVariance=3.0d+0

contains

  function gaussianRegressionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily gaussianRegression} posterior sampling likelihood class which builds the object
    from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodGaussianRegression)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (posteriorSampleLikelihoodClass             ), pointer       :: posteriorSampleLikelihood_
    class           (variogramClass                             ), pointer       :: variogram_
    integer                                                                      :: emulatorRebuildCount       , polynomialOrder    , &
         &                                                                          reportCount
    double precision                                                             :: sigmaBuffer                , logLikelihoodBuffer, &
         &                                                                          logLikelihoodErrorTolerance
    logical                                                                      :: emulateOutliers            , dummyEmulator      , &
         &                                                                          assumeZeroVarianceAtZeroLag
    type            (varying_string                             )                :: dumpEmulatorFileRoot

    !![
    <inputParameter>
      <name>emulatorRebuildCount</name>
      <description>The number of steps between rebuilds of the emulator.</description>
      <defaultValue>100</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>polynomialOrder</name>
      <description>The order of the polynomial to fit to the likelihood surface.</description>
      <defaultValue>2</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sigmaBuffer</name>
      <description>The buffer size in units of the likelihood error to use when deciding whether to emulate the likelihood.</description>
      <defaultValue>3.0d0</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>logLikelihoodBuffer</name>
      <description>The buffer size in log-likelihood to use when deciding whether to emulate the likelihood.</description>
      <defaultValue>10.0d0</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>logLikelihoodErrorTolerance</name>
      <description>The tolerance on the likelihood error to accept when deciding whether to emulate the likelihood.</description>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reportCount</name>
      <description>The number of steps between reports of emulator performance.</description>
      <defaultValue>10</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>emulateOutliers</name>
      <description>If true, then outlier chains are always emulated once the simulation is converged.</description>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>assumeZeroVarianceAtZeroLag</name>
      <description>If true, the variogram model is forced to go to zero for states with zero separation (as expected if the likelihood model being emulated is fully deterministic). Otherwise, the variance at zero separation is treated as a free parameter.</description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (parameters%isPresent('dumpEmulatorFileRoot')) then
       !![
       <inputParameter>
	 <name>dumpEmulatorFileRoot</name>
	 <description>The name of a file to which emulator internal state will be dumped. (If empty, no dump occurs.)</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    else
       dumpEmulatorFileRoot=var_str('')
    end if
    !![
    <inputParameter>
      <name>dummyEmulator</name>
      <description>If true, then the emulator is constructed, and performance measured, but likelihoods are always simulated directly.</description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="posteriorSampleLikelihood" name="posteriorSampleLikelihood_" source="parameters"/>
    <objectBuilder class="variogram"                 name="variogram_"                 source="parameters"/>
    !!]
    self=posteriorSampleLikelihoodGaussianRegression(emulatorRebuildCount,polynomialOrder,sigmaBuffer,logLikelihoodBuffer,logLikelihoodErrorTolerance,reportCount,emulateOutliers,assumeZeroVarianceAtZeroLag,char(dumpEmulatorFileRoot),dummyEmulator,posteriorSampleLikelihood_,variogram_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="posteriorSampleLikelihood_"/>
    <objectDestructor name="variogram_"                />
    !!]
    return
  end function gaussianRegressionConstructorParameters

  function gaussianRegressionConstructorInternal(emulatorRebuildCount,polynomialOrder,sigmaBuffer,logLikelihoodBuffer,logLikelihoodErrorTolerance,reportCount,emulateOutliers,assumeZeroVarianceAtZeroLag,dumpEmulatorFileRoot,dummyEmulator,posteriorSampleLikelihood_,variogram_) result(self)
    !!{
    Constructor for {\normalfont \ttfamily gaussianRegression} posterior sampling likelihood class.
    !!}
    implicit none
    type            (posteriorSampleLikelihoodGaussianRegression)                        :: self
    class           (posteriorSampleLikelihoodClass             ), intent(in   ), target :: posteriorSampleLikelihood_
    class           (variogramClass                             ), intent(in   ), target :: variogram_
    integer                                                      , intent(in   )         :: emulatorRebuildCount       , polynomialOrder    , &
         &                                                                                  reportCount
    double precision                                             , intent(in   )         :: sigmaBuffer                , logLikelihoodBuffer, &
         &                                                                                  logLikelihoodErrorTolerance
    logical                                                      , intent(in   )         :: emulateOutliers            , dummyEmulator      , &
         &                                                                                  assumeZeroVarianceAtZeroLag
    character       (len=*                                      ), intent(in   )         :: dumpEmulatorFileRoot
    !![
    <constructorAssign variables="emulatorRebuildCount, polynomialOrder, sigmaBuffer, logLikelihoodBuffer, logLikelihoodErrorTolerance, reportCount, emulateOutliers, assumeZeroVarianceAtZeroLag, dummyEmulator, dumpEmulatorFileRoot, *posteriorSampleLikelihood_, *variogram_"/>
    !!]

    ! Initialize state and counters.
    self%dumpEmulatorCount           =  0
    self%dumpEmulator                =  (self%dumpEmulatorFileRoot /= "")
    self%initialized                 =  .false.
    self%isGood                      =  .false.
    self%accumulatedStateCount       =  0
    self%simulationCount             =  0
    self%evaluationCount             =  0
    self%emulatorCheckCount          =  0
    self%emulatorFailCount           =  0
    return
  end function gaussianRegressionConstructorInternal

  subroutine gaussianRegressionDestructor(self)
    !!{
    Destructor for Gaussian regression likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodGaussianRegression), intent(inout) :: self

    !![
    <objectDestructor name="self%posteriorSampleLikelihood_"/>
    <objectDestructor name="self%variogram_"                />
    !!]
    return
  end subroutine gaussianRegressionDestructor

  double precision function gaussianRegressionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for a Gaussian regression likelihood function.
    !!}
    use :: Dates_and_Times               , only : Formatted_Date_and_Time
    use :: Display                       , only : displayIndent                  , displayMessage                 , displayUnindent, displayVerbosity, &
          &                                       verbosityLevelStandard         , displayMagenta                 , displayReset
    use :: Error_Functions               , only : Error_Function
    use :: Error                         , only : Error_Report
    use :: Linear_Algebra                , only : assignment(=)                  , matrix                         , vector
    use :: MPI_Utilities                 , only : mpiSelf
    use :: Models_Likelihoods_Constants  , only : logImpossible                  , logImprobable
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass      , posteriorSampleStateCorrelation
    use :: String_Handling               , only : operator(//)
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(inout), target           :: self
    class           (posteriorSampleStateClass                  ), intent(inout)                   :: simulationState
    type            (modelParameterList                         ), intent(inout), dimension(:)     :: modelParametersActive_                    , modelParametersInactive_
    class           (posteriorSampleConvergenceClass            ), intent(inout)                   :: simulationConvergence
    double precision                                             , intent(in   )                   :: temperature                               , logLikelihoodCurrent    , &
         &                                                                                            logPriorCurrent                           , logPriorProposed
    real                                                         , intent(inout)                   :: timeEvaluate
    double precision                                             , intent(  out), optional         :: logLikelihoodVariance
    logical                                                      , intent(inout), optional         :: forceAcceptance
    double precision                                             , allocatable  , dimension(:,:,:) :: states
    double precision                                             , allocatable  , dimension(  :,:) :: likelihoods                               , workspace2D
    double precision                                             , allocatable  , dimension(    :) :: likelihoodsCombined                       , workspace                , &
         &                                                                                            likelihoodsFitted                         , separations              , &
         &                                                                                            semiVariances
    integer                                                      , allocatable  , dimension(    :) :: stateCount
    double precision                                             , parameter                       :: emulatorFailureSignificance        =  3.0d+0
    double precision                                             , parameter                       :: emulatorFailureSignificanceAbsolute=  1.0d-3
    double precision                                             , parameter                       :: emulatorFailureLimit               =  5.0d+0
    integer                                                      , parameter                       :: checksTotalMinimum                 =100
    type            (polynomialIterator                         )                                  :: iterator1                                   , iterator2
    type            (vector                                     )                                  :: likelihoodSums                              , coefficients
    type            (matrix                                     )                                  :: stateSums                                   , regressionMatrix
    integer                                                                                        :: accumulatedStateCount                       , i                    , &
         &                                                                                            j                                           , k                    , &
         &                                                                                            evaluationsTotal                            , simulationsTotal     , &
         &                                                                                            stateCountAccept                            , stateCountAcceptChain, &
         &                                                                                            determinantSign                             , checksTotal          , &
         &                                                                                            failuresTotal                               , fileUnit             , &
         &                                                                                            correlationLength                           , l
    double precision                                                                               :: separation                                  , likelihoodError      , &
         &                                                                                            likelihoodEmulated                          , failureRateExpected
    logical                                                                                        :: simulatorWasCalled                          , isDuplicate          , &
         &                                                                                            useLikelihoodEmulated
    character       (len=12                                     )                                  :: label
    type            (varying_string                             )                                  :: message                                     , fileName
    !$GLC attributes unused :: forceAcceptance
    
    ! Report on emulation efficiency.
    if (mod(self%evaluationCount,self%reportCount) == 0 .and. self%evaluationCount > 0 .and. displayVerbosity() >= verbosityLevelStandard) then
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
          call displayMessage(message)
       end if
    end if
    ! Rebuild emulator if necessary.
    allocate(stateCount(0:mpiSelf%count()-1))
    stateCount                           =0
    stateCount           (mpiSelf%rank())=self   %accumulatedStateCount
    stateCount                           =mpiSelf%sum                  (stateCount)
    accumulatedStateCount                =        sum                  (stateCount)
    select type (simulationState)
    class is (posteriorSampleStateCorrelation)
       correlationLength=simulationState%correlationLength()
    class default
       correlationLength=1
    end select    
    if (accumulatedStateCount >= self%emulatorRebuildCount .and. (any(stateCount == self%emulatorRebuildCount) .or. .not.self%initialized) .and. correlationLength > 0) then
       if (mpiSelf%isMaster()) then
          call displayIndent('Rebuilding Gaussian process emulator ['//char(Formatted_Date_and_Time())//']')
          call flush(0)
       end if
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
          if (self%polynomialCoefficientCount > self%emulatorRebuildCount) then
             message='number of points used in emulator ['
             message=message                                              // &
                  &  self%emulatorRebuildCount                            // &
                  &  '] is too few to constrain global trend polynomial ['// &
                  &  self%polynomialCoefficientCount                      // &
                  &  ']'
             call Error_Report(message//{introspection:location})
          end if
          allocate(self%polynomialCoefficient(self%polynomialCoefficientCount                                ))
          allocate(self%likelihoodSums       (self%polynomialCoefficientCount                                ))
          allocate(self%stateSums            (self%polynomialCoefficientCount,self%polynomialCoefficientCount))
          allocate(self%coefficients         (self%polynomialCoefficientCount                                ))
          allocate(self%regressionMatrix     (self%regressionMatrixSize      ,self%regressionMatrixSize      ))
          allocate(self%stateOffset          (self%regressionMatrixSize                                      ))
          allocate(self%weight               (self%regressionMatrixSize                                      ))
          allocate(self%likelihoodResiduals  (self%emulatorRebuildCount                                      ))
          allocate(self%statesCombined       (simulationState%dimension()    ,self%emulatorRebuildCount      ))
          allocate(self%stateScales          (simulationState%dimension()                                    ))
          allocate(self%stateMeans           (simulationState%dimension()                                    ))
          self%initialized=.true.
       end if
       ! Gather the simulator state and likelihood from other chains.
       if (mpiSelf%isMaster()) call displayMessage('gathering states ['//char(Formatted_Date_and_Time())//']')
       allocate(states     (simulationState%dimension(),self%emulatorRebuildCount,0:mpiSelf%count()-1))
       allocate(likelihoods(                            self%emulatorRebuildCount,0:mpiSelf%count()-1))
       states                    =mpiSelf%gather(self%simulationState    )
       likelihoods               =mpiSelf%gather(self%simulatorLikelihood)
       stateCount                =0
       stateCount(mpiSelf%rank())=self   %accumulatedStateCount
       stateCount                =mpiSelf%sum   (stateCount              )
       ! Allocate workspace arrays.
       allocate(likelihoodsCombined( self%emulatorRebuildCount                                 ))
       allocate(likelihoodsFitted  ( self%emulatorRebuildCount                                 ))
       allocate(workspace          ( self%emulatorRebuildCount                                 ))
       allocate(workspace2D        ( self%emulatorRebuildCount, self%polynomialCoefficientCount))
       allocate(separations        ((self%emulatorRebuildCount*(self%emulatorRebuildCount-1))/2))
       allocate(semiVariances      ((self%emulatorRebuildCount*(self%emulatorRebuildCount-1))/2))
       ! Extract likelihoods and states. We may have more states stored than we actually need. The following algorithm combines
       ! states from chains such that we balance the states used across chains as much as possible. That is, we find the minimum
       ! number of states available from all chains and take that number of states from each chain. If more states are needed, we
       ! find the minimum number of remaining states available from all chains and take that number of states from each chain. The
       ! process is repeated until we have accumulated enough states.
       if (mpiSelf%isMaster()) call displayMessage('extracting states ['//char(Formatted_Date_and_Time())//']')
       j=0
       do while (j < self%emulatorRebuildCount)
          stateCountAccept=min(minval(stateCount,mask=stateCount > 0),max(self%emulatorRebuildCount/mpiSelf%count(),1))
          do i=0,mpiSelf%count()-1
             stateCountAcceptChain=min(stateCountAccept,self%emulatorRebuildCount-j,stateCount(i))
             k=stateCount(i)
             l=0
             do while (l < stateCountAcceptChain .and. k > 0)
                j                       =j                 +1
                l                       =l                 +1
                likelihoodsCombined(  j)=likelihoods(  k,i)
                self%statesCombined(:,j)=states     (:,k,i)
                if (k < stateCount(i)) states(:,k:stateCount(i)-1,i)=states(:,k+1:stateCount(i),i)
                stateCount         (  i)=stateCount (    i)-1
                k                       =k                 -correlationLength
             end do
          end do
       end do
       ! Deallocate workspace.
       deallocate(likelihoods)
       deallocate(stateCount )
       deallocate(states     )
       ! Check for duplicated states.
       do i=1,size(self%statesCombined,dim=2)
          do j=1,size(self%statesCombined,dim=2)
             if (i==j) cycle
             if (all(self%statesCombined(:,i) == self%statesCombined(:,j))) then
                message=""
                do k=1,size(self%statesCombined(:,i))
                   write (label,'(e12.6)') self%statesCombined(k,i)
                   message=message//label
                   if (k < size(self%statesCombined(:,i))) message=message//" "
                end do
                call displayIndent  ("duplicated states:")
                call displayMessage (message             )
                call displayUnindent(""                  )
                call Error_Report('duplicated states detected'//{introspection:location})
             end if
          end do
       end do
       ! Evaluate mean state.
       self%stateMeans=sum(self%statesCombined,dim=2)/size(self%statesCombined,dim=2)
       ! Subtract means from states.
       do i=1,size(self%statesCombined,dim=2)
          self%statesCombined(:,i)=self%statesCombined(:,i)-self%stateMeans
       end do
       ! Evaluate sums needed in fitting polynomial trend model.
       if (mpiSelf%isMaster()) call displayMessage('fitting global trends (step #1) ['//char(Formatted_Date_and_Time())//']')
       iterator1=polynomialIterator(self%polynomialOrder,simulationState%dimension())
       iterator2=polynomialIterator(self%polynomialOrder,simulationState%dimension())
       ! Pre-compute products over states.
       do while (iterator1%iterate())
          workspace2D(:,iterator1%counter())=1.0d0
          do j=1,iterator1%currentOrder()
             workspace2D(:,iterator1%counter())=workspace2D(:,iterator1%counter())*self%statesCombined(iterator1%index(j),:)
          end do
       end do
       call iterator1%reset()
       if (mpiSelf%isMaster()) call displayMessage('fitting global trends (step #2) ['//char(Formatted_Date_and_Time())//']')
       ! Compute sums over states.
       do while (iterator1%iterate())
          self%likelihoodSums(iterator1%counter())=sum(workspace2D(:,iterator1%counter())*likelihoodsCombined)
          call iterator2%reset()
          do while (iterator2%iterate())
             self%stateSums(iterator1%counter(),iterator2%counter())=sum(workspace2D(:,iterator1%counter())*workspace2D(:,iterator2%counter()))
          end do
       end do
       if (mpiSelf%isMaster()) call displayMessage('fitting global trends (step #3) ['//char(Formatted_Date_and_Time())//']')
       ! Assign vector and matrix in our linear system.
       likelihoodSums   =self%likelihoodSums
       stateSums        =self%stateSums
       ! Solve for regression coefficients.
       coefficients     =stateSums%linearSystemSolve(likelihoodSums)
       self%coefficients=coefficients
       ! Compute fitted likelihoods.
       if (mpiSelf%isMaster()) call displayMessage('fitting global trends (step #4) ['//char(Formatted_Date_and_Time())//']')
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
       if (mpiSelf%isMaster()) call displayMessage('computing variogram ['//char(Formatted_Date_and_Time())//']')
       k=0
       do i=1,self%emulatorRebuildCount-1
          do j=i+1,self%emulatorRebuildCount
             k=k+1
             separations  (k)=self%separation(self%statesCombined(:,i),self%statesCombined(:,j))
             semiVariances(k)=+0.5d0                         &
                  &           *(                             &
                  &             +self%likelihoodResiduals(i) &
                  &             -self%likelihoodResiduals(j) &
                  &            )**2
          end do
       end do       
       ! Fit the variogram model.
       if (mpiSelf%isMaster()) call displayMessage('fitting variogram ['//char(Formatted_Date_and_Time())//']')
       call self%variogram_%fit(separations,semiVariances)
       ! Compute regression matrix.
       k=0
       if (self%dumpEmulator) then
          self%dumpEmulatorCount=self%dumpEmulatorCount+1
          fileName=self%dumpEmulatorFileRoot//self%dumpEmulatorCount//".log"
          open(newUnit=fileUnit,file=char(fileName),form='formatted',status='unknown')
       end if
       do i=1,self%emulatorRebuildCount
          do j=1,self%emulatorRebuildCount
             separation                =self           %separation (self%statesCombined(:,i),self%statesCombined(:,j))
             self%regressionMatrix(i,j)=self%variogram_%correlation(separation)
             if (self%dumpEmulator) then
                k=k+1
                write (fileUnit,*) i,j,separation,semiVariances(k),self%regressionMatrix(i,j)
             end if
          end do
       end do
       if (self%dumpEmulator) close(fileUnit)
       ! Find the LU decomposition of the regression matrix for later use.
       if (mpiSelf%isMaster()) call displayMessage('computing regression matrix ['//char(Formatted_Date_and_Time())//']')
       self%regressionMatrix(1:self%emulatorRebuildCount  ,  self%emulatorRebuildCount+1)=1.0d0
       self%regressionMatrix(  self%emulatorRebuildCount+1,1:self%emulatorRebuildCount  )=1.0d0
       self%regressionMatrix(  self%emulatorRebuildCount+1,  self%emulatorRebuildCount+1)=0.0d0
       if (allocated(self%regressionMatrixLU)) deallocate(self%regressionMatrixLU)       
       regressionMatrix=matrix(self%regressionMatrix)
       self%regressionMatrixLU        =matrixLU(regressionMatrix)
       determinantSign                =regressionMatrix%signDeterminant()
       self%regressionMatrixIsSingular=(determinantSign == 0)
       if (mpiSelf%isMaster().and.self%regressionMatrixIsSingular) call displayMessage('   ==> regression matrix is singular')
       ! Retain (the most recent) 50% of the required number of points.
       self%simulatorLikelihood  (  1:self%accumulatedStateCount/2)=self%simulatorLikelihood(  self%accumulatedStateCount/2+1:2*(self%accumulatedStateCount/2))
       self%simulationState      (:,1:self%accumulatedStateCount/2)=self%simulationState    (:,self%accumulatedStateCount/2+1:2*(self%accumulatedStateCount/2))
       self%accumulatedStateCount                                  =                           self%accumulatedStateCount/2
       ! Deallocate workspace.
       deallocate(likelihoodsCombined)
       deallocate(likelihoodsFitted  )
       deallocate(workspace          )
       deallocate(workspace2D        )
       deallocate(semiVariances      )
       deallocate(separations        )
       ! Finished.
       if (mpiSelf%isMaster()) call displayUnindent('done ['//char(Formatted_Date_and_Time())//']')
    end if    
    ! Count evaluations.
    self%evaluationCount=self%evaluationCount+1
    ! Ensure arrays are allocated.
    if (.not.allocated(self%simulatorLikelihood)) then
       allocate(self%simulatorLikelihood(                            self%emulatorRebuildCount))
       allocate(self%simulationState    (simulationState%dimension(),self%emulatorRebuildCount))
    end if
    simulatorWasCalled   =.false.
    useLikelihoodEmulated=.false.
    if (self%initialized.and..not.self%regressionMatrixIsSingular) then       
       ! Perform the emulation.
       call self%emulate(simulationState,likelihoodEmulated,likelihoodError)
       if (present(logLikelihoodVariance)) logLikelihoodVariance=likelihoodError**2
       ! Test likelihood emulation. We do this whenever the emulator is currently rated "not good", and periodically otherwise to
       ! monitor emulator behavior. When first initialized, the emulator is rated "not good" such that it has to prove that it is
       ! valid before we actually begin using it.
       if (.not.self%isGood .or. (mod(self%evaluationCount,self%reportCount) == 0 .and. self%evaluationCount > 0)) then
          ! Every so many steps we evaluate the simulated likelihood and check that our emulator is reliable.
          if (likelihoodError <= sqrt(4.0d0*self%variogram_%variogram(0.0d0))*sigmaSimulatorVariance) then
             ! Evaluate the simulator for synchronization purposes only.
             call evaluateSimulator(synchronizeOnly=.true.,logLikelihood_=gaussianRegressionEvaluate,logLikelihoodVariance_=logLikelihoodVariance)
             ! Emulator variance is comparable to that of the simulator (which is the best we can do), so consider this to be a successful "check".
             self%emulatorCheckCount=self%emulatorCheckCount+1
             useLikelihoodEmulated  =.true.
          else
             ! Evaluate the simulator for comparison to our emulator.
             call evaluateSimulator(synchronizeOnly=.false.,logLikelihood_=gaussianRegressionEvaluate,logLikelihoodVariance_=logLikelihoodVariance)
             ! Check that a non-impossible likelihood was returned.
             if (gaussianRegressionEvaluate > logImpossible) then
                ! Count number of emulator checks and the failure rate. We include a final condition here that accounts for the
                ! variance in the simulator itself - we can't expect the emulator to perform better than the simulator.
                self%emulatorCheckCount=self%emulatorCheckCount+1
                if     (                                                                                                          &
                     &   abs(likelihoodEmulated-gaussianRegressionEvaluate) > emulatorFailureSignificance        *likelihoodError &
                     &  .and.                                                                                                     &
                     &   abs(likelihoodEmulated-gaussianRegressionEvaluate) > emulatorFailureSignificanceAbsolute                 &
                     & ) self%emulatorFailCount=self%emulatorFailCount+1            
             end if
          end if          
          checksTotal  =mpiSelf%sum(self%emulatorCheckCount)
          failuresTotal=mpiSelf%sum(self%emulatorFailCount )
          ! Determine if the emulator is sufficiently good to use or not.
          failureRateExpected=1.0d0-Error_Function(emulatorFailureSignificance/sqrt(2.0d0))
          self%isGood= checksTotal         > checksTotalMinimum                                         &
               &      .and.                                                                             &
               &       dble(failuresTotal) < emulatorFailureLimit*failureRateExpected*dble(checksTotal) &
               &      .and.                                                                             &
               &       .not.self%dummyEmulator
          ! Report on emulator failure rate.
          if (mpiSelf%isMaster() .and. displayVerbosity() >= verbosityLevelStandard .and. mod(self%evaluationCount,self%reportCount) == 0 .and. self%evaluationCount > 0) then
             write (label,'(i8)') failuresTotal
             message='Emulator failed '//trim(adjustl(label))//' times out of '
             write (label,'(i8)') checksTotal
             message=message//trim(adjustl(label))//' checks ('
             write (label,'(f6.2)') 100.0d0*dble(failuresTotal)/dble(checksTotal)
             message=message//trim(adjustl(label))//'% - expect '
             write (label,'(f6.2)') 100.0d0*failureRateExpected
             message=message//trim(adjustl(label))//'% for perfect emulator)'
             call displayMessage(message)
             if (.not.self%isGood) call displayMessage(displayMagenta()//'WARNING:'//displayReset()//' emulator failure rate is too high - emulator will not be used')
          end if
       else
          ! Apply logical to decide if we can use the emulated likelihood. If we can, simply return that value. Otherwise, we fall through to the simulation step below.
          !! If the simulation is converged and the chain is an outlier, we may be able to accept emulation always.
          if (simulationConvergence%isConverged().and.simulationConvergence%stateIsOutlier(simulationState%chainIndex()).and.self%emulateOutliers                                                              ) useLikelihoodEmulated=.true.
          !! If the emulated posterior is sufficiently smaller (even accounting for uncertainties) than the current posterior than
          !! we can use the emulated likelihood - the precise likelihood of such unlikely states should not mattter.
          if (likelihoodEmulated+logPriorProposed+self%sigmaBuffer*likelihoodError < logLikelihoodCurrent+logPriorCurrent-          self            %logLikelihoodBuffer                *temperature           ) useLikelihoodEmulated=.true.
          !! If the uncertainty in the emulated likelihood is below a threshold, accept it.
          if (                                                     likelihoodError <                                                self            %logLikelihoodErrorTolerance        *temperature           ) useLikelihoodEmulated=.true.
          !! If the uncertainty in the emulated likelihood is comparable to the variance in the simulator itself, we may as well
          !! just use the emulated likelihood. The variogram returns the semi-variance of the simulator, and we care about
          !! differences between the emulator and simulator, so we scale the simulator variance by a factor of 4.
          if (                                                     likelihoodError <                                      sqrt(4.0d0*self%variogram_%variogram                  (0.0d0))*sigmaSimulatorVariance) useLikelihoodEmulated=.true.
       end if
    end if
    ! Evaluate the simulator now if we have not already done so.
    if (.not.simulatorWasCalled) call evaluateSimulator(synchronizeOnly=useLikelihoodEmulated,logLikelihood_=gaussianRegressionEvaluate,logLikelihoodVariance_=logLikelihoodVariance)
    ! Check if the emulated likelihood is to be used.
    if (useLikelihoodEmulated) then
       ! Return the emulated likelihood.
       gaussianRegressionEvaluate=likelihoodEmulated
       logLikelihoodVariance     =likelihoodError   **2
       return
    else
       ! Increment the count of simulation evaluations.
       self%simulationCount=self%simulationCount+1
    end if
    return
    
  contains
    
    subroutine evaluateSimulator(synchronizeOnly,logLikelihood_,logLikelihoodVariance_)
      !!{
      Call the {\normalfont \ttfamily evaluate} method of the simulator. If {\normalfont \ttfamily synchronizeOnly} is true then
      no evaluation is actually needed, but we must still call the method to allow for possible MPI synchronization. In this case
      call with an impossible proposed prior such that the simulator can choose to not evaluate.
      !!}
      implicit none
      logical         , intent(in   ) :: synchronizeOnly
      double precision, intent(  out) :: logLikelihood_ , logLikelihoodVariance_

      if (synchronizeOnly) then
         logLikelihood_=self%posteriorSampleLikelihood_%evaluate(simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logImpossible   ,timeEvaluate,logLikelihoodVariance_)
      else
         logLikelihood_=self%posteriorSampleLikelihood_%evaluate(simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance_)
         ! Store the likelihood and state.
         if (logLikelihood_ > logImprobable .and. self%accumulatedStateCount < self%emulatorRebuildCount) then
            isDuplicate=.false.
            if (self%accumulatedStateCount > 0) then
               do i=1,self%accumulatedStateCount
                  isDuplicate=all(self%simulationState(:,i) == simulationState%get())
                  if (isDuplicate) exit
               end do
            end if
            if (.not.isDuplicate) then
               self%accumulatedStateCount                            =self%accumulatedStateCount+1
               self%simulatorLikelihood(  self%accumulatedStateCount)=logLikelihood_
               self%simulationState    (:,self%accumulatedStateCount)=simulationState%get()
            end if
         end if
      end if
      simulatorWasCalled=.true.
      return
    end subroutine evaluateSimulator
    
  end function gaussianRegressionEvaluate

  double precision function gaussianRegressionSeparation(self,state1,state2)
    !!{
    Determine the separation between two state vectors.
    !!}
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression)              , intent(in   ) :: self
    double precision                                             , dimension(:), intent(in   ) :: state1, state2

    gaussianRegressionSeparation=+sqrt(                         &
         &                             +sum(                    &
         &                                  +(                  &
         &                                    +(                &
         &                                      +state1         &
         &                                      -state2         &
         &                                     )                &
         &                                    /self%stateScales &
         &                                   )**2               &
         &                                 )                    &
         &                            )
    return
  end function gaussianRegressionSeparation

  subroutine gaussianRegressionFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodGaussianRegression), intent(inout) :: self

    ! Reset emulator state.
    self%initialized          =.false.
    self%accumulatedStateCount=0
    self%simulationCount      =0
    self%evaluationCount      =0
    self%emulatorCheckCount   =0
    self%emulatorFailCount    =0
    ! Free allocated space.
    if (allocated(self%polynomialCoefficient)) deallocate(self%polynomialCoefficient)
    if (allocated(self%likelihoodSums       )) deallocate(self%likelihoodSums       )
    if (allocated(self%stateSums            )) deallocate(self%stateSums            )
    if (allocated(self%coefficients         )) deallocate(self%coefficients         )
    if (allocated(self%regressionMatrix     )) deallocate(self%regressionMatrix     )
    if (allocated(self%stateOffset          )) deallocate(self%stateOffset          )
    if (allocated(self%weight               )) deallocate(self%weight               )
    if (allocated(self%likelihoodResiduals  )) deallocate(self%likelihoodResiduals  )
    if (allocated(self%statesCombined       )) deallocate(self%statesCombined       )
    if (allocated(self%stateScales          )) deallocate(self%stateScales          )
    if (allocated(self%stateMeans           )) deallocate(self%stateMeans           )
    ! Let the simulator know that the likelihood function may have changed.
    call self%posteriorSampleLikelihood_%functionChanged()
    return
  end subroutine gaussianRegressionFunctionChanged

  logical function gaussianRegressionWillEvaluate(self,simulationState,modelParameters_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
    !!{
    Return true if the log-likelihood will be evaluated.
    !!}
    use :: Models_Likelihoods_Constants  , only : logImpossible
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(inout)               :: self
    class           (posteriorSampleStateClass                  ), intent(inout)               :: simulationState
    type            (modelParameterList                         ), intent(in   ), dimension(:) :: modelParameters_
    class           (posteriorSampleConvergenceClass            ), intent(inout)               :: simulationConvergence
    double precision                                             , intent(in   )               :: temperature          , logLikelihoodCurrent, &
         &                                                                                        logPriorCurrent      , logPriorProposed
    double precision                                                                           :: likelihoodEmulated   , likelihoodEmulatedError
    !$GLC attributes unused :: modelParameters_

    if (logPriorProposed <= logImpossible) then
       ! Prior is impossible, no need to evaluate.
       gaussianRegressionWillEvaluate=.false.
       return
    else
       gaussianRegressionWillEvaluate=.true.
    end if
    if (self%initialized.and..not.self%regressionMatrixIsSingular) then
       if (mod(self%evaluationCount,self%reportCount) == 0 .and. self%evaluationCount > 0) then
          ! Emulation will be tested, so the simulator is always run.
          gaussianRegressionWillEvaluate=.true.
          return
       else
          ! Evaluate the emulator.
          call self%emulate(simulationState,likelihoodEmulated,likelihoodEmulatedError)
          ! If simulation is converged, this is an outlier chain, and we're told to emulate such outliers, then return.
          if (simulationConvergence%isConverged().and.simulationConvergence%stateIsOutlier(simulationState%chainIndex()).and.self%emulateOutliers) gaussianRegressionWillEvaluate=.false.
          ! If likelihood is well below current likelihood (and this should be sufficient that changes in priors won't affect the
          ! conclusion), then use the emulated likelihood. We scale the likelihood buffer by the current temperature to reflect the
          ! fact that transitions between states are easier when the temperature is high.
          if (likelihoodEmulated+logPriorProposed+self%sigmaBuffer*likelihoodEmulatedError < logLikelihoodCurrent+logPriorCurrent-           self          %logLikelihoodBuffer                *temperature            ) gaussianRegressionWillEvaluate=.false.
          ! Return if the error is below the tolerance. We increase the tolerance value in proportion to temperature since the
          ! likelihoods will be divided by this amount when evaluating transition probabilities.
          if (likelihoodEmulatedError                                                      <                                                 self          %logLikelihoodErrorTolerance        *temperature            ) gaussianRegressionWillEvaluate=.false.
          ! Return if the uncertainty in the emulated likelihood is comparable to the variance in the simulator itself, we may as
          ! well just use the emulated likelihood. The variogram returns the semi-variance of the simulator, and we care about
          ! differences between the emulator and simulator, so we scale the simulator variance by a factor of 4.
          if (likelihoodEmulatedError                                                      <                                      sqrt(4.0d0*self%variogram_%variogram                  (0.0d0))*sigmaSimulatorVariance) gaussianRegressionWillEvaluate=.false.
       end if
       return
    else
       ! Evaluate the simulator.
       gaussianRegressionWillEvaluate=.true.
    end if
    return
  end function gaussianRegressionWillEvaluate

  subroutine gaussianRegressionEmulate(self,simulationState,likelihoodEmulated,likelihoodEmulatedError)
    !!{
    Evaluate the model emulator.
    !!}
    use :: Linear_Algebra, only : assignment(=), vector
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(inout)               :: self
    class           (posteriorSampleStateClass                  ), intent(inout)               :: simulationState
    double precision                                             , intent(  out)               :: likelihoodEmulated, likelihoodEmulatedError
    double precision                                             , allocatable  , dimension(:) :: stateCurrent
    integer                                                                                    :: i                 , j
    double precision                                                                           :: separation        , likelihoodFit
    type            (vector                                     )                              :: stateOffsetVector , weightVector
    type            (polynomialIterator                         )                              :: iterator1

    ! Compute vector D.
    allocate(stateCurrent(simulationState%dimension()))
    stateCurrent=simulationState%get()-self%stateMeans
    do i=1,self%emulatorRebuildCount
       separation         =self           %separation (stateCurrent,self%statesCombined(:,i))
       self%stateOffset(i)=self%variogram_%correlation(separation                           )
    end do
    self%stateOffset(self%emulatorRebuildCount+1)=1.0d0
    ! Solve the linear system.
    stateOffsetVector=vector(self%stateOffset)
    weightVector     =self%regressionMatrixLU%squareSystemSolve(stateOffsetVector)
    self%weight      =weightVector
    ! Compute the likelihood and variance.
    likelihoodEmulated     =sum(self%likelihoodResiduals*self%weight(1:self%emulatorRebuildCount))
    likelihoodEmulatedError=self%variogram_%variogram()-sum(self%weight*self%stateOffset)
    if (likelihoodEmulatedError >= 0.0d0) then
       likelihoodEmulatedError=sqrt(likelihoodEmulatedError)
    else
       likelihoodEmulatedError=0.0d0
    end if
    iterator1=polynomialIterator(self%polynomialOrder,simulationState%dimension())
    do while (iterator1%iterate())
       likelihoodFit=self%coefficients(iterator1%counter())
       do j=1,iterator1%currentOrder()
          likelihoodFit=likelihoodFit*stateCurrent(iterator1%index(j))
       end do
       likelihoodEmulated=likelihoodEmulated+likelihoodFit
    end do
    return
  end subroutine gaussianRegressionEmulate

  subroutine gaussianRegressionRestore(self,simulationState,logLikelihood)
    !!{
    Process a previous state to restore likelihood function.
    !!}
    use :: Models_Likelihoods_Constants, only : logImpossible
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(inout)               :: self
    double precision                                             , intent(in   ), dimension(:) :: simulationState
    double precision                                             , intent(in   )               :: logLikelihood
    logical                                                                                    :: storeState

    if (logLikelihood > logImpossible) then
       if (.not.allocated(self%simulatorLikelihood)) then
          allocate(self%simulatorLikelihood(self%emulatorRebuildCount))
          allocate(self%simulationState    (size(simulationState),self%emulatorRebuildCount))
       end if
       if (self%accumulatedStateCount == self%emulatorRebuildCount) then
          ! Discard the oldest state.
          self%simulatorLikelihood  (  1:self%emulatorRebuildCount-1)=self%simulatorLikelihood  (  2:self%emulatorRebuildCount)
          self%simulationState      (:,1:self%emulatorRebuildCount-1)=self%simulationState      (:,2:self%emulatorRebuildCount)
          self%accumulatedStateCount                                 =self%accumulatedStateCount                               -1
       end if
       ! Store the state unless it is identical to the previous state.
       storeState=self%accumulatedStateCount == 0
       if (.not.storeState) storeState=any(simulationState /= self%simulationState(:,self%accumulatedStateCount))
       if (storeState) then
          self%accumulatedStateCount=min(self%accumulatedStateCount+1,self%emulatorRebuildCount)
          self%simulatorLikelihood(  self%accumulatedStateCount)=logLikelihood
          self%simulationState    (:,self%accumulatedStateCount)=simulationState
       end if
    end if
    return
  end subroutine gaussianRegressionRestore

  function polynomialIteratorConstructor(order,rank) result(self)
    !!{
    Create a polynomial iterator for a polynomial of specified {\normalfont \ttfamily order} and {\normalfont \ttfamily rank}.
    !!}
    implicit none
    type   (polynomialIterator)                :: self
    integer                    , intent(in   ) :: order, rank
    !![
    <constructorAssign variables="order, rank"/>
    !!]
    
    allocate(self%indices(order))
    call self%reset()
    return
  end function polynomialIteratorConstructor
  
  integer function polynomialCoefficientCount(n,d)
    !!{
    Return the number of coefficients at {\normalfont \ttfamily n}$^\mathrm{th}$ order in polynomial of dimension {\normalfont \ttfamily d}.
    !!}
    use :: Factorials, only : Factorial
    implicit none
    integer, intent(in   ) :: n, d

    if (n == 2) then
       polynomialCoefficientCount=(d*(1+d))/2
    else
       polynomialCoefficientCount=nint(Factorial(d+n-1)/Factorial(n)/Factorial(d-1))
    end if
    return
  end function polynomialCoefficientCount

  subroutine polynomialIteratorReset(self)
    !!{
    Reset a polynomial iterator.
    !!}
    implicit none
    class(polynomialIterator), intent(inout) :: self

    self%orderCurrent=-1
    self%stateCurrent=-1
    self%count       = 0
    self%indices     = 0
    return
  end subroutine polynomialIteratorReset

  logical function polynomialIteratorIterate(self)
    !!{
    Iterate over polynomial coefficients.
    !!}
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
    !!{
    Return the requested index of a polynomial iterator.
    !!}
    implicit none
    class  (polynomialIterator), intent(inout) :: self
    integer                    , intent(in   ) :: i

    polynomialIteratorIndex=self%indices(i)
    return
  end function polynomialIteratorIndex

  integer function polynomialIteratorCurrentOrder(self)
    !!{
    Return the current order of a polynomial iterator.
    !!}
    implicit none
    class  (polynomialIterator), intent(inout) :: self

    polynomialIteratorCurrentOrder=self%orderCurrent
    return
  end function polynomialIteratorCurrentOrder

  integer function polynomialIteratorCounter(self)
    !!{
    Return the current count of a polynomial iterator.
    !!}
    implicit none
    class  (polynomialIterator), intent(inout) :: self

    polynomialIteratorCounter=self%count
    return
  end function polynomialIteratorCounter
