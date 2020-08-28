!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implementation of a posterior sampling likelihood class which implements a likelihood using Gaussian regression to emulate
  !% another likelihood.

  use            :: Linear_Algebra, only : matrixLU
  use, intrinsic :: ISO_C_Binding , only : c_size_t

  !# <posteriorSampleLikelihood name="posteriorSampleLikelihoodGaussianRegression">
  !#  <description>
  !#   The likelihood is computed either using another likelihood function (the ``simulator'') or via Gaussian regression emulation of
  !#   that simulator. The details of the emulation algorithm are specified by the following sub-parameters:
  !#   \begin{description}
  !#   \item[{\normalfont \ttfamily emulatorRebuildCount}] The number of simulator evaluations from which the emulator is built;
  !#   \item[{\normalfont \ttfamily polynomialOrder}] The order of the polynomial fitted to the simulator likelihoods prior to Gaussian regression;
  !#   \item[{\normalfont \ttfamily sigmaBuffer}] See below;
  !#   \item[{\normalfont \ttfamily logLikelihoodBuffer}] See below;
  !#   \item[{\normalfont \ttfamily logLikelihoodErrorTolerance}] See below;
  !#   \item[{\normalfont \ttfamily reportCount}] The number of likelihood evaluations between successive reports on the status of the emulator;
  !#   \item[{\normalfont \ttfamily emulateOutliers}] If true, then outlier chains are always emulated post-convergence (this is safe if such chains are not
  !#    used in constructing proposals for non-outlier chains);
  !#   \item[{\normalfont \ttfamily simulatorLikelihood}] Contains another likelihood function definition which will be used to construct the simulator.
  !#   \end{description}
  !#   
  !#   In detail, this likelihood function first collects {\normalfont \ttfamily emulatorRebuildCount} likelihood evaluations from the
  !#   simulator. It then fits a polynomial of order {\normalfont \ttfamily polynomialOrder} and of dimension equal to the dimension of
  !#   the state vector to the simulated likelihoods. Gaussian regression is performed on the residuals of the simulated likelihoods
  !#   after this polynomial fit is removed. Once the emulator has been built in this way every second simulated state is discarded, and
  !#   accumulation of new simulated states continues. Once {\normalfont \ttfamily emulatorRebuildCount} simulated states have once again
  !#   been accumulated a new simulator is built. This ensures that the emulator does not lose all information used in building the
  !#   previous emulator\footnote{This would be unfortunate as the second emulator to be built would then contain information on only
  !#     those regions of the state space that were poorly emulated before.}, instead information from older emulators decays
  !#   exponentially.
  !#   
  !#   Once an emulator has been built, on each successive likelihood evaluation the emulated log-likelihood $\log\mathcal{L}_\mathrm{e}$
  !#   and its error estimate $\sigma_{\log\mathcal{L}_\mathrm{e}}$ are computed. The emulated likelihood is then returned if:
  !#   \begin{equation}
  !#   \log\mathcal{P}^\prime + \log\mathcal{L}_\mathrm{e} + N \sigma_{\log\mathcal{L}_\mathrm{e}} &lt; \log\mathcal{P} + \log\mathcal{L} - T \Delta\log\mathcal{L},
  !#   \end{equation}
  !#   where $N=${\normalfont \ttfamily sigmaBuffer}, $\Delta\log\mathcal{L}=${\normalfont \ttfamily logLikelihoodBuffer}, $T$ is the
  !#   temperature, $\log\mathcal{L}$ is the current log-likelhood, $\log\mathcal{P}$ is the current log-prior probability, and
  !#   $\log\mathcal{P}^\prime$ is the proposed log-prior probability, or if
  !#   \begin{equation}
  !#   \sigma_{\log\mathcal{L}_\mathrm{e}} &lt; T \sigma_{\log\mathcal{L}},
  !#   \end{equation}
  !#   where $\sigma_{\log\mathcal{L}}=${\normalfont \ttfamily logLikelihoodErrorTolerance}, otherwise the simulator is used to compute
  !#   the exact likelihood. In this way, the emulated likelihood is used if it is sufficiently below the current likelihood that, even
  !#   accounting for the emulation error, transition to the new state is highly unlikely, or if the error on the likelihood emulation is
  !#   sufficiently small that it will not have a significant effect on the transition probabilty to the proposed state.
  !#   
  !#   If verbosity is set to 2 or greater than a report will be issued every {\normalfont \ttfamily reportCount} evaluations. The report
  !#   will give the proportions of simulated vs. emulated evaluations. Additionally, during the evaluation where the report is issued,
  !#   both the emulated and simulated log-likelihoods are evaluated and are tested to see if they lie within
  !#   $3 \sigma_{\log\mathcal{L}_\mathrm{e}}$ of each other. The rate of failures (i.e. where the two differ by more than this amount)
  !#   is then reported.
  !#  </description>
  !# </posteriorSampleLikelihood>
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodGaussianRegression
     !% Implementation of a posterior sampling likelihood class which implements a likelihood using Gaussian regression to emulate
     !% another likelihood.
     private
     class           (posteriorSampleLikelihoodClass), pointer                     :: posteriorSampleLikelihood_ => null()
     integer                                                                       :: accumulatedStateCount     , emulatorRebuildCount       , &
          &                                                                           polynomialOrder           , polynomialCoefficientCount , &
          &                                                                           reportCount               , simulationCount            , &
          &                                                                           evaluationCount           , emulatorCheckCount         , &
          &                                                                           emulatorFailCount         , dumpEmulatorCount
     integer         (c_size_t                      )                              :: regressionMatrixSize
     double precision                                , allocatable, dimension(:  ) :: simulatorLikelihood       , polynomialCoefficient      , &
          &                                                                           likelihoodSums            , coefficients               , &
          &                                                                           stateOffset               , weight                     , &
          &                                                                           likelihoodResiduals       , stateScales                , &
          &                                                                           stateMeans
     double precision                                , allocatable, dimension(:,:) :: simulationState           , stateSums                  , &
          &                                                                           regressionMatrix          , statesCombined
     logical                                                                       :: initialized               , regressionMatrixIsSingular , &
          &                                                                           isGood
     type            (matrixLU                      ), allocatable                 :: regressionMatrixLU
     double precision                                                              :: C0                        , C1                         , &
          &                                                                           CR                        , sigmaBuffer                , &
          &                                                                           logLikelihoodBuffer       , logLikelihoodErrorTolerance
     logical                                                                       :: emulateOutliers           , dumpEmulator               , &
          &                                                                           dummyEmulator
     type            (varying_string                 )                             :: dumpEmulatorFileRoot
     ! Workspaces used when fitting the semi-variogram.
     double precision                                 , allocatable, dimension(:) :: separationsNormalized      , semiVariancesNormalized    , &
          &                                                                          separationsBinned          , semiVariancesBinned        , &
          &                                                                          separationsLimited
     double precision                                                             :: separationNormalization    , semiVarianceNormalization
     integer                                                                      :: binCount
   contains
     final     ::                    gaussianRegressionDestructor
     procedure :: evaluate        => gaussianRegressionEvaluate
     procedure :: functionChanged => gaussianRegressionFunctionChanged
     procedure :: willEvaluate    => gaussianRegressionWillEvaluate
     procedure :: restore         => gaussianRegressionRestore
  end type posteriorSampleLikelihoodGaussianRegression

  type(posteriorSampleLikelihoodGaussianRegression), pointer :: gaussianRegressionSelf
  !$omp threadprivate(gaussianRegressionSelf)

  interface posteriorSampleLikelihoodGaussianRegression
     !% Constructors for the {\normalfont \ttfamily gaussianRegression} posterior sampling likelihood class.
     module procedure gaussianRegressionConstructorParameters
     module procedure gaussianRegressionConstructorInternal
  end interface posteriorSampleLikelihoodGaussianRegression

  type polynomialIterator
     !% An object used for iterating over coefficients of polynomials.
     private
     integer                            :: order       , rank        , &
          &                                orderCurrent, stateCurrent, &
          &                                count
     integer, allocatable, dimension(:) :: indices
   contains
     !@ <objectMethods>
     !@   <object>polynomialIterator</object>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Reset the iterator object to the start of its sequence.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>iterate</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Move to the next iteration of polynomial coefficient indices. Returns true if successful. If no more iterations are available, returns false.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>index</method>
     !@     <type>\intone</type>
     !@     <arguments>\intone\ i\argin</arguments>
     !@     <description>Return the $i^\mathrm{th}$ index of the polynomial coefficient.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>currentOrder</method>
     !@     <type>\intone</type>
     !@     <arguments></arguments>
     !@     <description>Return the current order of the polynomial coefficient.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>counter</method>
     !@     <type>\intone</type>
     !@     <arguments></arguments>
     !@     <description>Return an incremental counter (i.e. begins at $0$ and increases by $1$ on each iteration).</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: index        => polynomialIteratorIndex
     procedure :: currentOrder => polynomialIteratorCurrentOrder
     procedure :: counter      => polynomialIteratorCounter
     procedure :: iterate      => polynomialIteratorIterate
     procedure :: reset        => polynomialIteratorReset
  end type polynomialIterator

  interface polynomialIterator
     module procedure polynomialIteratorConstructor
  end interface polynomialIterator

contains

  function gaussianRegressionConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily gaussianRegression} posterior sampling likelihood class which builds the object
    !% from a parameter set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodGaussianRegression)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (posteriorSampleLikelihoodClass             ), pointer       :: posteriorSampleLikelihood_
    integer                                                                      :: emulatorRebuildCount       , polynomialOrder    , &
         &                                                                          reportCount
    double precision                                                             :: sigmaBuffer                , logLikelihoodBuffer, &
         &                                                                          logLikelihoodErrorTolerance
    logical                                                                      :: emulateOutliers            , dummyEmulator
    type            (varying_string                             )                :: dumpEmulatorFileRoot

    !# <inputParameter>
    !#   <name>emulatorRebuildCount</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The number of steps between rebuilds of the emulator.</description>
    !#   <defaultValue>100</defaultValue>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>polynomialOrder</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The order of the polynomial to fit to the likelihood surface.</description>
    !#   <defaultValue>2</defaultValue>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigmaBuffer</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The buffer size in units of the likelihood error to use when deciding whether to emulate the likelihood.</description>
    !#   <defaultValue>3.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>logLikelihoodBuffer</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The buffer size in log-likelihood to use when deciding whether to emulate the likelihood.</description>
    !#   <defaultValue>10.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>logLikelihoodErrorTolerance</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The tolerance on the likelihood error to accept when deciding whether to emulate the likelihood.</description>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>reportCount</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The number of steps between reports of emulator performance.</description>
    !#   <defaultValue>10</defaultValue>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>emulateOutliers</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, then outlier chains are always emulated once the simulation is converged.</description>
    !#   <defaultValue>.true.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>dumpEmulatorFileRoot</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of a file to which emulator internal state will be dumped. (If empty, no dump occurs.)</description>
    !#   <defaultValue>var_str('')</defaultValue>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>dummyEmulator</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, then the emulator is constructed, and performance measured, but likelihoods are always simulated directly.</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
     !# <objectBuilder class="posteriorSampleLikelihood" name="posteriorSampleLikelihood_" source="parameters"/>
    self=posteriorSampleLikelihoodGaussianRegression(emulatorRebuildCount,polynomialOrder,sigmaBuffer,logLikelihoodBuffer,logLikelihoodErrorTolerance,reportCount,emulateOutliers,char(dumpEmulatorFileRoot),dummyEmulator,posteriorSampleLikelihood_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="posteriorSampleLikelihood_"/>
    return
  end function gaussianRegressionConstructorParameters

  function gaussianRegressionConstructorInternal(emulatorRebuildCount,polynomialOrder,sigmaBuffer,logLikelihoodBuffer,logLikelihoodErrorTolerance,reportCount,emulateOutliers,dumpEmulatorFileRoot,dummyEmulator,posteriorSampleLikelihood_) result(self)
    !% Constructor for ``gaussianRegression'' posterior sampling likelihood class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (posteriorSampleLikelihoodGaussianRegression)                        :: self
    class           (posteriorSampleLikelihoodClass             ), intent(in   ), target :: posteriorSampleLikelihood_
    integer                                                      , intent(in   )         :: emulatorRebuildCount       , polynomialOrder    , &
         &                                                                                  reportCount
    double precision                                             , intent(in   )         :: sigmaBuffer                , logLikelihoodBuffer, &
         &                                                                                  logLikelihoodErrorTolerance
    logical                                                      , intent(in   )         :: emulateOutliers            , dummyEmulator
    character       (len=*                                      ), intent(in   )         :: dumpEmulatorFileRoot
    !# <constructorAssign variables="emulatorRebuildCount, polynomialOrder, sigmaBuffer, logLikelihoodBuffer, logLikelihoodErrorTolerance, reportCount, emulateOutliers, dummyEmulator, dumpEmulatorFileRoot, *posteriorSampleLikelihood_"/>

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
    !% Destructor for Gaussian regression likelihood class.
    implicit none
    type(posteriorSampleLikelihoodGaussianRegression), intent(inout) :: self

    !# <objectDestructor name="self%posteriorSampleLikelihood_"/>
    return
  end subroutine gaussianRegressionDestructor

  double precision function gaussianRegressionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !% Return the log-likelihood for a Gaussian regression likelihood function.
    use :: Dates_and_Times               , only : Formatted_Date_and_Time
    use :: Error_Functions               , only : Error_Function
    use :: Galacticus_Display            , only : Galacticus_Display_Indent      , Galacticus_Display_Message     , Galacticus_Display_Unindent, Galacticus_Verbosity_Level, &
          &                                       verbosityInfo
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Linear_Algebra                , only : vector                         , matrix                         , assignment(=)
    use :: MPI_Utilities                 , only : mpiSelf
    use :: Memory_Management             , only : allocateArray
    use :: Models_Likelihoods_Constants  , only : logImpossible
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass      , posteriorSampleStateCorrelation
    use :: String_Handling               , only : operator(//)
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(inout)                   :: self
    class           (posteriorSampleStateClass                  ), intent(inout)                   :: simulationState
    type            (modelParameterList                         ), intent(in   ), dimension(:)     :: modelParametersActive_                    , modelParametersInactive_
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
    logical                                                                                        :: likelihoodIsSimulated
    character       (len=8                                      )                                  :: label
    type            (varying_string                             )                                  :: message                            , fileName
    !$GLC attributes unused :: forceAcceptance

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
          call Galacticus_Display_Indent('Rebuilding Gaussian process emulator ['//char(Formatted_Date_and_Time())//']')
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
             call Galacticus_Error_Report(message//{introspection:location})
          end if
          call allocateArray(self%polynomialCoefficient,    [self%polynomialCoefficientCount                                ] )
          call allocateArray(self%likelihoodSums       ,    [self%polynomialCoefficientCount                                ] )
          call allocateArray(self%stateSums            ,    [self%polynomialCoefficientCount,self%polynomialCoefficientCount] )
          call allocateArray(self%coefficients         ,    [self%polynomialCoefficientCount                                ] )
          call allocateArray(self%regressionMatrix     ,int([self%regressionMatrixSize      ,self%regressionMatrixSize      ]))
          call allocateArray(self%stateOffset          ,int([self%regressionMatrixSize                                      ]))
          call allocateArray(self%weight               ,int([self%regressionMatrixSize                                      ]))
          call allocateArray(self%likelihoodResiduals  ,    [self%emulatorRebuildCount                                      ] )
          call allocateArray(self%statesCombined       ,    [simulationState%dimension()    ,self%emulatorRebuildCount      ] )
          call allocateArray(self%stateScales          ,    [simulationState%dimension()                                    ] )
          call allocateArray(self%stateMeans           ,    [simulationState%dimension()                                    ] )
          self%initialized=.true.
       end if
       ! Gather the simulator state and likelihood from other chains.
       if (mpiSelf%isMaster()) call Galacticus_Display_Message('gathering states ['//char(Formatted_Date_and_Time())//']')
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
       if (mpiSelf%isMaster()) call Galacticus_Display_Message('extracting states ['//char(Formatted_Date_and_Time())//']')
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
       ! Evaluate mean state.
       self%stateMeans=sum(self%statesCombined,dim=2)/size(self%statesCombined,dim=2)
       ! Subtract means from states.
       do i=1,size(self%statesCombined,dim=2)
          self%statesCombined(:,i)=self%statesCombined(:,i)-self%stateMeans
       end do
       ! Evaluate sums needed in fitting polynomial trend model.
       if (mpiSelf%isMaster()) call Galacticus_Display_Message('fitting global trends (step #1) ['//char(Formatted_Date_and_Time())//']')
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
       if (mpiSelf%isMaster()) call Galacticus_Display_Message('fitting global trends (step #2) ['//char(Formatted_Date_and_Time())//']')
       ! Compute sums over states.
       do while (iterator1%iterate())
          self%likelihoodSums(iterator1%counter())=sum(workspace2D(:,iterator1%counter())*likelihoodsCombined)
          call iterator2%reset()
          do while (iterator2%iterate())
             self%stateSums(iterator1%counter(),iterator2%counter())=sum(workspace2D(:,iterator1%counter())*workspace2D(:,iterator2%counter()))
          end do
       end do
       if (mpiSelf%isMaster()) call Galacticus_Display_Message('fitting global trends (step #3) ['//char(Formatted_Date_and_Time())//']')
       ! Assign vector and matrix in our linear system.
       likelihoodSums   =self%likelihoodSums
       stateSums        =self%stateSums
       ! Solve for regression coefficients.
       coefficients     =stateSums%linearSystemSolve(likelihoodSums)
       self%coefficients=coefficients
       ! Compute fitted likelihoods.
       if (mpiSelf%isMaster()) call Galacticus_Display_Message('fitting global trends (step #4) ['//char(Formatted_Date_and_Time())//']')
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
       if (mpiSelf%isMaster()) call Galacticus_Display_Message('computing variogram ['//char(Formatted_Date_and_Time())//']')
       k=0
       do i=1,self%emulatorRebuildCount-1
          do j=i+1,self%emulatorRebuildCount
             k=k+1
             separations  (k)=gaussianRegressionSeparation(self,self%statesCombined(:,i),self%statesCombined(:,j))
             semiVariances(k)= 0.5d0                         &
                  &           *(                             &
                  &             +self%likelihoodResiduals(i) &
                  &             -self%likelihoodResiduals(j) &
                  &            )**2
          end do
       end do
       ! Fit the variogram model.
       if (mpiSelf%isMaster()) call Galacticus_Display_Message('fitting variogram ['//char(Formatted_Date_and_Time())//']')
       call gaussianRegressionFitVariogram(self,separations,semiVariances)
       ! Compute regression matrix.
       k=0
       if (self%dumpEmulator) then
          self%dumpEmulatorCount=self%dumpEmulatorCount+1
          fileName=self%dumpEmulatorFileRoot//self%dumpEmulatorCount//".log"
          open(newUnit=fileUnit,file=char(fileName),form='formatted',status='unknown')
       end if
       do i=1,self%emulatorRebuildCount
          do j=1,self%emulatorRebuildCount
             separation                =gaussianRegressionSeparation(self,self%statesCombined(:,i),self%statesCombined(:,j))
             self%regressionMatrix(i,j)=gaussianRegressionCorrelation(self,separation)
             if (self%dumpEmulator) then
                k=k+1
                write (fileUnit,*) i,j,separation,semiVariances(k),self%regressionMatrix(i,j)
             end if
          end do
       end do
       if (self%dumpEmulator) close(fileUnit)
       ! Find the LU decomposition of the regression matrix for later use.
       if (mpiSelf%isMaster()) call Galacticus_Display_Message('computing regression matrix ['//char(Formatted_Date_and_Time())//']')
       self%regressionMatrix(1:self%emulatorRebuildCount  ,  self%emulatorRebuildCount+1)=1.0d0
       self%regressionMatrix(  self%emulatorRebuildCount+1,1:self%emulatorRebuildCount  )=1.0d0
       self%regressionMatrix(  self%emulatorRebuildCount+1,  self%emulatorRebuildCount+1)=0.0d0
       if (allocated(self%regressionMatrixLU)) deallocate(self%regressionMatrixLU)
       regressionMatrix=matrix(self%regressionMatrix)
       self%regressionMatrixLU        =matrixLU(regressionMatrix)
       determinantSign                =regressionMatrix%signDeterminant()
       self%regressionMatrixIsSingular=(determinantSign == 0)
       if (mpiSelf%isMaster().and.self%regressionMatrixIsSingular) call Galacticus_Display_Message('   ==> regression matrix is singular')
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
       if (mpiSelf%isMaster()) call Galacticus_Display_Unindent('done ['//char(Formatted_Date_and_Time())//']')
    end if
    ! Count evaluations.
    self%evaluationCount=self%evaluationCount+1
    ! Ensure arrays are allocated.
    if (.not.allocated(self%simulatorLikelihood)) then
       call allocateArray(self%simulatorLikelihood,[                            self%emulatorRebuildCount])
       call allocateArray(self%simulationState    ,[simulationState%dimension(),self%emulatorRebuildCount])
    end if
    likelihoodIsSimulated=.false.
    if (self%initialized.and..not.self%regressionMatrixIsSingular) then
       ! Perform the emulation.
       call gaussianRegressionEmulate(self,simulationState,gaussianRegressionEvaluate,likelihoodError)
       if (present(logLikelihoodVariance)) logLikelihoodVariance=likelihoodError**2
       ! Test likelihood emulation. We do this whenever the emulator is currently rated "not good", and periodically otherwise to
       ! monitor emulator behavior. When first initialized, the emulator is rated "not good" such that it has to prove that it is
       ! valid before we actually begin using it.
       if (.not.self%isGood .or. (mod(self%evaluationCount,self%reportCount) == 0 .and. self%evaluationCount > 0)) then
          ! Every so many steps we evaluate the simulated likelihood and check that our emulator is reliable.
          likelihoodEmulated                  =gaussianRegressionEvaluate
          gaussianRegressionEvaluate=self%posteriorSampleLikelihood_%evaluate(simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance)
          likelihoodIsSimulated               =.true.
          ! Check that a non-impossible likelihood was returned.
          if (gaussianRegressionEvaluate > logImpossible) then
             ! Count number of emulator checks and the failure rate.
             self%emulatorCheckCount=self%emulatorCheckCount+1
             if     (                                                                                                                    &
                  &   abs(likelihoodEmulated-gaussianRegressionEvaluate) > emulatorFailureSignificance        *likelihoodError &
                  &  .and.                                                                                                               &
                  &   abs(likelihoodEmulated-gaussianRegressionEvaluate) > emulatorFailureSignificanceAbsolute                 &
                  & ) self%emulatorFailCount=self%emulatorFailCount+1
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
          if (mpiSelf%isMaster() .and. Galacticus_Verbosity_Level() >= verbosityInfo .and. mod(self%evaluationCount,self%reportCount) == 0 .and. self%evaluationCount > 0) then
             write (label,'(i8)') failuresTotal
             message='Emulator failed '//trim(adjustl(label))//' times out of '
             write (label,'(i8)') checksTotal
             message=message//trim(adjustl(label))//' checks ('
             write (label,'(f6.2)') 100.0d0*dble(failuresTotal)/dble(checksTotal)
             message=message//trim(adjustl(label))//'% - expect '
             write (label,'(f6.2)') 100.0d0*failureRateExpected
             message=message//trim(adjustl(label))//'% for perfect emulator)'
             call Galacticus_Display_Message(message)
             if (.not.self%isGood) call Galacticus_Display_Message('WARNING: emulator failure rate is too high - emulator will not be used')
          end if
       else
          if (simulationConvergence%isConverged().and.simulationConvergence%stateIsOutlier(simulationState%chainIndex()).and.self%emulateOutliers) return
          if (gaussianRegressionEvaluate+logPriorProposed+self%sigmaBuffer*likelihoodError < logLikelihoodCurrent+logPriorCurrent-self%logLikelihoodBuffer        *temperature) return
          if (                                                             likelihoodError <                                      self%logLikelihoodErrorTolerance*temperature) return
       end if
    end if
    ! Evaluate the likelihood using the simulator, unless we have already done so.
    self%simulationCount=self%simulationCount+1
    if (.not.likelihoodIsSimulated) then
       ! If prior probability is impossible, then don't even try to simulate.
       if (logPriorProposed <= logImpossible) then
          gaussianRegressionEvaluate=logImpossible
          if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
       else
          gaussianRegressionEvaluate=self%posteriorSampleLikelihood_%evaluate(simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance)
          ! Store the likelihood and state.
          if (gaussianRegressionEvaluate > logImpossible .and. self%accumulatedStateCount < self%emulatorRebuildCount) then
             self%accumulatedStateCount                            =self%accumulatedStateCount+1
             self%simulatorLikelihood(  self%accumulatedStateCount)=gaussianRegressionEvaluate
             self%simulationState    (:,self%accumulatedStateCount)=simulationState%get()
          end if
       end if
    end if
    return
  end function gaussianRegressionEvaluate

  elemental double precision function gaussianRegressionEvaluateVariogram(self,h)
    !% Compute the variogram at separation {\normalfont \ttfamily h}.
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(in   ) :: self
    double precision                                             , intent(in   ) :: h

    if (h <= 0.0d0) then
       gaussianRegressionEvaluateVariogram=self%C0
    else if (h < self%CR) then
       gaussianRegressionEvaluateVariogram=self%C0+self%C1*(1.5d0*h/self%CR-0.5d0*(h/self%CR)**3)
    else
       gaussianRegressionEvaluateVariogram=self%C0+self%C1
    end if
    return
  end function gaussianRegressionEvaluateVariogram

  elemental double precision function gaussianRegressionCorrelation(self,h)
    !% Compute correlation of the variogram model.
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(in   ) :: self
    double precision                                             , intent(in   ) :: h

    gaussianRegressionCorrelation=self%C0+self%C1-gaussianRegressionEvaluateVariogram(self,h)
    return
  end function gaussianRegressionCorrelation

  double precision function gaussianRegressionSeparation(self,state1,state2)
    !% Determine the separation between two state vectors.
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression)              , intent(in   ) :: self
    double precision                                             , dimension(:), intent(in   ) :: state1, state2

    gaussianRegressionSeparation=sqrt(sum(((state1-state2)/self%stateScales)**2))
    return
  end function gaussianRegressionSeparation

  subroutine gaussianRegressionFitVariogram(self,separations,semiVariances)
    !% Compute best fit coefficients for the variogram model.
    use            :: Multidimensional_Minimizer, only : multiDMinimizer
    use            :: Galacticus_Error          , only : Galacticus_Error_Report
    use            :: Interface_GSL             , only : GSL_Success            , GSL_ENoProg, GSL_Continue
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    use            :: Sorting                   , only : sortIndex
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(inout), target       :: self
    double precision                                             , intent(in   ), dimension(:) :: separations             , semiVariances
    double precision                                                            , dimension(3) :: C
    integer         (c_size_t                                   ), allocatable  , dimension(:) :: rank
    integer                                                      , parameter                   :: iterationsMaximum=10000
    double precision                                             , parameter                   :: gradientTolerance=1.0d-2
    double precision                                             , parameter                   :: binWidthMaximum  =2.0d0
    integer                                                      , parameter                   :: binCountMaximum  =100000
    type            (multiDMinimizer                            )                              :: minimizer_
    integer                                                                                    :: status                  , iteration
    integer         (c_size_t                                   )                              :: k                       , j            , &
         &                                                                                        i
    double precision                                                                           :: currentMinimum

    ! Allocate workspace.
    gaussianRegressionSelf => self
    allocate(self%separationsNormalized  (size(separations)))
    allocate(self%separationsLimited     (size(separations)))
    allocate(self%semiVariancesNormalized(size(separations)))
    allocate(self%semiVariancesBinned    (size(separations)))
    allocate(self%separationsBinned      (size(separations)))
    allocate(rank                        (size(separations)))
    ! Compute normalized separations and variances.
    self%separationNormalization  =sum(separations  )/dble(size(separations))
    self%semiVarianceNormalization=sum(semiVariances)/dble(size(separations))
    self%separationsNormalized    =separations       /self%separationNormalization
    self%semiVariancesNormalized  =semiVariances     /self%semiVarianceNormalization
    ! Get rank ordering by separation.
    rank=sortIndex(self%separationsNormalized)
    ! Compute binned estimates of the mean semi-variances.
    self%binCount=0
    j       =0
    k       =1
    do while (.true.)
       j=j+1
       if (j-k == binCountMaximum .or. self%separationsNormalized(rank(j)) > binWidthMaximum*self%separationsNormalized(rank(k)) .or. j == size(separations)) then
          self%binCount=self%binCount+1
          self%separationsBinned  (self%binCount)=0.0d0
          self%semiVariancesBinned(self%binCount)=0.0d0
          do i=k,j
             self%separationsBinned  (self%binCount)=self%separationsBinned  (self%binCount)+self%separationsNormalized  (rank(i))
             self%semiVariancesBinned(self%binCount)=self%semiVariancesBinned(self%binCount)+self%semiVariancesNormalized(rank(i))
          end do
          self%separationsBinned  (self%binCount)=self%separationsBinned  (self%binCount)/dble(j-k+1)
          self%semiVariancesBinned(self%binCount)=self%semiVariancesBinned(self%binCount)/dble(j-k+1)
          k=j+1
          if (j == size(separations)) exit
       end if
    end do
    ! Build the minimizer.
    minimizer_=multiDMinimizer(3_c_size_t,gaussianRegressionVariogramModelF,gaussianRegressionVariogramModelD,gaussianRegressionVariogramModelFD)
    C         =[self%semiVariancesBinned(1),self%semiVariancesBinned(self%binCount),self%separationsBinned(self%binCount/2)] ! Initial guess for the parameters.
    call minimizer_%set(x=C,stepSize=0.01d0,tolerance=0.1d0)
    ! Iterate the minimizer until a sufficiently good solution is found.
    currentMinimum=0.0d0
    iteration     =0
    do while (                                                                             &
         &     minimizer_%testGradient(toleranceAbsolute=gradientTolerance*currentMinimum) &
         &    .and.                                                                        &
         &     iteration <  iterationsMaximum                                              &
         &   )
       iteration     =iteration+1
       call minimizer_%iterate(status)
       currentMinimum=minimizer_%minimum()
       if (status == GSL_ENoProg) exit
       if (status /= GSL_Success) call Galacticus_Error_Report('failed to iterate minimizer'//{introspection:location})
    end do
    ! Extract the best fit parameters.
    C=minimizer_%x()
    self%C0=C(1)*self%semiVarianceNormalization
    self%C1=C(2)*self%semiVarianceNormalization
    self%CR=C(3)*self%separationNormalization
    ! Clean up.
    deallocate(self%separationsNormalized  )
    deallocate(self%separationsLimited     )
    deallocate(self%semiVariancesNormalized)
    deallocate(self%separationsBinned      )
    deallocate(self%semiVariancesBinned    )
    deallocate(rank                        )
    return
  end subroutine gaussianRegressionFitVariogram

  double precision function gaussianRegressionVariogramModelF(x)
    !% Function to be minimized when fitting the variogram.
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    where (gaussianRegressionSelf%separationsBinned(1:gaussianRegressionSelf%binCount) > x(3))
       gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)=x(3)
    elsewhere
       gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)=gaussianRegressionSelf%separationsBinned(1:gaussianRegressionSelf%binCount)
    end where
    gaussianRegressionVariogramModelF=sum(((x(1)+x(2)*(1.5d0*gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3)-0.5d0*(gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3))**3))/gaussianRegressionSelf%semiVariancesBinned(1:gaussianRegressionSelf%binCount)-1.0d0)**2)
    return
  end function gaussianRegressionVariogramModelF

  function gaussianRegressionVariogramModelD(x) result(df)
    !% Derivatives of the function to be minimized when fitting the variogram.
    implicit none
    double precision, intent(in   ), dimension(     : ) :: x
    double precision               , dimension(size(x)) :: df

    where (gaussianRegressionSelf%separationsBinned > x(3))
       gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)=x(3)
    elsewhere
       gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)=gaussianRegressionSelf%separationsBinned(1:gaussianRegressionSelf%binCount)
    end where
    df(1)=2.0d0*sum(((x(1)+x(2)*(1.5d0*gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3)-0.5d0*(gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3))**3))/gaussianRegressionSelf%semiVariancesBinned(1:gaussianRegressionSelf%binCount)-1.0d0)/gaussianRegressionSelf%semiVariancesBinned(1:gaussianRegressionSelf%binCount))
    df(2)=2.0d0*sum(((x(1)+x(2)*(1.5d0*gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3)-0.5d0*(gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3))**3))/gaussianRegressionSelf%semiVariancesBinned(1:gaussianRegressionSelf%binCount)-1.0d0)*(1.5d0*gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3)-0.5d0*(gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3))**3)/gaussianRegressionSelf%semiVariancesBinned(1:gaussianRegressionSelf%binCount))
    df(3)=2.0d0*sum(((x(1)+x(2)*(1.5d0*gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3)-0.5d0*(gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3))**3))/gaussianRegressionSelf%semiVariancesBinned(1:gaussianRegressionSelf%binCount)-1.0d0)*x(2)*(-1.5d0*gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3)+1.5d0*(gaussianRegressionSelf%separationsLimited(1:gaussianRegressionSelf%binCount)/x(3))**3)/x(3)/gaussianRegressionSelf%semiVariancesBinned(1:gaussianRegressionSelf%binCount))
    where(abs(df) < 1.0d-30)
       df=1.0d-30
    end where
    return
  end function gaussianRegressionVariogramModelD

  subroutine gaussianRegressionVariogramModelFD(x,f,df)
    !% Computes both function and derivatives to be minimized when fitting the variogram.
    implicit none
    double precision, intent(in   ), dimension(     : ) :: x
    double precision, intent(  out)                     :: f
    double precision, intent(  out), dimension(size(x)) :: df

    f =gaussianRegressionVariogramModelF(x)
    df=gaussianRegressionVariogramModelD(x)
    return
  end subroutine gaussianRegressionVariogramModelFD

  integer function polynomialCoefficientCount(n,d)
    !% Return the number of coefficients at {\normalfont \ttfamily n}$^\mathrm{th}$ order in polynomial of dimension {\normalfont \ttfamily d}.
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

  function polynomialIteratorConstructor(order,rank)
    !% Create a polynomial iterator for a poynomial of specified {\normalfont \ttfamily order} and {\normalfont \ttfamily rank}.
    use :: Memory_Management, only : allocateArray
    implicit none
    type   (polynomialIterator)                :: polynomialIteratorConstructor
    integer                    , intent(in   ) :: order                        , rank

    call allocateArray(polynomialIteratorConstructor%indices,[order])
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

  subroutine gaussianRegressionFunctionChanged(self)
    !% Respond to possible changes in the likelihood function.
    use :: Memory_Management, only : deallocateArray
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
    if (allocated(self%polynomialCoefficient)) call deallocateArray(self%polynomialCoefficient)
    if (allocated(self%likelihoodSums       )) call deallocateArray(self%likelihoodSums       )
    if (allocated(self%stateSums            )) call deallocateArray(self%stateSums            )
    if (allocated(self%coefficients         )) call deallocateArray(self%coefficients         )
    if (allocated(self%regressionMatrix     )) call deallocateArray(self%regressionMatrix     )
    if (allocated(self%stateOffset          )) call deallocateArray(self%stateOffset          )
    if (allocated(self%weight               )) call deallocateArray(self%weight               )
    if (allocated(self%likelihoodResiduals  )) call deallocateArray(self%likelihoodResiduals  )
    if (allocated(self%statesCombined       )) call deallocateArray(self%statesCombined       )
    if (allocated(self%stateScales          )) call deallocateArray(self%stateScales          )
    if (allocated(self%stateMeans           )) call deallocateArray(self%stateMeans           )
    ! Let the simulator know that the likelihood function may have changed.
    call self%posteriorSampleLikelihood_%functionChanged()
    return
  end subroutine gaussianRegressionFunctionChanged

  logical function gaussianRegressionWillEvaluate(self,simulationState,modelParameters_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
    !% Return true if the log-likelihood will be evaluated.
    use :: Galacticus_Display            , only : Galacticus_Verbosity_Level     , verbosityInfo
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
       if (mod(self%evaluationCount,self%reportCount) == 0 .and. self%evaluationCount > 0 .and. Galacticus_Verbosity_Level() >= verbosityInfo) then
          ! Emulation will be tested, so the simulator is always run.
          gaussianRegressionWillEvaluate=.true.
          return
       else
          ! Evaluate the emulator.
          call gaussianRegressionEmulate(self,simulationState,likelihoodEmulated,likelihoodEmulatedError)
          ! If simulation is converged, this is an outlier chain, and we're told to emulate such outliers, then return.
          if (simulationConvergence%isConverged().and.simulationConvergence%stateIsOutlier(simulationState%chainIndex()).and.self%emulateOutliers) gaussianRegressionWillEvaluate=.false.
          ! If likelihood is well below current likelihood (and this should be sufficient that changes in priors won't affect the
          ! conclusion), then use the emulated likelihood. We scale the likelihood buffer by the current temperature to reflect the
          ! fact that transitions between states are easier when the temperature is high.
          if (likelihoodEmulated+logPriorProposed+self%sigmaBuffer*likelihoodEmulatedError < logLikelihoodCurrent+logPriorCurrent-self%logLikelihoodBuffer        *temperature) gaussianRegressionWillEvaluate=.false.
          ! Return if the error is below the tolerance. We increase the tolerance value in proportion to temperature since the
          ! likelihoods will be divided by this amount when evaluating transition probabilties.
          if (likelihoodEmulatedError                                                      <                                      self%logLikelihoodErrorTolerance*temperature) gaussianRegressionWillEvaluate=.false.
       end if
       return
    else
       ! Evaluate the simulator.
       gaussianRegressionWillEvaluate=.true.
    end if
    return
  end function gaussianRegressionWillEvaluate

  subroutine gaussianRegressionEmulate(self,simulationState,likelihoodEmulated,likelihoodEmulatedError)
    !% Evaluate the model emulator.
    use :: Linear_Algebra, only : vector, assignment(=)
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(inout)               :: self
    class           (posteriorSampleStateClass                  ), intent(inout)               :: simulationState
    double precision                                             , intent(  out)               :: likelihoodEmulated        , likelihoodEmulatedError
    double precision                                             , allocatable  , dimension(:) :: stateCurrent
    double precision                                             , parameter                   :: likelihoodErrorLarge=1.0d6
    integer                                                                                    :: i                         , j
    double precision                                                                           :: separation                , likelihoodFit
    type            (vector                                     )                              :: stateOffsetVector         , weightVector
    type            (polynomialIterator                         )                              :: iterator1

    ! Compute vector D.
    allocate(stateCurrent(simulationState%dimension()))
    stateCurrent=simulationState%get()-self%stateMeans
    do i=1,self%emulatorRebuildCount
       separation         =gaussianRegressionSeparation (self,stateCurrent,self%statesCombined(:,i))
       self%stateOffset(i)=gaussianRegressionCorrelation(self,separation                           )
    end do
    self%stateOffset(self%emulatorRebuildCount+1)=1.0d0
    ! Solve the linear system.
    stateOffsetVector=vector(self%stateOffset)
    weightVector     =self%regressionMatrixLU%squareSystemSolve(stateOffsetVector)
    self%weight      =weightVector
    ! Compute the likelihood and variance.
    likelihoodEmulated     =sum(self%likelihoodResiduals*self%weight(1:self%emulatorRebuildCount))
    likelihoodEmulatedError=gaussianRegressionCorrelation(self,0.0d0)-sum(self%weight*self%stateOffset)
    if (likelihoodEmulatedError >= 0.0d0) then
       likelihoodEmulatedError=sqrt(likelihoodEmulatedError)
    else
       likelihoodEmulatedError=likelihoodErrorLarge
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
    !% Process a previous state to restore likelihood function.
    use :: Memory_Management           , only : allocateArray
    use :: Models_Likelihoods_Constants, only : logImpossible
    implicit none
    class           (posteriorSampleLikelihoodGaussianRegression), intent(inout)               :: self
    double precision                                             , intent(in   ), dimension(:) :: simulationState
    double precision                                             , intent(in   )               :: logLikelihood
    logical                                                                                    :: storeState

    if (logLikelihood > logImpossible) then
       if (.not.allocated(self%simulatorLikelihood)) then
          call allocateArray(self%simulatorLikelihood,[                      self%emulatorRebuildCount])
          call allocateArray(self%simulationState    ,[size(simulationState),self%emulatorRebuildCount])
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

