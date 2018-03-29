!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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
  
  !% Implementation of a posterior sampling likelihood class which implements a likelihood for \glc\ models.
  
  !# <posteriorSampleLikelihood name="posteriorSampleLikelihoodGalacticus">
  !#  <description>A posterior sampling likelihood class which implements a likelihood for \glc\ models.</description>
  !# </posteriorSampleLikelihood>
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodGalacticus
     !% Implementation of a posterior sampling likelihood class which implements a likelihood for \glc\ models.
     private
     type            (varying_string)                            :: timeCommand        , name                 , &
          &                                                         executable         , compilation          , &
          &                                                         baseParameters     , workPath             , &
          &                                                         scratchPath        , failPath
     type            (varying_string), allocatable, dimension(:) :: environment
     logical                                                     :: report             , randomize            , &
          &                                                         saveState          , cleanUp              , &
          &                                                         coreDump           , useFixedTrees        , &
          &                                                         fixedTreesInScratch, adjustMasses         , &
          &                                                         adjustOutputs
     integer                                                     :: storeCountPrevious , threads              , &
          &                                                         cpuLimit           , treesPerDecadeMinimum
     double precision                                            :: delayInterval
   contains
     procedure :: evaluate        => galacticusEvaluate
     procedure :: functionChanged => galacticusFunctionChanged
     procedure :: willEvaluate    => galacticusWillEvaluate
  end type posteriorSampleLikelihoodGalacticus

  interface posteriorSampleLikelihoodGalacticus
     !% Constructors for the {\normalfont \ttfamily galacticus} posterior sampling likelihood class.
     module procedure galacticusConstructorParameters
     module procedure galacticusConstructorInternal
  end interface posteriorSampleLikelihoodGalacticus

contains

  function galacticusConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily galacticus} posterior sampling likelihood class which builds the object
    !% from a parameter set.
    use Input_Parameters
    implicit none
    type            (posteriorSampleLikelihoodGalacticus)                            :: self
    type            (inputParameters                    ), intent(inout)             :: parameters
    type            (varying_string)                                                 :: failPath             , name         , &
         &                                                                              executable           , compilation  , &
         &                                                                              baseParameters       , workPath     , &
         &                                                                              scratchPath
    type            (varying_string                     ), allocatable, dimension(:) :: environment
    logical                                                                          :: report               , randomize    , &
         &                                                                              saveState            , cleanUp      , &
         &                                                                              coreDump             , useFixedTrees, &
         &                                                                              fixedTreesInScratch  , adjustMasses , &
         &                                                                              adjustOutputs
    integer                                                                          :: treesPerDecadeMinimum, threads      , &
         &                                                                              cpuLimit
    double precision                                                                 :: delayInterval

    !# <inputParameter>
    !#   <name>delayInterval</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The time (in seconds) to sleep between launching successive \glc\ models.</description>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>name</name>
    !#   <cardinality>1</cardinality>
    !#   <description>A name for the \glc\ model.</description>
    !#   <defaultValue>var_str('galacticus')</defaultValue>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>executable</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of the executable to run.</description>
    !#   <defaultValue>var_str('Galacticus.exe')</defaultValue>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>compilation</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The constraint compilation file to use.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>baseParameters</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The base set of parameters to use.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>workPath</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The path to the work directory.</description>
    !#   <defaultValue>var_str('./')</defaultValue>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>failPath</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The path to the failure directory.</description>
    !#   <defaultValue>var_str('./')</defaultValue>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>scratchPath</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The path to the scratch directory.</description>
    !#   <defaultValue>var_str('./')</defaultValue>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>report</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, report on progress.</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomize</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, randomize models (i.e. change the random seed).</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>saveState</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, save the state of models as they are running.</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>cleanUp</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, clean up after likelihood evaluation is complete by removing all files generated.</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>coreDump</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, create a core dump when a model fails.</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>useFixedTrees</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, used a fixed set of merger trees for the models.</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>fixedTreesInScratch</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, any fixed set of merger trees used will be stored in the scratch directory.</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>adjustMasses</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, allow constraints to adjust the range of halo masses run in the model.</description>
    !#   <defaultValue>.true.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>adjustOutputs</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, allow constraints to adjust the set of output redshifts used in the model.</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>treesPerDecadeMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The minimum number of trees per decade to run in any model.</description>
    !#   <defaultValue>10</defaultValue>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>threads</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The number of OpenMP threads to use for each model.</description>
    !#   <defaultValue>1</defaultValue>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>cpuLimit</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The maximum average CPU time per thread for which a model should be allowed to run.</description>
    !#   <defaultValue>31557600</defaultValue>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    allocate(environment(parameters%copiesCount('environment',zeroIfNotPresent=.true.)))
    do i=1,size(environment)
       !# <inputParameter>
       !#   <name>environment</name>
       !#   <variable>environment(i)</variable>
       !#   <instance>i</instance>
       !#   <cardinality>1</cardinality>
       !#   <description>.</description>
       !#   <source>parameters</source>
       !#   <type>string</type>
       !# </inputParameter>
    end do
    self=posteriorSampleLikelihoodGalacticus(delayInterval,name,executable,compilation,baseParameters,workPath,scratchPath,failPath,environment,report,randomize,saveState,cleanUp,coreDump,useFixedTrees,fixedTreesInScratch,adjustMasses,adjustOutputs,treesPerDecadeMinimum,threads,cpuLimit)
    !# <inputParametersValidate source="parameters"/>
    return
  end function galacticusConstructorParameters

  function galacticusConstructorInternal(delayInterval,name,executable,compilation,baseParameters,workPath,scratchPath,failPath,environment,report,randomize,saveState,cleanUp,coreDump,useFixedTrees,fixedTreesInScratch,adjustMasses,adjustOutputs,treesPerDecadeMinimum,threads,cpuLimit) result(self)
    !% Constructor for ``galacticus'' posterior sampling likelihood class.
    use File_Utilities
    use Galacticus_Error
    implicit none
    type            (posteriorSampleLikelihoodGalacticus)                                           :: self
    double precision                                     , intent(in   )                            :: delayInterval
    type            (varying_string                     ), intent(in   )                            :: failPath             , name         , &
         &                                                                                             executable           , compilation  , &
         &                                                                                             baseParameters       , workPath     , &
         &                                                                                             scratchPath
    type            (varying_string                     ), intent(in   ), allocatable, dimension(:) :: environment
    logical                                              , intent(in   )                            :: report               , randomize    , &
         &                                                                                             saveState            , cleanUp      , &
         &                                                                                             coreDump             , useFixedTrees, &
         &                                                                                             fixedTreesInScratch  , adjustMasses , &
         &                                                                                             adjustOutputs
    integer                                              , intent(in   )                            :: treesPerDecadeMinimum, threads      , &
         &                                                                                             cpuLimit
    !# <constructorAssign variables="delayInterval,name,executable,compilation,baseParameters,workPath,scratchPath,failPath,environment,report,randomize,saveState,cleanUp,coreDump,useFixedTrees,fixedTreesInScratch,adjustMasses,adjustOutputs,treesPerDecadeMinimum,threads,cpuLimit"/>

    self%storeCountPrevious=0
    ! Find the time command.
    self%timeCommand       =Executable_Find('time')
    if (self%timeCommand == "") call Galacticus_Error_Report('a working GNU time command is required'//{introspection:location})
    return
  end function galacticusConstructorInternal

  double precision function galacticusEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance)
    !% Return the log-likelihood for the \glc\ likelihood function. This function runs an external script that drives \glc\ and
    !% writes the resulting likelihood to file. This function then reads that likelihood from file. We make use of {\normalfont
    !% \ttfamily tmpfs} for this likelihood file so that no disk I/O is required.
    use Models_Likelihoods_Constants
    use Posterior_Sampling_State
    use Posterior_Sampling_Convergence
    use System_Command
    use String_Handling
    use MPI_Utilities
    use Galacticus_Error
    use File_Utilities
    implicit none
    class           (posteriorSampleLikelihoodGalacticus), intent(inout)               :: self
    class           (posteriorSampleStateClass          ), intent(inout)               :: simulationState
    type            (modelParameterList                 ), intent(in   ), dimension(:) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass    ), intent(inout)               :: simulationConvergence
    double precision                                     , intent(in   )               :: temperature           , logLikelihoodCurrent    , &
         &                                                                                logPriorCurrent       , logPriorProposed
    real                                                 , intent(inout)               :: timeEvaluate
    double precision                                     , intent(  out), optional     :: logLikelihoodVariance
    double precision                                     , allocatable  , dimension(:) :: stateVector
    type            (varying_string                     )                              :: wrapperCommand        , likelihoodFileName      , &
         &                                                                                timingFileName        , message
    character       (len=23                             )                              :: temperatureLabel      , reportLabel             , &
         &                                                                                randomizeLabel        , saveStateLabel          , &
         &                                                                                cleanUpLabel          , coreDumpLabel           , &
         &                                                                                useFixedTreesLabel    , fixedTreesInScratchLabel, &
         &                                                                                adjustMassesLabel     , adjustOutputsLabel      , &
         &                                                                                parameterValueLabel
    integer                                                                            :: i                     , likelihoodFileUnit      , &
         &                                                                                timingFileUnit        , status
    logical                                                                            :: storeModel
    real                                                                               :: timeSystem            , timeUser
    !GCC$ attributes unused :: logPriorCurrent, logLikelihoodCurrent

    ! If prior probability is impossible, then no need to waste time evaluating the likelihood.
    if (logPriorProposed <= logImpossible) then
       galacticusEvaluate=logImpossible
       if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
       return
    end if
    ! Generate the file name for the likelihood and timing files.
    likelihoodFileName='/dev/shm/galacticusLikelihood_'//mpiSelf%rankLabel()//'.dat'
    timingFileName    ='/dev/shm/glcTiming_'           //mpiSelf%rankLabel()//'.txt'
    ! Generate labels.
    write (temperatureLabel        ,'(e12.6)')      temperature 
    write (reportLabel             ,'(L1)'   ) self%report
    write (randomizeLabel          ,'(L1)'   ) self%randomize
    write (saveStateLabel          ,'(L1)'   ) self%saveState
    write (cleanUpLabel            ,'(L1)'   ) self%cleanUp
    write (coreDumpLabel           ,'(L1)'   ) self%coreDump
    write (useFixedTreesLabel      ,'(L1)'   ) self%useFixedTrees
    write (fixedTreesInScratchLabel,'(L1)'   ) self%fixedTreesInScratch
    write (adjustMassesLabel       ,'(L1)'   ) self%adjustMasses
    write (adjustOutputsLabel      ,'(L1)'   ) self%adjustOutputs
    ! Generate the command to run the Galacticus model.
    wrapperCommand=                                                                                 &
         & self%timeCommand                                                                 //' '// &
         & '--format "%S %U"'                                                               //' '// &
         & '--output '                            //             timingFileName             //' '// &
         & './constraints/galacticusLikelihood.pl'                                          //' '// &
         & '--mpiRank '                           //     mpiSelf%rankLabel               () //' '// &
         & '--likelihoodFile '                    //             likelihoodFileName         //' '// &
         & '--name '                              //     self   %name                       //' '// &
         & '--executable '                        //     self   %executable                 //' '// &
         & '--compilation '                       //     self   %compilation                //' '// &
         & '--baseParameters '                    //     self   %baseParameters             //' '// &
         & '--workPath '                          //     self   %workPath                   //' '// &
         & '--scratchPath '                       //     self   %scratchPath                //' '// &
         & '--failPath '                          //     self   %failPath                   //' '// &
         & '--report '                            //trim(        reportLabel               )//' '// &
         & '--randomize '                         //trim(        randomizeLabel            )//' '// &
         & '--saveState '                         //trim(        saveStateLabel            )//' '// &
         & '--cleanUp '                           //trim(        cleanUpLabel              )//' '// & 
         & '--coreDump '                          //trim(        coreDumpLabel             )//' '// &
         & '--useFixedTrees  '                    //trim(        useFixedTreesLabel        )//' '// &
         & '--fixedTreesInScratch '               //trim(        fixedTreesInScratchLabel  )//' '// &
         & '--adjustMasses  '                     //trim(        adjustMassesLabel         )//' '// &
         & '--adjustOutputs '                     //trim(        adjustOutputsLabel        )//' '// &
         & '--temperature '                       //trim(        temperatureLabel          )//' '// &
         & '--threads '                           //     self   %threads                    //' '// &
         & '--cpuLimit '                          //     self   %cpuLimit                   //' '// &
         & '--treesPerDecadeMinimum '             //     self   %treesPerDecadeMinimum
    do i=1,size(self%environment)
       wrapperCommand=wrapperCommand//'--environment '//self%environment(i)//' '
    end do
    ! Determine whether to store this model.
    select type (simulationState)
       class is (posteriorSampleStateCorrelation)
       storeModel=                                                                     &
            &       simulationConvergence%isConverged   (                            ) &
            & .and.                                                                    &
            &  .not.simulationConvergence%stateIsOutlier(simulationState%chainIndex()) &
            & .and.                                                                    &
            &       simulationState%count               (                            ) &
            &      >=                                                                  &
            &       self           %storeCountPrevious                                 &
            &      +                                                                   &
            &       simulationState%correlationLength   (                            )
       class default
       storeModel=.false.
    end select
    if (storeModel) self%storeCountPrevious=simulationState%count()
    if (storeModel) then
       wrapperCommand=wrapperCommand//' --storeModel '//simulationState%count()
    else
       wrapperCommand=wrapperCommand//' --storeModel none'
    end if
    ! Append the current state.
    stateVector=simulationState%get()
    do i=1,size(stateVector)
       stateVector(i)=modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
       write (parameterValueLabel,'(1x,e22.15)') stateVector(i)
       wrapperCommand=wrapperCommand//" --parameter "//modelParametersActive_(i)%modelParameter_%name()//"="//trim(adjustl(parameterValueLabel))
    end do
    deallocate(stateVector)
    ! Append any derived parameters
    do i=1,size(modelParametersInactive_)
       select type (modelParameter_ => modelParametersInactive_(i)%modelParameter_)
       class is (modelParameterDerived)
          wrapperCommand=wrapperCommand//" --parameter "//modelParameter_%name()//"=="//modelParameter_%definition()
       end select
    end do
    ! If requested, sleep for a short time to desynchronize launching of the Galacticus wrapper command. On some systems this is
    ! necessary to avoid resource conflicts.
    if (self%delayInterval > 0.0d0) call Sleep(int(self%delayInterval*dble(mpiSelf%rank())))
    ! Run the Galacticus model and evaluate likelihood.
    call System_Command_Do(wrapperCommand,status)
    if (status /= 0) then     
       message="Galacticus wrapper failed with status code "
       message=message//status
       if (File_Exists(likelihoodFileName)) then
          message=message//" (likelihood file '"//char(likelihoodFileName)//"' does exist)"
       else
          message=message//" (likelihood file '"//char(likelihoodFileName)//"' DOES NOT exist)"
       end if
       call Galacticus_Error_Report(message//{introspection:location})
    end if
    ! Read the likelihood back from file.
    open(newUnit=likelihoodFileUnit,file=char(likelihoodFileName),status='old',form='formatted')
    read                                     (likelihoodFileUnit,*) galacticusEvaluate
    if (present(logLikelihoodVariance)) read (likelihoodFileUnit,*) logLikelihoodVariance
    close(likelihoodFileUnit,status='delete')
    ! Read the timing data back from file.
    open(newUnit=timingFileUnit,file=char(timingFileName),status='old',form='formatted')
    read (timingFileUnit,*) timeSystem,timeUser
    close(timingFileUnit,status='delete')
    timeEvaluate=timeSystem+timeUser
    return
  end function galacticusEvaluate

  subroutine galacticusFunctionChanged(self)
    !% Respond to possible changes in the likelihood function.
    implicit none
    class(posteriorSampleLikelihoodGalacticus), intent(inout) :: self
    !GCC$ attributes unused :: self

    return
  end subroutine galacticusFunctionChanged

  logical function galacticusWillEvaluate(self,simulationState,modelParameters_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
    !% Return true if the log-likelihood will be evaluated.
    use Posterior_Sampling_State
    use Posterior_Sampling_Convergence
    use Models_Likelihoods_Constants
    implicit none
    class           (posteriorSampleLikelihoodGalacticus), intent(inout)               :: self
    class           (posteriorSampleStateClass          ), intent(inout)               :: simulationState
    type            (modelParameterList                 ), intent(in   ), dimension(:) :: modelParameters_
    class           (posteriorSampleConvergenceClass    ), intent(inout)               :: simulationConvergence
    double precision                                     , intent(in   )               :: temperature          , logLikelihoodCurrent, &
         &                                                                                logPriorCurrent      , logPriorProposed
    !GCC$ attributes unused :: self, simulationState, modelParameters_, simulationConvergence, temperature, logLikelihoodCurrent, logPriorCurrent

    ! Likelihood will not be evaluated if the proposed prior is impossible.
    galacticusWillEvaluate=(logPriorProposed > logImpossible)
    return
  end function galacticusWillEvaluate
