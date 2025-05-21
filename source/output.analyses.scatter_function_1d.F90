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
  Implements a generic 1D scatter function (i.e. the scatter of some property weighted by number density of
  objects binned by some property) output analysis class.
  !!}

  use :: ISO_Varying_String, only : varying_string

  !![
  <outputAnalysis name="outputAnalysisScatterFunction1D">
   <description>A generic 1D scatter function (i.e. the scatter of some property weighted by number density of objects binned by some property) output analysis class.</description>
   <deepCopy>
    <functionClass variables="meanFunction, meanSquaredFunction"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="meanFunction, meanSquaredFunction"/>
   </stateStorable>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisScatterFunction1D
     !!{
     A generic 1D scatter function (i.e. scatter of some property weighted by number density of objects binned by some property) output analysis class.
     !!}
     private
     class           (nodePropertyExtractorClass                  ), pointer                     :: nodePropertyExtractor_                => null(), outputAnalysisWeightPropertyExtractor_ => null()
     class           (outputAnalysisPropertyOperatorClass         ), pointer                     :: outputAnalysisPropertyOperator_       => null(), outputAnalysisPropertyUnoperator_      => null(), &
          &                                                                                         outputAnalysisWeightPropertyOperator_ => null()
     class           (outputAnalysisWeightOperatorClass           ), pointer                     :: outputAnalysisWeightOperator_         => null()
     class           (outputAnalysisDistributionOperatorClass     ), pointer                     :: outputAnalysisDistributionOperator_   => null()
     class           (galacticFilterClass                         ), pointer                     :: galacticFilter_                       => null()
     class           (outputTimesClass                            ), pointer                     :: outputTimes_                          => null()
     integer         (c_size_t                                    )                              :: bufferCount
     integer                                                                                     :: covarianceBinomialBinsPerDecade
     type            (enumerationOutputAnalysisCovarianceModelType)                              :: covarianceModel
     type            (varying_string                              )                              :: label                                          , comment                                         , &
          &                                                                                         propertyLabel                                  , propertyComment                                 , &
          &                                                                                         scatterLabel                                   , scatterComment                                  , &
          &                                                                                         propertyUnits                                  , scatterUnits                                    , &
          &                                                                                         xAxisLabel                                     , yAxisLabel                                      , &
          &                                                                                         targetLabel
     double precision                                                                            :: propertyUnitsInSI                              , scatterUnitsInSI                                , &
          &                                                                                         covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum
     type            (outputAnalysisMeanFunction1D                ), pointer                     :: meanFunction                          => null(), meanSquaredFunction                    => null()
     double precision                                              , allocatable, dimension(:  ) :: binCenter                                      , scatterValue                                    , &
          &                                                                                         scatterValueTarget                             , scatterCovarianceTarget1D                       , &
          &                                                                                         outputWeight
     double precision                                              , allocatable, dimension(:,:) :: scatterCovariance                              , scatterCovarianceTarget
     logical                                                                                     :: finalized                                      , likelihoodNormalize                             , &
          &                                                                                         xAxisIsLog                                     , yAxisIsLog
   contains
     !![
     <methods>
       <method description="Return the results of the scatter function operator." method="results"         />
       <method description="Finalize analysis of the scatter function operator."  method="finalizeAnalysis"/>
     </methods>
     !!]
     final     ::                     scatterFunction1DDestructor
     procedure :: analyze          => scatterFunction1DAnalyze
     procedure :: finalize         => scatterFunction1DFinalize
     procedure :: finalizeAnalysis => scatterFunction1DFinalizeAnalysis
     procedure :: reduce           => scatterFunction1DReduce
     procedure :: results          => scatterFunction1DResults
     procedure :: logLikelihood    => scatterFunction1DLogLikelihood
  end type outputAnalysisScatterFunction1D

  interface outputAnalysisScatterFunction1D
     !!{
     Constructors for the \refClass{outputAnalysisScatterFunction1D} output analysis class.
     !!}
     module procedure scatterFunction1DConstructorParameters
     module procedure scatterFunction1DConstructorInternal
  end interface outputAnalysisScatterFunction1D

contains

  function scatterFunction1DConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisScatterFunction1D} output analysis class which takes a parameter set as input.
    !!}
    use :: Error                  , only : Error_Report
    use :: Input_Parameters       , only : inputParameter                                , inputParameters
    use :: Output_Analyses_Options, only : enumerationOutputAnalysisCovarianceModelEncode
    implicit none
    type            (outputAnalysisScatterFunction1D        )                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    class           (nodePropertyExtractorClass             ), pointer                     :: nodePropertyExtractor_               , outputAnalysisWeightPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass    ), pointer                     :: outputAnalysisPropertyOperator_      , outputAnalysisPropertyUnoperator_     , &
         &                                                                                    outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisWeightOperatorClass      ), pointer                     :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass), pointer                     :: outputAnalysisDistributionOperator_
    class           (galacticFilterClass                    ), pointer                     :: galacticFilter_
    class           (outputTimesClass                       ), pointer                     :: outputTimes_
    double precision                                         , dimension(:  ), allocatable :: binCenter                            , outputWeight                          , &
         &                                                                                    scatterValueTarget                   , scatterCovarianceTarget1D
    double precision                                         , dimension(:,:), allocatable :: scatterCovarianceTarget
    integer         (c_size_t                               )                              :: bufferCount
    type            (varying_string                         )                              :: label                                , comment                               , &
         &                                                                                    propertyLabel                        , propertyComment                       , &
         &                                                                                    scatterLabel                         , scatterComment                        , &
         &                                                                                    propertyUnits                        , scatterUnits                          , &
         &                                                                                    covarianceModel                      , xAxisLabel                            , &
         &                                                                                    yAxisLabel                           , targetLabel
    integer                                                                                :: covarianceBinomialBinsPerDecade
    type            (inputParameters                        )                              :: unoperatorParameters
    type            (inputParameters                        )                              :: weightParameters
    double precision                                                                       :: propertyUnitsInSI                    , scatterUnitsInSI                      , &
         &                                                                                    covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum
    logical                                                                                :: likelihoodNormalize                  , xAxisIsLog                            , &
         &                                                                                    yAxisIsLog

    !![
    <objectBuilder class="nodePropertyExtractor"    name="nodePropertyExtractor_"       source="parameters"          />
    <objectBuilder class="nodePropertyExtractor"    name="outputAnalysisWeightPropertyExtractor_" source="weightParameters"    />
    <objectBuilder class="outputAnalysisPropertyOperator"     name="outputAnalysisPropertyOperator_"        source="parameters"          />
    <objectBuilder class="outputAnalysisPropertyOperator"     name="outputAnalysisWeightPropertyOperator_"  source="weightParameters"    />
    <objectBuilder class="outputAnalysisPropertyOperator"     name="outputAnalysisPropertyUnoperator_"      source="unoperatorParameters"/>
    <objectBuilder class="outputAnalysisWeightOperator"       name="outputAnalysisWeightOperator_"          source="parameters"          />
    <objectBuilder class="outputAnalysisDistributionOperator" name="outputAnalysisDistributionOperator_"    source="parameters"          />
    <objectBuilder class="galacticFilter"                     name="galacticFilter_"                        source="parameters"          />
    <objectBuilder class="outputTimes"                        name="outputTimes_"                           source="parameters"          />
    !!]
    unoperatorParameters=parameters%subParameters('unoperator',requireValue=.false.)
    weightParameters    =parameters%subParameters('weight'    ,requireValue=.false.)
    allocate(binCenter   (int(parameters%count('binCenter'))                     ))
    allocate(outputWeight(int(parameters%count('binCenter'))*outputTimes_%count()))
    if (parameters%count('outputWeight') /= parameters%count('binCenter')*outputTimes_%count()) &
         & call Error_Report('incorrect number of output weights provided'//{introspection:location})
    !![
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <variable>label</variable>
      <description>A label for the analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>xAxisLabel</name>
      <source>parameters</source>
      <description>A label for the $x$-axis in a plot of this analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>yAxisLabel</name>
      <source>parameters</source>
      <description>A label for the $y$-axis in a plot of this analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>xAxisIsLog</name>
      <source>parameters</source>
      <description>If true, indicates that the $x$-axis should be logarithmic in a plot of this analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>yAxisIsLog</name>
      <source>parameters</source>
      <description>If true, indicates that the $y$-axis should be logarithmic in a plot of this analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>comment</name>
      <source>parameters</source>
      <variable>comment</variable>
      <description>A descriptive comment for the analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>propertyLabel</name>
      <source>parameters</source>
      <variable>propertyLabel</variable>
      <description>A label for the property variable.</description>
    </inputParameter>
    <inputParameter>
      <name>propertyComment</name>
      <source>parameters</source>
      <variable>propertyComment</variable>
      <description>A descriptive comment for the property variable.</description>
    </inputParameter>
    <inputParameter>
      <name>propertyUnits</name>
      <source>parameters</source>
      <variable>propertyUnits</variable>
      <description>A human-readable description of the units for the property.</description>
    </inputParameter>
    <inputParameter>
      <name>propertyUnitsInSI</name>
      <source>parameters</source>
      <variable>propertyUnitsInSI</variable>
      <description>A units for the property in the SI system.</description>
    </inputParameter>
    <inputParameter>
      <name>scatterLabel</name>
      <source>parameters</source>
      <variable>scatterLabel</variable>
      <description>A label for the scatter.</description>
    </inputParameter>
    <inputParameter>
      <name>scatterComment</name>
      <source>parameters</source>
      <variable>scatterComment</variable>
      <description>A descriptive comment for the scatter.</description>
    </inputParameter>
    <inputParameter>
      <name>scatterUnits</name>
      <source>parameters</source>
      <variable>scatterUnits</variable>
      <description>A human-readable description of the units for the scatter.</description>
    </inputParameter>
    <inputParameter>
      <name>scatterUnitsInSI</name>
      <source>parameters</source>
      <variable>scatterUnitsInSI</variable>
      <description>A units for the scatter in the SI system.</description>
    </inputParameter>
    <inputParameter>
      <name>binCenter</name>
      <source>parameters</source>
      <variable>binCenter</variable>
      <description>The value of the property at the center of each bin.</description>
    </inputParameter>
    <inputParameter>
      <name>bufferCount</name>
      <source>parameters</source>
      <variable>bufferCount</variable>
      <description>The number of buffer bins to include below and above the range of actual bins.</description>
    </inputParameter>
    <inputParameter>
      <name>outputWeight</name>
      <source>parameters</source>
      <variable>outputWeight</variable>
      <description>The weight to assign to each bin at each output.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceModel</name>
      <source>parameters</source>
      <variable>covarianceModel</variable>
      <description>The model to use for computing covariances.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing volume function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>likelihoodNormalize</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true then normalize the likelihood to make it a probability density.</description>
    </inputParameter>
    !!]
    if (parameters%isPresent('scatterValueTarget')) then
       if (parameters%isPresent('scatterCovarianceTarget')) then
          !![
          <inputParameter>
            <name>scatterValueTarget</name>
            <source>parameters</source>
            <description>The target function for likelihood calculations.</description>
          </inputParameter>
          <inputParameter>
            <name>scatterCovarianceTarget</name>
            <source>parameters</source>
            <variable>scatterCovarianceTarget1D</variable>
            <description>The target function covariance for likelihood calculations.</description>
          </inputParameter>
          !!]
          if (size(scatterCovarianceTarget1D) == size(scatterValueTarget)**2) then
             allocate(scatterCovarianceTarget(size(scatterValueTarget),size(scatterValueTarget)))
             scatterCovarianceTarget=reshape(scatterCovarianceTarget1D,shape(scatterCovarianceTarget))
          else
             call Error_Report('scatterCovariance has wrong size'//{introspection:location})
          end if
       else
          call Error_Report('scatterCovariance must be specified if scatterTarget is present'//{introspection:location})
       end if
    else
       if (parameters%isPresent('scatterCovariance')) call Error_Report('scatterTarget must be specified if scatterCovariance is present'//{introspection:location})
    end if
    !![
    <inputParameter>
      <name>targetLabel</name>
      <source>parameters</source>
      <description>A label for the target dataset in a plot of this analysis.</description>
      <defaultValue>var_str('')</defaultValue>
    </inputParameter>
    !!]
    ! Build the object.
    !![
    <conditionalCall>
     <call>
      self=outputAnalysisScatterFunction1D(                                                                                            &amp;
           &amp;                        label                                                                                        , &amp;
           &amp;                        comment                                                                                      , &amp;
           &amp;                        propertyLabel                                                                                , &amp;
           &amp;                        propertyComment                                                                              , &amp;
           &amp;                        propertyUnits                                                                                , &amp;
           &amp;                        propertyUnitsInSI                                                                            , &amp;
           &amp;                        scatterLabel                                                                                 , &amp;
           &amp;                        scatterComment                                                                               , &amp;
           &amp;                        scatterUnits                                                                                 , &amp;
           &amp;                        scatterUnitsInSI                                                                             , &amp;
           &amp;                        binCenter                                                                                    , &amp;
           &amp;                        bufferCount                                                                                  , &amp;
           &amp;                        reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),outputTimes_%count()]), &amp;
           &amp;                        nodePropertyExtractor_                                                             , &amp;
           &amp;                        outputAnalysisWeightPropertyExtractor_                                                       , &amp;
           &amp;                        outputAnalysisPropertyOperator_                                                              , &amp;
           &amp;                        outputAnalysisWeightPropertyOperator_                                                        , &amp;
           &amp;                        outputAnalysisPropertyUnoperator_                                                            , &amp;
           &amp;                        outputAnalysisWeightOperator_                                                                , &amp;
           &amp;                        outputAnalysisDistributionOperator_                                                          , &amp;
           &amp;                        galacticFilter_                                                                              , &amp;
           &amp;                        outputTimes_                                                                                 , &amp;
           &amp;                        enumerationOutputAnalysisCovarianceModelEncode(char(covarianceModel),includesPrefix=.false.) , &amp;
           &amp;                        covarianceBinomialBinsPerDecade                                                              , &amp;
           &amp;                        covarianceBinomialMassHaloMinimum                                                            , &amp;
           &amp;                        covarianceBinomialMassHaloMaximum                                                            , &amp;
           &amp;                        likelihoodNormalize                                                                          , &amp;
           &amp;                        xAxisLabel                                                                                   , &amp;
           &amp;                        yAxisLabel                                                                                   , &amp;
           &amp;                        xAxisIsLog                                                                                   , &amp;
           &amp;                        yAxisIsLog                                                                                   , &amp;
           &amp;                        targetLabel                                                                                    &amp;
           &amp;                        {conditions}                                                                                   &amp;
           &amp;                       )
     </call>
     <argument name="scatterValueTarget"      value="scatterValueTarget"      parameterPresent="parameters"/>
     <argument name="scatterCovarianceTarget" value="scatterCovarianceTarget" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"      />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"/>
    <objectDestructor name="outputAnalysisPropertyOperator_"       />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_" />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"     />
    <objectDestructor name="outputAnalysisWeightOperator_"         />
    <objectDestructor name="outputAnalysisDistributionOperator_"   />
    <objectDestructor name="galacticFilter_"                       />
    <objectDestructor name="outputTimes_"                          />
    !!]
    return
  end function scatterFunction1DConstructorParameters

  function scatterFunction1DConstructorInternal(label,comment,propertyLabel,propertyComment,propertyUnits,propertyUnitsInSI,scatterLabel,scatterComment,scatterUnits,scatterUnitsInSI,binCenter,bufferCount,outputWeight,nodePropertyExtractor_,outputAnalysisWeightPropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisWeightPropertyOperator_,outputAnalysisPropertyUnoperator_,outputAnalysisWeightOperator_,outputAnalysisDistributionOperator_,galacticFilter_,outputTimes_,covarianceModel,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,likelihoodNormalize,xAxisLabel,yAxisLabel,xAxisIsLog,yAxisIsLog,targetLabel,scatterValueTarget,scatterCovarianceTarget) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisScatterFunction1D} output analysis class for internal use.
    !!}
    use :: Output_Analysis_Property_Operators, only : outputAnalysisPropertyOperatorClass         , outputAnalysisPropertyOperatorSequence, outputAnalysisPropertyOperatorSquare, propertyOperatorList
    use :: Output_Analysis_Weight_Operators  , only : outputAnalysisWeightOperatorClass           , weightOperatorList
    use :: Output_Analyses_Options           , only : enumerationOutputAnalysisCovarianceModelType
    implicit none
    type            (outputAnalysisScatterFunction1D             )                                          :: self
    type            (varying_string                              ), intent(in   )                           :: label                                       , comment                                      , &
         &                                                                                                     propertyLabel                               , propertyComment                              , &
         &                                                                                                     scatterLabel                                , scatterComment                               , &
         &                                                                                                     propertyUnits                               , scatterUnits
    type            (varying_string                              ), intent(in   ), optional                 :: xAxisLabel                                  , yAxisLabel                                   , &
         &                                                                                                     targetLabel
    double precision                                              , intent(in   )                           :: propertyUnitsInSI                           , scatterUnitsInSI
    double precision                                              , intent(in   )          , dimension(:  ) :: binCenter
    integer         (c_size_t                                    ), intent(in   )                           :: bufferCount
    double precision                                              , intent(in   )          , dimension(:,:) :: outputWeight
    logical                                                       , intent(in   ), optional                 :: xAxisIsLog                                  , yAxisIsLog                                   , &
         &                                                                                                     likelihoodNormalize
    class           (nodePropertyExtractorClass                  ), intent(inout), target                   :: nodePropertyExtractor_            , outputAnalysisWeightPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass         ), intent(inout), target                   :: outputAnalysisPropertyOperator_             , outputAnalysisPropertyUnoperator_            , &
         &                                                                                                     outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisWeightOperatorClass           ), intent(inout), target                   :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass     ), intent(inout), target                   :: outputAnalysisDistributionOperator_
    class           (galacticFilterClass                         ), intent(inout), target                   :: galacticFilter_
    class           (outputTimesClass                            ), intent(inout), target                   :: outputTimes_
    type            (enumerationOutputAnalysisCovarianceModelType), intent(in   )                           :: covarianceModel
    integer                                                       , intent(in   ), optional                 :: covarianceBinomialBinsPerDecade
    double precision                                              , intent(in   ), optional                 :: covarianceBinomialMassHaloMinimum           , covarianceBinomialMassHaloMaximum
    double precision                                              , intent(in   ), optional, dimension(:  ) :: scatterValueTarget
    double precision                                              , intent(in   ), optional, dimension(:,:) :: scatterCovarianceTarget
    type            (weightOperatorList                          ), pointer                                 :: weightOperatorWeight_                       , weightOperatorSquared_
    type            (propertyOperatorList                        ), pointer                                 :: propertyOperators_
    type            (outputAnalysisPropertyOperatorSequence      ), pointer                                 :: outputAnalysisWeightPropertyOperatorSquaring_
    type            (outputAnalysisPropertyOperatorSquare        ), pointer                                 :: outputAnalysisWeightPropertyOperatorSquare_
    !![
    <constructorAssign variables="label, comment, propertyLabel, propertyComment, propertyUnits, propertyUnitsInSI, scatterLabel, scatterComment, scatterUnits, scatterUnitsInSI, xAxisLabel, yAxisLabel, xAxisIsLog, yAxisIsLog, targetLabel, scatterValueTarget, scatterCovarianceTarget, *nodePropertyExtractor_, *outputAnalysisWeightPropertyExtractor_, *outputAnalysisPropertyOperator_, *outputAnalysisWeightPropertyOperator_, *outputAnalysisPropertyUnoperator_, *outputAnalysisWeightOperator_, *outputAnalysisDistributionOperator_, *galacticFilter_, *outputTimes_, bufferCount, covarianceModel, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum"/>
    !!]

    ! Set properties needed for descriptor.
    if (present(scatterCovarianceTarget))                                                                  &
         & self%scatterCovarianceTarget1D=reshape(scatterCovarianceTarget,[size(scatterCovarianceTarget)])
    self       %outputWeight             =reshape(outputWeight           ,[size(outputWeight           )])
    ! Mark as unfinalized.
    self%finalized=.false.
    ! Set normalization state for likelihood.
    self%likelihoodNormalize=.true.
    if (present(likelihoodNormalize)) self%likelihoodNormalize=likelihoodNormalize
    ! Build the squaring operator.
    allocate(propertyOperators_                                )
    allocate(propertyOperators_                           %next)
    allocate(outputAnalysisWeightPropertyOperatorSquaring_     )
    allocate(outputAnalysisWeightPropertyOperatorSquare_       )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSquare_"   constructor="outputAnalysisPropertyOperatorSquare  (                  )"/>
    !!]
    propertyOperators_     %operator_ => outputAnalysisWeightPropertyOperator_
    propertyOperators_%next%operator_ => outputAnalysisWeightPropertyOperatorSquare_
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSquaring_" constructor="outputAnalysisPropertyOperatorSequence(propertyOperators_)"/>
    !!]
    ! Build normal and squared mean function 1D objects.
    allocate(self%meanFunction       )
    allocate(self%meanSquaredFunction)
    !![
    <referenceConstruct isResult="yes" owner="self" object="meanFunction">
     <constructor>
      outputAnalysisMeanFunction1D(                                               &amp;
       &amp;                       label                                        , &amp;
       &amp;                       comment                                      , &amp;
       &amp;                       propertyLabel                                , &amp;
       &amp;                       propertyComment                              , &amp;
       &amp;                       propertyUnits                                , &amp;
       &amp;                       propertyUnitsInSI                            , &amp;
       &amp;                       scatterLabel                                 , &amp;
       &amp;                       scatterComment                               , &amp;
       &amp;                       scatterUnits                                 , &amp;
       &amp;                       scatterUnitsInSI                             , &amp;
       &amp;                       binCenter                                    , &amp;
       &amp;                       bufferCount                                  , &amp;
       &amp;                       outputWeight                                 , &amp;
       &amp;                       nodePropertyExtractor_             , &amp;
       &amp;                       outputAnalysisWeightPropertyExtractor_       , &amp;
       &amp;                       outputAnalysisPropertyOperator_              , &amp;
       &amp;                       outputAnalysisWeightPropertyOperator_        , &amp;
       &amp;                       outputAnalysisPropertyUnoperator_            , &amp;
       &amp;                       outputAnalysisWeightOperator_                , &amp;
       &amp;                       outputAnalysisDistributionOperator_          , &amp;
       &amp;                       galacticFilter_                              , &amp;
       &amp;                       outputTimes_                                 , &amp;
       &amp;                       covarianceModel                              , &amp;
       &amp;                       covarianceBinomialBinsPerDecade              , &amp;
       &amp;                       covarianceBinomialMassHaloMinimum            , &amp;
       &amp;                       covarianceBinomialMassHaloMaximum              &amp;
       &amp;                      )
     </constructor>
    </referenceConstruct>
    <referenceConstruct isResult="yes" owner="self" object="meanSquaredFunction">
     <constructor>
      outputAnalysisMeanFunction1D(                                               &amp;
       &amp;                       label                                        , &amp;
       &amp;                       comment                                      , &amp;
       &amp;                       propertyLabel                                , &amp;
       &amp;                       propertyComment                              , &amp;
       &amp;                       propertyUnits                                , &amp;
       &amp;                       propertyUnitsInSI                            , &amp;
       &amp;                       scatterLabel                                 , &amp;
       &amp;                       scatterComment                               , &amp;
       &amp;                       scatterUnits                                 , &amp;
       &amp;                       scatterUnitsInSI                             , &amp;
       &amp;                       binCenter                                    , &amp;
       &amp;                       bufferCount                                  , &amp;
       &amp;                       outputWeight                                 , &amp;
       &amp;                       nodePropertyExtractor_             , &amp;
       &amp;                       outputAnalysisWeightPropertyExtractor_       , &amp;
       &amp;                       outputAnalysisPropertyOperator_              , &amp;
       &amp;                       outputAnalysisWeightPropertyOperatorSquaring_, &amp;
       &amp;                       outputAnalysisPropertyUnoperator_            , &amp;
       &amp;                       outputAnalysisWeightOperator_                , &amp;
       &amp;                       outputAnalysisDistributionOperator_          , &amp;
       &amp;                       galacticFilter_                              , &amp;
       &amp;                       outputTimes_                                 , &amp;
       &amp;                       covarianceModel                              , &amp;
       &amp;                       covarianceBinomialBinsPerDecade              , &amp;
       &amp;                       covarianceBinomialMassHaloMinimum            , &amp;
       &amp;                       covarianceBinomialMassHaloMaximum              &amp;
       &amp;                      )
     </constructor>
    </referenceConstruct>
    !!]
    ! Clean up.
    !![
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSquare_"  />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSquaring_"/>
    !!]
    nullify(weightOperatorWeight_ )
    nullify(weightOperatorSquared_)
    return
  end function scatterFunction1DConstructorInternal

  subroutine scatterFunction1DDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisScatterFunction1D} output analysis class.
    !!}
    implicit none
    type(outputAnalysisScatterFunction1D), intent(inout) :: self

    !![
    <objectDestructor name="self%meanFunction"                          />
    <objectDestructor name="self%meanSquaredFunction"                   />
    <objectDestructor name="self%nodePropertyExtractor_"                />
    <objectDestructor name="self%outputAnalysisWeightPropertyExtractor_"/>
    <objectDestructor name="self%outputAnalysisPropertyOperator_"       />
    <objectDestructor name="self%outputAnalysisWeightPropertyOperator_" />
    <objectDestructor name="self%outputAnalysisPropertyUnoperator_"     />
    <objectDestructor name="self%outputAnalysisWeightOperator_"         />
    <objectDestructor name="self%outputAnalysisDistributionOperator_"   />
    <objectDestructor name="self%galacticFilter_"                       />
    <objectDestructor name="self%outputTimes_"                          />
    !!]
    return
  end subroutine scatterFunction1DDestructor

  subroutine scatterFunction1DAnalyze(self,node,iOutput)
    !!{
    Implement a scatterFunction1D output analysis.
    !!}
    implicit none
    class  (outputAnalysisScatterFunction1D), intent(inout) :: self
    type   (treeNode                       ), intent(inout) :: node
    integer(c_size_t                       ), intent(in   ) :: iOutput

    ! Analyze for all three volume functions.
    call self%meanFunction       %analyze(node,iOutput)
    call self%meanSquaredFunction%analyze(node,iOutput)
    return
  end subroutine scatterFunction1DAnalyze

  subroutine scatterFunction1DReduce(self,reduced)
    !!{
    Implement a scatterFunction1D output analysis reduction.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisScatterFunction1D), intent(inout) :: self
    class(outputAnalysisClass            ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisScatterFunction1D)
       call self%meanFunction       %reduce(reduced%meanFunction       )
       call self%meanSquaredFunction%reduce(reduced%meanSquaredFunction)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine scatterFunction1DReduce

  subroutine scatterFunction1DFinalizeAnalysis(self)
    !!{
    Finalize analysis of a {\normalfont \ttfamily scatterFunction1D} output analysis.
    !!}
    implicit none
    class           (outputAnalysisScatterFunction1D), intent(inout)                 :: self
    integer                                                                          :: i
    double precision                                 , allocatable  , dimension(:  ) :: meanValue     , squaredValue
    double precision                                 , allocatable  , dimension(:,:) :: meanCovariance

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
    ! Retrieve results from our 1-D mean functions.
    call self%meanFunction       %results(binCenter=self%binCenter,meanValue=meanValue   ,meanCovariance=meanCovariance)
    call self%meanSquaredFunction%results(                         meanValue=squaredValue                              )
    ! Compute the scatter. For the covariance we just assume this to be half of the covariance on the mean (e.g. "Concrete
    ! Mathematics: a foundation for computer science", by Graham, Knuth, Patashnik, Addison-Wesley, Reading, 1994, page 602).
    allocate(self%scatterValue     (size(self%binCenter)                     ))
    allocate(self%scatterCovariance(size(self%binCenter),size(self%binCenter)))
    self%scatterValue     =+squaredValue      &
         &                 -meanValue     **2
    self%scatterCovariance=+0.5d0             &
         &                 *meanCovariance
    do i=1,size(self%binCenter)
       if (self%scatterValue(i) < 0.0d0) then
          ! Scatter is ill-determined - set to zero.
          self%scatterValue     (i  )=0.0d0
          self%scatterCovariance(i,:)=0.0d0
          self%scatterCovariance(:,i)=0.0d0
       else
          ! Convert from variance to scatter.
          self%scatterValue(i)=sqrt(self%scatterValue(i))
       end if
    end do
    return
  end subroutine scatterFunction1DFinalizeAnalysis

  subroutine scatterFunction1DFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily scatterFunction1D} output analysis finalization.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(outputAnalysisScatterFunction1D), intent(inout)           :: self
    type (varying_string                 ), intent(in   ), optional :: groupName
    type (hdf5Object                     )               , target   :: analysesGroup, subGroup
    type (hdf5Object                     )               , pointer  :: inGroup
    type (hdf5Object                     )                          :: analysisGroup, dataset

    ! Finalize the analysis.
    call self%finalizeAnalysis()
    ! Output the resulting scatter function.
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup(char(self%label),char(self%comment))
    ! Write metadata describing this analysis.
    call    analysisGroup%writeAttribute(     char(self%   comment   )                    ,'description'                                                                                                      )
    call    analysisGroup%writeAttribute("function1D"                                     ,'type'                                                                                                             )
    call    analysisGroup%writeAttribute(     char(self%   xAxisLabel)                    ,'xAxisLabel'                                                                                                       )
    call    analysisGroup%writeAttribute(     char(self%   yAxisLabel)                    ,'yAxisLabel'                                                                                                       )
    call    analysisGroup%writeAttribute(          self%   xAxisIsLog                     ,'xAxisIsLog'                                                                                                       )
    call    analysisGroup%writeAttribute(          self%   yAxisIsLog                     ,'yAxisIsLog'                                                                                                       )
    call    analysisGroup%writeAttribute(     char(self%propertyLabel)                    ,'xDataset'                                                                                                         )
    call    analysisGroup%writeAttribute(     char(self% scatterLabel)                    ,'yDataset'                                                                                                         )
    call    analysisGroup%writeAttribute(     char(self% scatterLabel)//"Target"          ,'yDatasetTarget'                                                                                                   )
    call    analysisGroup%writeAttribute(     char(self% scatterLabel)//"Covariance"      ,'yCovariance'                                                                                                      )
    call    analysisGroup%writeAttribute(     char(self% scatterLabel)//"CovarianceTarget",'yCovarianceTarget'                                                                                                )
    ! Write computed datasets.
    call    analysisGroup%writeDataset  (          self%binCenter                         ,char(self%propertyLabel)                       ,char(self%propertyComment)                 ,datasetReturned=dataset)
    call    dataset      %writeAttribute(     char(self%propertyUnits       )             ,'units'                                                                                                            )
    call    dataset      %writeAttribute(          self%propertyUnitsInSI                 ,'unitsInSI'                                                                                                        )
    call    dataset      %close         (                                                                                                                                                                     )
    call    analysisGroup%writeDataset  (          self%scatterValue                      ,char(self% scatterLabel)                       ,char(self% scatterComment)                 ,datasetReturned=dataset)
    call    dataset      %writeAttribute(     char(self%    scatterUnits    )             ,'units'                                                                                                            )
    call    dataset      %writeAttribute(          self%scatterUnitsInSI                  ,'unitsInSI'                                                                                                        )
    call    dataset      %close         (                                                                                                                                                                     )
    call    analysisGroup%writeDataset  (          self%scatterCovariance                 ,char(self% scatterLabel)//"Covariance"         ,char(self% scatterComment)//" [covariance]",datasetReturned=dataset)
    call    dataset      %writeAttribute("["//char(self%    scatterUnits    )//"]²"       ,'units'                                                                                                            )
    call    dataset      %writeAttribute(          self%    scatterUnitsInSI   **2        ,'unitsInSI'                                                                                                        )
    call    dataset      %close         (                                                                                                                                                                     )
    ! If available, include the log-likelihood and target dataset.
    if (allocated(self%scatterValueTarget)) then
       call analysisGroup%writeAttribute(          self%logLikelihood()                   ,'logLikelihood'                                                                                                    )
       call analysisGroup%writeAttribute(     char(self%targetLabel         )             ,'targetLabel'                                                                                                      )
       call analysisGroup%writeDataset  (          self%scatterValueTarget                ,char(self%    scatterLabel)//"Target"          ,char(self% scatterComment)                 ,datasetReturned=dataset)
       call dataset      %writeAttribute(     char(self%    scatterUnits    )             ,'units'                                                                                                            )
       call dataset      %writeAttribute(          self%scatterUnitsInSI                  ,'unitsInSI'                                                                                                        )
       call dataset      %close         (                                                                                                                                                                     )
       call analysisGroup%writeDataset  (          self%scatterCovarianceTarget           ,char(self%    scatterLabel)//"CovarianceTarget",char(self% scatterComment)//" [covariance]",datasetReturned=dataset)
       call dataset      %writeAttribute("["//char(self%    scatterUnits    )//"]²"       ,'units'                                                                                                            )
       call dataset      %writeAttribute(          self%    scatterUnitsInSI   **2        ,'unitsInSI'                                                                                                        )
       call dataset      %close         (                                                                                                                                                                     )
    end if
    call    analysisGroup%close         (                                                                                                                                                                     )
    if (present(groupName)) &
         & call subGroup %close         (                                                                                                                                                                     )
    call    analysesGroup%close         (                                                                                                                                                                     )
    !$ call hdf5Access%unset()
    return
  end subroutine scatterFunction1DFinalize

  subroutine scatterFunction1DResults(self,binCenter,scatterValue,scatterCovariance)
    !!{
    Implement a scatterFunction1D output analysis finalization.
    !!}
    implicit none
    class           (outputAnalysisScatterFunction1D)                             , intent(inout)           :: self
    double precision                                 , allocatable, dimension(:  ), intent(inout), optional :: binCenter     , scatterValue
    double precision                                 , allocatable, dimension(:,:), intent(inout), optional :: scatterCovariance

    ! Finalize analysis.
    call self%finalizeAnalysis()
    ! Return results.
    if (present(binCenter        )) then
       if (allocated(binCenter        )) deallocate(binCenter        )
       allocate(binCenter(size(self%binCenter)))
       binCenter         =self%binCenter
    end if
    if (present(scatterValue     )) then
       if (allocated(scatterValue     )) deallocate(scatterValue     )
       allocate(scatterValue(size(self%scatterValue)))
       scatterValue     =self%scatterValue
    end if
    if (present(scatterCovariance)) then
       if (allocated(scatterCovariance)) deallocate(scatterCovariance)
       allocate(scatterCovariance(size(self%scatterCovariance,dim=1),size(self%scatterCovariance,dim=2)))
       scatterCovariance=self%scatterCovariance
    end if
    return
  end subroutine scatterFunction1DResults

  double precision function scatterFunction1DLogLikelihood(self)
    !!{
    Return the log-likelihood of a scatterFunction1D output analysis.
    !!}
    use :: Error                       , only : Error_Report
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Numerical_Constants_Math    , only : Pi
    use :: Interface_GSL               , only : GSL_Success
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisScatterFunction1D), intent(inout)                 :: self
    double precision                                 , allocatable  , dimension(:,:) :: scatterCovarianceCombined
    double precision                                 , allocatable  , dimension(:  ) :: scatterValueDifference
    type            (vector                         )                                :: residual
    type            (matrix                         )                                :: covariance
    integer                                                                          :: status

    ! Check for existence of a target distribution.
    if (allocated(self%scatterValueTarget)) then
       ! Finalize analysis.
       call self%finalizeAnalysis()
       ! Allocate workspaces.
       allocate(scatterCovarianceCombined(size(self%binCenter),size(self%binCenter)))
       allocate(scatterValueDifference   (size(self%binCenter)                     ))
       ! Find combined covariance and difference between model and target.
       scatterValueDifference   =+self%scatterValue            &
            &                    -self%scatterValueTarget
       scatterCovarianceCombined=+self%scatterCovariance       &
            &                    +self%scatterCovarianceTarget
       residual                 = vector(scatterValueDifference   )
       covariance               = matrix(scatterCovarianceCombined)
       ! Compute the log-likelihood.
       scatterFunction1DLogLikelihood          =-0.5d0*covariance%covarianceProduct(residual,status)
       if (status == GSL_Success) then
          if (self%likelihoodNormalize)                                                         &
               & scatterFunction1DLogLikelihood=+scatterFunction1DLogLikelihood                 &
               &                                -0.5d0*covariance%logarithmicDeterminant()      &
               &                                -0.5d0*dble(size(self%binCenter))*log(2.0d0*Pi)
       else
          scatterFunction1DLogLikelihood       =+logImprobable
       end if
    else
       scatterFunction1DLogLikelihood=0.0d0
       call Error_Report('no target distribution was provided for likelihood calculation'//{introspection:location})
    end if
    return
  end function scatterFunction1DLogLikelihood
