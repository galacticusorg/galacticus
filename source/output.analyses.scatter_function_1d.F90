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

  !!{RST
  Implements a generic 1D scatter function (i.e. the scatter of some property weighted by number density of objects binned by some property) output analysis class.
  !!}

  use :: ISO_Varying_String         , only : varying_string
  use :: Output_Analysis_Target_Data, only : outputAnalysisTargetDataClass, outputAnalysisTargetDataStandard

  !![
  <outputAnalysis name="outputAnalysisScatterFunction1D" docformat="rst">
   <description>
   A generic 1D scatter function (i.e. the scatter of some property weighted by number density of objects binned by some property) output analysis class.
   </description>
   <deepCopy>
    <functionClass variables="meanFunction, meanSquaredFunction"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="meanFunction, meanSquaredFunction"/>
   </stateStorable>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisScatterFunction1D
     !!{RST
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
          &                                                                                         propertyQuantity                               , scatterQuantity                                 , &
          &                                                                                         xAxisLabel                                     , yAxisLabel                                      , &
          &                                                                                         targetLabel
     ! Axis labels, log-scale flags, and the (optional) target dataset are also bundled into a
     ! single `outputAnalysisTargetDataStandard` instance so the wrapper-pipeline doesn't have
     ! to enumerate 2^N presence combinations for the otherwise individually-optional fields.
     ! The shadow `xAxisLabel` / `yAxisLabel` / `targetLabel` / `xAxisIsLog` / `yAxisIsLog` /
     ! `scatterValueTarget` / `scatterCovarianceTarget1D` fields are kept on the outer type so
     ! the auto-built descriptor (which walks the type definition, not contained types) can
     ! reconstruct a parameter block that recreates this object.
     type            (outputAnalysisTargetDataStandard            )                              :: targetData_
     double precision                                                                            :: propertyUnitsInSI                              , scatterUnitsInSI                                , &
          &                                                                                         covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum
     logical                                                                                     :: propertyIsComoving                             , scatterIsComoving
     type            (outputAnalysisMeanFunction1D                ), pointer                     :: meanFunction                          => null(), meanSquaredFunction                    => null()
     double precision                                              , allocatable, dimension(:  ) :: binCenter                                      , scatterValue                                    , &
          &                                                                                         scatterValueTarget                             , scatterCovarianceTarget1D                       , &
          &                                                                                         outputWeight
     double precision                                              , allocatable, dimension(:,:) :: scatterCovariance
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
     !!{RST
     Constructors for the ``outputAnalysisScatterFunction1D`` output analysis class.
     !!}
     module procedure scatterFunction1DConstructorParameters
     module procedure scatterFunction1DConstructorInternal
  end interface outputAnalysisScatterFunction1D

contains

  function scatterFunction1DConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``outputAnalysisScatterFunction1D`` output analysis class which takes a parameter set as input.
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
    type            (outputAnalysisTargetDataStandard       )                              :: targetData_
    double precision                                         , dimension(:  ), allocatable :: binCenter                            , outputWeight                          , &
         &                                                                                    scatterValueTarget                   , scatterCovarianceTarget1D
    double precision                                         , dimension(:,:), allocatable :: scatterCovarianceTarget
    integer         (c_size_t                               )                              :: bufferCount
    type            (varying_string                         )                              :: label                                , comment                               , &
         &                                                                                    propertyLabel                        , propertyComment                       , &
         &                                                                                    scatterLabel                         , scatterComment                        , &
         &                                                                                    propertyUnits                        , scatterUnits                          , &
         &                                                                                    propertyQuantity                     , scatterQuantity                       , &
         &                                                                                    covarianceModel                      , xAxisLabel                            , &
         &                                                                                    yAxisLabel                           , targetLabel
    integer                                                                                :: covarianceBinomialBinsPerDecade
    type            (inputParameters                        )                              :: unoperatorParameters
    type            (inputParameters                        )                              :: weightParameters
    double precision                                                                       :: propertyUnitsInSI                    , scatterUnitsInSI                      , &
         &                                                                                    covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum
    logical                                                                                :: xAxisIsLog                           , yAxisIsLog                            , &
         &                                                                                    propertyIsComoving                   , scatterIsComoving                     , &
         &                                                                                    likelihoodNormalize

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
    <inputParameter docformat="rst">
      <name>label</name>
      <source>parameters</source>
      <variable>label</variable>
      <description>
      A label for the analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>xAxisLabel</name>
      <source>parameters</source>
      <description>
      A label for the :math:`x`-axis in a plot of this analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>yAxisLabel</name>
      <source>parameters</source>
      <description>
      A label for the :math:`y`-axis in a plot of this analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>xAxisIsLog</name>
      <source>parameters</source>
      <description>
      If true, indicates that the :math:`x`-axis should be logarithmic in a plot of this analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>yAxisIsLog</name>
      <source>parameters</source>
      <description>
      If true, indicates that the :math:`y`-axis should be logarithmic in a plot of this analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>comment</name>
      <source>parameters</source>
      <variable>comment</variable>
      <description>
      A descriptive comment for the analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyLabel</name>
      <source>parameters</source>
      <variable>propertyLabel</variable>
      <description>
      A label for the property variable.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyComment</name>
      <source>parameters</source>
      <variable>propertyComment</variable>
      <description>
      A descriptive comment for the property variable.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyUnits</name>
      <source>parameters</source>
      <variable>propertyUnits</variable>
      <description>
      A human-readable description of the units for the property.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyQuantity</name>
      <source>parameters</source>
      <variable>propertyQuantity</variable>
      <description>
      An ``astropy.units``-parseable units string for the property.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyIsComoving</name>
      <source>parameters</source>
      <variable>propertyIsComoving</variable>
      <description>
      If true, the property is in comoving units.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>propertyUnitsInSI</name>
      <source>parameters</source>
      <variable>propertyUnitsInSI</variable>
      <description>
      A units for the property in the SI system.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>scatterLabel</name>
      <source>parameters</source>
      <variable>scatterLabel</variable>
      <description>
      A label for the scatter.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>scatterComment</name>
      <source>parameters</source>
      <variable>scatterComment</variable>
      <description>
      A descriptive comment for the scatter.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>scatterUnits</name>
      <source>parameters</source>
      <variable>scatterUnits</variable>
      <description>
      A human-readable description of the units for the scatter.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>scatterQuantity</name>
      <source>parameters</source>
      <variable>scatterQuantity</variable>
      <description>
      An ``astropy.units``-parseable units string for the scatter.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>scatterIsComoving</name>
      <source>parameters</source>
      <variable>scatterIsComoving</variable>
      <description>
      If true, the scatter is in comoving units.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>scatterUnitsInSI</name>
      <source>parameters</source>
      <variable>scatterUnitsInSI</variable>
      <description>
      A units for the scatter in the SI system.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>binCenter</name>
      <source>parameters</source>
      <variable>binCenter</variable>
      <description>
      The value of the property at the center of each bin.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>bufferCount</name>
      <source>parameters</source>
      <variable>bufferCount</variable>
      <description>
      The number of buffer bins to include below and above the range of actual bins.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputWeight</name>
      <source>parameters</source>
      <variable>outputWeight</variable>
      <description>
      The weight to assign to each bin at each output.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceModel</name>
      <source>parameters</source>
      <variable>covarianceModel</variable>
      <description>
      The model to use for computing covariances.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>
      The number of bins per decade of halo mass to use when constructing volume function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d8</defaultValue>
      <description>
      The minimum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d16</defaultValue>
      <description>
      The maximum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>likelihoodNormalize</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true then normalize the likelihood to make it a probability density.
      </description>
    </inputParameter>
    !!]
    if (parameters%isPresent('scatterValueTarget')) then
       if (parameters%isPresent('scatterCovarianceTarget')) then
          !![
          <inputParameter docformat="rst">
            <name>scatterValueTarget</name>
            <source>parameters</source>
            <description>
            The target function for likelihood calculations.
            </description>
          </inputParameter>
          <inputParameter docformat="rst">
            <name>scatterCovarianceTarget</name>
            <source>parameters</source>
            <variable>scatterCovarianceTarget1D</variable>
            <description>
            The target function covariance for likelihood calculations.
            </description>
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
    <inputParameter docformat="rst">
      <name>targetLabel</name>
      <source>parameters</source>
      <description>
      A label for the target dataset in a plot of this analysis.
      </description>
      <defaultValue>var_str('')</defaultValue>
    </inputParameter>
    !!]
    ! Bundle the (potentially partial) target data into a single object for the internal constructor.
    targetData_=outputAnalysisTargetDataStandard(                                         &
         &                                       xAxisLabel      =xAxisLabel            , &
         &                                       yAxisLabel      =yAxisLabel            , &
         &                                       targetLabel     =targetLabel           , &
         &                                       xAxisIsLog      =xAxisIsLog            , &
         &                                       yAxisIsLog      =yAxisIsLog            , &
         &                                       valueTarget     =scatterValueTarget    , &
         &                                       covarianceTarget=scatterCovarianceTarget &
         &                                      )
    ! Build the object.  No conditional arguments to dispatch on (scatter doesn't take a
    ! `binWidth`), so a plain Fortran call is enough — no need for a `<conditionalCall>`
    ! directive.
    self=outputAnalysisScatterFunction1D(                                                                                       &
         &                               label                                                                                , &
         &                               comment                                                                              , &
         &                               propertyLabel                                                                        , &
         &                               propertyComment                                                                      , &
         &                               propertyUnits                                                                        , &
         &                               propertyQuantity                                                                     , &
         &                               propertyIsComoving                                                                   , &
         &                               propertyUnitsInSI                                                                    , &
         &                               scatterLabel                                                                         , &
         &                               scatterComment                                                                       , &
         &                               scatterUnits                                                                         , &
         &                               scatterQuantity                                                                      , &
         &                               scatterIsComoving                                                                    , &
         &                               scatterUnitsInSI                                                                     , &
         &                               binCenter                                                                            , &
         &                               bufferCount                                                                          , &
         &                               reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),outputTimes_%count()]), &
         &                               nodePropertyExtractor_                                                               , &
         &                               outputAnalysisWeightPropertyExtractor_                                               , &
         &                               outputAnalysisPropertyOperator_                                                      , &
         &                               outputAnalysisWeightPropertyOperator_                                                , &
         &                               outputAnalysisPropertyUnoperator_                                                    , &
         &                               outputAnalysisWeightOperator_                                                        , &
         &                               outputAnalysisDistributionOperator_                                                  , &
         &                               galacticFilter_                                                                      , &
         &                               outputTimes_                                                                         , &
         &                               enumerationOutputAnalysisCovarianceModelEncode(char(covarianceModel),includesPrefix=.false.) , &
         &                               covarianceBinomialBinsPerDecade                                                      , &
         &                               covarianceBinomialMassHaloMinimum                                                    , &
         &                               covarianceBinomialMassHaloMaximum                                                    , &
         &                               likelihoodNormalize                                                                  , &
         &                               targetData_                                                                            &
         &                              )
    !![
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

  function scatterFunction1DConstructorInternal(label,comment,propertyLabel,propertyComment,propertyUnits,propertyQuantity,propertyIsComoving,propertyUnitsInSI,scatterLabel,scatterComment,scatterUnits,scatterQuantity,scatterIsComoving,scatterUnitsInSI,binCenter,bufferCount,outputWeight,nodePropertyExtractor_,outputAnalysisWeightPropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisWeightPropertyOperator_,outputAnalysisPropertyUnoperator_,outputAnalysisWeightOperator_,outputAnalysisDistributionOperator_,galacticFilter_,outputTimes_,covarianceModel,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,likelihoodNormalize,targetData_) result (self)
    !!{RST
    Constructor for the ``outputAnalysisScatterFunction1D`` output analysis class for internal use.
    !!}
    use :: Error                             , only : Error_Report
    use :: Output_Analysis_Property_Operators, only : outputAnalysisPropertyOperatorClass         , outputAnalysisPropertyOperatorSequence, outputAnalysisPropertyOperatorSquare, propertyOperatorList
    use :: Output_Analysis_Weight_Operators  , only : outputAnalysisWeightOperatorClass           , weightOperatorList
    use :: Output_Analyses_Options           , only : enumerationOutputAnalysisCovarianceModelType
    implicit none
    type            (outputAnalysisScatterFunction1D             )                                          :: self
    type            (varying_string                              ), intent(in   )                           :: label                                       , comment                                      , &
         &                                                                                                     propertyLabel                               , propertyComment                              , &
         &                                                                                                     scatterLabel                                , scatterComment                               , &
         &                                                                                                     propertyUnits                               , scatterUnits                                 , &
         &                                                                                                     propertyQuantity                            , scatterQuantity
    double precision                                              , intent(in   )                           :: propertyUnitsInSI                           , scatterUnitsInSI
    double precision                                              , intent(in   )          , dimension(:  ) :: binCenter
    integer         (c_size_t                                    ), intent(in   )                           :: bufferCount
    double precision                                              , intent(in   )          , dimension(:,:) :: outputWeight
    logical                                                       , intent(in   ), optional                 :: likelihoodNormalize
    logical                                                       , intent(in   )                           :: propertyIsComoving                          , scatterIsComoving
    class           (nodePropertyExtractorClass                  ), intent(inout), target                   :: nodePropertyExtractor_                      , outputAnalysisWeightPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass         ), intent(inout), target                   :: outputAnalysisPropertyOperator_             , outputAnalysisPropertyUnoperator_            , &
         &                                                                                                     outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisWeightOperatorClass           ), intent(inout), target                   :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass     ), intent(inout), target                   :: outputAnalysisDistributionOperator_
    class           (galacticFilterClass                         ), intent(inout), target                   :: galacticFilter_
    class           (outputTimesClass                            ), intent(inout), target                   :: outputTimes_
    type            (enumerationOutputAnalysisCovarianceModelType), intent(in   )                           :: covarianceModel
    integer                                                       , intent(in   ), optional                 :: covarianceBinomialBinsPerDecade
    double precision                                              , intent(in   ), optional                 :: covarianceBinomialMassHaloMinimum           , covarianceBinomialMassHaloMaximum
    class           (outputAnalysisTargetDataClass               ), intent(in   ), optional                 :: targetData_
    type            (weightOperatorList                          ), pointer                                 :: weightOperatorWeight_                       , weightOperatorSquared_
    type            (propertyOperatorList                        ), pointer                                 :: propertyOperators_
    type            (outputAnalysisPropertyOperatorSequence      ), pointer                                 :: outputAnalysisWeightPropertyOperatorSquaring_
    type            (outputAnalysisPropertyOperatorSquare        ), pointer                                 :: outputAnalysisWeightPropertyOperatorSquare_
    !![
    <constructorAssign variables="label, comment, propertyLabel, propertyComment, propertyUnits, propertyQuantity, propertyIsComoving, propertyUnitsInSI, scatterLabel, scatterComment, scatterUnits, scatterQuantity, scatterIsComoving, scatterUnitsInSI, *nodePropertyExtractor_, *outputAnalysisWeightPropertyExtractor_, *outputAnalysisPropertyOperator_, *outputAnalysisWeightPropertyOperator_, *outputAnalysisPropertyUnoperator_, *outputAnalysisWeightOperator_, *outputAnalysisDistributionOperator_, *galacticFilter_, *outputTimes_, bufferCount, covarianceModel, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum"/>
    !!]

    ! Initialise the bundled target-data fields.  An explicit `targetData_` must be of the
    ! concrete `outputAnalysisTargetDataStandard` type (the only impl in the project today);
    ! without one we default-construct, matching the per-arg defaults the previous signature
    ! exposed.
    if (present(targetData_)) then
       select type (targetData_)
       type is (outputAnalysisTargetDataStandard)
          self%targetData_=targetData_
       class default
          call Error_Report('targetData_ must be of type outputAnalysisTargetDataStandard'//{introspection:location})
       end select
    else
       self%targetData_=outputAnalysisTargetDataStandard()
    end if
    ! Mirror the bundled target-data fields onto the outer object so the auto-built descriptor
    ! (which walks the type definition, not contained types) can reconstruct a parameter block
    ! that recreates this object.  Reshape the 2D covariance into the 1D form the parameter
    ! reader produces.
    self%xAxisLabel =self%targetData_%xAxisLabel
    self%yAxisLabel =self%targetData_%yAxisLabel
    self%targetLabel=self%targetData_%targetLabel
    self%xAxisIsLog =self%targetData_%xAxisIsLog
    self%yAxisIsLog =self%targetData_%yAxisIsLog
    if (self%targetData_%hasTarget()) then
       self%scatterValueTarget       =self%targetData_%valueTarget
       self%scatterCovarianceTarget1D=reshape(self%targetData_%covarianceTarget,[size(self%targetData_%covarianceTarget)])
    end if
    self%outputWeight=reshape(outputWeight,[size(outputWeight)])
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
       &amp;                       propertyQuantity                             , &amp;
       &amp;                       propertyIsComoving                           , &amp;
       &amp;                       propertyUnitsInSI                            , &amp;
       &amp;                       scatterLabel                                 , &amp;
       &amp;                       scatterComment                               , &amp;
       &amp;                       scatterUnits                                 , &amp;
       &amp;                       scatterQuantity                              , &amp;
       &amp;                       scatterIsComoving                            , &amp;
       &amp;                       scatterUnitsInSI                             , &amp;
       &amp;                       binCenter                                    , &amp;
       &amp;                       bufferCount                                  , &amp;
       &amp;                       outputWeight                                 , &amp;
       &amp;                       nodePropertyExtractor_                       , &amp;
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
       &amp;                       propertyQuantity                             , &amp;
       &amp;                       propertyIsComoving                           , &amp;
       &amp;                       propertyUnitsInSI                            , &amp;
       &amp;                       scatterLabel                                 , &amp;
       &amp;                       scatterComment                               , &amp;
       &amp;                       scatterUnits                                 , &amp;
       &amp;                       scatterQuantity                              , &amp;
       &amp;                       scatterIsComoving                            , &amp;
       &amp;                       scatterUnitsInSI                             , &amp;
       &amp;                       binCenter                                    , &amp;
       &amp;                       bufferCount                                  , &amp;
       &amp;                       outputWeight                                 , &amp;
       &amp;                       nodePropertyExtractor_                       , &amp;
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
    !!{RST
    Destructor for the ``outputAnalysisScatterFunction1D`` output analysis class.
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
    !!{RST
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
    !!{RST
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
    !!{RST
    Finalize analysis of a ``scatterFunction1D`` output analysis.
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
    !!{RST
    Implement a ``scatterFunction1D`` output analysis finalization.
    !!}
    use :: Output_HDF5   , only : outputFile
    use :: HDF5_Access   , only : hdf5Access
    use :: IO_HDF5       , only : hdf5Object
    use :: Units_MetaData, only : unitType
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
    call    analysisGroup%writeAttribute(     char(self%             comment    )         ,'description'                                                                                                      )
    call    analysisGroup%writeAttribute("function1D"                                     ,'type'                                                                                                             )
    call    analysisGroup%writeAttribute(     char(self%targetData_%xAxisLabel  )         ,'xAxisLabel'                                                                                                       )
    call    analysisGroup%writeAttribute(     char(self%targetData_%yAxisLabel  )         ,'yAxisLabel'                                                                                                       )
    call    analysisGroup%writeAttribute(          self%targetData_%xAxisIsLog            ,'xAxisIsLog'                                                                                                       )
    call    analysisGroup%writeAttribute(          self%targetData_%yAxisIsLog            ,'yAxisIsLog'                                                                                                       )
    call    analysisGroup%writeAttribute(     char(self%propertyLabel)                    ,'xDataset'                                                                                                         )
    call    analysisGroup%writeAttribute(     char(self% scatterLabel)                    ,'yDataset'                                                                                                         )
    call    analysisGroup%writeAttribute(     char(self% scatterLabel)//"Target"          ,'yDatasetTarget'                                                                                                   )
    call    analysisGroup%writeAttribute(     char(self% scatterLabel)//"Covariance"      ,'yCovariance'                                                                                                      )
    call    analysisGroup%writeAttribute(     char(self% scatterLabel)//"CovarianceTarget",'yCovarianceTarget'                                                                                                )
    ! Write computed datasets.
    call    analysisGroup%writeDataset  (          self%binCenter                         ,char(self%propertyLabel)                       ,char(self%propertyComment)                 ,datasetReturned=dataset)
    call    dataset      %writeAttribute(unitType(self%propertyUnitsInSI  ,description=     char(self%propertyUnits)      ,quantity=     char(self%propertyQuantity)       ,isComoving=self%propertyIsComoving),'units')
    call    analysisGroup%writeDataset  (          self%scatterValue                      ,char(self% scatterLabel)                       ,char(self% scatterComment)                 ,datasetReturned=dataset)
    call    dataset      %writeAttribute(unitType(self%scatterUnitsInSI   ,description=     char(self%scatterUnits )      ,quantity=     char(self%scatterQuantity )       ,isComoving=self%scatterIsComoving ),'units')
    call    analysisGroup%writeDataset  (          self%scatterCovariance                 ,char(self% scatterLabel)//"Covariance"         ,char(self% scatterComment)//" [covariance]",datasetReturned=dataset)
    call    dataset      %writeAttribute(unitType(self%scatterUnitsInSI**2,description="["//char(self%scatterUnits )//"]²",quantity="("//char(self%scatterQuantity )//")^2",isComoving=self%scatterIsComoving ),'units')
    ! If available, include the log-likelihood and target dataset.
    if (self%targetData_%hasTarget()) then
       call analysisGroup%writeAttribute(          self%logLikelihood()                       ,'logLikelihood'                                                                                                )
       call analysisGroup%writeAttribute(     char(self%targetData_%targetLabel)              ,'targetLabel'                                                                                                  )
       call analysisGroup%writeDataset  (          self%targetData_%valueTarget               ,char(self%    scatterLabel)//"Target"          ,char(self% scatterComment)                 ,datasetReturned=dataset)
       call dataset      %writeAttribute(unitType(self%scatterUnitsInSI   ,description=     char(self%scatterUnits )      ,quantity=     char(self%scatterQuantity )       ,isComoving=self%scatterIsComoving ),'units')
       call analysisGroup%writeDataset  (          self%targetData_%covarianceTarget          ,char(self%    scatterLabel)//"CovarianceTarget",char(self% scatterComment)//" [covariance]",datasetReturned=dataset)
       call dataset      %writeAttribute(unitType(self%scatterUnitsInSI**2,description="["//char(self%scatterUnits )//"]²",quantity="("//char(self%scatterQuantity )//")^2",isComoving=self%scatterIsComoving ),'units')
    end if
    !$ call hdf5Access%unset()
    return
  end subroutine scatterFunction1DFinalize

  subroutine scatterFunction1DResults(self,binCenter,scatterValue,scatterCovariance)
    !!{RST
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
    !!{RST
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
    if (self%targetData_%hasTarget()) then
       ! Finalize analysis.
       call self%finalizeAnalysis()
       ! Allocate workspaces.
       allocate(scatterCovarianceCombined(size(self%binCenter),size(self%binCenter)))
       allocate(scatterValueDifference   (size(self%binCenter)                     ))
       ! Find combined covariance and difference between model and target.
       scatterValueDifference   =+self%scatterValue                  &
            &                    -self%targetData_%valueTarget
       scatterCovarianceCombined=+self%scatterCovariance             &
            &                    +self%targetData_%covarianceTarget
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
