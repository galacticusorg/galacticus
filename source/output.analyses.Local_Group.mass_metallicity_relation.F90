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

  !!{
  Implements an output analysis class that computes mass-metallicity relations for Local Group satellite
  galaxies.
  !!}

  !![
  <outputAnalysis name="outputAnalysisLocalGroupMassMetallicityRelation">
   <description>An output analysis class for Local Group satellite galaxy mass-metallicity relations.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisLocalGroupMassMetallicityRelation
     !!{
     An output analysis class for Local Group satellite galaxy mass-metallicity relations.
     !!}
     private
     class           (outputAnalysisClass        ), pointer                   :: outputAnalysis_                                 => null()
     class           (outputTimesClass           ), pointer                   :: outputTimes_                                    => null()
     double precision                             , allocatable, dimension(:) :: randomErrorPolynomialCoefficient                         , systematicErrorPolynomialCoefficient, &
          &                                                                      metallicitySystematicErrorPolynomialCoefficient
     integer                                                                  :: covarianceBinomialBinsPerDecade
     double precision                                                         :: covarianceBinomialMassHaloMinimum                        , covarianceBinomialMassHaloMaximum   , &
          &                                                                      randomErrorMinimum                                       , randomErrorMaximum
     type            (enumerationPositionTypeType)                            :: positionType
   contains
     final     ::                  localGroupMassMetallicityRelationDestructor
     procedure :: analyze       => localGroupMassMetallicityRelationAnalyze
     procedure :: finalize      => localGroupMassMetallicityRelationFinalize
     procedure :: reduce        => localGroupMassMetallicityRelationReduce
     procedure :: logLikelihood => localGroupMassMetallicityRelationLogLikelihood
  end type outputAnalysisLocalGroupMassMetallicityRelation

  interface outputAnalysisLocalGroupMassMetallicityRelation
     !!{
     Constructors for the \refClass{outputAnalysisLocalGroupMassMetallicityRelation} output analysis class.
     !!}
     module procedure localGroupMassMetallicityRelationConstructorParameters
     module procedure localGroupMassMetallicityRelationConstructorInternal
  end interface outputAnalysisLocalGroupMassMetallicityRelation

contains

  function localGroupMassMetallicityRelationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisLocalGroupMassMetallicityRelation} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters            , only : inputParameter               , inputParameters
    use :: Output_Times                , only : outputTimes                  , outputTimesClass
    use :: Galactic_Filters            , only : enumerationPositionTypeEncode
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    type            (outputAnalysisLocalGroupMassMetallicityRelation)                              :: self
    type            (inputParameters                                ), intent(inout)               :: parameters
    class           (outputTimesClass                               ), pointer                     :: outputTimes_
    double precision                                                 , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient               , systematicErrorPolynomialCoefficient, &
         &                                                                                            metallicitySystematicErrorPolynomialCoefficient
    integer                                                                                        :: covarianceBinomialBinsPerDecade
    double precision                                                                               :: covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum   , &
         &                                                                                            randomErrorMinimum                             , randomErrorMaximum
    type            (varying_string                                 )                              :: positionType

    ! Check and read parameters.
    allocate(metallicitySystematicErrorPolynomialCoefficient(max(1,parameters%count('metallicitySystematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(           systematicErrorPolynomialCoefficient(max(1,parameters%count(           'systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(               randomErrorPolynomialCoefficient(max(1,parameters%count(               'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    !![
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.1d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.1d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.07d0]</defaultValue>
      <description>The coefficients of the random error polynomial for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>metallicitySystematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>metallicitySystematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the metallicity systematic error polynomial.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <variable>covarianceBinomialBinsPerDecade</variable>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing Local Group stellar mass function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMinimum</variable>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing Local Group stellar mass function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMaximum</variable>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing Local Group stellar mass function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>positionType</name>
      <source>parameters</source>
      <defaultValue>var_str('orbital')</defaultValue>
      <description>The type of position to use in survey geometry filters.</description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=outputAnalysisLocalGroupMassMetallicityRelation(outputTimes_,enumerationPositionTypeEncode(positionType,includesPrefix=.false.),randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,metallicitySystematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function localGroupMassMetallicityRelationConstructorParameters

  function localGroupMassMetallicityRelationConstructorInternal(outputTimes_,positionType,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,metallicitySystematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLocalGroupMassMetallicityRelation} output analysis class for internal use.
    !!}
    use :: Atomic_Data                             , only : Atomic_Number
    use :: Abundances_Structure                    , only : Abundances_Index_From_Name
    use :: Galactic_Filters                        , only : filterList                                            , galacticFilterAll                      , galacticFilterHaloNotIsolated         , galacticFilterHostMassRange                    , &
          &                                                 galacticFilterSurveyGeometry                          , enumerationPositionTypeType
    use :: Geometry_Surveys                        , only : surveyGeometryFullSky
    use :: Interface_Local_Group_DB                , only : comparisonEquals                                      , comparisonLessThan                     , localGroupDB                          , setOperatorIntersection                        , &
          &                                                 setOperatorRelativeComplement                         , setOperatorUnion                       , attributeUncertainty
    use :: Node_Property_Extractors                , only : nodePropertyExtractorMassStellar                      , nodePropertyExtractorMetallicityStellar
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Numerical_Constants_Astronomical        , only : massSolar
    use :: Numerical_Ranges                        , only : Make_Range                                            , rangeTypeLinear
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerIdentity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10               , outputAnalysisPropertyOperatorLog10    , outputAnalysisPropertyOperatorSequence, outputAnalysisPropertyOperatorSystmtcPolynomial, &
          &                                                 outputAnalysisPropertyOperatorMetallicitySolarRelative, propertyOperatorList
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorSubsampling
    use :: Output_Times                            , only : outputTimesClass
    implicit none
    type            (outputAnalysisLocalGroupMassMetallicityRelation       )                                :: self
    integer                                                                 , intent(in   )                 :: covarianceBinomialBinsPerDecade
    double precision                                                        , intent(in   )                 :: covarianceBinomialMassHaloMinimum                          , covarianceBinomialMassHaloMaximum                                 , &
         &                                                                                                     randomErrorMinimum                                         , randomErrorMaximum
    double precision                                                        , intent(in   ), dimension(:  ) :: randomErrorPolynomialCoefficient                           , systematicErrorPolynomialCoefficient                              , &
         &                                                                                                     metallicitySystematicErrorPolynomialCoefficient
    type            (enumerationPositionTypeType                           ), intent(in   )                 :: positionType
    class           (outputTimesClass                                      ), intent(inout), target         :: outputTimes_
    type            (nodePropertyExtractorMassStellar                      )               , pointer        :: nodePropertyExtractor_
    type            (nodePropertyExtractorMetallicityStellar               )               , pointer        :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorMetallicitySolarRelative)               , pointer        :: outputAnalysisWeightPropertyOperatorMetallicity_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial       )               , pointer        :: outputAnalysisPropertyOperatorSystmtcPolynomial_           , outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisPropertyOperatorLog10                   )               , pointer        :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorSequence                )               , pointer        :: outputAnalysisPropertyOperator_                            , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10               )               , pointer        :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorSubsampling               )               , pointer        :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml   )               , pointer        :: outputAnalysisDistributionOperator_
    type            (surveyGeometryFullSky                                 )               , pointer        :: surveyGeometry_
    type            (galacticFilterHaloNotIsolated                         )               , pointer        :: galacticFilterHaloNotIsolated_
    type            (galacticFilterHostMassRange                           )               , pointer        :: galacticFilterHostMassRange_
    type            (galacticFilterSurveyGeometry                          )               , pointer        :: galacticFilterSurveyGeometry_
    type            (galacticFilterAll                                     )               , pointer        :: galacticFilter_
    type            (filterList                                            )               , pointer        :: filters_
    type            (propertyOperatorList                                  )               , pointer        :: operators_                                                  , weightPropertyOperators_
    integer                                                                 , allocatable  , dimension(:  ) :: countTarget
    logical                                                                 , allocatable  , dimension(:  ) :: isPresentMasses                                             , isPresentMetallicities                                           , &
         &                                                                                                     isPresentMetallicitiesUncertainties
    double precision                                                        , allocatable  , dimension(:  ) :: masses                                                      , massesTarget                                                     , &
         &                                                                                                     metallicitiesTarget                                         , functionValueTarget                                              , &
         &                                                                                                     metallicitiesUncertaintiesTarget                            , functionVarianceMeasurementTarget                                , &
         &                                                                                                     functionVarianceSampleTarget                                , functionValueTargetNonZero                                       , &
         &                                                                                                     massesNonZero
    double precision                                                        , allocatable  , dimension(:,:) :: outputWeight                                                , functionCovarianceTarget                                         , &
         &                                                                                                     functionCovarianceTargetNonZero
    double precision                                                        , parameter                     :: bufferWidthLogarithmic                          =+3.0d+0    , errorZeroPoint                                       =10.00d0    , &
         &                                                                                                     metallicityErrorPolynomialZeroPoint             =+0.0d+0
    integer         (c_size_t                                              ), parameter                     :: binCount                                        = 7_c_size_t, bufferCountMinimum                                   = 5_c_size_t
    double precision                                                        , parameter                     :: massMinimum                                     =+1.0d+3    , massMaximum                                          = 1.00d9    , &
         &                                                                                                     radiusOuter                                     =+3.0d-1    , metallicityUncertaintyDefault                        = 0.25d0
    logical                                                                 , parameter                     :: likelihoodNormalize                             =.false.
    integer         (c_size_t                                              )                                :: i                                                           , j                                                                , &
         &                                                                                                     bufferCount                                                 , binCountNonZero
    type            (localGroupDB                                          )                                :: localGroupDB_
    double precision                                                                                        :: massesWidthBin
    !![
    <constructorAssign variables="*outputTimes_, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, metallicitySystematicErrorPolynomialCoefficient, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum, randomErrorMinimum, randomErrorMaximum, positionType"/>
    !!]
    
    ! Construct mass bins.
    allocate(masses(binCount))
    masses        =Make_Range(log10(massMinimum),log10(massMaximum),int(binCount),rangeTypeLinear)
    massesWidthBin=+masses(2) &
         &         -masses(1)
    ! Construct the target distribution.
    !! Select galaxies within 300 kpc of the Milky Way, excluding the Milky Way itself, then retrieve stellar masses and
    !! metallicities for the selected galaxies.
    localGroupDB_=localGroupDB()
    call localGroupDB_%selectAll  (                                                                                                                       )
    call localGroupDB_%select     ('classification'    ,var_str('galaxy'    )           ,comparisonEquals                   ,setOperatorIntersection      )
    call localGroupDB_%select     ('distanceMilkyWay'  ,         radiusOuter            ,comparisonLessThan                 ,setOperatorIntersection      )
    call localGroupDB_%select     ('name'              ,var_str('The Galaxy')           ,comparisonEquals                   ,setOperatorRelativeComplement)
    call localGroupDB_%getProperty('massStellar'       ,massesTarget                    ,isPresentMasses                                                  )
    call localGroupDB_%getProperty('metallicityStellar',metallicitiesTarget             ,isPresentMetallicities                                           )
    call localGroupDB_%getProperty('metallicityStellar',metallicitiesUncertaintiesTarget,isPresentMetallicitiesUncertainties,attributeUncertainty         )
    ! Where metallicities have no uncertainty, assign a default one.
    where (isPresentMetallicities .and. .not.isPresentMetallicitiesUncertainties)
       metallicitiesUncertaintiesTarget=metallicityUncertaintyDefault
    end where
    ! Build the target dataset.
    allocate(countTarget                      (binCount         ))
    allocate(functionValueTarget              (binCount         ))
    allocate(functionVarianceMeasurementTarget(binCount         ))
    allocate(functionVarianceSampleTarget     (binCount         ))
    allocate(functionCovarianceTarget         (binCount,binCount))
    countTarget                      =0
    functionValueTarget              =0.0d0
    functionVarianceMeasurementTarget=0.0d0
    functionVarianceSampleTarget     =0.0d0
    functionCovarianceTarget         =0.0d0
    do i=1,size(massesTarget)
       if     (                                             &
            &   .not.isPresentMasses                    (i) &
            &  .or.                                         &
            &   .not.isPresentMetallicities             (i) &
            & ) cycle
       j=int((log10(massesTarget(i))-log10(massMinimum))/(masses(2)-masses(1))+0.5d0,kind=c_size_t)+1_c_size_t
       if (j > 0 .and. j <= binCount) then
          functionValueTarget              (j)=functionValueTarget              (j)+metallicitiesTarget             (i)
          functionVarianceMeasurementTarget(j)=functionVarianceMeasurementTarget(j)+metallicitiesUncertaintiesTarget(i)**2
          functionVarianceSampleTarget     (j)=functionVarianceSampleTarget     (j)+metallicitiesTarget             (i)**2
          countTarget                      (j)=countTarget                      (j)+1
       end if
    end do
    do j=1,binCount
       if (countTarget(j) > 0) then
          functionValueTarget              (j  )=  functionValueTarget              (j)/dble(countTarget(j))
          functionVarianceMeasurementTarget(j  )=  functionVarianceMeasurementTarget(j)/dble(countTarget(j))**2
          functionVarianceSampleTarget     (j  )=(                                                              &
               &                                  +functionVarianceSampleTarget     (j)/dble(countTarget(j))    &
               &                                  -functionValueTarget              (j)                     **2 &
               &                                 )                                     /dble(countTarget(j))
          functionCovarianceTarget         (j,j)=+ functionVarianceMeasurementTarget(j)                         &
               &                                 + functionVarianceSampleTarget     (j)
       end if
    end do
    ! Find non-empty bins.
    binCountNonZero=count(countTarget > 0)
    allocate(massesNonZero                  (binCountNonZero                ))
    allocate(functionValueTargetNonZero     (binCountNonZero                ))
    allocate(functionCovarianceTargetNonZero(binCountNonZero,binCountNonZero))
    functionValueTargetNonZero     =0.0d0
    functionCovarianceTargetNonZero=0.0d0
    j                              =0
    do i=1,binCount
       if (countTarget(i) == 0) cycle
       j                                   =j                            +1
       massesNonZero                  (j  )=masses                  (i  )
       functionValueTargetNonZero     (j  )=functionValueTarget     (i  )
       functionCovarianceTargetNonZero(j,j)=functionCovarianceTarget(i,i)
    end do
    ! Create a stellar mass property extractor.
    allocate(nodePropertyExtractor_                                )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                                 constructor="nodePropertyExtractorMassStellar               (                                                                                    )"/>
    !!]
    ! Create a stellar metallicity weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"                 constructor="nodePropertyExtractorMetallicityStellar         (Abundances_Index_From_Name('Fe')                                                   )"/>
    !!]
    ! Build a metallicity weight property operator.
    allocate(outputAnalysisWeightPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial (metallicityErrorPolynomialZeroPoint,metallicitySystematicErrorPolynomialCoefficient)"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorMetallicity_      )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorMetallicity_">
     <constructor>
      outputAnalysisPropertyOperatorMetallicitySolarRelative(                               &amp;
        &amp;                                                Atomic_Number(shortLabel="Fe") &amp;
        &amp;                                               )
     </constructor>
    </referenceConstruct>
    !!]
    allocate(weightPropertyOperators_                              )
    allocate(weightPropertyOperators_%next                         )
    weightPropertyOperators_     %operator_ => outputAnalysisWeightPropertyOperatorMetallicity_
    weightPropertyOperators_%next%operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"            constructor="outputAnalysisPropertyOperatorSequence         (weightPropertyOperators_                           )"/>
    !!]
    ! Create property operators and unoperators to perform conversion to/from logarithmic mass.
    allocate(outputAnalysisPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                   )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorZeroPoint,systematicErrorPolynomialCoefficient)"/>
    !!]
    allocate(operators_     )
    allocate(operators_%next)
    operators_     %operator_ => outputAnalysisPropertyOperatorLog10_
    operators_%next%operator_ => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence         (operators_                                         )"/>
    !!]
    allocate(outputAnalysisPropertyUnoperator_               )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                constructor="outputAnalysisPropertyOperatorAntiLog10        (                                                   )"/>
    !!]
    ! Create a subsampling weight operator.
    allocate(outputAnalysisWeightOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                    constructor="outputAnalysisWeightOperatorSubsampling        (                                                   )"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_             )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
    <constructor>
    outputAnalysisDistributionOperatorRandomErrorPlynml (                                  &amp;
         &amp;                                           randomErrorMinimum              , &amp;
         &amp;                                           randomErrorMaximum              , &amp;
         &amp;                                           errorZeroPoint                  , &amp;
         &amp;                                           randomErrorPolynomialCoefficient  &amp;
         &amp;                                          )
     </constructor>
    </referenceConstruct>
    !!]
    !!]
    ! Build survey geometry.
    !! Use a full sky geometry out to the outer radius considered (300 kpc).
    allocate(surveyGeometry_               )
    !![
    <referenceConstruct object="surveyGeometry_"                constructor="surveyGeometryFullSky            (distanceMinimum=0.0d0,distanceMaximum=radiusOuter)"/>
    !!]
    ! Build filters which select satellites in a specified range of host halo mass.
    allocate(galacticFilterHaloNotIsolated_)
    !![
    <referenceConstruct object="galacticFilterHaloNotIsolated_" constructor="galacticFilterHaloNotIsolated    (                                                 )"/>
    !!]
    allocate(galacticFilterSurveyGeometry_ )
    !![
    <referenceConstruct object="galacticFilterSurveyGeometry_"  constructor="galacticFilterSurveyGeometry     (positionType         ,surveyGeometry_            )"/>
    !!]
    allocate(galacticFilterHostMassRange_  )
    !![
    <referenceConstruct object="galacticFilterHostMassRange_">
     <constructor>
      galacticFilterHostMassRange(                      &amp;
        &amp;                     massMinimum =1.00d12, &amp;
        &amp;                     massMaximum =2.00d12, &amp;
        &amp;                     useFinalHost=.true.   &amp;
        &amp;                    )
     </constructor>
    </referenceConstruct>
    !!]
    allocate(filters_          )
    allocate(filters_%next     )
    allocate(filters_%next%next)
    filters_          %filter_ => galacticFilterHaloNotIsolated_ 
    filters_%next     %filter_ => galacticFilterHostMassRange_
    filters_%next%next%filter_ => galacticFilterSurveyGeometry_
    allocate(galacticFilter_)
    !![
    <referenceConstruct object="galacticFilter_" constructor="galacticFilterAll(filters_)"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(binCount,outputTimes_%count()))
    do i=1_c_size_t,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(i),0.0d0,absTol=1.0d-6)) then
          outputWeight(:,i)=1.0d0
       else
          outputWeight(:,i)=0.0d0
       end if
    end do
    if (any(sum(outputWeight,dim=2) /= 1.0d0)) call Error_Report('output weights do not equal unity'//{introspection:location})
    ! Compute the number of buffer bins to add to either side of the mass function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    bufferCount=max(int(bufferWidthLogarithmic/(masses(2)-masses(1)),kind=c_size_t)+1_c_size_t,bufferCountMinimum)
    ! Construct the analysis objects.
    allocate(outputAnalysisMeanFunction1D :: self%outputAnalysis_)
    select type (outputAnalysis_ => self%outputAnalysis_)
    type is (outputAnalysisMeanFunction1D)
       !![
       <referenceConstruct isResult="yes" object="outputAnalysis_">
	 <constructor>
	   outputAnalysisMeanFunction1D(                                                                &amp;
	   &amp;                        var_str('localGroupMassMetallicityRelation'                  ), &amp;
	   &amp;                        var_str('Mass-metallicity relation of Local Group satellites'), &amp;
	   &amp;                        var_str('massStellar'                                        ), &amp;
	   &amp;                        var_str('Stellar mass at the bin center'                     ), &amp;
	   &amp;                        var_str('M☉'                                                 ), &amp;
	   &amp;                        massSolar                                                     , &amp;
	   &amp;                        var_str('metallicityMean'                                    ), &amp;
           &amp;                        var_str('Mean stellar metallicity; ⟨[Fe/H]⟩'                 ), &amp;
           &amp;                        var_str('dimensionless'                                      ), &amp;
           &amp;                        0.0d0                                                         , &amp;
	   &amp;                        massesNonZero                                                 , &amp;
	   &amp;                        bufferCount                                                   , &amp;
	   &amp;                        outputWeight                                                  , &amp;
	   &amp;                        nodePropertyExtractor_                                        , &amp;
           &amp;                        outputAnalysisWeightPropertyExtractor_                        , &amp;
	   &amp;                        outputAnalysisPropertyOperator_                               , &amp;
           &amp;                        outputAnalysisWeightPropertyOperator_                         , &amp;
	   &amp;                        outputAnalysisPropertyUnoperator_                             , &amp;
	   &amp;                        outputAnalysisWeightOperator_                                 , &amp;
	   &amp;                        outputAnalysisDistributionOperator_                           , &amp;
	   &amp;                        galacticFilter_                                               , &amp;
	   &amp;                        outputTimes_                                                  , &amp;
	   &amp;                        outputAnalysisCovarianceModelBinomial                         , &amp;
	   &amp;                        covarianceBinomialBinsPerDecade                               , &amp;
	   &amp;                        covarianceBinomialMassHaloMinimum                             , &amp;
	   &amp;                        covarianceBinomialMassHaloMaximum                             , &amp;
           &amp;                        likelihoodNormalize                                           , &amp;
           &amp;                        var_str('$M_\star/\mathrm{M}_\odot$'                         ), &amp;
           &amp;                        var_str('$[\mathrm{Fe}/\mathrm{H}]$'                         ), &amp;
           &amp;                        .true.                                                        , &amp;
           &amp;                        .false.                                                       , &amp;
           &amp;                        var_str('Galacticus compilation'                             ), &amp;
           &amp;                        functionValueTargetNonZero                                    , &amp;
           &amp;                        functionCovarianceTargetNonZero                               , &amp;
	   &amp;                        massesWidthBin                                                  &amp; 
	   &amp;                       )
	 </constructor>
       </referenceConstruct>
       !!]         
    end select
    !![
    <objectDestructor name="nodePropertyExtractor_"                                />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"                />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisWeightPropertyOperatorMetallicity_"      />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisPropertyOperator_"                       />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"                  />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"      />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"                     />
    <objectDestructor name="outputAnalysisWeightOperator_"                         />
    <objectDestructor name="outputAnalysisDistributionOperator_"                   />
    <objectDestructor name="galacticFilterHaloNotIsolated_"                        />
    <objectDestructor name="galacticFilterHostMassRange_"                          />
    <objectDestructor name="galacticFilterSurveyGeometry_"                         />
    <objectDestructor name="galacticFilter_"                                       />
    <objectDestructor name="surveyGeometry_"                                       />
    !!]
    nullify(filters_                )
    nullify(operators_              )
    nullify(weightPropertyOperators_)
    return
  end function localGroupMassMetallicityRelationConstructorInternal

  subroutine localGroupMassMetallicityRelationDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisLocalGroupMassMetallicityRelation} output analysis class.
    !!}
    implicit none
    type(outputAnalysisLocalGroupMassMetallicityRelation), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysis_"/>
    <objectDestructor name="self%outputTimes_"   />
    !!]
    return
  end subroutine localGroupMassMetallicityRelationDestructor

  subroutine localGroupMassMetallicityRelationAnalyze(self,node,iOutput)
    !!{
    Implement a {\normalfont \ttfamily localGroupMassMetallicityRelation} output analysis.
    !!}
    implicit none
    class  (outputAnalysisLocalGroupMassMetallicityRelation), intent(inout) :: self
    type   (treeNode                                       ), intent(inout) :: node
    integer(c_size_t                                       ), intent(in   ) :: iOutput

    call self%outputAnalysis_%analyze(node,iOutput)
    return
  end subroutine localGroupMassMetallicityRelationAnalyze

  subroutine localGroupMassMetallicityRelationReduce(self,reduced)
    !!{
    Implement a {\normalfont \ttfamily localGroupMassMetallicityRelation} output analysis reduction.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisLocalGroupMassMetallicityRelation), intent(inout) :: self
    class(outputAnalysisClass                            ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisLocalGroupMassMetallicityRelation)
       call self%outputAnalysis_%reduce(reduced%outputAnalysis_)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine localGroupMassMetallicityRelationReduce

  subroutine localGroupMassMetallicityRelationFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily localGroupMassMetallicityRelation} output analysis finalization.
    !!}
    implicit none
    class(outputAnalysisLocalGroupMassMetallicityRelation), intent(inout)           :: self
    type (varying_string                                 ), intent(in   ), optional :: groupName

    call self%outputAnalysis_%finalize(groupName)
    return
  end subroutine localGroupMassMetallicityRelationFinalize

  double precision function localGroupMassMetallicityRelationLogLikelihood(self)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily localGroupMassMetallicityRelation} output analysis.
    !!}
    implicit none
    class(outputAnalysisLocalGroupMassMetallicityRelation), intent(inout) :: self

    localGroupMassMetallicityRelationLogLikelihood=self%outputAnalysis_%logLikelihood()
    return
  end function localGroupMassMetallicityRelationLogLikelihood
