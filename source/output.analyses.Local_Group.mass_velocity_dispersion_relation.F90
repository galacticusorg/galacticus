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
  Implements an output analysis class that computes mass-velocity dispersion relations for Local Group
  satellite galaxies.
  !!}

  !![
  <outputAnalysis name="outputAnalysisLocalGroupMassVelocityDispersionRelation">
   <description>An output analysis class for Local Group satellite galaxy mass-velocity dispersion relations.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisLocalGroupMassVelocityDispersionRelation
     !!{
     An output analysis class for Local Group satellite galaxy mass-velocity dispersion relations.
     !!}
     private
     class           (outputAnalysisClass        ), pointer                     :: outputAnalysis_                                        => null()
     class           (outputTimesClass           ), pointer                     :: outputTimes_                                           => null()
     class           (darkMatterHaloScaleClass   ), pointer                     :: darkMatterHaloScale_                                   => null()
     double precision                             , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient                                , systematicErrorPolynomialCoefficient, &
          &                                                                        velocityDispersionSystematicErrorPolynomialCoefficient
     integer                                                                    :: covarianceBinomialBinsPerDecade
     double precision                                                           :: covarianceBinomialMassHaloMinimum                               , covarianceBinomialMassHaloMaximum   , &
          &                                                                        randomErrorMinimum                                              , randomErrorMaximum
     type            (enumerationPositionTypeType)                              :: positionType
   contains
     final     ::                  localGroupMassVelocityDispersionRelationDestructor
     procedure :: analyze       => localGroupMassVelocityDispersionRelationAnalyze
     procedure :: finalize      => localGroupMassVelocityDispersionRelationFinalize
     procedure :: reduce        => localGroupMassVelocityDispersionRelationReduce
     procedure :: logLikelihood => localGroupMassVelocityDispersionRelationLogLikelihood
  end type outputAnalysisLocalGroupMassVelocityDispersionRelation

  interface outputAnalysisLocalGroupMassVelocityDispersionRelation
     !!{
     Constructors for the \refClass{outputAnalysisLocalGroupMassVelocityDispersionRelation} output analysis class.
     !!}
     module procedure localGroupMassVelocityDispersionRelationConstructorParameters
     module procedure localGroupMassVelocityDispersionRelationConstructorInternal
  end interface outputAnalysisLocalGroupMassVelocityDispersionRelation

contains

  function localGroupMassVelocityDispersionRelationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisLocalGroupMassVelocityDispersionRelation} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters            , only : inputParameter               , inputParameters
    use :: Output_Times                , only : outputTimes                  , outputTimesClass
    use :: Galactic_Filters            , only : enumerationPositionTypeEncode
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    type            (outputAnalysisLocalGroupMassVelocityDispersionRelation)                              :: self
    type            (inputParameters                                       ), intent(inout)               :: parameters
    class           (outputTimesClass                                      ), pointer                     :: outputTimes_
    class           (darkMatterHaloScaleClass                              ), pointer                     :: darkMatterHaloScale_
    double precision                                                        , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient                      , systematicErrorPolynomialCoefficient, &
         &                                                                                                   velocityDispersionSystematicErrorPolynomialCoefficient
    integer                                                                                               :: covarianceBinomialBinsPerDecade
    double precision                                                                                      :: covarianceBinomialMassHaloMinimum                     , covarianceBinomialMassHaloMaximum   , &
         &                                                                                                   randomErrorMinimum                                    , randomErrorMaximum
    type            (varying_string                                        )                              :: positionType

    ! Check and read parameters.
    allocate(velocityDispersionSystematicErrorPolynomialCoefficient(max(1,parameters%count('velocityDispersionSystematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(    systematicErrorPolynomialCoefficient              (max(1,parameters%count(                  'systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(        randomErrorPolynomialCoefficient              (max(1,parameters%count(                      'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
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
      <name>velocityDispersionSystematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>velocityDispersionSystematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the velocity dispersion systematic error polynomial.</description>
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
      <description>The number of bins per decade of halo mass to use when constructing Local Group mass-velocity dispersion covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMinimum</variable>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing Local Group mass-velocity dispersion covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMaximum</variable>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing Local Group mass-velocity dispersion covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>positionType</name>
      <source>parameters</source>
      <defaultValue>var_str('orbital')</defaultValue>
      <description>The type of position to use in survey geometry filters.</description>
    </inputParameter>
    <objectBuilder class="outputTimes"         name="outputTimes_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=outputAnalysisLocalGroupMassVelocityDispersionRelation(outputTimes_,darkMatterHaloScale_,enumerationPositionTypeEncode(positionType,includesPrefix=.false.),randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,velocityDispersionSystematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"        />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function localGroupMassVelocityDispersionRelationConstructorParameters

  function localGroupMassVelocityDispersionRelationConstructorInternal(outputTimes_,darkMatterHaloScale_,positionType,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,velocityDispersionSystematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLocalGroupMassVelocityDispersionRelation} output analysis class for internal use.
    !!}
    use :: Galactic_Filters                        , only : filterList                                          , galacticFilterAll                         , galacticFilterHaloNotIsolated         , galacticFilterHostMassRange                    , &
          &                                                 galacticFilterSurveyGeometry                        , galacticFilterHighPass                    , enumerationPositionTypeType
    use :: Geometry_Surveys                        , only : surveyGeometryFullSky
    use :: Galactic_Structure_Options              , only : componentTypeAll                                    , massTypeGalactic
    use :: Interface_Local_Group_DB                , only : comparisonEquals                                    , comparisonLessThan                        , localGroupDB                          , setOperatorIntersection                        , &
          &                                                 setOperatorRelativeComplement                       , setOperatorUnion                          , attributeUncertainty
    use :: Node_Property_Extractors                , only : nodePropertyExtractorMassStellar                    , nodePropertyExtractorVelocityDispersion   , nodePropertyExtractorScalarizer
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Numerical_Constants_Astronomical        , only : massSolar
    use :: Numerical_Ranges                        , only : Make_Range                                          , rangeTypeLinear
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerIdentity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10             , outputAnalysisPropertyOperatorLog10       , outputAnalysisPropertyOperatorSequence, outputAnalysisPropertyOperatorSystmtcPolynomial, &
          &                                                 propertyOperatorList
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorSubsampling
    use :: Output_Times                            , only : outputTimesClass
    implicit none
    type            (outputAnalysisLocalGroupMassVelocityDispersionRelation)                                :: self
    integer                                                                 , intent(in   )                 :: covarianceBinomialBinsPerDecade
    double precision                                                        , intent(in   )                 :: covarianceBinomialMassHaloMinimum                          , covarianceBinomialMassHaloMaximum                                 , &
         &                                                                                                     randomErrorMinimum                                         , randomErrorMaximum
    double precision                                                        , intent(in   ), dimension(:  ) :: randomErrorPolynomialCoefficient                           , systematicErrorPolynomialCoefficient                              , &
         &                                                                                                     velocityDispersionSystematicErrorPolynomialCoefficient
    type            (enumerationPositionTypeType                           ), intent(in   )                 :: positionType
    class           (outputTimesClass                                      ), intent(inout), target         :: outputTimes_
    class           (darkMatterHaloScaleClass                              ), intent(in   ), target         :: darkMatterHaloScale_
    type            (nodePropertyExtractorMassStellar                      )               , pointer        :: nodePropertyExtractor_
    type            (nodePropertyExtractorScalarizer                       )               , pointer        :: outputAnalysisWeightPropertyScalarizer_
    type            (nodePropertyExtractorVelocityDispersion               )               , pointer        :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial       )               , pointer        :: outputAnalysisPropertyOperatorSystmtcPolynomial_           , outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisPropertyOperatorLog10                   )               , pointer        :: outputAnalysisPropertyOperatorLog10_                       , outputAnalysisWeightPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorSequence                )               , pointer        :: outputAnalysisPropertyOperator_                            , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10               )               , pointer        :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorSubsampling               )               , pointer        :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml   )               , pointer        :: outputAnalysisDistributionOperator_
    type            (surveyGeometryFullSky                                 )               , pointer        :: surveyGeometry_
    type            (galacticFilterHaloNotIsolated                         )               , pointer        :: galacticFilterHaloNotIsolated_
    type            (galacticFilterHostMassRange                           )               , pointer        :: galacticFilterHostMassRange_
    type            (galacticFilterSurveyGeometry                          )               , pointer        :: galacticFilterSurveyGeometry_
    type            (galacticFilterHighPass                                )               , pointer        :: galacticFilterHighPass_
    type            (galacticFilterAll                                     )               , pointer        :: galacticFilter_
    type            (filterList                                            )               , pointer        :: filters_
    type            (propertyOperatorList                                  )               , pointer        :: operators_                                                  , weightPropertyOperators_
    integer                                                                 , allocatable  , dimension(:  ) :: countTarget
    logical                                                                 , allocatable  , dimension(:  ) :: isPresentMasses                                             , isPresentVelocityDispersions                                     , &
         &                                                                                                     isPresentVelocityDispersionsUncertainties
    double precision                                                        , allocatable  , dimension(:  ) :: masses                                                      , massesTarget                                                     , &
         &                                                                                                     velocityDispersionsTarget                                   , functionValueTarget                                              , &
         &                                                                                                     velocityDispersionsUncertaintiesTarget                      , functionVarianceMeasurementTarget                                , &
         &                                                                                                     functionVarianceSampleTarget                                , functionValueTargetNonZero                                       , &
         &                                                                                                     massesNonZero
    double precision                                                        , allocatable  , dimension(:,:) :: outputWeight                                                , functionCovarianceTarget                                         , &
         &                                                                                                     functionCovarianceTargetNonZero
    double precision                                                        , parameter                     :: bufferWidthLogarithmic                          =+3.0d+0    , errorZeroPoint                                       =10.00d+0   , &
         &                                                                                                     velocityDispersionErrorPolynomialZeroPoint      =+1.0d+0    , covarianceLarge                                      = 1.00d+0   , &
         &                                                                                                     velocityDispersionUndefined                     =+1.3d+0    , toleranceRelative                                    = 1.00d-2   , &
         &                                                                                                     velocityDispersionMinimum                       =+1.0d-3
    integer         (c_size_t                                              ), parameter                     :: binCount                                        = 7_c_size_t, bufferCountMinimum                                   = 5_c_size_t
    double precision                                                        , parameter                     :: massMinimum                                     =+1.0d+3    , massMaximum                                          = 1.00d+9   , &
         &                                                                                                     radiusOuter                                     =+3.0d-1    , sizeUncertaintyDefault                               = 0.25d+0
    logical                                                                 , parameter                     :: likelihoodNormalize                             =.false.
    integer         (c_size_t                                              )                                :: i                                                           , j                                                                , &
         &                                                                                                     bufferCount                                                 , binCountNonZero
    type            (localGroupDB                                          )                                :: localGroupDB_
    double precision                                                                                        :: massesWidthBin
    type            (varying_string                                        )               , dimension(1)   :: radiusSpecifier
    !![
    <constructorAssign variables="*outputTimes_, *darkMatterHaloScale_, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, velocityDispersionSystematicErrorPolynomialCoefficient, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum, randomErrorMinimum, randomErrorMaximum, positionType"/>
    !!]

    ! Construct mass bins.
    allocate(masses(binCount))
    masses        =Make_Range(log10(massMinimum),log10(massMaximum),int(binCount),rangeTypeLinear)
    massesWidthBin=+masses(2) &
         &         -masses(1)
    ! Construct the target distribution.
    !! Select galaxies within 300 kpc of the Milky Way, excluding the Milky Way itself, then retrieve stellar masses and
    !! size for the selected galaxies.
    localGroupDB_=localGroupDB()
    call localGroupDB_%selectAll  (                                                                                                                                          )
    call localGroupDB_%select     ('classification'           ,var_str('galaxy'    )                 ,comparisonEquals                         ,setOperatorIntersection      )
    call localGroupDB_%select     ('distanceMilkyWay'         ,         radiusOuter                  ,comparisonLessThan                       ,setOperatorIntersection      )
    call localGroupDB_%select     ('name'                     ,var_str('The Galaxy')                 ,comparisonEquals                         ,setOperatorRelativeComplement)
    call localGroupDB_%getProperty('massStellar'              ,massesTarget                          ,isPresentMasses                                                        )
    call localGroupDB_%getProperty('velocityDispersionStellar',velocityDispersionsTarget             ,isPresentVelocityDispersions                                           )
    call localGroupDB_%getProperty('velocityDispersionStellar',velocityDispersionsUncertaintiesTarget,isPresentVelocityDispersionsUncertainties,attributeUncertainty         )
    ! Convert velocity dispersions to logarithmic values.
    where (isPresentVelocityDispersions .and. isPresentVelocityDispersionsUncertainties)
       velocityDispersionsUncertaintiesTarget=velocityDispersionsUncertaintiesTarget/velocityDispersionsTarget/log(10.0d0)
    end where
    where (isPresentVelocityDispersions)
       velocityDispersionsTarget             =log10(velocityDispersionsTarget)
    end where
    ! Where velocityDispersions have no uncertainty, assign a default one.
    where (isPresentVelocityDispersions .and. .not.isPresentVelocityDispersionsUncertainties)
       velocityDispersionsUncertaintiesTarget=sizeUncertaintyDefault*log(10.0d0)
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
       if     (                                      &
            &   .not.isPresentMasses             (i) &
            &  .or.                                  &
            &   .not.isPresentVelocityDispersions(i) &
            & ) cycle
       j=int((log10(massesTarget(i))-log10(massMinimum))/(masses(2)-masses(1))+0.5d0,kind=c_size_t)+1_c_size_t
       if (j > 0 .and. j <= binCount) then
          functionValueTarget              (j)=functionValueTarget              (j)+velocityDispersionsTarget             (i)
          functionVarianceMeasurementTarget(j)=functionVarianceMeasurementTarget(j)+velocityDispersionsUncertaintiesTarget(i)**2
          functionVarianceSampleTarget     (j)=functionVarianceSampleTarget     (j)+velocityDispersionsTarget             (i)**2
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
    <referenceConstruct object="nodePropertyExtractor_"                                 constructor="nodePropertyExtractorMassStellar               (                                                                                                 )"/>
    !!]
    ! Create a velocity dispersion weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    radiusSpecifier(1)=var_str('stellarMassFraction{0.5}:all:galactic:lineOfSight:1.0')
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"                 constructor="nodePropertyExtractorVelocityDispersion        (radiusSpecifier,.false.,.false.,toleranceRelative,darkMatterHaloScale_                           )"/>
    !!]
    allocate(outputAnalysisWeightPropertyScalarizer_               )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyScalarizer_"                 constructor="nodePropertyExtractorScalarizer               (1,1,outputAnalysisWeightPropertyExtractor_                                                        )"/>
    !!]
    ! Build a size weight property operator.
    allocate(outputAnalysisWeightPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(velocityDispersionErrorPolynomialZeroPoint,velocityDispersionSystematicErrorPolynomialCoefficient )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                                                                  )"/>
    !!]
    allocate(weightPropertyOperators_                              )
    allocate(weightPropertyOperators_%next                         )
    weightPropertyOperators_     %operator_ => outputAnalysisWeightPropertyOperatorLog10_
    weightPropertyOperators_%next%operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence         (weightPropertyOperators_                                                                         )"/>
    !!]
    ! Create property operators and unoperators to perform conversion to/from logarithmic mass.
    allocate(outputAnalysisPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"                   constructor="outputAnalysisPropertyOperatorLog10            (                                                                                                 )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_"       constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorZeroPoint              ,systematicErrorPolynomialCoefficient                                )"/>
    !!]
    allocate(operators_     )
    allocate(operators_%next)
    operators_     %operator_ => outputAnalysisPropertyOperatorLog10_
    operators_%next%operator_ => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                        constructor="outputAnalysisPropertyOperatorSequence         (operators_                                                                                       )"/>
    !!]
    allocate(outputAnalysisPropertyUnoperator_               )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                      constructor="outputAnalysisPropertyOperatorAntiLog10        (                                                                                                                                           )"/>
    !!]
    ! Create a subsampling weight operator.
    allocate(outputAnalysisWeightOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                          constructor="outputAnalysisWeightOperatorSubsampling        (                                                                                                                                           )"/>
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
    allocate(galacticFilterHighPass_       )
    !![
    <referenceConstruct object="galacticFilterHighPass_"     >
     <constructor>
      galacticFilterHighPass     (                                                                &amp;
        &amp;                     threshold             =velocityDispersionMinimum              , &amp;
        &amp;                     nodePropertyExtractor_=outputAnalysisWeightPropertyScalarizer_  &amp;
        &amp;                    )
     </constructor>
    </referenceConstruct>
    !!]
    allocate(filters_               )
    allocate(filters_%next          )
    allocate(filters_%next%next     )
    allocate(filters_%next%next%next)
    filters_               %filter_ => galacticFilterHaloNotIsolated_ 
    filters_%next          %filter_ => galacticFilterHostMassRange_
    filters_%next%next     %filter_ => galacticFilterSurveyGeometry_
    filters_%next%next%next%filter_ => galacticFilterHighPass_
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
	   outputAnalysisMeanFunction1D(                                                                                                                  &amp;
	   &amp;                        var_str('localGroupMassVelocityDispersionRelation'                                                             ), &amp;
	   &amp;                        var_str('Mass-velocity dispersion relation of Local Group satellites'                                          ), &amp;
	   &amp;                        var_str('massStellar'                                                                                          ), &amp;
	   &amp;                        var_str('Stellar mass at the bin center'                                                                       ), &amp;
	   &amp;                        var_str('M☉'                                                                                                   ), &amp;
	   &amp;                        massSolar                                                                                                       , &amp;
	   &amp;                        var_str('velocityDispersionMean'                                                                               ), &amp;
           &amp;                        var_str('Mean logarithmic line-of-sight velocity dispersion at the half stellar-mass radius; ⟨log₁₀(σ/km s⁻¹)⟩'), &amp;
           &amp;                        var_str('dimensionless'                                                                                        ), &amp;
           &amp;                        0.0d0                                                                                                           , &amp;
	   &amp;                        massesNonZero                                                                                                   , &amp;
	   &amp;                        bufferCount                                                                                                     , &amp;
	   &amp;                        outputWeight                                                                                                    , &amp;
	   &amp;                        nodePropertyExtractor_                                                                                          , &amp;
           &amp;                        outputAnalysisWeightPropertyScalarizer_                                                                         , &amp;
	   &amp;                        outputAnalysisPropertyOperator_                                                                                 , &amp;
           &amp;                        outputAnalysisWeightPropertyOperator_                                                                           , &amp;
	   &amp;                        outputAnalysisPropertyUnoperator_                                                                               , &amp;
	   &amp;                        outputAnalysisWeightOperator_                                                                                   , &amp;
	   &amp;                        outputAnalysisDistributionOperator_                                                                             , &amp;
	   &amp;                        galacticFilter_                                                                                                 , &amp;
	   &amp;                        outputTimes_                                                                                                    , &amp;
	   &amp;                        outputAnalysisCovarianceModelBinomial                                                                           , &amp;
	   &amp;                        covarianceBinomialBinsPerDecade                                                                                 , &amp;
	   &amp;                        covarianceBinomialMassHaloMinimum                                                                               , &amp;
	   &amp;                        covarianceBinomialMassHaloMaximum                                                                               , &amp;
           &amp;                        likelihoodNormalize                                                                                             , &amp;
           &amp;                        var_str('$M_\star/\mathrm{M}_\odot$'                                                                           ), &amp;
           &amp;                        var_str('$\langle\log_{10}(\sigma_{\star, \mathrm{los}}/\mathrm{km}\,\mathrm{s}^{-1})\rangle$'                 ), &amp;
           &amp;                        .true.                                                                                                          , &amp;
           &amp;                        .false.                                                                                                         , &amp;
           &amp;                        var_str('Galacticus compilation'                                                                               ), &amp;
           &amp;                        functionValueTargetNonZero                                                                                      , &amp;
           &amp;                        functionCovarianceTargetNonZero                                                                                 , &amp;
	   &amp;                        massesWidthBin                                                                                                    &amp; 
	   &amp;                       )
	 </constructor>
       </referenceConstruct>
       !!]         
    end select
    !![
    <objectDestructor name="nodePropertyExtractor_"                                />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"                />
    <objectDestructor name="outputAnalysisWeightPropertyScalarizer_"               />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisWeightPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisPropertyOperator_"                       />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"                  />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"      />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"                     />
    <objectDestructor name="outputAnalysisWeightOperator_"                         />
    <objectDestructor name="outputAnalysisDistributionOperator_"                   />
    <objectDestructor name="galacticFilterHaloNotIsolated_"                        />
    <objectDestructor name="galacticFilterHostMassRange_"                          />
    <objectDestructor name="galacticFilterHighPass_"                               />
    <objectDestructor name="galacticFilterSurveyGeometry_"                         />
    <objectDestructor name="galacticFilter_"                                       />
    <objectDestructor name="surveyGeometry_"                                       />
    !!]
    nullify(filters_                )
    nullify(operators_              )
    nullify(weightPropertyOperators_)
    return
  end function localGroupMassVelocityDispersionRelationConstructorInternal

  subroutine localGroupMassVelocityDispersionRelationDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisLocalGroupMassVelocityDispersionRelation} output analysis class.
    !!}
    implicit none
    type(outputAnalysisLocalGroupMassVelocityDispersionRelation), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysis_"     />
    <objectDestructor name="self%outputTimes_"        />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine localGroupMassVelocityDispersionRelationDestructor

  subroutine localGroupMassVelocityDispersionRelationAnalyze(self,node,iOutput)
    !!{
    Implement a {\normalfont \ttfamily localGroupMassVelocityDispersionRelation} output analysis.
    !!}
    implicit none
    class  (outputAnalysisLocalGroupMassVelocityDispersionRelation), intent(inout) :: self
    type   (treeNode                                              ), intent(inout) :: node
    integer(c_size_t                                              ), intent(in   ) :: iOutput

    call self%outputAnalysis_%analyze(node,iOutput)
    return
  end subroutine localGroupMassVelocityDispersionRelationAnalyze

  subroutine localGroupMassVelocityDispersionRelationReduce(self,reduced)
    !!{
    Implement a {\normalfont \ttfamily localGroupMassVelocityDispersionRelation} output analysis reduction.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisLocalGroupMassVelocityDispersionRelation), intent(inout) :: self
    class(outputAnalysisClass                            ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisLocalGroupMassVelocityDispersionRelation)
       call self%outputAnalysis_%reduce(reduced%outputAnalysis_)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine localGroupMassVelocityDispersionRelationReduce

  subroutine localGroupMassVelocityDispersionRelationFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily localGroupMassVelocityDispersionRelation} output analysis finalization.
    !!}
    implicit none
    class(outputAnalysisLocalGroupMassVelocityDispersionRelation), intent(inout)           :: self
    type (varying_string                                        ), intent(in   ), optional :: groupName

    call self%outputAnalysis_%finalize(groupName)
    return
  end subroutine localGroupMassVelocityDispersionRelationFinalize

  double precision function localGroupMassVelocityDispersionRelationLogLikelihood(self)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily localGroupMassVelocityDispersionRelation} output analysis.
    !!}
    implicit none
    class(outputAnalysisLocalGroupMassVelocityDispersionRelation), intent(inout) :: self

    localGroupMassVelocityDispersionRelationLogLikelihood=self%outputAnalysis_%logLikelihood()
    return
  end function localGroupMassVelocityDispersionRelationLogLikelihood
