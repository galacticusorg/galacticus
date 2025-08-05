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
  Implements an output analysis class that computes mass functions for Local Group satellite galaxies.
  !!}

  !![
  <outputAnalysis name="outputAnalysisLocalGroupStellarMassFunction">
   <description>An output analysis class for Local Group satellite galaxy mass functions.</description>
   <deepCopy>
    <functionClass variables="volumeFunctionSatellites, volumeFunctionCentrals"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="volumeFunctionSatellites, volumeFunctionCentrals"/>
   </stateStorable>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisLocalGroupStellarMassFunction
     !!{
     An output analysis class for Local Group satellite galaxy mass functions.
     !!}
     private
     type            (outputAnalysisVolumeFunction1D), pointer                     :: volumeFunctionSatellites          => null(), volumeFunctionCentrals => null()
     class           (outputTimesClass              ), pointer                     :: outputTimes_                      => null()
     double precision                                , allocatable, dimension(:  ) :: randomErrorPolynomialCoefficient           , systematicErrorPolynomialCoefficient
     double precision                                , allocatable, dimension(:  ) :: masses                                     , massFunction                        , &
          &                                                                           massFunctionTarget
     double precision                                , allocatable, dimension(:,:) :: covariance
     logical                                                                       :: finalized
     integer                                                                       :: covarianceBinomialBinsPerDecade
     double precision                                                              :: covarianceBinomialMassHaloMinimum          , covarianceBinomialMassHaloMaximum   , &
          &                                                                           randomErrorMinimum                         , randomErrorMaximum                  , &
          &                                                                           negativeBinomialScatterFractional          , logLikelihoodZero                   , &
          &                                                                           countFailures
     type            (enumerationPositionTypeType   )                              :: positionType
   contains
     !![
     <methods>
       <method description="Finalize analysis." method="finalizeAnalysis" />
     </methods>
     !!]
     final     ::                     localGroupStellarMassFunctionDestructor
     procedure :: analyze          => localGroupStellarMassFunctionAnalyze
     procedure :: finalize         => localGroupStellarMassFunctionFinalize
     procedure :: finalizeAnalysis => localGroupStellarMassFunctionFinalizeAnalysis
     procedure :: reduce           => localGroupStellarMassFunctionReduce
     procedure :: logLikelihood    => localGroupStellarMassFunctionLogLikelihood
  end type outputAnalysisLocalGroupStellarMassFunction

  interface outputAnalysisLocalGroupStellarMassFunction
     !!{
     Constructors for the \refClass{outputAnalysisLocalGroupStellarMassFunction} output analysis class.
     !!}
     module procedure localGroupStellarMassFunctionConstructorParameters
     module procedure localGroupStellarMassFunctionConstructorInternal
  end interface outputAnalysisLocalGroupStellarMassFunction

contains

  function localGroupStellarMassFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisLocalGroupStellarMassFunction} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters            , only : inputParameter               , inputParameters
    use :: Output_Times                , only : outputTimes                  , outputTimesClass
    use :: Galactic_Filters            , only : enumerationPositionTypeEncode
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    type            (outputAnalysisLocalGroupStellarMassFunction)                              :: self
    type            (inputParameters                            ), intent(inout)               :: parameters
    class           (outputTimesClass                           ), pointer                     :: outputTimes_
    double precision                                             , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient , systematicErrorPolynomialCoefficient
    integer                                                                                    :: covarianceBinomialBinsPerDecade
    double precision                                                                           :: covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum   , &
         &                                                                                        randomErrorMinimum               , randomErrorMaximum                  , &
         &                                                                                        negativeBinomialScatterFractional, logLikelihoodZero
    type            (varying_string                             )                              :: positionType

    ! Check and read parameters.
    if (parameters%isPresent(    'randomErrorPolynomialCoefficient')) then
       allocate(    randomErrorPolynomialCoefficient(parameters%count(    'randomErrorPolynomialCoefficient')))
    else
       allocate(    randomErrorPolynomialCoefficient(1                                                   ))
    end if
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                   ))
    end if
    !![
    <inputParameter>
      <name>negativeBinomialScatterFractional</name>
      <source>parameters</source>
      <defaultValue>0.18d0</defaultValue>
      <defaultSource>\citep{boylan-kolchin_theres_2010}</defaultSource>
      <description>The fractional scatter (relative to the Poisson scatter) in the negative binomial distribution used in likelihood calculations.</description>
    </inputParameter>
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
    <inputParameter>
      <name>logLikelihoodZero</name>
      <source>parameters</source>
      <defaultValue>logImprobable</defaultValue>
      <description>The log-likelihood to assign to bins where the model expectation is zero.</description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=outputAnalysisLocalGroupStellarMassFunction(outputTimes_,enumerationPositionTypeEncode(positionType,includesPrefix=.false.),negativeBinomialScatterFractional,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,logLikelihoodZero)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function localGroupStellarMassFunctionConstructorParameters

  function localGroupStellarMassFunctionConstructorInternal(outputTimes_,positionType,negativeBinomialScatterFractional,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,logLikelihoodZero) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLocalGroupStellarMassFunction} output analysis class for internal use.
    !!}
    use :: Galactic_Filters                        , only : filterList                                         , galacticFilterAll                  , galacticFilterHaloIsolated            , galacticFilterHaloNotIsolated                  , &
          &                                                 galacticFilterHostMassRange                        , galacticFilterSurveyGeometry       , enumerationPositionTypeType
    use :: Geometry_Surveys                        , only : surveyGeometryCombined                             , surveyGeometryList                 , surveyGeometryLocalGroupClassical     , surveyGeometryLocalGroupDES                    , &
          &                                                 surveyGeometryLocalGroupSDSS
    use :: Interface_Local_Group_DB                , only : comparisonEquals                                   , comparisonLessThan                 , localGroupDB                          , setOperatorIntersection                        , &
          &                                                 setOperatorRelativeComplement                      , setOperatorUnion
    use :: Node_Property_Extractors                , only : nodePropertyExtractorMassStellar
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Numerical_Constants_Astronomical        , only : massSolar
    use :: Numerical_Ranges                        , only : Make_Range                                         , rangeTypeLinear
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerIdentity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10            , outputAnalysisPropertyOperatorLog10, outputAnalysisPropertyOperatorSequence, outputAnalysisPropertyOperatorSystmtcPolynomial, &
          &                                                 propertyOperatorList
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorSubsampling
    use :: Output_Times                            , only : outputTimesClass
    implicit none
    type            (outputAnalysisLocalGroupStellarMassFunction        )                                :: self
    integer                                                              , intent(in   )                 :: covarianceBinomialBinsPerDecade
    double precision                                                     , intent(in   )                 :: covarianceBinomialMassHaloMinimum                         , covarianceBinomialMassHaloMaximum              , &
         &                                                                                                  negativeBinomialScatterFractional                         , logLikelihoodZero
    double precision                                                     , intent(in   )                 :: randomErrorMinimum                                        , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:  ) :: randomErrorPolynomialCoefficient                          , systematicErrorPolynomialCoefficient
    type            (enumerationPositionTypeType                        ), intent(in   )                 :: positionType
    class           (outputTimesClass                                   ), intent(inout), target         :: outputTimes_
    type            (nodePropertyExtractorMassStellar                   )               , pointer        :: nodePropertyExtractor_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer        :: outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisPropertyOperatorLog10                )               , pointer        :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorSequence             )               , pointer        :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10            )               , pointer        :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorSubsampling            )               , pointer        :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerIdentity       )               , pointer        :: outputAnalysisDistributionNormalizerCentrals_             , outputAnalysisDistributionNormalizerSatellites_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer        :: outputAnalysisDistributionOperator_
    type            (surveyGeometryLocalGroupClassical                  )               , pointer        :: surveyGeometryClassical_
    type            (surveyGeometryLocalGroupSDSS                       )               , pointer        :: surveyGeometrySDSS_
    type            (surveyGeometryLocalGroupDES                        )               , pointer        :: surveyGeometryDES_
    type            (surveyGeometryCombined                             )               , pointer        :: surveyGeometry_
    type            (surveyGeometryList                                 )               , pointer        :: surveyGeometryList_
    type            (galacticFilterHaloIsolated                         )               , pointer        :: galacticFilterHaloIsolated_
    type            (galacticFilterHaloNotIsolated                      )               , pointer        :: galacticFilterHaloNotIsolated_
    type            (galacticFilterHostMassRange                        )               , pointer        :: galacticFilterHostMassRange_
    type            (galacticFilterSurveyGeometry                       )               , pointer        :: galacticFilterSurveyGeometry_
    type            (galacticFilterAll                                  )               , pointer        :: galacticFilterSatellites_                                 , galacticFilterCentrals_
    type            (filterList                                         )               , pointer        :: filtersSatellites_                                        , filtersCentrals_
    type            (propertyOperatorList                               )               , pointer        :: operators_
    double precision                                                     , allocatable  , dimension(:  ) :: massesSatellites                                          , massesCentrals                                             , &
         &                                                                                                  massesTarget
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeightSatellites                                    , outputWeightCentrals
    double precision                                                     , parameter                     :: bufferWidthLogarithmic                         =3.0d+0    , errorZeroPoint                                  =10.0d0
    integer         (c_size_t                                           ), parameter                     :: binCountSatellites                             =9_c_size_t, binCountCentrals                                =2_c_size_t, &
         &                                                                                                  bufferCountMinimum                             =5_c_size_t, bufferCountCentrals                             =0_c_size_t
    double precision                                                     , parameter                     :: massSatelliteMinimum                           =1.0d+2    , massSatelliteMaximum                            =1.0d10    , &
         &                                                                                                  massCentralMinimum                             =1.0d+0    , massCentralMaximum                              =1.0d20    , &
         &                                                                                                  radiusOuter                                    =3.0d-1    , massThresholdClassical                          =1.0d05
    integer         (c_size_t                                           )                                :: i                                                         , j                                                          , &
         &                                                                                                  bufferCountSatellites
    type            (localGroupDB                                       )                                :: localGroupDB_
    !![
    <constructorAssign variables="*outputTimes_, positionType, negativeBinomialScatterFractional, randomErrorMinimum, randomErrorMaximum, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum, logLikelihoodZero"/>
    !!]

    ! Initialize.
    self%finalized=.false.
    ! Compute failure count for the negative binomial distribution used in likelihood calculations.
    self%countFailures=1.0d0/negativeBinomialScatterFractional**2
    ! Construct mass bins.
    allocate(massesSatellites(binCountSatellites))
    allocate(massesCentrals  (binCountCentrals  ))
    massesSatellites=Make_Range(log10(massSatelliteMinimum),log10(massSatelliteMaximum),int(binCountSatellites),rangeTypeLinear)
    massesCentrals  =Make_Range(log10(massCentralMinimum  ),log10(massCentralMaximum  ),int(binCountCentrals  ),rangeTypeLinear)
    ! Construct the target distribution.
    !! Select Classical, SDSS, and DES Local Group galaxies, within 300 kpc of the Milky Way, and excluding the Milky Way itself,
    !! then retrieve stellar masses for the selected galaxies.
    localGroupDB_=localGroupDB()
    call localGroupDB_%select     ('discoverySurvey' ,var_str('classical' ),comparisonEquals  ,setOperatorUnion             )
    call localGroupDB_%select     ('discoverySurvey' ,var_str('SDSS'      ),comparisonEquals  ,setOperatorUnion             )
    call localGroupDB_%select     ('discoverySurvey' ,var_str('DES'       ),comparisonEquals  ,setOperatorUnion             )
    call localGroupDB_%select     ('classification'  ,var_str('galaxy'    ),comparisonEquals  ,setOperatorIntersection      )
    call localGroupDB_%select     ('distanceMilkyWay',         radiusOuter ,comparisonLessThan,setOperatorIntersection      )
    call localGroupDB_%select     ('name'            ,var_str('The Galaxy'),comparisonEquals  ,setOperatorRelativeComplement)
    call localGroupDB_%getProperty('massStellar'     ,massesTarget                                                          )
    allocate(self%massFunctionTarget(binCountSatellites))
    self%massFunctionTarget=0.0d0
    massesTarget      =log10(massesTarget)
    do i=1,size(massesTarget)
       j=int((massesTarget(i)-log10(massSatelliteMinimum))/(massesSatellites(2)-massesSatellites(1))+0.5d0,kind=c_size_t)+1_c_size_t
       if (j > 0 .and. j <= binCountSatellites) self%massFunctionTarget(j)=self%massFunctionTarget(j)+1.0d0
    end do
    ! Create a stellar mass property extractor.
    allocate(nodePropertyExtractor_                )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                           constructor="nodePropertyExtractorMassStellar               (                                                   )"/>
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
    ! Create an identity weight operator.
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
    ! Build survey geometry.
    !! For classical satellites use the corresponding geometry (which excludes the Zone of Avoidance) and assume they are detected out to the outer radius considered (300 kpc).
    allocate(surveyGeometryClassical_               )
    !![
    <referenceConstruct object="surveyGeometryClassical_"       constructor="surveyGeometryLocalGroupClassical(radiusOuter       ,massThresholdClassical)"/>
    !!]
    !! For SDSS satellites use the corresponding geometry and assume they are detected out to at most the outer radius considered (300 kpc).
    allocate(surveyGeometrySDSS_                    )
    !![
    <referenceConstruct object="surveyGeometrySDSS_"            constructor="surveyGeometryLocalGroupSDSS     (radiusOuter                              )"/>
    !!]
    allocate(surveyGeometryDES_                     )
    !![
    <referenceConstruct object="surveyGeometryDES_"             constructor="surveyGeometryLocalGroupDES      (radiusOuter                              )"/>
    !!]
    !! Combine the survey geometries.
    allocate(surveyGeometryList_            )
    allocate(surveyGeometryList_  %next     )
    allocate(surveyGeometryList_  %next%next)
    surveyGeometryList_          %surveyGeometry_ => surveyGeometryClassical_
    surveyGeometryList_%next     %surveyGeometry_ => surveyGeometrySDSS_
    surveyGeometryList_%next%next%surveyGeometry_ => surveyGeometryDES_
    allocate(surveyGeometry_                )
    !![
    <referenceConstruct object="surveyGeometry_"                constructor="surveyGeometryCombined           (surveyGeometryList_                       )"/>
    !!]
    ! Build filters which select satellites/centrals in a specified range of host halo mass, and applies a survey geometry.
    allocate(galacticFilterHaloIsolated_   )
    !![
    <referenceConstruct object="galacticFilterHaloIsolated_"    constructor="galacticFilterHaloIsolated       (                                          )"/>
    !!]
    allocate(galacticFilterHaloNotIsolated_)
    !![
    <referenceConstruct object="galacticFilterHaloNotIsolated_" constructor="galacticFilterHaloNotIsolated    (                                          )"/>
    !!]
    allocate(galacticFilterSurveyGeometry_ )
    !![
    <referenceConstruct object="galacticFilterSurveyGeometry_"  constructor="galacticFilterSurveyGeometry     (positionType,surveyGeometry_              )"/>
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
    allocate(filtersCentrals_            )
    allocate(filtersCentrals_  %next     )
    filtersCentrals_            %filter_ => galacticFilterHaloIsolated_
    filtersCentrals_       %next%filter_ => galacticFilterHostMassRange_
    allocate(galacticFilterCentrals_  )
    !![
    <referenceConstruct object="galacticFilterCentrals_"   constructor="galacticFilterAll(filtersCentrals_  )"/>
    !!]
    allocate(filtersSatellites_          )
    allocate(filtersSatellites_%next     )
    allocate(filtersSatellites_%next%next)
    filtersSatellites_          %filter_ => galacticFilterHaloNotIsolated_
    filtersSatellites_%next     %filter_ => galacticFilterHostMassRange_
    filtersSatellites_%next%next%filter_ => galacticFilterSurveyGeometry_
    allocate(galacticFilterSatellites_)
    !![
    <referenceConstruct object="galacticFilterSatellites_" constructor="galacticFilterAll(filtersSatellites_)"/>
    !!]
    ! Build an identity distribution normalizers for centrals and satellites.
    allocate(outputAnalysisDistributionNormalizerCentrals_  )
    allocate(outputAnalysisDistributionNormalizerSatellites_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerCentrals_"   constructor="outputAnalysisDistributionNormalizerIdentity()"/>
    <referenceConstruct object="outputAnalysisDistributionNormalizerSatellites_" constructor="outputAnalysisDistributionNormalizerIdentity()"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeightSatellites(binCountSatellites,outputTimes_%count()))
    allocate(outputWeightCentrals  (binCountCentrals  ,outputTimes_%count()))
    do i=1_c_size_t,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(i),0.0d0,absTol=1.0d-6)) then
          outputWeightSatellites(:,i)=1.0d0
          outputWeightCentrals  (:,i)=1.0d0
       else
          outputWeightSatellites(:,i)=0.0d0
          outputWeightCentrals  (:,i)=0.0d0
       end if
    end do
    if (any(sum(outputWeightSatellites,dim=2) /= 1.0d0)) call Error_Report('output weights do not equal unity for satellites'//{introspection:location})
    if (any(sum(outputWeightCentrals  ,dim=2) /= 1.0d0)) call Error_Report('output weights do not equal unity for centrals'  //{introspection:location})
    ! Compute the number of buffer bins to add to either side of the mass function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    bufferCountSatellites=max(int(bufferWidthLogarithmic/(massesSatellites(2)-massesSatellites(1)),kind=c_size_t)+1_c_size_t,bufferCountMinimum)
    ! Construct the volume function 1D objects.
    allocate(self%volumeFunctionSatellites)
    allocate(self%volumeFunctionCentrals  )
    !![
    <referenceConstruct isResult="yes" owner="self" object="volumeFunctionSatellites">
     <constructor>
      outputAnalysisVolumeFunction1D(                                                         &amp;
       &amp;                         var_str('localGroupStellarMassFunction')               , &amp;
       &amp;                         var_str('Mass function of Local Group satellites')     , &amp;
       &amp;                         var_str('massStellar')                                 , &amp;
       &amp;                         var_str('Stellar mass at the bin center')              , &amp;
       &amp;                         var_str('M☉')                                          , &amp;
       &amp;                         massSolar                                              , &amp;
       &amp;                         var_str('massFunction')                                , &amp;
       &amp;                         var_str('Differential satellite stellar mass function'), &amp;
       &amp;                         var_str(' ')                                           , &amp;
       &amp;                         0.0d0                                                  , &amp;
       &amp;                         massesSatellites                                       , &amp;
       &amp;                         bufferCountSatellites                                  , &amp;
       &amp;                         outputWeightSatellites                                 , &amp;
       &amp;                         nodePropertyExtractor_                                 , &amp;
       &amp;                         outputAnalysisPropertyOperator_                        , &amp;
       &amp;                         outputAnalysisPropertyUnoperator_                      , &amp;
       &amp;                         outputAnalysisWeightOperator_                          , &amp;
       &amp;                         outputAnalysisDistributionOperator_                    , &amp;
       &amp;                         outputAnalysisDistributionNormalizerSatellites_        , &amp;
       &amp;                         galacticFilterSatellites_                              , &amp;
       &amp;                         outputTimes_                                           , &amp;
       &amp;                         outputAnalysisCovarianceModelBinomial                  , &amp;
       &amp;                         covarianceBinomialBinsPerDecade                        , &amp;
       &amp;                         covarianceBinomialMassHaloMinimum                      , &amp;
       &amp;                         covarianceBinomialMassHaloMaximum                        &amp;
       &amp;                        )
     </constructor>
    </referenceConstruct>
    <referenceConstruct isResult="yes" owner="self" object="volumeFunctionCentrals">
     <constructor>
      outputAnalysisVolumeFunction1D(                                                        &amp;
       &amp;                         var_str(' ')                                          , &amp;
       &amp;                         var_str(' ')                                          , &amp;
       &amp;                         var_str(' ')                                          , &amp;
       &amp;                         var_str(' ')                                          , &amp;
       &amp;                         var_str(' ')                                          , &amp;
       &amp;                         0.0d0                                                 , &amp;
       &amp;                         var_str(' ')                                          , &amp;
       &amp;                         var_str(' ')                                          , &amp;
       &amp;                         var_str(' ')                                          , &amp;
       &amp;                         0.0d0                                                 , &amp;
       &amp;                         massesCentrals                                        , &amp;
       &amp;                         bufferCountCentrals                                   , &amp;
       &amp;                         outputWeightCentrals                                  , &amp;
       &amp;                         nodePropertyExtractor_                                , &amp;
       &amp;                         outputAnalysisPropertyOperator_                       , &amp;
       &amp;                         outputAnalysisPropertyUnoperator_                     , &amp;
       &amp;                         outputAnalysisWeightOperator_                         , &amp;
       &amp;                         outputAnalysisDistributionOperator_                   , &amp;
       &amp;                         outputAnalysisDistributionNormalizerCentrals_         , &amp;
       &amp;                         galacticFilterCentrals_                               , &amp;
       &amp;                         outputTimes_                                          , &amp;
       &amp;                         outputAnalysisCovarianceModelBinomial                 , &amp;
       &amp;                         covarianceBinomialBinsPerDecade                       , &amp;
       &amp;                         covarianceBinomialMassHaloMinimum                     , &amp;
       &amp;                         covarianceBinomialMassHaloMaximum                       &amp;
       &amp;                        )
     </constructor>
    </referenceConstruct>
    <objectDestructor name="nodePropertyExtractor_"                          />
    <objectDestructor name="outputAnalysisPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisPropertyUnoperator_"               />
    <objectDestructor name="outputAnalysisWeightOperator_"                   />
    <objectDestructor name="outputAnalysisDistributionOperator_"             />
    <objectDestructor name="galacticFilterHaloIsolated_"                     />
    <objectDestructor name="galacticFilterHaloNotIsolated_"                  />
    <objectDestructor name="galacticFilterHostMassRange_"                    />
    <objectDestructor name="galacticFilterCentrals_"                         />
    <objectDestructor name="galacticFilterSatellites_"                       />
    <objectDestructor name="galacticFilterSurveyGeometry_"                   />
    <objectDestructor name="outputAnalysisDistributionNormalizerSatellites_" />
    <objectDestructor name="outputAnalysisDistributionNormalizerCentrals_"   />
    <objectDestructor name="surveyGeometryClassical_"                        />
    <objectDestructor name="surveyGeometrySDSS_"                             />
    <objectDestructor name="surveyGeometryDES_"                              />
    <objectDestructor name="surveyGeometry_"                                 />
    !!]
    nullify(filtersSatellites_ )
    nullify(filtersCentrals_   )
    nullify(operators_         )
    nullify(surveyGeometryList_)
    return
  end function localGroupStellarMassFunctionConstructorInternal

  subroutine localGroupStellarMassFunctionDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisLocalGroupStellarMassFunction} output analysis class.
    !!}
    implicit none
    type(outputAnalysisLocalGroupStellarMassFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%volumeFunctionSatellites"/>
    <objectDestructor name="self%volumeFunctionCentrals"  />
    <objectDestructor name="self%outputTimes_"            />
    !!]
    return
  end subroutine localGroupStellarMassFunctionDestructor

  subroutine localGroupStellarMassFunctionAnalyze(self,node,iOutput)
    !!{
    Implement a {\normalfont \ttfamily localGroupStellarMassFunction} output analysis.
    !!}
    implicit none
    class  (outputAnalysisLocalGroupStellarMassFunction), intent(inout) :: self
    type   (treeNode                                   ), intent(inout) :: node
    integer(c_size_t                                   ), intent(in   ) :: iOutput

    ! Analyze for all three volume functions.
    call self%volumeFunctionSatellites%analyze(node,iOutput)
    call self%volumeFunctionCentrals  %analyze(node,iOutput)
    return
  end subroutine localGroupStellarMassFunctionAnalyze

  subroutine localGroupStellarMassFunctionReduce(self,reduced)
    !!{
    Implement a {\normalfont \ttfamily localGroupStellarMassFunction} output analysis reduction.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisLocalGroupStellarMassFunction), intent(inout) :: self
    class(outputAnalysisClass                        ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisLocalGroupStellarMassFunction)
       call self%volumeFunctionSatellites%reduce(reduced%volumeFunctionSatellites)
       call self%volumeFunctionCentrals  %reduce(reduced%volumeFunctionCentrals  )
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine localGroupStellarMassFunctionReduce

  subroutine localGroupStellarMassFunctionFinalizeAnalysis(self)
    !!{
    Finalize analysis of a {\normalfont \ttfamily localGroupStellarMassFunction} output analysis.
    !!}
    implicit none
    class           (outputAnalysisLocalGroupStellarMassFunction), intent(inout)                 :: self
    double precision                                             , allocatable  , dimension(:  ) :: massFunctionCentrals
    double precision                                             , allocatable  , dimension(:,:) :: covarianceCentrals
    double precision                                                                             :: weight              , weightVariance
    integer         (c_size_t                                   )                                :: i                   , j

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
    ! Retrieve results from our 1-D volume functions.
    call self%volumeFunctionSatellites%results(binCenter=self%masses,functionValue=self%massFunction        ,functionCovariance=self%covariance        )
    call self%volumeFunctionCentrals  %results(                      functionValue=     massFunctionCentrals,functionCovariance=     covarianceCentrals)
    ! Normalize the mass function.
    weight           = sum(massFunctionCentrals)
    weightVariance   = sum(covarianceCentrals  )
    if (weight > 0.0d0) then
       self%massFunction=+self%massFunction/weight
       self%covariance  =+self%covariance  /weight**2
       do    i=1_c_size_t,size(self%massFunction)
          do j=1_c_size_t,size(self%massFunction)
             self%covariance(i,j)=+self%covariance  (i,j)    &
                  &               +self%massFunction(i  )    &
                  &               *self%massFunction(  j)    &
                  &               *weightVariance            &
                  &               /weight                **2
          end do
       end do
    end if
    return
  end subroutine localGroupStellarMassFunctionFinalizeAnalysis

  subroutine localGroupStellarMassFunctionFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily localGroupStellarMassFunction} output analysis finalization.
    !!}
    use :: Output_HDF5                     , only : outputFile
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(outputAnalysisLocalGroupStellarMassFunction), intent(inout)           :: self
    type (varying_string                             ), intent(in   ), optional :: groupName
    type (hdf5Object                                 )               , target   :: analysesGroup, subGroup
    type (hdf5Object                                 )               , pointer  :: inGroup
    type (hdf5Object                                 )                          :: analysisGroup, dataset

    ! Finalize analysis.
    call self%finalizeAnalysis()
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'                         )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName)                    )
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('localGroupStellarMassFunction','Analysis of stellar mass functions of Local Group satellite galaxies')
    call analysisGroup%writeAttribute('Local Group stellar mass function','description'                                                                                  )
    call analysisGroup%writeAttribute('function1D'                       ,'type'                                                                                         )
    call analysisGroup%writeAttribute('$M_\star\,[\mathrm{M}_\odot]$'    ,'xAxisLabel'                                                                                   )
    call analysisGroup%writeAttribute('$N$'                              ,'yAxisLabel'                                                                                   )
    call analysisGroup%writeAttribute(.true.                             ,'xAxisIsLog'                                                                                   )
    call analysisGroup%writeAttribute(.true.                             ,'yAxisIsLog'                                                                                   )
    call analysisGroup%writeAttribute('massStellar'                      ,'xDataset'                                                                                     )
    call analysisGroup%writeAttribute('massFunction'                     ,'yDataset'                                                                                     )
    call analysisGroup%writeAttribute('massFunctionTarget'               ,'yDatasetTarget'                                                                               )
    call analysisGroup%writeAttribute('massFunctionCovariance'           ,'yCovariance'                                                                                  )
    call analysisGroup%writeAttribute('Observed'                         ,'targetLabel'                                                                                  )
    call analysisGroup%writeDataset  (self%masses                        ,'massStellar'           ,'Stellar mass at the bin center'              ,datasetReturned=dataset)
    call dataset      %writeAttribute('M☉'                               ,'units'                                                                                        )
    call dataset      %writeAttribute(massSolar                          ,'unitsInSI'                                                                                    )
    call analysisGroup%writeDataset  (self%massFunction                  ,'massFunction'          ,'Satellite number per bin [model]'                                    )
    call analysisGroup%writeDataset  (self%covariance                    ,'massFunctionCovariance','Satellite number per bin [model; covariance]'                        )
    call analysisGroup%writeDataset  (self%massFunctionTarget            ,'massFunctionTarget'    ,'Satellite number per bin [observed]'                                 )
    call analysisGroup%writeAttribute(self%logLikelihood     ()          ,'logLikelihood'                                                                                )
    !$ call hdf5Access%unset()
    return
  end subroutine localGroupStellarMassFunctionFinalize

  double precision function localGroupStellarMassFunctionLogLikelihood(self)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily localGroupStellarMassFunction} output analysis. The likelihood function
    assumes that the model prediction for the number of satellite galaxies in any given mass bin follows a negative binomial
    distribution as was found for dark matter subhalos \citep[][see also
    \protect\citealt{lu_connection_2016}]{boylan-kolchin_theres_2010}. This has been confirmed by examining the results of many
    tree realizations, although it in principal could be model-dependent.
    !!}
    use :: Numerical_Constants_Math         , only : Pi
    use :: Statistics_Distributions_Discrete, only : distributionFunctionDiscrete1DNegativeBinomial
    implicit none
    class           (outputAnalysisLocalGroupStellarMassFunction   ), intent(inout) :: self
    type            (distributionFunctionDiscrete1DNegativeBinomial)                :: distribution
    integer                                                                         :: i
    double precision                                                                :: negativeBinomialProbabilitySuccess, countEffective, &
         &                                                                             variance

    call self%finalizeAnalysis()
    localGroupStellarMassFunctionLogLikelihood=0.0d0
    do i=1,size(self%masses)
       if (self%massFunction(i) <= 0.0d0) then
          if (nint(self%massFunctionTarget(i)) > 0) localGroupStellarMassFunctionLogLikelihood=+localGroupStellarMassFunctionLogLikelihood &
               &                                                                        +self%logLikelihoodZero
       else
          negativeBinomialProbabilitySuccess =+  1.0d0                                                          &
               &                              /(                                                                &
               &                                +1.0d0                                                          &
               &                                +self%negativeBinomialScatterFractional**2*self%massFunction(i) &
               &                               )
          if (negativeBinomialProbabilitySuccess >= 1.0d0) then
             if (nint(self%massFunctionTarget(i)) > 0) localGroupStellarMassFunctionLogLikelihood=+localGroupStellarMassFunctionLogLikelihood &
                  &                                                                               +self%logLikelihoodZero
          else
             ! Compute the likelihood assuming a negative binomial distribution. Note that we "de-normalize" the likelihood by
             ! multiplying by √[2πσᵢ²] (the normalization term in the corresponding normal distribution). This is useful to allow
             ! (-logℒ) to be used as a metric for significant shifts in the model results, without changing the relative
             ! likelihood of models (as this de-normalization shift is a constant multiplicative factor).
             countEffective                            = dble(max(1.0d0,self%massFunctionTarget(i)))
             variance                                  =+       countEffective                       &
                  &                                     *(                                           &
                  &                                       +     1.0d0                                &
                  &                                       +self%negativeBinomialScatterFractional**2 &
                  &                                       *     countEffective                       &
                  &                                      )
             distribution                              = distributionFunctionDiscrete1DNegativeBinomial                (negativeBinomialProbabilitySuccess,     self%countFailures         )
             localGroupStellarMassFunctionLogLikelihood=+localGroupStellarMassFunctionLogLikelihood                                                                                          &
                  &                                     +distribution                                  %massLogarithmic(                                   nint(self%massFunctionTarget(i))) &
                  &                                     +0.50d0                                                                                                                              &
                  &                                     *log(                                                                                                                                &
                  &                                          +2.0d0                                                                                                                          &
                  &                                          *Pi                                                                                                                             &
                  &                                          *variance                                                                                                                       &
                  &                                         )
          end if
       end if
    end do
    return
  end function localGroupStellarMassFunctionLogLikelihood
