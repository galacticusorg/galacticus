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
  Implements an ICM X-ray luminosity-temperature relation output analysis class.
  !!}

  !![
  <outputAnalysis name="outputAnalysisICMXrayLuminosityTemperature">
   <description>An ICM X-ray luminosity-temperature relation output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisICMXrayLuminosityTemperature
     !!{
     An ICM X-ray luminosity-temperature relation output analysis class.
     !!}
     private
     double precision                          , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient          , randomErrorPolynomialCoefficient
     class           (darkMatterHaloScaleClass), pointer                     :: darkMatterHaloScale_                 => null()
     class           (coolingFunctionClass    ), pointer                     :: coolingFunction_                     => null()
     class           (cosmologyFunctionsClass ), pointer                     :: cosmologyFunctions_                  => null()
     double precision                                                        :: randomErrorMinimum                            , randomErrorMaximum
   contains
     final :: icmXrayLuminosityTemperatureDestructor
  end type outputAnalysisICMXrayLuminosityTemperature

  interface outputAnalysisICMXrayLuminosityTemperature
     !!{
     Constructors for the \refClass{outputAnalysisICMXrayLuminosityTemperature} output analysis class.
     !!}
     module procedure icmXrayLuminosityTemperatureConstructorParameters
     module procedure icmXrayLuminosityTemperatureConstructorInternal
  end interface outputAnalysisICMXrayLuminosityTemperature

contains

  function icmXrayLuminosityTemperatureConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisICMXrayLuminosityTemperature} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters       , only : inputParameter          , inputParameters
    use :: Cooling_Functions      , only : coolingFunctionClass
    use :: Cosmology_Functions    , only : cosmologyFunctionsClass
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
    implicit none
    type            (outputAnalysisICMXrayLuminosityTemperature)                              :: self
    type            (inputParameters                           ), intent(inout)               :: parameters
    double precision                                            , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient
    class           (outputTimesClass                          ), pointer                     :: outputTimes_
    class           (darkMatterHaloScaleClass                  ), pointer                     :: darkMatterHaloScale_
    class           (coolingFunctionClass                      ), pointer                     :: coolingFunction_
    class           (cosmologyFunctionsClass                   ), pointer                     :: cosmologyFunctions_
    double precision                                                                          :: randomErrorMinimum                  , randomErrorMaximum

    ! Check and read parameters.
    allocate(systematicErrorPolynomialCoefficient(max(1,parameters%count('systematicErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    allocate(    randomErrorPolynomialCoefficient(max(1,parameters%count(    'randomErrorPolynomialCoefficient',zeroIfNotPresent=.true.))))
    !![
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the random error polynomial.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.05d0</defaultValue>
      <description>The minimum random error for X-ray temperature.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.05d0</defaultValue>
      <description>The maximum random error for X-ray temperature.</description>
    </inputParameter>
    <objectBuilder class="outputTimes"         name="outputTimes_"         source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="coolingFunction"     name="coolingFunction_"     source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisICMXrayLuminosityTemperature(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,outputTimes_,cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"        />
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="coolingFunction_"    />
    !!]
    return
  end function icmXrayLuminosityTemperatureConstructorParameters

  function icmXrayLuminosityTemperatureConstructorInternal(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,outputTimes_,cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisICMXrayLuminosityTemperature} output analysis class for internal use.
    !!}
    use :: Galactic_Filters                      , only : filterList                                         , galacticFilterAll                      , galacticFilterBasicMass               , galacticFilterHaloIsolated
    use :: Error                                 , only : Error_Report
    use :: Node_Property_Extractors              , only : nodePropertyExtractorICMXRayLuminosity             , nodePropertyExtractorICMXRayTemperature
    use :: Numerical_Comparison                  , only : Values_Agree
    use :: Numerical_Constants_Prefixes          , only : kilo
    use :: Numerical_Constants_Units             , only : electronVolt                                       , ergs
    use :: Numerical_Ranges                      , only : Make_Range                                         , rangeTypeLinear
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10            , outputAnalysisPropertyOperatorLog10    , outputAnalysisPropertyOperatorSequence, outputAnalysisPropertyOperatorSystmtcPolynomial, &
          &                                               propertyOperatorList
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                          , only : outputTimesClass
    implicit none
    type            (outputAnalysisICMXrayLuminosityTemperature         )                                :: self
    double precision                                                     , intent(in   )                 :: randomErrorMinimum                                      , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                    , randomErrorPolynomialCoefficient
    class           (outputTimesClass                                   ), intent(inout), target         :: outputTimes_
    class           (cosmologyFunctionsClass                            ), intent(in   ), target         :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                           ), intent(in   ), target         :: darkMatterHaloScale_
    class           (coolingFunctionClass                               ), intent(in   ), target         :: coolingFunction_
    integer                                                              , parameter                     :: covarianceBinomialBinsPerDecade                 =10
    double precision                                                     , parameter                     :: covarianceBinomialMassHaloMinimum               = 1.0d08, covarianceBinomialMassHaloMaximum    =1.0d16
    double precision                                                     , allocatable  , dimension(:  ) :: temperatures
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeight
    type            (galacticFilterAll                                  ), pointer                       :: galacticFilter_
    type            (galacticFilterHaloIsolated                         ), pointer                       :: galacticFilterHaloIsolated_
    type            (galacticFilterBasicMass                            ), pointer                       :: galacticFilterBasicMass_
    type            (filterList                                         ), pointer                       :: filters_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity               ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorSequence             ), pointer                       :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: outputAnalysisPropertyOperatorLog10_                    , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10            ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (nodePropertyExtractorICMXRayTemperature            ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorICMXRayLuminosity             ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    ), pointer                       :: outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (propertyOperatorList                               ), pointer                       :: propertyOperators_
    double precision                                                     , parameter                     :: errorPolynomialZeroPoint                        =1.0d0
    integer         (c_size_t                                           ), parameter                     :: bufferCount                                     =10
    logical                                                              , parameter                     :: likelihoodNormalize                             =.false.
    double precision                                                     , parameter                     :: temperatureMinimum                              =0.1d0  , temperatureMaximum                   =1.0d01, &
         &                                                                                                  countTemperaturesPerDecade                      =5.0d0
    integer         (c_size_t                                           )                                :: iOutput                                                 , countTemperatures
    !![
    <constructorAssign variables="systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient, randomErrorMinimum, randomErrorMaximum, *cosmologyFunctions_, *darkMatterHaloScale_, *coolingFunction_"/>
    !!]
    
    ! Construct bins in temperature.
    countTemperatures=int(log10(temperatureMaximum/temperatureMinimum)*countTemperaturesPerDecade+1.0d0,kind=c_size_t)
    allocate(temperatures(countTemperatures))
    temperatures=Make_Range(log10(temperatureMinimum),log10(temperatureMaximum),int(countTemperatures),rangeTypeLinear)    
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(countTemperatures,outputTimes_%count()))
    outputWeight=0.0d0
    do iOutput=1,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(iOutput),0.0d0,absTol=1.0d-10)) outputWeight(:,iOutput)=1.0d0
    end do
    if (any(sum(outputWeight,dim=2) /= 1.0d0)) call Error_Report('zero redshift output is required'//{introspection:location})
    ! Build a filter which selects isolated nodes above some coarse lower mass limit suitable for this analysis.
    allocate(galacticFilterHaloIsolated_)
    !![
    <referenceConstruct object="galacticFilterHaloIsolated_"                      constructor="galacticFilterHaloIsolated                     (                                                                                                             )"/>
    !!]
    allocate(galacticFilterBasicMass_)
    !![
    <referenceConstruct object="galacticFilterBasicMass_"                         constructor="galacticFilterBasicMass                        (massThreshold=1.0d11                                                                                         )"/>
    !!]
    allocate(galacticFilter_                                 )
    allocate(filters_                                        )
    allocate(filters_%next                                   )
    filters_     %filter_ => galacticFilterHaloIsolated_
    filters_%next%filter_ => galacticFilterBasicMass_
    !![
    <referenceConstruct object="galacticFilter_"                                  constructor="galacticFilterAll                              (filters_                                                                                                     )"/>
    !!]
    ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                    constructor="outputAnalysisWeightOperatorIdentity           (                                                                                                             )"/>
    !!]
    ! Build systematic, and log10() property operators.
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient                                                )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                                                                             )"/>
    !!]
    allocate(propertyOperators_                              )
    allocate(propertyOperators_%next                         )
    propertyOperators_     %operator_ => outputAnalysisPropertyOperatorLog10_
    propertyOperators_%next%operator_ => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence         (propertyOperators_                                                                                           )"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_             )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
      outputAnalysisDistributionOperatorRandomErrorPlynml(                                  &amp;
        &amp;                                             randomErrorMinimum              , &amp;
        &amp;                                             randomErrorMaximum              , &amp;
        &amp;                                             errorPolynomialZeroPoint        , &amp;
        &amp;                                             randomErrorPolynomialCoefficient  &amp;
        &amp;                                            )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build weight property operators.
    allocate(outputAnalysisWeightPropertyOperator_           )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"            constructor="outputAnalysisPropertyOperatorLog10            (                                                                                                             )"/>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_               )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                constructor="outputAnalysisPropertyOperatorAntiLog10        (                                                                                                             )"/>
    !!]
    ! Create an X-ray temperature property extractor.
    allocate(nodePropertyExtractor_                          )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                           constructor="nodePropertyExtractorICMXRayTemperature        (cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_                                                    )"/>
    !!]
    ! Create an X-ray luminosity property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_          )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"           constructor="nodePropertyExtractorICMXRayLuminosity         (cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_                                                    )"/>
    !!]
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                                              &
         &                                                         var_str('icmXrayLuminosityTemperature'                                     ), &
         &                                                         var_str('ICM X-ray luminosity-temperature relation'                        ), &
         &                                                         var_str('temperatureICMXray'                                               ), &
         &                                                         var_str('X-ray temperature of the ICM'                                     ), &
         &                                                         var_str('keV'                                                              ), &
         &                                                         kilo*electronVolt                                                           , &
         &                                                         var_str('luminosityICMXray'                                                ), &
         &                                                         var_str('Mean logarithmic (base-10) X-ray luminosity of the ICM'           ), &
         &                                                         var_str('ergs/s'                                                           ), &
         &                                                         ergs                                                                        , &
         &                                                         temperatures                                                                , &
         &                                                         bufferCount                                                                 , &
         &                                                         outputWeight                                                                , &
         &                                                         nodePropertyExtractor_                                                      , &
         &                                                         outputAnalysisWeightPropertyExtractor_                                      , &
         &                                                         outputAnalysisPropertyOperator_                                             , &
         &                                                         outputAnalysisWeightPropertyOperator_                                       , &
         &                                                         outputAnalysisPropertyUnoperator_                                           , &
         &                                                         outputAnalysisWeightOperator_                                               , &
         &                                                         outputAnalysisDistributionOperator_                                         , &
         &                                                         galacticFilter_                                                             , &
         &                                                         outputTimes_                                                                , &
         &                                                         outputAnalysisCovarianceModelBinomial                                       , &
         &                                                         covarianceBinomialBinsPerDecade                                             , &
         &                                                         covarianceBinomialMassHaloMinimum                                           , &
         &                                                         covarianceBinomialMassHaloMaximum                                           , &
         &                                                         likelihoodNormalize                                                         , &
         &                                                         var_str('$T_\mathrm{ICM}$ [keV]'                                           ), &
         &                                                         var_str('$\langle \log_{10} L_\mathrm{ICM}/\mathrm{ergs\, s}^{-1} \rangle$'), &
         &                                                         .true.                                                                      , &
         &                                                         .false.                                                                       &
         &                                                        )
    ! Clean up.
    !![
    <objectDestructor name="galacticFilter_"                                 />
    <objectDestructor name="galacticFilterHaloIsolated_"                     />
    <objectDestructor name="galacticFilterBasicMass_"                        />
    <objectDestructor name="outputAnalysisDistributionOperator_"             />
    <objectDestructor name="outputAnalysisWeightOperator_"                   />
    <objectDestructor name="outputAnalysisPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"           />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisPropertyUnoperator_"               />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"           />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"          />
    <objectDestructor name="nodePropertyExtractor_"                          />
    !!]
    nullify(propertyOperators_)
    return
  end function icmXrayLuminosityTemperatureConstructorInternal


  subroutine icmXrayLuminosityTemperatureDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisICMXrayLuminosityTemperature} output analysis class.
    !!}
    implicit none
    type(outputAnalysisICMXrayLuminosityTemperature), intent(inout) :: self
    
    !![
    <objectDestructor name="self%outputTimes_"        />
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%coolingFunction_"    />
    !!]
    return
  end subroutine icmXrayLuminosityTemperatureDestructor
