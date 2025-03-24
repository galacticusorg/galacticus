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
Implements a correlation function output analysis class for the \cite{hearin_dark_2013} analysis.
!!}

  !![
  <outputAnalysis name="outputAnalysisCorrelationFunctionHearin2013SDSS">
   <description>A correlation function output analysis class for the \cite{hearin_dark_2013} analysis.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisCorrelationFunction) :: outputAnalysisCorrelationFunctionHearin2013SDSS
     !!{
     A correlation function function output analysis class for the \cite{hearin_dark_2013} analysis.
     !!}
     private
     double precision                        , allocatable, dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient
     double precision                                                    :: randomErrorMinimum                        , randomErrorMaximum
  end type outputAnalysisCorrelationFunctionHearin2013SDSS

  interface outputAnalysisCorrelationFunctionHearin2013SDSS
     !!{
     Constructors for the {\normalfont \ttfamily correlationFunctionHearin2013SDSS} output analysis class.
     !!}
     module procedure correlationFunctionHearin2013SDSSConstructorParameters
     module procedure correlationFunctionHearin2013SDSSConstructorInternal
  end interface outputAnalysisCorrelationFunctionHearin2013SDSS

contains

  function correlationFunctionHearin2013SDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily correlationFunctionHearin2013SDSS} output analysis class which takes a parameter set as input.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: Input_Parameters  , only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisCorrelationFunctionHearin2013SDSS)                              :: self
    type            (inputParameters                                ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass                        ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                               ), pointer                     :: outputTimes_
    class           (darkMatterProfileDMOClass                      ), pointer                     :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass                        ), pointer                     :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass                       ), pointer                     :: darkMatterHaloScale_
    class           (haloModelPowerSpectrumModifierClass            ), pointer                     :: haloModelPowerSpectrumModifier_
    class           (powerSpectrumClass                             ), pointer                     :: powerSpectrum_
    double precision                                                 , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient
    double precision                                                                               :: massHaloMinimum                 , massHaloMaximum                     , &
         &                                                                                            randomErrorMinimum              , randomErrorMaximum
    integer                                                                                        :: massHaloBinsPerDecade

    ! Check and read parameters.
    if (parameters%isPresent(    'randomErrorPolynomialCoefficient')) then
       allocate(    randomErrorPolynomialCoefficient(parameters%count(    'randomErrorPolynomialCoefficient')))
    else
       allocate(    randomErrorPolynomialCoefficient(1                                                       ))
    end if
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                       ))
    end if
    !![
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.07d0]</defaultValue>
      <description>The coefficients of the random error polynomial for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>massHaloBinsPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing the mass function covariance matrix for main branch galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massHaloMinimum</name>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massHaloMaximum</name>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"/>
    <objectBuilder class="outputTimes"                    name="outputTimes_"                    source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"           name="darkMatterProfileDMO_"           source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"             name="darkMatterHaloBias_"             source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"/>
    <objectBuilder class="powerSpectrum"                  name="powerSpectrum_"                  source="parameters"/>
    <objectBuilder class="haloModelPowerSpectrumModifier" name="haloModelPowerSpectrumModifier_" source="parameters"/>
    !!]
    self=outputAnalysisCorrelationFunctionHearin2013SDSS(massHaloBinsPerDecade,massHaloMinimum, massHaloMaximum,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,cosmologyFunctions_,outputTimes_,darkMatterProfileDMO_,darkMatterHaloBias_,darkMatterHaloScale_,haloModelPowerSpectrumModifier_,powerSpectrum_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="outputTimes_"                   />
    <objectDestructor name="darkMatterProfileDMO_"          />
    <objectDestructor name="darkMatterHaloBias_"            />
    <objectDestructor name="darkMatterHaloScale_"           />
    <objectDestructor name="haloModelPowerSpectrumModifier_"/>
    <objectDestructor name="powerSpectrum_"                 />
    !!]
    return
  end function correlationFunctionHearin2013SDSSConstructorParameters

  function correlationFunctionHearin2013SDSSConstructorInternal(massHaloBinsPerDecade,massHaloMinimum,massHaloMaximum,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,cosmologyFunctions_,outputTimes_,darkMatterProfileDMO_,darkMatterHaloBias_,darkMatterHaloScale_,haloModelPowerSpectrumModifier_,powerSpectrum_) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily correlationFunctionHearin2013SDSS} output analysis class for internal use.
    !!}
    use            :: Cosmology_Functions                   , only : cosmologyFunctionsClass                            , cosmologyFunctionsMatterLambda
    use            :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use            :: Galactic_Filters                      , only : galacticFilterStellarMass
    use            :: Input_Paths                           , only : inputPath                                          , pathTypeDataStatic
    use            :: Geometry_Surveys                      , only : surveyGeometryHearin2014SDSS
    use, intrinsic :: ISO_C_Binding                         , only : c_size_t
    use            :: Node_Property_Extractors              , only : nodePropertyExtractorMassStellar
    use            :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use            :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorCsmlgyAnglrDstnc     , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorLog10, outputAnalysisPropertyOperatorSequence, &
          &                                                          outputAnalysisPropertyOperatorSystmtcPolynomial    , propertyOperatorList
    implicit none
    type            (outputAnalysisCorrelationFunctionHearin2013SDSS    )                                :: self
    double precision                                                     , intent(in   )  , dimension(:) :: randomErrorPolynomialCoefficient              , systematicErrorPolynomialCoefficient
    double precision                                                     , intent(in   )                 :: massHaloMinimum                               , massHaloMaximum                     , &
         &                                                                                                  randomErrorMinimum                            , randomErrorMaximum
    integer                                                              , intent(in   )                 :: massHaloBinsPerDecade
    class           (cosmologyFunctionsClass                            ), intent(in   ), target         :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(in   ), target         :: outputTimes_
    class           (darkMatterProfileDMOClass                          ), intent(in   ), target         :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass                            ), intent(in   ), target         :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass                           ), intent(in   ), target         :: darkMatterHaloScale_
    class           (haloModelPowerSpectrumModifierClass                ), intent(in   ), target         :: haloModelPowerSpectrumModifier_
    class           (powerSpectrumClass                                 ), intent(in   ), target         :: powerSpectrum_
    type            (cosmologyParametersSimple                          ), pointer                       :: cosmologyParametersData_
    type            (cosmologyFunctionsMatterLambda                     ), pointer                       :: cosmologyFunctionsData_
    type            (galacticFilterStellarMass                          ), pointer                       :: galacticFilter_
    type            (surveyGeometryHearin2014SDSS                       ), pointer                       :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: massPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer        :: massPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer        :: massDistributionOperator_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc    ), pointer                       :: massPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc     ), pointer                       :: separationPropertyOperator_
    type            (nodePropertyExtractorMassStellar                   ), pointer                       :: massPropertyExtractor_
    type            (outputAnalysisPropertyOperatorSequence             ), pointer                       :: massPropertyOperator_
    type            (propertyOperatorList                               ), pointer                       :: propertyOperators_
    double precision                                                     , parameter                     :: errorPolynomialZeroPoint              =11.3d+0
    integer         (c_size_t                                           ), parameter                     :: wavenumberCount                       =60
    double precision                                                     , parameter                     :: wavenumberMinimum                     = 1.0d-3, wavenumberMaximum=1.0d+4
    logical                                                              , parameter                     :: halfIntegral                          =.false.
    !![
    <constructorAssign variables="randomErrorMinimum, randomErrorMaximum, randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient"/>
    !!]
    
    ! Build a filter which selects galaxies above some minimum stellar mass.
    allocate(galacticFilter_         )
    !![
    <referenceConstruct object="galacticFilter_" constructor="galacticFilterStellarMass   (                    1.0d8              )"/>
    !!]
    ! Build the SDSS survey geometry of Hearin et al. (2013) with their imposed redshift limits.
    allocate(surveyGeometry_         )
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryHearin2014SDSS(cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Create the data cosmology.
    allocate(cosmologyParametersData_)
    allocate(cosmologyFunctionsData_ )
    !![
    <referenceConstruct object="cosmologyParametersData_">
     <constructor>
       cosmologyParametersSimple     (                         &amp;
        &amp;                         OmegaMatter    = 0.27d0, &amp;
        &amp;                         OmegaDarkEnergy= 0.73d0, &amp;
        &amp;                         HubbleConstant =70.00d0, &amp;
        &amp;                         temperatureCMB = 0.00d0, &amp;
        &amp;                         OmegaBaryon    = 0.00d0  &amp;
        &amp;                        )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData_">
     <constructor>
       cosmologyFunctionsMatterLambda(                         &amp;
        &amp;                         cosmologyParametersData_ &amp;
        &amp;                        )
     </constructor>
    </referenceConstruct>
    !!]
    ! Stellar mass property extractor.
    allocate(massPropertyExtractor_                )
    !![
    <referenceConstruct object="massPropertyExtractor_" constructor="nodePropertyExtractorMassStellar                               (                                                                           )"/>
    !!]
    ! Sequence of property operators to correct for cosmological model, convert to logarithm, and apply systematic errors.
    allocate(massPropertyOperatorCsmlgyLmnstyDstnc_)
    !![
    <referenceConstruct object="massPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc(cosmologyFunctions_     ,cosmologyFunctionsData_              ,outputTimes_)"/>
    !!]
    allocate(massPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="massPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                                           )"/>
    !!]
    ! Systematic error model.
    allocate(massPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="massPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient              )"/>
    !!]
    allocate(propertyOperators_                    )
    allocate(propertyOperators_%next               )
    allocate(propertyOperators_%next%next          )
    propertyOperators_          %operator_ => massPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_%next     %operator_ => massPropertyOperatorLog10_
    propertyOperators_%next%next%operator_ => massPropertyOperatorSystmtcPolynomial_
    allocate(massPropertyOperator_                 )
    !![
    <referenceConstruct object="massPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence         (propertyOperators_                                                         )"/>
    !!]
    ! Build a random error distribution operator.
    allocate(massDistributionOperator_)
    !![
    <referenceConstruct object="massDistributionOperator_">
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
    ! Build an operator for separations which corrects for cosmological model.
    allocate(separationPropertyOperator_)
    !![
    <referenceConstruct object="separationPropertyOperator_"            constructor="outputAnalysisPropertyOperatorCsmlgyAnglrDstnc (cosmologyFunctions_     ,cosmologyFunctionsData_              ,outputTimes_)"/>
    !!]
    ! Build the object.
    self%outputAnalysisCorrelationFunction=                                                                                                                              &
         & outputAnalysisCorrelationFunction(                                                                                                                            &
         &                                   var_str('Hearin2013SDSS'                                                 )                                                , &
         &                                   var_str('Correlation function for the Hearin et al. (2013) SDSS analysis')                                                , &
         &                                   char(inputPath(pathTypeDataStatic)//'/observations/correlationFunctions/Projected_Correlation_Functions_Hearin_2013.hdf5'), &
         &                                   massHaloBinsPerDecade                                                                                                     , &
         &                                   massHaloMinimum                                                                                                           , &
         &                                   massHaloMaximum                                                                                                           , &
         &                                   wavenumberCount                                                                                                           , &
         &                                   wavenumberMinimum                                                                                                         , &
         &                                   wavenumberMaximum                                                                                                         , &
         &                                   halfIntegral                                                                                                              , &
         &                                   galacticFilter_                                                                                                           , &
         &                                   surveyGeometry_                                                                                                           , &
         &                                   cosmologyFunctions_                                                                                                       , &
         &                                   outputTimes_                                                                                                              , &
         &                                   darkMatterProfileDMO_                                                                                                     , &
         &                                   darkMatterHaloBias_                                                                                                       , &
         &                                   darkMatterHaloScale_                                                                                                      , &
         &                                   haloModelPowerSpectrumModifier_                                                                                           , &
         &                                   powerSpectrum_                                                                                                            , &
         &                                   massDistributionOperator_                                                                                                 , &
         &                                   massPropertyOperator_                                                                                                     , &
         &                                   separationPropertyOperator_                                                                                               , &
         &                                   massPropertyExtractor_                                                                                                      &
         &                                  )
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                       />
    <objectDestructor name="galacticFilter_"                       />
    <objectDestructor name="cosmologyParametersData_"              />
    <objectDestructor name="cosmologyFunctionsData_"               />
    <objectDestructor name="massPropertyExtractor_"                />
    <objectDestructor name="massPropertyOperatorCsmlgyLmnstyDstnc_"/>
    <objectDestructor name="massPropertyOperatorLog10_"            />
    <objectDestructor name="massPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="massPropertyOperator_"                 />
    <objectDestructor name="massDistributionOperator_"             />
    <objectDestructor name="separationPropertyOperator_"           />
    !!]
    nullify(propertyOperators_)
    return
  end function correlationFunctionHearin2013SDSSConstructorInternal
