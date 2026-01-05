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
  Implements a thermal Sunyaev-Zeldovich signal vs. stellar mass analysis class.
  !!}

  !![
  <outputAnalysis name="outputAnalysisSunyaevZeldovichPlanck2013">
   <description>A thermal Sunyaev-Zeldovich signal vs. stellar mass analysis class using the results of \cite{planck_collaboration_planck_2013}.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisSunyaevZeldovichPlanck2013
     !!{
     A thermal Sunyaev-Zeldovich signal vs stellar mass analysis class using the results of \cite{planck_collaboration_planck_2013}.
     !!}
     private
     double precision                          , allocatable, dimension(:) :: systematicErrorPolynomialCoefficient          , randomErrorPolynomialCoefficient
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_                 => null()
     class           (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_                  => null()
     class           (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_                 => null()
     class           (chemicalStateClass      ), pointer                   :: chemicalState_                       => null()
     double precision                                                      :: randomErrorMinimum                            , randomErrorMaximum
   contains
     final :: sunyaevZeldovichPlanck2013Destructor
  end type outputAnalysisSunyaevZeldovichPlanck2013
  
  interface outputAnalysisSunyaevZeldovichPlanck2013
     !!{
     Constructors for the \refClass{outputAnalysisSunyaevZeldovichPlanck2013} output analysis class.
     !!}
     module procedure sunyaevZeldovichPlanck2013ConstructorParameters
     module procedure sunyaevZeldovichPlanck2013ConstructorInternal
  end interface outputAnalysisSunyaevZeldovichPlanck2013

contains

  function sunyaevZeldovichPlanck2013ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisSunyaevZeldovichPlanck2013} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisSunyaevZeldovichPlanck2013)                              :: self
    type            (inputParameters                         ), intent(inout)               :: parameters
    double precision                                          , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient
    class           (cosmologyParametersClass                ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass                 ), pointer                     :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                ), pointer                     :: darkMatterHaloScale_
    class           (chemicalStateClass                      ), pointer                     :: chemicalState_
    class           (outputTimesClass                        ), pointer                     :: outputTimes_
    double precision                                                                        :: randomErrorMinimum                  , randomErrorMaximum

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
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for stellar masses.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"       name="cosmologyParameters_"       source="parameters"/>
    <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"       name="darkMatterHaloScale_"       source="parameters"/>
    <objectBuilder class="chemicalState"             name="chemicalState_"             source="parameters"/>
    <objectBuilder class="outputTimes"               name="outputTimes_"               source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisSunyaevZeldovichPlanck2013(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,chemicalState_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="chemicalState_"      />
    <objectDestructor name="outputTimes_"        />
    !!]
    return
  end function sunyaevZeldovichPlanck2013ConstructorParameters

  function sunyaevZeldovichPlanck2013ConstructorInternal(systematicErrorPolynomialCoefficient,randomErrorPolynomialCoefficient,randomErrorMinimum,randomErrorMaximum,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,chemicalState_,outputTimes_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisSunyaevZeldovichPlanck2013} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                            , cosmologyFunctionsMatterLambda                 , densityCosmologicalCritical
    use :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterStellarMass                          , galacticFilterAll                              , galacticFilterHaloIsolated               , &
         &                                                filterList
    use :: Geometry_Surveys                      , only : surveyGeometryLiWhite2009SDSS
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassStellar                   , nodePropertyExtractorICMSZ
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorRandomErrorPlynml
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10            , outputAnalysisPropertyOperatorLog10            , outputAnalysisPropertyOperatorCosmologySZ, &
          &                                               outputAnalysisPropertyOperatorSequence             , outputAnalysisPropertyOperatorSystmtcPolynomial, propertyOperatorList
    use :: Output_Analysis_Utilities             , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                          , only : outputTimesClass
    use :: Statistics_Distributions              , only : distributionFunction1DBeta
    implicit none
    type            (outputAnalysisSunyaevZeldovichPlanck2013           )                                :: self
    double precision                                                     , intent(in   )                 :: randomErrorMinimum                                        , randomErrorMaximum
    double precision                                                     , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient                      , randomErrorPolynomialCoefficient
    class           (cosmologyParametersClass                           ), intent(inout), target         :: cosmologyParameters_
    class           (cosmologyFunctionsClass                            ), intent(inout), target         :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                           ), intent(in   ), target         :: darkMatterHaloScale_
    class           (chemicalStateClass                                 ), intent(in   ), target         :: chemicalState_
    class           (outputTimesClass                                   ), intent(inout), target         :: outputTimes_
    integer                                                              , parameter                     :: covarianceBinomialBinsPerDecade                 =10
    double precision                                                     , parameter                     :: covarianceBinomialMassHaloMinimum               = 1.000d08, covarianceBinomialMassHaloMaximum=1.0d16
    double precision                                                     , allocatable  , dimension(:  ) :: masses                                                    , functionValueTarget                     , &
         &                                                                                                  functionErrorTarget
    double precision                                                     , allocatable  , dimension(:,:) :: outputWeight                                              , functionCovarianceTarget
    type            (galacticFilterStellarMass                          ), pointer                       :: galacticFilterStellarMass_
    type            (galacticFilterHaloIsolated                         ), pointer                       :: galacticFilterCentrals_
    type            (galacticFilterAll                                  ), pointer                       :: galacticFilter_
    type            (filterList                                         ), pointer                       :: filters_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity               ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorSequence             ), pointer                       :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorLog10                ), pointer                       :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10            ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorCosmologySZ          ), pointer                       :: outputAnalysisWeightPropertyOperator_
    type            (nodePropertyExtractorMassStellar                   ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorICMSZ                         ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    ), pointer                       :: outputAnalysisPropertyOperatorSystmtcPolynomial_
    type            (propertyOperatorList                               ), pointer                       :: propertyOperators_
    type            (surveyGeometryLiWhite2009SDSS                      ), pointer                       :: surveyGeometry_
    type            (cosmologyParametersSimple                          ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     ), pointer                       :: cosmologyFunctionsData
    logical                                                              , parameter                     :: likelihoodNormalize                             =.false.
    double precision                                                     , parameter                     :: errorPolynomialZeroPoint                        =11.300d00
    integer         (c_size_t                                           ), parameter                     :: bufferCount                                     =10
    integer         (c_size_t                                           )                                :: iBin                                                      , binCount
    !![
    <constructorAssign variables="systematicErrorPolynomialCoefficient, randomErrorPolynomialCoefficient, randomErrorMinimum, randomErrorMaximum, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_, *chemicalState_"/>
    !!]
    
    ! Construct the target data.
    binCount=20_c_size_t
    allocate(masses(binCount))
    allocate(functionValueTarget     (binCount         ))
    allocate(functionErrorTarget     (binCount         ))
    allocate(functionCovarianceTarget(binCount,binCount))
    ! Target data: Planck Intermediate Results XI (https://ui.adsabs.harvard.edu/abs/2013A%26A...557A..52P), Table 1. Errors are
    ! the bootstrap errors reported in that table.
    masses                  =[10.05d0, 10.15d0, 10.25d0, 10.35d0, 10.45d0, 10.55d0, 10.65d0, 10.75d0, 10.85d0, 10.95d0, 11.05d0, 11.15d0, 11.25d0, 11.35d0, 11.45d0, 11.55d0,  11.65d0,  11.75d0,  11.85d0,  11.95d0]
    functionValueTarget     =[ 0.47d0,  0.79d0,  0.44d0,  0.90d0,  0.05d0,  0.65d0,  0.80d0,  0.25d0, -0.05d0,  1.54d0,  1.27d0,  1.70d0,  5.20d0, 11.20d0, 29.00d0, 60.70d0, 123.00d0, 266.00d0, 445.00d0, 721.00d0]
    functionErrorTarget     =[ 0.44d0,  0.39d0,  0.37d0,  0.37d0,  0.34d0,  0.37d0,  0.40d0,  0.43d0,  0.75d0,  0.58d0,  0.78d0,  1.10d0,  1.80d0,  2.40d0,  3.80d0,  6.80d0,  16.00d0,  36.00d0,  84.00d0, 210.00d0]
    masses                  =1.0d+1**masses
    functionValueTarget     =1.0d-6* functionValueTarget
    functionErrorTarget     =1.0d-6* functionErrorTarget
    functionCovarianceTarget=0.0d0
    do iBin=1_c_size_t,binCount
       functionCovarianceTarget(iBin,iBin)=functionErrorTarget(iBin)**2
    end do
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
     <constructor>
      cosmologyParametersSimple(                            &amp;
        &amp;                   OmegaMatter    = 0.27200d0, &amp;
        &amp;                   OmegaDarkEnergy= 0.72800d0, &amp;
        &amp;                   HubbleConstant =70.40000d0, &amp;
        &amp;                   temperatureCMB = 2.72548d0, &amp;
        &amp;                   OmegaBaryon    = 0.00000d0  &amp;
        &amp;                  )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
    <constructor>
      cosmologyFunctionsMatterLambda(                        &amp;
        &amp;                        cosmologyParametersData &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Construct survey geometry. We use the Li & White (2009) SDSS geometry here, since Planck Intermediate Results XI uses the
    ! NYU VAGC catalog (as do Li & White), with a depth of r=17.7 (while Li & White use r=17.6 - sufficiently close). Planck
    ! Intermediate Results XI uses a minimum redshift of 0.03, while the maximum redshift of 0.3 here is chosen based on the
    ! distributions of redshift shown in Figure 1 of that paper.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryLiWhite2009SDSS(redshiftMinimum=0.03d0,redshiftMaximum=0.30d0,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(binCount,outputTimes_%count()))
    do iBin=1,binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,outputTimes_,masses(iBin),allowSingleEpoch=.true.)
    end do
    ! Build a filter which select galaxies with stellar mass above some coarse lower limit suitable for this sample.
    allocate(galacticFilterCentrals_        )
    allocate(galacticFilterStellarMass_     )
    allocate(galacticFilter_                )
    allocate(filters_                       )
    allocate(filters_                  %next)
    filters_     %filter_ => galacticFilterCentrals_
    filters_%next%filter_ => galacticFilterStellarMass_
    !![
    <referenceConstruct object="galacticFilterCentrals_"    constructor="galacticFilterHaloIsolated(                      )"/>
    <referenceConstruct object="galacticFilterStellarMass_" constructor="galacticFilterStellarMass (massThreshold=1.0d9   )"/>
    <referenceConstruct object="galacticFilter_"            constructor="galacticFilterAll         (              filters_)"/>
    !!]
    ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                         )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                    constructor="outputAnalysisWeightOperatorIdentity           (                                                             )"/>
    !!]
    ! Build log10() property operator.
    allocate(outputAnalysisPropertyOperatorSystmtcPolynomial_      )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)"/>
    !!]
    allocate(outputAnalysisPropertyOperatorLog10_                  )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                             )"/>
    !!]
    allocate(propertyOperators_                                    )
    allocate(propertyOperators_%next                               )
    propertyOperators_     %operator_  => outputAnalysisPropertyOperatorLog10_
    propertyOperators_%next%operator_  => outputAnalysisPropertyOperatorSystmtcPolynomial_
    allocate(outputAnalysisPropertyOperator_                       )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence         (propertyOperators_                                           )"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
     outputAnalysisDistributionOperatorRandomErrorPlynml (                                  &amp;
        &amp;                                             randomErrorMinimum              , &amp;
        &amp;                                             randomErrorMaximum              , &amp;
        &amp;                                             errorPolynomialZeroPoint        , &amp;
        &amp;                                             randomErrorPolynomialCoefficient  &amp;
        &amp;                                            )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build a weight property operator. We use an operator which performs the cosmological scaling (e.g. unnumbered equation at
    ! the end of page 2 of Planck Intermediate Results XI), using the cosmology employed in their analysis.
    allocate(outputAnalysisWeightPropertyOperator_ )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_">
     <constructor>
     outputAnalysisPropertyOperatorCosmologySZ(cosmologyParametersData,cosmologyFunctionsData,outputTimes_)
     </constructor>
    </referenceConstruct>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_     )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_" constructor="outputAnalysisPropertyOperatorAntiLog10(                  )"/>
    !!]
    ! Create a stellar mass property extractor.
    allocate(nodePropertyExtractor_                )
    !![
    <referenceConstruct object="nodePropertyExtractor_"            constructor="nodePropertyExtractorMassStellar       ()"/>
    !!]
    ! Create a thermal Sunyaev-Zeldovich property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_">
     <constructor>
      nodePropertyExtractorICMSZ(                                                       &amp;
        &amp;                                              cosmologyParameters_       , &amp;
        &amp;                                              cosmologyFunctions_        , &amp;
        &amp;                                              darkMatterHaloScale_       , &amp;
        &amp;                                              chemicalState_             , &amp;
        &amp;                    densityContrast          =500.0d0                    , &amp;
        &amp;                    densityContrastRelativeTo=densityCosmologicalCritical, &amp;
        &amp;                    distanceAngular          =500.0d0                      &amp;
        &amp;                   )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                                     &
         &                                                         var_str('sunyaevZeldovichPlanck2013'                              ), &
         &                                                         var_str('Sunyaev-Zeldovich signal vs. central galaxy stellar mass'), &
         &                                                         var_str('massStellar'                                             ), &
         &                                                         var_str('Stellar mass'                                            ), &
         &                                                         var_str('M☉'                                                      ), &
         &                                                         massSolar                                                          , &
         &                                                         var_str('thermalSZ'                                               ), &
         &                                                         var_str('Thermal S-Z signal (̃Y₅₀₀)'                              ), &
         &                                                         var_str(' '                                                       ), &
         &                                                         0.0d0                                                              , &
         &                                                         log10(masses)                                                      , &
         &                                                         bufferCount                                                        , &
         &                                                         outputWeight                                                       , &
         &                                                         nodePropertyExtractor_                                             , &
         &                                                         outputAnalysisWeightPropertyExtractor_                             , &
         &                                                         outputAnalysisPropertyOperator_                                    , &
         &                                                         outputAnalysisWeightPropertyOperator_                              , &
         &                                                         outputAnalysisPropertyUnoperator_                                  , &
         &                                                         outputAnalysisWeightOperator_                                      , &
         &                                                         outputAnalysisDistributionOperator_                                , &
         &                                                         galacticFilter_                                                    , &
         &                                                         outputTimes_                                                       , &
         &                                                         outputAnalysisCovarianceModelBinomial                              , &
         &                                                         covarianceBinomialBinsPerDecade                                    , &
         &                                                         covarianceBinomialMassHaloMinimum                                  , &
         &                                                         covarianceBinomialMassHaloMaximum                                  , &
         &                                                         likelihoodNormalize                                                , &
         &                                                         var_str('$M_\star/\mathrm{M}_\odot$'                              ), &
         &                                                         var_str('$\widetilde{Y}_{500}/\hbox{arcmin}^2$'                   ), &
         &                                                         .true.                                                             , &
         &                                                         .true.                                                             , &
         &                                                         var_str('Planck Intermediate Results XI (2013)'                   ), &
         &                                                         functionValueTarget                                                , &
         &                                                         functionCovarianceTarget                                             &
         &                                                        )
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                                 />
    <objectDestructor name="cosmologyParametersData"                         />
    <objectDestructor name="cosmologyFunctionsData"                          />
    <objectDestructor name="galacticFilter_"                                 />
    <objectDestructor name="galacticFilterStellarMass_"                      />
    <objectDestructor name="galacticFilterCentrals_"                         />
    <objectDestructor name="outputAnalysisDistributionOperator_"             />
    <objectDestructor name="outputAnalysisWeightOperator_"                   />
    <objectDestructor name="outputAnalysisPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisPropertyUnoperator_"               />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"           />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"          />
    <objectDestructor name="nodePropertyExtractor_"                          />
    !!]
    nullify(propertyOperators_)
    nullify(filters_          )
    return
  end function sunyaevZeldovichPlanck2013ConstructorInternal

  subroutine sunyaevZeldovichPlanck2013Destructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisSunyaevZeldovichPlanck2013} output analysis class.
    !!}
    implicit none
    type(outputAnalysisSunyaevZeldovichPlanck2013), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%chemicalState_"      />
    <objectDestructor name="self%outputTimes_"        />
    !!]
    return
  end subroutine sunyaevZeldovichPlanck2013Destructor
