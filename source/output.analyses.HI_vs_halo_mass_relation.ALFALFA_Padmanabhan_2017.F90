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
  Implements an HI vs halo mass relation analysis class.
  !!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <outputAnalysis name="outputAnalysisHIVsHaloMassRelationPadmanabhan2017">
   <description>An HI vs halo mass relation output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisHIVsHaloMassRelationPadmanabhan2017
     !!{
     An HI vs halo mass relation output analysis class.
     !!}
     private
     class           (cosmologyParametersClass         ), pointer                     :: cosmologyParameters_                 => null()
     class           (cosmologyFunctionsClass          ), pointer                     :: cosmologyFunctions_                  => null()
     class           (outputAnalysisMolecularRatioClass), pointer                     :: outputAnalysisMolecularRatio_        => null()
     class           (virialDensityContrastClass       ), pointer                     :: virialDensityContrast_               => null()
     class           (darkMatterProfileDMOClass        ), pointer                     :: darkMatterProfileDMO_                => null()
     double precision                                   , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient
     integer         (c_size_t                         )                              :: likelihoodBin
   contains
     final :: hiVsHaloMassRelationPadmanabhan2017Destructor
  end type outputAnalysisHIVsHaloMassRelationPadmanabhan2017
  
  interface outputAnalysisHIVsHaloMassRelationPadmanabhan2017
     !!{
     Constructors for the \refClass{outputAnalysisHIVsHaloMassRelationPadmanabhan2017} output analysis class.
     !!}
     module procedure hiVsHaloMassRelationPadmanabhan2017ConstructorParameters
     module procedure hiVsHaloMassRelationPadmanabhan2017ConstructorInternal
  end interface outputAnalysisHIVsHaloMassRelationPadmanabhan2017

contains

  function hiVsHaloMassRelationPadmanabhan2017ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisHIVsHaloMassRelationPadmanabhan2017} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions             , only : cosmologyFunctions                                      , cosmologyFunctionsClass
    use :: Cosmology_Parameters            , only : cosmologyParameters                                     , cosmologyParametersClass
    use :: Functions_Global                , only : Virial_Density_Contrast_Percolation_Objects_Constructor_
    use :: Input_Parameters                , only : inputParameter                                          , inputParameters
    use :: Output_Analysis_Molecular_Ratios, only : outputAnalysisMolecularRatio                            , outputAnalysisMolecularRatioClass
    implicit none
    type            (outputAnalysisHIVsHaloMassRelationPadmanabhan2017)                              :: self
    type            (inputParameters                                  ), intent(inout)               :: parameters
    double precision                                                   , allocatable  , dimension(:) :: systematicErrorPolynomialCoefficient
    class           (cosmologyFunctionsClass                          ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                                 ), pointer                     :: outputTimes_
    class           (cosmologyParametersClass                         ), pointer                     :: cosmologyParameters_
    class           (virialDensityContrastClass                       ), pointer                     :: virialDensityContrast_
    class           (outputAnalysisMolecularRatioClass                ), pointer                     :: outputAnalysisMolecularRatio_
    class           (darkMatterProfileDMOClass                        ), pointer                     :: darkMatterProfileDMO_
    class           (*                                                ), pointer                     :: percolationObjects_
    integer         (c_size_t                                         )                              :: likelihoodBin

    ! Check and read parameters.
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                      ))
    end if
    !![
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for HI vs halo mass relation.</description>
    </inputParameter>
    <inputParameter>
      <name>likelihoodBin</name>
      <source>parameters</source>
      <defaultValue>0_c_size_t</defaultValue>
      <description>If $>0$ then use only the mass bin given by this value in the likelihood calculation.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"          name="cosmologyParameters_"          source="parameters"/>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="outputTimes"                  name="outputTimes_"                  source="parameters"/>
    <objectBuilder class="outputAnalysisMolecularRatio" name="outputAnalysisMolecularRatio_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="virialDensityContrast"        name="virialDensityContrast_"        source="parameters"/>
    !!]
    percolationObjects_ => Virial_Density_Contrast_Percolation_Objects_Constructor_(parameters)
    self                =  outputAnalysisHIVsHaloMassRelationPadmanabhan2017(likelihoodBin,systematicErrorPolynomialCoefficient,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,outputAnalysisMolecularRatio_,outputTimes_,percolationObjects_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"         />
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="outputTimes_"                 />
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="outputAnalysisMolecularRatio_"/>
    <objectDestructor name="virialDensityContrast_"       />
    !!]
    return
  end function hiVsHaloMassRelationPadmanabhan2017ConstructorParameters

  function hiVsHaloMassRelationPadmanabhan2017ConstructorInternal(likelihoodBin,systematicErrorPolynomialCoefficient,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,outputAnalysisMolecularRatio_,outputTimes_,percolationObjects_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisHIVsHaloMassRelationPadmanabhan2017} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                           , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersClass                          , cosmologyParametersSimple                      , hubbleUnitsLittleH
    use :: Dark_Matter_Halo_Scales               , only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Galactic_Filters                      , only : filterList                                        , galacticFilterAll                              , galacticFilterHaloIsolated                     , galacticFilterISMMass
    use :: Error                                 , only : Error_Report
    use :: Galacticus_Nodes                      , only : nodeComponentBasic                                , treeNode
    use :: Geometry_Surveys                      , only : surveyGeometryMartin2010ALFALFA
    use :: ISO_Varying_String                    , only : var_str
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassHalo                     , nodePropertyExtractorMassISM
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Numerical_Ranges                      , only : Make_Range                                        , rangeTypeLinear
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorIdentity
    use :: Output_Analysis_Molecular_Ratios      , only : outputAnalysisMolecularRatioClass
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10           , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorFilterHighPass   , outputAnalysisPropertyOperatorHIMass, &
          &                                               outputAnalysisPropertyOperatorLog10               , outputAnalysisPropertyOperatorSequence         , outputAnalysisPropertyOperatorSystmtcPolynomial, propertyOperatorList
    use :: Output_Analysis_Utilities             , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                          , only : outputTimesClass
    use :: Virial_Density_Contrast               , only : virialDensityContrastBryanNorman1998              , virialDensityContrastPercolation
    implicit none
    type            (outputAnalysisHIVsHaloMassRelationPadmanabhan2017 )                                :: self
    integer         (c_size_t                                          ), intent(in   )                 :: likelihoodBin
    double precision                                                    , intent(in   ), dimension(:  ) :: systematicErrorPolynomialCoefficient
    class           (cosmologyParametersClass                          ), intent(inout), target         :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), intent(inout), target         :: cosmologyFunctions_
    class           (virialDensityContrastClass                        ), intent(in   ), target         :: virialDensityContrast_
    class           (outputTimesClass                                  ), intent(inout), target         :: outputTimes_
    class           (outputAnalysisMolecularRatioClass                 ), intent(in   ), target         :: outputAnalysisMolecularRatio_
    class           (darkMatterProfileDMOClass                         ), intent(inout), target         :: darkMatterProfileDMO_
    class           (*                                                 ), intent(in   ), target         :: percolationObjects_
    integer         (c_size_t                                          ), parameter                     :: massHaloCount                                         =26
    double precision                                                    , parameter                     :: massHaloMinimum                                       = 1.0d10, massHaloMaximum                                     =1.0d15
    integer                                                             , parameter                     :: covarianceBinomialBinsPerDecade                       =10
    double precision                                                    , parameter                     :: covarianceBinomialMassHaloMinimum                     = 1.0d08        , covarianceBinomialMassHaloMaximum           =1.0d16
    double precision                                                    , allocatable  , dimension(:  ) :: massHalo                                                              , massHILogarithmicTarget                              , &
         &                                                                                                 massHaloLogarithmic
    double precision                                                    , allocatable  , dimension(:,:) :: outputWeight                                                          , massHILogarithmicCovarianceTarget
    type            (galacticFilterISMMass                             ), pointer                       :: galacticFilterISMMass_
    type            (galacticFilterHaloIsolated                        ), pointer                       :: galacticFilterHaloIsolated_
    type            (galacticFilterAll                                 ), pointer                       :: galacticFilterAll_
    type            (filterList                                        ), pointer                       :: filters_
    type            (outputAnalysisDistributionOperatorIdentity        ), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity              ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorHIMass              ), pointer                       :: outputAnalysisWeightPropertyOperatorHIMass_
    type            (outputAnalysisPropertyOperatorLog10               ), pointer                       :: outputAnalysisPropertyOperator_                                      , outputAnalysisWeightPropertyOperatorLog10_            , &
         &                                                                                                 outputAnalysisWeightPropertyOperatorLog10Second_
    type            (outputAnalysisPropertyOperatorAntiLog10           ), pointer                       :: outputAnalysisPropertyUnoperator_                                    , outputAnalysisWeightPropertyOperatorAntiLog10_
    type            (outputAnalysisPropertyOperatorSequence            ), pointer                       :: outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc   ), pointer                       :: outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial   ), pointer                       :: outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    type            (outputAnalysisPropertyOperatorFilterHighPass      ), pointer                       :: outputAnalysisWeightPropertyOperatorFilterHighPass_
    type            (nodePropertyExtractorMassHalo                     ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorMassISM                      ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (propertyOperatorList                              ), pointer                       :: propertyOperators_
    type            (cosmologyParametersSimple                         ), pointer                       :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                    ), pointer                       :: cosmologyFunctionsData
    type            (virialDensityContrastPercolation                  ), pointer                       :: virialDensityContrastDefinition_
    type            (virialDensityContrastBryanNorman1998              ), pointer                       :: virialDensityContrastData
    type            (darkMatterHaloScaleVirialDensityContrastDefinition), pointer                       :: darkMatterHaloScaleData
    type            (surveyGeometryMartin2010ALFALFA                   ), pointer                       :: surveyGeometry_
    type            (treeNode                                          ), pointer                       :: nodeWork
    class           (nodeComponentBasic                                ), pointer                       :: basicWork
    double precision                                                    , parameter                     :: errorPolynomialZeroPoint                              =10.0d00
    double precision                                                    , parameter                     :: massHILimit                                           = 1.0d06
    logical                                                             , parameter                     :: likelihoodNormalize                                   =.false.
    double precision                                                    , parameter                     :: alphaFit                                              =0.09d0        , betaFit                                      =-0.58d0, &
         &                                                                                                 log10Velocity0Fit                                     =1.56d0        , log10Velocity1Fit                            =+4.64d0, &
         &                                                                                                 massReferenceFit                                      =1.0d11        , errorAlphaFit                                =+0.02d0, &
         &                                                                                                 errorBetaFit                                          =0.12d0        , errorLog10Velocity0Fit                       =+0.03d0, &
         &                                                                                                 errorLog10Velocity1Fit                                =0.75d0
    double precision                                                                                    :: velocityVirial                                                       , fractionHydrogenCosmic                               , &
         &                                                                                                 velocity0Fit                                                         , velocity1Fit                                         , &
         &                                                                                                 jacobianAlpha                                                        , jacobianBeta                                         , &
         &                                                                                                 jacobianVelocity0                                                    , jacobianVelocity1
    integer         (c_size_t                                          )                                :: iBin
    !![
    <constructorAssign variables="systematicErrorPolynomialCoefficient, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *virialDensityContrast_, *outputAnalysisMolecularRatio_"/>
    !!]
    
    ! Construct survey geometry.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryMartin2010ALFALFA(cosmologyParameters_)"/>
    !!]
    ! Create output time weights.
    allocate(outputWeight(massHaloCount,outputTimes_%count()))
    outputWeight(1,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,outputTimes_,massHILimit,allowSingleEpoch=.true.)
    forall(iBin=2:massHaloCount)
       outputWeight(iBin,:)=outputWeight(1,:)
    end forall
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
     <constructor>
      cosmologyParametersSimple     (                            &amp;
         &amp;                       OmegaMatter    = 0.28100d0, &amp;
         &amp;                       OmegaDarkEnergy= 0.71900d0, &amp;
         &amp;                       HubbleConstant =71.00000d0, &amp;
         &amp;                       temperatureCMB = 2.72548d0, &amp;
         &amp;                       OmegaBaryon    = 0.04620d0  &amp;
         &amp;                      )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
     <constructor>
      cosmologyFunctionsMatterLambda(                            &amp;
         &amp;                       cosmologyParametersData     &amp;
         &amp;                      )
     </constructor>
    </referenceConstruct>
    !!]
    ! Create bins in halo mass.
    massHaloLogarithmic=Make_Range(log10(massHaloMinimum),log10(massHaloMaximum),int(massHaloCount),rangeType=rangeTypeLinear)
    massHalo           =10.0d0**massHaloLogarithmic
    ! Build a filter which select central galaxies with HI mass above some coarse lower limit suitable for this sample.
    allocate(galacticFilterISMMass_          )
    allocate(galacticFilterHaloIsolated_     )
    allocate(galacticFilterAll_              )
    allocate(filters_                        )
    allocate(filters_                   %next)
    filters_     %filter_ => galacticFilterHaloIsolated_
    filters_%next%filter_ => galacticFilterISMMass_
    !![
    <referenceConstruct object="galacticFilterISMMass_"      constructor="galacticFilterISMMass      (massThreshold=1.0d4   )"/>
    <referenceConstruct object="galacticFilterHaloIsolated_" constructor="galacticFilterHaloIsolated (                      )"/>
    <referenceConstruct object="galacticFilterAll_"          constructor="galacticFilterAll          (              filters_)"/>
    !!]
    ! Build identity distribution operator.
    allocate(outputAnalysisDistributionOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_"                    constructor="outputAnalysisDistributionOperatorIdentity            (                                                                                  )"/>
    !!]
    ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_                         )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                          constructor="outputAnalysisWeightOperatorIdentity                  (                                                                                  )"/>
    !!]
    ! Build log10() property operator.
    allocate(outputAnalysisPropertyOperator_                       )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"                        constructor="outputAnalysisPropertyOperatorLog10                   (                                                                                  )"/>
    !!]
    ! Build a sequence (HI mass, log10, polynomial systematic, anti-log10, cosmological luminosity distance, high-pass filter) of weight property operators.
    allocate(outputAnalysisWeightPropertyOperatorFilterHighPass_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorFilterHighPass_"    constructor="outputAnalysisPropertyOperatorFilterHighPass          (log10(massHILimit)                                                                )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc       (cosmologyFunctions_          ,cosmologyFunctionsData              ,outputTimes_   )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorSystmtcPolynomial_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial       (errorPolynomialZeroPoint     ,systematicErrorPolynomialCoefficient                )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorHIMass_           )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorHIMass_"            constructor="outputAnalysisPropertyOperatorHIMass                  (outputAnalysisMolecularRatio_                                                     )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10                   (                                                                                  )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorLog10Second_      )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorLog10Second_"       constructor="outputAnalysisPropertyOperatorLog10                   (                                                                                  )"/>
    !!]
    allocate(outputAnalysisWeightPropertyOperatorAntiLog10_        )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorAntiLog10_"         constructor="outputAnalysisPropertyOperatorAntiLog10               (                                                                                  )"/>
    !!]
    allocate(propertyOperators_                              )
    allocate(propertyOperators_%next                         )
    allocate(propertyOperators_%next%next                    )
    allocate(propertyOperators_%next%next%next               )
    allocate(propertyOperators_%next%next%next%next          )
    allocate(propertyOperators_%next%next%next%next%next     )
    allocate(propertyOperators_%next%next%next%next%next%next)
    propertyOperators_                              %operator_ => outputAnalysisWeightPropertyOperatorHIMass_
    propertyOperators_%next                         %operator_ => outputAnalysisWeightPropertyOperatorLog10_
    propertyOperators_%next%next                    %operator_ => outputAnalysisWeightPropertyOperatorSystmtcPolynomial_
    propertyOperators_%next%next%next               %operator_ => outputAnalysisWeightPropertyOperatorAntiLog10_
    propertyOperators_%next%next%next%next          %operator_ => outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_
    propertyOperators_%next%next%next%next%next     %operator_ => outputAnalysisWeightPropertyOperatorLog10Second_
    propertyOperators_%next%next%next%next%next%next%operator_ => outputAnalysisWeightPropertyOperatorFilterHighPass_
    allocate(outputAnalysisWeightPropertyOperator_                 )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_"                  constructor="outputAnalysisPropertyOperatorSequence                (propertyOperators_                                                                                                             )"/>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_                     )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_"                      constructor="outputAnalysisPropertyOperatorAntiLog10               (                                                                                                                               )"/>
    !!]
    ! Create an HI mass weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_                )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"                 constructor="nodePropertyExtractorMassISM                          (                                                                                                                               )"/>
    !!]
    ! Create a halo mass weight property extractor. The virial density contrast is chosen to equal that expected for a
    ! friends-of-friends algorithm with linking length parameter b=0.2 since that is what was used by Sheth, Mo & Tormen (2001) in
    ! their original calibration of their halo mass function (as used by Padmanabhan & Refregier 2017).
    allocate(virialDensityContrastDefinition_                                )
    !![
    <referenceConstruct object="virialDensityContrastDefinition_"                       constructor="virialDensityContrastPercolation                      (0.2d0                        ,cosmologyFunctions_       ,percolationObjects_                                                  )"/>
    !!]
    allocate(nodePropertyExtractor_                       )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                                 constructor="nodePropertyExtractorMassHalo                         (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_)"/>
    !!]
    ! Create a halo scale object from which to compute virial velocities. Padmanabhan & Refrigier use the Bryan & Norman (1998)
    ! virial density contrast definition. However (Padmanabhan, private communication), they assume it gives the density contrast
    ! relative to mean density, instead of critical density (as was used by Bryan & Norman). This means that virial velocities are
    ! lower by a factor Ωₘ^⅙ than they should be. We explicitly account for this factor when computing virial velocity below.
    allocate(virialDensityContrastData                             )
    !![
    <referenceConstruct object="virialDensityContrastData"                              constructor="virialDensityContrastBryanNorman1998                  (.false.,cosmologyParametersData      ,cosmologyFunctionsData                              )"/>
    !!]
    allocate(darkMatterHaloScaleData)
    !![
    <referenceConstruct object="darkMatterHaloScaleData"                                constructor="darkMatterHaloScaleVirialDensityContrastDefinition    (       cosmologyParametersData      ,cosmologyFunctionsData    ,virialDensityContrastData)"/>
    !!]
    ! Generate the target dataset.
    allocate(massHILogarithmicTarget          (massHaloCount              ))
    allocate(massHILogarithmicCovarianceTarget(massHaloCount,massHaloCount))
    nodeWork  => treeNode      (                 )
    basicWork => nodeWork%basic(autoCreate=.true.)
    call basicWork%timeSet            (cosmologyFunctionsData%cosmicTime(expansionFactor=1.0d0))
    call basicWork%timeLastIsolatedSet(cosmologyFunctionsData%cosmicTime(expansionFactor=1.0d0))
    fractionHydrogenCosmic           =+cosmologyParametersData%OmegaBaryon() &
         &                            /cosmologyParametersData%OmegaMatter()
    massHILogarithmicCovarianceTarget=+0.0d0
    velocity0Fit                     =+10.0d0**log10Velocity0Fit
    velocity1Fit                     =+10.0d0**log10Velocity1Fit
    do iBin=1,massHaloCount
       call basicWork%massSet(massHalo(iBin))
       ! Compute virial velocity, including the extra factor of Ωₘ^⅙ to correct for incorrect definition of virial density
       ! contrast used by Padmanabhan & Refrigier (2017).
       velocityVirial               =+darkMatterHaloScaleData%velocityVirial(nodeWork)                &
            &                        *cosmologyParametersData%OmegaMatter   (        )**(1.0d0/6.0d0)
       massHILogarithmicTarget(iBin)=+alphaFit                                                       &
            &                        *fractionHydrogenCosmic                                         &
            &                        *  massHalo(iBin)                                               &
            &                        *(                                                              &
            &                          +  massHalo(iBin)                                             &
            &                          /(                                                            &
            &                            +massReferenceFit                                           &
            &                            /cosmologyParametersData%HubbleConstant(hubbleUnitsLittleH) &
            &                           )                                                            &
            &                         )**betaFit                                                     &
            &                        *exp(-(velocity0Fit  /velocityVirial)**3)                       &
            &                        *exp(-(velocityVirial/velocity1Fit  )**3)
       jacobianAlpha    =+massHILogarithmicTarget(iBin)/alphaFit
       jacobianBeta     =+massHILogarithmicTarget(iBin)*log(massHalo(iBin)/(massReferenceFit/cosmologyParametersData%HubbleConstant(hubbleUnitsLittleH)))
       jacobianVelocity0=-massHILogarithmicTarget(iBin)*3.0d0*log(10.0d0)*velocity0Fit  **3/velocityVirial**3
       jacobianVelocity1=+massHILogarithmicTarget(iBin)*3.0d0*log(10.0d0)*velocityVirial**3/velocity1Fit  **3
       massHILogarithmicCovarianceTarget(iBin,iBin)=+jacobianAlpha    **2*errorAlphaFit         **2 &
            &                                       +jacobianBeta     **2*errorBetaFit          **2 &
            &                                       +jacobianVelocity0**2*errorLog10Velocity0Fit**2 &
            &                                       +jacobianVelocity1**2*errorLog10Velocity1Fit**2
       ! Convert to log10.
       massHILogarithmicCovarianceTarget(iBin,iBin)=+      massHILogarithmicCovarianceTarget(iBin,iBin)     &
            &                                       /      massHILogarithmicTarget          (iBin     ) **2 &
            &                                       /log  (10.0d0                                      )**2
       massHILogarithmicTarget          (iBin     )=+log10(massHILogarithmicTarget          (iBin     ))
    end do
    call nodeWork%destroy()
    deallocate(nodeWork)
    self%likelihoodBin=likelihoodBin
    if (self%likelihoodBin > 0_c_size_t) then
       ! Assume that only a single bin of the relation is to be populated. Set the target dataset in all other bins to zero so
       ! that they do not contribute to the likelihood.
       if (self%likelihoodBin > massHaloCount) call Error_Report('likelihoodBin is out of range'//{introspection:location})
       do iBin=1,massHaloCount
          if (iBin /= self%likelihoodBin) massHILogarithmicTarget(iBin)=0.0d0
       end do
    end if
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                        &
         &                                                         var_str('hiHaloMassRelationPadmanabhan2017'          ), &
         &                                                         var_str('HI vs. halo mass relation'                  ), &
         &                                                         var_str('massHalo'                                   ), &
         &                                                         var_str('Halo mass'                                  ), &
         &                                                         var_str('M☉'                                        ), &
         &                                                         massSolar                                             , &
         &                                                         var_str('massHILog10'                                ), &
         &                                                         var_str('⟨log₁₀(HI mass/M☉)⟩'                        ), &
         &                                                         var_str(' '                                          ), &
         &                                                         0.0d0                                                 , &
         &                                                         massHaloLogarithmic                                   , &
         &                                                         0_c_size_t                                            , &
         &                                                         outputWeight                                          , &
         &                                                         nodePropertyExtractor_                                , &
         &                                                         outputAnalysisWeightPropertyExtractor_                , &
         &                                                         outputAnalysisPropertyOperator_                       , &
         &                                                         outputAnalysisWeightPropertyOperator_                 , &
         &                                                         outputAnalysisPropertyUnoperator_                     , &
         &                                                         outputAnalysisWeightOperator_                         , &
         &                                                         outputAnalysisDistributionOperator_                   , &
         &                                                         galacticFilterAll_                                    , &
         &                                                         outputTimes_                                          , &
         &                                                         outputAnalysisCovarianceModelBinomial                 , &
         &                                                         covarianceBinomialBinsPerDecade                       , &
         &                                                         covarianceBinomialMassHaloMinimum                     , &
         &                                                         covarianceBinomialMassHaloMaximum                     , &
         &                                                         likelihoodNormalize                                   , &
         &                                                         var_str('$M_\mathrm{halo}/\mathrm{M}_\odot$'         ), &
         &                                                         var_str('$\log_{10}(M_\mathrm{HI}/\mathrm{M}_\odot)$'), &
         &                                                         .true.                                                , &
         &                                                         .false.                                               , &
         &                                                         var_str('Padmanabhan \\& Refrigier (2017)'           ), &
         &                                                         massHILogarithmicTarget                               , &
         &                                                         massHILogarithmicCovarianceTarget                       &
         &                                                        )
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                                       />
    <objectDestructor name="galacticFilterAll_"                                    />
    <objectDestructor name="galacticFilterISMMass_"                                />
    <objectDestructor name="galacticFilterHaloIsolated_"                           />
    <objectDestructor name="outputAnalysisDistributionOperator_"                   />
    <objectDestructor name="outputAnalysisWeightOperator_"                         />
    <objectDestructor name="outputAnalysisPropertyOperator_"                       />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"                     />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"                 />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"                />
    <objectDestructor name="nodePropertyExtractor_"                                />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorFilterHighPass_"   />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorCsmlgyLmnstyDstnc_"/>
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSystmtcPolynomial_"/>
    <objectDestructor name="outputAnalysisWeightPropertyOperatorHIMass_"           />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorLog10Second_"      />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorAntiLog10_"        />
    <objectDestructor name="cosmologyParametersData"                               />
    <objectDestructor name="cosmologyFunctionsData"                                />
    <objectDestructor name="virialDensityContrastDefinition_"                      />
    <objectDestructor name="virialDensityContrastData"                             />
    <objectDestructor name="darkMatterHaloScaleData"                               />
    !!]
    nullify(filters_          )
    nullify(propertyOperators_)
    return
  end function hiVsHaloMassRelationPadmanabhan2017ConstructorInternal

  subroutine hiVsHaloMassRelationPadmanabhan2017Destructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisHIVsHaloMassRelationPadmanabhan2017} output analysis class.
    !!}
    implicit none
    type(outputAnalysisHIVsHaloMassRelationPadmanabhan2017), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"         />
    <objectDestructor name="self%cosmologyFunctions_"          />
    <objectDestructor name="self%outputAnalysisMolecularRatio_"/>
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    <objectDestructor name="self%virialDensityContrast_"       />
    !!]
    return
  end subroutine hiVsHaloMassRelationPadmanabhan2017Destructor
