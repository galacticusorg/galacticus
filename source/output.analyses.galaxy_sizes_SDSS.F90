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
Implements a galaxy size output analysis class for SDSS data.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <outputAnalysis name="outputAnalysisGalaxySizesSDSS">
   <description>
    An output analysis class which computes the mass-dependent $z\approx 0.07$ galaxy size distribution of \cite{shen_size_2003} from
    the \gls{sdss}. The size function reported by \cite{shen_size_2003} is converted to the appropriate cosmology for the given \glc\
    model (assuming that sizes scale as the angular diameter distance, and masses as the square of the luminosity distance). The model
    sizes and masses are then used to construct a mass-dependent radius function by binning into a 2-D histogram using the size and
    mass bins reported by \cite{shen_size_2003} (modified as described above) as the centers of the bins (with bin boundaries placed
    at the geometric means of consecutive bin centers). Distributions are computed for both late-type and early-type galaxies,
    classified on the basis of the stellar mass spheroid-to-total ratio, with the division at a ratio given by {\normalfont \ttfamily
    [massStellarRatio]}.
   </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisGalaxySizesSDSS
     !!{
     A galaxySizesSDSS output analysis class.
     !!}
     private
     integer                                              :: distributionNumber
     double precision                                     :: massStellarRatio               , sizeSourceLensing
     class           (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_   => null()
     class           (gravitationalLensingClass), pointer :: gravitationalLensing_ => null()
   contains
     final :: galaxySizesSDSSDestructor
  end type outputAnalysisGalaxySizesSDSS

  interface outputAnalysisGalaxySizesSDSS
     !!{
     Constructors for the \refClass{outputAnalysisGalaxySizesSDSS} output analysis class.
     !!}
     module procedure galaxySizesSDSSConstructorParameters
     module procedure galaxySizesSDSSConstructorInternal
  end interface outputAnalysisGalaxySizesSDSS

contains

  function galaxySizesSDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisGalaxySizesSDSS} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisGalaxySizesSDSS)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass             ), pointer       :: outputTimes_
    class           (gravitationalLensingClass    ), pointer       :: gravitationalLensing_
    double precision                                               :: massStellarRatio     , sizeSourceLensing
    integer                                                        :: distributionNumber

    !![
    <inputParameter>
      <name>distributionNumber</name>
      <source>parameters</source>
      <description>The number (1-34) of the distribution to compute.</description>
    </inputParameter>
    <inputParameter>
      <name>massStellarRatio</name>
      <source>parameters</source>
      <defaultValue>0.3d0</defaultValue>
      <description>The stellar mass bulge-to-total ratio used to discriminate late-type vs. early-type galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>sizeSourceLensing</name>
      <source>parameters</source>
      <variable>sizeSourceLensing</variable>
      <defaultValue>2.0d-3</defaultValue>
      <description>The characteristic source size for gravitational lensing calculations.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="outputTimes"          name="outputTimes_"          source="parameters"/>
    <objectBuilder class="gravitationalLensing" name="gravitationalLensing_" source="parameters"/>
    !!]
    self=outputAnalysisGalaxySizesSDSS(distributionNumber,massStellarRatio,sizeSourceLensing,cosmologyFunctions_,outputTimes_,gravitationalLensing_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="outputTimes_"         />
    <objectDestructor name="gravitationalLensing_"/>
    !!]
    return
  end function galaxySizesSDSSConstructorParameters

  function galaxySizesSDSSConstructorInternal(distributionNumber,massStellarRatio,sizeSourceLensing,cosmologyFunctions_,outputTimes_,gravitationalLensing_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisGalaxySizesSDSS} output analysis class.
    !!}
    use :: Cosmology_Functions                     , only : cosmologyFunctionsClass                      , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                    , only : cosmologyParametersSimple
    use :: Galactic_Filters                        , only : filterList                                   , galacticFilterAll                             , galacticFilterNot                              , galacticFilterStellarMass                        , &
          &                                                 galacticFilterStellarMassMorphology
    use :: Error                                   , only : Error_Report
    use :: Input_Paths                             , only : inputPath                                    , pathTypeDataStatic
    use :: Geometry_Surveys                        , only : surveyGeometryLiWhite2009SDSS
    use :: Gravitational_Lensing                   , only : gravitationalLensingClass
    use :: HDF5_Access                             , only : hdf5Access
    use :: IO_HDF5                                 , only : hdf5Object
    use :: ISO_Varying_String                      , only : var_str                                      , varying_string
    use :: Node_Property_Extractors                , only : nodePropertyExtractorRadiusHalfMassStellar   , nodePropertyExtractorMassStellar
    use :: Numerical_Constants_Astronomical        , only : megaParsec
    use :: Numerical_Constants_Prefixes            , only : kilo                                         , milli
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Normalizers, only : normalizerList                               , outputAnalysisDistributionNormalizerBinWidth  , outputAnalysisDistributionNormalizerSequence   , outputAnalysisDistributionNormalizerUnitarity
    use :: Output_Analysis_Distribution_Operators  , only : distributionOperatorList                     , lensedPropertySize                            , outputAnalysisDistributionOperatorClass        , outputAnalysisDistributionOperatorDiskSizeInclntn, &
          &                                                 outputAnalysisDistributionOperatorGrvtnlLnsng, outputAnalysisDistributionOperatorIdentity    , outputAnalysisDistributionOperatorSequence
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10      , outputAnalysisPropertyOperatorCsmlgyAnglrDstnc, outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorLog10              , &
          &                                                 outputAnalysisPropertyOperatorMultiply       , outputAnalysisPropertyOperatorSequence        , propertyOperatorList
    use :: Output_Analysis_Utilities               , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorNormal
    use :: Output_Times                            , only : outputTimesClass
    implicit none
    type            (outputAnalysisGalaxySizesSDSS                  )                              :: self
    integer                                                                       , intent(in   )  :: distributionNumber
    double precision                                                              , intent(in   )  :: massStellarRatio                               , sizeSourceLensing
    class           (cosmologyFunctionsClass                        ), target     , intent(in   )  :: cosmologyFunctions_
    class           (outputTimesClass                               ), target     , intent(inout)  :: outputTimes_
    class           (gravitationalLensingClass                      ), target     , intent(in   )  :: gravitationalLensing_
    type            (cosmologyParametersSimple                      ), pointer                     :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                 ), pointer                     :: cosmologyFunctionsData
    type            (nodePropertyExtractorRadiusHalfMassStellar     ), pointer                     :: nodePropertyExtractor_
    type            (nodePropertyExtractorMassStellar               ), pointer                     :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc ), pointer                     :: outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc), pointer                     :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorLog10            ), pointer                     :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorMultiply         ), pointer                     :: outputAnalysisPropertyOperatorMultiply_
    type            (outputAnalysisPropertyOperatorSequence         ), pointer                     :: outputAnalysisPropertyOperatorSequence_         , outputAnalysisWeightPropertyOperatorSequence_
    type            (outputAnalysisWeightOperatorNormal             ), pointer                     :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerSequence   ), pointer                     :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisPropertyOperatorAntiLog10        ), pointer                     :: outputAnalysisPropertyOperatorAntiLog10_
    type            (outputAnalysisDistributionOperatorSequence     ), pointer                     :: outputAnalysisDistributionOperator_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng  ), pointer                     :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    class           (outputAnalysisDistributionOperatorClass        ), pointer                     :: outputAnalysisDistributionOperatorProjection_
    type            (distributionOperatorList                       ), pointer                     :: distributionOperatorSequence
    type            (surveyGeometryLiWhite2009SDSS                  ), pointer                     :: surveyGeometry_
    type            (outputAnalysisDistributionNormalizerBinWidth   ), pointer                     :: outputAnalysisDistributionNormalizerBinWidth_
    type            (outputAnalysisDistributionNormalizerUnitarity  ), pointer                     :: outputAnalysisDistributionNormalizerUnitarity_
    type            (normalizerList                                 ), pointer                     :: normalizerSequence
    type            (propertyOperatorList                           ), pointer                     :: propertyOperatorSequence                        , weightPropertyOperatorSequence
    type            (galacticFilterStellarMass                      ), pointer                     :: galacticFilterMassStellarMinimum_               , galacticFilterMassStellarMaximum_
    type            (galacticFilterStellarMassMorphology            ), pointer                     :: galacticFilterMorphology_
    type            (galacticFilterNot                              ), pointer                     :: galacticFilterMassStellarMaximumInverted_       , galacticFilterMorphologyInverted_
    type            (galacticFilterAll                              ), pointer                     :: galacticFilterAll_
    type            (filterList                                     ), pointer                     :: filters_                                        , filter_
    double precision                                                 , allocatable, dimension(:  ) :: radii                                           , functionValueTarget                                , &
         &                                                                                            functionErrorTarget
    double precision                                                 , allocatable, dimension(:,:) :: outputWeight                                    , functionCovarianceTarget
    integer                                                          , parameter                   :: covarianceBinomialBinsPerDecade         =2
    double precision                                                 , parameter                   :: covarianceBinomialMassHaloMinimum       =3.00d11, covarianceBinomialMassHaloMaximum           =1.0d15
    !  Random error (in dex) on galaxy stellar masses in the Shen et al. (2003) sample. Shen et al. (2003) quote 95% confidence
    !  interval on masses of +/-40%, which corresponds to a standard deviation of 0.0806 dex assuming a normal distribution in
    !  log10(stellar mass).
    double precision                                                 , parameter                   :: massStellarErrorDex                     =8.06d-2
    integer         (c_size_t                                       )                              :: iBin
    type            (hdf5Object                                     )                              :: dataFile                                       , distribution
    double precision                                                                               :: massStellarMinimum                             , massStellarMaximum                                  , &
         &                                                                                            indexSersicMinimum                             , indexSersicMaximum
    character       (len=16                                         )                              :: distributionName                               , massStellarMinimumLogarithmic                       , &
         &                                                                                            massStellarMaximumLogarithmic
    type            (varying_string                                 )                              :: description
    logical                                                                                        :: isLateType
    !![
    <constructorAssign variables="distributionNumber, massStellarRatio, sizeSourceLensing, *cosmologyFunctions_, *gravitationalLensing_"/>
    !!]

    ! Validate input.
    if (distributionNumber < 1 .or. distributionNumber > 34) call Error_Report('distributionNumber ∈ [1..34] is required'//{introspection:location})
    ! Construct sizes matched to those used by  Shen et al. (2003). Also read stellar mass and Sersic index ranges.
    write (distributionName,'(a,i2.2)') 'distribution',distributionNumber
    !$ call hdf5Access%set()
    dataFile    =hdf5Object          (char(inputPath(pathTypeDataStatic)//'observations/galaxySizes/Galaxy_Sizes_By_Mass_SDSS_Shen_2003.hdf5'),readOnly=.true.)
    distribution=dataFile  %openGroup(distributionName                                                                                                        )
    call distribution%readDataset  ('radius'             ,radii              )
    call distribution%readDataset  ('radiusFunction'     ,functionValueTarget)
    call distribution%readDataset  ('radiusFunctionError',functionErrorTarget)
    call distribution%readAttribute('massMinimum'        ,massStellarMinimum )
    call distribution%readAttribute('massMaximum'        ,massStellarMaximum )
    call distribution%readAttribute('sersicIndexMinimum' ,indexSersicMinimum )
    call distribution%readAttribute('sersicIndexMaximum' ,indexSersicMaximum )
    !$ call hdf5Access%unset()
    self %binCount=size(radii)
    allocate(functionCovarianceTarget(self%binCount,self%binCount))
    functionCovarianceTarget=0.0d0
    do iBin=1,self%binCount
       ! Negative values indicate upper limits, which mean no galaxies were observed in the given bin.
       functionValueTarget     (iBin     )=max(0.0d0,functionValueTarget(iBin))
       functionCovarianceTarget(iBin,iBin)=          functionErrorTarget(iBin) **2
    end do
    ! Determine if this distribution is for late-type galaxies.
    isLateType=(indexSersicMinimum == 0.0d0)
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
     <constructor>
     cosmologyParametersSimple     (                            &amp;
          &amp;                     OmegaMatter    = 0.30000d0, &amp;
          &amp;                     OmegaDarkEnergy= 0.70000d0, &amp;
          &amp;                     HubbleConstant =70.00000d0, &amp;
          &amp;                     temperatureCMB = 2.72548d0, &amp;
          &amp;                     OmegaBaryon    = 0.04550d0  &amp;
          &amp;                    )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
     <constructor>
     cosmologyFunctionsMatterLambda(                            &amp;
          &amp;                     cosmologyParametersData     &amp;
          &amp;                    )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the SDSS survey geometry of Shen et al. (2013) with their imposed redshift limits.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryLiWhite2009SDSS(redshiftMinimum=1.0d-3,redshiftMaximum=huge(0.0d0),cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(self%binCount,outputTimes_%count()))
    do iBin=1,self%binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,self%cosmologyFunctions_,outputTimes_,massStellarMinimum)
    end do
    ! Create a half-mass radius property extractor.
    allocate(nodePropertyExtractor_        )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                           constructor="nodePropertyExtractorRadiusHalfMassStellar       (                                                                                                                                                             )"/>
    !!]
    ! Create a stellar mass property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_        )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_"           constructor="nodePropertyExtractorMassStellar                 (                                                                                                                                                             )"/>
    !!]
    ! Create multiply, log10, cosmological angular distance, and cosmological luminosity distance property operators.
    allocate(outputAnalysisPropertyOperatorMultiply_         )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorMultiply_"          constructor="outputAnalysisPropertyOperatorMultiply            (kilo                                                                                                                                                        )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10               (                                                                                                                                                            )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_"  constructor="outputAnalysisPropertyOperatorCsmlgyAnglrDstnc    (cosmologyFunctions_             ,cosmologyFunctionsData,outputTimes_                                                                                        )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc   (cosmologyFunctions_             ,cosmologyFunctionsData,outputTimes_                                                                                        )"/>
    !!]
    allocate(propertyOperatorSequence          )
    allocate(propertyOperatorSequence%next     )
    allocate(propertyOperatorSequence%next%next)
    propertyOperatorSequence          %operator_ => outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_
    propertyOperatorSequence%next     %operator_ => outputAnalysisPropertyOperatorMultiply_
    propertyOperatorSequence%next%next%operator_ => outputAnalysisPropertyOperatorLog10_
    allocate(outputAnalysisPropertyOperatorSequence_ )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorSequence_"          constructor="outputAnalysisPropertyOperatorSequence             (propertyOperatorSequence                                                                                                                                   )"/>
    !!]
    allocate(weightPropertyOperatorSequence          )
    allocate(weightPropertyOperatorSequence%next     )
    weightPropertyOperatorSequence     %operator_ => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    weightPropertyOperatorSequence%next%operator_ => outputAnalysisPropertyOperatorLog10_
    allocate(outputAnalysisWeightPropertyOperatorSequence_ )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSequence_"    constructor="outputAnalysisPropertyOperatorSequence             (weightPropertyOperatorSequence                                                                                                                              )"/>
    !!]
    ! Create a normal-weight weight operator.
    allocate(outputAnalysisWeightOperator_           )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                    constructor="outputAnalysisWeightOperatorNormal                 (log10(massStellarMinimum),log10(massStellarMaximum),massStellarErrorDex,outputAnalysisWeightPropertyExtractor_,outputAnalysisWeightPropertyOperatorSequence_)"/>
    !!]
    ! Create an projection distribution operator.
    if (isLateType) then
       allocate(outputAnalysisDistributionOperatorDiskSizeInclntn :: outputAnalysisDistributionOperatorProjection_)
       select type (outputAnalysisDistributionOperatorProjection_)
       type is (outputAnalysisDistributionOperatorDiskSizeInclntn)
          !![
          <referenceConstruct object="outputAnalysisDistributionOperatorProjection_" constructor="outputAnalysisDistributionOperatorDiskSizeInclntn(                                                                                                                                                           )"/>
          !!]
       end select
    else
       allocate(outputAnalysisDistributionOperatorIdentity        :: outputAnalysisDistributionOperatorProjection_)
       select type (outputAnalysisDistributionOperatorProjection_)
       type is (outputAnalysisDistributionOperatorIdentity)
          !![
          <referenceConstruct object="outputAnalysisDistributionOperatorProjection_" constructor="outputAnalysisDistributionOperatorIdentity       (                                                                                                                                                           )"/>
          !!]
       end select
    end if
    ! Create a gravitational lensing distribution operator.
    allocate(outputAnalysisDistributionOperatorGrvtnlLnsng_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperatorGrvtnlLnsng_">
    <constructor>
    outputAnalysisDistributionOperatorGrvtnlLnsng       (                                  &amp;
         &amp;                                           gravitationalLensing_           , &amp;
         &amp;                                           outputTimes_                    , &amp;
         &amp;                                           sizeSourceLensing               , &amp;
         &amp;                                           lensedPropertySize                &amp;
         &amp;                                          )
     </constructor>
    </referenceConstruct>
    !!]
    allocate(outputAnalysisDistributionOperator_     )
    allocate(distributionOperatorSequence            )
    allocate(distributionOperatorSequence       %next)
    distributionOperatorSequence     %operator_ => outputAnalysisDistributionOperatorProjection_
    distributionOperatorSequence%next%operator_ => outputAnalysisDistributionOperatorGrvtnlLnsng_
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
    <constructor>
    outputAnalysisDistributionOperatorSequence(                             &amp;
         &amp;                                 distributionOperatorSequence &amp;
         &amp;                                )
     </constructor>
    </referenceConstruct>
    !!]
    ! Create anti-log10 operator.
    allocate(outputAnalysisPropertyOperatorAntiLog10_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorAntiLog10_"         constructor="outputAnalysisPropertyOperatorAntiLog10          (                                                                                                                                                             )"/>
    !!]
    ! Create a filter to select galaxies in the required stellar mass and morphology range.
    allocate   (galacticFilterMassStellarMinimum_        )
    !![
    <referenceConstruct object="galacticFilterMassStellarMinimum_"                 constructor="galacticFilterStellarMass                       (massStellarMinimum*1.0d-1**(5.0d0*massStellarErrorDex)                                                                                                       )"/>
    !!]
    allocate   (galacticFilterMassStellarMaximum_        )
    !![
    <referenceConstruct object="galacticFilterMassStellarMaximum_"                 constructor="galacticFilterStellarMass                       (massStellarMaximum*1.0d+1**(5.0d0*massStellarErrorDex)                                                                                                       )"/>
    !!]
    allocate   (galacticFilterMassStellarMaximumInverted_)
    !![
    <referenceConstruct object="galacticFilterMassStellarMaximumInverted_"         constructor="galacticFilterNot                               (galacticFilterMassStellarMaximum_                                                                                                                            )"/>
    !!]
    allocate   (galacticFilterMorphology_                )
    !![
    <referenceConstruct object="galacticFilterMorphology_"                         constructor="galacticFilterStellarMassMorphology             (massStellarRatio                                                                                                                                             )"/>
    !!]
    if (isLateType) then
       allocate(galacticFilterMorphologyInverted_        )
       !![
       <referenceConstruct object="galacticFilterMorphologyInverted_"              constructor="galacticFilterNot                               (galacticFilterMorphology_                                                                                                                                    )"/>
       !!]
    else
       nullify (galacticFilterMorphologyInverted_        )
    end if
    allocate(filters_                                )
    filter_            => filters_
    filter_   %filter_ => galacticFilterMassStellarMinimum_
    allocate(filter_%next)
    filter_            => filter_%next
    filter_   %filter_ => galacticFilterMassStellarMaximumInverted_
    allocate(filter_%next)
    filter_            => filter_%next
    if (isLateType) then
       filter_%filter_ => galacticFilterMorphologyInverted_
    else
       filter_%filter_ => galacticFilterMorphology_
    end if
    allocate   (galacticFilterAll_                       )
    !![
    <referenceConstruct object="galacticFilterAll_"                                constructor="galacticFilterAll                               (filters_                                                                                                                                                  )"/>
    !!]
    ! Create a distribution normalizer which normalizes to bin width.
    allocate(outputAnalysisDistributionNormalizerBinWidth_ )
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerBinWidth_"  constructor="outputAnalysisDistributionNormalizerBinWidth ()"/>
    !!]
    allocate(outputAnalysisDistributionNormalizerUnitarity_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerUnitarity_" constructor="outputAnalysisDistributionNormalizerUnitarity()"/>
    !!]
    allocate(normalizerSequence     )
    allocate(normalizerSequence%next)
    normalizerSequence     %normalizer_ => outputAnalysisDistributionNormalizerUnitarity_
    normalizerSequence%next%normalizer_ => outputAnalysisDistributionNormalizerBinWidth_
    allocate(outputAnalysisDistributionNormalizer_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizer_" constructor="outputAnalysisDistributionNormalizerSequence(normalizerSequence)"/>
    !!]
    ! Construct the object. We convert radii to log10(radii) here.
    write (distributionName,'(i2.2)') distributionNumber
    description="Distribution of half-mass radii; "
    if (isLateType) then
       description=description//"late-type; "
    else
       description=description//"early-type; "
    end if
    write (massStellarMinimumLogarithmic,'(f5.2)') log10(massStellarMinimum)
    write (massStellarMaximumLogarithmic,'(f5.2)') log10(massStellarMaximum)
    description=description//"$"//trim(adjustl(massStellarMinimumLogarithmic))//" < \log_{10}(M_\star/M_\odot) < "//trim(adjustl(massStellarMaximumLogarithmic))//"$"
    self%outputAnalysisVolumeFunction1D=                                                          &
         & outputAnalysisVolumeFunction1D(                                                        &
         &                                var_str('galaxySizesSDSS')//trim(distributionName)    , &
         &                                description                                           , &
         &                                var_str('radius'                                     ), &
         &                                var_str('Radius at the bin center'                   ), &
         &                                var_str('kpc'                                        ), &
         &                                milli*megaParsec                                      , &
         &                                var_str('galaxySizesSDSSFunction'                    ), &
         &                                var_str('Galaxy size function averaged over each bin'), &
         &                                var_str('dimensionless'                              ), &
         &                                0.0d0                                                 , &
         &                                log10(radii)                                          , &
         &                                0_c_size_t                                            , &
         &                                outputWeight                                          , &
         &                                nodePropertyExtractor_                                , &
         &                                outputAnalysisPropertyOperatorSequence_               , &
         &                                outputAnalysisPropertyOperatorAntiLog10_              , &
         &                                outputAnalysisWeightOperator_                         , &
         &                                outputAnalysisDistributionOperator_                   , &
         &                                outputAnalysisDistributionNormalizer_                 , &
         &                                galacticFilterAll_                                    , &
         &                                outputTimes_                                          , &
         &                                outputAnalysisCovarianceModelPoisson                  , &
         &                                covarianceBinomialBinsPerDecade                       , &
         &                                covarianceBinomialMassHaloMinimum                     , &
         &                                covarianceBinomialMassHaloMaximum                     , &
         &                                .false.                                               , &
         &                                var_str('$r_{1/2}/\mathrm{kpc}$'                   )  , &
         &                                var_str('$\mathrm{d}p/\mathrm{d}\log_{10} r_{1/2}$')  , &
         &                                .true.                                                , &
         &                                .false.                                               , &
         &                                var_str('Shen et al. (2003)')                         , &
         &                                functionValueTarget                                   , &
         &                                functionCovarianceTarget                                &
         &                               )
    !![
    <objectDestructor name="surveyGeometry_"                                 />
    <objectDestructor name="nodePropertyExtractor_"                          />
    <objectDestructor name="outputAnalysisPropertyOperatorSequence_"         />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSequence_"   />
    <objectDestructor name="outputAnalysisPropertyOperatorMultiply_"         />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_" />
    <objectDestructor name="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_"/>
    <objectDestructor name="outputAnalysisPropertyOperatorAntiLog10_"        />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"           />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"          />
    <objectDestructor name="outputAnalysisWeightOperator_"                   />
    <objectDestructor name="outputAnalysisDistributionOperator_"             />
    <objectDestructor name="outputAnalysisDistributionOperatorGrvtnlLnsng_"  />
    <objectDestructor name="outputAnalysisDistributionOperatorProjection_"   />
    <objectDestructor name="galacticFilterAll_"                              />
    <objectDestructor name="galacticFilterMassStellarMinimum_"               />
    <objectDestructor name="galacticFilterMassStellarMaximum_"               />
    <objectDestructor name="galacticFilterMassStellarMaximumInverted_"       />
    <objectDestructor name="galacticFilterMorphology_"                       />
    <objectDestructor name="galacticFilterMorphologyInverted_"               />
    <objectDestructor name="cosmologyParametersData"                         />
    <objectDestructor name="cosmologyFunctionsData"                          />
    <objectDestructor name="outputAnalysisDistributionNormalizerBinWidth_"   />
    <objectDestructor name="outputAnalysisDistributionNormalizerUnitarity_"  />
    !!]
    nullify(propertyOperatorSequence      )
    nullify(weightPropertyOperatorSequence)
    nullify(normalizerSequence            )
    nullify(distributionOperatorSequence  )
    nullify(filters_                      )
    return
  end function galaxySizesSDSSConstructorInternal

  subroutine galaxySizesSDSSDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisGalaxySizesSDSS} output analysis class.
    !!}
    implicit none
    type(outputAnalysisGalaxySizesSDSS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%gravitationalLensing_"/>
    !!]
    return
  end subroutine galaxySizesSDSSDestructor
