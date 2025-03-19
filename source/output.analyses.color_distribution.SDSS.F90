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
Implements a color distribution output analysis class for SDSS data.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <outputAnalysis name="outputAnalysisColorDistributionSDSS">
   <description>An SDSS color distribution function output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisColorDistributionSDSS
     !!{
     An SDSS color distribution output analysis class.
     !!}
     private
     integer                                   :: distributionNumber
     class  (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
   contains
     final :: colorDistributionSDSSDestructor
  end type outputAnalysisColorDistributionSDSS

  interface outputAnalysisColorDistributionSDSS
     !!{
     Constructors for the ``colorDistributionSDSS'' output analysis class.
     !!}
     module procedure colorDistributionSDSSConstructorParameters
     module procedure colorDistributionSDSSConstructorInternal
  end interface outputAnalysisColorDistributionSDSS

contains

  function colorDistributionSDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``colorDistributionSDSS'' output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisColorDistributionSDSS)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass                   ), pointer       :: outputTimes_
    integer                                                              :: distributionNumber

    !![
    <inputParameter>
      <name>distributionNumber</name>
      <source>parameters</source>
      <description>The number (1-16) of the distribution to compute.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="outputTimes"        name="outputTimes_"        source="parameters"/>
    !!]
    self=outputAnalysisColorDistributionSDSS(distributionNumber,cosmologyFunctions_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="outputTimes_"       />
    !!]
    return
  end function colorDistributionSDSSConstructorParameters

  function colorDistributionSDSSConstructorInternal(distributionNumber,cosmologyFunctions_,outputTimes_) result(self)
    !!{
    Internal constructor for the ``colorDistributionSDSS'' output analysis class.
    !!}
    use :: Cosmology_Functions                     , only : cosmologyFunctionsClass                           , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                    , only : cosmologyParametersSimple
    use :: Galactic_Filters                        , only : galacticFilterStellarMass
    use :: Error                                   , only : Error_Report
    use :: Input_Paths                             , only : inputPath                                         , pathTypeDataStatic
    use :: Geometry_Surveys                        , only : surveyGeometryMonteroDorta2009SDSS
    use :: HDF5_Access                             , only : hdf5Access
    use :: IO_HDF5                                 , only : hdf5Object
    use :: ISO_Varying_String                      , only : var_str                                           , varying_string
    use :: Node_Property_Extractors                , only : nodePropertyExtractorLmnstyStllrCF2000            , nodePropertyExtractorRatio
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Normalizers, only : normalizerList                                    , outputAnalysisDistributionNormalizerBinWidth, outputAnalysisDistributionNormalizerSequence, outputAnalysisDistributionNormalizerUnitarity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorRandomErrorFixed
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc   , outputAnalysisPropertyOperatorIdentity      , outputAnalysisPropertyOperatorMagnitude     , outputAnalysisPropertyOperatorSequence       , &
          &                                                 propertyOperatorList
    use :: Output_Analysis_Utilities               , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorNormal
    use :: Output_Times                            , only : outputTimesClass
    implicit none
    type            (outputAnalysisColorDistributionSDSS               )                              :: self
    integer                                                                          , intent(in   )  :: distributionNumber
    class           (cosmologyFunctionsClass                           ), target     , intent(in   )  :: cosmologyFunctions_
    class           (outputTimesClass                                  ), target     , intent(inout)  :: outputTimes_
    type            (cosmologyParametersSimple                         ), pointer                     :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                    ), pointer                     :: cosmologyFunctionsData
    type            (nodePropertyExtractorRatio                        ), pointer                     :: nodePropertyExtractorRatio_
    type            (nodePropertyExtractorLmnstyStllrCF2000            ), pointer                     :: nodePropertyExtractorBandR_                     , nodePropertyExtractorBandU_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc   ), pointer                     :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorMagnitude           ), pointer                     :: outputAnalysisPropertyOperatorMagnitude_
    type            (outputAnalysisPropertyOperatorSequence            ), pointer                     :: outputAnalysisWeightPropertyOperatorSequence_
    type            (outputAnalysisPropertyOperatorIdentity            ), pointer                     :: outputAnalysisPropertyOperatorIdentity_
    type            (outputAnalysisWeightOperatorNormal                ), pointer                     :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerSequence      ), pointer                     :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionOperatorRandomErrorFixed), pointer                     :: outputAnalysisDistributionOperator_
    type            (surveyGeometryMonteroDorta2009SDSS                )               , pointer      :: surveyGeometry_
    type            (normalizerList                                    ), pointer                     :: normalizerSequence                               , normalizer_
    type            (propertyOperatorList                              ), pointer                     :: propertyOperatorSequence                         , weightPropertyOperatorSequence
    type            (galacticFilterStellarMass                         ), pointer                     :: galacticFilter_
    double precision                                                    , allocatable, dimension(:  ) :: colors                                           , functionValueTarget                                , &
         &                                                                                               functionErrorTarget
    double precision                                                    , allocatable, dimension(:,:) :: outputWeight                                    , functionCovarianceTarget
    integer                                                             , parameter                   :: covarianceBinomialBinsPerDecade         =2
    double precision                                                    , parameter                   :: covarianceBinomialMassHaloMinimum       =3.00d+11, covarianceBinomialMassHaloMaximum    =1.0d15
    double precision                                                    , parameter                   :: redshiftBand                            =1.00d-01
    double precision                                                    , parameter                   :: magnitudeErrorR                         =1.50d-02
    double precision                                                    , parameter                   :: magnitudeErrorU                         =3.50d-02
    double precision                                                    , parameter                   :: massStellarMinimum                      =1.00d+06
    integer         (c_size_t                                          )                              :: iBin                                             , bufferCount
    type            (hdf5Object                                        )                              :: dataFile                                         , distribution
    double precision                                                                                  :: magnitudeMinimum                                 , magnitudeMaximum
    character       (len=16                                            )                              :: distributionName                                 , magnitudeMinimumLabel                               , &
         &                                                                                               magnitudeMaximumLabel
    type            (varying_string                                    )                              :: description
    !![
    <constructorAssign variables="distributionNumber, *cosmologyFunctions_"/>
    !!]

    ! Validate input.
    if (distributionNumber < 1 .or. distributionNumber > 16) call Error_Report('distributionNumber ∈ [1..16] is required'//{introspection:location})
    ! Construct colors matched to those used by Baldry et al. (2004). Also read magnitude range.
    write (distributionName,'(a,i2.2)') 'distribution',distributionNumber
    !$ call hdf5Access%set()
    call dataFile    %openFile     (char(inputPath(pathTypeDataStatic)//'observations/galaxyColors/colorDistributionsBaldry2004.hdf5'),readOnly=.true.             )
    distribution=dataFile%openGroup(distributionName)
    call distribution%readDataset  (                                         'color'                                                  ,         colors             )
    call distribution%readDataset  (                                         'distribution'                                           ,         functionValueTarget)
    call distribution%readDataset  (                                         'distributionError'                                      ,         functionErrorTarget)
    call distribution%readAttribute(                                         'magnitudeMinimum'                                       ,         magnitudeMinimum   )
    call distribution%readAttribute(                                         'magnitudeMaximum'                                       ,         magnitudeMaximum   )
    call distribution%close        (                                                                                                                               )
    call dataFile    %close        (                                                                                                                               )
    !$ call hdf5Access%unset()
    self %binCount=size(colors)
   allocate(functionCovarianceTarget(self%binCount,self%binCount))
    functionCovarianceTarget=0.0d0
    do iBin=1,self%binCount
       ! Negative values indicate upper limits, which mean no galaxies were observed in the given bin.
       functionValueTarget     (iBin     )=max(0.0d0,functionValueTarget(iBin))
       functionCovarianceTarget(iBin,iBin)=          functionErrorTarget(iBin) **2
    end do
     ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
     <constructor>
      cosmologyParametersSimple     (                            &amp;
        &amp;                        OmegaMatter    = 0.30000d0, &amp;
        &amp;                        OmegaDarkEnergy= 0.70000d0, &amp;
        &amp;                        HubbleConstant =70.00000d0, &amp;
        &amp;                        temperatureCMB = 2.72548d0, &amp;
        &amp;                        OmegaBaryon    = 0.04550d0  &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
     <constructor>
      cosmologyFunctionsMatterLambda(                            &amp;
        &amp;                        cosmologyParametersData     &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the SDSS survey geometry of Baldry et al. (2004) with their imposed redshift limits.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_" constructor="surveyGeometryMonteroDorta2009SDSS(band='r',redshiftMinimum=4.0d-3,redshiftMaximum=8.0d-2,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(self%binCount,outputTimes_%count()))
    do iBin=1,self%binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,self%cosmologyFunctions_,outputTimes_,magnitudeAbsoluteLimit=magnitudeMaximum)
    end do
    ! Create stellar luminosity property extractors.
    allocate(nodePropertyExtractorBandR_           )
    allocate(nodePropertyExtractorBandU_           )
    !![
    <referenceConstruct object="nodePropertyExtractorBandR_"            constructor="nodePropertyExtractorLmnstyStllrCF2000  ('SDSS_r','observed',depthOpticalISMCoefficient=1.0d0,depthOpticalCloudsCoefficient=1.0d0,wavelengthExponent=0.7d0,outputTimes_=outputTimes_,redshiftBand=redshiftBand,outputMask=sum(outputWeight,dim=1) > 0.0d0)"/>
    <referenceConstruct object="nodePropertyExtractorBandU_"            constructor="nodePropertyExtractorLmnstyStllrCF2000  ('SDSS_u','observed',depthOpticalISMCoefficient=1.0d0,depthOpticalCloudsCoefficient=1.0d0,wavelengthExponent=0.7d0,outputTimes_=outputTimes_,redshiftBand=redshiftBand,outputMask=sum(outputWeight,dim=1) > 0.0d0)"/>
    !!]
    ! Create a ratio property extractor.
    allocate(nodePropertyExtractorRatio_        )
    !![
    <referenceConstruct object="nodePropertyExtractorRatio_"            constructor="nodePropertyExtractorRatio              ('color','SDSS u-r color',nodePropertyExtractorBandU_,nodePropertyExtractorBandR_                                 )"/>
    !!]
    ! Create magnitude, and cosmological luminosity distance property operators.
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    allocate(outputAnalysisPropertyOperatorMagnitude_        )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc   (cosmologyFunctions_,cosmologyFunctionsData,outputTimes_                                                                              )"/>
    <referenceConstruct object="outputAnalysisPropertyOperatorMagnitude_"         constructor="outputAnalysisPropertyOperatorMagnitude()"/>
    !!]
    allocate(weightPropertyOperatorSequence                  )
    allocate(weightPropertyOperatorSequence%next             )
    weightPropertyOperatorSequence     %operator_ => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    weightPropertyOperatorSequence%next%operator_ => outputAnalysisPropertyOperatorMagnitude_
    allocate(outputAnalysisWeightPropertyOperatorSequence_   )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperatorSequence_"    constructor="outputAnalysisPropertyOperatorSequence            (weightPropertyOperatorSequence                                                                                                       )"/>
    !!]
    ! Create a normal-weight weight operator.
    allocate(outputAnalysisWeightOperator_                   )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                    constructor="outputAnalysisWeightOperatorNormal                (magnitudeMinimum,magnitudeMaximum,magnitudeErrorR,nodePropertyExtractorBandR_,outputAnalysisWeightPropertyOperatorSequence_)"/>
    !!]
    ! Create random error distribution operator.
    allocate(outputAnalysisDistributionOperator_             )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_"              constructor="outputAnalysisDistributionOperatorRandomErrorFixed(sqrt(magnitudeErrorR**2+magnitudeErrorU**2))"/>
    !!]
    ! Create identity operator.
    allocate(outputAnalysisPropertyOperatorIdentity_         )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorIdentity_"          constructor="outputAnalysisPropertyOperatorIdentity            (                                                                                                                                     )"/>
    !!]
    ! Create a filter to select galaxies above some minimum stellar mass.
    allocate(galacticFilter_                                 )
    !![
    <referenceConstruct object="galacticFilter_"                                  constructor="galacticFilterStellarMass                         (massStellarMinimum                                                                                                                   )"/>
    !!]
    ! Create a distribution normalizer which normalizes to bin width and unitarity.
    allocate(normalizerSequence)
    normalizer_ => normalizerSequence
    allocate(outputAnalysisDistributionNormalizerUnitarity :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerUnitarity)
       !![
       <referenceConstruct object="normalizer_" constructor="outputAnalysisDistributionNormalizerUnitarity ()"/>
       !!]
    end select
    allocate(normalizer_%next)
    normalizer_ => normalizer_%next
    allocate(outputAnalysisDistributionNormalizerBinWidth  :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerBinWidth )
       !![
       <referenceConstruct object="normalizer_" constructor="outputAnalysisDistributionNormalizerBinWidth  ()"/>
       !!]
    end select
    allocate(outputAnalysisDistributionNormalizer_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizer_"            constructor="outputAnalysisDistributionNormalizerSequence      (normalizerSequence                                                                                                                   )"/>
    !!]
    ! Determine number of buffer bins.
    bufferCount=int(3.0d0*(colors(2)-colors(1))/sqrt(magnitudeErrorR**2+magnitudeErrorU**2),kind=c_size_t)+1_c_size_t
    ! Construct the object.
    write (distributionName,'(i2.2)') distributionNumber
    description="Distribution of SDSS u-r color; "
    write (magnitudeMinimumLabel,'(f6.2)') magnitudeMinimum
    write (magnitudeMaximumLabel,'(f6.2)') magnitudeMaximum
    description=description//"$"//trim(adjustl(magnitudeMinimumLabel))//" < r < "//trim(adjustl(magnitudeMaximumLabel))//"$"
    self%outputAnalysisVolumeFunction1D=                                                            &
         & outputAnalysisVolumeFunction1D(                                                          &
         &                                var_str('colorDistributionSDSS')//trim(distributionName), &
         &                                description                                             , &
         &                                var_str('color'                                        ), &
         &                                var_str('Color at the bin center'                      ), &
         &                                var_str('dimensionless'                                ), &
         &                                0.0d0                                                   , &
         &                                var_str('colorDistributionSDSSFunction'                ), &
         &                                var_str('Color distribution averaged over each bin'    ), &
         &                                var_str('dimensionless'                                ), &
         &                                0.0d0                                                   , &
         &                                colors                                                  , &
         &                                bufferCount                                             , &
         &                                outputWeight                                            , &
         &                                nodePropertyExtractorRatio_                   , &
         &                                outputAnalysisPropertyOperatorMagnitude_                , &
         &                                outputAnalysisPropertyOperatorIdentity_                 , &
         &                                outputAnalysisWeightOperator_                           , &
         &                                outputAnalysisDistributionOperator_                     , &
         &                                outputAnalysisDistributionNormalizer_                   , &
         &                                galacticFilter_                                         , &
         &                                outputTimes_                                            , &
         &                                outputAnalysisCovarianceModelPoisson                    , &
         &                                covarianceBinomialBinsPerDecade                         , &
         &                                covarianceBinomialMassHaloMinimum                       , &
         &                                covarianceBinomialMassHaloMaximum                       , &
         &                                .false.                                                 , &
         &                                var_str('$u-r$'                        )                , &
         &                                var_str('$\mathrm{d}p/\mathrm{d}(u-r)$')                , &
         &                                .false.                                                 , &
         &                                .false.                                                 , &
         &                                var_str('Baldry et al. (2004)')                         , &
         &                                functionValueTarget                                     , &
         &                                functionCovarianceTarget                                  &
         &                               )
    ! Clean up.
    !![
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSequence_"   />
    <objectDestructor name="outputAnalysisPropertyOperatorMagnitude_"        />
    <objectDestructor name="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_"/>
    <objectDestructor name="outputAnalysisPropertyOperatorIdentity_"         />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"           />
    <objectDestructor name="outputAnalysisDistributionOperator_"             />
    <objectDestructor name="nodePropertyExtractorBandR_"                     />
    <objectDestructor name="nodePropertyExtractorBandU_"                     />
    <objectDestructor name="nodePropertyExtractorRatio_"                     />
    <objectDestructor name="outputAnalysisWeightOperator_"                   />
    <objectDestructor name="galacticFilter_"                                 />
    <objectDestructor name="cosmologyParametersData"                         />
    <objectDestructor name="cosmologyFunctionsData"                          />
    !!]
    nullify(propertyOperatorSequence      )
    nullify(weightPropertyOperatorSequence)
    nullify(normalizerSequence            )
    return
  end function colorDistributionSDSSConstructorInternal

  subroutine colorDistributionSDSSDestructor(self)
    !!{
    Destructor for the ``colorDistributionSDSS'' output analysis class.
    !!}
    implicit none
    type(outputAnalysisColorDistributionSDSS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine colorDistributionSDSSDestructor
