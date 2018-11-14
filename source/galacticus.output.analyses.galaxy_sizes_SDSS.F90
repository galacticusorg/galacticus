!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a galaxy size output analysis class for SDSS data.
  
  use Cosmology_Functions

  !# <outputAnalysis name="outputAnalysisGalaxySizesSDSS" defaultThreadPrivate="yes">
  !#  <description>A stellar mass function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisGalaxySizesSDSS
     !% A galaxySizesSDSS output analysis class.
     private
     integer                                            :: distributionNumber
     double precision                                   :: massStellarRatio
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_
   contains
     final :: galaxySizesSDSSDestructor
  end type outputAnalysisGalaxySizesSDSS

  interface outputAnalysisGalaxySizesSDSS
     !% Constructors for the ``galaxySizesSDSS'' output analysis class.
     module procedure galaxySizesSDSSConstructorParameters
     module procedure galaxySizesSDSSConstructorInternal
  end interface outputAnalysisGalaxySizesSDSS

contains

  function galaxySizesSDSSConstructorParameters(parameters) result (self)
    !% Constructor for the ``galaxySizesSDSS'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisGalaxySizesSDSS)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass             ), pointer       :: outputTimes_
    double precision                                               :: massStellarRatio
    integer                                                        :: distributionNumber

    !# <inputParameter>
    !#   <name>distributionNumber</name>
    !#   <source>parameters</source>
    !#   <description>The number (1-34) of the distribution to compute.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massStellarRatio</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.3d0</defaultValue>
    !#   <description>The stellar mass bulge-to-total ratio used to discriminate late-type vs. early-type galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !# <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    self=outputAnalysisGalaxySizesSDSS(distributionNumber,massStellarRatio,cosmologyFunctions_,outputTimes_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function galaxySizesSDSSConstructorParameters

  function galaxySizesSDSSConstructorInternal(distributionNumber,massStellarRatio,cosmologyFunctions_,outputTimes_) result(self)
    !% Internal constructor for the ``galaxySizesSDSS'' output analysis class.
    use ISO_Varying_String
    use Output_Times
    use Output_Analyses_Options
    use Output_Analysis_Utilities
    use Galacticus_Error
    use Galacticus_Paths
    use IO_HDF5
    use Memory_Management
    use Numerical_Comparison
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    type            (outputAnalysisGalaxySizesSDSS                  )                              :: self
    integer                                                                       , intent(in   )  :: distributionNumber
    double precision                                                              , intent(in   )  :: massStellarRatio
    class           (cosmologyFunctionsClass                        ), target     , intent(in   )  :: cosmologyFunctions_
    class           (outputTimesClass                               ), target     , intent(inout)  :: outputTimes_
    type            (cosmologyParametersSimple                      ), pointer                     :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                 ), pointer                     :: cosmologyFunctionsData
    type            (outputAnalysisPropertyExtractorHalfMassRadius  ), pointer                     :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyExtractorMassStellar     ), pointer                     :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc ), pointer                     :: outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc), pointer                     :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorLog10            ), pointer                     :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorMultiply         ), pointer                     :: outputAnalysisPropertyOperatorMultiply_
    type            (outputAnalysisPropertyOperatorSequence         ), pointer                     :: outputAnalysisPropertyOperatorSequence_         , outputAnalysisWeightPropertyOperatorSequence_
    type            (outputAnalysisWeightOperatorNormal             ), pointer                     :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerSequence   ), pointer                     :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisPropertyOperatorAntiLog10        ), pointer                     :: outputAnalysisPropertyOperatorAntiLog10_
    class           (outputAnalysisDistributionOperatorClass        ), pointer                     :: outputAnalysisDistributionOperator_
    type            (surveyGeometryLiWhite2009SDSS                  )               , pointer      :: surveyGeometry_
    type            (normalizerList                                 ), pointer                     :: normalizerSequence                              , normalizer_
    type            (propertyOperatorList                           ), pointer                     :: propertyOperatorSequence                        , weightPropertyOperatorSequence
    type            (galacticFilterStellarMass                      ), pointer                     :: galacticFilterMassStellarMinimum_               , galacticFilterMassStellarMaximum_
    type            (galacticFilterStellarMassMorphology            ), pointer                     :: galacticFilterMorphology_
    type            (galacticFilterNot                              ), pointer                     :: galacticFilterMassStellarMaximumInverted_       , galacticFilterMorphologyInverted_
    type            (galacticFilterAll                              ), pointer                     :: galacticFilterAll_
    type            (filterList                                     ), pointer                     :: filters_                                        , filter_
    double precision                                                 , allocatable, dimension(:  ) :: radii
    double precision                                                 , allocatable, dimension(:,:) :: outputWeight
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
    character       (len=16                                         )                              :: distributionName
    logical                                                                                        :: isLateType
    !# <constructorAssign variables="distributionNumber, massStellarRatio, *cosmologyFunctions_"/>

    ! Validate input.
    if (distributionNumber < 1 .or. distributionNumber > 34) call Galacticus_Error_Report('distributionNumber ∈ [1..34] is required'//{introspection:location})
    ! Construct sizes matched to those used by  Shen et al. (2003). Also read stellar mass and Sersic index ranges.
    write (distributionName,'(a,i2.2)') 'distribution',distributionNumber
    !$ call hdf5Access%set()
    call dataFile    %openFile     (char(galacticusPath(pathTypeDataStatic)//'observations/galaxySizes/Galaxy_Sizes_By_Mass_SDSS_Shen_2003.hdf5'),readOnly=.true.            )
    distribution=dataFile%openGroup(distributionName)
    call distribution%readDataset  (                                         'radius'                                                            ,         radii             )
    call distribution%readAttribute(                                         'massMinimum'                                                       ,         massStellarMinimum)
    call distribution%readAttribute(                                         'massMaximum'                                                       ,         massStellarMaximum)
    call distribution%readAttribute(                                         'sersicIndexMinimum'                                                ,         indexSersicMinimum)
    call distribution%readAttribute(                                         'sersicIndexMaximum'                                                ,         indexSersicMaximum)
    call distribution%close        (                                                                                                                                         )
    call dataFile    %close        (                                                                                                                                         )
    !$ call hdf5Access%unset()
    self %binCount=size(radii)
    ! Determine if this distribution is for late-type galaxies.
    isLateType=(indexSersicMinimum == 0.0d0)
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    cosmologyParametersData=cosmologyParametersSimple     (                            &
         &                                                 OmegaMatter    = 0.30000d0, &
         &                                                 OmegaDarkEnergy= 0.70000d0, &
         &                                                 HubbleConstant =70.00000d0, &
         &                                                 temperatureCMB = 2.72548d0, &
         &                                                 OmegaBaryon    = 0.04550d0  &
         &                                                )
    cosmologyFunctionsData =cosmologyFunctionsMatterLambda(                            &
         &                                                 cosmologyParametersData     &
         &                                                )
    ! Build the SDSS survey geometry of Shen et al. (2013) with their imposed redshift limits.
    allocate(surveyGeometry_)
    surveyGeometry_=surveyGeometryLiWhite2009SDSS(redshiftMinimum=1.0d-3,redshiftMaximum=huge(0.0d0),cosmologyFunctions_=cosmologyFunctions_)
    ! Compute weights that apply to each output redshift.
    call allocateArray(outputWeight,[self%binCount,outputTimes_%count()])
    do iBin=1,self%binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,self%cosmologyFunctions_,outputTimes_,massStellarMinimum)
    end do
    ! Create a half-mass radius property extractor.
    allocate(outputAnalysisPropertyExtractor_        )
    outputAnalysisPropertyExtractor_                =outputAnalysisPropertyExtractorHalfMassRadius     (                                                                                                                                                            )
    ! Create a stellar mass property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_        )
    outputAnalysisWeightPropertyExtractor_          =outputAnalysisPropertyExtractorMassStellar        (                                                                                                                                                            )
    ! Create multiply, log10, cosmological angular distance, and cosmologyical luminosity distance property operators.
    allocate(outputAnalysisPropertyOperatorMultiply_         )
    outputAnalysisPropertyOperatorMultiply_         =outputAnalysisPropertyOperatorMultiply            (kilo                                                                                                                                                        )
    allocate(outputAnalysisPropertyOperatorLog10_            )
    outputAnalysisPropertyOperatorLog10_            =outputAnalysisPropertyOperatorLog10               (                                                                                                                                                            )
    allocate(outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_)
    outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_ =outputAnalysisPropertyOperatorCsmlgyAnglrDstnc    (cosmologyFunctions_             ,cosmologyFunctionsData,outputTimes_                                                                                        )
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_=outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc   (cosmologyFunctions_             ,cosmologyFunctionsData,outputTimes_                                                                                        )
    allocate(propertyOperatorSequence          )
    allocate(propertyOperatorSequence%next     )
    allocate(propertyOperatorSequence%next%next)
    propertyOperatorSequence          %operator_ => outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_
    propertyOperatorSequence%next     %operator_ => outputAnalysisPropertyOperatorMultiply_
    propertyOperatorSequence%next%next%operator_ => outputAnalysisPropertyOperatorLog10_
    allocate(outputAnalysisPropertyOperatorSequence_ )
    outputAnalysisPropertyOperatorSequence_         =outputAnalysisPropertyOperatorSequence            (propertyOperatorSequence                                                                                                                                   )
    allocate(weightPropertyOperatorSequence          )
    allocate(weightPropertyOperatorSequence%next     )
    weightPropertyOperatorSequence          %operator_ => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    weightPropertyOperatorSequence%next%operator_ => outputAnalysisPropertyOperatorLog10_
    allocate(outputAnalysisWeightPropertyOperatorSequence_ )
    outputAnalysisWeightPropertyOperatorSequence_   =outputAnalysisPropertyOperatorSequence           (weightPropertyOperatorSequence                                                                                                                              )
    ! Create a normal-weight weight operator.
    allocate(outputAnalysisWeightOperator_           )
    outputAnalysisWeightOperator_                   =outputAnalysisWeightOperatorNormal               (log10(massStellarMinimum),log10(massStellarMaximum),massStellarErrorDex,outputAnalysisWeightPropertyExtractor_,outputAnalysisWeightPropertyOperatorSequence_)
    ! Create an inclination distribution operator.
    if (isLateType) then
       allocate(outputAnalysisDistributionOperatorDiskSizeInclntn :: outputAnalysisDistributionOperator_)
       select type (outputAnalysisDistributionOperator_)
       type is (outputAnalysisDistributionOperatorDiskSizeInclntn)
          outputAnalysisDistributionOperator_          =outputAnalysisDistributionOperatorDiskSizeInclntn(                                                                                                                                                             )
       end select
    else
       allocate(outputAnalysisDistributionOperatorIdentity        :: outputAnalysisDistributionOperator_)
       select type (outputAnalysisDistributionOperator_)
       type is (outputAnalysisDistributionOperatorIdentity)
          outputAnalysisDistributionOperator_          =outputAnalysisDistributionOperatorIdentity       (                                                                                                                                                             )
       end select
    end if
    ! Create anti-log10 operator.
    allocate(outputAnalysisPropertyOperatorAntiLog10_)
    outputAnalysisPropertyOperatorAntiLog10_        =outputAnalysisPropertyOperatorAntiLog10          (                                                                                                                                                             )
    ! Create a filter to select galaxies in the required stellar mass and morphology range.
    allocate   (galacticFilterMassStellarMinimum_        )
    galacticFilterMassStellarMinimum_                =galacticFilterStellarMass                       (massStellarMinimum*1.0d-1**(5.0d0*massStellarErrorDex)                                                                                                       )
    allocate   (galacticFilterMassStellarMaximum_        )
    galacticFilterMassStellarMaximum_                =galacticFilterStellarMass                       (massStellarMaximum*1.0d+1**(5.0d0*massStellarErrorDex)                                                                                                       )
    allocate   (galacticFilterMassStellarMaximumInverted_)
    galacticFilterMassStellarMaximumInverted_        =galacticFilterNot                               (galacticFilterMassStellarMaximum_                                                                                                                            )
    allocate   (galacticFilterMorphology_                )
    galacticFilterMorphology_                        =galacticFilterStellarMassMorphology             (massStellarRatio                                                                                                                                             )
    if (isLateType) then
       allocate(galacticFilterMorphologyInverted_        )
       galacticFilterMorphologyInverted_             =galacticFilterNot                               (galacticFilterMorphology_                                                                                                                                    )
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
    galacticFilterAll_                               =galacticFilterAll                               (filters_                                                                                                                                                  )
    ! Create a distribution normalizer which normalizes to bin width.
    allocate(normalizerSequence)
    normalizer_ => normalizerSequence
    allocate(outputAnalysisDistributionNormalizerUnitarity  :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerUnitarity  )
       normalizer_=outputAnalysisDistributionNormalizerUnitarity ()
    end select
    allocate(normalizer_%next)
    normalizer_ => normalizer_%next
    allocate(outputAnalysisDistributionNormalizerBinWidth   :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerBinWidth  )
       normalizer_=outputAnalysisDistributionNormalizerBinWidth  ()
    end select
    allocate(outputAnalysisDistributionNormalizer_)
    outputAnalysisDistributionNormalizer_=outputAnalysisDistributionNormalizerSequence(normalizerSequence)
    ! Construct the object. We convert radii to log10(radii) here.
    write (distributionName,'(i2.2)') distributionNumber
    self%outputAnalysisVolumeFunction1D=                                                          &
         & outputAnalysisVolumeFunction1D(                                                        &
         &                                var_str('galaxySizesSDSS')//trim(distributionName)    , &
         &                                var_str('Distribution of galaxy half-mass radii'     ), &
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
         &                                outputAnalysisPropertyExtractor_                      , &
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
         &                                covarianceBinomialMassHaloMaximum                       &
         &                               )
    ! Clean up.
    nullify(outputAnalysisPropertyExtractor_                )
    nullify(outputAnalysisPropertyOperatorSequence_         )
    nullify(outputAnalysisWeightPropertyOperatorSequence_   )
    nullify(outputAnalysisPropertyOperatorLog10_            )
    nullify(outputAnalysisPropertyOperatorCsmlgyAnglrDstnc_ )
    nullify(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    nullify(outputAnalysisPropertyOperatorAntiLog10_        )
    nullify(outputAnalysisDistributionNormalizer_           )
    nullify(outputAnalysisWeightPropertyExtractor_          )
    nullify(outputAnalysisWeightOperator_                   )
    nullify(propertyOperatorSequence                        )
    nullify(weightPropertyOperatorSequence                  )
    nullify(normalizerSequence                              )
    nullify(filters_                                        )
    nullify(galacticFilterAll_                              )
    nullify(galacticFilterMassStellarMinimum_               )
    nullify(galacticFilterMassStellarMaximum_               )
    nullify(galacticFilterMassStellarMaximumInverted_       )
    nullify(galacticFilterMorphology_                       )
    nullify(galacticFilterMorphologyInverted_               )
    nullify(cosmologyParametersData                         )
    nullify(cosmologyFunctionsData                          )
    return
  end function galaxySizesSDSSConstructorInternal

  subroutine galaxySizesSDSSDestructor(self)
    !% Destructor for the ``galaxySizesSDSS'' output analysis class.
    implicit none
    type(outputAnalysisGalaxySizesSDSS), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologyFunctions_" />    
    return
  end subroutine galaxySizesSDSSDestructor
