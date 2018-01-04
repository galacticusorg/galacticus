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

!% Contains a module which implements a color distribution output analysis class for SDSS data.
  
  use Cosmology_Functions

  !# <outputAnalysis name="outputAnalysisColorDistributionSDSS" defaultThreadPrivate="yes">
  !#  <description>An SDSS color distribution function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisColorDistributionSDSS
     !% An SDSS color distribution output analysis class.
     private
     integer                                   :: distributionNumber
     class  (cosmologyFunctionsClass), pointer :: cosmologyFunctions_
   contains
     final :: colorDistributionSDSSDestructor
  end type outputAnalysisColorDistributionSDSS

  interface outputAnalysisColorDistributionSDSS
     !% Constructors for the ``colorDistributionSDSS'' output analysis class.
     module procedure colorDistributionSDSSConstructorParameters
     module procedure colorDistributionSDSSConstructorInternal
  end interface outputAnalysisColorDistributionSDSS

contains

  function colorDistributionSDSSConstructorParameters(parameters) result (self)
    !% Constructor for the ``colorDistributionSDSS'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisColorDistributionSDSS)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    integer                                                              :: distributionNumber

    !# <inputParameter>
    !#   <name>distributionNumber</name>
    !#   <source>parameters</source>
    !#   <description>The number (1-16) of the distribution to compute.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    self=outputAnalysisColorDistributionSDSS(distributionNumber,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function colorDistributionSDSSConstructorParameters

  function colorDistributionSDSSConstructorInternal(distributionNumber,cosmologyFunctions_) result(self)
    !% Internal constructor for the ``colorDistributionSDSS'' output analysis class.
    use ISO_Varying_String
    use Galacticus_Output_Times
    use Output_Analyses_Options
    use Output_Analysis_Utilities
    use Galacticus_Error
    use Galacticus_Input_Paths
    use IO_HDF5
    use Memory_Management
    use Numerical_Comparison
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    type            (outputAnalysisColorDistributionSDSS               )                              :: self
    integer                                                                          , intent(in   )  :: distributionNumber
    class           (cosmologyFunctionsClass                           ), target     , intent(in   )  :: cosmologyFunctions_
    type            (cosmologyParametersSimple                         ), pointer                     :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                    ), pointer                     :: cosmologyFunctionsData
    type            (outputAnalysisPropertyExtractorRatio              ), pointer                     :: outputAnalysisPropertyExtractorRatio_
    type            (outputAnalysisPropertyExtractorLmnstyStllrCF2000  ), pointer                     :: outputAnalysisPropertyExtractorBandR_            , outputAnalysisPropertyExtractorBandU_
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
    double precision                                                    , allocatable, dimension(:  ) :: colors
    double precision                                                    , allocatable, dimension(:,:) :: outputWeight
    integer                                                             , parameter                   :: covarianceBinomialBinsPerDecade         =2
    double precision                                                    , parameter                   :: covarianceBinomialMassHaloMinimum       =3.00d+11, covarianceBinomialMassHaloMaximum    =1.0d15
    double precision                                                    , parameter                   :: redshiftBand                            =1.00d-01
    double precision                                                    , parameter                   :: magnitudeErrorR                         =1.50d-02
    double precision                                                    , parameter                   :: magnitudeErrorU                         =3.50d-02
    double precision                                                    , parameter                   :: massStellarMinimum                      =1.00d+06
    integer         (c_size_t                                          )                              :: iBin                                             , bufferCount
    type            (hdf5Object                                        )                              :: dataFile                                         , distribution
    double precision                                                                                  :: magnitudeMinimum                                 , magnitudeMaximum
    character       (len=16                                            )                              :: distributionName
    !# <constructorAssign variables="distributionNumber, *cosmologyFunctions_"/>

    ! Validate input.
    if (distributionNumber < 1 .or. distributionNumber > 16) call Galacticus_Error_Report('distributionNumber ∈ [1..16] is required'//{introspection:location})
    ! Construct colors matched to those used by Baldry et al. (2004). Also read magnitude range.
    write (distributionName,'(a,i2.2)') 'distribution',distributionNumber
    !$omp critical(HDF5_Access)
    call dataFile    %openFile     (char(Galacticus_Input_Path()//'data/observations/galaxyColors/colorDistributionsBaldry2004.hdf5'),readOnly=.true.          )
    distribution=dataFile%openGroup(distributionName)
    call distribution%readDataset  (                              'color'                                                            ,         colors          )
    call distribution%readAttribute(                              'magnitudeMinimum'                                                 ,         magnitudeMinimum)
    call distribution%readAttribute(                              'magnitudeMaximum'                                                 ,         magnitudeMaximum)
    call distribution%close        (                                                                                                                           )
    call dataFile    %close        (                                                                                                                           )
    !$omp end critical(HDF5_Access)
    self %binCount=size(colors)
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
    ! Build the SDSS survey geometry of Baldry et al. (2004) with their imposed redshift limits.
    allocate(surveyGeometry_)
    surveyGeometry_=surveyGeometryMonteroDorta2009SDSS(band='r',redshiftMinimum=4.0d-3,redshiftMaximum=8.0d-2,cosmologyFunctions_=cosmologyFunctions_)
    ! Compute weights that apply to each output redshift.
    call allocateArray(outputWeight,[self%binCount,Galacticus_Output_Time_Count()])
    do iBin=1,self%binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,self%cosmologyFunctions_,magnitudeAbsoluteLimit=magnitudeMaximum)
    end do
    ! Create stellar luminosity property extractors.
    allocate(outputAnalysisPropertyExtractorBandR_           )
    allocate(outputAnalysisPropertyExtractorBandU_           )
    outputAnalysisPropertyExtractorBandR_           =outputAnalysisPropertyExtractorLmnstyStllrCF2000  ('SDSS_r','observed',depthOpticalISMCoefficient=1.0d0,depthOpticalCloudsCoefficient=1.0d0,wavelengthExponent=0.7d0,redshiftBand=redshiftBand,outputMask=sum(outputWeight,dim=1) > 0.0d0)
    outputAnalysisPropertyExtractorBandU_           =outputAnalysisPropertyExtractorLmnstyStllrCF2000  ('SDSS_u','observed',depthOpticalISMCoefficient=1.0d0,depthOpticalCloudsCoefficient=1.0d0,wavelengthExponent=0.7d0,redshiftBand=redshiftBand,outputMask=sum(outputWeight,dim=1) > 0.0d0)
    ! Create a ratio property extractor.
    allocate(outputAnalysisPropertyExtractorRatio_        )
    outputAnalysisPropertyExtractorRatio_           =outputAnalysisPropertyExtractorRatio              (outputAnalysisPropertyExtractorBandU_,outputAnalysisPropertyExtractorBandR_                                                          )
    ! Creat magnitude, and cosmological luminosity distance property operators.
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    allocate(outputAnalysisPropertyOperatorMagnitude_        )
    outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_=outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc   (cosmologyFunctions_,cosmologyFunctionsData                                                                                           )
    outputAnalysisPropertyOperatorMagnitude_        =outputAnalysisPropertyOperatorMagnitude()
    allocate(weightPropertyOperatorSequence                  )
    allocate(weightPropertyOperatorSequence%next             )
    weightPropertyOperatorSequence     %operator_ => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    weightPropertyOperatorSequence%next%operator_ => outputAnalysisPropertyOperatorMagnitude_
    allocate(outputAnalysisWeightPropertyOperatorSequence_   )
    outputAnalysisWeightPropertyOperatorSequence_   =outputAnalysisPropertyOperatorSequence            (weightPropertyOperatorSequence                                                                                                       )
    ! Create a normal-weight weight operator.
    allocate(outputAnalysisWeightOperator_                   )
    outputAnalysisWeightOperator_                   =outputAnalysisWeightOperatorNormal                (magnitudeMinimum,magnitudeMaximum,magnitudeErrorR,outputAnalysisPropertyExtractorBandR_,outputAnalysisWeightPropertyOperatorSequence_)
    ! Create random error distribution operator.
    allocate(outputAnalysisDistributionOperator_             )
    outputAnalysisDistributionOperator_             =outputAnalysisDistributionOperatorRandomErrorFixed(sqrt(magnitudeErrorR**2+magnitudeErrorU**2))
    ! Create identity operator.
    allocate(outputAnalysisPropertyOperatorIdentity_         )
    outputAnalysisPropertyOperatorIdentity_         =outputAnalysisPropertyOperatorIdentity            (                                                                                                                                     )
    ! Create a filter to select galaxies above some minimum stellar mass.
    allocate(galacticFilter_                                 )
    galacticFilter_                                 =galacticFilterStellarMass                         (massStellarMinimum                                                                                                                   )
    ! Create a distribution normalizer which normalizes to bin width and unitarity.
    allocate(normalizerSequence)
    normalizer_ => normalizerSequence
    allocate(outputAnalysisDistributionNormalizerUnitarity :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerUnitarity)
       normalizer_=outputAnalysisDistributionNormalizerUnitarity ()
    end select
    allocate(normalizer_%next)
    normalizer_ => normalizer_%next
    allocate(outputAnalysisDistributionNormalizerBinWidth  :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerBinWidth )
       normalizer_=outputAnalysisDistributionNormalizerBinWidth  ()
    end select
    allocate(outputAnalysisDistributionNormalizer_)
    outputAnalysisDistributionNormalizer_           =outputAnalysisDistributionNormalizerSequence      (normalizerSequence                                                                                                                   )
    ! Determine number of buffer bins.
    bufferCount=int(3.0d0*(colors(2)-colors(1))/sqrt(magnitudeErrorR**2+magnitudeErrorU**2),kind=c_size_t)+1_c_size_t
    ! Construct the object.
    write (distributionName,'(i2.2)') distributionNumber
    self%outputAnalysisVolumeFunction1D=                                                            &
         & outputAnalysisVolumeFunction1D(                                                          &
         &                                var_str('colorDistributionSDSS')//trim(distributionName), &
         &                                var_str('Distribution of SDSS u-r color'               ), &
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
         &                                outputAnalysisPropertyExtractorRatio_                   , &
         &                                outputAnalysisPropertyOperatorMagnitude_                , &
         &                                outputAnalysisPropertyOperatorIdentity_                 , &
         &                                outputAnalysisWeightOperator_                           , &
         &                                outputAnalysisDistributionOperator_                     , &
         &                                outputAnalysisDistributionNormalizer_                   , &
         &                                galacticFilter_                                         , &
         &                                outputAnalysisCovarianceModelPoisson                    , &
         &                                covarianceBinomialBinsPerDecade                         , &
         &                                covarianceBinomialMassHaloMinimum                       , &
         &                                covarianceBinomialMassHaloMaximum                         &
         &                               )
    ! Clean up.
    nullify(outputAnalysisWeightPropertyOperatorSequence_   )
    nullify(outputAnalysisPropertyOperatorMagnitude_        )
    nullify(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    nullify(outputAnalysisPropertyOperatorIdentity_         )
    nullify(outputAnalysisDistributionNormalizer_           )
    nullify(outputAnalysisDistributionOperator_             )
    nullify(outputAnalysisPropertyExtractorBandR_           )
    nullify(outputAnalysisPropertyExtractorBandU_           )
    nullify(outputAnalysisPropertyExtractorRatio_           )
    nullify(outputAnalysisWeightOperator_                   )
    nullify(propertyOperatorSequence                        )
    nullify(weightPropertyOperatorSequence                  )
    nullify(normalizerSequence                              )
    nullify(galacticFilter_                                 )
    nullify(cosmologyParametersData                         )
    nullify(cosmologyFunctionsData                          )
    return
  end function colorDistributionSDSSConstructorInternal

  subroutine colorDistributionSDSSDestructor(self)
    !% Destructor for the ``colorDistributionSDSS'' output analysis class.
    implicit none
    type(outputAnalysisColorDistributionSDSS), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologyFunctions_" />    
    return
  end subroutine colorDistributionSDSSDestructor
