!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements an HI mass function output analysis class.

  use Geometry_Surveys
  use Cosmology_Functions

  !# <outputAnalysis name="outputAnalysisMassFunctionHI" defaultThreadPrivate="yes">
  !#  <description>An HI mass function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisMassFunctionHI
     !% An HI mass function output analysis class.
     private
     class(surveyGeometryClass    ), pointer :: surveyGeometry_     => null()
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null(), cosmologyFunctionsData => null()
   contains
     final :: massFunctionHIDestructor     
  end type outputAnalysisMassFunctionHI

  interface outputAnalysisMassFunctionHI
     !% Constructors for the ``massFunctionHI'' output analysis class.
     module procedure massFunctionHIConstructorParameters
     module procedure massFunctionHIConstructorInternal
     module procedure massFunctionHIConstructorFile
  end interface outputAnalysisMassFunctionHI

contains

  function massFunctionHIConstructorParameters(parameters) result (self)
    !% Constructor for the ``massFunctionHI'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    use Output_Analysis_Molecular_Ratios
    implicit none
    type            (outputAnalysisMassFunctionHI           )                             :: self
    type            (inputParameters                        ), intent(inout)              :: parameters
    class           (galacticFilterClass                    ), pointer                    :: galacticFilter_
    class           (surveyGeometryClass                    ), pointer                    :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), pointer                    :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputTimesClass                       ), pointer                    :: outputTimes_
    class           (outputAnalysisDistributionOperatorClass), pointer                    :: outputAnalysisDistributionOperator_
    class           (outputAnalysisPropertyOperatorClass    ), pointer                    :: outputAnalysisPropertyOperator_
    class           (outputAnalysisMolecularRatioClass      ), pointer                    :: outputAnalysisMolecularRatio_
    double precision                                         , dimension(:) , allocatable :: masses
    integer                                                                               :: covarianceBinomialBinsPerDecade
    double precision                                                                      :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum
    type            (inputParameters                        )                             :: dataAnalysisParameters
    type            (varying_string                         )                             :: label                              , comment
    
    ! Check and read parameters.
    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    allocate(masses(parameters%count('masses')))
    !# <inputParameter>
    !#   <name>label</name>
    !#   <source>parameters</source>
    !#   <variable>label</variable>
    !#   <description>A label for the mass function.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>comment</name>
    !#   <source>parameters</source>
    !#   <variable>comment</variable>
    !#   <description>A descriptive comment for the mass function.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>masses</name>
    !#   <source>parameters</source>
    !#   <variable>masses</variable>
    !#   <description>The masses corresponding to bin centers.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialBinsPerDecade</name>
    !#   <source>parameters</source>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins per decade of halo mass to use when constructing HI mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing HI mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing HI mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="galacticFilter"                     name="galacticFilter_"                     source="parameters"            />
    !# <objectBuilder class="outputTimes"                        name="outputTimes_"                        source="parameters"            />
    !# <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"            />
    !# <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctionsData"              source="dataAnalysisParameters"/>
    !# <objectBuilder class="outputAnalysisPropertyOperator"     name="outputAnalysisPropertyOperator_"     source="parameters"            />
    !# <objectBuilder class="outputAnalysisDistributionOperator" name="outputAnalysisDistributionOperator_" source="parameters"            />
    !# <objectBuilder class="outputAnalysisMolecularRatio"       name="outputAnalysisMolecularRatio_"       source="parameters"            />
    !# <objectBuilder class="surveyGeometry"                     name="surveyGeometry_"                     source="parameters"            />
    self=outputAnalysisMassFunctionHI(label,comment,masses,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputAnalysisMolecularRatio_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum)
    !# <inputParametersValidate source="parameters"/>
    return
  end function massFunctionHIConstructorParameters

  function massFunctionHIConstructorFile(label,comment,fileName,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputAnalysisMolecularRatio_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result (self)
    !% Constructor for the ``massFunctionHI'' output analysis class which reads bin information from a standard format file.
    use IO_HDF5
    use Output_Analysis_Molecular_Ratios
    implicit none
    type            (outputAnalysisMassFunctionHI           )                             :: self
    type            (varying_string                         ), intent(in   )              :: label                              , comment
    character       (len=*                                  ), intent(in   )              :: fileName
    class           (galacticFilterClass                    ), intent(in   ), target      :: galacticFilter_
    class           (surveyGeometryClass                    ), intent(in   ), target      :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), intent(in   ), target      :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputTimesClass                       ), intent(inout), target      :: outputTimes_
    class           (outputAnalysisPropertyOperatorClass    ), intent(in   ), target      :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass), intent(in   ), target      :: outputAnalysisDistributionOperator_
    class           (outputAnalysisMolecularRatioClass      ), intent(in   ), target      :: outputAnalysisMolecularRatio_
    double precision                                         , dimension(:) , allocatable :: masses
    integer                                                                               :: covarianceBinomialBinsPerDecade
    double precision                                                                      :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum
    type            (hdf5Object                             )                             :: dataFile
    
    !$ call hdf5Access%set()
    call dataFile%openFile   (fileName,readOnly=.true.)
    call dataFile%readDataset('mass'  ,masses         )
    call dataFile%close      (                        )
    !$ call hdf5Access%unset()
    ! Construct the object.
    self=outputAnalysisMassFunctionHI(label,comment,masses,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputAnalysisMolecularRatio_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum)
    return
  end function massFunctionHIConstructorFile

  function massFunctionHIConstructorInternal(label,comment,masses,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputAnalysisMolecularRatio_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result(self)
    !% Constructor for the ``massFunctionHI'' output analysis class which takes a parameter set as input.
    use ISO_Varying_String
    use Memory_Management
    use Output_Times
    use String_Handling
    use Galacticus_Error
    use Output_Analyses_Options
    use Output_Analysis_Utilities
    use Output_Analysis_Molecular_Ratios
    use Numerical_Constants_Astronomical
    implicit none
    type            (outputAnalysisMassFunctionHI                   )                                :: self
    type            (varying_string                                 ), intent(in   )                 :: label                                                 , comment
    double precision                                                 , intent(in   ), dimension(:  ) :: masses
    class           (galacticFilterClass                            ), intent(in   ), target         :: galacticFilter_
    class           (surveyGeometryClass                            ), intent(in   ), target         :: surveyGeometry_
    class           (cosmologyFunctionsClass                        ), intent(in   ), target         :: cosmologyFunctions_                                   , cosmologyFunctionsData
    class           (outputTimesClass                               ), intent(inout), target         :: outputTimes_
    class           (outputAnalysisPropertyOperatorClass            ), intent(in   ), target         :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass        ), intent(in   ), target         :: outputAnalysisDistributionOperator_
    class           (outputAnalysisMolecularRatioClass              ), intent(in   ), target         :: outputAnalysisMolecularRatio_
    integer                                                          , intent(in   )                 :: covarianceBinomialBinsPerDecade
    double precision                                                 , intent(in   )                 :: covarianceBinomialMassHaloMinimum                     , covarianceBinomialMassHaloMaximum
    type            (outputAnalysisPropertyExtractorMassISM         )               , pointer        :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyOperatorHIMass           )               , pointer        :: outputAnalysisPropertyOperatorHIMass_
    type            (outputAnalysisPropertyOperatorLog10            )               , pointer        :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10        )               , pointer        :: outputAnalysisPropertyOperatorAntiLog10_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc)               , pointer        :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSequence         )               , pointer        :: outputAnalysisPropertyOperatorSequence_
    type            (outputAnalysisWeightOperatorCsmlgyVolume       )               , pointer        :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerSequence   )               , pointer        :: outputAnalysisDistributionNormalizer_
    type            (normalizerList                                 )               , pointer        :: normalizerSequence                                    , normalizer_
    type            (propertyOperatorList                           )               , pointer        :: propertyOperatorSequence
    double precision                                                 , allocatable  , dimension(:,:) :: outputWeight
    double precision                                                 , parameter                     :: bufferWidthLogarithmic                          =3.0d0
    integer         (c_size_t                                       ), parameter                     :: bufferCountMinimum                              =5
    integer         (c_size_t                                       )                                :: iBin                                                  , bufferCount
    !# <constructorAssign variables="*surveyGeometry_, *cosmologyFunctions_, *cosmologyFunctionsData"/>

    ! Compute weights that apply to each output redshift.
    self%binCount=size(masses,kind=c_size_t)
    call allocateArray(outputWeight,[self%binCount,outputTimes_%count()])
    do iBin=1,self%binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(self%surveyGeometry_,self%cosmologyFunctions_,outputTimes_,masses(iBin))
    end do
    ! Create a HI mass property extractor.
    allocate(outputAnalysisPropertyExtractor_)
    outputAnalysisPropertyExtractor_                =outputAnalysisPropertyExtractorMassISM         (                                                                 )
    ! Prepend log10, cosmological luminosity distance, and HI mass property operators.
    allocate(outputAnalysisPropertyOperatorHIMass_           )
    outputAnalysisPropertyOperatorHIMass_           =outputAnalysisPropertyOperatorHIMass           (outputAnalysisMolecularRatio_                                    )
    allocate(outputAnalysisPropertyOperatorLog10_            )
    outputAnalysisPropertyOperatorLog10_            =outputAnalysisPropertyOperatorLog10            (                                                                 )
    allocate(outputAnalysisPropertyOperatorAntiLog10_        )
    outputAnalysisPropertyOperatorAntiLog10_        =outputAnalysisPropertyOperatorAntiLog10        (                                                                 )
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_=outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc(cosmologyFunctions_          ,cosmologyFunctionsData,outputTimes_)
    select type (outputAnalysisPropertyOperator_)
    type is (outputAnalysisPropertyOperatorSequence)
       ! Existing property operator is a sequence operator - simply prepend our log10, and cosmological luminosity distance operators to it.
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorLog10_            )
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorHIMass_           )
       outputAnalysisPropertyOperatorSequence_ => outputAnalysisPropertyOperator_
    class default
       ! Existing operator is some other type - combine with our operators into a sequence operator.
       allocate(propertyOperatorSequence               )
       allocate(propertyOperatorSequence%next          )
       allocate(propertyOperatorSequence%next%next     )
       allocate(propertyOperatorSequence%next%next%next)
       propertyOperatorSequence               %operator_ => outputAnalysisPropertyOperatorHIMass_
       propertyOperatorSequence%next          %operator_ => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
       propertyOperatorSequence%next%next     %operator_ => outputAnalysisPropertyOperatorLog10_
       propertyOperatorSequence%next%next%next%operator_ => outputAnalysisPropertyOperator_
       allocate(outputAnalysisPropertyOperatorSequence_)
       outputAnalysisPropertyOperatorSequence_=outputAnalysisPropertyOperatorSequence(propertyOperatorSequence)
    end select
    ! Create a cosmological volume correction weight operator.
    allocate(outputAnalysisWeightOperator_)
    outputAnalysisWeightOperator_=outputAnalysisWeightOperatorCsmlgyVolume(cosmologyFunctions_,cosmologyFunctionsData,surveyGeometry_)
    ! Create a bin width distribution normalizer.
    allocate(normalizerSequence)
    normalizer_ => normalizerSequence
    allocate(outputAnalysisDistributionNormalizerBinWidth   :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerBinWidth  )
       normalizer_=outputAnalysisDistributionNormalizerBinWidth  ()
    end select
    allocate(normalizer_%next)
    normalizer_ => normalizer_%next
    allocate(outputAnalysisDistributionNormalizerLog10ToLog :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerLog10ToLog)
       normalizer_=outputAnalysisDistributionNormalizerLog10ToLog()
    end select
    allocate(outputAnalysisDistributionNormalizer_)
    outputAnalysisDistributionNormalizer_=outputAnalysisDistributionNormalizerSequence(normalizerSequence)
    ! Compute the number of buffer bins to add to either side of the mass function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    bufferCount=max(int(bufferWidthLogarithmic/log10(masses(2)/masses(1)))+1,bufferCountMinimum)
    ! Construct the object. We convert masses to log10(masses) here.
    self%outputAnalysisVolumeFunction1D=                                                           &
         & outputAnalysisVolumeFunction1D(                                                         &
         &                                'massFunctionHI'//label                                , &
         &                                comment                                                , &
         &                                var_str('massHI'                                      ), &
         &                                var_str('HI mass at the bin center'                   ), &
         &                                var_str('M☉'                                         ), &
         &                                massSolar                                              , &
         &                                var_str('massFunction'                                ), &
         &                                var_str('HI mass function averaged over each bin'     ), &
         &                                var_str('ᵪMpc⁻³'                                      ), &
         &                                megaParsec**(-3)                                       , &
         &                                log10(masses)                                          , &
         &                                bufferCount                                            , &
         &                                outputWeight                                           , &
         &                                outputAnalysisPropertyExtractor_                       , &
         &                                outputAnalysisPropertyOperatorSequence_                , &
         &                                outputAnalysisPropertyOperatorAntiLog10_               , &
         &                                outputAnalysisWeightOperator_                          , &
         &                                outputAnalysisDistributionOperator_                    , &
         &                                outputAnalysisDistributionNormalizer_                  , &
         &                                galacticFilter_                                        , &
         &                                outputTimes_                                           , &
         &                                outputAnalysisCovarianceModelBinomial                  , &
         &                                covarianceBinomialBinsPerDecade                        , &
         &                                covarianceBinomialMassHaloMinimum                      , &
         &                                covarianceBinomialMassHaloMaximum                        &
         &                               )
    ! Clean up.
    nullify(outputAnalysisPropertyExtractor_                )
    nullify(outputAnalysisPropertyOperatorLog10_            )
    nullify(outputAnalysisPropertyOperatorAntiLog10_        )
    nullify(outputAnalysisPropertyOperatorSequence_         )
    nullify(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    nullify(outputAnalysisDistributionNormalizer_           )
    nullify(outputAnalysisWeightOperator_                   )
    nullify(propertyOperatorSequence                        )
    nullify(normalizerSequence                              )
    return
  end function massFunctionHIConstructorInternal

  subroutine massFunctionHIDestructor(self)
    !% Destructor for  the ``massFunctionHI'' output analysis class.
    type(outputAnalysisMassFunctionHI), intent(inout) :: self
    
    !# <objectDestructor name="self%surveyGeometry_"       />
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%cosmologyFunctionsData"/>
    return
  end subroutine massFunctionHIDestructor
