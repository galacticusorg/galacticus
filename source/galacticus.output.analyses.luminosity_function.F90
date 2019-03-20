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

!% Contains a module which implements a luminosity function output analysis class.

  use Geometry_Surveys
  use Cosmology_Functions

  !# <outputAnalysis name="outputAnalysisLuminosityFunction">
  !#  <description>A luminosity function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisLuminosityFunction
     !% A luminosity function output analysis class.
     private
     class(surveyGeometryClass    ), pointer :: surveyGeometry_     => null()
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null(), cosmologyFunctionsData => null()
   contains
     final :: luminosityFunctionDestructor
  end type outputAnalysisLuminosityFunction

  interface outputAnalysisLuminosityFunction
     !% Constructors for the ``luminosityFunction'' output analysis class.
     module procedure luminosityFunctionConstructorParameters
     module procedure luminosityFunctionConstructorInternal
     module procedure luminosityFunctionConstructorFile
  end interface outputAnalysisLuminosityFunction

contains

  function luminosityFunctionConstructorParameters(parameters) result (self)
    !% Constructor for the ``luminosityFunction'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisLuminosityFunction       )                             :: self
    type            (inputParameters                        ), intent(inout)              :: parameters
    class           (outputTimesClass                       ), pointer                    :: outputTimes_
    class           (galacticFilterClass                    ), pointer                    :: galacticFilter_
    class           (surveyGeometryClass                    ), pointer                    :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), pointer                    :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputAnalysisDistributionOperatorClass), pointer                    :: outputAnalysisDistributionOperator_
    class           (outputAnalysisPropertyOperatorClass    ), pointer                    :: outputAnalysisPropertyOperator_
    double precision                                         , dimension(:) , allocatable :: magnitudesAbsolute
    integer                                                                               :: covarianceBinomialBinsPerDecade
    double precision                                                                      :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum, &
         &                                                                                   redshiftBand
    type            (inputParameters                        )                             :: dataAnalysisParameters
    type            (varying_string                         )                             :: label                              , comment                          , &
         &                                                                                   filterName                         , filterType

    ! Check and read parameters.
    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    allocate(magnitudesAbsolute(parameters%count('magnitudesAbsolute')))
    !# <inputParameter>
    !#   <name>label</name>
    !#   <source>parameters</source>
    !#   <variable>label</variable>
    !#   <description>A label for the luminosity function.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>comment</name>
    !#   <source>parameters</source>
    !#   <variable>comment</variable>
    !#   <description>A descriptive comment for the luminosity function.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>magnitudesAbsolute</name>
    !#   <source>parameters</source>
    !#   <variable>masses</variable>
    !#   <description>The absolute magnitudes corresponding to bin centers.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialBinsPerDecade</name>
    !#   <source>parameters</source>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins per decade of halo mass to use when constructing luminosity function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing luminosity function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing luminosity function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="galacticFilter"                     name="galacticFilter_"                     source="parameters"            />
    !# <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"            />
    !# <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctionsData"              source="dataAnalysisParameters"/>
    !# <objectBuilder class="outputAnalysisPropertyOperator"     name="outputAnalysisPropertyOperator_"     source="parameters"            />
    !# <objectBuilder class="outputAnalysisDistributionOperator" name="outputAnalysisDistributionOperator_" source="parameters"            />
    !# <objectBuilder class="surveyGeometry"                     name="surveyGeometry_"                     source="parameters"            />
    !# <objectBuilder class="outputTimes"                        name="outputTimes_"                        source="parameters"            />
    self=outputAnalysisLuminosityFunction(label,comment,magnitudesAbsolute,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,char(filterName),char(filterType),redshiftBand)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="galacticFilter_"                    />
    !# <objectDestructor name="cosmologyFunctions_"                />
    !# <objectDestructor name="cosmologyFunctionsData"             />
    !# <objectDestructor name="outputAnalysisPropertyOperator_"    />
    !# <objectDestructor name="outputAnalysisDistributionOperator_"/>
    !# <objectDestructor name="surveyGeometry_"                    />
    !# <objectDestructor name="outputTimes_"                       />
    return
  end function luminosityFunctionConstructorParameters

  function luminosityFunctionConstructorFile(label,comment,fileName,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,filterName,filterType,redshiftBand) result (self)
    !% Constructor for the ``luminosityFunction'' output analysis class which reads bin information from a standard format file.
    use IO_HDF5
    implicit none
    type            (outputAnalysisLuminosityFunction       )                             :: self
    type            (varying_string                         ), intent(in   )              :: label                              , comment
    character       (len=*                                  ), intent(in   )              :: fileName
    class           (galacticFilterClass                    ), intent(in   ), target      :: galacticFilter_
    class           (surveyGeometryClass                    ), intent(in   ), target      :: surveyGeometry_
    class           (outputTimesClass                       ), intent(in   ), target      :: outputTimes_
    class           (cosmologyFunctionsClass                ), intent(in   ), target      :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputAnalysisPropertyOperatorClass    ), intent(in   ), target      :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass), intent(in   ), target      :: outputAnalysisDistributionOperator_
    character       (len=*                                  ), intent(in   )              :: filterName                         , filterType
    double precision                                         , intent(in   ), optional    :: redshiftBand
    double precision                                         , dimension(:) , allocatable :: magnitudesAbsolute
    integer                                                                               :: covarianceBinomialBinsPerDecade
    double precision                                                                      :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum
    type            (hdf5Object                             )                             :: dataFile
    
    !$ call hdf5Access%set()
    call dataFile%openFile   (fileName           ,readOnly=.true.   )
    call dataFile%readDataset('magnitudeAbsolute',magnitudesAbsolute)
    call dataFile%close      (                                      )
    !$ call hdf5Access%unset()
    ! Construct the object.
    self=outputAnalysisLuminosityFunction(label,comment,magnitudesAbsolute,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,filterName,filterType,redshiftBand)
    return
  end function luminosityFunctionConstructorFile

  function luminosityFunctionConstructorInternal(label,comment,magnitudesAbsolute,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,filterName,filterType,redshiftBand) result(self)
    !% Constructor for the ``luminosityFunction'' output analysis class which takes a parameter set as input.
    use ISO_Varying_String
    use Memory_Management
    use String_Handling
    use Galacticus_Error
    use Numerical_Constants_Astronomical
    use Output_Analyses_Options
    use Output_Analysis_Utilities
    implicit none
    type            (outputAnalysisLuminosityFunction                )                                :: self
    type            (varying_string                                  ), intent(in   )                 :: label                                                 , comment
    double precision                                                  , intent(in   ), dimension(:  ) :: magnitudesAbsolute
    class           (galacticFilterClass                             ), intent(in   ), target         :: galacticFilter_
    class           (surveyGeometryClass                             ), intent(in   ), target         :: surveyGeometry_
    class           (outputTimesClass                                ), intent(in   ), target         :: outputTimes_
    class           (cosmologyFunctionsClass                         ), intent(in   ), target         :: cosmologyFunctions_                                   , cosmologyFunctionsData
    class           (outputAnalysisPropertyOperatorClass             ), intent(in   ), target         :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass         ), intent(in   ), target         :: outputAnalysisDistributionOperator_
    integer                                                           , intent(in   )                 :: covarianceBinomialBinsPerDecade
    double precision                                                  , intent(in   )                 :: covarianceBinomialMassHaloMinimum                     , covarianceBinomialMassHaloMaximum
    character       (len=*                                           ), intent(in   )                 :: filterName                                            , filterType
    double precision                                                  , intent(in   ), optional       :: redshiftBand
    type            (outputAnalysisPropertyExtractorLmnstyStllrCF2000)               , pointer        :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyOperatorMagnitude         )               , pointer        :: outputAnalysisPropertyOperatorMagnitude_
    type            (outputAnalysisPropertyOperatorIdentity          )               , pointer        :: outputAnalysisPropertyOperatorIdentity_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc )               , pointer        :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSequence          )               , pointer        :: outputAnalysisPropertyOperatorSequence_
    type            (outputAnalysisWeightOperatorCsmlgyVolume        )               , pointer        :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerBinWidth    )               , pointer        :: outputAnalysisDistributionNormalizer_
    type            (propertyOperatorList                            )               , pointer        :: propertyOperatorSequence
    double precision                                                  , allocatable  , dimension(:,:) :: outputWeight
    double precision                                                  , parameter                     :: bufferWidth                                     =7.5d0
    integer         (c_size_t                                        ), parameter                     :: bufferCountMinimum                              =5
    integer         (c_size_t                                        )                                :: iBin                                                  , bufferCount
    !# <constructorAssign variables="*surveyGeometry_, *cosmologyFunctions_, *cosmologyFunctionsData, *outputTimes_"/>

    ! Compute weights that apply to each output redshift.
    self%binCount=size(magnitudesAbsolute,kind=c_size_t)
    call allocateArray(outputWeight,[self%binCount,self%outputTimes_%count()])
    do iBin=1,self%binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(self%surveyGeometry_,self%cosmologyFunctions_,self%outputTimes_,magnitudeAbsoluteLimit=magnitudesAbsolute(iBin))
    end do
    ! Create a luminosity property extractor.
    allocate(outputAnalysisPropertyExtractor_)
    outputAnalysisPropertyExtractor_                =outputAnalysisPropertyExtractorLmnstyStllrCF2000(filterName         ,filterType  ,depthOpticalISMCoefficient=1.0d0,depthOpticalCloudsCoefficient=1.0d0,wavelengthExponent          =0.7d0,outputTimes_=outputTimes_,redshiftBand=redshiftBand,outputMask=sum(outputWeight,dim=1) > 0.0d0)
    ! Prepend magnitude and cosmological luminosity distance property operators.
    allocate(outputAnalysisPropertyOperatorMagnitude_        )
    outputAnalysisPropertyOperatorMagnitude_        =outputAnalysisPropertyOperatorMagnitude         (                                                                                                  )
    allocate(outputAnalysisPropertyOperatorIdentity_         )
    outputAnalysisPropertyOperatorIdentity_         =outputAnalysisPropertyOperatorIdentity          (                                                                                                  )
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_=outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc (cosmologyFunctions_,cosmologyFunctionsData,outputTimes_                                           )
    select type (outputAnalysisPropertyOperator_)
    type is (outputAnalysisPropertyOperatorSequence)
       ! Existing property operator is a sequence operator - simply prepend our magnitude and cosmological luminosity distance operators to it.
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorMagnitude_        )
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
       outputAnalysisPropertyOperatorSequence_ => outputAnalysisPropertyOperator_
    class default
       ! Existing operator is some other type - combine with our operators into a sequence operator.
       allocate(propertyOperatorSequence          )
       allocate(propertyOperatorSequence%next     )
       allocate(propertyOperatorSequence%next%next)
       propertyOperatorSequence          %operator_ => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
       propertyOperatorSequence%next     %operator_ => outputAnalysisPropertyOperatorMagnitude_
       propertyOperatorSequence%next%next%operator_ => outputAnalysisPropertyOperator_
       allocate(outputAnalysisPropertyOperatorSequence_)
       outputAnalysisPropertyOperatorSequence_=outputAnalysisPropertyOperatorSequence(propertyOperatorSequence)
    end select
    ! Create a cosmological volume correction weight operator.
    allocate(outputAnalysisWeightOperator_)
    outputAnalysisWeightOperator_=outputAnalysisWeightOperatorCsmlgyVolume(cosmologyFunctions_,cosmologyFunctionsData,surveyGeometry_)
    ! Create a bin width distribution normalizer.
    allocate(outputAnalysisDistributionNormalizer_)
    outputAnalysisDistributionNormalizer_=outputAnalysisDistributionNormalizerBinWidth()
    ! Compute the number of buffer bins to add to either side of the luminosity function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    bufferCount=max(int(bufferWidth/(magnitudesAbsolute(2)-magnitudesAbsolute(1)))+1,bufferCountMinimum)
    ! Construct the object.
    self%outputAnalysisVolumeFunction1D=                                                         &
         & outputAnalysisVolumeFunction1D(                                                       &
         &                                'luminosityFunction'//label                          , &
         &                                comment                                              , &
         &                                var_str('magnitudeAbsolute'                         ), &
         &                                var_str('absolute magnitude at the bin center'      ), &
         &                                var_str(' '                                         ), &
         &                                0.0d0                                                , &
         &                                var_str('luminosityFunction'                        ), &
         &                                var_str('luminosity function averaged over each bin'), &
         &                                var_str('ᵪMpc⁻³'                                    ), &
         &                                megaParsec**(-3)                                     , &
         &                                magnitudesAbsolute                                   , &
         &                                bufferCount                                          , &
         &                                outputWeight                                         , &
         &                                outputAnalysisPropertyExtractor_                     , &
         &                                outputAnalysisPropertyOperatorSequence_              , &
         &                                outputAnalysisPropertyOperatorIdentity_              , &
         &                                outputAnalysisWeightOperator_                        , &
         &                                outputAnalysisDistributionOperator_                  , &
         &                                outputAnalysisDistributionNormalizer_                , &
         &                                galacticFilter_                                      , &
         &                                outputTimes_                                         , &
         &                                outputAnalysisCovarianceModelBinomial                , &
         &                                covarianceBinomialBinsPerDecade                      , &
         &                                covarianceBinomialMassHaloMinimum                    , &
         &                                covarianceBinomialMassHaloMaximum                      &
         &                               )
    ! Clean up.
    nullify(outputAnalysisPropertyExtractor_                )
    nullify(outputAnalysisPropertyOperatorMagnitude_        )
    nullify(outputAnalysisPropertyOperatorIdentity_         )
    nullify(outputAnalysisPropertyOperatorSequence_         )
    nullify(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    nullify(outputAnalysisDistributionNormalizer_           )
    nullify(outputAnalysisWeightOperator_                   )
    nullify(propertyOperatorSequence                        )
    return
  end function luminosityFunctionConstructorInternal

  subroutine luminosityFunctionDestructor(self)
    !% Destructor for  the ``luminosityFunction'' output analysis class.
    type(outputAnalysisLuminosityFunction), intent(inout) :: self
    
    !# <objectDestructor name="self%surveyGeometry_"       />
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%cosmologyFunctionsData"/>
    return
  end subroutine luminosityFunctionDestructor
  
