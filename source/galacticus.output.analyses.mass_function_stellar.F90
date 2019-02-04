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

!% Contains a module which implements a stellar mass function output analysis class.

  use Geometry_Surveys
  use Cosmology_Functions

  !# <outputAnalysis name="outputAnalysisMassFunctionStellar">
  !#  <description>A stellar mass function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisMassFunctionStellar
     !% A massFunctionStellar output analysis class.
     private
     class(surveyGeometryClass    ), pointer :: surveyGeometry_     => null()
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null(), cosmologyFunctionsData => null()
   contains
     final :: massFunctionStellarDestructor     
  end type outputAnalysisMassFunctionStellar

  interface outputAnalysisMassFunctionStellar
     !% Constructors for the ``massFunctionStellar'' output analysis class.
     module procedure massFunctionStellarConstructorParameters
     module procedure massFunctionStellarConstructorInternal
     module procedure massFunctionStellarConstructorFile
  end interface outputAnalysisMassFunctionStellar

contains

  function massFunctionStellarConstructorParameters(parameters) result (self)
    !% Constructor for the ``massFunctionStellar'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type            (outputAnalysisMassFunctionStellar      )                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    class           (galacticFilterClass                    ), pointer                     :: galacticFilter_
    class           (surveyGeometryClass                    ), pointer                     :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), pointer                     :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputAnalysisDistributionOperatorClass), pointer                     :: outputAnalysisDistributionOperator_
    class           (outputAnalysisPropertyOperatorClass    ), pointer                     :: outputAnalysisPropertyOperator_
    class           (outputTimesClass                       ), pointer                     :: outputTimes_ 
    double precision                                         , dimension(:  ), allocatable :: masses                             , functionValueTarget              , &
         &                                                                                    functionCovarianceTarget1D
    double precision                                         , dimension(:,:), allocatable :: functionCovarianceTarget
    integer                                                                                :: covarianceBinomialBinsPerDecade
    double precision                                                                       :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum
    type            (inputParameters                        )                              :: dataAnalysisParameters
    type            (varying_string                         )                              :: label                              , comment
    
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
    !#   <description>The number of bins per decade of halo mass to use when constructing stellar mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing stellar mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing stellar mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    if (parameters%isPresent('functionValueTarget')) then
       if (parameters%isPresent('functionCovarianceTarget')) then
          !# <inputParameter>
          !#   <name>functionValueTarget</name>
          !#   <source>parameters</source>
          !#   <description>The target function for likelihood calculations.</description>
          !#   <type>real</type>
          !#   <cardinality>0..1</cardinality>
          !# </inputParameter> 
          !# <inputParameter>
          !#   <name>functionCovarianceTarget</name>
          !#   <source>parameters</source>
          !#   <variable>functionCovarianceTarget1D</variable>
          !#   <description>The target function covariance for likelihood calculations.</description>
          !#   <type>real</type>
          !#   <cardinality>0..1</cardinality>
          !# </inputParameter>
          if (size(functionCovarianceTarget1D) == size(functionValueTarget)**2) then
             allocate(functionCovarianceTarget(size(functionValueTarget),size(functionValueTarget)))
             functionCovarianceTarget=reshape(functionCovarianceTarget1D,shape(functionCovarianceTarget))
          else
             call Galacticus_Error_Report('functionCovariance has wrong size'//{introspection:location})
          end if
       else
          call Galacticus_Error_Report('functionCovariance must be specified if functionTarget is present'//{introspection:location})
       end if
    else
       if (parameters%isPresent('functionCovariance')) call Galacticus_Error_Report('functionTarget must be specified if functionCovariance is present'//{introspection:location})
    end if
    !# <objectBuilder class="galacticFilter"                     name="galacticFilter_"                     source="parameters"            />
    !# <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"            />
    !# <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctionsData"              source="dataAnalysisParameters"/>
    !# <objectBuilder class="outputAnalysisPropertyOperator"     name="outputAnalysisPropertyOperator_"     source="parameters"            />
    !# <objectBuilder class="outputAnalysisDistributionOperator" name="outputAnalysisDistributionOperator_" source="parameters"            />
    !# <objectBuilder class="surveyGeometry"                     name="surveyGeometry_"                     source="parameters"            />
    !# <objectBuilder class="outputTimes"                        name="outputTimes_"                        source="parameters"          />
    !# <conditionalCall>
    !#  <call>self=outputAnalysisMassFunctionStellar(label,comment,masses,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum{conditions})</call>
    !#  <argument name="functionValueTarget"      value="functionValueTarget"      parameterPresent="parameters"/>
    !#  <argument name="functionCovarianceTarget" value="functionCovarianceTarget" parameterPresent="parameters"/>
    !# </conditionalCall>
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="galacticFilter_"                    />
    !# <objectDestructor name="cosmologyFunctions_"                />
    !# <objectDestructor name="cosmologyFunctionsData"             />
    !# <objectDestructor name="outputAnalysisPropertyOperator_"    />
    !# <objectDestructor name="outputAnalysisDistributionOperator_"/>
    !# <objectDestructor name="surveyGeometry_"                    />
    !# <objectDestructor name="outputTimes_"                       />
    return
  end function massFunctionStellarConstructorParameters

  function massFunctionStellarConstructorFile(label,comment,fileName,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result (self)
    !% Constructor for the ``massFunctionStellar'' output analysis class which reads bin information from a standard format file.
    use IO_HDF5
    implicit none
    type            (outputAnalysisMassFunctionStellar      )                              :: self
    type            (varying_string                         ), intent(in   )               :: label                              , comment
    character       (len=*                                  ), intent(in   )               :: fileName
    class           (galacticFilterClass                    ), intent(in   ) , target      :: galacticFilter_
    class           (surveyGeometryClass                    ), intent(in   ) , target      :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), intent(in   ) , target      :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputAnalysisPropertyOperatorClass    ), intent(in   ) , target      :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass), intent(in   ) , target      :: outputAnalysisDistributionOperator_
    class           (outputTimesClass                       ), intent(in   ) , target      :: outputTimes_
    double precision                                         , dimension(:  ), allocatable :: masses                             , functionValueTarget
    double precision                                         , dimension(:,:), allocatable :: functionCovarianceTarget
    integer                                                                                :: covarianceBinomialBinsPerDecade
    double precision                                                                       :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum
    type            (hdf5Object                             )                              :: dataFile
    logical                                                                                :: haveTarget
    
    !$ call hdf5Access%set()
    call dataFile%openFile   (fileName,readOnly=.true.)
    call dataFile%readDataset('mass'  ,masses         )
    haveTarget=dataFile%hasDataset('massFunctionObserved').and.dataFile%hasDataset('covariance')
    if (haveTarget) then
       call dataFile%readDataset('massFunctionObserved',functionValueTarget     )
       call dataFile%readDataset('covariance'          ,functionCovarianceTarget)
    end if
    call dataFile%close      (                        )
    !$ call hdf5Access%unset()
    ! Construct the object.
    !# <conditionalCall>
    !#  <call>self=outputAnalysisMassFunctionStellar(label,comment,masses,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum{conditions})</call>
    !#  <argument name="functionValueTarget"      value="functionValueTarget"      condition="haveTarget"/>
    !#  <argument name="functionCovarianceTarget" value="functionCovarianceTarget" condition="haveTarget"/>
    !# </conditionalCall>
    return
  end function massFunctionStellarConstructorFile

  function massFunctionStellarConstructorInternal(label,comment,masses,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,functionValueTarget,functionCovarianceTarget) result(self)
    !% Constructor for the ``massFunctionStellar'' output analysis class which takes a parameter set as input.
    use ISO_Varying_String
    use Memory_Management
    use String_Handling
    use Galacticus_Error
    use Numerical_Constants_Astronomical
    use Output_Analyses_Options
    use Output_Analysis_Utilities
    implicit none
    type            (outputAnalysisMassFunctionStellar              )                                          :: self
    type            (varying_string                                 ), intent(in   )                           :: label                                                 , comment
    double precision                                                 , intent(in   )          , dimension(:  ) :: masses
    class           (galacticFilterClass                            ), intent(in   ), target                   :: galacticFilter_
    class           (surveyGeometryClass                            ), intent(in   ), target                   :: surveyGeometry_
    class           (cosmologyFunctionsClass                        ), intent(in   ), target                   :: cosmologyFunctions_                                   , cosmologyFunctionsData
    class           (outputAnalysisPropertyOperatorClass            ), intent(in   ), target                   :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass        ), intent(in   ), target                   :: outputAnalysisDistributionOperator_
    class           (outputTimesClass                               ), intent(in   ), target                   :: outputTimes_
    integer                                                          , intent(in   )                           :: covarianceBinomialBinsPerDecade
    double precision                                                 , intent(in   )                           :: covarianceBinomialMassHaloMinimum                     , covarianceBinomialMassHaloMaximum
    double precision                                                 , intent(in   ), optional, dimension(:  ) :: functionValueTarget
    double precision                                                 , intent(in   ), optional, dimension(:,:) :: functionCovarianceTarget
    type            (outputAnalysisPropertyExtractorMassStellar     )               , pointer                  :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyOperatorLog10            )               , pointer                  :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10        )               , pointer                  :: outputAnalysisPropertyOperatorAntiLog10_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc)               , pointer                  :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSequence         )               , pointer                  :: outputAnalysisPropertyOperatorSequence_
    type            (outputAnalysisWeightOperatorCsmlgyVolume       )               , pointer                  :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerSequence   )               , pointer                  :: outputAnalysisDistributionNormalizer_
    type            (normalizerList                                 )               , pointer                  :: normalizerSequence                                    , normalizer_
    type            (propertyOperatorList                           )               , pointer                  :: propertyOperatorSequence
    double precision                                                 , allocatable            , dimension(:,:) :: outputWeight
    double precision                                                 , parameter                               :: bufferWidthLogarithmic                          =3.0d0
    integer         (c_size_t                                       ), parameter                               :: bufferCountMinimum                              =5
    integer         (c_size_t                                       )                                          :: iBin                                                  , bufferCount
    !# <constructorAssign variables="*surveyGeometry_, *cosmologyFunctions_, *cosmologyFunctionsData, *outputTimes_"/>

    ! Compute weights that apply to each output redshift.
    self%binCount=size(masses,kind=c_size_t)
    call allocateArray(outputWeight,[self%binCount,self%outputTimes_%count()])
    do iBin=1,self%binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(self%surveyGeometry_,self%cosmologyFunctions_,self%outputTimes_,masses(iBin))
    end do
    ! Create a stellar mass property extractor.
    allocate(outputAnalysisPropertyExtractor_)
    outputAnalysisPropertyExtractor_                =outputAnalysisPropertyExtractorMassStellar     (                                                       )
    ! Prepend log10 and cosmological luminosity distance property operators.
    allocate(outputAnalysisPropertyOperatorLog10_            )
    outputAnalysisPropertyOperatorLog10_            =outputAnalysisPropertyOperatorLog10            (                                                       )
    allocate(outputAnalysisPropertyOperatorAntiLog10_        )
    outputAnalysisPropertyOperatorAntiLog10_        =outputAnalysisPropertyOperatorAntiLog10        (                                                       )
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_=outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc(cosmologyFunctions_,cosmologyFunctionsData,outputTimes_)
    select type (outputAnalysisPropertyOperator_)
    type is (outputAnalysisPropertyOperatorSequence)
       ! Existing property operator is a sequence operator - simply prepend our log10 and cosmological luminosity distance operators to it.
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorLog10_            )
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
       outputAnalysisPropertyOperatorSequence_ => outputAnalysisPropertyOperator_
    class default
       ! Existing operator is some other type - combine with our operators into a sequence operator.
       allocate(propertyOperatorSequence          )
       allocate(propertyOperatorSequence%next     )
       allocate(propertyOperatorSequence%next%next)
       propertyOperatorSequence          %operator_ => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
       propertyOperatorSequence%next     %operator_ => outputAnalysisPropertyOperatorLog10_
       propertyOperatorSequence%next%next%operator_ => outputAnalysisPropertyOperator_
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
         &                                'massFunctionStellar'//label                           , &
         &                                comment                                                , &
         &                                var_str('massStellar'                                 ), &
         &                                var_str('Stellar mass at the bin center'              ), &
         &                                var_str('M☉'                                         ), &
         &                                massSolar                                              , &
         &                                var_str('massFunction'                                ), &
         &                                var_str('Stellar mass function averaged over each bin'), &
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
         &                                covarianceBinomialMassHaloMaximum                      , &
         &                                functionValueTarget                                    , &
         &                                functionCovarianceTarget                                 &
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
  end function massFunctionStellarConstructorInternal

  subroutine massFunctionStellarDestructor(self)
    !% Destructor for  the ``massFunctionStellar'' output analysis class.
    type(outputAnalysisMassFunctionStellar), intent(inout) :: self
    
    !# <objectDestructor name="self%surveyGeometry_"       />
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%cosmologyFunctionsData"/>
    return
  end subroutine massFunctionStellarDestructor
