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
Implements a luminosity function output analysis class.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Geometry_Surveys   , only : surveyGeometryClass

  !![
  <outputAnalysis name="outputAnalysisLuminosityFunction">
   <description>A luminosity function output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisLuminosityFunction
     !!{
     A luminosity function output analysis class.
     !!}
     private
     class           (surveyGeometryClass    ), pointer                   :: surveyGeometry_     => null()
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_ => null(), cosmologyFunctionsData => null()
     double precision                         , allocatable, dimension(:) :: magnitudesAbsolute
   contains
     final :: luminosityFunctionDestructor
  end type outputAnalysisLuminosityFunction

  interface outputAnalysisLuminosityFunction
     !!{
     Constructors for the \refClass{outputAnalysisLuminosityFunction} output analysis class.
     !!}
     module procedure luminosityFunctionConstructorParameters
     module procedure luminosityFunctionConstructorInternal
     module procedure luminosityFunctionConstructorFile
  end interface outputAnalysisLuminosityFunction

contains

  function luminosityFunctionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLuminosityFunction} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisLuminosityFunction       )                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    class           (outputTimesClass                       ), pointer                     :: outputTimes_
    class           (galacticFilterClass                    ), pointer                     :: galacticFilter_
    class           (surveyGeometryClass                    ), pointer                     :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), pointer                     :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputAnalysisDistributionOperatorClass), pointer                     :: outputAnalysisDistributionOperator_
    class           (outputAnalysisPropertyOperatorClass    ), pointer                     :: outputAnalysisPropertyOperator_
    double precision                                         , dimension(:  ), allocatable :: magnitudesAbsolute                 , functionValueTarget              , &
         &                                                                                    functionCovarianceTarget1D
    double precision                                         , dimension(:,:), allocatable :: functionCovarianceTarget
    integer                                                                                :: covarianceBinomialBinsPerDecade
    double precision                                                                       :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum, &
         &                                                                                    redshiftBand
    type            (inputParameters                        )                              :: dataAnalysisParameters
    type            (varying_string                         )                              :: label                              , comment                          , &
         &                                                                                    filterName                         , filterType                       , &
         &                                                                                    targetLabel

    ! Check and read parameters.
    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    allocate(magnitudesAbsolute(parameters%count('magnitudesAbsolute')))
    !![
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <variable>label</variable>
      <description>A label for the luminosity function.</description>
    </inputParameter>
    <inputParameter>
      <name>comment</name>
      <source>parameters</source>
      <variable>comment</variable>
      <description>A descriptive comment for the luminosity function.</description>
    </inputParameter>
    <inputParameter>
      <name>magnitudesAbsolute</name>
      <source>parameters</source>
      <variable>magnitudesAbsolute</variable>
      <description>The absolute magnitudes corresponding to bin centers.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing luminosity function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing luminosity function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing luminosity function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    !!]
    if (parameters%isPresent('targetLabel')) then
       !![
       <inputParameter>
         <name>targetLabel</name>
         <source>parameters</source>
         <description>Label for the target dataset.</description>
       </inputParameter>
       !!]
    end if
    if (parameters%isPresent('functionValueTarget')) then
       if (parameters%isPresent('functionCovarianceTarget')) then
          !![
          <inputParameter>
            <name>functionValueTarget</name>
            <source>parameters</source>
            <description>The target function for likelihood calculations.</description>
          </inputParameter>
          <inputParameter>
            <name>functionCovarianceTarget</name>
            <source>parameters</source>
            <variable>functionCovarianceTarget1D</variable>
            <description>The target function covariance for likelihood calculations.</description>
          </inputParameter>
          !!]
          if (size(functionCovarianceTarget1D) == size(functionValueTarget)**2) then
             allocate(functionCovarianceTarget(size(functionValueTarget),size(functionValueTarget)))
             functionCovarianceTarget=reshape(functionCovarianceTarget1D,shape(functionCovarianceTarget))
          else
             call Error_Report('functionCovariance has wrong size'//{introspection:location})
          end if
       else
          call Error_Report('functionCovariance must be specified if functionTarget is present'//{introspection:location})
       end if
    else
       if (parameters%isPresent('functionCovariance')) call Error_Report('functionTarget must be specified if functionCovariance is present'//{introspection:location})
    end if
    !![
    <objectBuilder class="galacticFilter"                     name="galacticFilter_"                     source="parameters"            />
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"            />
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctionsData"              source="dataAnalysisParameters"/>
    <objectBuilder class="outputAnalysisPropertyOperator"     name="outputAnalysisPropertyOperator_"     source="parameters"            />
    <objectBuilder class="outputAnalysisDistributionOperator" name="outputAnalysisDistributionOperator_" source="parameters"            />
    <objectBuilder class="surveyGeometry"                     name="surveyGeometry_"                     source="parameters"            />
    <objectBuilder class="outputTimes"                        name="outputTimes_"                        source="parameters"            />
    <conditionalCall>
     <call>self=outputAnalysisLuminosityFunction(label,comment,magnitudesAbsolute,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,char(filterName),char(filterType),redshiftBand{conditions})</call>
     <argument name="targetLabel"              value="targetLabel"              parameterPresent="parameters"/>
     <argument name="functionValueTarget"      value="functionValueTarget"      parameterPresent="parameters"/>
     <argument name="functionCovarianceTarget" value="functionCovarianceTarget" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"                    />
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="cosmologyFunctionsData"             />
    <objectDestructor name="outputAnalysisPropertyOperator_"    />
    <objectDestructor name="outputAnalysisDistributionOperator_"/>
    <objectDestructor name="surveyGeometry_"                    />
    <objectDestructor name="outputTimes_"                       />
    !!]
    return
  end function luminosityFunctionConstructorParameters

  function luminosityFunctionConstructorFile(label,comment,fileName,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,filterName,filterType,redshiftBand) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisLuminosityFunction} output analysis class which reads bin information from a standard format file.
    !!}
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    type            (outputAnalysisLuminosityFunction       )                              :: self
    type            (varying_string                         ), intent(in   )               :: label                              , comment
    character       (len=*                                  ), intent(in   )               :: fileName
    class           (galacticFilterClass                    ), intent(in   ) , target      :: galacticFilter_
    class           (surveyGeometryClass                    ), intent(in   ) , target      :: surveyGeometry_
    class           (outputTimesClass                       ), intent(inout) , target      :: outputTimes_
    class           (cosmologyFunctionsClass                ), intent(in   ) , target      :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputAnalysisPropertyOperatorClass    ), intent(inout) , target      :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass), intent(in   ) , target      :: outputAnalysisDistributionOperator_
    character       (len=*                                  ), intent(in   )               :: filterName                         , filterType
    double precision                                         , intent(in   ) , optional    :: redshiftBand
    double precision                                         , dimension(:  ), allocatable :: magnitudesAbsolute                 , functionValueTarget              , &
         &                                                                                    functionErrorTarget
    double precision                                         , dimension(:,:), allocatable :: functionCovarianceTarget
    integer                                                  , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                         , intent(in   )               :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum
    integer                                                                                :: i
    type            (hdf5Object                             )                              :: dataFile
    type            (varying_string                         )                              :: targetLabel
    logical                                                                                :: haveTarget

    !$ call hdf5Access%set()
    call dataFile%openFile   (fileName           ,readOnly=.true.   )
    call dataFile%readDataset('magnitudeAbsolute',magnitudesAbsolute)
    haveTarget=dataFile%hasDataset('luminosityFunction').and.dataFile%hasDataset('luminosityFunctionError')
    if (haveTarget) then
       call dataFile%readAttribute('label'                  ,targetLabel        )
       call dataFile%readDataset  ('luminosityFunction'     ,functionValueTarget)
       call dataFile%readDataset  ('luminosityFunctionError',functionErrorTarget)
    end if
    call dataFile%close      (                                      )
    !$ call hdf5Access%unset()
    if (haveTarget) then
       allocate(functionCovarianceTarget(size(functionErrorTarget),size(functionErrorTarget)))
       functionCovarianceTarget=0.0d0
       do i=1,size(functionErrorTarget)
          functionCovarianceTarget(i,i)=functionErrorTarget(i)**2
       end do
    end if
    ! Construct the object.
    !![
    <conditionalCall>
     <call>self=outputAnalysisLuminosityFunction(label,comment,magnitudesAbsolute,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,filterName,filterType,redshiftBand{conditions})</call>
     <argument name="targetLabel"              value="targetLabel"              condition="haveTarget"/>
     <argument name="functionValueTarget"      value="functionValueTarget"      condition="haveTarget"/>
     <argument name="functionCovarianceTarget" value="functionCovarianceTarget" condition="haveTarget"/>
    </conditionalCall>
    !!]
    return
  end function luminosityFunctionConstructorFile

  function luminosityFunctionConstructorInternal(label,comment,magnitudesAbsolute,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,filterName,filterType,redshiftBand,targetLabel,functionValueTarget,functionCovarianceTarget) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisLuminosityFunction} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions                     , only : cosmologyFunctionsClass
    use :: Galactic_Filters                        , only : galacticFilterClass
    use :: Geometry_Surveys                        , only : surveyGeometryClass
    use :: ISO_Varying_String                      , only : var_str                                     , varying_string
    use :: Node_Property_Extractors                , only : nodePropertyExtractorLmnstyStllrCF2000
    use :: Numerical_Constants_Astronomical        , only : megaParsec
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerBinWidth
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorClass
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorClass         , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorIdentity, outputAnalysisPropertyOperatorMagnitude, &
          &                                                 outputAnalysisPropertyOperatorSequence      , propertyOperatorList
    use :: Output_Analysis_Utilities               , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorCsmlgyVolume
    implicit none
    type            (outputAnalysisLuminosityFunction                )                                          :: self
    type            (varying_string                                  ), intent(in   )                           :: label                                                 , comment
    double precision                                                  , intent(in   )          , dimension(:  ) :: magnitudesAbsolute
    class           (galacticFilterClass                             ), intent(in   ), target                   :: galacticFilter_
    class           (surveyGeometryClass                             ), intent(in   ), target                   :: surveyGeometry_
    class           (outputTimesClass                                ), intent(inout), target                   :: outputTimes_
    class           (cosmologyFunctionsClass                         ), intent(in   ), target                   :: cosmologyFunctions_                                   , cosmologyFunctionsData
    class           (outputAnalysisPropertyOperatorClass             ), intent(inout), target                   :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass         ), intent(in   ), target                   :: outputAnalysisDistributionOperator_
    integer                                                           , intent(in   )                           :: covarianceBinomialBinsPerDecade
    double precision                                                  , intent(in   )                           :: covarianceBinomialMassHaloMinimum                     , covarianceBinomialMassHaloMaximum
    character       (len=*                                           ), intent(in   )                           :: filterName                                            , filterType
    double precision                                                  , intent(in   ), optional                 :: redshiftBand
    type            (varying_string                                  ), intent(in   ), optional                 :: targetLabel
    double precision                                                  , intent(in   ), optional, dimension(:  ) :: functionValueTarget
    double precision                                                  , intent(in   ), optional, dimension(:,:) :: functionCovarianceTarget
    type            (nodePropertyExtractorLmnstyStllrCF2000          )               , pointer                  :: nodePropertyExtractor_
    type            (outputAnalysisPropertyOperatorMagnitude         )               , pointer                  :: outputAnalysisPropertyOperatorMagnitude_
    type            (outputAnalysisPropertyOperatorIdentity          )               , pointer                  :: outputAnalysisPropertyOperatorIdentity_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc )               , pointer                  :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSequence          )               , pointer                  :: outputAnalysisPropertyOperatorSequence_
    type            (outputAnalysisWeightOperatorCsmlgyVolume        )               , pointer                  :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerBinWidth    )               , pointer                  :: outputAnalysisDistributionNormalizer_
    type            (propertyOperatorList                            )               , pointer                  :: propertyOperatorSequence
    double precision                                                  , allocatable            , dimension(:,:) :: outputWeight
    double precision                                                  , parameter                               :: bufferWidth                                     =7.5d0
    integer         (c_size_t                                        ), parameter                               :: bufferCountMinimum                              =5
    integer         (c_size_t                                        )                                          :: iBin                                                  , bufferCount
    !![
    <constructorAssign variables="magnitudesAbsolute, *surveyGeometry_, *cosmologyFunctions_, *cosmologyFunctionsData"/>
    !!]

    ! Compute weights that apply to each output redshift.
    self%binCount=size(magnitudesAbsolute,kind=c_size_t)
    allocate(outputWeight(self%binCount,outputTimes_%count()))
    do iBin=1,self%binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(self%surveyGeometry_,self%cosmologyFunctions_,outputTimes_,magnitudeAbsoluteLimit=magnitudesAbsolute(iBin))
    end do
    ! Create a luminosity property extractor.
    allocate(nodePropertyExtractor_)
    !![
    <referenceConstruct object="nodePropertyExtractor_"                           constructor="nodePropertyExtractorLmnstyStllrCF2000(filterName         ,filterType  ,depthOpticalISMCoefficient=1.0d0,depthOpticalCloudsCoefficient=1.0d0,wavelengthExponent          =0.7d0,outputTimes_=outputTimes_,redshiftBand=redshiftBand,outputMask=sum(outputWeight,dim=1) > 0.0d0)"/>
    !!]
    ! Prepend magnitude and cosmological luminosity distance property operators.
    allocate(outputAnalysisPropertyOperatorMagnitude_        )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorMagnitude_"         constructor="outputAnalysisPropertyOperatorMagnitude         (                                                                                                  )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorIdentity_         )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorIdentity_"          constructor="outputAnalysisPropertyOperatorIdentity          (                                                                                                  )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc (cosmologyFunctions_,cosmologyFunctionsData,outputTimes_                                           )"/>
    !!]
    select type (outputAnalysisPropertyOperator_)
    type is (outputAnalysisPropertyOperatorSequence)
       ! Existing property operator is a sequence operator - simply prepend our magnitude and cosmological luminosity distance operators to it.
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorMagnitude_        )
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
       !![
       <referenceAcquire target="outputAnalysisPropertyOperatorSequence_" source="outputAnalysisPropertyOperator_"/>
       !!]
    class default
       ! Existing operator is some other type - combine with our operators into a sequence operator.
       allocate(propertyOperatorSequence          )
       allocate(propertyOperatorSequence%next     )
       allocate(propertyOperatorSequence%next%next)
       propertyOperatorSequence          %operator_ => outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
       propertyOperatorSequence%next     %operator_ => outputAnalysisPropertyOperatorMagnitude_
       propertyOperatorSequence%next%next%operator_ => outputAnalysisPropertyOperator_
       allocate(outputAnalysisPropertyOperatorSequence_)
       !![
       <referenceConstruct object="outputAnalysisPropertyOperatorSequence_" constructor="outputAnalysisPropertyOperatorSequence(propertyOperatorSequence)"/>
       !!]
    end select
    ! Create a cosmological volume correction weight operator.
    allocate(outputAnalysisWeightOperator_)
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_" constructor="outputAnalysisWeightOperatorCsmlgyVolume(cosmologyFunctions_,cosmologyFunctionsData,surveyGeometry_)"/>
    !!]
    ! Create a bin width distribution normalizer.
    allocate(outputAnalysisDistributionNormalizer_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizer_" constructor="outputAnalysisDistributionNormalizerBinWidth()"/>
    !!]
    ! Compute the number of buffer bins to add to either side of the luminosity function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    bufferCount=max(int(bufferWidth/(magnitudesAbsolute(2)-magnitudesAbsolute(1)))+1,bufferCountMinimum)
    ! Construct the object.
    self%outputAnalysisVolumeFunction1D=                                                            &
         & outputAnalysisVolumeFunction1D(                                                          &
         &                                'luminosityFunction'//label                             , &
         &                                comment                                                 , &
         &                                var_str('magnitudeAbsolute'                            ), &
         &                                var_str('absolute magnitude at the bin center'         ), &
         &                                var_str(' '                                            ), &
         &                                0.0d0                                                   , &
         &                                var_str('luminosityFunction'                           ), &
         &                                var_str('luminosity function averaged over each bin'   ), &
         &                                var_str('ᵪMpc⁻³'                                       ), &
         &                                megaParsec**(-3)                                        , &
         &                                magnitudesAbsolute                                      , &
         &                                bufferCount                                             , &
         &                                outputWeight                                            , &
         &                                nodePropertyExtractor_                                  , &
         &                                outputAnalysisPropertyOperatorSequence_                 , &
         &                                outputAnalysisPropertyOperatorIdentity_                 , &
         &                                outputAnalysisWeightOperator_                           , &
         &                                outputAnalysisDistributionOperator_                     , &
         &                                outputAnalysisDistributionNormalizer_                   , &
         &                                galacticFilter_                                         , &
         &                                outputTimes_                                            , &
         &                                outputAnalysisCovarianceModelBinomial                   , &
         &                                covarianceBinomialBinsPerDecade                         , &
         &                                covarianceBinomialMassHaloMinimum                       , &
         &                                covarianceBinomialMassHaloMaximum                       , &
         &                                .false.                                                 , &
         &                                var_str('$M$'                                          ), &
         &                                var_str('$\mathrm{d}n/\mathrm{d}M$ [$_\chi$Mpc$^{-3}$]'), &
         &                                .false.                                                 , &
         &                                .true.                                                  , &
         &                                targetLabel                                             , &
         &                                functionValueTarget                                     , &
         &                                functionCovarianceTarget                                  &
         &                               )
    ! Clean up.
    !![
    <objectDestructor name="nodePropertyExtractor_"                          />
    <objectDestructor name="outputAnalysisPropertyOperatorMagnitude_"        />
    <objectDestructor name="outputAnalysisPropertyOperatorIdentity_"         />
    <objectDestructor name="outputAnalysisPropertyOperatorSequence_"         />
    <objectDestructor name="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_"/>
    <objectDestructor name="outputAnalysisDistributionNormalizer_"           />
    <objectDestructor name="outputAnalysisWeightOperator_"                   />
    !!]
    nullify(propertyOperatorSequence)
    return
  end function luminosityFunctionConstructorInternal

  subroutine luminosityFunctionDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisLuminosityFunction} output analysis class.
    !!}
    type(outputAnalysisLuminosityFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%surveyGeometry_"       />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%cosmologyFunctionsData"/>
    !!]
    return
  end subroutine luminosityFunctionDestructor

