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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Implements a luminosity function output analysis class.
!!}

  use :: Cosmology_Functions             , only : cosmologyFunctionsClass
  use :: Geometry_Surveys                , only : surveyGeometryClass
  use :: Dust_Attenuations               , only : dustAttenuationClass
  use :: Star_Formation_Histories        , only : starFormationHistoryClass        , starFormationHistoryNull
  use :: HII_Region_Luminosity_Functions , only : hiiRegionLuminosityFunctionClass
  use :: HII_Region_Mass_Functions       , only : hiiRegionMassFunctionClass
  use :: HII_Region_Density_Distributions, only : hiiRegionDensityDistributionClass
  use :: HII_Region_Escape_Fraction      , only : hiiRegionEscapeFractionClass

  !![
  <outputAnalysis name="outputAnalysisLuminosityFunctionHalpha" docformat="rst">
   <description>
   Computes the H\ :math:`\alpha` emission line galaxy luminosity function in bins of line luminosity, optionally including [NII] doublet contamination (``includeNitrogenII``), applying ISM dust attenuation, and with configurable binomial covariance matrix parameters for halo mass range and optional target dataset.

   Line luminosities are computed by :galacticus-class:`nodePropertyExtractorLuminosityEmissionLine`, which convolves
   the star formation history of each galaxy with tabulated Cloudy models of the H**ii** regions described by the
   ``[hiiRegionLuminosityFunction]``, ``[hiiRegionMassFunction]``, ``[hiiRegionDensityDistribution]``, and
   ``[hiiRegionEscapeFraction]`` objects. A star formation history must therefore be stored (see
   ``[starFormationHistory]``); the ``null`` default stores none and no line emission results.

   Note that the tabulation selected by ``[cloudyTableFileName]`` may already include the effect of dust *within* the
   H**ii** region, depending on its dust-to-metals ratio---all of the currently distributed tabulations do. The
   ``[dustAttenuation]`` model given here is applied on top, and so should describe only the dust *outside* the
   region: configuring a birth cloud component in addition would count the same dust twice.
   </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisLuminosityFunctionHalpha
     !!{RST
     A luminosity function output analysis class.
     !!}
     private
     class           (surveyGeometryClass              ), pointer                   :: surveyGeometry_               => null()
     class           (cosmologyFunctionsClass          ), pointer                   :: cosmologyFunctions_           => null(), cosmologyFunctionsData => null()
     class           (dustAttenuationClass             ), pointer                   :: dustAttenuation_              => null()
     class           (starFormationHistoryClass        ), pointer                   :: starFormationHistory_         => null()
     class           (hiiRegionLuminosityFunctionClass ), pointer                   :: hiiRegionLuminosityFunction_  => null()
     class           (hiiRegionMassFunctionClass       ), pointer                   :: hiiRegionMassFunction_        => null()
     class           (hiiRegionDensityDistributionClass), pointer                   :: hiiRegionDensityDistribution_ => null()
     class           (hiiRegionEscapeFractionClass     ), pointer                   :: hiiRegionEscapeFraction_      => null()
     double precision                                   , allocatable, dimension(:) :: luminosities
     type            (varying_string                   )                            :: cloudyTableFileName
     double precision                                                               :: toleranceRelative
     logical                                                                        :: includeNitrogenII
   contains
     final :: luminosityFunctionHalphaDestructor
  end type outputAnalysisLuminosityFunctionHalpha

  interface outputAnalysisLuminosityFunctionHalpha
     !!{RST
     Constructors for the :galacticus-class:`outputAnalysisLuminosityFunctionHalpha` output analysis class.
     !!}
     module procedure luminosityFunctionHalphaConstructorParameters
     module procedure luminosityFunctionHalphaConstructorInternal
     module procedure luminosityFunctionHalphaConstructorFile
  end interface outputAnalysisLuminosityFunctionHalpha

contains

  function luminosityFunctionHalphaConstructorParameters(parameters) result (self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisLuminosityFunctionHalpha` output analysis class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisLuminosityFunctionHalpha )                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    class           (galacticFilterClass                    ), pointer                     :: galacticFilter_
    class           (surveyGeometryClass                    ), pointer                     :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), pointer                     :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputTimesClass                       ), pointer                     :: outputTimes_
    class           (outputAnalysisDistributionOperatorClass), pointer                     :: outputAnalysisDistributionOperator_
    class           (outputAnalysisPropertyOperatorClass    ), pointer                     :: outputAnalysisPropertyOperator_
    class           (starFormationHistoryClass              ), pointer                     :: starFormationHistory_
    class           (hiiRegionLuminosityFunctionClass       ), pointer                     :: hiiRegionLuminosityFunction_
    class           (hiiRegionMassFunctionClass             ), pointer                     :: hiiRegionMassFunction_
    class           (hiiRegionDensityDistributionClass      ), pointer                     :: hiiRegionDensityDistribution_
    class           (hiiRegionEscapeFractionClass           ), pointer                     :: hiiRegionEscapeFraction_
    class           (dustAttenuationClass                   ), pointer                     :: dustAttenuation_
    double precision                                         , dimension(:  ), allocatable :: luminosities                       , functionValueTarget              , &
         &                                                                                    functionCovarianceTarget1D
    double precision                                         , dimension(:,:), allocatable :: functionCovarianceTarget
    integer                                                                                :: covarianceBinomialBinsPerDecade
    double precision                                                                       :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum, &
         &                                                                                    toleranceRelative
    type            (inputParameters                        )                              :: dataAnalysisParameters
    type            (varying_string                         )                              :: label                              , comment                          , &
         &                                                                                    targetLabel                        , cloudyTableFileName
    logical                                                                                :: includeNitrogenII

    ! Check and read parameters.
    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    allocate(luminosities(parameters%count('luminosities')))
    !![
    <inputParameter docformat="rst">
      <name>label</name>
      <source>parameters</source>
      <description>
      A label for the luminosity function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>comment</name>
      <source>parameters</source>
      <description>
      A descriptive comment for the luminosity function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>luminosities</name>
      <source>parameters</source>
      <description>
      The luminosities corresponding to bin centers.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>
      The number of bins per decade of halo mass to use when constructing luminosity function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d8</defaultValue>
      <description>
      The minimum halo mass to consider when constructing luminosity function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d16</defaultValue>
      <description>
      The maximum halo mass to consider when constructing luminosity function covariance matrices for main branch galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>includeNitrogenII</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, include contamination by the [NII] (6548Å :math:`+` 6584Å) doublet.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>cloudyTableFileName</name>
      <source>parameters</source>
      <defaultValue>var_str('%DATASTATICPATH%/hiiRegions/emissionLineLuminosities_BC2003_highResolution_imfChabrier.hdf5')</defaultValue>
      <description>
      The file of tabulated emission line luminosities to use.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>toleranceRelative</name>
      <source>parameters</source>
      <defaultValue>1.0d-3</defaultValue>
      <description>
      The relative tolerance used in integration over stellar population spectra when computing line luminosities.
      </description>
    </inputParameter>
    !!]
    if (parameters%isPresent('targetLabel')) then
       !![
       <inputParameter docformat="rst">
         <name>targetLabel</name>
         <source>parameters</source>
         <description>
         Label for the target dataset.
         </description>
       </inputParameter>
       !!]
    end if
    if (parameters%isPresent('functionValueTarget')) then
       if (parameters%isPresent('functionCovarianceTarget')) then
          !![
          <inputParameter docformat="rst">
            <name>functionValueTarget</name>
            <source>parameters</source>
            <description>
            The target function for likelihood calculations.
            </description>
          </inputParameter>
          <inputParameter docformat="rst">
            <name>functionCovarianceTarget</name>
            <source>parameters</source>
            <variable>functionCovarianceTarget1D</variable>
            <description>
            The target function covariance for likelihood calculations.
            </description>
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
    <objectBuilder class="outputTimes"                        name="outputTimes_"                        source="parameters"            />
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"            />
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctionsData"              source="dataAnalysisParameters"/>
    <objectBuilder class="outputAnalysisPropertyOperator"     name="outputAnalysisPropertyOperator_"     source="parameters"            />
    <objectBuilder class="outputAnalysisDistributionOperator" name="outputAnalysisDistributionOperator_" source="parameters"            />
    <objectBuilder class="surveyGeometry"                     name="surveyGeometry_"                     source="parameters"            />
    <objectBuilder class="starFormationHistory"               name="starFormationHistory_"               source="parameters"            />
    <objectBuilder class="hiiRegionLuminosityFunction"        name="hiiRegionLuminosityFunction_"        source="parameters"            />
    <objectBuilder class="hiiRegionMassFunction"              name="hiiRegionMassFunction_"              source="parameters"            />
    <objectBuilder class="hiiRegionDensityDistribution"       name="hiiRegionDensityDistribution_"       source="parameters"            />
    <objectBuilder class="hiiRegionEscapeFraction"            name="hiiRegionEscapeFraction_"            source="parameters"            />
    <objectBuilder class="dustAttenuation"                    name="dustAttenuation_"                    source="parameters"            />
    <conditionalCall>
     <call>self=outputAnalysisLuminosityFunctionHalpha(label,comment,luminosities,includeNitrogenII,cloudyTableFileName,toleranceRelative,galacticFilter_,surveyGeometry_,dustAttenuation_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,starFormationHistory_,hiiRegionLuminosityFunction_,hiiRegionMassFunction_,hiiRegionDensityDistribution_,hiiRegionEscapeFraction_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum{conditions})</call>
     <argument name="targetLabel"              value="targetLabel"              parameterPresent="parameters"/>
     <argument name="functionValueTarget"      value="functionValueTarget"      parameterPresent="parameters"/>
     <argument name="functionCovarianceTarget" value="functionCovarianceTarget" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"                    />
    <objectDestructor name="outputTimes_"                       />
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="cosmologyFunctionsData"             />
    <objectDestructor name="outputAnalysisPropertyOperator_"    />
    <objectDestructor name="outputAnalysisDistributionOperator_"/>
    <objectDestructor name="surveyGeometry_"                    />
    <objectDestructor name="starFormationHistory_"              />
    <objectDestructor name="hiiRegionLuminosityFunction_"       />
    <objectDestructor name="hiiRegionMassFunction_"             />
    <objectDestructor name="hiiRegionDensityDistribution_"      />
    <objectDestructor name="hiiRegionEscapeFraction_"           />
    <objectDestructor name="dustAttenuation_"                   />
    !!]
    return
  end function luminosityFunctionHalphaConstructorParameters

  function luminosityFunctionHalphaConstructorFile(label,comment,fileName,includeNitrogenII,cloudyTableFileName,toleranceRelative,galacticFilter_,surveyGeometry_,dustAttenuation_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,starFormationHistory_,hiiRegionLuminosityFunction_,hiiRegionMassFunction_,hiiRegionDensityDistribution_,hiiRegionEscapeFraction_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result (self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisLuminosityFunctionHalpha` output analysis class which reads bin information from a standard format file.
    !!}
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5File
    implicit none
    type            (outputAnalysisLuminosityFunctionHalpha )                              :: self
    type            (varying_string                         ), intent(in   )               :: label                              , comment
    character       (len=*                                  ), intent(in   )               :: fileName
    logical                                                  , intent(in   )               :: includeNitrogenII
    type            (varying_string                         ), intent(in   )               :: cloudyTableFileName
    double precision                                         , intent(in   )               :: toleranceRelative
    class           (galacticFilterClass                    ), intent(in   ) , target      :: galacticFilter_
    class           (surveyGeometryClass                    ), intent(in   ) , target      :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), intent(in   ) , target      :: cosmologyFunctions_                , cosmologyFunctionsData
    class           (outputTimesClass                       ), intent(inout) , target      :: outputTimes_
    class           (outputAnalysisPropertyOperatorClass    ), intent(inout) , target      :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass), intent(in   ) , target      :: outputAnalysisDistributionOperator_
    class           (starFormationHistoryClass              ), intent(in   ) , target      :: starFormationHistory_
    class           (hiiRegionLuminosityFunctionClass       ), intent(in   ) , target      :: hiiRegionLuminosityFunction_
    class           (hiiRegionMassFunctionClass             ), intent(in   ) , target      :: hiiRegionMassFunction_
    class           (hiiRegionDensityDistributionClass      ), intent(in   ) , target      :: hiiRegionDensityDistribution_
    class           (hiiRegionEscapeFractionClass           ), intent(in   ) , target      :: hiiRegionEscapeFraction_
    class           (dustAttenuationClass                   ), intent(in   ) , target      :: dustAttenuation_
    double precision                                         , dimension(:  ), allocatable :: luminosities                       , functionValueTarget              , &
         &                                                                                    functionErrorTarget
    double precision                                         , dimension(:,:), allocatable :: functionCovarianceTarget
    integer                                                  , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                         , intent(in   )               :: covarianceBinomialMassHaloMinimum  , covarianceBinomialMassHaloMaximum
    integer                                                                                :: i
    type            (hdf5File                               )                              :: dataFile
    type            (varying_string                         )                              :: targetLabel
    logical                                                                                :: haveTarget

    !$ call hdf5Access%set()
    dataFile=hdf5File(fileName,readOnly=.true.)
    call    dataFile%readDataset  ('luminosity'             ,luminosities       )
    haveTarget=dataFile%hasDataset('luminosityFunction').and.dataFile%hasDataset('luminosityFunctionError')
    if (haveTarget) then
       call dataFile%readAttribute('label'                  ,targetLabel        )
       call dataFile%readDataset  ('luminosityFunction'     ,functionValueTarget)
       call dataFile%readDataset  ('luminosityFunctionError',functionErrorTarget)
    end if
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
     <call>self=outputAnalysisLuminosityFunctionHalpha(label,comment,luminosities,includeNitrogenII,cloudyTableFileName,toleranceRelative,galacticFilter_,surveyGeometry_,dustAttenuation_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,starFormationHistory_,hiiRegionLuminosityFunction_,hiiRegionMassFunction_,hiiRegionDensityDistribution_,hiiRegionEscapeFraction_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum{conditions})</call>
     <argument name="targetLabel"              value="targetLabel"              condition="haveTarget"/>
     <argument name="functionValueTarget"      value="functionValueTarget"      condition="haveTarget"/>
     <argument name="functionCovarianceTarget" value="functionCovarianceTarget" condition="haveTarget"/>
    </conditionalCall>
    !!]
    return
  end function luminosityFunctionHalphaConstructorFile

  function luminosityFunctionHalphaConstructorInternal(label,comment,luminosities,includeNitrogenII,cloudyTableFileName,toleranceRelative,galacticFilter_,surveyGeometry_,dustAttenuation_,cosmologyFunctions_,cosmologyFunctionsData,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputTimes_,starFormationHistory_,hiiRegionLuminosityFunction_,hiiRegionMassFunction_,hiiRegionDensityDistribution_,hiiRegionEscapeFraction_,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,targetLabel,functionValueTarget,functionCovarianceTarget) result(self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisLuminosityFunctionHalpha` output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions                     , only : cosmologyFunctionsClass
    use :: Galactic_Filters                        , only : galacticFilterClass
    use :: Galactic_Structure_Options              , only : componentTypeAll
    use :: Error                                   , only : Error_Report
    use :: Geometry_Surveys                        , only : surveyGeometryClass
    use :: ISO_Varying_String                      , only : var_str                                        , varying_string
    use :: Dust_Attenuations                       , only : dustAttenuationClass
    use :: Node_Property_Extractors                , only : multiExtractorList                             , nodePropertyExtractorDustAttenuation        , nodePropertyExtractorLuminosityEmissionLine    , nodePropertyExtractorScalarizer
    use :: Numerical_Constants_Astronomical        , only : megaParsec
    use :: Numerical_Constants_Units               , only : ergs
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Normalizers, only : normalizerList                                 , outputAnalysisDistributionNormalizerBinWidth, outputAnalysisDistributionNormalizerLog10ToLog , outputAnalysisDistributionNormalizerSequence
    use :: Output_Analysis_Target_Data             , only : outputAnalysisTargetDataStandard
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorClass
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10        , outputAnalysisPropertyOperatorClass         , outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorLog10         , &
          &                                                 outputAnalysisPropertyOperatorSequence         , propertyOperatorList
    use :: Output_Analysis_Utilities               , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorCosmologyVolume
    use :: Output_Times                            , only : outputTimesClass
    implicit none
    type            (outputAnalysisLuminosityFunctionHalpha         )                                          :: self
    type            (varying_string                                 ), intent(in   )                           :: label                                                 , comment
    logical                                                          , intent(in   )                           :: includeNitrogenII
    type            (varying_string                                 ), intent(in   )                           :: cloudyTableFileName
    double precision                                                 , intent(in   )                           :: toleranceRelative
    double precision                                                 , intent(in   )          , dimension(:  ) :: luminosities
    class           (galacticFilterClass                            ), intent(in   ), target                   :: galacticFilter_
    class           (surveyGeometryClass                            ), intent(in   ), target                   :: surveyGeometry_
    class           (cosmologyFunctionsClass                        ), intent(in   ), target                   :: cosmologyFunctions_                                   , cosmologyFunctionsData
    class           (outputTimesClass                               ), intent(inout), target                   :: outputTimes_
    class           (outputAnalysisPropertyOperatorClass            ), intent(inout), target                   :: outputAnalysisPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass        ), intent(in   ), target                   :: outputAnalysisDistributionOperator_
    class           (dustAttenuationClass                           ), intent(in   ), target                   :: dustAttenuation_
    class           (starFormationHistoryClass                      ), intent(in   ), target                   :: starFormationHistory_
    class           (hiiRegionLuminosityFunctionClass               ), intent(in   ), target                   :: hiiRegionLuminosityFunction_
    class           (hiiRegionMassFunctionClass                     ), intent(in   ), target                   :: hiiRegionMassFunction_
    class           (hiiRegionDensityDistributionClass              ), intent(in   ), target                   :: hiiRegionDensityDistribution_
    class           (hiiRegionEscapeFractionClass                   ), intent(in   ), target                   :: hiiRegionEscapeFraction_
    integer                                                          , intent(in   )                           :: covarianceBinomialBinsPerDecade
    double precision                                                 , intent(in   )                           :: covarianceBinomialMassHaloMinimum                     , covarianceBinomialMassHaloMaximum
    type            (varying_string                                 ), intent(in   ), optional                 :: targetLabel
    double precision                                                 , intent(in   ), optional, dimension(:  ) :: functionValueTarget
    double precision                                                 , intent(in   ), optional, dimension(:,:) :: functionCovarianceTarget
    type            (nodePropertyExtractorLuminosityEmissionLine    )               , pointer                  :: nodePropertyExtractorLines_
    type            (nodePropertyExtractorDustAttenuation           )               , pointer                  :: nodePropertyExtractorAttenuated_
    type            (nodePropertyExtractorScalarizer                )               , pointer                  :: nodePropertyExtractor_
    type            (multiExtractorList                             )               , pointer                  :: extractors                                            , extractor_
    type            (multiExtractorList                             ), allocatable            , dimension(:  ) :: extractorsLines
    type            (outputAnalysisPropertyOperatorLog10            )               , pointer                  :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10        )               , pointer                  :: outputAnalysisPropertyOperatorAntiLog10_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc)               , pointer                  :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSequence         )               , pointer                  :: outputAnalysisPropertyOperatorSequence_
    type            (outputAnalysisWeightOperatorCosmologyVolume       )               , pointer                  :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerSequence   )               , pointer                  :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionNormalizerBinWidth   )               , pointer                  :: outputAnalysisDistributionNormalizerBinWidth_
    type            (outputAnalysisDistributionNormalizerLog10ToLog )               , pointer                  :: outputAnalysisDistributionNormalizerLog10ToLog_
    type            (normalizerList                                 )               , pointer                  :: normalizerSequence
    type            (propertyOperatorList                           )               , pointer                  :: propertyOperatorSequence
    double precision                                                 , allocatable            , dimension(:,:) :: outputWeight
    type            (varying_string                                 ), allocatable            , dimension(:  ) :: lineNames
    double precision                                                 , parameter                               :: bufferWidth                                     =1.0d0
    integer         (c_size_t                                       ), parameter                               :: bufferCountMinimum                              =5
    integer         (c_size_t                                       )                                          :: iBin                                                  , bufferCount
    integer                                                                                                   :: iLine
    type            (outputAnalysisTargetDataStandard)                              :: outputAnalysisTargetData_
    !![
    <constructorAssign variables="luminosities, includeNitrogenII, cloudyTableFileName, toleranceRelative, *surveyGeometry_, *cosmologyFunctions_, *cosmologyFunctionsData, *dustAttenuation_, *starFormationHistory_, *hiiRegionLuminosityFunction_, *hiiRegionMassFunction_, *hiiRegionDensityDistribution_, *hiiRegionEscapeFraction_"/>
    !!]

    ! Compute weights that apply to each output redshift.
    self%binCount=size(luminosities,kind=c_size_t)
    allocate(outputWeight(self%binCount,outputTimes_%count()))
    do iBin=1,self%binCount
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(self%surveyGeometry_,self%cosmologyFunctions_,outputTimes_,luminosity=luminosities(iBin))
    end do
    ! Create a luminosity property extractor. The line luminosities themselves are unattenuated; dust is applied by
    ! wrapping them, which lets any `dustAttenuation` model be used and keeps a single implementation of the physics.
    ! The wrapper emits only the sum, formed elementwise over its children, so one child is built per line, each
    ! emitting a single element. The sum is then a single value --- the summed, attenuated luminosity of the lines ---
    ! which a scalarizer presents as the scalar this analysis requires.
    if (includeNitrogenII) then
       allocate(lineNames(3))
       lineNames(1)=var_str('balmerAlpha6565')
       lineNames(2)=var_str('nitrogenII6550' )
       lineNames(3)=var_str('nitrogenII6585' )
    else
       allocate(lineNames(1))
       lineNames(1)=var_str('balmerAlpha6565')
    end if
    ! Line luminosities are computed from the star formation history, so without one this analysis would run happily
    ! and report a luminosity function of zero in every bin. Since `starFormationHistory` defaults to `null`, that is
    ! what an otherwise unchanged parameter file would get, so refuse it here rather than return an empty analysis.
    select type (starFormationHistory_)
    type is (starFormationHistoryNull)
       call Error_Report(                                                                                        &
            &            'this analysis computes emission line luminosities from the star formation history,'//  &
            &            ' so a [starFormationHistory] which stores one must be specified - the default of'  //  &
            &            ' `null` stores none, which would give a luminosity function of zero in every bin'  //  &
            &            {introspection:location}                                                                &
            &           )
    class default
       ! A star formation history is stored, so line luminosities can be computed.
    end select
    ! `extractorsLines` keeps a handle on each line extractor which does not depend on the linked list, since that
    ! list belongs to the wrapper once it has been passed to it.
    allocate(extractorsLines(size(lineNames)))
    extractors => null()
    extractor_ => null()
    do iLine=1,size(lineNames)
       if (associated(extractor_)) then
          allocate(extractor_%next)
          extractor_ => extractor_%next
       else
          allocate(extractors      )
          extractor_ => extractors
       end if
       allocate(nodePropertyExtractorLines_)
       !![
       <referenceConstruct object="nodePropertyExtractorLines_" constructor="nodePropertyExtractorLuminosityEmissionLine(cloudyTableFileName,componentTypeAll,lineNames(iLine:iLine),toleranceRelative,starFormationHistory_,outputTimes_,hiiRegionLuminosityFunction_,hiiRegionMassFunction_,hiiRegionDensityDistribution_,hiiRegionEscapeFraction_)"/>
       !!]
       extractor_            %extractor_ => nodePropertyExtractorLines_
       extractorsLines(iLine)%extractor_ => nodePropertyExtractorLines_
    end do
    allocate(nodePropertyExtractorAttenuated_)
    allocate(nodePropertyExtractor_          )
    !![
    <referenceConstruct object="nodePropertyExtractorAttenuated_" constructor="nodePropertyExtractorDustAttenuation(dustAttenuation_,.false.,.true.,.true.,var_str('luminosityEmissionLine'),extractors)"/>
    <referenceConstruct object="nodePropertyExtractor_"           constructor="nodePropertyExtractorScalarizer     (1,1,nodePropertyExtractorAttenuated_)"/>
    !!]
    ! Release our own references to the line extractors --- the wrapper holds its own.
    do iLine=1,size(lineNames)
       !![
       <objectDestructor name="extractorsLines(iLine)%extractor_"/>
       !!]
    end do
    ! Prepend log10 and cosmological luminosity distance property operators.
    allocate(outputAnalysisPropertyOperatorLog10_            )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                                                                                                      )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorAntiLog10_        )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorAntiLog10_"         constructor="outputAnalysisPropertyOperatorAntiLog10        (                                                                                                                                      )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc(cosmologyFunctions_           ,cosmologyFunctionsData,outputTimes_                                                                    )"/>
    !!]
    select type (outputAnalysisPropertyOperator_)
    type is (outputAnalysisPropertyOperatorSequence)
       ! Existing property operator is a sequence operator - simply prepend our magnitude and cosmological luminosity distance operators to it.
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorLog10_            )
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
       propertyOperatorSequence%next     %operator_ => outputAnalysisPropertyOperatorLog10_
       propertyOperatorSequence%next%next%operator_ => outputAnalysisPropertyOperator_
       allocate(outputAnalysisPropertyOperatorSequence_)
       !![
       <referenceConstruct object="outputAnalysisPropertyOperatorSequence_" constructor="outputAnalysisPropertyOperatorSequence(propertyOperatorSequence)"/>
       !!]
    end select
    ! Create a cosmological volume correction weight operator.
    allocate(outputAnalysisWeightOperator_)
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"                    constructor="outputAnalysisWeightOperatorCosmologyVolume       (cosmologyFunctions_,cosmologyFunctionsData                    ,surveyGeometry_)"/>
    !!]
    ! Create a bin width distribution normalizer.
    allocate(outputAnalysisDistributionNormalizerBinWidth_  )
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerBinWidth_"   constructor="outputAnalysisDistributionNormalizerBinWidth  ()"/>
    !!]
    allocate(outputAnalysisDistributionNormalizerLog10ToLog_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerLog10ToLog_" constructor="outputAnalysisDistributionNormalizerLog10ToLog()"/>
    !!]
    allocate(normalizerSequence     )
    allocate(normalizerSequence%next)
    normalizerSequence     %normalizer_ => outputAnalysisDistributionNormalizerBinWidth_
    normalizerSequence%next%normalizer_ => outputAnalysisDistributionNormalizerLog10ToLog_
    allocate(outputAnalysisDistributionNormalizer_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizer_"            constructor="outputAnalysisDistributionNormalizerSequence(normalizerSequence)"/>
    !!]
    ! Compute the number of buffer bins to add to either side of the luminosity function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    bufferCount=max(int(bufferWidth/log10(luminosities(2)/luminosities(1)))+1,bufferCountMinimum)
    ! Construct the object.
    outputAnalysisTargetData_=outputAnalysisTargetDataStandard(                                                                                                              &
         &                                                     xAxisLabel      =var_str('$L_{\mathrm{H}\alpha}$ [ergs/s]'                                                 ), &
         &                                                     yAxisLabel      =var_str('$\mathrm{d}n/\mathrm{d}\log_\mathrm{e} L_{\mathrm{H}\alpha}$ [$_\chi$Mpc$^{-3}$]'), &
         &                                                     xAxisIsLog      =.true.                                                                                     , &
         &                                                     yAxisIsLog      =.true.                                                                                     , &
         &                                                     targetLabel     =targetLabel                                                                                , &
         &                                                     valueTarget     =functionValueTarget                                                                        , &
         &                                                     covarianceTarget=functionCovarianceTarget                                                                     &
         &                                                    )
    self%outputAnalysisVolumeFunction1D=                                                         &
         & outputAnalysisVolumeFunction1D(                                                       &
         &                                'luminosityFunctionHalpha'//label                    , &
         &                                comment                                              , &
         &                                var_str('luminosity'                                ), &
         &                                var_str('Hα luminosity at the bin center'           ), &
         &                                var_str('ergs/s'                                    ), &
         &                                var_str('erg/s'                                     ), &
         &                                .false.                                              , &
         &                                ergs                                                 , &
         &                                var_str('luminosityFunction'                        ), &
         &                                var_str('luminosity function averaged over each bin'), &
         &                                var_str('ᵪMpc⁻³'                                    ), &
         &                                var_str('Mpc^-3'                                    ), &
         &                                .true.                                               , &
         &                                megaParsec**(-3)                                     , &
         &                                log10(luminosities)                                  , &
         &                                bufferCount                                          , &
         &                                outputWeight                                         , &
         &                                nodePropertyExtractor_                               , &
         &                                outputAnalysisPropertyOperatorSequence_              , &
         &                                outputAnalysisPropertyOperatorAntiLog10_             , &
         &                                outputAnalysisWeightOperator_                        , &
         &                                outputAnalysisDistributionOperator_                  , &
         &                                outputAnalysisDistributionNormalizer_                , &
         &                                galacticFilter_                                      , &
         &                                outputTimes_                                         , &
         &                                outputAnalysisCovarianceModelBinomial                , &
         &                                covarianceBinomialBinsPerDecade                      , &
         &                                covarianceBinomialMassHaloMinimum                    , &
         &                                covarianceBinomialMassHaloMaximum                    , &
         &                                .false.                                              , &
         &                                outputAnalysisTargetData_                              &
         &                               )
    ! Clean up.
    !![
    <objectDestructor name="nodePropertyExtractorAttenuated_"                />
    <objectDestructor name="nodePropertyExtractor_"                          />
    <objectDestructor name="outputAnalysisPropertyOperatorLog10_"            />
    <objectDestructor name="outputAnalysisPropertyOperatorAntiLog10_"        />
    <objectDestructor name="outputAnalysisPropertyOperatorSequence_"         />
    <objectDestructor name="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_"/>
    <objectDestructor name="outputAnalysisDistributionNormalizer_"           />
    <objectDestructor name="outputAnalysisWeightOperator_"                   />
    <objectDestructor name="outputAnalysisDistributionNormalizerBinWidth_"   />
    <objectDestructor name="outputAnalysisDistributionNormalizerLog10ToLog_" />
    !!]
    nullify(propertyOperatorSequence)
    nullify(normalizerSequence      )
    nullify(extractors              )
    nullify(extractor_              )
    return
  end function luminosityFunctionHalphaConstructorInternal

  subroutine luminosityFunctionHalphaDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`outputAnalysisLuminosityFunctionHalpha` output analysis class.
    !!}
    type(outputAnalysisLuminosityFunctionHalpha), intent(inout) :: self

    !![
    <objectDestructor name="self%surveyGeometry_"              />
    <objectDestructor name="self%dustAttenuation_"             />
    <objectDestructor name="self%cosmologyFunctions_"          />
    <objectDestructor name="self%cosmologyFunctionsData"       />
    <objectDestructor name="self%starFormationHistory_"        />
    <objectDestructor name="self%hiiRegionLuminosityFunction_" />
    <objectDestructor name="self%hiiRegionMassFunction_"       />
    <objectDestructor name="self%hiiRegionDensityDistribution_"/>
    <objectDestructor name="self%hiiRegionEscapeFraction_"     />
    !!]
    return
  end subroutine luminosityFunctionHalphaDestructor

