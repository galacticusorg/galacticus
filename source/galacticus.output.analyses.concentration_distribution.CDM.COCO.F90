!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  Contains a module which implements a concentration distribution output analysis class for COCO CDM data.
  !!}

  !![
  <outputAnalysis name="outputAnalysisConcentrationDistributionCDMCOCO">
   <description>A concentration distribution function output analysis class for COCO CDM data.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisConcentrationDistributionCDMCOCO
     !!{
     A concentration distribution output analysis class for COCO CDM data.
     !!}
     private
     integer :: distributionNumber
   contains
  end type outputAnalysisConcentrationDistributionCDMCOCO

  interface outputAnalysisConcentrationDistributionCDMCOCO
     !!{
     Constructors for the ``concentrationDistributionCDMCOCO'' output analysis class.
     !!}
     module procedure concentrationDistributionCDMCOCOConstructorParameters
     module procedure concentrationDistributionCDMCOCOConstructorInternal
  end interface outputAnalysisConcentrationDistributionCDMCOCO

contains

  function concentrationDistributionCDMCOCOConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``concentrationDistributionCDMCOCO'' output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctions , cosmologyFunctionsClass
    use :: Cosmology_Parameters             , only : cosmologyParameters, cosmologyParametersClass
    use :: Input_Parameters                 , only : inputParameter     , inputParameters
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassError , nbodyHaloMassErrorClass
    implicit none
    type            (outputAnalysisConcentrationDistributionCDMCOCO)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (cosmologyParametersClass                      ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass                              ), pointer       :: outputTimes_
    class           (nbodyHaloMassErrorClass                       ), pointer       :: nbodyHaloMassError_
    integer                                                                         :: distributionNumber
    double precision                                                                :: formationTimeRecent

    !![
    <inputParameter>
      <name>distributionNumber</name>
      <source>parameters</source>
      <description>The number (1-7) of the distribution to compute.</description>
    </inputParameter>
    <inputParameter>
      <name>formationTimeRecent</name>
      <source>parameters</source>
      <description>Halos which ``formed'' more recently than this time in the past will be excluded from the analysis.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="outputTimes"         name="outputTimes_"         source="parameters"/>
    <objectBuilder class="nbodyHaloMassError"  name="nbodyHaloMassError_"  source="parameters"/>
    !!]
    self=outputAnalysisConcentrationDistributionCDMCOCO(distributionNumber,formationTimeRecent,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,outputTimes_)
    !![
    <inputParametersValidate source="parameters" />
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="outputTimes_"        />
    <objectDestructor name="nbodyHaloMassError_" />
    !!]
    return
  end function concentrationDistributionCDMCOCOConstructorParameters

  function concentrationDistributionCDMCOCOConstructorInternal(distributionNumber,formationTimeRecent,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,outputTimes_) result(self)
    !!{
    Internal constructor for the ``concentrationDistributionCDMCOCO'' output analysis class.
    !!}
    use :: Cosmology_Functions                     , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters                    , only : cosmologyParametersClass
    use :: Galactic_Filters                        , only : filterList                                      , galacticFilterAll                           , galacticFilterFormationTime                   , galacticFilterHaloIsolated
    use :: Galacticus_Error                        , only : Galacticus_Error_Report
    use :: Galacticus_Paths                        , only : galacticusPath                                  , pathTypeDataStatic
    use :: HDF5_Access                             , only : hdf5Access
    use :: IO_HDF5                                 , only : hdf5Object
    use :: ISO_Varying_String                      , only : var_str
    use :: Memory_Management                       , only : allocateArray
    use :: Node_Property_Extractors                , only : nodePropertyExtractorConcentration              , nodePropertyExtractorMassHalo
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Normalizers, only : normalizerList                                  , outputAnalysisDistributionNormalizerBinWidth, outputAnalysisDistributionNormalizerLog10ToLog, outputAnalysisDistributionNormalizerSequence, &
          &                                                 outputAnalysisDistributionNormalizerUnitarity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorRndmErrNbdyCnc
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10         , outputAnalysisPropertyOperatorIdentity      , outputAnalysisPropertyOperatorLog10
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorNbodyMass
    use :: Output_Times                            , only : outputTimesClass
    use :: Statistics_NBody_Halo_Mass_Errors       , only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast                 , only : fixedDensityTypeCritical                        , virialDensityContrastFixed
    implicit none
    type            (outputAnalysisConcentrationDistributionCDMCOCO  )                                :: self
    class           (cosmologyParametersClass                        ), target     , intent(in   )    :: cosmologyParameters_
    class           (cosmologyFunctionsClass                         ), target     , intent(in   )    :: cosmologyFunctions_
    class           (outputTimesClass                                ), target     , intent(inout)    :: outputTimes_
    class           (nbodyHaloMassErrorClass                         ), target     , intent(in   )    :: nbodyHaloMassError_
    integer                                                                        , intent(in   )    :: distributionNumber
    double precision                                                               , intent(in   )    :: formationTimeRecent
    type            (nodePropertyExtractorConcentration              ), pointer                       :: nodePropertyExtractor_
    type            (outputAnalysisPropertyOperatorLog10             ), pointer                       :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10         ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorNbodyMass           ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorIdentity          ), pointer                       :: outputAnalysisWeightPropertyOperator_
    type            (nodePropertyExtractorMassHalo                   ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisDistributionNormalizerSequence    ), pointer                       :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionOperatorRndmErrNbdyCnc), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisDistributionNormalizerUnitarity   ), pointer                       :: outputAnalysisDistributionNormalizerUnitarity_
    type            (outputAnalysisDistributionNormalizerBinWidth    ), pointer                       :: outputAnalysisDistributionNormalizerBinWidth_
    type            (outputAnalysisDistributionNormalizerLog10ToLog  ), pointer                       :: outputAnalysisDistributionNormalizerLog10ToLog_
    type            (normalizerList                                  ), pointer                       :: normalizer_
    type            (galacticFilterHaloIsolated                      ), pointer                       :: galacticFilterHaloIsolated_
    type            (galacticFilterFormationTime                     ), pointer                       :: galacticFilterFormationTime_
    type            (galacticFilterAll                               ), pointer                       :: galacticFilter_
    type            (filterList                                      ), pointer                       :: filters_
    type            (virialDensityContrastFixed                      ), pointer                       :: virialDensityContrast_
    double precision                                                  , allocatable, dimension(:    ) :: concentrations                                      , masses                                  , &
         &                                                                                               functionTarget
    double precision                                                  , allocatable, dimension(:,:  ) :: outputWeight                                        , functionTargets                         , &
         &                                                                                               functionCovarianceTarget
    double precision                                                  , allocatable, dimension(:,:,:) :: functionCovarianceTargets
    double precision                                                  , parameter                     :: massParticle                            = +1.612d+05
    double precision                                                  , parameter                     :: haloDensityContrast                     = +2.000d+02
    double precision                                                  , parameter   , dimension(3)    :: concentrationFitA                       =[                                                      &
         &                                                                                                                                         -0.270d+00,                                           &
         &                                                                                                                                         +1.780d+00,                                           &
         &                                                                                                                                         -0.490d+00                                            &
         &                                                                                                                                        ]
    double precision                                                  , parameter                     :: concentrationFitB                       = -0.550d+00
    integer                                                           , parameter                     :: covarianceBinomialBinsPerDecade         =  2
    double precision                                                  , parameter                     :: covarianceBinomialMassHaloMinimum       = +3.000d+11, covarianceBinomialMassHaloMaximum=1.0d15
    integer         (c_size_t                                        )                                :: iOutput                                             , bufferCount
    type            (hdf5Object                                      )                                :: dataFile
    double precision                                                                                  :: massMinimum                                         , massMaximum
    character       (len=16                                          )                                :: distributionName
    character       (len= 5                                          )                                :: massMinimumLabel                                    , massMaximumLabel
    !![
    <constructorAssign variables="distributionNumber"/>
    !!]

    ! Validate input.
    if (distributionNumber < 1 .or. distributionNumber > 7) call Galacticus_Error_Report('distributionNumber ∈ [1..7] is required'//{introspection:location})
    !$ call hdf5Access%set()
    call dataFile%openFile   (char(galacticusPath(pathTypeDataStatic)//'darkMatter/concentrationDistributionCocoCDM.hdf5'),readOnly=.true.                   )
    call dataFile%readDataset(                                         'concentration'                                    ,         concentrations           )
    call dataFile%readDataset(                                         'mass'                                             ,         masses                   )
    call dataFile%readDataset(                                         'distribution'                                     ,         functionTargets          )
    call dataFile%readDataset(                                         'distributionCovariance'                           ,         functionCovarianceTargets)
    call dataFile%close      (                                                                                                                               )
    !$ call hdf5Access%unset()
    allocate(functionTarget          (size(concentrations)                     ))
    allocate(functionCovarianceTarget(size(concentrations),size(concentrations)))
    functionTarget          =functionTargets          (  :,distributionNumber)
    functionCovarianceTarget=functionCovarianceTargets(:,:,distributionNumber)
    deallocate(functionTargets          )
    deallocate(functionCovarianceTargets)
    self%binCount=size(concentrations)
    ! Determine minimum and maximum halo masses for this distribution.
    massMinimum=masses(distributionNumber)/sqrt(masses(2)/masses(1))
    massMaximum=masses(distributionNumber)*sqrt(masses(2)/masses(1))
    write (massMinimumLabel,'(f5.2)') log10(massMinimum)
    write (massMaximumLabel,'(f5.2)') log10(massMaximum)
    ! Compute weights that apply to each output redshift.
    call allocateArray(outputWeight,[self%binCount,outputTimes_%count()])
    outputWeight=0.0d0
    do iOutput=1,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(iOutput),0.0d0,absTol=1.0d-10)) outputWeight(:,iOutput)=1.0d0
    end do
    ! Build a filter which selects isolated halos, and rejects halos which formed too recently.
    allocate(galacticFilter_                  )
    allocate(galacticFilterHaloIsolated_      )
    allocate(galacticFilterFormationTime_     )
    allocate(filters_                         )
    allocate(filters_                    %next)
    filters_     %filter_ => galacticFilterHaloIsolated_
    filters_%next%filter_ => galacticFilterFormationTime_
    !![
    <referenceConstruct object="galacticFilterHaloIsolated_"  constructor="galacticFilterHaloIsolated  (                   )"/>
    <referenceConstruct object="galacticFilterFormationTime_" constructor="galacticFilterFormationTime (formationTimeRecent)"/>
    <referenceConstruct object="galacticFilter_"              constructor="galacticFilterAll           (filters_           )"/>
    !!]
    ! Create a distribution normalizer which normalizes to bin width and unitarity.
    allocate(normalizer_          )
    allocate(normalizer_%next     )
    allocate(normalizer_%next%next)
    allocate(outputAnalysisDistributionNormalizerUnitarity_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerUnitarity_"  constructor="outputAnalysisDistributionNormalizerUnitarity ()"/>
    !!]
    allocate(outputAnalysisDistributionNormalizerBinWidth_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerBinWidth_"   constructor="outputAnalysisDistributionNormalizerBinWidth  ()"/>
    !!]
    allocate(outputAnalysisDistributionNormalizerLog10ToLog_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerLog10ToLog_" constructor="outputAnalysisDistributionNormalizerLog10ToLog()"/>
    !!]
    normalizer_          %normalizer_ => outputAnalysisDistributionNormalizerUnitarity_
    normalizer_%next     %normalizer_ => outputAnalysisDistributionNormalizerBinWidth_
    normalizer_%next%next%normalizer_ => outputAnalysisDistributionNormalizerLog10ToLog_
    allocate(outputAnalysisDistributionNormalizer_ )
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizer_">
     <constructor>
      outputAnalysisDistributionNormalizerSequence    (                                        &amp;
          &amp;                                        normalizer_                             &amp;
          &amp;                                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build log10() property operator.
    allocate(outputAnalysisPropertyOperator_       )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_">
     <constructor>
      outputAnalysisPropertyOperatorLog10             (                                        &amp;
        &amp;                                         )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_     )
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_">
     <constructor>
      outputAnalysisPropertyOperatorAntiLog10         (                                        &amp;
        &amp;                                         )
     </constructor>
    </referenceConstruct>
    !!]
    ! Create a virial density contrast object matched to the definition used by Ludlow et al. (2016).
    allocate(virialDensityContrast_                )
    !![
    <referenceConstruct object="virialDensityContrast_">
     <constructor>
      virialDensityContrastFixed                      (                                        &amp;
         &amp;                                         haloDensityContrast                   , &amp;
         &amp;                                         fixedDensityTypeCritical              , &amp;
         &amp;                                         2.0d0                                 , &amp;
         &amp;                                         cosmologyParameters_                  , &amp;
         &amp;                                         cosmologyFunctions_                     &amp;
         &amp;                                        )
     </constructor>
    </referenceConstruct>
    !!]
    ! Create a concentration property extractor.
    allocate(nodePropertyExtractor_      )
    !![
    <referenceConstruct object="nodePropertyExtractor_">
     <constructor>
      nodePropertyExtractorConcentration              (                                        &amp;
          &amp;                                        virialDensityContrast_                  &amp;
          &amp;                                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Create a halo mass property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_">
     <constructor>
      nodePropertyExtractorMassHalo                   (                                        &amp;
         &amp;                                         virialDensityContrast_                  &amp;
         &amp;                                        )
     </constructor>
    </referenceConstruct>
    !!]
    ! Create an identity property operator.
    allocate(outputAnalysisWeightPropertyOperator_ )
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyOperator_">
     <constructor>
      outputAnalysisPropertyOperatorIdentity          (                                        &amp;
        &amp;                                         )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build error distribution operator.
    allocate(outputAnalysisDistributionOperator_   )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
      outputAnalysisDistributionOperatorRndmErrNbdyCnc(                                        &amp;
         &amp;                                         concentrationFitA                     , &amp;
         &amp;                                         concentrationFitB                     , &amp;
         &amp;                                         massParticle                          , &amp;
         &amp;                                         outputAnalysisWeightPropertyExtractor_  &amp;
         &amp;                                        )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build N-body mass distribution weight operator.
    allocate(outputAnalysisWeightOperator_         )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_">
     <constructor>
      outputAnalysisWeightOperatorNbodyMass          (                                        &amp;
         &amp;                                        massMinimum                           , &amp;
         &amp;                                        massMaximum                           , &amp;
         &amp;                                        outputAnalysisWeightPropertyExtractor_, &amp;
         &amp;                                        outputAnalysisWeightPropertyOperator_ , &amp;
         &amp;                                        nbodyHaloMassError_                     &amp;
         &amp;                                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Determine number of buffer bins.
    bufferCount=0
    ! Construct the object.
    write (distributionName,'(i2.2)') distributionNumber
    self%outputAnalysisVolumeFunction1D=                                                                                                                                               &
         & outputAnalysisVolumeFunction1D(                                                                                                                                             &
         &                                var_str('concentrationDistributionCDMCOCO')//trim(distributionName)                                                                        , &
         &                                var_str('Concentration distribution for $')//massMinimumLabel//' \le \log_{10} M_\mathrm{200c}/\mathrm{M}_\odot < '//massMaximumLabel//'$' , &
         &                                var_str('concentration'                                                                                                                   ), &
         &                                var_str('Concentration at the bin center'                                                                                                 ), &
         &                                var_str('dimensionless'                                                                                                                   ), &
         &                                0.0d0                                                                                                                                      , &
         &                                var_str('concentrationFunction'                                                                                                           ), &
         &                                var_str('Concentration distribution averaged over each bin'                                                                               ), &
         &                                var_str('dimensionless'                                                                                                                   ), &
         &                                0.0d0                                                                                                                                      , &
         &                                log10(concentrations)                                                                                                                      , &
         &                                bufferCount                                                                                                                                , &
         &                                outputWeight                                                                                                                               , &
         &                                nodePropertyExtractor_                                                                                                                     , &
         &                                outputAnalysisPropertyOperator_                                                                                                            , &
         &                                outputAnalysisPropertyUnoperator_                                                                                                          , &
         &                                outputAnalysisWeightOperator_                                                                                                              , &
         &                                outputAnalysisDistributionOperator_                                                                                                        , &
         &                                outputAnalysisDistributionNormalizer_                                                                                                      , &
         &                                galacticFilter_                                                                                                                            , &
         &                                outputTimes_                                                                                                                               , &
         &                                outputAnalysisCovarianceModelPoisson                                                                                                       , &
         &                                covarianceBinomialBinsPerDecade                                                                                                            , &
         &                                covarianceBinomialMassHaloMinimum                                                                                                          , &
         &                                covarianceBinomialMassHaloMaximum                                                                                                          , &
         &                                .false.                                                                                                                                    , &
         &                                var_str('$c_\mathrm{200c}$'                                                                                                               ), &
         &                                var_str('$\mathrm{d}p/\mathrm{d}\log_{10}c_\mathrm{200c}$'                                                                                ), &
         &                                .true.                                                                                                                                     , &
         &                                .false.                                                                                                                                    , &
         &                                var_str('Benson et al. (2019)'                                                                                                            ), &
         &                                functionTarget                                                                                                                             , &
         &                                functionCovarianceTarget                                                                                                                     &
         &                               )
    !![
    <objectDestructor name="galacticFilterHaloIsolated_"                    />
    <objectDestructor name="galacticFilterFormationTime_"                   />
    <objectDestructor name="galacticFilter_"                                />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"          />
    <objectDestructor name="outputAnalysisPropertyOperator_"                />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"              />
    <objectDestructor name="virialDensityContrast_"                         />
    <objectDestructor name="nodePropertyExtractor_"               />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"         />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"          />
    <objectDestructor name="outputAnalysisDistributionOperator_"            />
    <objectDestructor name="outputAnalysisWeightOperator_"                  />
    <objectDestructor name="outputAnalysisDistributionNormalizerUnitarity_" />
    <objectDestructor name="outputAnalysisDistributionNormalizerBinWidth_"  />
    <objectDestructor name="outputAnalysisDistributionNormalizerLog10ToLog_"/>
    !!]
    nullify(normalizer_)
    nullify(filters_   )
  return
  end function concentrationDistributionCDMCOCOConstructorInternal
