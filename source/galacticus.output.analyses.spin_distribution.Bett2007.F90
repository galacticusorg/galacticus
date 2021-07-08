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
Contains a module which implements a spin parameter distribution output analysis class.
!!}

  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass

  !![
  <outputAnalysis name="outputAnalysisSpinDistributionBett2007">
   <description>A stellar mass function output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisSpinDistributionBett2007
     !!{
     A spinDistributionBett2007 output analysis class.
     !!}
     private
  end type outputAnalysisSpinDistributionBett2007

  interface outputAnalysisSpinDistributionBett2007
     !!{
     Constructors for the ``spinDistributionBett2007'' output analysis class.
     !!}
     module procedure spinDistributionBett2007ConstructorParameters
     module procedure spinDistributionBett2007ConstructorInternal
  end interface outputAnalysisSpinDistributionBett2007

contains

  function spinDistributionBett2007ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``spinDistributionBett2007'' output analysis class which takes a parameter set as input.
    !!}
    use :: Functions_Global, only : Virial_Density_Contrast_Percolation_Objects_Constructor_
    use :: Input_Parameters, only : inputParameter                                          , inputParameters
    implicit none
    type            (outputAnalysisSpinDistributionBett2007)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass                      ), pointer       :: outputTimes_
    class           (nbodyHaloMassErrorClass               ), pointer       :: nbodyHaloMassError_
    class           (haloMassFunctionClass                 ), pointer       :: haloMassFunction_
    class           (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass             ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterProfileScaleRadiusClass     ), pointer       :: darkMatterProfileScaleRadius_
    class           (*                                     ), pointer       :: percolationObjects_
    double precision                                                        :: timeRecent                   , logNormalRange
    logical                                                                 :: errorTolerant

    !![
    <inputParameter>
      <name>errorTolerant</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>Error tolerance for the N-body spin distribution operator.</description>
    </inputParameter>
    <inputParameter>
      <name>timeRecent</name>
      <source>parameters</source>
      <description>Halos which experienced a major node merger within a time $\Delta t=${\normalfont \ttfamily [timeRecent]} of the analysis time will be excluded from the analysis.</description>
    </inputParameter>
    <inputParameter>
      <name>logNormalRange</name>
      <source>parameters</source>
      <defaultValue>100.0d0</defaultValue>
      <defaultSource>Approximately the range expected for the \cite{bett_spin_2007} ``QE'' cut.</defaultSource>
      <description>The multiplicative range of the log-normal distribution used to model the distribution of the mass and energy terms in the spin parameter. Specifically, the lognormal distribution is truncated outside the range $(\lambda_\mathrm{m}/R,\lambda_\mathrm{m} R$, where $\lambda_\mathrm{m}$ is the measured spin, and $R=${\normalfont \ttfamily [logNormalRange]}</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="outputTimes"                  name="outputTimes_"                  source="parameters"/>
    <objectBuilder class="nbodyHaloMassError"           name="nbodyHaloMassError_"           source="parameters"/>
    <objectBuilder class="haloMassFunction"             name="haloMassFunction_"             source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    percolationObjects_ => Virial_Density_Contrast_Percolation_Objects_Constructor_(parameters)
    self                =  outputAnalysisSpinDistributionBett2007(timeRecent,logNormalRange,errorTolerant,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfileDMO_,darkMatterProfileScaleRadius_,outputTimes_,percolationObjects_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="outputTimes_"                 />
    <objectDestructor name="nbodyHaloMassError_"          />
    <objectDestructor name="haloMassFunction_"            />
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function spinDistributionBett2007ConstructorParameters

  function spinDistributionBett2007ConstructorInternal(timeRecent,logNormalRange,errorTolerant,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfileDMO_,darkMatterProfileScaleRadius_,outputTimes_,percolationObjects_) result(self)
    !!{
    Internal constructor for the ``spinDistributionBett2007'' output analysis class.
    !!}
    use :: Cosmology_Functions                     , only : cosmologyFunctionsClass
    use :: Dark_Matter_Halo_Scales                 , only : darkMatterHaloScaleClass
    use :: Galactic_Filters                        , only : filterList                                       , galacticFilterAll                           , galacticFilterHaloIsolated                    , galacticFilterHaloMass                      , &
          &                                                 galacticFilterNodeMajorMergerRecent              , galacticFilterNot
    use :: Galacticus_Error                        , only : Galacticus_Error_Report
    use :: Galacticus_Paths                        , only : galacticusPath                                   , pathTypeDataStatic
    use :: Halo_Mass_Functions                     , only : haloMassFunctionClass
    use :: Halo_Spin_Distributions                 , only : haloSpinDistributionDeltaFunction                , haloSpinDistributionNbodyErrors
    use :: IO_HDF5                                 , only : hdf5Access                                       , hdf5Object
    use :: ISO_Varying_String                      , only : var_str
    use :: Memory_Management                       , only : allocateArray
    use :: Node_Property_Extractors                , only : nodePropertyExtractorSpin
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Normalizers, only : normalizerList                                   , outputAnalysisDistributionNormalizerBinWidth, outputAnalysisDistributionNormalizerLog10ToLog, outputAnalysisDistributionNormalizerSequence, &
          &                                                 outputAnalysisDistributionNormalizerUnitarity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorSpinNBodyErrors
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10          , outputAnalysisPropertyOperatorLog10
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                            , only : outputTimesClass
    use :: Statistics_NBody_Halo_Mass_Errors       , only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast                 , only : virialDensityContrastPercolation
    implicit none
    type            (outputAnalysisSpinDistributionBett2007           )                              :: self
    double precision                                                                , intent(in   )  :: timeRecent                                       , logNormalRange
    logical                                                                         , intent(in   )  :: errorTolerant
    class           (cosmologyFunctionsClass                          ), target     , intent(inout)  :: cosmologyFunctions_
    class           (outputTimesClass                                 ), target     , intent(inout)  :: outputTimes_
    class           (nbodyHaloMassErrorClass                          ), target     , intent(in   )  :: nbodyHaloMassError_
    class           (haloMassFunctionClass                            ), target     , intent(in   )  :: haloMassFunction_
    class           (darkMatterHaloScaleClass                         ), target     , intent(in   )  :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                        ), target     , intent(in   )  :: darkMatterProfileDMO_
    class           (darkMatterProfileScaleRadiusClass                ), target     , intent(in   )  :: darkMatterProfileScaleRadius_
    class           (*                                                ), target     , intent(in   )  :: percolationObjects_
    type            (nodePropertyExtractorSpin                        ), pointer                     :: nodePropertyExtractor_
    type            (outputAnalysisPropertyOperatorLog10              ), pointer                     :: outputAnalysisPropertyOperator_
    type            (outputAnalysisWeightOperatorIdentity             ), pointer                     :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerSequence     ), pointer                     :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionNormalizerUnitarity    ), pointer                     :: outputAnalysisDistributionNormalizerUnitarity_
    type            (outputAnalysisDistributionNormalizerBinWidth     ), pointer                     :: outputAnalysisDistributionNormalizerBinWidth_
    type            (outputAnalysisDistributionNormalizerLog10ToLog   ), pointer                     :: outputAnalysisDistributionNormalizerLog10ToLog_
    type            (outputAnalysisPropertyOperatorAntiLog10          ), pointer                     :: outputAnalysisPropertyOperatorAntiLog10_
    type            (outputAnalysisDistributionOperatorSpinNBodyErrors), pointer                     :: outputAnalysisDistributionOperator_
    type            (normalizerList                                   ), pointer                     :: normalizerSequence
    type            (virialDensityContrastPercolation                 ), pointer                     :: virialDensityContrast_
    type            (galacticFilterHaloIsolated                       ), pointer                     :: galacticFilterHaloIsolated_
    type            (galacticFilterNodeMajorMergerRecent              ), pointer                     :: galacticFilterNodeMajorMergerRecent_
    type            (galacticFilterHaloMass                           ), pointer                     :: galacticFilterHaloMass_
    type            (galacticFilterNot                                ), pointer                     :: galacticFilterNot_
    type            (galacticFilterAll                                ), pointer                     :: galacticFilterAll_
    type            (filterList                                       ), pointer                     :: filters_                                         , filter_
    type            (haloSpinDistributionNbodyErrors                  ), pointer                     :: haloSpinDistribution_
    type            (haloSpinDistributionDeltaFunction                ), pointer                     :: haloSpinDistributionDeltaFunction_
    double precision                                                   , allocatable, dimension(:  ) :: spins                                            , functionTarget                          , &
         &                                                                                              functionErrorTarget
    double precision                                                   , allocatable, dimension(:,:) :: outputWeight                                     , functionCovarianceTarget
    double precision                                                   , parameter                   :: massParticleMillennium                  =1.178d09
    double precision                                                   , parameter                   :: countParticleMininum                    =3.000d02
    integer         (c_size_t                                         ), parameter                   :: bufferCountMinimum                      =5
    double precision                                                   , parameter                   :: bufferWidthLogarithmic                  =0.500d00
    integer                                                            , parameter                   :: covarianceBinomialBinsPerDecade         =2
    double precision                                                   , parameter                   :: covarianceBinomialMassHaloMinimum       =3.000d11, covarianceBinomialMassHaloMaximum=1.0d15
    integer         (c_size_t                                         )                              :: i                                                , bufferCount
    type            (hdf5Object                                       )                              :: dataFile

    ! Construct spins matched to those used by Bett et al. (2007).
    !$ call hdf5Access%set()
    call dataFile%openFile   (char(galacticusPath(pathTypeDataStatic)//'darkMatter/bett2007HaloSpinDistribution.hdf5'),readOnly=.true.             )
    call dataFile%readDataset(                                         'spinParameter'                                ,         spins              )
    call dataFile%readDataset(                                         'distribution'                                 ,         functionTarget     )
    call dataFile%readDataset(                                         'distributionError'                            ,         functionErrorTarget)
    call dataFile%close      (                                                                                                                     )
    !$ call hdf5Access%unset()
    ! Convert target distribution from per log₁₀(λ) to per ln(λ), and construct covariance matrix.
    self%binCount      =size(spins)
    functionTarget     =functionTarget     /log(10.0d0)
    functionErrorTarget=functionErrorTarget/log(10.0d0)
    allocate(functionCovarianceTarget(self%binCount,self%binCount))
    functionCovarianceTarget=0.0d0
    do i=1,self%binCount
       functionCovarianceTarget(i,i)=functionErrorTarget(i)**2
    end do
    ! Compute weights that apply to each output redshift.
    call allocateArray(outputWeight,[self%binCount,outputTimes_%count()])
    outputWeight=0.0d0
    do i=1,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(i),0.0d0,absTol=1.0d-10)) outputWeight(:,i)=1.0d0
    end do
    if (any(sum(outputWeight,dim=2) /= 1.0d0)) call Galacticus_Error_Report('zero redshift output is required'//{introspection:location})
    ! Build an N-body halo spin distribution class.
    allocate(haloSpinDistributionDeltaFunction_)
    allocate(haloSpinDistribution_             )
    !![
    <referenceConstruct object="haloSpinDistributionDeltaFunction_">
     <constructor>
      haloSpinDistributionDeltaFunction(                                                                                              &amp;
        &amp;                           spin                              =                                                 0.0d0     &amp;
        &amp;                          )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="haloSpinDistribution_">
     <constructor>
      haloSpinDistributionNbodyErrors  (                                                                                              &amp;
        &amp;                                                              haloSpinDistributionDeltaFunction_                       , &amp;
        &amp;                           massParticle                      =massParticleMillennium                                   , &amp;
        &amp;                           particleCountMinimum              =                                               300       , &amp;
        &amp;                           energyEstimateParticleCountMaximum=                                              1000.000d0 , &amp;
        &amp;                           logNormalRange                    =logNormalRange                                           , &amp;
        &amp;                           time                              =cosmologyFunctions_               %cosmicTime(   1.000d0), &amp;
        &amp;                           nbodyHaloMassError_               =nbodyHaloMassError_                                      , &amp;
        &amp;                           haloMassFunction_                 =haloMassFunction_                                        , &amp;
        &amp;                           darkMatterHaloScale_              =darkMatterHaloScale_                                     , &amp;
        &amp;                           darkMatterProfileDMO_             =darkMatterProfileDMO_                                    , &amp;
        &amp;                           darkMatterProfileScaleRadius_     =darkMatterProfileScaleRadius_                              &amp;
        &amp;                          )
     </constructor>
    </referenceConstruct>
    !!]
    ! Create a spin parameter property extractor.
    allocate(nodePropertyExtractor_        )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                   constructor="nodePropertyExtractorSpin                       (                                                                                       )"/>
    !!]
    ! Create a log10 property operator.
    allocate(outputAnalysisPropertyOperator_         )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"          constructor="outputAnalysisPropertyOperatorLog10              (                                                                                      )"/>
    !!]
    ! Create an identity weight operator.
    allocate(outputAnalysisWeightOperator_           )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"            constructor="outputAnalysisWeightOperatorIdentity             (                                                                                      )"/>
    !!]
    ! Create an N-body spin error distribution operator.
    allocate(outputAnalysisDistributionOperator_     )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_"      constructor="outputAnalysisDistributionOperatorSpinNBodyErrors(errorTolerant                              , haloSpinDistribution_                    )"/>
    !!]
    ! Create anit-log10 operator.
    allocate(outputAnalysisPropertyOperatorAntiLog10_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorAntiLog10_" constructor="outputAnalysisPropertyOperatorAntiLog10          (                                                                                      )"/>
    !!]
    ! Create a virial density contrast class to match the friends-of-friends halo definition used by Bett et al. (2007).
    allocate(virialDensityContrast_)
    !![
    <referenceConstruct object="virialDensityContrast_"                   constructor="virialDensityContrastPercolation                 (0.2d0                                      ,cosmologyFunctions_   ,percolationObjects_)"/>
    !!]
    ! Create a filter to select isolated halos with no recent major merger.
    allocate(galacticFilterHaloIsolated_             )
    !![
    <referenceConstruct object="galacticFilterHaloIsolated_"              constructor="galacticFilterHaloIsolated                       (                                                                                      )"/>
    !!]
    allocate(galacticFilterHaloMass_                )
    !![
    <referenceConstruct object="galacticFilterHaloMass_"                  constructor="galacticFilterHaloMass                           (countParticleMininum*massParticleMillennium,virialDensityContrast_                    )"/>
    !!]
    allocate(galacticFilterNodeMajorMergerRecent_    )
    !![
    <referenceConstruct object="galacticFilterNodeMajorMergerRecent_"     constructor="galacticFilterNodeMajorMergerRecent              (timeRecent                                                                            )"/>
    !!]
    allocate(galacticFilterNot_                      )
    !![
    <referenceConstruct object="galacticFilterNot_"                       constructor="galacticFilterNot                                (galacticFilterNodeMajorMergerRecent_                                                  )"/>
    !!]
    allocate(filters_                                )
    filter_ => filters_
    filter_%filter_ => galacticFilterHaloIsolated_
    allocate(filter_%next)
    filter_ => filter_%next
    filter_%filter_ => galacticFilterHaloMass_
    allocate(filter_%next)
    filter_ => filter_%next
    filter_%filter_ => galacticFilterNot_
    allocate(galacticFilterAll_                      )
    !![
    <referenceConstruct object="galacticFilterAll_"                       constructor="galacticFilterAll                                (filters_                                                                              )"/>
    !!]
    ! Create a distribution normalizer which normalizes to unit integral, and then to bin width.
    allocate(outputAnalysisDistributionNormalizerUnitarity_ )
    allocate(outputAnalysisDistributionNormalizerBinWidth_  )
    allocate(outputAnalysisDistributionNormalizerLog10ToLog_)
    allocate(outputAnalysisDistributionNormalizer_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerUnitarity_"  constructor="outputAnalysisDistributionNormalizerUnitarity ()"/>
    <referenceConstruct object="outputAnalysisDistributionNormalizerBinWidth_"   constructor="outputAnalysisDistributionNormalizerBinWidth  ()"/>
    <referenceConstruct object="outputAnalysisDistributionNormalizerLog10ToLog_" constructor="outputAnalysisDistributionNormalizerLog10ToLog()"/>
    !!]
    allocate(normalizerSequence          )
    allocate(normalizerSequence%next     )
    allocate(normalizerSequence%next%next)
    normalizerSequence          %normalizer_ => outputAnalysisDistributionNormalizerUnitarity_
    normalizerSequence%next     %normalizer_ => outputAnalysisDistributionNormalizerBinWidth_
    normalizerSequence%next%next%normalizer_ => outputAnalysisDistributionNormalizerLog10ToLog_
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizer_" constructor="outputAnalysisDistributionNormalizerSequence(normalizerSequence)"/>
    !!]
    ! Compute the number of buffer bins to add to either side of the mass function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    bufferCount=max(int(bufferWidthLogarithmic/log10(spins(2)/spins(1)))+1,bufferCountMinimum)
    ! Construct the object. We convert spins to log10(spins) here.
    self%outputAnalysisVolumeFunction1D=outputAnalysisVolumeFunction1D(                                                     &
         &                                                             var_str('spinDistributionBett2007')                , &
         &                                                             var_str('Distribution of halo spin parameters'    ), &
         &                                                             var_str('spin'                                    ), &
         &                                                             var_str('Spin at the bin center'                  ), &
         &                                                             var_str('dimensionless'                           ), &
         &                                                             0.0d0                                              , &
         &                                                             var_str('spinDistributionFunction'                ), &
         &                                                             var_str('Spin distribution averaged over each bin'), &
         &                                                             var_str('dimensionless'                           ), &
         &                                                             0.0d0                                              , &
         &                                                             log10(spins)                                       , &
         &                                                             bufferCount                                        , &
         &                                                             outputWeight                                       , &
         &                                                             nodePropertyExtractor_                             , &
         &                                                             outputAnalysisPropertyOperator_                    , &
         &                                                             outputAnalysisPropertyOperatorAntiLog10_           , &
         &                                                             outputAnalysisWeightOperator_                      , &
         &                                                             outputAnalysisDistributionOperator_                , &
         &                                                             outputAnalysisDistributionNormalizer_              , &
         &                                                             galacticFilterAll_                                 , &
         &                                                             outputTimes_                                       , &
         &                                                             outputAnalysisCovarianceModelPoisson               , &
         &                                                             covarianceBinomialBinsPerDecade                    , &
         &                                                             covarianceBinomialMassHaloMinimum                  , &
         &                                                             covarianceBinomialMassHaloMaximum                  , &
         &                                                             .false.                                            , &
         &                                                             var_str('$\lambda$'                               ), &
         &                                                             var_str('$\mathrm{d}p/\mathrm{d}\ln\lambda$'      ), &
         &                                                             .true.                                             , &
         &                                                             .true.                                             , &
         &                                                             var_str('Bett et al. (2007)'                      ), &
         &                                                             functionTarget                                     , &
         &                                                             functionCovarianceTarget                             &
         &                                                            )
    !![
    <objectDestructor name="haloSpinDistributionDeltaFunction_"             />
    <objectDestructor name="haloSpinDistribution_"                          />
    <objectDestructor name="nodePropertyExtractor_"               />
    <objectDestructor name="outputAnalysisPropertyOperator_"                />
    <objectDestructor name="outputAnalysisWeightOperator_"                  />
    <objectDestructor name="outputAnalysisDistributionOperator_"            />
    <objectDestructor name="outputAnalysisPropertyOperatorAntiLog10_"       />
    <objectDestructor name="galacticFilterHaloIsolated_"                    />
    <objectDestructor name="galacticFilterHaloMass_"                        />
    <objectDestructor name="galacticFilterNodeMajorMergerRecent_"           />
    <objectDestructor name="galacticFilterNot_"                             />
    <objectDestructor name="galacticFilterAll_"                             />
    <objectDestructor name="virialDensityContrast_"                         />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"          />
    <objectDestructor name="outputAnalysisDistributionNormalizerUnitarity_" />
    <objectDestructor name="outputAnalysisDistributionNormalizerBinWidth_"  />
    <objectDestructor name="outputAnalysisDistributionNormalizerLog10ToLog_"/>
    !!]
    nullify(filters_          )
    nullify(filter_           )
    nullify(normalizerSequence)
    return
  end function spinDistributionBett2007ConstructorInternal
