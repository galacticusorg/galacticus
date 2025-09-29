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
  Implements a spin parameter distribution output analysis class.
  !!}
  
  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass

  !![
  <outputAnalysis name="outputAnalysisSpinDistribution">
   <description>A stellar mass function output analysis class.</description>
   <runTimeFileDependencies paths="fileName"/>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisSpinDistribution
     !!{
     A spinDistribution output analysis class.
     !!}
     private
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_          => null()
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_           => null()
     class           (nbodyHaloMassErrorClass          ), pointer :: nbodyHaloMassError_           => null()
     class           (haloMassFunctionClass            ), pointer :: haloMassFunction_             => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_          => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     class           (virialDensityContrastClass       ), pointer :: virialDensityContrast_        => null(), virialDensityContrastDefinition_   => null()
     double precision                                             :: timeRecent                             , logNormalRange                              , &
          &                                                          massMinimum                            , massMaximum                                 , &
          &                                                          spinMinimum                            , spinMaximum                                 , &
          &                                                          countSpinsPerDecade                    , redshift                                    , &
          &                                                          time                                   , energyEstimateParticleCountMaximum          , &
          &                                                          massParticle
     integer         (c_size_t                         )          :: countSpins
     integer                                                      :: particleCountMinimum
     logical                                                      :: errorTolerant
     type            (varying_string                   )          :: fileName
   contains
     final     ::                  spinDistributionDestructor
     procedure :: logLikelihood => spinDistributionLogLikelihood
  end type outputAnalysisSpinDistribution

  interface outputAnalysisSpinDistribution
     !!{
     Constructors for the \refClass{outputAnalysisSpinDistribution} output analysis class.
     !!}
     module procedure spinDistributionConstructorParameters
       module procedure spinDistributionConstructorFile
         module procedure spinDistributionConstructorInternal
        end interface outputAnalysisSpinDistribution

contains

  function spinDistributionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisSpinDistribution} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisSpinDistribution   )                              :: self
    type            (inputParameters                  ), intent(inout)               :: parameters
    class           (cosmologyParametersClass         ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                 ), pointer                     :: outputTimes_
    class           (nbodyHaloMassErrorClass          ), pointer                     :: nbodyHaloMassError_
    class           (haloMassFunctionClass            ), pointer                     :: haloMassFunction_
    class           (darkMatterHaloScaleClass         ), pointer                     :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass), pointer                     :: darkMatterProfileScaleRadius_
    class           (virialDensityContrastClass       ), pointer                     :: virialDensityContrast_       , virialDensityContrastDefinition_
    double precision                                   , dimension(:  ), allocatable :: functionValueTarget          , functionCovarianceTarget1D
    double precision                                   , dimension(:,:), allocatable :: functionCovarianceTarget
    double precision                                                                 :: timeRecent                   , logNormalRange                    , &
         &                                                                              massMinimum                  , massMaximum                       , &
         &                                                                              spinMinimum                  , spinMaximum                       , &
         &                                                                              countSpinsPerDecade          , redshift                          , &
         &                                                                              time                         , energyEstimateParticleCountMaximum, &
         &                                                                              massParticle
    integer         (c_size_t                         )                              :: countSpins
    integer                                                                          :: particleCountMinimum
    logical                                                                          :: errorTolerant
    type            (varying_string                   )                              :: targetLabel                  , fileName                          , &
         &                                                                              label                        , comment
    
    !![
    <inputParameter>
      <name>errorTolerant</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>Error tolerance for the N-body spin distribution operator.</description>
    </inputParameter>
    <inputParameter>
      <name>logNormalRange</name>
      <source>parameters</source>
      <defaultValue>100.0d0</defaultValue>
      <defaultSource>Approximately the range expected for the \cite{bett_spin_2007} ``QE'' cut.</defaultSource>
      <description>The multiplicative range of the log-normal distribution used to model the distribution of the mass and energy terms in the spin parameter. Specifically, the lognormal distribution is truncated outside the range $(\lambda_\mathrm{m}/R,\lambda_\mathrm{m} R$, where $\lambda_\mathrm{m}$ is the measured spin, and $R=${\normalfont \ttfamily [logNormalRange]}</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"          name="cosmologyParameters_"             source="parameters"                                                />
    <objectBuilder class="outputTimes"                  name="outputTimes_"                     source="parameters"                                                />
    <objectBuilder class="nbodyHaloMassError"           name="nbodyHaloMassError_"              source="parameters"                                                />
    <objectBuilder class="haloMassFunction"             name="haloMassFunction_"                source="parameters"                                                />
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"             source="parameters"                                                />
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_"    source="parameters"                                                />
    <objectBuilder class="virialDensityContrast"        name="virialDensityContrast_"           source="parameters"                                                />
    <objectBuilder class="virialDensityContrast"        name="virialDensityContrastDefinition_" source="parameters" parameterName="virialDensityContrastDefinition"/>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"              source="parameters"                                                />
    !!]
    if (parameters%isPresent('fileName')) then
       !![
       <inputParameter>
         <name>fileName</name>
         <source>parameters</source>
         <description>The name of the file from which to read spin distribution function parameters.</description>
       </inputParameter>
       <inputParameter>
         <name>comment</name>
         <source>parameters</source>
         <description>A comment describing this analysis.</description>
       </inputParameter>
       <inputParameter>
         <name>label</name>
         <source>parameters</source>
         <description>A label for this analysis.</description>
       </inputParameter>
       !!]
       self=outputAnalysisSpinDistribution(char(fileName),label,comment,logNormalRange,errorTolerant,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfileScaleRadius_,outputTimes_,virialDensityContrast_,virialDensityContrastDefinition_)
    else
       !![
       <inputParameter>
         <name>label</name>
         <source>parameters</source>
         <variable>label</variable>
         <description>A label for the spin distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>comment</name>
         <source>parameters</source>
         <variable>comment</variable>
         <description>A descriptive comment for the spin distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>redshift</name>
         <source>parameters</source>
         <description>The redshift at which to compute the spin distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>massMinimum</name>
         <source>parameters</source>
         <description>Minimum halo mass for the spin distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>massMaximum</name>
         <source>parameters</source>
         <description>Maximum halo mass for the spin distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>spinMinimum</name>
         <source>parameters</source>
         <description>Minimum spin for the spin distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>spinMaximum</name>
         <source>parameters</source>
         <description>Maximum spin for the spin distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>countSpinsPerDecade</name>
         <source>parameters</source>
         <description>Number of spins per decade at which to compute the spin distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>timeRecent</name>
         <source>parameters</source>
         <description>Halos which experienced a major node merger within a time $\Delta t=${\normalfont \ttfamily [timeRecent]} of the analysis time will be excluded from the analysis.</description>
       </inputParameter>
       <inputParameter>
         <name>particleCountMinimum</name>
         <source>parameters</source>
         <description>The minimum particle count to assume when computing N-body errors on spins.</description>
       </inputParameter>
       <inputParameter>
         <name>massParticle</name>
         <source>parameters</source>
         <description>The mass of the particle used in the N-body simulation from which spins were measured.</description>
       </inputParameter>
       <inputParameter>
         <name>energyEstimateParticleCountMaximum</name>
         <source>parameters</source>
         <description>The maximum number of particles used in estimating halo energies when measuring spins from the N-body simulation.</description>
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
       time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
       countSpins=int(log10(spinMaximum/spinMinimum)*countSpinsPerDecade,kind=c_size_t)+1_c_size_t
       !![
       <conditionalCall>
        <call>
         self=outputAnalysisSpinDistribution(                                    &amp;
          &amp;                              label                             , &amp;
          &amp;                              comment                           , &amp;
          &amp;                              time                              , &amp;
          &amp;                              massMinimum                       , &amp;
          &amp;                              massMaximum                       , &amp;
          &amp;                              spinMinimum                       , &amp;
          &amp;                              spinMaximum                       , &amp;
          &amp;                              countSpins                        , &amp;
          &amp;                              timeRecent                        , &amp;
          &amp;                              massParticle                      , &amp;
          &amp;                              particleCountMinimum              , &amp;
          &amp;                              energyEstimateParticleCountMaximum, &amp;
          &amp;                              logNormalRange                    , &amp;
          &amp;                              errorTolerant                     , &amp;
          &amp;                              cosmologyParameters_              , &amp;
          &amp;                              cosmologyFunctions_               , &amp;
          &amp;                              nbodyHaloMassError_               , &amp;
          &amp;                              haloMassFunction_                 , &amp;
          &amp;                              darkMatterHaloScale_              , &amp;
          &amp;                              darkMatterProfileScaleRadius_     , &amp;
          &amp;                              outputTimes_                      , &amp;
          &amp;                              virialDensityContrast_            , &amp;
          &amp;                              virialDensityContrastDefinition_    &amp;
          &amp;                              {conditions}                        &amp;
          &amp;                             )
        </call>
        <argument name="targetLabel"              value="targetLabel"              parameterPresent="parameters"/>
        <argument name="functionValueTarget"      value="functionValueTarget"      parameterPresent="parameters"/>
        <argument name="functionCovarianceTarget" value="functionCovarianceTarget" parameterPresent="parameters"/>
       </conditionalCall>
       !!]
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"                    />
    <objectDestructor name="nbodyHaloMassError_"             />
    <objectDestructor name="haloMassFunction_"               />
    <objectDestructor name="darkMatterHaloScale_"            />
    <objectDestructor name="darkMatterProfileScaleRadius_"   />
    <objectDestructor name="virialDensityContrast_"          />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    <objectDestructor name="cosmologyFunctions_"             />
    <objectDestructor name="cosmologyParameters_"            />
    !!]
    return
  end function spinDistributionConstructorParameters

  function spinDistributionConstructorFile(fileName,label,comment,logNormalRange,errorTolerant,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfileScaleRadius_,outputTimes_,virialDensityContrast_,virialDensityContrastDefinition_) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisSpinDistribution} output analysis class which reads all required properties from file.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctionsClass
    use :: IO_HDF5                          , only : hdf5Object
    use :: HDF5_Access                      , only : hdf5Access
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast          , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisSpinDistribution   )                                :: self
    character       (len=*                            ), intent(in   )                 :: fileName
    type            (varying_string                   ), intent(in   )                 :: label                        , comment
    logical                                            , intent(in   )                 :: errorTolerant
    double precision                                   , intent(in   )                 :: logNormalRange
    class           (outputTimesClass                 ), intent(inout)                 :: outputTimes_
    class           (cosmologyParametersClass         ), intent(in   )                 :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), intent(inout)                 :: cosmologyFunctions_
    class           (nbodyHaloMassErrorClass          ), intent(in   )                 :: nbodyHaloMassError_
    class           (haloMassFunctionClass            ), intent(in   )                 :: haloMassFunction_
    class           (darkMatterHaloScaleClass         ), intent(in   )                 :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass), intent(in   )                 :: darkMatterProfileScaleRadius_
    class           (virialDensityContrastClass       ), intent(in   )                 :: virialDensityContrast_       , virialDensityContrastDefinition_
    double precision                                   , allocatable  , dimension(:  ) :: functionValueTarget          , spin
    integer         (c_size_t                         ), allocatable  , dimension(:  ) :: functionCountTarget
    double precision                                   , allocatable  , dimension(:,:) :: functionCovarianceTarget
    type            (varying_string                   )                                :: targetLabel
    type            (hdf5Object                       )                                :: dataFile                     , simulationGroup                   , &
         &                                                                                attributesGroup
    integer                                                                            :: i                            , particleCountMinimum
    double precision                                                                   :: massParticle                 , redshift                          , &
         &                                                                                timeRecent                   , massMinimum                       , &
         &                                                                                massMaximum                  , energyEstimateParticleCountMaximum, &
         &                                                                                time

    !$ call hdf5Access%set  ()
    dataFile       =hdf5Object(fileName,readOnly=.true.)
    simulationGroup=dataFile       %openGroup('simulation0001')
    attributesGroup=simulationGroup%openGroup('simulation'    )
    call simulationGroup%readDataset  ('spin'                              ,spin                              )
    call simulationGroup%readDataset  ('spinDistributionFunction'          ,functionValueTarget               )
    call simulationGroup%readDataset  ('count'                             ,functionCountTarget               )
    call attributesGroup%readAttribute('labelTarget'                       ,targetLabel                       )
    call attributesGroup%readAttribute('massMinimum'                       ,massMinimum                       )
    call attributesGroup%readAttribute('massMaximum'                       ,massMaximum                       )
    call attributesGroup%readAttribute('massParticle'                      ,massParticle                      )
    call attributesGroup%readAttribute('redshift'                          ,redshift                          )
    call attributesGroup%readAttribute('timeRecent'                        ,timeRecent                        )
    call attributesGroup%readAttribute('particleCountMinimum'              ,particleCountMinimum              )
    call attributesGroup%readAttribute('energyEstimateParticleCountMaximum',energyEstimateParticleCountMaximum)
    !$ call hdf5Access%unset()
    ! Compute a (diagonal) covariance matrix from the counts.
    allocate(functionCovarianceTarget(size(functionValueTarget),size(functionValueTarget)))
    functionCovarianceTarget=0.0d0
    do i=1,size(functionValueTarget)
       if (functionCountTarget(i) > 0_c_size_t)                                                    &
            & functionCovarianceTarget(i,i)=functionValueTarget(i)**2/dble(functionCountTarget(i))
    end do
    ! Convert redshift to time.
    time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    ! Build the object.
    self=outputAnalysisSpinDistribution(label,comment,time,massMinimum,massMaximum,spin(1),spin(size(spin)),size(spin,kind=c_size_t),timeRecent,massParticle,particleCountMinimum,energyEstimateParticleCountMaximum,logNormalRange,errorTolerant,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfileScaleRadius_,outputTimes_,virialDensityContrast_,virialDensityContrastDefinition_,targetLabel,functionValueTarget,functionCovarianceTarget)
    return
  end function spinDistributionConstructorFile

  function spinDistributionConstructorInternal(label,comment,time,massMinimum,massMaximum,spinMinimum,spinMaximum,countSpins,timeRecent,massParticle,particleCountMinimum,energyEstimateParticleCountMaximum,logNormalRange,errorTolerant,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfileScaleRadius_,outputTimes_,virialDensityContrast_,virialDensityContrastDefinition_,targetLabel,functionValueTarget,functionCovarianceTarget) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisSpinDistribution} output analysis class.
    !!}
    use :: Dark_Matter_Halo_Scales                 , only : darkMatterHaloScaleClass
    use :: Galactic_Filters                        , only : filterList                                       , galacticFilterAll                           , galacticFilterHaloIsolated                    , galacticFilterHaloMassRange                 , &
          &                                                 galacticFilterNodeMajorMergerRecent              , galacticFilterNot
    use :: Error                                   , only : Error_Report
    use :: Halo_Mass_Functions                     , only : haloMassFunctionClass
    use :: Halo_Spin_Distributions                 , only : haloSpinDistributionDeltaFunction                , haloSpinDistributionNbodyErrors
    use :: ISO_Varying_String                      , only : var_str
    use :: Node_Property_Extractors                , only : nodePropertyExtractorSpin
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Numerical_Ranges                        , only : Make_Range                                       , rangeTypeLogarithmic
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Normalizers, only : normalizerList                                   , outputAnalysisDistributionNormalizerBinWidth, outputAnalysisDistributionNormalizerLog10ToLog, outputAnalysisDistributionNormalizerSequence, &
          &                                                 outputAnalysisDistributionNormalizerUnitarity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorSpinNBodyErrors
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10          , outputAnalysisPropertyOperatorLog10
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                            , only : outputTimesClass
    use :: Statistics_NBody_Halo_Mass_Errors       , only : nbodyHaloMassErrorClass
    implicit none
    type            (outputAnalysisSpinDistribution                   )                                             :: self
    type            (varying_string                                   )                             , intent(in   ) :: label                                            , comment
    double precision                                                                                , intent(in   ) :: timeRecent                                       , logNormalRange                          , &
         &                                                                                                             massMinimum                                      , massMaximum                             , &
         &                                                                                                             spinMinimum                                      , spinMaximum                             , &
         &                                                                                                             energyEstimateParticleCountMaximum               , time                                    , &
         &                                                                                                             massParticle
    integer         (c_size_t                                         )                             , intent(in   ) :: countSpins
    integer                                                                                         , intent(in   ) :: particleCountMinimum
    logical                                                                                         , intent(in   ) :: errorTolerant
    class           (cosmologyParametersClass                         ), target                     , intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass                          ), target                     , intent(in   ) :: cosmologyFunctions_
    class           (outputTimesClass                                 ), target                     , intent(inout) :: outputTimes_
    class           (nbodyHaloMassErrorClass                          ), target                     , intent(in   ) :: nbodyHaloMassError_
    class           (haloMassFunctionClass                            ), target                     , intent(in   ) :: haloMassFunction_
    class           (darkMatterHaloScaleClass                         ), target                     , intent(in   ) :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass                ), target                     , intent(in   ) :: darkMatterProfileScaleRadius_
    class           (virialDensityContrastClass                       ), target                     , intent(in   ) :: virialDensityContrast_                           , virialDensityContrastDefinition_
    type            (varying_string                                   ), optional                   , intent(in   ) :: targetLabel
    double precision                                                   , optional   , dimension(:  ), intent(in   ) :: functionValueTarget
    double precision                                                   , optional   , dimension(:,:), intent(in   ) :: functionCovarianceTarget
    type            (nodePropertyExtractorSpin                        ), pointer                                    :: nodePropertyExtractor_
    type            (outputAnalysisPropertyOperatorLog10              ), pointer                                    :: outputAnalysisPropertyOperator_
    type            (outputAnalysisWeightOperatorIdentity             ), pointer                                    :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerSequence     ), pointer                                    :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionNormalizerUnitarity    ), pointer                                    :: outputAnalysisDistributionNormalizerUnitarity_
    type            (outputAnalysisDistributionNormalizerBinWidth     ), pointer                                    :: outputAnalysisDistributionNormalizerBinWidth_
    type            (outputAnalysisDistributionNormalizerLog10ToLog   ), pointer                                    :: outputAnalysisDistributionNormalizerLog10ToLog_
    type            (outputAnalysisPropertyOperatorAntiLog10          ), pointer                                    :: outputAnalysisPropertyOperatorAntiLog10_
    type            (outputAnalysisDistributionOperatorSpinNBodyErrors), pointer                                    :: outputAnalysisDistributionOperator_
    type            (normalizerList                                   ), pointer                                    :: normalizerSequence
    type            (galacticFilterHaloIsolated                       ), pointer                                    :: galacticFilterHaloIsolated_
    type            (galacticFilterNodeMajorMergerRecent              ), pointer                                    :: galacticFilterNodeMajorMergerRecent_
    type            (galacticFilterHaloMassRange                      ), pointer                                    :: galacticFilterHaloMassRange_
    type            (galacticFilterNot                                ), pointer                                    :: galacticFilterNot_
    type            (galacticFilterAll                                ), pointer                                    :: galacticFilterAll_
    type            (filterList                                       ), pointer                                    :: filters_                                         , filter_
    type            (haloSpinDistributionNbodyErrors                  ), pointer                                    :: haloSpinDistribution_
    type            (haloSpinDistributionDeltaFunction                ), pointer                                    :: haloSpinDistributionDeltaFunction_
    double precision                                                   , allocatable, dimension(:  )                :: spins
    double precision                                                   , allocatable, dimension(:,:)                :: outputWeight
    integer         (c_size_t                                         ), parameter                                  :: bufferCountMinimum                      =5
    double precision                                                   , parameter                                  :: bufferWidthLogarithmic                  =0.500d00
    integer                                                            , parameter                                  :: covarianceBinomialBinsPerDecade         =2
    double precision                                                   , parameter                                  :: covarianceBinomialMassHaloMinimum       =3.000d11, covarianceBinomialMassHaloMaximum=1.0d15
    integer         (c_size_t                                         )                                             :: i                                                , bufferCount
    !![
    <constructorAssign variables="label, comment, time, massMinimum, massMaximum, spinMinimum, spinMaximum, countSpins, timeRecent, massParticle, particleCountMinimum, energyEstimateParticleCountMaximum, logNormalRange, errorTolerant, *cosmologyParameters_, *cosmologyFunctions_, *nbodyHaloMassError_, *haloMassFunction_, *darkMatterHaloScale_, *darkMatterProfileScaleRadius_, *outputTimes_, *virialDensityContrast_, *virialDensityContrastDefinition_, targetLabel, functionValueTarget, functionCovarianceTarget"/>
    !!]
    
    ! Build grid of spins.
    allocate(spins(countSpins))
    spins=Make_Range(spinMinimum,spinMaximum,int(countSpins),rangeType=rangeTypeLogarithmic)
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(countSpins,outputTimes_%count()))
    outputWeight=0.0d0
    do i=1,outputTimes_%count()
       if (Values_Agree(outputTimes_%time(i),time,absTol=1.0d-10)) outputWeight(:,i)=1.0d0
    end do
    if (any(sum(outputWeight,dim=2) /= 1.0d0)) call Error_Report('required output time is not available'//{introspection:location})
    ! Build an N-body halo spin distribution class.
    allocate(haloSpinDistributionDeltaFunction_)
    allocate(haloSpinDistribution_             )
    !![
    <referenceConstruct object="haloSpinDistributionDeltaFunction_">
     <constructor>
      haloSpinDistributionDeltaFunction(                                                                       &amp;
        &amp;                           spin                              =0.0d0                               &amp;
        &amp;                          )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="haloSpinDistribution_">
     <constructor>
      haloSpinDistributionNbodyErrors  (                                                                       &amp;
        &amp;                                                              haloSpinDistributionDeltaFunction_, &amp;
        &amp;                           massParticle                      =massParticle                      , &amp;
        &amp;                           particleCountMinimum              =particleCountMinimum              , &amp;
        &amp;                           energyEstimateParticleCountMaximum=energyEstimateParticleCountMaximum, &amp;
        &amp;                           logNormalRange                    =logNormalRange                    , &amp;
        &amp;                           time                              =time                              , &amp;
        &amp;                           nbodyHaloMassError_               =nbodyHaloMassError_               , &amp;
        &amp;                           cosmologyFunctions_               =cosmologyFunctions_               , &amp;
        &amp;                           haloMassFunction_                 =haloMassFunction_                 , &amp;
        &amp;                           darkMatterHaloScale_              =darkMatterHaloScale_              , &amp;
        &amp;                           darkMatterProfileScaleRadius_     =darkMatterProfileScaleRadius_       &amp;
        &amp;                          )
     </constructor>
    </referenceConstruct>
    !!]
    ! Create a spin parameter property extractor.
    allocate(nodePropertyExtractor_        )
    !![
    <referenceConstruct object="nodePropertyExtractor_"                   constructor="nodePropertyExtractorSpin                        (darkMatterHaloScale_                                                   )"/>
    !!]
    ! Create a log10 property operator.
    allocate(outputAnalysisPropertyOperator_         )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"          constructor="outputAnalysisPropertyOperatorLog10              (                                                                       )"/>
    !!]
    ! Create an identity weight operator.
    allocate(outputAnalysisWeightOperator_           )
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_"            constructor="outputAnalysisWeightOperatorIdentity             (                                                                       )"/>
    !!]
    ! Create an N-body spin error distribution operator.
    allocate(outputAnalysisDistributionOperator_     )
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_"      constructor="outputAnalysisDistributionOperatorSpinNBodyErrors(errorTolerant                                   ,haloSpinDistribution_ )"/>
    !!]
    ! Create anti-log10 operator.
    allocate(outputAnalysisPropertyOperatorAntiLog10_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorAntiLog10_" constructor="outputAnalysisPropertyOperatorAntiLog10          (                                                                       )"/>
    !!]
    ! Create a filter to select isolated halos with no recent major merger.
    allocate(galacticFilterHaloIsolated_             )
    !![
    <referenceConstruct object="galacticFilterHaloIsolated_"              constructor="galacticFilterHaloIsolated                       (                                                                       )"/>
    !!]
    allocate(galacticFilterHaloMassRange_            )
    !![
    <referenceConstruct object="galacticFilterHaloMassRange_"             constructor="galacticFilterHaloMassRange                      (massMinimum                         ,massMaximum,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_)"/>
    !!]
    allocate(galacticFilterNodeMajorMergerRecent_    )
    !![
    <referenceConstruct object="galacticFilterNodeMajorMergerRecent_"     constructor="galacticFilterNodeMajorMergerRecent              (timeRecent                                                             )"/>
    !!]
    allocate(galacticFilterNot_                      )
    !![
    <referenceConstruct object="galacticFilterNot_"                       constructor="galacticFilterNot                                (galacticFilterNodeMajorMergerRecent_                                   )"/>
    !!]
    allocate(filters_                                )
    filter_ => filters_
    filter_%filter_ => galacticFilterHaloIsolated_
    allocate(filter_%next)
    filter_ => filter_%next
    filter_%filter_ => galacticFilterHaloMassRange_
    allocate(filter_%next)
    filter_ => filter_%next
    filter_%filter_ => galacticFilterNot_
    allocate(galacticFilterAll_                      )
    !![
    <referenceConstruct object="galacticFilterAll_"                       constructor="galacticFilterAll                                (filters_                                                               )"/>
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
         &                                                             var_str('spinDistribution')//label                 , &
         &                                                             comment                                            , &
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
         &                                                             var_str('$\mathrm{d}p/\mathrm{d}\log\lambda$'     ), &
         &                                                             .true.                                             , &
         &                                                             .true.                                             , &
         &                                                             targetLabel                                        , &
         &                                                             functionValueTarget                                , &
         &                                                             functionCovarianceTarget                             &
         &                                                            )
    !![
    <objectDestructor name="haloSpinDistributionDeltaFunction_"             />
    <objectDestructor name="haloSpinDistribution_"                          />
    <objectDestructor name="nodePropertyExtractor_"                         />
    <objectDestructor name="outputAnalysisPropertyOperator_"                />
    <objectDestructor name="outputAnalysisWeightOperator_"                  />
    <objectDestructor name="outputAnalysisDistributionOperator_"            />
    <objectDestructor name="outputAnalysisPropertyOperatorAntiLog10_"       />
    <objectDestructor name="galacticFilterHaloIsolated_"                    />
    <objectDestructor name="galacticFilterHaloMassRange_"                   />
    <objectDestructor name="galacticFilterNodeMajorMergerRecent_"           />
    <objectDestructor name="galacticFilterNot_"                             />
    <objectDestructor name="galacticFilterAll_"                             />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"          />
    <objectDestructor name="outputAnalysisDistributionNormalizerUnitarity_" />
    <objectDestructor name="outputAnalysisDistributionNormalizerBinWidth_"  />
    <objectDestructor name="outputAnalysisDistributionNormalizerLog10ToLog_"/>
    !!]
    nullify(filters_          )
    nullify(filter_           )
    nullify(normalizerSequence)
    return
  end function spinDistributionConstructorInternal

  double precision function spinDistributionLogLikelihood(self)
    !!{
    Return the log-likelihood of the spin distribution function.
    !!}
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use            :: Error                       , only : Error_Report
    use            :: Interface_GSL               , only : GSL_Success
    use            :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisSpinDistribution), intent(inout)                 :: self
    double precision                                , allocatable  , dimension(:,:) :: functionCovarianceCombined
    double precision                                , allocatable  , dimension(:  ) :: functionValueDifference
    logical                                         , allocatable  , dimension(:  ) :: mask
    type            (vector                        )                                :: residual
    type            (matrix                        )                                :: covariance
    integer         (c_size_t                      )                                :: i                         , j , &
         &                                                                             ii                        , jj
    integer                                                                         :: status
    
    ! Check for existence of a target distribution.
    if (allocated(self%functionValueTarget)) then
       ! Finalize analysis.
       call self%finalizeAnalysis()
       ! Find bins which have a measured target value.
       mask=self%functionValueTarget > 0.0d0
       if (count(mask) > 0) then
          ! Allocate workspaces.
          allocate(functionCovarianceCombined(count(mask),count(mask)))
          allocate(functionValueDifference   (count(mask)            ))
          ! Find combined covariance and difference between model and target.
          ii=0
          do i=1,self%binCount
             if (mask(i)) then
                ii=ii+1
                functionValueDifference(ii)=+self%functionValue      (i) &
                     &                      -self%functionValueTarget(i)
                jj=0
                do j=1,self%binCount
                   if (mask(j)) then
                      jj=jj+1
                      ! Compute total covariance.
                      functionCovarianceCombined(ii,jj)=+self%functionCovarianceTarget(i,j) &
                           &                            +self%functionCovariance      (i,j)
                   end if
                end do
             end if
          end do
          residual  =vector(functionValueDifference   )
          covariance=matrix(functionCovarianceCombined)
          ! Compute the log-likelihood.
          spinDistributionLogLikelihood=-0.5d0*covariance%covarianceProduct(residual,status)
          if (status /= GSL_Success) spinDistributionLogLikelihood=logImprobable
       else
          spinDistributionLogLikelihood=+0.0d0
       end if
    else
       spinDistributionLogLikelihood   =+0.0d0
       call Error_Report('no target distribution was provided for likelihood calculation'//{introspection:location})
    end if
    return
  end function spinDistributionLogLikelihood

  subroutine spinDistributionDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisSpinDistribution} output analysis class.
    !!}
    implicit none
    type(outputAnalysisSpinDistribution), intent(inout) :: self

    !![
    <objectDestructor name="self%nbodyHaloMassError_"             />
    <objectDestructor name="self%haloMassFunction_"               />
    <objectDestructor name="self%darkMatterHaloScale_"            />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"   />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%cosmologyParameters_"            />
    !!]
    return
  end subroutine spinDistributionDestructor
