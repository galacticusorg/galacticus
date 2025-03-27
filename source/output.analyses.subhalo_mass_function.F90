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
  Implements an output analysis class that computes subhalo mass functions.
  !!}
  
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <outputAnalysis name="outputAnalysisSubhaloMassFunction">
   <description>An output analysis class for subhalo mass functions.</description>
   <deepCopy>
    <functionClass variables="volumeFunctionsSubHalos, volumeFunctionsHostHalos"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="volumeFunctionsSubHalos, volumeFunctionsHostHalos"/>
   </stateStorable>
   <runTimeFileDependencies paths="fileName"/>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisSubhaloMassFunction
     !!{
     An output analysis class for subhalo mass functions.
     !!}
     private
     class           (cosmologyFunctionsClass       ), pointer                     :: cosmologyFunctions_               => null()
     class           (outputTimesClass              ), pointer                     :: outputTimes_                      => null()
     class           (virialDensityContrastClass    ), pointer                     :: virialDensityContrast_            => null(), virialDensityContrastDefinition_ => null()
     class           (cosmologyParametersClass      ), pointer                     :: cosmologyParameters_              => null()
     class           (darkMatterProfileDMOClass     ), pointer                     :: darkMatterProfileDMO_             => null()
     type            (outputAnalysisVolumeFunction1D), pointer                     :: volumeFunctionsSubHalos           => null(), volumeFunctionsHostHalos         => null()
     double precision                                , allocatable, dimension(:  ) :: massRatios                                 , massFunction                              , &
          &                                                                           massFunctionTarget
     double precision                                , allocatable, dimension(:,:) :: covariance                                 , massFunctionCovarianceTarget
     type            (varying_string                )                              :: labelTarget                                , fileName
     double precision                                                              :: negativeBinomialScatterFractional          , countFailures                             , &
          &                                                                           massRatioMinimum                           , massRatioMaximum                          , &
          &                                                                           time                                       , redshift
     integer         (c_size_t                      )                              :: countMassRatios 
     logical                                                                       :: finalized
   contains
     !![
     <methods>
       <method description="Finalize analysis." method="finalizeAnalysis" />
     </methods>
     !!]
     final     ::                     subhaloMassFunctionDestructor
     procedure :: analyze          => subhaloMassFunctionAnalyze
     procedure :: finalize         => subhaloMassFunctionFinalize
     procedure :: finalizeAnalysis => subhaloMassFunctionFinalizeAnalysis
     procedure :: reduce           => subhaloMassFunctionReduce
     procedure :: logLikelihood    => subhaloMassFunctionLogLikelihood
  end type outputAnalysisSubhaloMassFunction

  interface outputAnalysisSubhaloMassFunction
     !!{
     Constructors for the {\normalfont \ttfamily subhaloMassFunction} output analysis class.
     !!}
     module procedure subhaloMassFunctionConstructorParameters
     module procedure subhaloMassFunctionConstructorFile
     module procedure subhaloMassFunctionConstructorInternal
  end interface outputAnalysisSubhaloMassFunction

contains

  function subhaloMassFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily subhaloMassFunction} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters       , only : inputParameter            , inputParameters
    use :: Output_Times           , only : outputTimesClass
    use :: Cosmology_Functions    , only : cosmologyFunctionsClass
    use :: Virial_Density_Contrast, only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisSubhaloMassFunction)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass        ), pointer       :: darkMatterProfileDMO_
    class           (outputTimesClass                 ), pointer       :: outputTimes_
    class           (virialDensityContrastClass       ), pointer       :: virialDensityContrast_, virialDensityContrastDefinition_
    double precision                                                   :: massRatioMinimum      , massRatioMaximum                 , &
         &                                                                redshift              , negativeBinomialScatterFractional
    integer         (c_size_t                         )                :: countMassRatios
    type            (varying_string                   )                :: fileName

    if (parameters%isPresent('fileName')) then
       !![
       <inputParameter>
         <name>fileName</name>
         <source>parameters</source>
         <description>The name of the file from which to read the target dataset.</description>
       </inputParameter>
       !!]
    else
       !![
       <inputParameter>
         <name>massRatioMinimum</name>
         <source>parameters</source>
         <defaultValue>1.0d-4</defaultValue>
         <description>The minimum mass ratio to consider.</description>
       </inputParameter>
       <inputParameter>
         <name>massRatioMaximum</name>
         <source>parameters</source>
         <defaultValue>1.0d0</defaultValue>
         <description>The maximum mass ratio to consider.</description>
       </inputParameter>
       <inputParameter>
         <name>countMassRatios</name>
         <source>parameters</source>
         <defaultValue>10_c_size_t</defaultValue>
         <description>The number of bins in mass ratio to use.</description>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift at which to compute the subhalo mass function.</description>
    </inputParameter>
    <inputParameter>
      <name>negativeBinomialScatterFractional</name>
      <source>parameters</source>
      <defaultValue>0.18d0</defaultValue>
      <defaultSource>\citep{boylan-kolchin_theres_2010}</defaultSource>
      <description>The fractional scatter (relative to the Poisson scatter) in the negative binomial distribution used in likelihood calculations.</description>
    </inputParameter>
    <objectBuilder class="outputTimes"           name="outputTimes_"                     source="parameters"                                                />
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"             source="parameters"                                                />
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"              source="parameters"                                                />
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"            source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_"           source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrastDefinition_" source="parameters" parameterName="virialDensityContrastDefinition"/>
    !!]
    if (parameters%isPresent('fileName')) then
       !![
       <conditionalCall>
        <call>self=outputAnalysisSubhaloMassFunction(darkMatterProfileDMO_,outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,fileName,negativeBinomialScatterFractional{conditions})</call>
         <argument name="redshift" value="redshift" parameterPresent="parameters"/>
       </conditionalCall>
       !!]
    else
       self=outputAnalysisSubhaloMassFunction(darkMatterProfileDMO_,outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),massRatioMinimum,massRatioMaximum,countMassRatios,negativeBinomialScatterFractional)
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"                    />
    <objectDestructor name="cosmologyParameters_"            />
    <objectDestructor name="cosmologyFunctions_"             />
    <objectDestructor name="darkMatterProfileDMO_"           />
    <objectDestructor name="virialDensityContrast_"          />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function subhaloMassFunctionConstructorParameters
  
  function subhaloMassFunctionConstructorFile(darkMatterProfileDMO_,outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,fileName,negativeBinomialScatterFractional,redshift) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily subhaloMassFunction} output analysis class for internal use.
    !!}
    use :: HDF5_Access            , only : hdf5Access
    use :: IO_HDF5                , only : hdf5Object
    use :: Output_Times           , only : outputTimesClass
    use :: Cosmology_Functions    , only : cosmologyFunctionsClass
    use :: File_Utilities         , only : File_Name_Expand
    use :: Virial_Density_Contrast, only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisSubhaloMassFunction)                                :: self
    type            (varying_string                   ), intent(in   )                 :: fileName
    double precision                                   , intent(in   )                 :: negativeBinomialScatterFractional
    class           (outputTimesClass                 ), intent(inout)                 :: outputTimes_
    class           (darkMatterProfileDMOClass        ), intent(inout)                 :: darkMatterProfileDMO_
    class           (virialDensityContrastClass       ), intent(in   )                 :: virialDensityContrast_           , virialDensityContrastDefinition_
    class           (cosmologyParametersClass         ), intent(inout)                 :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), intent(inout), target         :: cosmologyFunctions_
    double precision                                   , intent(in   ), optional       :: redshift
    double precision                                   , allocatable  , dimension(:  ) :: massRatiosTarget                 , massFunctionTarget               , &
         &                                                                                massFunctionErrorTarget
    double precision                                   , allocatable  , dimension(:,:) :: massFunctionCovarianceTarget
    double precision                                                                   :: massRatioMinimum                 , massRatioMaximum                 , &
         &                                                                                time                             , redshift_
    integer         (c_size_t                         )                                :: countMassRatios                  , i
    type            (varying_string                   )                                :: labelTarget
    type            (hdf5Object                       )                                :: file                             , massFunctionGroup

    ! Read properties from the file.
    !$ call hdf5Access%set()
    call file             %openFile     (char(File_Name_Expand(char(fileName))),readOnly=.true.                 )
    call file             %readAttribute('label'                               ,         labelTarget            )
    call file             %readAttribute('redshift'                            ,         redshift_              )
    massFunctionGroup=file%openGroup('massFunction')
    call massFunctionGroup%readDataset  ('massRatio'                           ,         massRatiosTarget       )
    call massFunctionGroup%readDataset  ('massFunction'                        ,         massFunctionTarget     )
    call massFunctionGroup%readDataset  ('massFunctionError'                   ,         massFunctionErrorTarget)
    call massFunctionGroup%close        (                                                                       )
    call file             %close        (                                                                       )
    !$ call hdf5Access%unset()
    ! Override the redshift if one is provided.
    if (present(redshift)) redshift_=redshift
    ! Construct the mass function.
    time            =cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift_))
    countMassRatios =size(massRatiosTarget)
    massRatioMinimum=massRatiosTarget(              1)
    massRatioMaximum=massRatiosTarget(countMassRatios)
    allocate(massFunctionCovarianceTarget(countMassRatios,countMassRatios))
    massFunctionCovarianceTarget=0.0d0
    do i=1_c_size_t,countMassRatios
       massFunctionCovarianceTarget(i,i)=massFunctionErrorTarget(i)**2
    end do
    self=outputAnalysisSubhaloMassFunction(darkMatterProfileDMO_,outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,time,massRatioMinimum,massRatioMaximum,countMassRatios,negativeBinomialScatterFractional,massFunctionTarget,massFunctionCovarianceTarget,labelTarget)
    !![
    <constructorAssign variables="fileName, redshift"/>
    !!]
    return
  end function subhaloMassFunctionConstructorFile

  function subhaloMassFunctionConstructorInternal(darkMatterProfileDMO_,outputTimes_,virialDensityContrastDefinition_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,time,massRatioMinimum,massRatioMaximum,countMassRatios,negativeBinomialScatterFractional,massFunctionTarget,massFunctionCovarianceTarget,labelTarget) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily subhaloMassFunction} output analysis class for internal use.
    !!}
    use :: Galactic_Filters                        , only : filterList                                  , galacticFilterAll                  , galacticFilterHaloIsolated   , galacticFilterHaloNotIsolated     , &
          &                                                 galacticFilterLowPass
    use :: Node_Property_Extractors                , only : nodePropertyExtractorHostNode               , nodePropertyExtractorMassBound     , nodePropertyExtractorMassHalo, nodePropertyExtractorRadiusOrbital, &
          &                                                 nodePropertyExtractorRadiusVirial           , nodePropertyExtractorRatio
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Numerical_Ranges                        , only : Make_Range                                  , rangeTypeLinear
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorIdentity
    use :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerIdentity
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10     , outputAnalysisPropertyOperatorLog10
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorSubsampling
    use :: Output_Times                            , only : outputTimesClass
    use :: Virial_Density_Contrast                 , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisSubhaloMassFunction           )                                          :: self
    double precision                                              , intent(in   )                           :: negativeBinomialScatterFractional               , massRatioMinimum                      , &
         &                                                                                                     massRatioMaximum                                , time
    integer         (c_size_t                                    ), intent(in   )                           :: countMassRatios
    class           (outputTimesClass                            ), intent(inout), target                   :: outputTimes_
    class           (virialDensityContrastClass                  ), intent(in   ), target                   :: virialDensityContrast_                          , virialDensityContrastDefinition_
    class           (cosmologyParametersClass                    ), intent(in   ), target                   :: cosmologyParameters_
    class           (cosmologyFunctionsClass                     ), intent(in   ), target                   :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass                   ), intent(inout), target                   :: darkMatterProfileDMO_
    double precision                                              , intent(in   ), dimension(:)  , optional :: massFunctionTarget
    double precision                                              , intent(in   ), dimension(:,:), optional :: massFunctionCovarianceTarget
    type            (varying_string                              ), intent(in   )                , optional :: labelTarget
    type            (nodePropertyExtractorMassBound              )               , pointer                  :: nodePropertyExtractorMassBound_
    type            (nodePropertyExtractorHostNode               )               , pointer                  :: nodePropertyExtractorMassHost_                  , nodePropertyExtractorRadiusVirialHost_
    type            (nodePropertyExtractorMassHalo               )               , pointer                  :: nodePropertyExtractorMassHalo_
    type            (nodePropertyExtractorRadiusOrbital          )               , pointer                  :: nodePropertyExtractorRadiusOrbital_
    type            (nodePropertyExtractorRadiusVirial           )               , pointer                  :: nodePropertyExtractorRadiusVirial_
    type            (nodePropertyExtractorRatio                  )               , pointer                  :: nodePropertyExtractor_                          , nodePropertyExtractorRadiusFractional_
    type            (outputAnalysisPropertyOperatorLog10         )               , pointer                  :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10     )               , pointer                  :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorSubsampling     )               , pointer                  :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerIdentity)               , pointer                  :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionOperatorIdentity  )               , pointer                  :: outputAnalysisDistributionOperator_
    type            (galacticFilterHaloIsolated                  )               , pointer                  :: galacticFilterHosts_
    type            (galacticFilterHaloNotIsolated               )               , pointer                  :: galacticFilterIsSubhalo_
    type            (galacticFilterLowPass                       ), pointer                                 :: galacticFilterVirialRadius_
    type            (galacticFilterAll                           ), pointer                                 :: galacticFilterSubhalos_
    type            (filterList                                  ), pointer                                 :: filters_
    double precision                                              , allocatable  , dimension(:  )           :: massRatios                                      , massesHosts
    double precision                                              , allocatable  , dimension(:,:)           :: outputWeightSubhalos                            , outputWeightHosts
    integer         (c_size_t                                    ), parameter                               :: binCountHosts                        =2_c_size_t
    double precision                                              , parameter                               :: massHostLogarithmicMaximum           =1.0d2
    integer         (c_size_t                                    )                                          :: i
    !![
    <constructorAssign variables="negativeBinomialScatterFractional, countMassRatios, massRatioMinimum, massRatioMaximum, massFunctionTarget, massFunctionCovarianceTarget, labelTarget, *cosmologyFunctions_, *outputTimes_, *cosmologyParameters_, *darkMatterProfileDMO_, *virialDensityContrast_, *virialDensityContrastDefinition_"/>
    !!]

    ! Initialize.
    self%finalized=.false.
    self%redshift =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    ! Compute failure count for the negative binomial distribution used in likelihood calculations.
    self%countFailures=1.0d0/negativeBinomialScatterFractional**2
    ! Construct mass bins.
    allocate(massRatios (self%countMassRatios))
    allocate(massesHosts(                   2))
    massRatios =Make_Range(log10(self%massRatioMinimum),log10(massRatioMaximum)   ,int(self%countMassRatios),rangeTypeLinear)
    massesHosts=Make_Range(                      0.0d0 ,massHostLogarithmicMaximum,int(     binCountHosts  ),rangeTypeLinear)
    ! Create a mass ratio property extractor.
    allocate(nodePropertyExtractorMassBound_       )
    allocate(nodePropertyExtractorRadiusOrbital_   )
    allocate(nodePropertyExtractorRadiusVirial_    )
    allocate(nodePropertyExtractorMassHalo_        )
    allocate(nodePropertyExtractorMassHost_        )
    allocate(nodePropertyExtractor_                )
    allocate(nodePropertyExtractorRadiusVirialHost_)
    allocate(nodePropertyExtractorRadiusFractional_)
    !![
    <referenceConstruct object="nodePropertyExtractorMassBound_"        constructor="nodePropertyExtractorMassBound    (                                                                                                                                                   )"/>
    <referenceConstruct object="nodePropertyExtractorRadiusOrbital_"    constructor="nodePropertyExtractorRadiusOrbital(                                                                                                                                                   )"/>
    <referenceConstruct object="nodePropertyExtractorMassHalo_"         constructor="nodePropertyExtractorMassHalo     (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_                     )"/>
    <referenceConstruct object="nodePropertyExtractorRadiusVirial_"     constructor="nodePropertyExtractorRadiusVirial (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_                     )"/>
    <referenceConstruct object="nodePropertyExtractorMassHost_"         constructor="nodePropertyExtractorHostNode     (nodePropertyExtractorMassHalo_                                                                                                                     )"/>
    <referenceConstruct object="nodePropertyExtractorRadiusVirialHost_" constructor="nodePropertyExtractorHostNode     (nodePropertyExtractorRadiusVirial_                                                                                                                 )"/>
    <referenceConstruct object="nodePropertyExtractor_"                 constructor="nodePropertyExtractorRatio        ('massRatio'     ,'Ratio of subhalo to host mass'                        ,nodePropertyExtractorMassBound_    ,nodePropertyExtractorMassHost_        )"/>
    <referenceConstruct object="nodePropertyExtractorRadiusFractional_" constructor="nodePropertyExtractorRatio        ('radiusFraction','Ratio of subhalo orbital radius to host virial radius',nodePropertyExtractorRadiusOrbital_,nodePropertyExtractorRadiusVirialHost_)"/>
    !!]
    ! Create property operators and unoperators to perform conversion to/from logarithmic mass ratio.
    allocate(outputAnalysisPropertyOperator_  )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_"   constructor="outputAnalysisPropertyOperatorLog10    ()"/>
    !!]
    allocate(outputAnalysisPropertyUnoperator_)
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_" constructor="outputAnalysisPropertyOperatorAntiLog10()"/>
    !!]
    ! Create an identity weight operator.
    allocate(outputAnalysisWeightOperator_)
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_" constructor="outputAnalysisWeightOperatorSubsampling    ()"/>
    !!]
    ! Build filters which select subhalos/hosts.
    allocate(galacticFilterHosts_        )
    !![
    <referenceConstruct object="galacticFilterHosts_"        constructor="galacticFilterHaloIsolated   (                                            )"/>
    !!]
    allocate(galacticFilterIsSubhalo_    )
    !![
    <referenceConstruct object="galacticFilterIsSubhalo_"    constructor="galacticFilterHaloNotIsolated(                                            )"/>
    !!]
    allocate(galacticFilterVirialRadius_ )
    !![
    <referenceConstruct object="galacticFilterVirialRadius_" constructor="galacticFilterLowPass        (1.0d0,nodePropertyExtractorRadiusFractional_)"/>
    !!]
    allocate(galacticFilterSubhalos_     )
    allocate(filters_                    )
    allocate(filters_               %next)
    filters_     %filter_ => galacticFilterIsSubhalo_
    filters_%next%filter_ => galacticFilterVirialRadius_
    !![
    <referenceConstruct object="galacticFilterSubhalos_"     constructor="galacticFilterAll            (filters_                                    )"/>
    !!]
    ! Build an identity distribution operator.
    allocate(outputAnalysisDistributionOperator_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_" constructor="outputAnalysisDistributionOperatorIdentity()"/>
    !!]
    ! Build an identity distribution normalizers for hosts and subhalos.
    allocate(outputAnalysisDistributionNormalizer_)
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizer_" constructor="outputAnalysisDistributionNormalizerIdentity()"/>
    !!]
    ! Compute weights that apply to each output redshift.
    allocate(outputWeightSubhalos(self%countMassRatios,outputTimes_%count()))
    allocate(outputWeightHosts   (     binCountHosts  ,outputTimes_%count()))
    do i=1_c_size_t,outputTimes_%count()
       if (Values_Agree(outputTimes_%time(i),time,absTol=1.0d-6)) then
          outputWeightSubhalos(:,i)=1.0d0
          outputWeightHosts   (:,i)=1.0d0
       else
          outputWeightSubhalos(:,i)=0.0d0
          outputWeightHosts   (:,i)=0.0d0
       end if
    end do
    ! Construct the volume function 1D objects.
    allocate(self%volumeFunctionsSubHalos )
    allocate(self%volumeFunctionsHostHalos)
    !![
    <referenceConstruct isResult="yes" owner="self" object="volumeFunctionsSubHalos">
     <constructor>
      outputAnalysisVolumeFunction1D(                                                     &amp;
       &amp;                         var_str('subhaloMassFunction')                     , &amp;
       &amp;                         var_str('Subhalo mass function')                   , &amp;
       &amp;                         var_str('massBoundRatio')                          , &amp;
       &amp;                         var_str('Ratio of subhalo bound mass to host mass'), &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         1.0d0                                              , &amp;
       &amp;                         var_str('massFunction')                            , &amp;
       &amp;                         var_str('Differential subhalo mass function')      , &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         0.0d0                                              , &amp;
       &amp;                         massRatios                                         , &amp;
       &amp;                         0_c_size_t                                         , &amp;
       &amp;                         outputWeightSubhalos                               , &amp;
       &amp;                         nodePropertyExtractor_                             , &amp;
       &amp;                         outputAnalysisPropertyOperator_                    , &amp;
       &amp;                         outputAnalysisPropertyUnoperator_                  , &amp;
       &amp;                         outputAnalysisWeightOperator_                      , &amp;
       &amp;                         outputAnalysisDistributionOperator_                , &amp;
       &amp;                         outputAnalysisDistributionNormalizer_              , &amp;
       &amp;                         galacticFilterSubhalos_                            , &amp;
       &amp;                         outputTimes_                                       , &amp;
       &amp;                         outputAnalysisCovarianceModelPoisson                 &amp;
       &amp;                        )
     </constructor>
    </referenceConstruct>
    <referenceConstruct isResult="yes" owner="self" object="volumeFunctionsHostHalos">
     <constructor>
      outputAnalysisVolumeFunction1D(                                                     &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         0.0d0                                              , &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         var_str(' ')                                       , &amp;
       &amp;                         0.0d0                                              , &amp;
       &amp;                         massesHosts                                        , &amp;
       &amp;                         0_c_size_t                                         , &amp;
       &amp;                         outputWeightHosts                                  , &amp;
       &amp;                         nodePropertyExtractorMassBound_                    , &amp;
       &amp;                         outputAnalysisPropertyOperator_                    , &amp;
       &amp;                         outputAnalysisPropertyUnoperator_                  , &amp;
       &amp;                         outputAnalysisWeightOperator_                      , &amp;
       &amp;                         outputAnalysisDistributionOperator_                , &amp;
       &amp;                         outputAnalysisDistributionNormalizer_              , &amp;
       &amp;                         galacticFilterHosts_                               , &amp;
       &amp;                         outputTimes_                                       , &amp;
       &amp;                         outputAnalysisCovarianceModelPoisson                 &amp;
       &amp;                        )
     </constructor>
    </referenceConstruct>
    <objectDestructor name="nodePropertyExtractorMassBound_"       />
    <objectDestructor name="nodePropertyExtractorMassHost_"        />
    <objectDestructor name="nodePropertyExtractorMassHalo_"        />
    <objectDestructor name="nodePropertyExtractorRadiusOrbital_"   />
    <objectDestructor name="nodePropertyExtractorRadiusVirial_"    />
    <objectDestructor name="nodePropertyExtractorRadiusVirialHost_"/>
    <objectDestructor name="nodePropertyExtractorRadiusFractional_"/>
    <objectDestructor name="nodePropertyExtractor_"                />
    <objectDestructor name="outputAnalysisPropertyOperator_"       />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"     />
    <objectDestructor name="outputAnalysisWeightOperator_"         />
    <objectDestructor name="outputAnalysisDistributionOperator_"   />
    <objectDestructor name="galacticFilterHosts_"                  />
    <objectDestructor name="galacticFilterSubhalos_"               />
    <objectDestructor name="galacticFilterIsSubhalo_"              />
    <objectDestructor name="galacticFilterVirialRadius_"           />
    <objectDestructor name="outputAnalysisDistributionNormalizer_" />
    !!]
    nullify(filters_)
    return
  end function subhaloMassFunctionConstructorInternal

  subroutine subhaloMassFunctionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily subhaloMassFunction} output analysis class.
    !!}
    implicit none
    type(outputAnalysisSubhaloMassFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%volumeFunctionsSubHalos"         />
    <objectDestructor name="self%volumeFunctionsHostHalos"        />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%outputTimes_"                    />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%darkMatterProfileDMO_"           />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
   !!]
    return
  end subroutine subhaloMassFunctionDestructor

  subroutine subhaloMassFunctionAnalyze(self,node,iOutput)
    !!{
    Implement a {\normalfont \ttfamily subhaloMassFunction} output analysis.
    !!}
    implicit none
    class  (outputAnalysisSubhaloMassFunction), intent(inout) :: self
    type   (treeNode                         ), intent(inout) :: node
    integer(c_size_t                         ), intent(in   ) :: iOutput

    ! Analyze for all three volume functions.
    call self%volumeFunctionsSubHalos %analyze(node,iOutput)
    call self%volumeFunctionsHostHalos%analyze(node,iOutput)
    return
  end subroutine subhaloMassFunctionAnalyze

  subroutine subhaloMassFunctionReduce(self,reduced)
    !!{
    Implement a {\normalfont \ttfamily subhaloMassFunction} output analysis reduction.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisSubhaloMassFunction), intent(inout) :: self
    class(outputAnalysisClass              ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisSubhaloMassFunction)
       call self%volumeFunctionsSubHalos %reduce(reduced%volumeFunctionsSubHalos )
       call self%volumeFunctionsHostHalos%reduce(reduced%volumeFunctionsHostHalos)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine subhaloMassFunctionReduce

  subroutine subhaloMassFunctionFinalizeAnalysis(self)
    !!{
    Finalize analysis of a {\normalfont \ttfamily subhaloMassFunction} output analysis.
    !!}
    implicit none
    class           (outputAnalysisSubhaloMassFunction), intent(inout)               :: self
    double precision                                   , allocatable  , dimension(:) :: massFunctionHosts
    double precision                                                                 :: weight

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
    ! Retrieve results from our 1-D volume functions.
    call self%volumeFunctionsSubHalos %results(binCenter=self%massRatios,functionValue=self%massFunction     ,functionCovariance=self%covariance)
    call self%volumeFunctionsHostHalos%results(                          functionValue=     massFunctionHosts                                   )
    ! Normalize the mass function.
    weight=sum(massFunctionHosts)
    if (weight > 0.0d0) then
       self%massFunction=+self%massFunction/weight
       self%covariance  =+self%covariance  /weight**2
    else
       self%massFunction=0.0d0
       self%covariance  =0.0d0
    end if
    return
  end subroutine subhaloMassFunctionFinalizeAnalysis

  subroutine subhaloMassFunctionFinalize(self,groupName)
    !!{
    Implement a {\normalfont \ttfamily subhaloMassFunction} output analysis finalization.
    !!}
    use :: Output_HDF5                     , only : outputFile
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(outputAnalysisSubhaloMassFunction), intent(inout)           :: self
    type (varying_string                   ), intent(in   ), optional :: groupName
    type (hdf5Object                       )               , target   :: analysesGroup, subGroup
    type (hdf5Object                       )               , pointer  :: inGroup
    type (hdf5Object                       )                          :: analysisGroup, dataset

    ! Finalize analysis.
    call self%finalizeAnalysis()
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'                         )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName)                    )
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('subhaloMassFunction','Analysis of subhalo mass functions')
    call analysisGroup   %writeAttribute('Subhalo mass function'              ,'description'                                                                                         )
    call analysisGroup   %writeAttribute('function1D'                         ,'type'                                                                                                )
    call analysisGroup   %writeAttribute('$M_\mathrm{sub}/M_\mathrm{host}$'   ,'xAxisLabel'                                                                                          )
    call analysisGroup   %writeAttribute('$N(M_\mathrm{sub}/M_\mathrm{host})$','yAxisLabel'                                                                                          )
    call analysisGroup   %writeAttribute(.true.                               ,'xAxisIsLog'                                                                                          )
    call analysisGroup   %writeAttribute(.true.                               ,'yAxisIsLog'                                                                                          )
    call analysisGroup   %writeAttribute('massRatio'                          ,'xDataset'                                                                                            )
    call analysisGroup   %writeAttribute('massFunction'                       ,'yDataset'                                                                                            )
    call analysisGroup   %writeAttribute('massFunctionCovariance'             ,'yCovariance'                                                                                         )
    call analysisGroup   %writeDataset  (self%massRatios                      ,'massRatio'                   ,'Mass ratio at the bin center'                 ,datasetReturned=dataset)
    call dataset         %writeAttribute(' '                                  ,'units'                                                                                               )
    call dataset         %writeAttribute(1.0d0                                ,'unitsInSI'                                                                                           )
    call dataset         %close         (                                                                                                                                            )
    call analysisGroup   %writeDataset  (self%massFunction                    ,'massFunction'                ,'Subhalo number per bin [model]'                                       )
    call analysisGroup   %writeDataset  (self%covariance                      ,'massFunctionCovariance'      ,'Subhalo number per bin [model; covariance]'                           )
    if (allocated(self%massFunctionTarget)) then
       call analysisGroup%writeAttribute(self%logLikelihood               ()  ,'logLikelihood'                                                                                       )
       call analysisGroup%writeAttribute('massFunctionTarget'                 ,'yDatasetTarget'                                                                                      )
       call analysisGroup%writeAttribute('massFunctionCovarianceTarget'       ,'yCovarianceTarget'                                                                                   )
       call analysisGroup%writeAttribute(char(self%labelTarget)               ,'targetLabel'                                                                                         )
       call analysisGroup%writeDataset  (self%massFunctionTarget              ,'massFunctionTarget'          ,'Subhalo number per bin [observed]'                                    )
       call analysisGroup%writeDataset  (self%massFunctionCovarianceTarget    ,'massFunctionCovarianceTarget','Subhalo number per bin [observed; covariance]'                        )
    end if
    call analysisGroup   %close         (                                                                                                                                            )
    if (present(groupName)) &
         & call subGroup %close         (                                                                                                                                            )
    call analysesGroup   %close         (                                                                                                                                            )
    !$ call hdf5Access%unset()
    return
  end subroutine subhaloMassFunctionFinalize

  double precision function subhaloMassFunctionLogLikelihood(self)
    !!{
    Return the log-likelihood of a {\normalfont \ttfamily subhaloMassFunction} output analysis. The likelihood function
    assumes that the model prediction for the number of subhalos in any given mass bin follows a negative binomial
    distribution as was found for dark matter subhalos \citep[][see also
    \protect\citealt{lu_connection_2016}]{boylan-kolchin_theres_2010}. This has been confirmed by examining the results of many
    tree realizations, although it in principal could be model-dependent.
    !!}
    use :: Numerical_Constants_Math         , only : Pi
    use :: Models_Likelihoods_Constants     , only : logImpossible
    use :: Statistics_Distributions_Discrete, only : distributionFunctionDiscrete1DNegativeBinomial
    implicit none
    class           (outputAnalysisSubhaloMassFunction             ), intent(inout) :: self
    type            (distributionFunctionDiscrete1DNegativeBinomial)                :: distribution
    integer                                                                         :: i
    double precision                                                                :: negativeBinomialProbabilitySuccess, countEffective, &
         &                                                                             variance

    call self%finalizeAnalysis()
    subhaloMassFunctionLogLikelihood=0.0d0
    do i=1,size(self%massRatios)
       if (self%massFunction(i) <= 0.0d0) then
          if (nint(self%massFunctionTarget(i)) > 0) then
             subhaloMassFunctionLogLikelihood=logImpossible
             return
          end if
       else
          ! Compute the likelihood assuming a negative binomial distribution. Note that we "de-normalize" the likelihood by
          ! multiplying by √[2πσᵢ²] (the normalization term in the corresponding normal distribution). This is useful to allow
          ! (-logℒ) to be used as a metric for significant shifts in the model results, without changing the relative
          ! likelihood of models (as this de-normalization shift is a constant multiplicative factor).
          countEffective                    = dble(max(1.0d0,self%massFunctionTarget(i)))
          variance                          =+       countEffective                       &
               &                             *(                                           &
               &                               +     1.0d0                                &
               &                               +self%negativeBinomialScatterFractional**2 &
               &                               *     countEffective                       &
               &                              )
          negativeBinomialProbabilitySuccess=+  1.0d0                                                          &
               &                             /(                                                                &
               &                               +1.0d0                                                          &
               &                               +self%negativeBinomialScatterFractional**2*self%massFunction(i) &
               &                              )
          distribution                      = distributionFunctionDiscrete1DNegativeBinomial                (negativeBinomialProbabilitySuccess,     self%countFailures         )
          subhaloMassFunctionLogLikelihood  =+subhaloMassFunctionLogLikelihood                                                                                                    &
               &                             +distribution                                  %massLogarithmic(                                   nint(self%massFunctionTarget(i))) &
               &                             +0.50d0                                                                                                                              &
               &                             *log(                                                                                                                                &
               &                                  +2.0d0                                                                                                                          &
               &                                  *Pi                                                                                                                             &
               &                                  *variance                                                                                                                       &
               &                                 )
       end if
    end do
    return
  end function subhaloMassFunctionLogLikelihood
