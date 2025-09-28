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

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !!{
  Implements a concentration distribution output analysis class.
  !!}

  !![
  <outputAnalysis name="outputAnalysisConcentrationDistribution">
   <description>A concentration distribution function output analysis class.</description>
   <runTimeFileDependencies paths="fileName"/>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisConcentrationDistribution
     !!{
     A concentration distribution output analysis class.
     !!}
     private
     class           (cosmologyParametersClass  ), pointer :: cosmologyParameters_             => null()
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_              => null()
     class           (nbodyHaloMassErrorClass   ), pointer :: nbodyHaloMassError_              => null()
     class           (darkMatterProfileDMOClass ), pointer :: darkMatterProfileDMO_            => null()
     class           (virialDensityContrastClass), pointer :: virialDensityContrastDefinition_ => null(), virialDensityContrast_ => null()
     double precision                                      :: rootVarianceFractionalMinimum             , redshift                        , &
          &                                                   massMinimum                               , massMaximum                     , &
          &                                                   concentrationMinimum                      , concentrationMaximum            , &
          &                                                   countConcentrationsPerDecade              , timeRecent                      , &
          &                                                   massParticle
     type            (varying_string            )          :: fileName
   contains
     final     ::                  concentrationDistributionDestructor
     procedure :: logLikelihood => concentrationDistributionLogLikelihood
  end type outputAnalysisConcentrationDistribution

  interface outputAnalysisConcentrationDistribution
     !!{
     Constructors for the \refClass{outputAnalysisConcentrationDistribution} output analysis class.
     !!}
     module procedure concentrationDistributionConstructorParameters
     module procedure concentrationDistributionConstructorFile
     module procedure concentrationDistributionConstructorInternal
  end interface outputAnalysisConcentrationDistribution

contains

  function concentrationDistributionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisConcentrationDistribution} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctions, cosmologyFunctionsClass
    use :: Input_Parameters                 , only : inputParameter    , inputParameters
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassError, nbodyHaloMassErrorClass
    implicit none
    type            (outputAnalysisConcentrationDistribution)                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    class           (cosmologyParametersClass               ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                       ), pointer                     :: outputTimes_
    class           (nbodyHaloMassErrorClass                ), pointer                     :: nbodyHaloMassError_
    class           (darkMatterProfileDMOClass              ), pointer                     :: darkMatterProfileDMO_
    class           (virialDensityContrastClass             ), pointer                     :: virialDensityContrastDefinition_, virialDensityContrast_
    double precision                                         , dimension(:  ), allocatable :: functionValueTarget             , functionCovarianceTarget1D
    double precision                                         , dimension(:,:), allocatable :: functionCovarianceTarget
    double precision                                                                       :: time                            , timeRecent                   , &
         &                                                                                    massMinimum                     , massMaximum                  , &
         &                                                                                    concentrationMinimum            , concentrationMaximum         , &
         &                                                                                    countConcentrationsPerDecade    , redshift                     , &
         &                                                                                    massParticle                    , rootVarianceFractionalMinimum
    integer         (c_size_t                               )                              :: countConcentrations
    type            (varying_string                         )                              :: targetLabel                     , fileName                     , &
         &                                                                                    label                           , comment

    !![
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"             source="parameters"                                                />
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"              source="parameters"                                                />
    <objectBuilder class="outputTimes"           name="outputTimes_"                     source="parameters"                                                />
    <objectBuilder class="nbodyHaloMassError"    name="nbodyHaloMassError_"              source="parameters"                                                />
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"            source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_"           source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrastDefinition_" source="parameters" parameterName="virialDensityContrastDefinition"/>
    <inputParameter>
      <name>rootVarianceFractionalMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum fractional root variance (relative to the target dataset).</description>
    </inputParameter>
    !!]
    if (parameters%isPresent('fileName')) then
       !![
       <inputParameter>
         <name>fileName</name>
         <source>parameters</source>
         <description>The name of the file from which to read concentration distribution function parameters.</description>
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
       self=outputAnalysisConcentrationDistribution(char(fileName),label,comment,rootVarianceFractionalMinimum,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,virialDensityContrast_,virialDensityContrastDefinition_,outputTimes_)
    else
       !![
       <inputParameter>
         <name>label</name>
         <source>parameters</source>
         <variable>label</variable>
         <description>A label for the concentration distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>comment</name>
         <source>parameters</source>
         <variable>comment</variable>
         <description>A descriptive comment for the concentration distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>redshift</name>
         <source>parameters</source>
         <description>The redshift at which to compute the concentration distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>massMinimum</name>
         <source>parameters</source>
         <description>Minimum halo mass for the concentration distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>massMaximum</name>
         <source>parameters</source>
         <description>Maximum halo mass for the concentration distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>concentrationMinimum</name>
         <source>parameters</source>
         <description>Minimum concentration for the concentration distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>concentrationMaximum</name>
         <source>parameters</source>
         <description>Maximum concentration for the concentration distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>countConcentrationsPerDecade</name>
         <source>parameters</source>
         <description>Number of concentrations per decade at which to compute the concentration distribution function.</description>
       </inputParameter>
       <inputParameter>
         <name>timeRecent</name>
         <source>parameters</source>
         <description>Halos which experienced a major node merger within a time $\Delta t=${\normalfont \ttfamily [timeRecent]} of the analysis time will be excluded from the analysis.</description>
       </inputParameter>
       <inputParameter>
         <name>massParticle</name>
         <source>parameters</source>
         <description>The particle mass in the source N-body simulation.</description>
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
       countConcentrations=int(log10(concentrationMaximum/concentrationMinimum)*countConcentrationsPerDecade,kind=c_size_t)+1_c_size_t
       !![
       <conditionalCall>
        <call>
         self=outputAnalysisConcentrationDistribution(                                  &amp;
          &amp;                                       label                           , &amp;
          &amp;                                       comment                         , &amp;
          &amp;                                       time                            , &amp;
          &amp;                                       massMinimum                     , &amp;
          &amp;                                       massMaximum                     , &amp;
          &amp;                                       concentrationMinimum            , &amp;
          &amp;                                       concentrationMaximum            , &amp;
          &amp;                                       countConcentrations             , &amp;
          &amp;                                       timeRecent                      , &amp;
          &amp;                                       massParticle                    , &amp;
	  &amp;                                       rootVarianceFractionalMinimum   , &amp;
	  &amp;                                       darkMatterProfileDMO_           , &amp;
          &amp;                                       cosmologyParameters_            , &amp;
          &amp;                                       cosmologyFunctions_             , &amp;
          &amp;                                       nbodyHaloMassError_             , &amp;
          &amp;                                       virialDensityContrast_          , &amp;
          &amp;                                       virialDensityContrastDefinition_, &amp;
          &amp;                                       outputTimes_                      &amp;
          &amp;                                       {conditions}                      &amp;
          &amp;                                      )
        </call>
        <argument name="targetLabel"              value="targetLabel"              parameterPresent="parameters"/>
        <argument name="functionValueTarget"      value="functionValueTarget"      parameterPresent="parameters"/>
        <argument name="functionCovarianceTarget" value="functionCovarianceTarget" parameterPresent="parameters"/>
       </conditionalCall>
       !!]
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"            />
    <objectDestructor name="cosmologyFunctions_"             />
    <objectDestructor name="outputTimes_"                    />
    <objectDestructor name="nbodyHaloMassError_"             />
    <objectDestructor name="darkMatterProfileDMO_"           />
    <objectDestructor name="virialDensityContrast_"          />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function concentrationDistributionConstructorParameters

  function concentrationDistributionConstructorFile(fileName,label,comment,rootVarianceFractionalMinimum,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,virialDensityContrast_,virialDensityContrastDefinition_,outputTimes_) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisConcentrationDistribution} output analysis class which reads all required properties from file.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctionsClass
    use :: IO_HDF5                          , only : hdf5Object
    use :: HDF5_Access                      , only : hdf5Access
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast          , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisConcentrationDistribution)                                :: self
    character       (len=*                                  ), intent(in   )                 :: fileName
    type            (varying_string                         ), intent(in   )                 :: label                           , comment
    double precision                                         , intent(in   )                 :: rootVarianceFractionalMinimum
    class           (outputTimesClass                       ), intent(inout)                 :: outputTimes_
    class           (cosmologyParametersClass               ), intent(in   )                 :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), intent(inout)                 :: cosmologyFunctions_
    class           (nbodyHaloMassErrorClass                ), intent(in   )                 :: nbodyHaloMassError_
    class           (darkMatterProfileDMOClass              ), intent(in   )                 :: darkMatterProfileDMO_
    class           (virialDensityContrastClass             ), intent(in   )                 :: virialDensityContrastDefinition_, virialDensityContrast_
    double precision                                         , allocatable  , dimension(:  ) :: functionValueTarget             , concentration
    integer         (c_size_t                               ), allocatable  , dimension(:  ) :: functionCountTarget
    double precision                                         , allocatable  , dimension(:,:) :: functionCovarianceTarget
    type            (varying_string                         )                                :: targetLabel
    type            (hdf5Object                             )                                :: dataFile                        , simulationGroup      , &
         &                                                                                      attributesGroup
    integer                                                                                  :: i
    double precision                                                                         :: time                            , redshift             , &
         &                                                                                      timeRecent                      , massMinimum          , &
         &                                                                                      massMaximum                     , massParticle

    !$ call hdf5Access%set  ()
    call dataFile%openFile(fileName,readOnly=.true.)
    simulationGroup=dataFile       %openGroup('simulation0001')
    attributesGroup=simulationGroup%openGroup('simulation'    )
    call simulationGroup   %readDataset  ('concentration'                              ,concentration           )
    call simulationGroup   %readDataset  ('concentrationDistributionFunction'          ,functionValueTarget     )
    if (simulationGroup%hasDataset('concentrationDistributionFunctionCovariance')) then
       call simulationGroup%readDataset  ('concentrationDistributionFunctionCovariance',functionCovarianceTarget)
    else
       call simulationGroup%readDataset  ('count'                                      ,functionCountTarget     )
    end if
    call attributesGroup   %readAttribute('labelTarget'                                ,targetLabel             )
    call attributesGroup   %readAttribute('massMinimum'                                ,massMinimum             )
    call attributesGroup   %readAttribute('massMaximum'                                ,massMaximum             )
    call attributesGroup   %readAttribute('redshift'                                   ,redshift                )
    call attributesGroup   %readAttribute('massParticle'                               ,massParticle            )
    call attributesGroup   %readAttribute('timeRecent'                                 ,timeRecent              )
    call attributesGroup   %close        (                                                                      )
    call simulationGroup   %close        (                                                                      )
    call dataFile          %close        (                                                                      )
    !$ call hdf5Access%unset()
    ! Compute a (diagonal) covariance matrix from the counts if necessary.
    if (.not.allocated(functionCovarianceTarget)) then
       allocate(functionCovarianceTarget(size(functionValueTarget),size(functionValueTarget)))
       functionCovarianceTarget=0.0d0
       do i=1,size(functionValueTarget)
          if (functionCountTarget(i) > 0_c_size_t)                                                    &
               & functionCovarianceTarget(i,i)=functionValueTarget(i)**2/dble(functionCountTarget(i))
       end do
    end if
    ! Convert redshift to time.
    time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    ! Build the object.
    self=outputAnalysisConcentrationDistribution(label,comment,time,massMinimum,massMaximum,concentration(1),concentration(size(concentration)),size(concentration,kind=c_size_t),timeRecent,massParticle,rootVarianceFractionalMinimum,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,virialDensityContrast_,virialDensityContrastDefinition_,outputTimes_,targetLabel,functionValueTarget,functionCovarianceTarget)
    !![
    <constructorAssign variables="fileName"/>
    !!]
    return
  end function concentrationDistributionConstructorFile

  function concentrationDistributionConstructorInternal(label,comment,time,massMinimum,massMaximum,concentrationMinimum,concentrationMaximum,countConcentrations,timeRecent,massParticle,rootVarianceFractionalMinimum,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,virialDensityContrast_,virialDensityContrastDefinition_,outputTimes_,targetLabel,functionValueTarget,functionCovarianceTarget) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisConcentrationDistribution} output analysis class.
    !!}
    use :: Cosmology_Functions                     , only : cosmologyFunctionsClass
    use :: Galactic_Filters                        , only : filterList                                      , galacticFilterAll                           , galacticFilterFormationTime                   , galacticFilterHaloIsolated
    use :: ISO_Varying_String                      , only : var_str
    use :: Node_Property_Extractors                , only : nodePropertyExtractorConcentration              , nodePropertyExtractorMassHalo
    use :: Numerical_Ranges                        , only : Make_Range                                      , rangeTypeLogarithmic
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelPoisson
    use :: Output_Analysis_Distribution_Normalizers, only : normalizerList                                  , outputAnalysisDistributionNormalizerBinWidth, outputAnalysisDistributionNormalizerLog10ToLog, outputAnalysisDistributionNormalizerSequence, &
          &                                                 outputAnalysisDistributionNormalizerUnitarity
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorRndmErrNbdyCnc
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10         , outputAnalysisPropertyOperatorIdentity      , outputAnalysisPropertyOperatorLog10
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorNbodyMass
    use :: Output_Times                            , only : outputTimesClass
    use :: Statistics_NBody_Halo_Mass_Errors       , only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast                 , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisConcentrationDistribution         )                                             :: self
    type            (varying_string                                  )                             , intent(in   ) :: label                                               , comment
    class           (cosmologyParametersClass                        ), target                     , intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass                         ), target                     , intent(in   ) :: cosmologyFunctions_
    class           (outputTimesClass                                ), target                     , intent(inout) :: outputTimes_
    class           (nbodyHaloMassErrorClass                         ), target                     , intent(in   ) :: nbodyHaloMassError_
    class           (virialDensityContrastClass                      ), target                     , intent(in   ) :: virialDensityContrastDefinition_                    , virialDensityContrast_
    class           (darkMatterProfileDMOClass                       ), target                     , intent(in   ) :: darkMatterProfileDMO_
    double precision                                                                               , intent(in   ) :: massMinimum                                         , massMaximum                             , &
         &                                                                                                            concentrationMinimum                                , concentrationMaximum                    , &
         &                                                                                                            timeRecent                                          , time                                    , &
         &                                                                                                            massParticle                                        , rootVarianceFractionalMinimum
    integer         (c_size_t                                        )                             , intent(in   ) :: countConcentrations
    type            (varying_string                                  ), optional                   , intent(in   ) :: targetLabel
    double precision                                                  , optional   , dimension(:  ), intent(in   ) :: functionValueTarget
    double precision                                                  , optional   , dimension(:,:), intent(in   ) :: functionCovarianceTarget
    type            (nodePropertyExtractorConcentration              ), pointer                                    :: nodePropertyExtractor_
    type            (outputAnalysisPropertyOperatorLog10             ), pointer                                    :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10         ), pointer                                    :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorNbodyMass           ), pointer                                    :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorIdentity          ), pointer                                    :: outputAnalysisWeightPropertyOperator_
    type            (nodePropertyExtractorMassHalo                   ), pointer                                    :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisDistributionNormalizerSequence    ), pointer                                    :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionOperatorRndmErrNbdyCnc), pointer                                    :: outputAnalysisDistributionOperator_
    type            (outputAnalysisDistributionNormalizerUnitarity   ), pointer                                    :: outputAnalysisDistributionNormalizerUnitarity_
    type            (outputAnalysisDistributionNormalizerBinWidth    ), pointer                                    :: outputAnalysisDistributionNormalizerBinWidth_
    type            (outputAnalysisDistributionNormalizerLog10ToLog  ), pointer                                    :: outputAnalysisDistributionNormalizerLog10ToLog_
    type            (normalizerList                                  ), pointer                                    :: normalizer_
    type            (galacticFilterHaloIsolated                      ), pointer                                    :: galacticFilterHaloIsolated_
    type            (galacticFilterFormationTime                     ), pointer                                    :: galacticFilterFormationTime_
    type            (galacticFilterAll                               ), pointer                                    :: galacticFilter_
    type            (filterList                                      ), pointer                                    :: filters_
    double precision                                                  , allocatable, dimension(:  )                :: concentrations
    double precision                                                  , allocatable, dimension(:,:)                :: outputWeight
    double precision                                                  , parameter  , dimension(3  )                :: concentrationFitA                       =[                                                      &
         &                                                                                                                                                      -0.270d+00,                                           &
         &                                                                                                                                                      +1.780d+00,                                           &
         &                                                                                                                                                      -0.490d+00                                            &
         &                                                                                                                                                     ]
    double precision                                                  , parameter                                  :: concentrationFitB                       = -0.550d+00
    integer                                                           , parameter                                  :: covarianceBinomialBinsPerDecade         =  2
    double precision                                                  , parameter                                  :: covarianceBinomialMassHaloMinimum       = +3.000d+11, covarianceBinomialMassHaloMaximum=1.0d15
    integer         (c_size_t                                        )                                             :: iOutput                                             , bufferCount
    !![
    <constructorAssign variables="rootVarianceFractionalMinimum, massMinimum, massMaximum, concentrationMinimum, concentrationMaximum, timeRecent, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *nbodyHaloMassError_, *virialDensityContrastDefinition_, *virialDensityContrast_"/>
    !!]

    ! Set parameters needed for descriptor.
    self%redshift                    =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    self%countConcentrationsPerDecade=dble(countConcentrations-1_c_size_t)/log10(concentrationMaximum/concentrationMinimum)
    ! Build grid of concentrations.
    allocate(concentrations(countConcentrations))
    concentrations=Make_Range(concentrationMinimum,concentrationMaximum,int(countConcentrations),rangeType=rangeTypeLogarithmic)
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(countConcentrations,outputTimes_%count()))
    outputWeight=0.0d0
    do iOutput=1,outputTimes_%count()
       if (Values_Agree(outputTimes_%time(iOutput),time,absTol=1.0d-10)) outputWeight(:,iOutput)=1.0d0
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
    <referenceConstruct object="galacticFilterHaloIsolated_"  constructor="galacticFilterHaloIsolated  (          )"/>
    <referenceConstruct object="galacticFilterFormationTime_" constructor="galacticFilterFormationTime (timeRecent)"/>
    <referenceConstruct object="galacticFilter_"              constructor="galacticFilterAll           (filters_  )"/>
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
    ! Create a concentration property extractor.
    allocate(nodePropertyExtractor_      )
    !![
    <referenceConstruct object="nodePropertyExtractor_">
     <constructor>
      nodePropertyExtractorConcentration              (                                        &amp;
      &amp;                                            .false.                               , &amp;
      &amp;                                            cosmologyParameters_                  , &amp;
      &amp;                                            cosmologyFunctions_                   , &amp;
      &amp;                                            darkMatterProfileDMO_                 , &amp;
      &amp;                                            virialDensityContrast_                , &amp;
      &amp;                                            virialDensityContrastDefinition_        &amp;
      &amp;                                           )
     </constructor>
    </referenceConstruct>
    !!]
    ! Create a halo mass property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_">
     <constructor>
      nodePropertyExtractorMassHalo                   (                                        &amp;
      &amp;                                            .false.                               , &amp;
      &amp;                                            cosmologyFunctions_                   , &amp;
      &amp;                                            cosmologyParameters_                  , &amp;
      &amp;                                            darkMatterProfileDMO_                 , &amp;
      &amp;                                            virialDensityContrast_                , &amp;
      &amp;                                            virialDensityContrastDefinition_        &amp;
      &amp;                                           )
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
    self%outputAnalysisVolumeFunction1D=                                                                &
         & outputAnalysisVolumeFunction1D(                                                              &
         &                                var_str('concentrationDistribution')//label                 , &
         &                                comment                                                     , &
         &                                var_str('concentration'                                    ), &
         &                                var_str('Concentration at the bin center'                  ), &
         &                                var_str('dimensionless'                                    ), &
         &                                0.0d0                                                       , &
         &                                var_str('concentrationFunction'                            ), &
         &                                var_str('Concentration distribution averaged over each bin'), &
         &                                var_str('dimensionless'                                    ), &
         &                                0.0d0                                                       , &
         &                                log10(concentrations)                                       , &
         &                                bufferCount                                                 , &
         &                                outputWeight                                                , &
         &                                nodePropertyExtractor_                                      , &
         &                                outputAnalysisPropertyOperator_                             , &
         &                                outputAnalysisPropertyUnoperator_                           , &
         &                                outputAnalysisWeightOperator_                               , &
         &                                outputAnalysisDistributionOperator_                         , &
         &                                outputAnalysisDistributionNormalizer_                       , &
         &                                galacticFilter_                                             , &
         &                                outputTimes_                                                , &
         &                                outputAnalysisCovarianceModelPoisson                        , &
         &                                covarianceBinomialBinsPerDecade                             , &
         &                                covarianceBinomialMassHaloMinimum                           , &
         &                                covarianceBinomialMassHaloMaximum                           , &
         &                                .false.                                                     , &
         &                                var_str('$c$'                                              ), &
         &                                var_str('$\mathrm{d}p/\mathrm{d}\log c$'                   ), &
         &                                .true.                                                      , &
         &                                .false.                                                     , &
         &                                targetLabel                                                 , &
         &                                functionValueTarget                                         , &
         &                                functionCovarianceTarget                                      &
         &                               )
    !![
    <objectDestructor name="galacticFilterHaloIsolated_"                    />
    <objectDestructor name="galacticFilterFormationTime_"                   />
    <objectDestructor name="galacticFilter_"                                />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"          />
    <objectDestructor name="outputAnalysisPropertyOperator_"                />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"              />
    <objectDestructor name="nodePropertyExtractor_"                         />
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
  end function concentrationDistributionConstructorInternal

  double precision function concentrationDistributionLogLikelihood(self)
    !!{
    Return the log-likelihood of the concentration distribution function.
    !!}
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use            :: Error                       , only : Error_Report
    use            :: Interface_GSL               , only : GSL_Success
    use            :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisConcentrationDistribution), intent(inout)                 :: self
    double precision                                         , allocatable  , dimension(:,:) :: functionCovarianceCombined
    double precision                                         , allocatable  , dimension(:  ) :: functionValueDifference
    logical                                                  , allocatable  , dimension(:  ) :: mask
    type            (vector                                 )                                :: residual
    type            (matrix                                 )                                :: covariance
    integer         (c_size_t                               )                                :: i                         , j , &
         &                                                                                      ii                        , jj
    integer                                                                                  :: status
    
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
                      functionCovarianceCombined       (ii,jj)=    +self%functionCovarianceTarget     (i,j)       &
                           &                                       +self%functionCovariance           (i,j)
                      if (ii == jj) &
                           & functionCovarianceCombined(ii,jj)=max(                                               &
                           &                                       +functionCovarianceCombined        (ii,jj)   , &
                           &                                       +self%functionValueTarget          ( i   )     &
                           &                                       *self%functionValueTarget          (    j)     &
                           &                                       *self%rootVarianceFractionalMinimum       **2  &
                           &                                      )
                    end if
                end do
             end if
          end do
          residual  =vector(functionValueDifference   )
          covariance=matrix(functionCovarianceCombined)
          ! Compute the log-likelihood.
          concentrationDistributionLogLikelihood=-0.5d0*covariance%covarianceProduct(residual,status)
          if (status /= GSL_Success) concentrationDistributionLogLikelihood=logImprobable
       else
          concentrationDistributionLogLikelihood=+0.0d0
       end if
    else
       concentrationDistributionLogLikelihood   =+0.0d0
       call Error_Report('no target distribution was provided for likelihood calculation'//{introspection:location})
    end if
    return
  end function concentrationDistributionLogLikelihood

  subroutine concentrationDistributionDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisConcentrationDistribution} output analysis class.
    !!}
    type(outputAnalysisConcentrationDistribution), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%nbodyHaloMassError_"             />
    <objectDestructor name="self%darkMatterProfileDMO_"           />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine concentrationDistributionDestructor
