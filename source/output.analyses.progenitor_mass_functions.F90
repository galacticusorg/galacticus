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
  Implements a concentration distribution output analysis class for dark matter halo progenitor mass functions.
  !!}
  
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOClass
  use :: Galactic_Filters                , only : galacticFilterAll
  use :: Node_Property_Extractors        , only : nodePropertyExtractorMassHalo
  use :: Output_Analysis_Weight_Operators, only : outputAnalysisWeightOperatorNbodyMass

  !![
  <outputAnalysis name="outputAnalysisProgenitorMassFunction">
   <description>A dark matter halo progenitor mass function output analysis class.</description>
   <deepCopy>
    <functionClass variables="galacticFilterParentMass_, outputAnalysisWeightOperatorNbodyMass_, nodePropertyExtractorMassParent_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="galacticFilterParentMass_, outputAnalysisWeightOperatorNbodyMass_, nodePropertyExtractorMassParent_"/>
   </stateStorable>
   <runTimeFileDependencies paths="fileName"/>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisProgenitorMassFunction
     !!{
     A dark matter halo progenitor mass function output analysis class.
     !!}
     private
     class           (cosmologyFunctionsClass              ), pointer                   :: cosmologyFunctions_                    => null()
     class           (cosmologyParametersClass             ), pointer                   :: cosmologyParameters_                   => null()
     class           (nbodyHaloMassErrorClass              ), pointer                   :: nbodyHaloMassError_                    => null()
     class           (virialDensityContrastClass           ), pointer                   :: virialDensityContrastDefinition_       => null(), virialDensityContrast_     => null()
     class           (darkMatterProfileDMOClass            ), pointer                   :: darkMatterProfileDMO_                  => null()
     type            (galacticFilterAll                    ), pointer                   :: galacticFilterParentMass_              => null()
     type            (outputAnalysisWeightOperatorNbodyMass), pointer                   :: outputAnalysisWeightOperatorNbodyMass_ => null()
     type            (nodePropertyExtractorMassHalo        ), pointer                   :: nodePropertyExtractorMassParent_       => null()
     double precision                                       , allocatable, dimension(:) :: rootVarianceTargetFractional
     double precision                                                                   :: massRatioMinimum                                , massRatioMaximum                    , &
          &                                                                                massParentMinimum                               , massParentMaximum                   , &
          &                                                                                timeProgenitor                                  , timeParent                          , &
          &                                                                                redshiftProgenitor                              , redshiftParent                      , &
          &                                                                                weightParents                                   , massRatioLikelihoodMinimum          , &
          &                                                                                massRatioLikelihoodMaximum
     integer                                                                            :: indexParent                                     , indexRedshift
     integer         (c_size_t                             )                            :: countMassRatio                                  , indexOutput
     logical                                                                            :: alwaysIsolatedOnly                              , covarianceDiagonalize               , &
          &                                                                                covarianceTargetOnly                            , likelihoodInLog                     , &
          &                                                                                weightsFinalized
     type            (varying_string                       )                            :: fileName
  contains
     final     ::                     progenitorMassFunctionDestructor
     procedure :: newTree          => progenitorMassFunctionNewTree
     procedure :: reduce           => progenitorMassFunctionReduce
     procedure :: finalizeAnalysis => progenitorMassFunctionFinalizeAnalysis
     procedure :: finalize         => progenitorMassFunctionFinalize
     procedure :: logLikelihood    => progenitorMassFunctionLogLikelihood
  end type outputAnalysisProgenitorMassFunction

  interface outputAnalysisProgenitorMassFunction
     !!{
     Constructors for the {\normalfont \ttfamily progenitorMassFunction} output analysis class.
     !!}
     module procedure progenitorMassFunctionConstructorParameters
     module procedure progenitorMassFunctionConstructorFile
     module procedure progenitorMassFunctionConstructorInternal
  end interface outputAnalysisProgenitorMassFunction

contains
  
  function progenitorMassFunctionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily progenitorMassFunction} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctionsClass
    use :: Input_Parameters                 , only : inputParameter            , inputParameters
    use :: ISO_Varying_String               , only : char
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast          , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisProgenitorMassFunction)                              :: self
    type            (inputParameters                     ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass             ), pointer                     :: cosmologyFunctions_
    class           (cosmologyParametersClass            ), pointer                     :: cosmologyParameters_
    class           (nbodyHaloMassErrorClass             ), pointer                     :: nbodyHaloMassError_
    class           (outputTimesClass                    ), pointer                     :: outputTimes_
    class           (virialDensityContrastClass          ), pointer                     :: virialDensityContrastDefinition_, virialDensityContrast_
    class           (darkMatterProfileDMOClass           ), pointer                     :: darkMatterProfileDMO_
    double precision                                      , dimension(:  ), allocatable :: functionValueTarget             , functionCovarianceTarget1D, &
         &                                                                                 rootVarianceTargetFractional
    double precision                                      , dimension(:,:), allocatable :: functionCovarianceTarget
    double precision                                                                    :: massRatioMinimum                , massRatioMaximum          , &
         &                                                                                 massParentMinimum               , massParentMaximum         , &
         &                                                                                 redshiftProgenitor              , redshiftParent            , &
         &                                                                                 massRatioLikelihoodMinimum      , massRatioLikelihoodMaximum
    integer         (c_size_t                            )                              :: countMassRatio
    integer                                                                             :: indexParent                     , indexRedshift
    type            (varying_string                      )                              :: label                           , comment                   , &
         &                                                                                 targetLabel                     , fileName
    logical                                                                             :: alwaysIsolatedOnly              , covarianceDiagonalize     , &
          &                                                                                covarianceTargetOnly            , likelihoodInLog
    
    allocate(rootVarianceTargetFractional(max(1,parameters%count('rootVarianceTargetFractional',zeroIfNotPresent=.true.))))
    !![
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"             source="parameters"                                                />
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"              source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_"           source="parameters"                                                />
    <objectBuilder class="nbodyHaloMassError"    name="nbodyHaloMassError_"              source="parameters"                                                /> 
    <objectBuilder class="outputTimes"           name="outputTimes_"                     source="parameters"                                                />
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"            source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrastDefinition_" source="parameters" parameterName="virialDensityContrastDefinition"/>
    <inputParameter>
      <name>covarianceDiagonalize</name>
      <source>parameters</source>
      <description>If true, all off-diagonal elements of the covariance matrix are set to zero.</description>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>covarianceTargetOnly</name>
      <source>parameters</source>
      <description>If true, only the covariance of the target dataset is accounted for (otherwise the model covariance is added).</description>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>rootVarianceTargetFractional</name>
      <source>parameters</source>
      <description>The diagonal of the covariance matrix is forced to be at least equal to this fraction multiplied by the target dataset squared. This may be a list of values corresponding to each element of the target dataset. If the list is shorter than the target dataset the final value in the list is applied to all remaining elements in the target dataset.</description>
      <defaultValue>[0.0d0]</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>likelihoodInLog</name>
      <source>parameters</source>
      <description>If true, the likelihood is computed in $\log\phi$ instead of in $\phi$.</description>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>massRatioLikelihoodMinimum</name>
      <source>parameters</source>
      <description>The minimum mass ratio to include in likelihood calculations.</description>
      <defaultValue>0.0d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>massRatioLikelihoodMaximum</name>
      <source>parameters</source>
      <description>The maximum mass ratio to include in likelihood calculations.</description>
      <defaultValue>huge(0.0d0)</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>alwaysIsolatedOnly</name>
      <source>parameters</source>
      <description>If true, include only progenitors which have been always isolated halos.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftParent</name>
      <source>parameters</source>
      <description>Redshift of the parent halos.</description>
    </inputParameter>
    !!]
    if (parameters%isPresent('fileName')) then
       !![
       <inputParameter>
         <name>fileName</name>
         <source>parameters</source>
         <description>The name of the file from which to read progenitor mass function parameters.</description>
       </inputParameter>
       <inputParameter>
         <name>indexParent</name>
         <source>parameters</source>
         <description>The parent mass index to use from the file.</description>
       </inputParameter>
       <inputParameter>
         <name>indexRedshift</name>
         <source>parameters</source>
         <description>The progenitor redshift index to use from the file.</description>
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
       <inputParameter>
         <name>targetLabel</name>
         <source>parameters</source>
         <description>Label for the target dataset.</description>
       </inputParameter>
       !!]
       self=outputAnalysisProgenitorMassFunction(char(fileName),label,comment,targetLabel,indexParent,indexRedshift,redshiftParent,massRatioLikelihoodMinimum,massRatioLikelihoodMaximum,covarianceDiagonalize,covarianceTargetOnly,rootVarianceTargetFractional,likelihoodInLog,alwaysIsolatedOnly,darkMatterProfileDMO_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_,nbodyHaloMassError_,outputTimes_)
    else
       !![
       <inputParameter>
         <name>label</name>
         <source>parameters</source>
         <variable>label</variable>
         <description>A label for the progenitor mass function.</description>
       </inputParameter>
       <inputParameter>
         <name>comment</name>
         <source>parameters</source>
         <variable>comment</variable>
         <description>A descriptive comment for the progenitor mass function.</description>
       </inputParameter>
       <inputParameter>
         <name>massRatioMinimum</name>
         <source>parameters</source>
         <description>Minimum mass ratio for the progenitor mass function.</description>
       </inputParameter>
       <inputParameter>
         <name>massRatioMaximum</name>
         <source>parameters</source>
         <description>Maximum mass ratio for the progenitor mass function.</description>
       </inputParameter>
       <inputParameter>
         <name>countMassRatio</name>
         <source>parameters</source>
         <description>Number of mass ratios at which to compute the progenitor mass function.</description>
       </inputParameter>
       <inputParameter>
         <name>massParentMinimum</name>
         <source>parameters</source>
         <description>Minimum mass of the parent halo for the progenitor mass function.</description>
       </inputParameter>
       <inputParameter>
         <name>massParentMaximum</name>
         <source>parameters</source>
         <description>Maximum mass of the parent halo for the progenitor mass function.</description>
       </inputParameter>
       <inputParameter>
         <name>redshiftProgenitor</name>
         <source>parameters</source>
         <description>Redshift of the progenitor halos.</description>
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
       <conditionalCall>
        <call>
         self=outputAnalysisProgenitorMassFunction(                                                                                                     &amp;
          &amp;                                    label                                                                                              , &amp;
          &amp;                                    comment                                                                                            , &amp;
          &amp;                                    massRatioMinimum                                                                                   , &amp;
          &amp;                                    massRatioMaximum                                                                                   , &amp;
          &amp;                                    countMassRatio                                                                                     , &amp;
          &amp;                                    massParentMinimum                                                                                  , &amp;
          &amp;                                    massParentMaximum                                                                                  , &amp;
          &amp;                                    cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftProgenitor)), &amp;
          &amp;                                    cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftParent    )), &amp;
          &amp;                                    alwaysIsolatedOnly                                                                                 , &amp;
          &amp;                                    massRatioLikelihoodMinimum                                                                         , &amp;
          &amp;                                    massRatioLikelihoodMaximum                                                                         , &amp;
          &amp;                                    covarianceDiagonalize                                                                              , &amp;
          &amp;                                    covarianceTargetOnly                                                                               , &amp;
          &amp;                                    rootVarianceTargetFractional                                                                       , &amp;
          &amp;                                    likelihoodInLog                                                                                    , &amp;
	  &amp;                                    darkMatterProfileDMO_                                                                              , &amp;
          &amp;                                    cosmologyParameters_                                                                               , &amp;
          &amp;                                    cosmologyFunctions_                                                                                , &amp;
          &amp;                                    virialDensityContrast_                                                                             , &amp;
          &amp;                                    virialDensityContrastDefinition_                                                                   , &amp;
          &amp;                                    nbodyHaloMassError_                                                                                , &amp;
          &amp;                                    outputTimes_                                                                                         &amp;
          &amp;                                    {conditions}                                                                                         &amp;
          &amp;                                   )
        </call>
        <argument name="targetLabel"              value="targetLabel"              parameterPresent="parameters"/>
        <argument name="functionValueTarget"      value="functionValueTarget"      parameterPresent="parameters"/>
        <argument name="functionCovarianceTarget" value="functionCovarianceTarget" parameterPresent="parameters"/>
       </conditionalCall>
       <inputParametersValidate source="parameters"/>
       !!]
    end if
    !![
    <objectDestructor name="cosmologyFunctions_"             />
    <objectDestructor name="cosmologyParameters_"            />
    <objectDestructor name="virialDensityContrast_"          />
    <objectDestructor name="outputTimes_"                    />
    <objectDestructor name="nbodyHaloMassError_"             />
    <objectDestructor name="darkMatterProfileDMO_"           />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function progenitorMassFunctionConstructorParameters
  
  function progenitorMassFunctionConstructorFile(fileName,label,comment,targetLabel,indexParent,indexRedshift,redshiftParent,massRatioLikelihoodMinimum,massRatioLikelihoodMaximum,covarianceDiagonalize,covarianceTargetOnly,rootVarianceTargetFractional,likelihoodInLog,alwaysIsolatedOnly,darkMatterProfileDMO_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_,nbodyHaloMassError_,outputTimes_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily progenitorMassFunction} output analysis class which reads all required properties from file.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctionsClass
    use :: HDF5_Access                      , only : hdf5Access
    use :: IO_HDF5                          , only : hdf5Object
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass
    use :: File_Utilities                   , only : File_Name_Expand
    use :: Virial_Density_Contrast          , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisProgenitorMassFunction)                                               :: self
    character       (len=*                               ), intent(in   )                                :: fileName
    type            (varying_string                      ), intent(in   )                                :: label                          , targetLabel               , &
         &                                                                                                  comment
    integer                                               , intent(in   )                                :: indexParent                    , indexRedshift
    double precision                                      , intent(in   )                                :: massRatioLikelihoodMinimum     , massRatioLikelihoodMaximum, &
         &                                                                                                  redshiftParent
    double precision                                      , intent(in   ), allocatable, dimension(:    ) :: rootVarianceTargetFractional
    logical                                               , intent(in   )                                :: covarianceDiagonalize          , covarianceTargetOnly      , &
         &                                                                                                  likelihoodInLog                , alwaysIsolatedOnly
    class           (cosmologyFunctionsClass             ), intent(inout), target                        :: cosmologyFunctions_
    class           (cosmologyParametersClass            ), intent(inout), target                        :: cosmologyParameters_
    class           (outputTimesClass                    ), intent(inout), target                        :: outputTimes_
    class           (virialDensityContrastClass          ), intent(in   ), target                        :: virialDensityContrastDefinition_, virialDensityContrast_
    class           (nbodyHaloMassErrorClass             ), intent(in   ), target                        :: nbodyHaloMassError_
    class           (darkMatterProfileDMOClass           ), intent(in   ), target                        :: darkMatterProfileDMO_
    double precision                                                                                     :: massParentMinimum               , massParentMaximum        , &
         &                                                                                                  timeProgenitor                  , timeParent               , &
         &                                                                                                  redshiftProgenitor
    double precision                                                     , allocatable, dimension(:    ) :: functionValueTarget             , massRatio                , &
         &                                                                                                  redshiftProgenitors             , massParents              , &
         &                                                                                                  massParentsMinimum              , massParentsMaximum
    double precision                                                     , allocatable, dimension(:,:  ) :: functionCovarianceTarget
    double precision                                                     , allocatable, dimension(:,:,:) :: functionValuesTarget
    integer         (c_size_t                            )               , allocatable, dimension(:,:,:) :: functionCountsTarget
    type            (hdf5Object                          )                                               :: dataFile                        , simulationGroup
    integer                                                                                              :: i
    logical                                                                                              :: haveBoundaries
    
    !$ call hdf5Access%set  ()
    call dataFile%openFile(char(File_Name_Expand(fileName)),readOnly=.true.)
    simulationGroup=dataFile       %openGroup ('simulation0001'   )
    haveBoundaries =simulationGroup%hasDataset('massParentMinimum')
    call    simulationGroup%readDataset('massRatioProgenitor'   ,massRatio           )
    call    simulationGroup%readDataset('progenitorMassFunction',functionValuesTarget)
    call    simulationGroup%readDataset('count'                 ,functionCountsTarget)
    call    simulationGroup%readDataset('massParent'            ,massParents         )
    call    simulationGroup%readDataset('redshiftProgenitor'    ,redshiftProgenitors )
    if (haveBoundaries) then
       call simulationGroup%readDataset('massParentMinimum'     ,massParentsMinimum  )
       call simulationGroup%readDataset('massParentMaximum'     ,massParentsMaximum  )
    end if
    call    simulationGroup%close      (                                             )
    call    dataFile       %close      (                                             )
    !$ call hdf5Access%unset()
    ! Extract parent mass range and progenitor redshift.
    if (haveBoundaries) then
       massParentMinimum=massParentsMinimum (indexParent  +1)
       massParentMaximum=massParentsMaximum (indexParent  +1)
    else
       massParentMinimum=massParents        (indexParent  +1)/sqrt(massParents(2)/massParents(1))
       massParentMaximum=massParents        (indexParent  +1)*sqrt(massParents(2)/massParents(1))
    end if
    redshiftProgenitor  =redshiftProgenitors(indexRedshift+1)
    ! Extract the target function values.
    allocate(functionValueTarget(size(massRatio)))
    functionValueTarget=functionValuesTarget(indexParent+1,:,indexRedshift+1)
    ! Compute a (diagonal) covariance matrix from the counts.
    allocate(functionCovarianceTarget(size(functionValueTarget),size(functionValueTarget)))
    functionCovarianceTarget=0.0d0
    do i=1,size(functionValueTarget)
       if (functionCountsTarget(indexParent+1,i,indexRedshift+1) > 0_c_size_t) &
            & functionCovarianceTarget(i,i)=functionValueTarget(i)**2/dble(functionCountsTarget(indexParent+1,i,indexRedshift+1))
    end do
    ! Convert redshifts to times.
    timeProgenitor    =cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftProgenitor))
    timeParent        =cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftParent    ))
    ! Build the object.
    self              =outputAnalysisProgenitorMassFunction(label,comment,massRatio(1),massRatio(size(massRatio)),size(massRatio,kind=c_size_t),massParentMinimum,massParentMaximum,timeProgenitor,timeParent,alwaysIsolatedOnly,massRatioLikelihoodMinimum,massRatioLikelihoodMaximum,covarianceDiagonalize,covarianceTargetOnly,rootVarianceTargetFractional,likelihoodInLog,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,virialDensityContrastDefinition_,nbodyHaloMassError_,outputTimes_,targetLabel,functionValueTarget,functionCovarianceTarget)
    !![
    <constructorAssign variables="fileName, indexParent, indexRedshift"/>
    !!]
    return
  end function progenitorMassFunctionConstructorFile

  function progenitorMassFunctionConstructorInternal(label,comment,massRatioMinimum,massRatioMaximum,countMassRatio,massParentMinimum,massParentMaximum,timeProgenitor,timeParent,alwaysIsolatedOnly,massRatioLikelihoodMinimum,massRatioLikelihoodMaximum,covarianceDiagonalize,covarianceTargetOnly,rootVarianceTargetFractional,likelihoodInLog,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,virialDensityContrastDefinition_,nbodyHaloMassError_,outputTimes_,targetLabel,functionValueTarget,functionCovarianceTarget) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily progenitorMassFunction} output analysis class.
    !!}
    use :: Galactic_Filters                        , only : filterList                                      , galacticFilterDescendantNode                , galacticFilterHaloAlwaysIsolated              , galacticFilterHaloIsolated                  , &
          &                                                 galacticFilterHaloMass                          , galacticFilterNot
    use :: Node_Property_Extractors                , only : nodePropertyExtractorDescendantNode             , nodePropertyExtractorRatio
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Numerical_Ranges                        , only : Make_Range                                      , rangeTypeLogarithmic
    use :: Cosmology_Functions                     , only : cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                    , only : cosmologyParametersSimple
    use :: Output_Analysis_Distribution_Normalizers, only : normalizerList                                  , outputAnalysisDistributionNormalizerBinWidth, outputAnalysisDistributionNormalizerLog10ToLog, outputAnalysisDistributionNormalizerSequence
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorProperty            , outputAnalysisWeightOperatorSequence        , outputAnalysisWeightOperatorSubsampling       , weightOperatorList
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorMassRatioNBody
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorAntiLog10         , outputAnalysisPropertyOperatorIdentity      , outputAnalysisPropertyOperatorLog10
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelPoisson
    use :: Statistics_NBody_Halo_Mass_Errors       , only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast                 , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisProgenitorMassFunction            )                                          :: self
    type            (varying_string                                  ), intent(in   )                           :: label                                                  , comment
    double precision                                                  , intent(in   )                           :: massRatioMinimum                                       , massRatioMaximum                        , &
         &                                                                                                         massParentMinimum                                      , massParentMaximum                       , &
         &                                                                                                         timeProgenitor                                         , timeParent
    double precision                                                  , intent(in   )          , dimension(:  ) :: rootVarianceTargetFractional
    integer         (c_size_t                                        ), intent(in   )                           :: countMassRatio
    logical                                                           , intent(in   )                           :: alwaysIsolatedOnly                                     , covarianceDiagonalize                   , &
         &                                                                                                         covarianceTargetOnly                                   , likelihoodInLog
    double precision                                                  , intent(in   )                           :: massRatioLikelihoodMinimum                             , massRatioLikelihoodMaximum
    class           (cosmologyParametersClass                        ), intent(inout), target                   :: cosmologyParameters_
    class           (cosmologyFunctionsClass                         ), intent(inout), target                   :: cosmologyFunctions_
    class           (outputTimesClass                                ), intent(inout), target                   :: outputTimes_
    class           (nbodyHaloMassErrorClass                         ), intent(in   ), target                   :: nbodyHaloMassError_
    class           (virialDensityContrastClass                      ), intent(in   ), target                   :: virialDensityContrastDefinition_                       , virialDensityContrast_
    class           (darkMatterProfileDMOClass                       ), intent(in   ), target                   :: darkMatterProfileDMO_
    type            (varying_string                                  ), intent(in   ), optional                 :: targetLabel
    double precision                                                  , intent(in   ), optional, dimension(:  ) :: functionValueTarget
    double precision                                                  , intent(in   ), optional, dimension(:,:) :: functionCovarianceTarget
    double precision                                                  , parameter                               :: timeTolerance                                  =1.0d-04
    double precision                                                  , parameter                               :: massRatioBuffer                                =1.0d-01
    integer                                                           , parameter                               :: covarianceBinomialBinsPerDecade                =2
    double precision                                                  , parameter                               :: covarianceBinomialMassHaloMinimum              =3.0d+11, covarianceBinomialMassHaloMaximum=1.0d15
    logical                                                           , parameter                               :: allowSelf                                      =.false.
    double precision                                                  , allocatable            , dimension(:  ) :: massRatios
    double precision                                                  , allocatable            , dimension(:,:) :: outputWeight
    type            (galacticFilterAll                               ), pointer                                 :: galacticFilter_
    type            (galacticFilterHaloIsolated                      ), pointer                                 :: galacticFilterHaloIsolated_
    type            (galacticFilterDescendantNode                    ), pointer                                 :: galacticFilterParentNode_
    type            (galacticFilterNot                               ), pointer                                 :: galacticFilterNot_
    type            (galacticFilterHaloMass                          ), pointer                                 :: galacticFilterProgenitorMass_                          , galacticFilterParentMassMinimum_        , &
         &                                                                                                         galacticFilterParentMassMaximum_
    type            (galacticFilterHaloAlwaysIsolated                ), pointer                                 :: galacticFilterHaloAlwaysIsolated_
    type            (filterList                                      ), pointer                                 :: filters_                                               , filtersParent_
    type            (nodePropertyExtractorMassHalo                   ), pointer                                 :: nodePropertyExtractorMassProgenitor_
    type            (nodePropertyExtractorRatio                      ), pointer                                 :: nodePropertyExtractorMassRatio_
    type            (nodePropertyExtractorDescendantNode             ), pointer                                 :: nodePropertyExtractorParentNode_
    type            (outputAnalysisDistributionNormalizerSequence    ), pointer                                 :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionNormalizerBinWidth    ), pointer                                 :: outputAnalysisDistributionNormalizerBinWidth_
    type            (outputAnalysisDistributionNormalizerLog10ToLog  ), pointer                                 :: outputAnalysisDistributionNormalizerLog10ToLog_
    type            (normalizerList                                  ), pointer                                 :: normalizer_
    type            (outputAnalysisWeightOperatorSequence            ), pointer                                 :: outputAnalysisWeightOperator_
    type            (weightOperatorList                              ), pointer                                 :: weightOperator_
    type            (outputAnalysisWeightOperatorProperty            ), pointer                                 :: outputAnalysisWeightOperatorMassRatio_
    type            (outputAnalysisWeightOperatorSubsampling         ), pointer                                 :: outputAnalysisWeightOperatorSubsampling_
    type            (outputAnalysisDistributionOperatorMassRatioNBody), pointer                                 :: outputAnalysisDistributionOperator_
    type            (outputAnalysisPropertyOperatorLog10             ), pointer                                 :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10         ), pointer                                 :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorIdentity          ), pointer                                 :: outputAnalysisPropertyIdentity_
    integer         (c_size_t                                        )                                          :: iOutput                                                , bufferCount
    !![
    <constructorAssign variables="massRatioMinimum, massRatioMaximum, countMassRatio, massParentMinimum, massParentMaximum, timeProgenitor, timeParent, alwaysIsolatedOnly, massRatioLikelihoodMinimum, massRatioLikelihoodMaximum, covarianceDiagonalize, covarianceTargetOnly, rootVarianceTargetFractional, likelihoodInLog, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *virialDensityContrast_, *virialDensityContrastDefinition_, *nbodyHaloMassError_, *outputTimes_"/>
    !!]

    ! Initialize state.
    self%weightsFinalized=.false.
    ! Set redshifts.
    self%redshiftProgenitor=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeProgenitor))
    self%redshiftParent    =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeParent    ))
    ! Build grid of mass ratios.
    allocate(massRatios(countMassRatio))
    massRatios=Make_Range(massRatioMinimum,massRatioMaximum,int(countMassRatio),rangeType=rangeTypeLogarithmic)
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(countMassRatio,outputTimes_%count()))
    outputWeight=0.0d0
    do iOutput=1,outputTimes_%count()
       if (Values_Agree(outputTimes_%time(iOutput),timeProgenitor,absTol=timeTolerance)) then
          outputWeight(:,iOutput)=1.0d0
          self%indexOutput=iOutput
       end if
    end do
    ! Initialize accumulated weight of parent nodes.
    self%weightParents=0.0d0
    ! Build a filter which selects isolated halos, above a suitable lower mass, and with parents in the correct mass range.
    allocate(     galacticFilter_                           )
    allocate(     galacticFilterHaloIsolated_               )
    allocate(     galacticFilterProgenitorMass_             )
    allocate(     galacticFilterParentMassMinimum_          )
    allocate(     galacticFilterParentMassMaximum_          )
    allocate(self%galacticFilterParentMass_                 )
    allocate(     galacticFilterNot_                        )
    allocate(     galacticFilterParentNode_                 )
    allocate(     filters_                                  )
    allocate(     filters_                        %next     )
    allocate(     filters_                        %next%next)
    allocate(     filtersParent_                            )
    allocate(     filtersParent_                  %next     )
    filters_                %filter_ => galacticFilterHaloIsolated_
    filters_      %next     %filter_ => galacticFilterProgenitorMass_
    filters_      %next%next%filter_ => galacticFilterParentNode_
    filtersParent_          %filter_ => galacticFilterParentMassMinimum_
    filtersParent_%next     %filter_ => galacticFilterNot_
    if (self%alwaysIsolatedOnly) then
       allocate(galacticFilterHaloAlwaysIsolated_               )
       allocate(filters_                         %next%next%next)
       allocate(filtersParent_                   %next%next     )
       filters_      %next%next%next%filter_ => galacticFilterHaloAlwaysIsolated_
       filtersParent_%next%next     %filter_ => galacticFilterHaloAlwaysIsolated_
       !![
       <referenceConstruct object="galacticFilterHaloAlwaysIsolated_" constructor="galacticFilterHaloAlwaysIsolated()"/>
       !!]
    else
       nullify(galacticFilterHaloAlwaysIsolated_)
    end if
    !![
    <referenceConstruct                             object="galacticFilterHaloIsolated_"      constructor="galacticFilterHaloIsolated  (                                                                                           )"/>
    <referenceConstruct                             object="galacticFilterProgenitorMass_"    constructor="galacticFilterHaloMass      (massParentMinimum*massRatioMinimum*massRatioBuffer,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_)"/>
    <referenceConstruct                             object="galacticFilterParentMassMinimum_" constructor="galacticFilterHaloMass      (massParentMinimum                                 ,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_)"/>
    <referenceConstruct                             object="galacticFilterParentMassMaximum_" constructor="galacticFilterHaloMass      (massParentMaximum                                 ,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_)"/>
    <referenceConstruct                             object="galacticFilterNot_"               constructor="galacticFilterNot           (galacticFilterParentMassMaximum_                                                                               )"/>
    <referenceConstruct isResult="yes" owner="self" object="galacticFilterParentMass_"        constructor="galacticFilterAll           (filtersParent_                                                                                                 )"/>
    <referenceConstruct                             object="galacticFilterParentNode_"        constructor="galacticFilterDescendantNode(timeParent                                        ,allowSelf,cosmologyFunctions_,self%galacticFilterParentMass_)"/>
    <referenceConstruct                             object="galacticFilter_"                  constructor="galacticFilterAll           (filters_                                                                                                       )"/>
    !!]
    ! Build a node property extractor which gives the ratio of the progenitor and parent halo masses.
    allocate(     nodePropertyExtractorMassProgenitor_)
    allocate(self%nodePropertyExtractorMassParent_    )
    allocate(     nodePropertyExtractorMassRatio_     )
    allocate(     nodePropertyExtractorParentNode_    )
    !![
    <referenceConstruct                             object="nodePropertyExtractorMassProgenitor_" constructor="nodePropertyExtractorMassHalo      (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_                           )"/>
    <referenceConstruct isResult="yes" owner="self" object="nodePropertyExtractorMassParent_"     constructor="nodePropertyExtractorMassHalo      (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_                           )"/>
    <referenceConstruct                             object="nodePropertyExtractorParentNode_"     constructor="nodePropertyExtractorDescendantNode(                                                            timeParent,cosmologyFunctions_   ,self%nodePropertyExtractorMassParent_                      )"/>
    <referenceConstruct                             object="nodePropertyExtractorMassRatio_"      constructor="nodePropertyExtractorRatio         ('massRatio','Ratio of progenitor and parent masses.',nodePropertyExtractorMassProgenitor_,     nodePropertyExtractorParentNode_                          )"/>
    !!]
    ! Create a distribution normalizer which normalizes to bin width.
    allocate(outputAnalysisDistributionNormalizerBinWidth_       )
    allocate(outputAnalysisDistributionNormalizerLog10ToLog_     )
    allocate(outputAnalysisDistributionNormalizer_               )
    allocate(normalizer_                                         )
    allocate(normalizer_                                    %next)
    normalizer_     %normalizer_ => outputAnalysisDistributionNormalizerBinWidth_
    normalizer_%next%normalizer_ => outputAnalysisDistributionNormalizerLog10ToLog_
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerBinWidth_"   constructor="outputAnalysisDistributionNormalizerBinWidth  (           )"/>
    <referenceConstruct object="outputAnalysisDistributionNormalizerLog10ToLog_" constructor="outputAnalysisDistributionNormalizerLog10ToLog(           )"/>
    <referenceConstruct object="outputAnalysisDistributionNormalizer_"           constructor="outputAnalysisDistributionNormalizerSequence  (normalizer_)"/>
    !!]
    ! Build log10() property operator.
    allocate(     outputAnalysisPropertyOperator_    )
    !![
    <referenceConstruct                             object="outputAnalysisPropertyOperator_"          constructor="outputAnalysisPropertyOperatorLog10             (                                                                                                                                        )"/>
    !!]
    ! Build anti-log10() property operator.
    allocate(     outputAnalysisPropertyUnoperator_  )
    !![
    <referenceConstruct                             object="outputAnalysisPropertyUnoperator_"        constructor="outputAnalysisPropertyOperatorAntiLog10         (                                                                                                                                        )"/>
    !!]
    ! Build an identity property operator.
    allocate(     outputAnalysisPropertyIdentity_    )
    !![
    <referenceConstruct                             object="outputAnalysisPropertyIdentity_"          constructor="outputAnalysisPropertyOperatorIdentity          (                                                                                                                                        )"/>
    !!]
    ! Build a weight operator which weights by the mass ratio.
    allocate(     outputAnalysisWeightOperatorMassRatio_)
    !![
    <referenceConstruct                             object="outputAnalysisWeightOperatorMassRatio_"   constructor="outputAnalysisWeightOperatorProperty            (                                         nodePropertyExtractorMassRatio_ ,outputAnalysisPropertyIdentity_                               )"/>
    !!]
    ! Build a weight operator for node subsampling.
    allocate(     outputAnalysisWeightOperatorSubsampling_)
    !![
    <referenceConstruct                             object="outputAnalysisWeightOperatorSubsampling_" constructor="outputAnalysisWeightOperatorSubsampling         (                                                                                                                                        )"/>
    !!]
    allocate(outputAnalysisWeightOperator_     )
    allocate(weightOperator_                   )
    allocate(weightOperator_              %next)
    weightOperator_     %operator_ => outputAnalysisWeightOperatorMassRatio_
    weightOperator_%next%operator_ => outputAnalysisWeightOperatorSubsampling_
    !![
    <referenceConstruct                             object="outputAnalysisWeightOperator_"            constructor="outputAnalysisWeightOperatorSequence            (weightOperator_                                                                                                                         )"/>
    !!]   
    ! Build a weight operator for the parent node mass uncertainty.
    allocate(self%outputAnalysisWeightOperatorNbodyMass_      )
    !![
    <referenceConstruct isResult="yes" owner="self" object="outputAnalysisWeightOperatorNbodyMass_"   constructor="outputAnalysisWeightOperatorNbodyMass           (massParentMinimum,massParentMaximum,self%nodePropertyExtractorMassParent_,outputAnalysisPropertyIdentity_,     nbodyHaloMassError_      )"/>
    !!]
    ! Build an identity distribution operator.
    allocate(     outputAnalysisDistributionOperator_)
    !![
    <referenceConstruct                             object="outputAnalysisDistributionOperator_"      constructor="outputAnalysisDistributionOperatorMassRatioNBody(massParentMinimum,massParentMaximum,     timeParent                      ,nbodyHaloMassError_            ,self%galacticFilterParentMass_)"/>
    !!]
    ! Determine number of buffer bins.
    bufferCount=0_c_size_t
    ! Construct the object.
    self%outputAnalysisVolumeFunction1D=                                                              &
         & outputAnalysisVolumeFunction1D(                                                            &
         &                                var_str('progenitorMassFunction')//label                  , &
         &                                comment                                                   , &
         &                                var_str('massRatio'                                      ), &
         &                                var_str('Mass ratio at the bin center'                   ), &
         &                                var_str('dimensionless'                                  ), &
         &                                0.0d0                                                     , &
         &                                var_str('progenitorMassFunction'                         ), &
         &                                var_str('Progenitor mass function averaged over each bin'), &
         &                                var_str('dimensionless'                                  ), &
         &                                0.0d0                                                     , &
         &                                log10(massRatios)                                         , &
         &                                bufferCount                                               , &
         &                                outputWeight                                              , &
         &                                nodePropertyExtractorMassRatio_                           , &
         &                                outputAnalysisPropertyOperator_                           , &
         &                                outputAnalysisPropertyUnoperator_                         , &
         &                                outputAnalysisWeightOperator_                             , &
         &                                outputAnalysisDistributionOperator_                       , &
         &                                outputAnalysisDistributionNormalizer_                     , &
         &                                galacticFilter_                                           , &
         &                                outputTimes_                                              , &
         &                                outputAnalysisCovarianceModelPoisson                      , &
         &                                covarianceBinomialBinsPerDecade                           , &
         &                                covarianceBinomialMassHaloMinimum                         , &
         &                                covarianceBinomialMassHaloMaximum                         , &
         &                                .false.                                                   , &
         &                                var_str('$x=M_\mathrm{progenitor}/M_\mathrm{parent}$'    ), &
         &                                var_str('$\mathrm{d}f/\mathrm{d}\log_\mathrm{e}x$'       ), &
         &                                .true.                                                    , &
         &                                .true.                                                    , &
         &                                targetLabel                                               , &
         &                                functionValueTarget                                       , &
         &                                functionCovarianceTarget                                    &
         &                               )
    !![
    <objectDestructor name="galacticFilterHaloIsolated_"                    />
    <objectDestructor name="galacticFilterProgenitorMass_"                  />
    <objectDestructor name="galacticFilterParentMassMinimum_"               />
    <objectDestructor name="galacticFilterParentMassMaximum_"               />
    <objectDestructor name="galacticFilterParentNode_"                      />
    <objectDestructor name="galacticFilterNot_"                             />
    <objectDestructor name="galacticFilter_"                                />
    <objectDestructor name="nodePropertyExtractorMassProgenitor_"           />
    <objectDestructor name="nodePropertyExtractorParentNode_"               />
    <objectDestructor name="nodePropertyExtractorMassRatio_"                />
    <objectDestructor name="outputAnalysisDistributionNormalizerBinWidth_"  />
    <objectDestructor name="outputAnalysisDistributionNormalizerLog10ToLog_"/>
    <objectDestructor name="outputAnalysisDistributionNormalizer_"          />
    <objectDestructor name="outputAnalysisPropertyIdentity_"                />
    <objectDestructor name="outputAnalysisPropertyOperator_"                />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"              />
    <objectDestructor name="outputAnalysisWeightOperatorMassRatio_"         />
    <objectDestructor name="outputAnalysisWeightOperatorSubsampling_"       />
    <objectDestructor name="outputAnalysisWeightOperator_"                  />
    <objectDestructor name="outputAnalysisDistributionOperator_"            />
    !!]
    if (self%alwaysIsolatedOnly) then
       !![
       <objectDestructor name="galacticFilterHaloAlwaysIsolated_"/>
       !!]
    end if
    nullify(filters_       )
    nullify(filtersParent_ )
    nullify(normalizer_    )
    nullify(weightOperator_)
    return
  end function progenitorMassFunctionConstructorInternal

  subroutine progenitorMassFunctionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily progenitorMassFunction} output analysis class.
    !!}
    implicit none
    type(outputAnalysisProgenitorMassFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilterParentMass_"             />
    <objectDestructor name="self%outputAnalysisWeightOperatorNbodyMass_"/>
    <objectDestructor name="self%nodePropertyExtractorMassParent_"      />
    <objectDestructor name="self%cosmologyFunctions_"                   />
    <objectDestructor name="self%cosmologyParameters_"                  />
    <objectDestructor name="self%virialDensityContrast_"                />
    <objectDestructor name="self%nbodyHaloMassError_"                   />
    <objectDestructor name="self%darkMatterProfileDMO_"                 />
    <objectDestructor name="self%virialDensityContrastDefinition_"      />
    !!]
    return
  end subroutine progenitorMassFunctionDestructor
  
  subroutine progenitorMassFunctionNewTree(self,tree,iOutput)
    !!{
    Record the weight of parent nodes.
    !!}
    use :: Galacticus_Nodes       , only : nodeComponentBasic
    use :: Merger_Tree_Walkers    , only : mergerTreeWalkerIsolatedNodes
    use :: Numerical_Comparison   , only : Values_Agree
    use :: Output_Analyses_Options, only : enumerationOutputAnalysisPropertyTypeType, enumerationOutputAnalysisPropertyQuantityType
    implicit none
    class           (outputAnalysisProgenitorMassFunction         ), intent(inout) :: self
    type            (mergerTree                                   ), intent(inout) :: tree
    integer         (c_size_t                                     ), intent(in   ) :: iOutput
    double precision                                               , parameter     :: timeTolerance=1.0d-4
    type            (treeNode                                     ), pointer       :: node
    class           (nodeComponentBasic                           ), pointer       :: basic
    type            (mergerTreeWalkerIsolatedNodes                )                :: treeWalker
    double precision                                                               :: weight                , mass
    type            (enumerationOutputAnalysisPropertyTypeType    )                :: propertyType 
    type            (enumerationOutputAnalysisPropertyQuantityType)                :: propertyQuantity

    ! Only accumulate tree weight if the output corresponds to the one for which we are constructing the conditional mass
    ! function.
    if (iOutput /= self%indexOutput) return
    ! Walk the tree, applying a filter to find parent nodes, and accumulate their weights.
    treeWalker=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       basic => node%basic()
              if (Values_Agree(basic%time(),self%timeParent,absTol=timeTolerance) .and. self%galacticFilterParentMass_%passes(node)) then
          weight            =+node%hostTree%volumeWeight
          mass              =+self%nodePropertyExtractorMassParent_      %extract      (       node                                                )
          propertyType      = self%nodePropertyExtractorMassParent_      %type         (                                                           )
          propertyQuantity  = self%nodePropertyExtractorMassParent_      %quantity     (                                                           )
          weight            =+self%outputAnalysisWeightOperatorNbodyMass_%operate      (weight,node,mass,mass,propertyType,propertyQuantity,iOutput)
          self%weightParents=+self                                       %weightParents                                                              &
               &             +                                            weight
       end if
    end do
    return
  end subroutine progenitorMassFunctionNewTree

  subroutine progenitorMassFunctionReduce(self,reduced)
    !!{
    Implement reduction over progenitor mass functions.
    !!}
    use    :: Error  , only : Error_Report
    !$ use :: OMP_Lib, only : OMP_Set_Lock, OMP_Unset_Lock
    implicit none
    class(outputAnalysisProgenitorMassFunction), intent(inout) :: self
    class(outputAnalysisClass                 ), intent(inout) :: reduced
    
    select type (reduced)
    class is (outputAnalysisProgenitorMassFunction)
       !$ call OMP_Set_Lock(reduced%accumulateLock)
       reduced%weightParents=+reduced%weightParents &
            &                +self   %weightParents
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    call self%outputAnalysisVolumeFunction1D%reduce(reduced)
    return
  end subroutine progenitorMassFunctionReduce

  subroutine progenitorMassFunctionFinalizeAnalysis(self)
    !!{
    Implement analysis finalization for progenitor mass functions. We simply normalize the accumulated weight of parent nodes.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
#endif
    implicit none
    class(outputAnalysisProgenitorMassFunction), intent(inout) :: self

    call self%outputAnalysisVolumeFunction1D%finalizeAnalysis()
    ! If already finalized, no need to do anything.
    if (self%weightsFinalized) return
    self%weightsFinalized=.true.
#ifdef USEMPI
    ! If running under MPI, perform a summation reduction of the parent weights across all processes.
    self%weightParents=mpiSelf%sum(self%weightParents)
#endif
    if (self%weightParents > 0.0d0) then
       self%functionValue     =+self%functionValue         &
            &                  /self%weightParents
       self%functionCovariance=+self%functionCovariance    &
            &                  /self%weightParents     **2
    end if
    return
  end subroutine progenitorMassFunctionFinalizeAnalysis

  subroutine progenitorMassFunctionFinalize(self,groupName)
    !!{
    Implement analysis finalization for progenitor mass functions. We simply normalize the accumulated weight of parent nodes.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(outputAnalysisProgenitorMassFunction), intent(inout)           :: self
    type (varying_string                      ), intent(in   ), optional :: groupName
    type (hdf5Object                          )               , target   :: analysesGroup, subGroup
    type (hdf5Object                          )               , pointer  :: inGroup
    type (hdf5Object                          )                          :: analysisGroup

    call self                               %finalizeAnalysis(         )
    call self%outputAnalysisVolumeFunction1D%finalize        (groupName)
    ! Add attributes giving the range of mass ratios considered. Also re-write the log-likelihood attribute here as it will have
    ! been written as part of the "outputAnalysisVolumeFunction1D" parent class, but we need to do our own calculation.    
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=analysesGroup%openGroup(char(self%label),char(self%comment))
    call analysisGroup%writeAttribute(self%logLikelihood             (),'logLikelihood'             )
    call analysisGroup%writeAttribute(self%massRatioLikelihoodMinimum  ,'massRatioLikelihoodMinimum')
    call analysisGroup%writeAttribute(self%massRatioLikelihoodMaximum  ,'massRatioLikelihoodMaximum')
    call analysisGroup%close         (                                                              )
    if (present(groupName)) &
         & call subGroup%close       (                                                              )
    call analysesGroup%close         (                                                              )
    !$ call hdf5Access%unset()
    return
  end subroutine progenitorMassFunctionFinalize

  double precision function progenitorMassFunctionLogLikelihood(self)
    !!{
    Return the log-likelihood of the progenitor mass function.
    !!}
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use            :: Interface_GSL               , only : GSL_Success
    use            :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (outputAnalysisProgenitorMassFunction), intent(inout)                 :: self
    double precision                                      , allocatable  , dimension(:,:) :: functionCovarianceCombined
    double precision                                      , allocatable  , dimension(:  ) :: functionValueDifference
    logical                                               , allocatable  , dimension(:  ) :: mask
    double precision                                      , parameter                     :: logRatioZero              =7.0d0
    type            (vector                              )                                :: residual
    type            (matrix                              )                                :: covariance
    integer         (c_size_t                            )                                :: i                               , j             , &
         &                                                                                   ii                              , jj
    integer                                                                               :: status
    double precision                                                                      :: covarianceTermTarget            , covarianceTerm
    
    ! Check for existence of a target distribution.
    if (allocated(self%functionValueTarget)) then
       ! Finalize analysis.
       call self%finalizeAnalysis()
       ! Find bins which satisfy the mass ratio limits and have a measured target value.
       mask   = self%binCenter           > self%massRatioLikelihoodMinimum &
            &  .and.                                                       &
            &   self%binCenter           < self%massRatioLikelihoodMaximum &
            &  .and.                                                       &
            &   self%functionValueTarget > 0.0d0
       if (count(mask) > 0) then
          ! Allocate workspaces.
          allocate(functionCovarianceCombined(count(mask),count(mask)))
          allocate(functionValueDifference   (count(mask)            ))
          ! Find combined covariance and difference between model and target.
          ii=0
          do i=1,self%binCount
             if (mask(i)) then
                ii=ii+1
                if (self%likelihoodInLog) then
                   !  Compute difference between model and target.
                   if (self%functionValue(i) > 0.0d0) then
                      ! Map values to compute difference in log(φ).
                      functionValueDifference(ii)=log(+self%functionValue      (i) &
                           &                          /self%functionValueTarget(i) &
                           &                         )
                   else
                      functionValueDifference(ii)=logRatioZero
                   end if
                else
                   ! Compute difference in φ.
                   functionValueDifference   (ii)=    +self%functionValue      (i) &
                        &                             -self%functionValueTarget(i)
                end if
                jj=0
                do j=1,self%binCount
                   if (mask(j)) then
                      jj=jj+1
                      ! Compute covariance terms for model and target
                      covarianceTermTarget=self%functionCovarianceTarget(i,j)
                      if (self%covarianceTargetOnly) then
                         covarianceTerm   =0.0d0
                      else
                         covarianceTerm   =self%functionCovariance      (i,j)
                      end if
                      ! Map to log(φ) if requested.
                      if (self%likelihoodInLog) then
                         covarianceTermTarget=+covarianceTermTarget        &
                              &               /self%functionValueTarget(i) &
                              &               /self%functionValueTarget(j)
                         if (self%functionValue(i) > 0.0d0 .and. self%functionValue(j) > 0.0d0) then
                            covarianceTerm   =+covarianceTerm              &
                                 &            /self%functionValue      (i) &
                                 &            /self%functionValue      (j)
                         else
                            covarianceTerm   =+0.0d0
                         end if
                      end if
                      ! Compute total covariance.
                      functionCovarianceCombined(ii,jj)=+covarianceTermTarget &
                           &                            +covarianceTerm
                      ! Set a floor in covariance.
                      if (self%likelihoodInLog) then
                         functionCovarianceCombined(ii,jj)=max(                                                                                                                           &
                              &                                +self%rootVarianceTargetFractional(min(i,size(self%rootVarianceTargetFractional)))                                         &
                              &                                *self%rootVarianceTargetFractional(min(j,size(self%rootVarianceTargetFractional)))                                       , &
                              &                                +                                                                                       functionCovarianceCombined(ii,jj)  &
                              &                               )
                      else
                         functionCovarianceCombined(ii,jj)=max(                                                                                                                           &
                              &                                +self%rootVarianceTargetFractional(min(i,size(self%rootVarianceTargetFractional)))*self%functionValueTarget       (i    )  &
                              &                                *self%rootVarianceTargetFractional(min(j,size(self%rootVarianceTargetFractional)))*self%functionValueTarget       (    j), &
                              &                                +                                                                                       functionCovarianceCombined(ii,jj)  &
                              &                               )
                      end if
                      ! Zero off-diagonal terms if requested.
                      if (self%covarianceDiagonalize .and. ii /= jj) &
                           & functionCovarianceCombined(ii,jj)=0.0d0
                   end if
                end do
             end if
          end do
          residual  =vector(functionValueDifference   )
          covariance=matrix(functionCovarianceCombined)
          ! Compute the log-likelihood.
          progenitorMassFunctionLogLikelihood=-0.5d0*covariance%covarianceProduct(residual,status)
          if (status /= GSL_Success) progenitorMassFunctionLogLikelihood=logImprobable
       else
          progenitorMassFunctionLogLikelihood=+0.0d0
       end if
    else
       progenitorMassFunctionLogLikelihood   =+0.0d0
    end if
    return
  end function progenitorMassFunctionLogLikelihood
