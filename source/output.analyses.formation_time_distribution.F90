!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Contains a module which implements a concentration distribution output analysis class for dark matter halo time formation distribution.
  !!}
  
  use :: Galactic_Filters                , only : galacticFilterAll
  use :: Node_Property_Extractors        , only : nodePropertyExtractorMassHalo

  !![
  <outputAnalysis name="outputAnalysisFormationTimeDistribution">
   <description>A dark matter halo time formation distribution  output analysis class.</description>
   <deepCopy>
    <functionClass variables="galacticFilterParentMass_, nodePropertyExtractorMassParent_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="galacticFilterParentMass_, nodePropertyExtractorMassParent_"/>
   </stateStorable>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisFormationTimeDistribution
     !!{
     A dark matter halo time formation distribution output analysis class.
     !!}
     private
     class           (cosmologyFunctionsClass              ), pointer                   :: cosmologyFunctions_                    => null()
     class           (cosmologyParametersClass             ), pointer                   :: cosmologyParameters_                   => null()
     class           (nbodyHaloMassErrorClass              ), pointer                   :: nbodyHaloMassError_                    => null()
     class           (virialDensityContrastClass           ), pointer                   :: virialDensityContrastDefinition_       => null(), virialDensityContrast_     => null()
     class           (darkMatterProfileDMOClass            ), pointer                   :: darkMatterProfileDMO_                  => null()
     type            (galacticFilterAll                    ), pointer                   :: galacticFilterParentMass_              => null()
     type            (nodePropertyExtractorMassHalo        ), pointer                   :: nodePropertyExtractorMassParent_       => null()
     double precision                                       , allocatable, dimension(:) :: rootVarianceTargetFractional
     double precision                                                                   :: massParentMinimum                               , massParentMaximum                   , &
          &                                                                                timeProgenitor                                  , timeParent                          , &
          &                                                                                redshiftProgenitor                              , redshiftParent                      , &
          &                                                                                weightParents                                   , redshiftMinimum                     , &
          &                                                                                redshiftMaximum
     integer                                                                            :: indexParent                                     , indexRedshift
     integer         (c_size_t                             )                            :: countRedshiftProgenitor                         , indexOutput
     logical                                                                            :: alwaysIsolatedOnly                              , covarianceDiagonalize               , &
          &                                                                                covarianceTargetOnly                            , weightsFinalized
     type            (varying_string                       )                            :: fileName
  contains
     final     ::                     formationTimeDistributionDestructor
     procedure :: reduce           => formationTimeDistributionReduce
     procedure :: finalizeAnalysis => formationTimeDistributionFinalizeAnalysis
  end type outputAnalysisFormationTimeDistribution

  interface outputAnalysisFormationTimeDistribution
     !!{
     Constructors for the ``formationTimeDistribution'' output analysis class.
     !!}
     module procedure formationTimeDistributionConstructorParameters
     module procedure formationTimeDistributionConstructorFile
     module procedure formationTimeDistributionConstructorInternal
  end interface outputAnalysisFormationTimeDistribution

contains
  
  function formationTimeDistributionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``formationTimeDistribution'' output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctionsClass
    use :: Input_Parameters                 , only : inputParameter            , inputParameters
    use :: ISO_Varying_String               , only : char
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast          , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisFormationTimeDistribution)                              :: self
    type            (inputParameters                     ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass             ), pointer                     :: cosmologyFunctions_
    class           (cosmologyParametersClass            ), pointer                     :: cosmologyParameters_
    class           (nbodyHaloMassErrorClass             ), pointer                     :: nbodyHaloMassError_
    class           (outputTimesClass                    ), pointer                     :: outputTimes_
    class           (darkMatterProfileDMOClass           ), pointer                     :: darkMatterProfileDMO_
    class           (virialDensityContrastClass          ), pointer                     :: virialDensityContrastDefinition_, virialDensityContrast_
    double precision                                      , dimension(:  ), allocatable :: functionValueTarget             , functionCovarianceTarget1D, &
         &                                                                                 rootVarianceTargetFractional
    double precision                                      , dimension(:,:), allocatable :: functionCovarianceTarget
    double precision                                                                    :: massParentMinimum               , massParentMaximum         , &                                                 
         &                                                                                 redshiftProgenitor              , redshiftParent            , &
         &                                                                                 redshiftMinimum                 , redshiftMaximum
    integer         (c_size_t                            )                              :: countRedshiftProgenitor
    integer                                                                             :: indexParent                     , indexRedshift
    type            (varying_string                      )                              :: label                           , comment                   , &
         &                                                                                 targetLabel                     , fileName
    logical                                                                             :: alwaysIsolatedOnly              , covarianceDiagonalize     , &
          &                                                                                covarianceTargetOnly
    
    allocate(rootVarianceTargetFractional(max(1,parameters%count('rootVarianceTargetFractional',zeroIfNotPresent=.true.))))
    !![
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"             source="parameters"                                                />
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"              source="parameters"                                                />
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"            source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_"           source="parameters"                                                />
    <objectBuilder class="nbodyHaloMassError"    name="nbodyHaloMassError_"              source="parameters"                                                /> 
    <objectBuilder class="outputTimes"           name="outputTimes_"                     source="parameters"                                                />
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
    !!]
    if (parameters%isPresent('rootVarianceTargetFractional')) then
       if (parameters%count('rootVarianceTargetFractional') == 1) then
          ! For a single value, read as a scalar to allow processing of math in the parameter value.
          !![
	  <inputParameter>	    
	    <name>rootVarianceTargetFractional</name>
	    <variable>rootVarianceTargetFractional(1)</variable>
	    <source>parameters</source>
	    <description>The diagonal of the covariance matrix is forced to be at least equal to this fraction multiplied by the target dataset squared. This may be a list of values corresponding to each element of the target dataset. If the list is shorter than the target dataset the final value in the list is applied to all remaining elements in the target dataset.</description>
	  </inputParameter>
          !!]
       else
          !![
	  <inputParameter>
	    <name>rootVarianceTargetFractional</name>
	    <source>parameters</source>
	    <description>The diagonal of the covariance matrix is forced to be at least equal to this fraction multiplied by the target dataset squared. This may be a list of values corresponding to each element of the target dataset. If the list is shorter than the target dataset the final value in the list is applied to all remaining elements in the target dataset.</description>
	  </inputParameter>
          !!]
       end if
    else
       rootVarianceTargetFractional=0.0d0
    end if
    !![
    <inputParameter>
      <name>redshiftMinimum</name>
      <source>parameters</source>
      <description>The minimum redshift to include in redshift formation calculations.</description>
      <defaultValue>0.0d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <source>parameters</source>
      <description>The maximum redshift to include in redshift formation calculations.</description>
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
         <description>The name of the file from which to read redshift formation time  parameters.</description>
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
       self=outputAnalysisFormationTimeDistribution(char(fileName),label,comment,targetLabel,redshiftParent,indexParent,indexRedshift,covarianceDiagonalize,covarianceTargetOnly,rootVarianceTargetFractional,alwaysIsolatedOnly,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_,nbodyHaloMassError_,outputTimes_)
    else
       !![
       <inputParameter>
         <name>label</name>
         <source>parameters</source>
         <variable>label</variable>
         <description>A label for the redshift formation time distribution.</description>
       </inputParameter>
       <inputParameter>
         <name>comment</name>
         <source>parameters</source>
         <variable>comment</variable>
         <description>A descriptive comment for the redshift formation time distribution.</description>
       </inputParameter>
       <inputParameter>
         <name>redshiftMinimum</name>
         <source>parameters</source>
         <description>Minimum redshift for the redshift formation time distribution.</description>
       </inputParameter>
       <inputParameter>
         <name>redshiftMaximum</name>
         <source>parameters</source>
         <description>Maximum redshift for the redshift formation time distribution.</description>
       </inputParameter>
       <inputParameter>
         <name>countRedshiftProgenitor</name>
         <source>parameters</source>
         <description>Number of redshift of progenitors at which to compute the redshift formation time distribution.</description>
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
          self=outputAnalysisFormationTimeDistribution(                                                                                                 &amp;
          &amp;                                    label                                                                                              , &amp;
          &amp;                                    comment                                                                                            , &amp;
          &amp;                                    redshiftMinimum                                                                                    , &amp;
          &amp;                                    redshiftMaximum                                                                                    , &amp;
          &amp;                                    countRedshiftProgenitor                                                                            , &amp;
          &amp;                                    massParentMinimum                                                                                  , &amp;
          &amp;                                    massParentMaximum                                                                                  , &amp;
          &amp;                                    cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftProgenitor)), &amp;
          &amp;                                    cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftParent    )), &amp;
          &amp;                                    alwaysIsolatedOnly                                                                                 , &amp;
          &amp;                                    covarianceDiagonalize                                                                              , &amp;
          &amp;                                    covarianceTargetOnly                                                                               , &amp;
          &amp;                                    rootVarianceTargetFractional                                                                       , &amp;
          &amp;                                    cosmologyParameters_                                                                               , &amp;
          &amp;                                    cosmologyFunctions_                                                                                , &amp;
          &amp;                                    darkMatterProfileDMO_                                                                              , &amp;
          &amp;                                    virialDensityContrast_                                                                             , &amp;
          &amp;                                    virialDensityContrastDefinition_                                                                   , &amp;
          &amp;                                    nbodyHaloMassError_                                                                                , &amp;
          &amp;                                    outputTimes_                                                                                         &amp;
          &amp;                                    {conditions}                                                                                         &amp;
          <!-- &amp;                                    targetLabel                                                                                        , &amp; -->
         <!-- &amp;                                    functionValueTarget                                                                                , &amp; -->
         <!-- &amp;                                    functionCovarianceTarget                                                                             &amp; -->
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
    <objectDestructor name="darkMatterProfileDMO_"           />
    <objectDestructor name="virialDensityContrast_"          />
    <objectDestructor name="outputTimes_"                    />
    <objectDestructor name="nbodyHaloMassError_"             />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function formationTimeDistributionConstructorParameters
  
  function formationTimeDistributionConstructorFile(fileName,label,comment,targetLabel,redshiftParent,indexParent,indexRedshift,covarianceDiagonalize,covarianceTargetOnly,rootVarianceTargetFractional,alwaysIsolatedOnly,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_,nbodyHaloMassError_,outputTimes_) result(self)
    !!{
    Constructor for the ``formationTimeDistribution'' output analysis class which reads all required properties from file.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctionsClass
    use :: HDF5_Access                      , only : hdf5Access
    use :: IO_HDF5                          , only : hdf5Object
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass
    use :: File_Utilities                   , only : File_Name_Expand
    use :: Virial_Density_Contrast          , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisFormationTimeDistribution)                                               :: self
    character       (len=*                               ), intent(in   )                                :: fileName
    type            (varying_string                      ), intent(in   )                                :: label                          , targetLabel               , &
         &                                                                                                  comment
    integer                                               , intent(in   )                                :: indexParent                    , indexRedshift
    double precision                                      , intent(in   )                                :: redshiftParent
    !--double precision                                      , intent(in   )                                :: redshiftMinimum                , redshiftMaximum           , redshiftParent --/>
    double precision                                      , intent(in   ), allocatable, dimension(:    ) :: rootVarianceTargetFractional
    logical                                               , intent(in   )                                :: covarianceDiagonalize          , covarianceTargetOnly      , &
         &                                                                                                  alwaysIsolatedOnly
    class           (cosmologyFunctionsClass             ), intent(inout), target                        :: cosmologyFunctions_
    class           (cosmologyParametersClass            ), intent(inout), target                        :: cosmologyParameters_
    class           (darkMatterProfileDMOClass           ), intent(in   ), target                        :: darkMatterProfileDMO_
    class           (outputTimesClass                    ), intent(inout), target                        :: outputTimes_
    class           (virialDensityContrastClass          ), intent(in   ), target                        :: virialDensityContrastDefinition_, virialDensityContrast_
    class           (nbodyHaloMassErrorClass             ), intent(in   ), target                        :: nbodyHaloMassError_
    double precision                                                                                     :: massParentMinimum               , massParentMaximum        , &
         &                                                                                                  timeProgenitor                  , timeParent               , &
         &                                                                                                  redshiftMinimum                 , redshiftMaximum
    double precision                                                     , allocatable, dimension(:    ) :: functionValueTarget             , massRatio                , &
         &                                                                                                  redshiftProgenitor_val          , massParents              , &
         &                                                                                                  massParentsMinimum              , massParentsMaximum
    double precision                                                                                     :: redshiftProgenitor
    double precision                                                     , allocatable, dimension(:,:  ) :: functionCovarianceTarget
    double precision                                                     , allocatable, dimension(:,:,:) :: functionValuesTarget
    integer         (c_size_t                            )               , allocatable, dimension(:,:,:) :: functionCountsTarget
    type            (hdf5Object                          )                                               :: dataFile                        , simulationGroup
    integer                                                                                              :: i
    logical                                                                                              :: haveBoundaries
    
    !$ call hdf5Access%set  ()
    call dataFile%openFile(char(File_Name_Expand(fileName)),readOnly=.true.)
    simulationGroup=dataFile       %openGroup ('simulation0001/timeFormation'   )
    haveBoundaries =simulationGroup%hasDataset('massParentMinimum')
    call    simulationGroup%readDataset('redshift'              ,redshiftProgenitor_val)
    call    simulationGroup%readDataset('distribution'          ,functionValuesTarget)
    call    simulationGroup%readDataset('count'                 ,functionCountsTarget)
    call    simulationGroup%readDataset('massParent'            ,massParents         )
    !call    simulationGroup%readDataset('redshiftProgenitor'    ,redshiftProgenitors )
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
    !Extracting redshift min and max
    redshiftMinimum =redshiftProgenitor_val(1)
    redshiftMaximum =redshiftProgenitor_val(size(redshiftProgenitor_val))
    redshiftProgenitor  =redshiftProgenitor_val(indexRedshift+1)
    !redshiftParent  =redshiftProgenitor_val(indexParent+1)
    ! Extract the target function values.
    allocate(functionValueTarget(size(redshiftProgenitor_val)))
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
    self              =outputAnalysisFormationTimeDistribution(label,comment,redshiftMinimum,redshiftMaximum,size(redshiftProgenitor_val,kind=c_size_t),massParentMinimum,massParentMaximum,timeProgenitor,timeParent,alwaysIsolatedOnly,covarianceDiagonalize,covarianceTargetOnly,rootVarianceTargetFractional,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_,nbodyHaloMassError_,outputTimes_,targetLabel,functionValueTarget,functionCovarianceTarget)
    !![
    <constructorAssign variables="fileName, indexParent, indexRedshift"/>
    !!]
    return
  end function formationTimeDistributionConstructorFile

  function formationTimeDistributionConstructorInternal(label,comment,redshiftMinimum,redshiftMaximum,countRedshiftProgenitor,massParentMinimum,massParentMaximum,timeProgenitor,timeParent,alwaysIsolatedOnly,covarianceDiagonalize,covarianceTargetOnly,rootVarianceTargetFractional,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_,nbodyHaloMassError_,outputTimes_,targetLabel,functionValueTarget,functionCovarianceTarget) result(self)
    !!{
    Internal constructor for the ``formationTimeDistribution'' output analysis class.
    !!}
    use :: HDF5_Access                             , only : hdf5Access
    use :: Galactic_Filters                        , only : filterList                                      , galacticFilterDescendantNode                , galacticFilterHaloAlwaysIsolated             ,&
      &                                                     galacticFilterHaloIsolated                      , galacticFilterHaloMass                      , galacticFilterNot                            ,&
      &                                                     galacticFilterHaloMassRange
    use :: Node_Property_Extractors                , only : nodePropertyExtractorDescendantNode             , nodePropertyExtractorNodeFormationTime
    use :: Numerical_Comparison                    , only : Values_Agree
    use :: Numerical_Ranges                        , only : Make_Range                                      , rangeTypeLinear
    use :: Cosmology_Functions                     , only : cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                    , only : cosmologyParametersSimple
    use :: Dark_Matter_Profiles_DMO                , only : darkMatterProfileDMOClass
    use :: Error                                   , only : Error_Report
    use :: Output_Analysis_Distribution_Normalizers, only : normalizerList                                  , outputAnalysisDistributionNormalizerUnitarity , outputAnalysisDistributionNormalizerBinWidth,&
      &                                                     outputAnalysisDistributionNormalizerSequence
    use :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorIdentity
    use :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorSubsampling         , weightOperatorList
    use :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorIdentity
    use :: Output_Analyses_Options                 , only : outputAnalysisCovarianceModelPoisson
    use :: Statistics_NBody_Halo_Mass_Errors       , only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast                 , only : virialDensityContrastClass
    implicit none
    type            (outputAnalysisFormationTimeDistribution         )                                          :: self
    type            (varying_string                                  ), intent(in   )                           :: label                                                  , comment
    double precision                                                  , intent(in   )                           :: massParentMinimum                                      , massParentMaximum     , &
         &                                                                                                         timeProgenitor                                         , timeParent
    double precision                                                  , intent(in   )          , dimension(:  ) :: rootVarianceTargetFractional
    integer         (c_size_t                                        ), intent(in   )                           :: countRedshiftProgenitor
    logical                                                           , intent(in   )                           :: alwaysIsolatedOnly                                     , covarianceDiagonalize                   , &
         &                                                                                                         covarianceTargetOnly
    double precision                                                  , intent(in   )                           :: redshiftMinimum                                        , redshiftMaximum
    class           (cosmologyParametersClass                        ), intent(inout), target                   :: cosmologyParameters_
    class           (cosmologyFunctionsClass                         ), intent(inout), target                   :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass                       ), intent(in   ), target                   :: darkMatterProfileDMO_
    class           (outputTimesClass                                ), intent(inout), target                   :: outputTimes_
    class           (nbodyHaloMassErrorClass                         ), intent(in   ), target                   :: nbodyHaloMassError_
    class           (virialDensityContrastClass                      ), intent(in   ), target                   :: virialDensityContrastDefinition_                       , virialDensityContrast_
    type            (varying_string                                  ), intent(in   ), optional                 :: targetLabel
    double precision                                                  , intent(in   ), optional, dimension(:  ) :: functionValueTarget
    double precision                                                  , intent(in   ), optional, dimension(:,:) :: functionCovarianceTarget
    double precision                                                  , parameter                               :: timeTolerance                                  =1.0d-04
    !double precision                                                  , parameter                               :: massRatioBuffer                                =1.0d-01
    integer                                                           , parameter                               :: covarianceBinomialBinsPerDecade                =2
    double precision                                                  , parameter                               :: covarianceBinomialMassHaloMinimum              =3.0d+11, covarianceBinomialMassHaloMaximum=1.0d15
    logical                                                           , parameter                               :: allowSelf                                      =.false.
    double precision                                                  , allocatable            , dimension(:  ) :: redshiftProgenitors
    double precision                                                  , allocatable            , dimension(:,:) :: outputWeight
    type            (galacticFilterAll                               ), pointer                                 :: galacticFilter_
    type            (galacticFilterHaloIsolated                      ), pointer                                 :: galacticFilterHaloIsolated_
    type            (galacticFilterDescendantNode                    ), pointer                                 :: galacticFilterParentNode_
    type            (galacticFilterNot                               ), pointer                                 :: galacticFilterNot_
    !type            (galacticFilterHaloMass                          ), pointer                                 :: galacticFilterParentMassMinimum_                      , galacticFilterParentMassMaximum_
    type            (galacticFilterHaloMassRange                     ), pointer                                 :: galacticFilterParentMassMinimum_                     ,galacticFilterParentMassMaximum_,&
       &                                                                                                           galacticFilterHaloMassRange_ 
    type            (galacticFilterHaloAlwaysIsolated                ), pointer                                 :: galacticFilterHaloAlwaysIsolated_
    type            (filterList                                      ), pointer                                 :: filters_                                               , filtersParent_
    type            (nodePropertyExtractorMassHalo                   ), pointer                                 :: nodePropertyExtractorMassProgenitor_
    type            (nodePropertyExtractorDescendantNode             ), pointer                                 :: nodePropertyExtractorParentNode_
    type            (nodePropertyExtractorNodeFormationTime          ), pointer                                 :: nodePropertyExtractorNodeFormationTime_
    type            (outputAnalysisDistributionNormalizerSequence    ), pointer                                 :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionNormalizerUnitarity   ), pointer                                 :: outputAnalysisDistributionNormalizerUnitarity_
    type            (outputAnalysisDistributionNormalizerBinWidth    ), pointer                                 :: outputAnalysisDistributionNormalizerBinWidth_
    type            (normalizerList                                  ), pointer                                 :: normalizer_
    type            (weightOperatorList                              ), pointer                                 :: weightOperator_
    type            (outputAnalysisWeightOperatorSubsampling         ), pointer                                 :: outputAnalysisWeightOperatorSubsampling_
    type            (outputAnalysisPropertyOperatorIdentity          ), pointer                                 :: outputAnalysisPropertyOperatorIdentity_
    type            (outputAnalysisDistributionOperatorIdentity      ), pointer                                 :: outputAnalysisDistributionOperatorIdentity_
    integer         (c_size_t                                        )                                          :: iOutput                                                , bufferCount
    type            (varying_string                                  )                                          :: message
    character       (len=10                                          )                                          :: timeLabel
    !![
    <constructorAssign variables="redshiftMinimum, redshiftMaximum, countRedshiftProgenitor, massParentMinimum, massParentMaximum, timeProgenitor, timeParent, alwaysIsolatedOnly, covarianceDiagonalize, covarianceTargetOnly, rootVarianceTargetFractional, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *virialDensityContrast_, *virialDensityContrastDefinition_, *nbodyHaloMassError_, *outputTimes_"/>
    !!]

    ! Initialize state.
    self%weightsFinalized=.false.
    ! Set redshifts.
    self%redshiftProgenitor=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeProgenitor))
    self%redshiftParent    =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeParent    ))
    ! Build a grid of formation redshifts
    allocate(redshiftProgenitors(countRedshiftProgenitor))
    redshiftProgenitors=Make_Range(redshiftMinimum, redshiftMaximum, int(countRedshiftProgenitor),rangeType=rangeTypeLinear) 
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(countRedshiftProgenitor,outputTimes_%count()))
    outputWeight=0.0d0
    self%indexOutput=-1_c_size_t
    do iOutput=1_c_size_t,outputTimes_%count()
       if (Values_Agree(outputTimes_%time(iOutput),timeProgenitor,absTol=timeTolerance)) then
          outputWeight(:,iOutput)=1.0d0
          self%indexOutput=iOutput
       end if
    end do
    if (self%indexOutput < 0_c_size_t) then
       message='no matching output time found - seeking '
       write (timeLabel,'(f9.3)') timeProgenitor
       message=message//trim(adjustl(timeLabel))//' in ['
       do iOutput=1_c_size_t,outputTimes_%count()
          write (timeLabel,'(f9.3)') outputTimes_%time(iOutput)
          message=message//' '//trim(adjustl(timeLabel))
          if (iOutput < outputTimes_%count()) message=message//','
       end do
       message=message//' ]'
       call Error_Report(message//{introspection:location})
    end if
    ! Initialize accumulated weight of parent nodes.
    self%weightParents=0.0d0
    ! Build a filter which selects isolated halos, above a suitable lower mass, and with parents in the correct mass range.
    allocate(     galacticFilter_                           )
    allocate(     galacticFilterHaloIsolated_               )
    allocate(     galacticFilterParentMassMinimum_          )
    allocate(     galacticFilterParentMassMaximum_          )
    allocate(self%galacticFilterParentMass_                 )
    allocate(     galacticFilterNot_                        )
    allocate(     galacticFilterParentNode_                 )
    allocate(     filters_                                  )
    allocate(     filters_                        %next     )
    allocate(     filtersParent_                            )
    allocate(     filtersParent_                  %next     )
    filters_                %filter_ => galacticFilterHaloIsolated_  
    filters_           %next%filter_ => galacticFilterParentNode_
    filtersParent_          %filter_ => galacticFilterParentMassMinimum_
    filtersParent_     %next%filter_ => galacticFilterNot_
    if (self%alwaysIsolatedOnly) then
       allocate(galacticFilterHaloAlwaysIsolated_               )
       allocate(filters_                         %next%next     )
       allocate(filtersParent_                   %next%next     )
       filters_           %next%next%filter_ => galacticFilterHaloAlwaysIsolated_
       filtersParent_     %next%next%filter_ => galacticFilterHaloAlwaysIsolated_
       !![
       <referenceConstruct object="galacticFilterHaloAlwaysIsolated_" constructor="galacticFilterHaloAlwaysIsolated()"/>
       !!]
    else
       nullify(galacticFilterHaloAlwaysIsolated_)
    end if
    
    ! Build a filter which selects isolated halos, above a suitable lower mass for formation time distribution. 
    allocate(     galacticFilterHaloIsolated_               )
    allocate(     galacticFilterHaloMassRange_              )
    allocate(     filters_                                  )
    allocate(     filters_                        %next     )
    filters_                %filter_ => galacticFilterHaloIsolated_
    filters_      %next     %filter_ => galacticFilterHaloMassRange_
    !![
    <referenceConstruct                             object="galacticFilter_"                  constructor="galacticFilterAll(filters_)"/>
    !!]

    !![
    <referenceConstruct                             object="galacticFilterHaloIsolated_"      constructor="galacticFilterHaloIsolated  (                                                                                           )"/>
    <referenceConstruct                             object="galacticFilterHaloMassRange_"     constructor="galacticFilterHaloMassRange (massParentMinimum,massParentMaximum,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_)"/>
    <referenceConstruct isResult="yes" owner="self" object="galacticFilterParentMass_"        constructor="galacticFilterAll           (filtersParent_                                                                                                 )"/>
    <referenceConstruct                             object="galacticFilterParentNode_"        constructor="galacticFilterDescendantNode(timeParent                                        ,allowSelf,cosmologyFunctions_,self%galacticFilterParentMass_)"/>
    !!]
    ! Build a node property extractor which gives the ratio of the progenitor and parent halo masses.
    allocate(     nodePropertyExtractorMassProgenitor_)
    allocate(self%nodePropertyExtractorMassParent_    )
    allocate(     nodePropertyExtractorParentNode_    )
    !![
    <referenceConstruct                             object="nodePropertyExtractorMassProgenitor_" constructor="nodePropertyExtractorMassHalo      (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_         )"/>
    <referenceConstruct isResult="yes" owner="self" object="nodePropertyExtractorMassParent_"     constructor="nodePropertyExtractorMassHalo      (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_         )"/>
    <referenceConstruct                             object="nodePropertyExtractorParentNode_"     constructor="nodePropertyExtractorDescendantNode(                                                            timeParent,cosmologyFunctions_   ,self%nodePropertyExtractorMassParent_    )"/>
    !!]
    ! Build a node property extractor for distribution of formation time
    allocate(nodePropertyExtractorNodeFormationTime_)
    !![
    <referenceConstruct                             object="nodePropertyExtractorNodeFormationTime_"               constructor="nodePropertyExtractorNodeFormationTime()"/>
    !!]
    ! Create a distribution normalizer which normalizes to bin width.
    allocate(outputAnalysisDistributionNormalizerUnitarity_       )
    allocate(outputAnalysisDistributionNormalizerBinWidth_       )
    allocate(outputAnalysisDistributionNormalizer_               )
    allocate(normalizer_                                         )
    allocate(normalizer_                                    %next)
    normalizer_     %normalizer_ => outputAnalysisDistributionNormalizerUnitarity_
    normalizer_%next%normalizer_ => outputAnalysisDistributionNormalizerBinWidth_
    !![
    <referenceConstruct object="outputAnalysisDistributionNormalizerUnitarity_"  constructor="outputAnalysisDistributionNormalizerUnitarity (            )"/>
    <referenceConstruct object="outputAnalysisDistributionNormalizerBinWidth_"   constructor="outputAnalysisDistributionNormalizerBinWidth  (           )"/>
    <referenceConstruct object="outputAnalysisDistributionNormalizer_"           constructor="outputAnalysisDistributionNormalizerSequence  (normalizer_)"/>
    !!]
    ! Build log10() property operator identity --keep unchanged
    allocate(     outputAnalysisPropertyOperatorIdentity_  )
    !![
    <referenceConstruct                             object="outputAnalysisPropertyOperatorIdentity_"  constructor="outputAnalysisPropertyOperatorIdentity          (                                                                                                                                        )"/>
    !!]
    ! Build outputAnalysisDistributionOperatorIdentity
    allocate(     outputAnalysisDistributionOperatorIdentity_  )
    !![
    <referenceConstruct                             object="outputAnalysisDistributionOperatorIdentity_"  constructor="outputAnalysisDistributionOperatorIdentity( )"/>
    !!]
    ! Build a weight operator that accounts for subsampling weights.
    allocate(     outputAnalysisWeightOperatorSubsampling_ )
    !![
    <referenceConstruct                             object="outputAnalysisWeightOperatorSubsampling_" constructor="outputAnalysisWeightOperatorSubsampling         (                                                                                                                                        )"/>
    !!]
    ! Determine number of buffer bins.
    bufferCount=0_c_size_t
    ! Construct the object.
    self%outputAnalysisVolumeFunction1D=                                                                 &
         & outputAnalysisVolumeFunction1D(                                                               &
         &                                var_str('formationTimeDistribution')//label                  , &
         &                                comment                                                      , &
         &                                var_str('redshiftFormation'                                 ), &
         &                                var_str('Formation redshift at the bin center'              ), &
         &                                var_str('dimensionless'                                     ), &
         &                                0.0d0                                                        , &
         &                                var_str('formationTimeDistribution'                         ), &
         &                                var_str('Distribution formation time averaged over each bin'), &
         &                                var_str('dimensionless'                                     ), &
         &                                0.0d0                                                        , &
         &                                redshiftProgenitors                                          , &
         &                                bufferCount                                                  , &
         &                                outputWeight                                                 , &
         &                                nodePropertyExtractorNodeFormationTime_                      , &
         &                                outputAnalysisPropertyOperatorIdentity_                      , &
         &                                outputAnalysisPropertyOperatorIdentity_                      , &
         &                                outputAnalysisWeightOperatorSubsampling_                     , &
         &                                outputAnalysisDistributionOperatorIdentity_                  , &
         &                                outputAnalysisDistributionNormalizer_                        , &
         &                                galacticFilter_                                              , &
         &                                outputTimes_                                                 , &
         &                                outputAnalysisCovarianceModelPoisson                         , &
         &                                covarianceBinomialBinsPerDecade                              , &
         &                                covarianceBinomialMassHaloMinimum                            , &
         &                                covarianceBinomialMassHaloMaximum                            , &
         &                                .false.                                                      , &
         &                                var_str('$x= z $'                                           ), &
         &                                var_str('Distribution'                                      ), &
         &                                .false.                                                      , &
         &                                .true.                                                       ,                                 targetLabel                                                  , &
         &                                functionValueTarget                                          , &
         &                                functionCovarianceTarget                                       &
         &                               )
    !![
    <objectDestructor name="galacticFilterHaloIsolated_"                    />
    <objectDestructor name="galacticFilterHaloMassRange_"                    />
    <objectDestructor name="galacticFilterParentMassMinimum_"               />
    <objectDestructor name="galacticFilterParentMassMaximum_"               />
    <objectDestructor name="galacticFilterParentNode_"                      />
    <objectDestructor name="galacticFilterNot_"                             />
    <objectDestructor name="galacticFilter_"                                />
    <objectDestructor name="nodePropertyExtractorNodeFormationTime_"         />
    <objectDestructor name="nodePropertyExtractorMassProgenitor_"           />
    <objectDestructor name="nodePropertyExtractorParentNode_"               />
    <objectDestructor name="outputAnalysisDistributionNormalizerUnitarity_" />
    <objectDestructor name="outputAnalysisDistributionNormalizerBinWidth_"  />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"          />
    <objectDestructor name="outputAnalysisDistributionOperatorIdentity_"    />
    <objectDestructor name="outputAnalysisPropertyOperatorIdentity_"        />
    <objectDestructor name="outputAnalysisWeightOperatorSubsampling_"       />
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
  end function formationTimeDistributionConstructorInternal

  subroutine formationTimeDistributionDestructor(self)
    !!{
    Destructor for the ``formationTimeDistribution'' output analysis class.
    !!}
    implicit none
    type(outputAnalysisFormationTimeDistribution), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilterParentMass_"             />
    <objectDestructor name="self%nodePropertyExtractorMassParent_"      />
    <objectDestructor name="self%cosmologyFunctions_"                   />
    <objectDestructor name="self%cosmologyParameters_"                  />
    <objectDestructor name="self%darkMatterProfileDMO_"                 />
    <objectDestructor name="self%virialDensityContrast_"                />
    <objectDestructor name="self%nbodyHaloMassError_"                   />
    <objectDestructor name="self%virialDensityContrastDefinition_"      />
    !!]
    return
  end subroutine formationTimeDistributionDestructor
  
  subroutine formationTimeDistributionReduce(self,reduced)
    !!{
    Implement reduction over progenitor mass functions.
    !!}
    use    :: Error  , only : Error_Report
    !$ use :: OMP_Lib, only : OMP_Set_Lock, OMP_Unset_Lock
    implicit none
    class(outputAnalysisFormationTimeDistribution), intent(inout) :: self
    class(outputAnalysisClass                    ), intent(inout) :: reduced
    
    select type (reduced)
    class is (outputAnalysisFormationTimeDistribution)
       !$ call OMP_Set_Lock(reduced%accumulateLock)
       reduced%weightParents=+reduced%weightParents &
            &                +self   %weightParents
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    call self%outputAnalysisVolumeFunction1D%reduce(reduced)
    return
  end subroutine formationTimeDistributionReduce

  subroutine formationTimeDistributionFinalizeAnalysis(self)
    !!{
    Implement analysis finalization for formation time distribution. We simply normalize the accumulated weight of parent nodes.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
#endif
    implicit none
    class(outputAnalysisFormationTimeDistribution), intent(inout) :: self

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
  end subroutine formationTimeDistributionFinalizeAnalysis
