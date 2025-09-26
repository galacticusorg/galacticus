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
  Implements a star forming main sequence output analysis class.
  !!}

  use :: Cosmology_Functions                       , only : cosmologyFunctionsClass
  use :: Geometry_Surveys                          , only : surveyGeometryClass
  use :: Star_Formation_Rates_Disks                , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids            , only : starFormationRateSpheroidsClass
  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass

  !![
  <outputAnalysis name="outputAnalysisStarFormingMainSequence">
   <description>A star forming main sequence output analysis class.</description>
   <runTimeFileDependencies paths="fileName"/>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisStarFormingMainSequence
     !!{
     A star forming main sequence output analysis class.
     !!}
     private
     class           (surveyGeometryClass                      ), pointer :: surveyGeometry_                       => null()
     class           (cosmologyFunctionsClass                  ), pointer :: cosmologyFunctions_                   => null(), cosmologyFunctionsData => null()
     class           (starFormationRateDisksClass              ), pointer :: starFormationRateDisks_               => null()
     class           (starFormationRateSpheroidsClass          ), pointer :: starFormationRateSpheroids_           => null()
     class           (starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     type            (varying_string                           )          :: fileName
     double precision                                                     :: massMinimum                                    , massMaximum                     , &
          &                                                                  countMassesPerDecade     
  end type outputAnalysisStarFormingMainSequence

  interface outputAnalysisStarFormingMainSequence
     !!{
     Constructors for the \refClass{outputAnalysisStarFormingMainSequence} output analysis class.
     !!}
     module procedure starFormingMainSequenceConstructorParameters
     module procedure starFormingMainSequenceConstructorFile
     module procedure starFormingMainSequenceConstructorInternal
  end interface outputAnalysisStarFormingMainSequence

contains

  function starFormingMainSequenceConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisStarFormingMainSequence} output analysis class which takes a parameter set as input.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: Numerical_Ranges  , only : Make_Range    , rangeTypeLogarithmic
    implicit none
    type            (outputAnalysisStarFormingMainSequence    )                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    class           (galacticFilterClass                      ), pointer                     :: galacticFilter_
    class           (surveyGeometryClass                      ), pointer                     :: surveyGeometry_
    class           (cosmologyFunctionsClass                  ), pointer                     :: cosmologyFunctions_                  , cosmologyFunctionsData
    class           (outputTimesClass                         ), pointer                     :: outputTimes_
    class           (starFormationRateDisksClass              ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass          ), pointer                     :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass), pointer                     :: starFormationRateNuclearStarClusters_
    class           (outputAnalysisPropertyOperatorClass      ), pointer                     :: outputAnalysisPropertyOperator_      , outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass  ), pointer                     :: outputAnalysisDistributionOperator_
    double precision                                           , dimension(:  ), allocatable :: meanValueTarget                      , meanCovarianceTarget1D               , &
         &                                                                                      massesStellar
    double precision                                           , dimension(:,:), allocatable :: meanCovarianceTarget
    double precision                                                                         :: massMinimum                          , massMaximum                          , &
         &                                                                                      countMassesPerDecade
    type            (inputParameters                          )                              :: dataAnalysisParameters
    integer         (c_size_t                                 )                              :: countMasses
    type            (varying_string                           )                              :: targetLabel                          , fileName                             , &
         &                                                                                      label                                , comment

    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    !![
    <objectBuilder class="galacticFilter"                       name="galacticFilter_"                       source="parameters"                                                                 />
    <objectBuilder class="cosmologyFunctions"                   name="cosmologyFunctions_"                   source="parameters"                                                                 />
    <objectBuilder class="cosmologyFunctions"                   name="cosmologyFunctionsData"                source="dataAnalysisParameters"                                                     />
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"                                                                 />
    <objectBuilder class="surveyGeometry"                       name="surveyGeometry_"                       source="parameters"                                                                 />
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"                                                                 />
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"                                                                 />
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"                                                                 />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"       source="parameters"                                                                 />
    <objectBuilder class="outputAnalysisDistributionOperator"   name="outputAnalysisDistributionOperator_"   source="parameters"                                                                 />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisWeightPropertyOperator_" source="parameters"             parameterName="outputAnalysisWeightPropertyOperator"/>
    !!]
    if (parameters%isPresent('fileName')) then
       !![
       <inputParameter>
         <name>fileName</name>
         <source>parameters</source>
         <description>The name of the file from which to read star forming main sequence function parameters.</description>
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
       self=outputAnalysisStarFormingMainSequence(char(fileName),label,comment,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputTimes_,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputAnalysisWeightPropertyOperator_,starFormationRateDisks_,starFormationRateSpheroids_, starFormationRateNuclearStarClusters_)
    else
       !![
       <inputParameter>
         <name>label</name>
         <source>parameters</source>
         <variable>label</variable>
         <description>A label for the star forming main sequence function.</description>
       </inputParameter>
       <inputParameter>
         <name>comment</name>
         <source>parameters</source>
         <variable>comment</variable>
         <description>A descriptive comment for the star forming main sequence function.</description>
       </inputParameter>
       <inputParameter>
         <name>massMinimum</name>
         <source>parameters</source>
         <description>Minimum stellar mass for the star forming main sequence function.</description>
       </inputParameter>
       <inputParameter>
         <name>massMaximum</name>
         <source>parameters</source>
         <description>Maximum stellar mass for the star forming main sequence function.</description>
       </inputParameter>
       <inputParameter>
         <name>countMassesPerDecade</name>
         <source>parameters</source>
         <description>Number of masses per decade at which to compute the star forming main sequence function.</description>
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
       if (parameters%isPresent('meanValueTarget')) then
          if (parameters%isPresent('meanCovarianceTarget')) then
             !![
	     <inputParameter>
               <name>meanValueTarget</name>
	       <source>parameters</source>
	       <description>The target function for likelihood calculations.</description>
             </inputParameter>
             <inputParameter>
               <name>meanCovarianceTarget</name>
	       <source>parameters</source>
	       <variable>meanCovarianceTarget1D</variable>
	       <description>The target function covariance for likelihood calculations.</description>
             </inputParameter>
             !!]
             if (size(meanCovarianceTarget1D) == size(meanValueTarget)**2) then
                allocate(meanCovarianceTarget(size(meanValueTarget),size(meanValueTarget)))
                meanCovarianceTarget=reshape(meanCovarianceTarget1D,shape(meanCovarianceTarget))
             else
                call Error_Report('functionCovariance has wrong size'//{introspection:location})
             end if
          else
             call Error_Report('functionCovariance must be specified if functionTarget is present'//{introspection:location})
          end if
       else
          if (parameters%isPresent('functionCovariance')) call Error_Report('functionTarget must be specified if functionCovariance is present'//{introspection:location})
       end if
       ! Build grid of stellar masses.
       countMasses=int(log10(massMaximum/massMinimum)*countMassesPerDecade,kind=c_size_t)+1_c_size_t
       if (countMasses <= 1_c_size_t) call Error_Report('>1 mass bins is required'//{introspection:location})
       allocate(massesStellar(countMasses))
       massesStellar=Make_Range(massMinimum,massMaximum,int(countMasses),rangeType=rangeTypeLogarithmic)
       !![
       <conditionalCall>
        <call>
         self=outputAnalysisStarFormingMainSequence(                                       &amp;
          &amp;                                     label                                , &amp;
          &amp;                                     comment                              , &amp;
          &amp;                                     massesStellar                        , &amp;
  	  &amp;                                     galacticFilter_                      , &amp;
	  &amp;                                     surveyGeometry_                      , &amp;
          &amp;                                     cosmologyFunctions_                  , &amp;
	  &amp;                                     cosmologyFunctionsData               , &amp;
          &amp;                                     outputTimes_                         , &amp;
	  &amp;                                     outputAnalysisPropertyOperator_      , &amp;
	  &amp;                                     outputAnalysisDistributionOperator_  , &amp;
	  &amp;                                     outputAnalysisWeightPropertyOperator_, &amp;
	  &amp;                                     starFormationRateDisks_              , &amp;
	  &amp;                                     starFormationRateSpheroids_          , &amp;
          &amp;                                     starFormationRateNuclearStarClusters_  &amp;
          &amp;                                     {conditions}                           &amp;
          &amp;                                    )
        </call>
        <argument name="targetLabel"          value="targetLabel"          parameterPresent="parameters"/>
        <argument name="meanValueTarget"      value="meanValueTarget"      parameterPresent="parameters"/>
        <argument name="meanCovarianceTarget" value="meanCovarianceTarget" parameterPresent="parameters"/>
       </conditionalCall>
       !!]
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"                      />
    <objectDestructor name="surveyGeometry_"                      />
    <objectDestructor name="cosmologyFunctions_"                  />
    <objectDestructor name="cosmologyFunctionsData"               />
    <objectDestructor name="outputTimes_"                         />
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    <objectDestructor name="outputAnalysisPropertyOperator_"      />
    <objectDestructor name="outputAnalysisDistributionOperator_"  />
    <objectDestructor name="outputAnalysisWeightPropertyOperator_"/>
    !!]
    return
  end function starFormingMainSequenceConstructorParameters

  function starFormingMainSequenceConstructorFile(fileName,label,comment,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputTimes_,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputAnalysisWeightPropertyOperator_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisStarFormingMainSequence} output analysis class which reads all required properties from file.
    !!}
    use :: IO_HDF5    , only : hdf5Object
    use :: HDF5_Access, only : hdf5Access
    implicit none
    type            (outputAnalysisStarFormingMainSequence    )                                :: self
    character       (len=*                                    ), intent(in   )                 :: fileName
    type            (varying_string                           ), intent(in   )                 :: label                                , comment
    class           (outputTimesClass                         ), intent(inout)                 :: outputTimes_
    class           (galacticFilterClass                      ), intent(inout), target         :: galacticFilter_
    class           (surveyGeometryClass                      ), intent(in   ), target         :: surveyGeometry_
    class           (cosmologyFunctionsClass                  ), intent(in   ), target         :: cosmologyFunctions_                  , cosmologyFunctionsData
    class           (starFormationRateDisksClass              ), intent(in   ), target         :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass          ), intent(in   ), target         :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass), intent(in   ), target         :: starFormationRateNuclearStarClusters_
    class           (outputAnalysisPropertyOperatorClass      ), intent(inout), target         :: outputAnalysisPropertyOperator_
    class           (outputAnalysisPropertyOperatorClass      ), intent(inout), target         :: outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass  ), intent(inout), target         :: outputAnalysisDistributionOperator_
    double precision                                           , allocatable  , dimension(:  ) :: meanValueTarget                      , massesStellar
    double precision                                           , allocatable  , dimension(:,:) :: meanCovarianceTarget
    double precision                                                                           :: massesStellarBinWidthLogarithmic
    type            (varying_string                           )                                :: targetLabel
    type            (hdf5Object                               )                                :: dataFile
    integer                                                                                    :: i                                    , j
    
    !$ call hdf5Access%set  ()
    call        dataFile%openFile     (fileName                                   ,readOnly=.true.                          )
    call        dataFile%readDataset  ('massStellar'                              ,         massesStellar                   )
    if (size(massesStellar) == 1)                                                                                             &
         & call dataFile%readAttribute('binWidth'                                 ,         massesStellarBinWidthLogarithmic)
    call        dataFile%readDataset  ('starFormingMainSequenceFunction'          ,         meanValueTarget                 )
    call        dataFile%readDataset  ('starFormingMainSequenceFunctionCovariance',         meanCovarianceTarget            )
    call        dataFile%readAttribute('labelTarget'                              ,         targetLabel                     )
    call        dataFile%close        (                                                                                     )
    !$ call hdf5Access%unset()
    ! Convert to logarithmic specific star formation rate.
    do i=1,size(meanCovarianceTarget,dim=1)
       do j=1,size(meanCovarianceTarget,dim=2)
          meanCovarianceTarget(i,j)=+meanCovarianceTarget(i,j) &
               &                    /meanValueTarget     (i  ) &
               &                    /meanValueTarget     (  j) &
               &                    /log(10.0d0)**2
       end do
    end do
    meanValueTarget=log10(meanValueTarget)
    ! Build the object.
    !![
    <conditionalCall>
      <call>self=starFormingMainSequenceConstructorInternal(label,comment,massesStellar,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputTimes_,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputAnalysisWeightPropertyOperator_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_,targetLabel,meanValueTarget,meanCovarianceTarget{conditions})</call>
      <argument name="massesStellarBinWidthLogarithmic" value="massesStellarBinWidthLogarithmic" condition="size(massesStellar) == 1"/>
    </conditionalCall>
    <constructorAssign variables="fileName"/>
    !!]
    return
  end function starFormingMainSequenceConstructorFile

  function starFormingMainSequenceConstructorInternal(label,comment,massesStellar,galacticFilter_,surveyGeometry_,cosmologyFunctions_,cosmologyFunctionsData,outputTimes_,outputAnalysisPropertyOperator_,outputAnalysisDistributionOperator_,outputAnalysisWeightPropertyOperator_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_,targetLabel,meanValueTarget,meanCovarianceTarget,massesStellarBinWidthLogarithmic) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisStarFormingMainSequence} output analysis class.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass
    use :: Galactic_Filters                      , only : galacticFilterClass
    use :: Error                                 , only : Error_Report
    use :: ISO_Varying_String                    , only : var_str
    use :: Node_Property_Extractors              , only : nodePropertyExtractorMassStellar           , nodePropertyExtractorRatio         , nodePropertyExtractorStarFormationRate
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorClass
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10    , outputAnalysisPropertyOperatorClass, outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc, outputAnalysisPropertyOperatorLog10, &
          &                                               outputAnalysisPropertyOperatorSequence     , propertyOperatorList
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Analysis_Utilities             , only : Output_Analysis_Output_Weight_Survey_Volume
    use :: Output_Times                          , only : outputTimesClass
    implicit none
    type            (outputAnalysisStarFormingMainSequence          )                                               :: self
    type            (varying_string                                 )                               , intent(in   ) :: label                                                       , comment
    class           (galacticFilterClass                            ), intent(inout), target                        :: galacticFilter_
    class           (surveyGeometryClass                            ), intent(in   ), target                        :: surveyGeometry_
    class           (cosmologyFunctionsClass                        ), intent(in   ), target                        :: cosmologyFunctions_                                         , cosmologyFunctionsData
    class           (outputTimesClass                               ), target                       , intent(inout) :: outputTimes_
    class           (outputAnalysisPropertyOperatorClass            ), intent(inout), target                        :: outputAnalysisPropertyOperator_                             , outputAnalysisWeightPropertyOperator_
    class           (outputAnalysisDistributionOperatorClass        ), intent(inout), target                        :: outputAnalysisDistributionOperator_
    class           (starFormationRateDisksClass                    ), intent(in   ), target                        :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                ), intent(in   ), target                        :: starFormationRateSpheroids_
    class           (starFormationRateNuclearStarClustersClass      ), intent(in   ), target                        :: starFormationRateNuclearStarClusters_
    type            (varying_string                                 ), optional                     , intent(in   ) :: targetLabel
    double precision                                                 , optional     , dimension(:  ), intent(in   ) :: meanValueTarget                                             , massesStellar
    double precision                                                 , optional     , dimension(:,:), intent(in   ) :: meanCovarianceTarget
    double precision                                                 , optional                     , intent(in   ) :: massesStellarBinWidthLogarithmic
    type            (nodePropertyExtractorMassStellar               )               , pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorStarFormationRate         )               , pointer                       :: nodePropertyExtractorStarFormationRate_
    type            (nodePropertyExtractorRatio                     )               , pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisPropertyOperatorLog10            )               , pointer                       :: outputAnalysisPropertyOperatorLog10_
    type            (outputAnalysisPropertyOperatorAntiLog10        )               , pointer                       :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc)               , pointer                       :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_
    type            (outputAnalysisPropertyOperatorSequence         )               , pointer                       :: outputAnalysisPropertyOperatorSequence_                     , outputAnalysisWeightPropertyOperatorSequence_
    type            (outputAnalysisWeightOperatorIdentity           )               , pointer                       :: outputAnalysisWeightOperator_
    type            (propertyOperatorList                           )               , pointer                       :: propertyOperatorSequence                                    , weightPropertyOperatorSequence
    double precision                                                 , allocatable, dimension(:,:)                  :: outputWeight
    integer                                                          , parameter                                    :: covarianceBinomialBinsPerDecade                 =  2
    double precision                                                 , parameter                                    :: covarianceBinomialMassHaloMinimum               = +3.000d+11, covarianceBinomialMassHaloMaximum            =1.0d15
    double precision                                                 , parameter                                    :: bufferWidthLogarithmic                          =3.0d0
    integer         (c_size_t                                       ), parameter                                    :: bufferCountMinimum                              =5
    integer         (c_size_t                                       )                                               :: iBin                                                        , bufferCount                                         , &
         &                                                                                                             countMasses
    !![
    <constructorAssign variables="*surveyGeometry_, *cosmologyFunctions_, *cosmologyFunctionsData, *starFormationRateDisks_, *starFormationRateSpheroids_, *starFormationRateNuclearStarClusters_"/>
    !!]

    ! Set properties needed for descriptor.
    countMasses              =size(massesStellar             )
    self%massMinimum         =     massesStellar(          1)
    self%massMaximum         =     massesStellar(countMasses)
    self%countMassesPerDecade=dble(countMasses-1_c_size_t)*log10(self%massMaximum/self%massMinimum)
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(countMasses,outputTimes_%count()))
    outputWeight=0.0d0
    do iBin=1,countMasses
       outputWeight(iBin,:)=Output_Analysis_Output_Weight_Survey_Volume(self%surveyGeometry_,self%cosmologyFunctions_,outputTimes_,massesStellar(iBin))
    end do
    ! Build a property operator that converts stellar mass to logarithmic form and applies a correction for cosmological luminosity distance.
    allocate(outputAnalysisPropertyOperatorLog10_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorLog10_"             constructor="outputAnalysisPropertyOperatorLog10            (                                                       )"/>
    !!]
    allocate(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
    !![
    <referenceConstruct object="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_" constructor="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc(cosmologyFunctions_,cosmologyFunctionsData,outputTimes_)"/>
    !!]
    select type (outputAnalysisPropertyOperator_)
    type is (outputAnalysisPropertyOperatorSequence)
       ! Existing property operator is a sequence operator - simply prepend our log10 and cosmological luminosity distance operators to it.
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorLog10_            )
       call outputAnalysisPropertyOperator_%prepend(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_)
       outputAnalysisPropertyOperatorSequence_ => outputAnalysisPropertyOperator_
    class default
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
    select type (outputAnalysisWeightPropertyOperator_)
    type is (outputAnalysisPropertyOperatorSequence)
       ! Existing weight property operator is a sequence operator - simply prepend our log10 operator to it.
       call outputAnalysisWeightPropertyOperator_%prepend(outputAnalysisPropertyOperatorLog10_)
       outputAnalysisWeightPropertyOperatorSequence_ => outputAnalysisWeightPropertyOperator_
    class default
       allocate(weightPropertyOperatorSequence     )
       allocate(weightPropertyOperatorSequence%next)
       weightPropertyOperatorSequence     %operator_ => outputAnalysisPropertyOperatorLog10_
       weightPropertyOperatorSequence%next%operator_ => outputAnalysisWeightPropertyOperator_
       allocate(outputAnalysisWeightPropertyOperatorSequence_)
       !![
       <referenceConstruct object="outputAnalysisWeightPropertyOperatorSequence_" constructor="outputAnalysisPropertyOperatorSequence(weightPropertyOperatorSequence)"/>
       !!]
    end select
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_)
    !![
    <referenceConstruct object="outputAnalysisPropertyUnoperator_">
     <constructor>
      outputAnalysisPropertyOperatorAntiLog10()
     </constructor>
    </referenceConstruct>
    !!]
    ! Create a stellar mass property extractor.
    allocate(nodePropertyExtractor_)
    !![
    <referenceConstruct object="nodePropertyExtractor_">
      <constructor>
	nodePropertyExtractorMassStellar()
      </constructor>
    </referenceConstruct>
    !!]
    ! Create a specific star formation rate property extractor.
    allocate(nodePropertyExtractorStarFormationRate_)
    !![
    <referenceConstruct object="nodePropertyExtractorStarFormationRate_">
      <constructor>
	nodePropertyExtractorStarFormationRate(                                                                             &amp;
	&amp;                                  starFormationRateDisks_              =starFormationRateDisks_              , &amp;
	&amp;                                  starFormationRateSpheroids_          =starFormationRateSpheroids_          , &amp;
        &amp;                                  starFormationRateNuclearStarClusters_=starFormationRateNuclearStarClusters_  &amp; 
        &amp;                                 )
      </constructor>
    </referenceConstruct>
    !!]
    allocate(outputAnalysisWeightPropertyExtractor_)
    !![
    <referenceConstruct object="outputAnalysisWeightPropertyExtractor_">
      <constructor>
	nodePropertyExtractorRatio(                                                              &amp;
	&amp;                      name                ="specificStarFormationRate"            , &amp;
	&amp;                      description         ="Specific star formation rate."        , &amp;
	&amp;                      propertyNumerator_  =nodePropertyExtractorStarFormationRate_, &amp;
	&amp;                      propertyDenominator_=nodePropertyExtractor_                   &amp;
        &amp;                     )
      </constructor>
    </referenceConstruct>
    !!]
    ! Create an identity weight operator.
    allocate(outputAnalysisWeightOperator_)
    !![
    <referenceConstruct object="outputAnalysisWeightOperator_" constructor="outputAnalysisWeightOperatorIdentity()"/>
    !!]
    ! Compute the number of buffer bins to add to either side of the mass function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    if (size(massesStellar) == 1) then
       if (.not.present(massesStellarBinWidthLogarithmic)) call Error_Report('bin width must be specified'//{introspection:location})
       bufferCount=max(int(bufferWidthLogarithmic/massesStellarBinWidthLogarithmic        )+1,bufferCountMinimum)
    else
       bufferCount=max(int(bufferWidthLogarithmic/log10(massesStellar(2)/massesStellar(1)))+1,bufferCountMinimum)
    end if
    ! Construct the object.
    self%outputAnalysisMeanFunction1D=                                                                                      &
         & outputAnalysisMeanFunction1D(                                                                                    &
         &                              var_str('starFormingMainSequence')//label                                         , &
         &                              comment                                                                           , &
         &                              var_str('massStellar'                                                            ), &
         &                              var_str('Stellar mass at the bin center'                                         ), &
         &                              var_str('$\mathrm{M}_\odot$'                                                     ), &
         &                              massSolar                                                                         , &
         &                              var_str('rateStarFormationSpecific'                                              ), &
         &                              var_str('Logarithmic specific star formation rate averaged over each bin'        ), &
         &                              var_str(' '                                                                      ), &
         &                              1.0d0                                                                             , &
         &                              log10(massesStellar)                                                              , &
         &                              bufferCount                                                                       , &
         &                              outputWeight                                                                      , &
         &                              nodePropertyExtractor_                                                            , &
         &                              outputAnalysisWeightPropertyExtractor_                                            , &
         &                              outputAnalysisPropertyOperatorSequence_                                           , &
         &                              outputAnalysisWeightPropertyOperatorSequence_                                     , &
         &                              outputAnalysisPropertyUnoperator_                                                 , &
         &                              outputAnalysisWeightOperator_                                                     , &
         &                              outputAnalysisDistributionOperator_                                               , &
         &                              galacticFilter_                                                                   , &
         &                              outputTimes_                                                                      , &
         &                              outputAnalysisCovarianceModelBinomial                                             , &
         &                              covarianceBinomialBinsPerDecade                                                   , &
         &                              covarianceBinomialMassHaloMinimum                                                 , &
         &                              covarianceBinomialMassHaloMaximum                                                 , &
         &                              .false.                                                                           , &
         &                              var_str('$M_\star\, [\mathrm{M}_\odot]$'                                         ), &
         &                              var_str('$\langle \log_{10} (\dot{M}_\star/M_\star / \mathrm{Gyr}^{-1}) \rangle$'), &
         &                              .true.                                                                            , &
         &                              .false.                                                                           , &
         &                              targetLabel                                                                       , &
         &                              meanValueTarget                                                                   , &
         &                              meanCovarianceTarget                                                              , &
         &                              massesStellarBinWidthLogarithmic                                                    &
         &                             )
    !![
    <objectDestructor name="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc_"/>
    <objectDestructor name="outputAnalysisPropertyOperatorSequence_"         />
    <objectDestructor name="outputAnalysisWeightPropertyOperatorSequence_"   />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"               />
    <objectDestructor name="nodePropertyExtractor_"                          />
    <objectDestructor name="nodePropertyExtractorStarFormationRate_"         />
    <objectDestructor name="outputAnalysisWeightPropertyExtractor_"          />   
    !!]
    nullify(propertyOperatorSequence)
    return
  end function starFormingMainSequenceConstructorInternal
