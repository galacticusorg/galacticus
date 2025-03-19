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
Implements a generic 1D volume function (i.e. number density of objects binned by some property, e.g. a
mass function) output analysis class.
!!}

  use               :: Galactic_Filters                        , only : galacticFilterClass
  use   , intrinsic :: ISO_C_Binding                           , only : c_size_t
  use               :: ISO_Varying_String                      , only : varying_string
  use               :: Node_Property_Extractors                , only : nodePropertyExtractorClass
  !$ use            :: OMP_Lib                                 , only : omp_lock_kind
  use               :: Output_Analysis_Distribution_Normalizers, only : outputAnalysisDistributionNormalizerClass
  use               :: Output_Analysis_Distribution_Operators  , only : outputAnalysisDistributionOperatorClass
  use               :: Output_Analysis_Property_Operators      , only : outputAnalysisPropertyOperatorClass
  use               :: Output_Analysis_Weight_Operators        , only : outputAnalysisWeightOperatorClass
  use               :: Output_Times                            , only : outputTimesClass
  use               :: Output_Analyses_Options                 , only : enumerationOutputAnalysisCovarianceModelType

  !![
  <outputAnalysis name="outputAnalysisCrossCorrelator1D">
   <description>
     A generic 1D cross-correlator (i.e. the cross-correlation of two weights binned by some property, e.g. a mass function)
     output analysis class.
  
     The assumptions used when constructing the covariance matrix are controlled by the parameter {\normalfont \ttfamily
     [covarianceModel]}, and follow the method described for the \refClass{outputAnalysisVolumeFunction1D} class.
   </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisCrossCorrelator1D
     !!{
     A generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output
     analysis class.
     !!}
     private
     class           (nodePropertyExtractorClass                  ), pointer                     :: nodePropertyExtractor_                => null()
     class           (outputAnalysisPropertyOperatorClass         ), pointer                     :: outputAnalysisPropertyOperator_       => null()                                                             , &
          &                                                                                         outputAnalysisPropertyUnoperator_     => null()
     class           (outputAnalysisWeightOperatorClass           ), pointer                     :: outputAnalysisWeightOperator1_        => null(), outputAnalysisWeightOperator2_                    => null()
     class           (outputAnalysisDistributionOperatorClass     ), pointer                     :: outputAnalysisDistributionOperator_   => null()
     class           (outputAnalysisDistributionNormalizerClass   ), pointer                     :: outputAnalysisDistributionNormalizer_ => null()
     class           (galacticFilterClass                         ), pointer                     :: galacticFilter_                       => null()
     class           (outputTimesClass                            ), pointer                     :: outputTimes_                          => null()
     double precision                                              , dimension(:,:), allocatable :: outputWeight                                   , functionCovariance                                         , &
          &                                                                                         weightMainBranch
     double precision                                              , dimension(:  ), allocatable :: binCenter                                      , weightMainBranchCross
     double precision                                              , dimension(:  ), allocatable :: binMinimum                                     , binMaximum
     integer         (c_size_t                                    )                              :: binCount                                       , bufferCount                                                , &
          &                                                                                         binCountTotal                                  , covarianceModelBinomialBinCount
     integer                                                                                     :: covarianceBinomialBinsPerDecade
     type            (enumerationOutputAnalysisCovarianceModelType)                              :: covarianceModel
     double precision                                                                            :: covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum                          , &
          &                                                                                         covarianceModelHaloMassMinimumLogarithmic      , covarianceModelHaloMassIntervalLogarithmicInverse          , &
          &                                                                                         binWidth
     logical                                                                                     :: finalized
     !$ integer      (omp_lock_kind                               )                              :: accumulateLock
   contains
     !![
     <methods>
       <method description="Return the results of the volume function operator." method="results"         />
       <method description="Finalize the analysis of this function."             method="finalizeAnalysis"/>
     </methods>
     !!]
     final     ::                     crossCorrelator1DDestructor
     procedure :: analyze          => crossCorrelator1DAnalyze
     procedure :: finalize         => crossCorrelator1DFinalize
     procedure :: results          => crossCorrelator1DResults
     procedure :: reduce           => crossCorrelator1DReduce
     procedure :: finalizeAnalysis => crossCorrelator1DFinalizeAnalysis
  end type outputAnalysisCrossCorrelator1D

  interface outputAnalysisCrossCorrelator1D
     !!{
     Constructors for the ``crossCorrelator1D'' output analysis class.
     !!}
     module procedure crossCorrelator1DConstructorParameters
     module procedure crossCorrelator1DConstructorInternal
  end interface outputAnalysisCrossCorrelator1D

contains

  function crossCorrelator1DConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``crossCorrelator1D'' output analysis class which takes a parameter set as input.
    !!}
    use :: Error                  , only : Error_Report
    use :: Input_Parameters       , only : inputParameter                                , inputParameters
    use :: Output_Analyses_Options, only : enumerationOutputAnalysisCovarianceModelEncode
    implicit none
    type            (outputAnalysisCrossCorrelator1D          )                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    class           (nodePropertyExtractorClass               ), pointer                     :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass      ), pointer                     :: outputAnalysisPropertyOperator_          , outputAnalysisPropertyUnoperator_
    class           (outputAnalysisWeightOperatorClass        ), pointer                     :: outputAnalysisWeightOperator1_           , outputAnalysisWeightOperator2_
    class           (outputAnalysisDistributionOperatorClass  ), pointer                     :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass), pointer                     :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                      ), pointer                     :: galacticFilter_
    class           (outputTimesClass                         ), pointer                     :: outputTimes_
    double precision                                           , dimension(:  ), allocatable :: binCenter                                , outputWeight
    integer         (c_size_t                                 )                              :: bufferCount
    type            (varying_string                           )                              :: covarianceModel
    integer                                                                                  :: covarianceBinomialBinsPerDecade
    type            (inputParameters                          )                              :: unoperatorParameters
    double precision                                                                         :: covarianceBinomialMassHaloMinimum        , covarianceBinomialMassHaloMaximum , &
         &                                                                                      binWidth

         ! Check and read parameters.
    unoperatorParameters=parameters%subParameters('unoperatorParameters',requireValue=.false.)
    !![
    <objectBuilder class="nodePropertyExtractor"                name="nodePropertyExtractor_"                source="parameters"                                                       />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"       source="parameters"                                                       />
    <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyUnoperator_"     source="unoperatorParameters"                                             />
    <objectBuilder class="outputAnalysisWeightOperator"         name="outputAnalysisWeightOperator1_"        source="parameters"          parameterName="outputAnalysisWeightOperator1"/>
    <objectBuilder class="outputAnalysisWeightOperator"         name="outputAnalysisWeightOperator2_"        source="parameters"          parameterName="outputAnalysisWeightOperator2"/>
    <objectBuilder class="outputAnalysisDistributionOperator"   name="outputAnalysisDistributionOperator_"   source="parameters"                                                       />
    <objectBuilder class="outputAnalysisDistributionNormalizer" name="outputAnalysisDistributionNormalizer_" source="parameters"                                                       />
    <objectBuilder class="galacticFilter"                       name="galacticFilter_"                       source="parameters"                                                       />
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"                                                       />
    !!]
    allocate(binCenter   (int(parameters%count('binCenter'))                          ))
    allocate(outputWeight(int(parameters%count('binCenter'))*self%outputTimes_%count()))
    if (parameters%count('outputWeight') /= parameters%count('binCenter')*self%outputTimes_%count()) &
         & call Error_Report('incorrect number of output weights provided'//{introspection:location})
    !![
    <inputParameter>
      <name>binCenter</name>
      <source>parameters</source>
      <variable>binCenter</variable>
      <description>The value of the property at the center of each bin.</description>
    </inputParameter>
    !!]
    if (size(binCenter) == 1) then
       !![
       <inputParameter>
	 <name>binWidth</name>
	 <source>parameters</source>
	 <variable>binWidth</variable>
	 <description>The width of the bins.</description>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>bufferCount</name>
      <source>parameters</source>
      <variable>bufferCount</variable>
      <description>The number of buffer bins to include below and above the range of actual bins.</description>
    </inputParameter>
    <inputParameter>
      <name>outputWeight</name>
      <source>parameters</source>
      <variable>outputWeight</variable>
      <description>The weight to assign to each bin at each output.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceModel</name>
      <source>parameters</source>
      <variable>covarianceModel</variable>
      <description>The model to use for computing covariances.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing volume function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <conditionalCall>
      <call>
	self=outputAnalysisCrossCorrelator1D(                                                                                                    &amp;
        &amp;                                binCenter                                                                                         , &amp;
        &amp;                                bufferCount                                                                                       , &amp;
        &amp;                                reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),self%outputTimes_%count()]), &amp;
        &amp;                                nodePropertyExtractor_                                                                            , &amp;
        &amp;                                outputAnalysisPropertyOperator_                                                                   , &amp;
        &amp;                                outputAnalysisPropertyUnoperator_                                                                 , &amp;
        &amp;                                outputAnalysisWeightOperator1_                                                                    , &amp;
        &amp;                                outputAnalysisWeightOperator2_                                                                    , &amp;
        &amp;                                outputAnalysisDistributionOperator_                                                               , &amp;
        &amp;                                outputAnalysisDistributionNormalizer_                                                             , &amp;
        &amp;                                galacticFilter_                                                                                   , &amp;
        &amp;                                outputTimes_                                                                                      , &amp;
        &amp;                                enumerationOutputAnalysisCovarianceModelEncode(char(covarianceModel),includesPrefix=.false.)      , &amp;
        &amp;                                covarianceBinomialBinsPerDecade                                                                   , &amp;
        &amp;                                covarianceBinomialMassHaloMinimum                                                                 , &amp;
        &amp;                                covarianceBinomialMassHaloMaximum                                                                   &amp;
        &amp;                                {conditions}                                                                                        &amp;
        &amp;                               )
      </call>
      <argument name="binWidth" value="binWidth" condition="size(binCenter) == 1"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"               />
    <objectDestructor name="outputAnalysisPropertyOperator_"      />
    <objectDestructor name="outputAnalysisPropertyUnoperator_"    />
    <objectDestructor name="outputAnalysisWeightOperator1_"       />
    <objectDestructor name="outputAnalysisWeightOperator2_"       />
    <objectDestructor name="outputAnalysisDistributionOperator_"  />
    <objectDestructor name="outputAnalysisDistributionNormalizer_"/>
    <objectDestructor name="galacticFilter_"                      />
    <objectDestructor name="outputTimes_"                         />
    !!]
    return
  end function crossCorrelator1DConstructorParameters

  function crossCorrelator1DConstructorInternal(binCenter,bufferCount,outputWeight,nodePropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisPropertyUnoperator_,outputAnalysisWeightOperator1_,outputAnalysisWeightOperator2_,outputAnalysisDistributionOperator_,outputAnalysisDistributionNormalizer_,galacticFilter_,outputTimes_,covarianceModel,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,binWidth) result (self)
    !!{
    Constructor for the ``crossCorrelator1D'' output analysis class for internal use.
    !!}
    use    :: Error                   , only : Error_Report
    use    :: Node_Property_Extractors, only : nodePropertyExtractorClass           , nodePropertyExtractorScalar
    use    :: Output_Analyses_Options , only : outputAnalysisCovarianceModelBinomial
    !$ use :: OMP_Lib                 , only : OMP_Init_Lock
    implicit none
    type            (outputAnalysisCrossCorrelator1D             )                                          :: self
    double precision                                              , intent(in   )          , dimension(:  ) :: binCenter
    integer         (c_size_t                                    ), intent(in   )                           :: bufferCount
    double precision                                              , intent(in   )          , dimension(:,:) :: outputWeight
    class           (nodePropertyExtractorClass                  ), intent(in   ), target                   :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass         ), intent(in   ), target                   :: outputAnalysisPropertyOperator_      , outputAnalysisPropertyUnoperator_
    class           (outputAnalysisWeightOperatorClass           ), intent(in   ), target                   :: outputAnalysisWeightOperator1_       , outputAnalysisWeightOperator2_
    class           (outputAnalysisDistributionOperatorClass     ), intent(in   ), target                   :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass   ), intent(in   ), target                   :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                         ), intent(in   ), target                   :: galacticFilter_
    class           (outputTimesClass                            ), intent(in   ), target                   :: outputTimes_
    type            (enumerationOutputAnalysisCovarianceModelType), intent(in   )                           :: covarianceModel
    integer                                                       , intent(in   ), optional                 :: covarianceBinomialBinsPerDecade
    double precision                                              , intent(in   ), optional                 :: covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum, &
         &                                                                                                     binWidth
    integer         (c_size_t                                    )                                          :: i
    !![
    <constructorAssign variables="binCenter, bufferCount, outputWeight, *nodePropertyExtractor_, *outputAnalysisPropertyOperator_, *outputAnalysisPropertyUnoperator_, *outputAnalysisWeightOperator1_, *outputAnalysisWeightOperator2_, *outputAnalysisDistributionOperator_, *outputAnalysisDistributionNormalizer_, *galacticFilter_, *outputTimes_, covarianceModel, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum, binWidth"/>
    !!]

    ! Validate.
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       ! This is acceptable.
    class default
       call Error_Report('property extrator must be of scalar class'//{introspection:location})
    end select
    ! Count bins.
    self%binCount     =size(binCenter,kind=c_size_t)
    self%binCountTotal=self%binCount+2*bufferCount
    ! Determine bin minima and maxima. Allocate arrays for bins that include buffer regions.
    allocate(self%binMinimum(-bufferCount+1:-bufferCount+self%binCountTotal))
    allocate(self%binMaximum(-bufferCount+1:-bufferCount+self%binCountTotal))
    if (present(binWidth)) then
       self%binMinimum(1:self%binCount)=+binCenter-0.5d0*binWidth
       self%binMaximum(1:self%binCount)=+binCenter+0.5d0*binWidth
       do i=1,bufferCount
          self%binMinimum(            1-i)=self%binMinimum(            1)-dble(i)*binWidth
          self%binMaximum(            1-i)=self%binMaximum(            1)-dble(i)*binWidth
          self%binMinimum(self%binCount+i)=self%binMinimum(self%binCount)+dble(i)*binWidth
          self%binMaximum(self%binCount+i)=self%binMaximum(self%binCount)+dble(i)*binWidth
       end do
    else
       do i=1,self%binCount
          if (i == 1) then
             self%binMinimum(i)=+binCenter(i)+0.5d0*(self%binCenter(i  )-self%binCenter(i+1))
          else
             self%binMinimum(i)=             +0.5d0*(self%binCenter(i-1)+self%binCenter(i  ))
          end if
          if (i == self%binCount) then
             self%binMaximum(i)=+binCenter(i)+0.5d0*(self%binCenter(i  )-self%binCenter(i-1))
          else
             self%binMaximum(i)=             +0.5d0*(self%binCenter(i  )+self%binCenter(i+1))
          end if
       end do
       do i=1,bufferCount
          self%binMinimum(            1-i)=self%binMinimum(            1)-dble(i)*(self%binMinimum(            2)-self%binMinimum(              1))
          self%binMaximum(            1-i)=self%binMaximum(            1)-dble(i)*(self%binMinimum(            2)-self%binMinimum(              1))
          self%binMinimum(self%binCount+i)=self%binMinimum(self%binCount)+dble(i)*(self%binMaximum(self%binCount)-self%binMaximum(self%binCount-1))
          self%binMaximum(self%binCount+i)=self%binMaximum(self%binCount)+dble(i)*(self%binMaximum(self%binCount)-self%binMaximum(self%binCount-1))
       end do
    end if
    ! Allocate and initialize function covariance.
    allocate(self%functionCovariance(self%binCount,self%binCount))
    self%functionCovariance=0.0d0
    ! Allocate and initialize binomial covariance model halo arrays if necessary.
    if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       self%covarianceModelBinomialBinCount                  =int(log10(self%covarianceBinomialMassHaloMaximum/self%covarianceBinomialMassHaloMinimum)*dble(self%covarianceBinomialBinsPerDecade)+0.5d0)
       self%covarianceModelHaloMassMinimumLogarithmic        =log10(self%covarianceBinomialMassHaloMinimum)
       self%covarianceModelHaloMassIntervalLogarithmicInverse=dble(self%covarianceModelBinomialBinCount)/log10(self%covarianceBinomialMassHaloMaximum/self%covarianceBinomialMassHaloMinimum)
       allocate(self%weightMainBranch     (self%binCount,self%covarianceModelBinomialBinCount))
       allocate(self%weightMainBranchCross(self%covarianceModelBinomialBinCount))
       self%weightMainBranch     =0.0d0
       self%weightMainBranchCross=0.0d0
    end if
    ! Initialize finalization status.
    self%finalized=.false.
    ! Initialize OpenMP accumulation lock.
    !$ call OMP_Init_Lock(self%accumulateLock)
   return
  end function crossCorrelator1DConstructorInternal

  subroutine crossCorrelator1DDestructor(self)
    !!{
    Destructor for  the ``crossCorrelator1D'' output analysis class.
    !!}
    !$ use :: OMP_Lib, only : OMP_Destroy_Lock
    implicit none
    type(outputAnalysisCrossCorrelator1D), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"               />
    <objectDestructor name="self%outputAnalysisPropertyOperator_"      />
    <objectDestructor name="self%outputAnalysisPropertyUnoperator_"    />
    <objectDestructor name="self%outputAnalysisWeightOperator1_"       />
    <objectDestructor name="self%outputAnalysisWeightOperator2_"       />
    <objectDestructor name="self%outputAnalysisDistributionOperator_"  />
    <objectDestructor name="self%outputAnalysisDistributionNormalizer_"/>
    <objectDestructor name="self%galacticFilter_"                      />
    <objectDestructor name="self%outputTimes_"                         />
    !!]
    ! Destroy OpenMP lock.
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine crossCorrelator1DDestructor

  subroutine crossCorrelator1DAnalyze(self,node,iOutput)
    !!{
    Implement a crossCorrelator1D output analysis.
    !!}
    use    :: Galacticus_Nodes        , only : nodeComponentBasic                   , treeNode
    use    :: Node_Property_Extractors, only : nodePropertyExtractorScalar
    use    :: Output_Analyses_Options , only : outputAnalysisCovarianceModelBinomial, enumerationOutputAnalysisPropertyTypeType, enumerationOutputAnalysisPropertyQuantityType
    !$ use :: OMP_Lib                 , only : OMP_Set_Lock                         , OMP_Unset_Lock
    implicit none
    class           (outputAnalysisCrossCorrelator1D              ), intent(inout)                 :: self
    type            (treeNode                                     ), intent(inout)                 :: node
    integer         (c_size_t                                     ), intent(in   )                 :: iOutput
    double precision                                               , allocatable  , dimension(:  ) :: distribution
    double precision                                               , allocatable  , dimension(:,:) :: covariance
    class           (nodeComponentBasic                           ), pointer                       :: basic
    double precision                                                                               :: propertyValue         , weightValue1, &
         &                                                                                            propertyValueIntrinsic, weightValue2
    type            (enumerationOutputAnalysisPropertyTypeType    )                                :: propertyType 
    type            (enumerationOutputAnalysisPropertyQuantityType)                                :: propertyQuantity
    integer         (c_size_t                                     )                                :: j                     , k           , &
         &                                                                                            indexHaloMass

    ! If weights for this output are all zero, we can skip analysis.
    if (all(self%outputWeight(:,iOutput) == 0.0d0)) return
    ! Filter this node.
    if (.not.self%galacticFilter_%passes(node)) return
    ! Allocate work arrays.
    allocate(distribution(-self%bufferCount+1:self%binCount+self%bufferCount))
    ! Extract the property from the node.
    propertyType             =self%nodePropertyExtractor_%type    (    )
    propertyQuantity         =self%nodePropertyExtractor_%quantity(    )
    select type (extractor_ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       propertyValue         =                 extractor_%extract (node)
    class default
       propertyValue         =0.0d0
    end select
    propertyValueIntrinsic   =propertyValue
    ! Apply property operators.
    propertyValue=self%outputAnalysisPropertyOperator_%operate(propertyValue,node,propertyType,iOutput)
    ! Apply distribution operators.
    distribution=self%outputAnalysisDistributionOperator_%operateScalar(propertyValue,propertyType,self%binMinimum,self%binMaximum,iOutput,node)
    ! Apply output weights.
    distribution(1:self%binCount)=+distribution              (1:self%binCount        ) &
         &                        *self        %outputWeight ( :             ,iOutput)
    ! Compute the weights.
    weightValue1=node%hostTree%volumeWeight
    weightValue2=node%hostTree%volumeWeight
    ! Apply weight operators.
    weightValue1=self%outputAnalysisWeightOperator1_%operate(weightValue1,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,iOutput)
    weightValue2=self%outputAnalysisWeightOperator2_%operate(weightValue2,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,iOutput)
    ! Accumulate covariance. If using the binomial model for main branch galaxies, handle them separately.
    if (node%isOnMainBranch() .and. self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       ! Find the bin to which this halo mass belongs.
       basic         => node%basic()
       indexHaloMass =  floor((log10(basic%mass())-self%covarianceModelHaloMassMinimumLogarithmic)*self%covarianceModelHaloMassIntervalLogarithmicInverse)+1
       ! Accumulate weights to halo mass arrays.
       if (indexHaloMass >= 1 .and. indexHaloMass <= self%covarianceModelBinomialBinCount) then
          self             %weightMainBranch     ( :             ,indexHaloMass)=     &
               &  +self    %weightMainBranch     ( :             ,indexHaloMass)      &
               &  +         distribution         (1:self%binCount)                    &
               &  *         weightValue1
          self             %weightMainBranchCross(                indexHaloMass)=     &
               &  +    self%weightMainBranchCross(                indexHaloMass)      &
               &  +sum(     distribution          (1:self%binCount              ))**2 &
               &  *         weightValue1                                              &
               &  *         weightValue2
       end if
    else
       ! Construct contribution to the covariance matrix assuming Poisson statistics.
       allocate(covariance(self%binCount,self%binCount))
       forall(j=1:self%binCount)
          forall(k=j:self%binCount)
             covariance(j,k)=+distribution(j  )*weightValue1 &
                  &          *distribution(  k)*weightValue2
             covariance(k,j)=+covariance  (j,k)
          end forall
       end forall
       ! Accumulate covariance.
       !$ call OMP_Set_Lock(self%accumulateLock)
       self        %functionCovariance= &
            & +self%functionCovariance  &
            & +             covariance
       !$ call OMP_Unset_Lock(self%accumulateLock)
       deallocate(covariance)
    end if
    ! Deallocate workspace.
    deallocate(distribution)
    return
  end subroutine crossCorrelator1DAnalyze

  subroutine crossCorrelator1DReduce(self,reduced)
    !!{
    Implement a crossCorrelator1D output analysis reduction.
    !!}
    use    :: Error                  , only : Error_Report
    use    :: Output_Analyses_Options, only : outputAnalysisCovarianceModelBinomial
    !$ use :: OMP_Lib                , only : OMP_Set_Lock                         , OMP_Unset_Lock
    implicit none
    class(outputAnalysisCrossCorrelator1D), intent(inout) :: self
    class(outputAnalysisClass            ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisCrossCorrelator1D)
       !$ call OMP_Set_Lock(reduced%accumulateLock)
       if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
          reduced%weightMainBranch     =reduced%weightMainBranch     +self%weightMainBranch
          reduced%weightMainBranchCross=reduced%weightMainBranchCross+self%weightMainBranchCross
       end if
       reduced%functionCovariance      =reduced%functionCovariance   +self%functionCovariance
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine crossCorrelator1DReduce

  subroutine crossCorrelator1DFinalize(self,groupName)
    !!{
    Implement a crossCorrelator1D output analysis finalization.
    !!}
    implicit none
    class(outputAnalysisCrossCorrelator1D), intent(inout)           :: self
    type (varying_string                 ), intent(in   ), optional :: groupName

    ! Finalize analysis.
    call self%finalizeAnalysis()
    return
  end subroutine crossCorrelator1DFinalize

  subroutine crossCorrelator1DFinalizeAnalysis(self)
    !!{
    Compute final covariances and normalize.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities          , only : mpiSelf
#endif
    use :: Output_Analyses_Options, only : outputAnalysisCovarianceModelBinomial
    implicit none
    class           (outputAnalysisCrossCorrelator1D), intent(inout) :: self
    integer         (c_size_t                       )                :: i                    , j, &
         &                                                              m
    double precision                                                 :: weightMainBranchTotal

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
    ! If using the binomial covariance model, add on the contribution to covariance from main branch galaxies.
    if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       do m=1,self%covarianceModelBinomialBinCount
          weightMainBranchTotal=sum(self%weightMainBranch(:,m))
          if (weightMainBranchTotal > 0.0d0) then
             do i=1,self%binCount
                ! For on-diagonal terms the variance is just n pᵢ (1-pᵢ). Here, pᵢ is the probability of contributing to bin i in
                ! our function. We find this by dividing the weight of main branch galaxies in bin i by the sum of weights over
                ! all bins. Finally, n is just the sum of squared weights for this halo mass bin.
                self               %functionCovariance   (i,i)=                       &
                     &         self%functionCovariance   (i,i)                        &
                     & +       self%weightMainBranchCross(  m)                        & ! +   n
                     & *       self%weightMainBranch     (i,m)/weightMainBranchTotal  & !     pᵢ
                     & *(1.0d0-self%weightMainBranch     (i,m)/weightMainBranchTotal)   !  (1-pᵢ)
                do j=1,self%binCount
                   if (i == j) cycle
                   ! For off-diagonal terms the covariance is - n pᵢ pⱼ.
                   self            %functionCovariance   (i,j)=                       &
                        &  +   self%functionCovariance   (i,j)                        &
                        &  -   self%weightMainBranchCross(  m)                        & ! -   n
                        &  *   self%weightMainBranch     (i,m)/weightMainBranchTotal  & !     pᵢ
                        &  *   self%weightMainBranch     (j,m)/weightMainBranchTotal    !     pⱼ
                end do
             end do
          end if
       end do
    end if
#ifdef USEMPI
    ! If running under MPI, perform a summation reduction across all processes.
    self%functionCovariance=mpiSelf%sum(self%functionCovariance)
#endif
    ! Apply final distribution operators - pass only the non-buffer bin values here.
    call self%outputAnalysisDistributionNormalizer_%normalize(                                                               &
         &                                                    covariance          =self%functionCovariance                 , &
         &                                                    propertyValueMinimum=self%binMinimum        (1:self%binCount), &
         &                                                    propertyValueMaximum=self%binMaximum        (1:self%binCount)  &
         &                                                   )
    ! Apply any "unoperator" to output bin values. This can be used to reverse transformations (e.g. if masses were converted to
    ! log10 for analysis, that can be undone here).
    do i=1,self%binCount
       self%binCenter(i)=self%outputAnalysisPropertyUnoperator_%operate(self%binCenter(i))
    end do
    return
  end subroutine crossCorrelator1DFinalizeAnalysis

  subroutine crossCorrelator1DResults(self,binCenter,functionCovariance)
    !!{
    Implement a crossCorrelator1D output analysis finalization.
    !!}
    implicit none
    class           (outputAnalysisCrossCorrelator1D)                             , intent(inout)           :: self
    double precision                                 , allocatable, dimension(:  ), intent(inout), optional :: binCenter
    double precision                                 , allocatable, dimension(:,:), intent(inout), optional :: functionCovariance

    ! Finalize analysis.
    call self%finalizeAnalysis()
    ! Return results.
    if (present(binCenter         )) then
       if (allocated(binCenter         )) deallocate(binCenter         )
       allocate(binCenter(size(self%binCenter)))
       binCenter         =self%binCenter
    end if
    if (present(functionCovariance)) then
       if (allocated(functionCovariance)) deallocate(functionCovariance)
       allocate(functionCovariance(size(self%functionCovariance,dim=1),size(self%functionCovariance,dim=2)))
       functionCovariance=self%functionCovariance
    end if
    return
  end subroutine crossCorrelator1DResults
