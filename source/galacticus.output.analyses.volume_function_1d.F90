!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a generic 1D volume function (i.e. number density of objects binned by some property, e.g. a
!% mass function) output analysis class.

  !$ use            OMP_Lib
  use, intrinsic :: ISO_C_Binding
  use               ISO_Varying_String
  use               Output_Analysis_Property_Extractions
  use               Output_Analysis_Property_Operators
  use               Output_Analysis_Weight_Operators
  use               Output_Analysis_Distribution_Operators
  use               Output_Analysis_Distribution_Normalizers
  use               Output_Analyses_Options
  use               Galactic_Filters

  !# <outputAnalysis name="outputAnalysisVolumeFunction1D" defaultThreadPrivate="yes">
  !#  <description>A generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisClass) :: outputAnalysisVolumeFunction1D
     !% A generic 1D volume function (i.e. number density of objects binned by some property, e.g. a mass function) output
     !% analysis class.
     private
     type            (varying_string                           )                              :: label                                          , comment                                          , &
          &                                                                                      propertyLabel                                  , propertyComment                                  , &
          &                                                                                      distributionLabel                              , distributionComment                              , &
          &                                                                                      propertyUnits                                  , distributionUnits
     double precision                                                                         :: propertyUnitsInSI                              , distributionUnitsInSI
     class           (outputAnalysisPropertyExtractorClass     ), pointer                     :: outputAnalysisPropertyExtractor_      => null()
     class           (outputAnalysisPropertyOperatorClass      ), pointer                     :: outputAnalysisPropertyOperator_       => null()                                                   , &
          &                                                                                      outputAnalysisPropertyUnoperator_     => null()
     class           (outputAnalysisWeightOperatorClass        ), pointer                     :: outputAnalysisWeightOperator_         => null()
     class           (outputAnalysisDistributionOperatorClass  ), pointer                     :: outputAnalysisDistributionOperator_   => null()
     class           (outputAnalysisDistributionNormalizerClass), pointer                     :: outputAnalysisDistributionNormalizer_ => null()
     class           (galacticFilterClass                      ), pointer                     :: galacticFilter_                       => null()
     double precision                                           , dimension(:,:), allocatable :: outputWeight                                   , functionCovariance                               , &
          &                                                                                      weightMainBranch                               , weightSquaredMainBranch
     double precision                                           , dimension(:  ), allocatable :: binCenter                                      , functionValue
     double precision                                           , dimension(:  ), allocatable :: binMinimum                                     , binMaximum
     integer         (c_size_t                                 )                              :: binCount                                       , bufferCount                                      , &
          &                                                                                      binCountTotal                                  , covarianceModelBinomialBinCount
     integer                                                                                  :: covarianceModel                                , covarianceBinomialBinsPerDecade
     double precision                                                                         :: covarianceBinomialMassHaloMinimum              , covarianceBinomialMassHaloMaximum                , &
          &                                                                                      covarianceModelHaloMassMinimumLogarithmic      , covarianceModelHaloMassIntervalLogarithmicInverse
     logical                                                                                  :: finalized
     !$ integer      (omp_lock_kind                            )                              :: accumulateLock
   contains
     !@ <objectMethods>
     !@   <object>outputAnalysisVolumeFunction1D</object>
     !@   <objectMethod>
     !@     <method>results</method>
     !@     <arguments>\doubleone\ [binCenter]\arginout, \doubletwo\ [functionValue]\arginout, \doubletwo\ [functionCovariance]\arginout</arguments>
     !@     <type>\void</type>
     !@     <description>Return the results of the volume function operator.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::             volumeFunction1DDestructor
     procedure :: analyze  => volumeFunction1DAnalyze
     procedure :: finalize => volumeFunction1DFinalize
     procedure :: results  => volumeFunction1DResults
  end type outputAnalysisVolumeFunction1D

  interface outputAnalysisVolumeFunction1D
     !% Constructors for the ``volumeFunction1D'' output analysis class.
     module procedure volumeFunction1DConstructorParameters
     module procedure volumeFunction1DConstructorInternal
  end interface outputAnalysisVolumeFunction1D

contains

  function volumeFunction1DConstructorParameters(parameters) result(self)
    !% Constructor for the ``volumeFunction1D'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    use Galacticus_Output_Times
    use Memory_Management
    use Galacticus_Error
    implicit none
    type            (outputAnalysisVolumeFunction1D           )                              :: self
    type            (inputParameters                          ), intent(inout)               :: parameters
    class           (outputAnalysisPropertyExtractorClass     ), pointer                     :: outputAnalysisPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass      ), pointer                     :: outputAnalysisPropertyOperator_      , outputAnalysisPropertyUnoperator_
    class           (outputAnalysisWeightOperatorClass        ), pointer                     :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass  ), pointer                     :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass), pointer                     :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                      ), pointer                     :: galacticFilter_
    double precision                                           , dimension(:  ), allocatable :: binCenter                            , outputWeight
    integer         (c_size_t                                 )                              :: bufferCount
    type            (varying_string                           )                              :: label                                , comment                          , &
         &                                                                                      propertyLabel                        , propertyComment                  , &
         &                                                                                      distributionLabel                    , distributionComment              , &
         &                                                                                      propertyUnits                        , distributionUnits                , &
         &                                                                                      covarianceModel
    integer                                                                                  :: covarianceBinomialBinsPerDecade
    type            (inputParameters                          )                              :: unoperatorParameters
    double precision                                                                         :: propertyUnitsInSI                    , distributionUnitsInSI            , &
         &                                                                                      covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum
    
    ! Check and read parameters.
    unoperatorParameters=parameters%subParameters('unoperatorParameters',requireValue=.false.)
    call allocateArray(binCenter   ,[int(parameters%count('binCenter'),kind=c_size_t)                               ])
    call allocateArray(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t)*Galacticus_Output_Time_Count()])
    if (parameters%count('outputWeight') /= parameters%count('binCenter')*Galacticus_Output_Time_Count()) &
         & call Galacticus_Error_Report('incorrect number of output weights provided'//{introspection:location})
    !# <inputParameter>
    !#   <name>label</name>
    !#   <source>parameters</source>
    !#   <variable>label</variable>
    !#   <description>A label for the analysis.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>comment</name>
    !#   <source>parameters</source>
    !#   <variable>comment</variable>
    !#   <description>A descriptive comment for the analysis.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>propertyLabel</name>
    !#   <source>parameters</source>
    !#   <variable>propertyLabel</variable>
    !#   <description>A label for the property variable.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>propertyComment</name>
    !#   <source>parameters</source>
    !#   <variable>propertyComment</variable>
    !#   <description>A descriptive comment for the property variable.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>propertyUnits</name>
    !#   <source>parameters</source>
    !#   <variable>propertyUnits</variable>
    !#   <description>A human-readable description of the units for the property.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>propertyUnitsInSI</name>
    !#   <source>parameters</source>
    !#   <variable>propertyUnitsInSI</variable>
    !#   <description>A units for the property in the SI system.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>distributionLabel</name>
    !#   <source>parameters</source>
    !#   <variable>distributionLabel</variable>
    !#   <description>A label for the distribution.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>distributionComment</name>
    !#   <source>parameters</source>
    !#   <variable>distributionComment</variable>
    !#   <description>A descriptive comment for the distribution.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>distributionUnits</name>
    !#   <source>parameters</source>
    !#   <variable>distributionUnits</variable>
    !#   <description>A human-readable description of the units for the distribution.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>distributionUnitsInSI</name>
    !#   <source>parameters</source>
    !#   <variable>distributionUnitsInSI</variable>
    !#   <description>A units for the distribution in the SI system.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>    
    !# <inputParameter>
    !#   <name>binCenter</name>
    !#   <source>parameters</source>
    !#   <variable>binCenter</variable>
    !#   <description>The value of the property at the center of each bin.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bufferCount</name>
    !#   <source>parameters</source>
    !#   <variable>bufferCount</variable>
    !#   <description>The number of buffer bins to include below and above the range of actual bins.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputWeight</name>
    !#   <source>parameters</source>
    !#   <variable>outputWeight</variable>
    !#   <description>The weight to assign to each bin at each output.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceModel</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceModel</variable>
    !#   <description>The model to use for computing covariances.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialBinsPerDecade</name>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins per decade of halo mass to use when constructing volume function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing volume function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="outputAnalysisPropertyExtractor"      name="outputAnalysisPropertyExtractor_"      source="parameters"          />
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"       source="parameters"          />
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyUnoperator_"     source="unoperatorParameters"/>
    !# <objectBuilder class="outputAnalysisWeightOperator"         name="outputAnalysisWeightOperator_"         source="parameters"          />
    !# <objectBuilder class="outputAnalysisDistributionOperator"   name="outputAnalysisDistributionOperator_"   source="parameters"          />
    !# <objectBuilder class="outputAnalysisDistributionNormalizer" name="outputAnalysisDistributionNormalizer_" source="parameters"          />
    !# <objectBuilder class="galacticFilter"                       name="galacticFilter_"                       source="parameters"          />
    ! Build the object.
    self=outputAnalysisVolumeFunction1D(                                                                                                         &
         &                              label                                                                                                  , &
         &                              comment                                                                                                , &
         &                              propertyLabel                                                                                          , &
         &                              propertyComment                                                                                        , &
         &                              propertyUnits                                                                                          , &
         &                              propertyUnitsInSI                                                                                      , &
         &                              distributionLabel                                                                                      , &
         &                              distributionComment                                                                                    , &
         &                              distributionUnits                                                                                      , &
         &                              distributionUnitsInSI                                                                                  , &
         &                              binCenter                                                                                              , &
         &                              bufferCount                                                                                            , &
         &                              reshape(outputWeight,[int(parameters%count('binCenter'),kind=c_size_t),Galacticus_Output_Time_Count()]), &
         &                              outputAnalysisPropertyExtractor_                                                                       , &
         &                              outputAnalysisPropertyOperator_                                                                        , &
         &                              outputAnalysisPropertyUnoperator_                                                                      , &
         &                              outputAnalysisWeightOperator_                                                                          , &
         &                              outputAnalysisDistributionOperator_                                                                    , &
         &                              outputAnalysisDistributionNormalizer_                                                                  , &
         &                              galacticFilter_                                                                                        , &
         &                              enumerationOutputAnalysisCovarianceModelEncode(char(covarianceModel),includesPrefix=.false.)           , &
         &                              covarianceBinomialBinsPerDecade                                                                        , &
         &                              covarianceBinomialMassHaloMinimum                                                                      , &
         &                              covarianceBinomialMassHaloMaximum                                                                        &                      
         &                             )
    !# <inputParametersValidate source="parameters"/>
    return
  end function volumeFunction1DConstructorParameters

  function volumeFunction1DConstructorInternal(label,comment,propertyLabel,propertyComment,propertyUnits,propertyUnitsInSI,distributionLabel,distributionComment,distributionUnits,distributionUnitsInSI,binCenter,bufferCount,outputWeight,outputAnalysisPropertyExtractor_,outputAnalysisPropertyOperator_,outputAnalysisPropertyUnoperator_,outputAnalysisWeightOperator_,outputAnalysisDistributionOperator_,outputAnalysisDistributionNormalizer_,galacticFilter_,covarianceModel,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum) result (self)
    !% Constructor for the ``volumeFunction1D'' output analysis class for internal use.
    use Memory_Management
    implicit none
    type            (outputAnalysisVolumeFunction1D           )                                :: self
    type            (varying_string                           ), intent(in   )                 :: label                                , comment                          , &
         &                                                                                        propertyLabel                        , propertyComment                  , &
         &                                                                                        distributionLabel                    , distributionComment              , &
         &                                                                                        propertyUnits                        , distributionUnits
    double precision                                                                           :: propertyUnitsInSI                    , distributionUnitsInSI
    double precision                                           , intent(in   ), dimension(:  ) :: binCenter
    integer         (c_size_t                                 ), intent(in   )                 :: bufferCount
    double precision                                           , intent(in   ), dimension(:,:) :: outputWeight
    class           (outputAnalysisPropertyExtractorClass     ), intent(in   ), target         :: outputAnalysisPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass      ), intent(in   ), target         :: outputAnalysisPropertyOperator_      , outputAnalysisPropertyUnoperator_
    class           (outputAnalysisWeightOperatorClass        ), intent(in   ), target         :: outputAnalysisWeightOperator_
    class           (outputAnalysisDistributionOperatorClass  ), intent(in   ), target         :: outputAnalysisDistributionOperator_
    class           (outputAnalysisDistributionNormalizerClass), intent(in   ), target         :: outputAnalysisDistributionNormalizer_
    class           (galacticFilterClass                      ), intent(in   ), target         :: galacticFilter_
    integer                                                    , intent(in   )                 :: covarianceModel
    integer                                                    , intent(in   ), optional       :: covarianceBinomialBinsPerDecade
    double precision                                           , intent(in   ), optional       :: covarianceBinomialMassHaloMinimum    , covarianceBinomialMassHaloMaximum
    integer         (c_size_t                                 )                                :: i
    !# <constructorAssign variables="label, comment, propertyLabel, propertyComment, propertyUnits, propertyUnitsInSI, distributionLabel, distributionComment, distributionUnits, distributionUnitsInSI, binCenter, bufferCount, outputWeight, *outputAnalysisPropertyExtractor_, *outputAnalysisPropertyOperator_, *outputAnalysisPropertyUnoperator_, *outputAnalysisWeightOperator_, *outputAnalysisDistributionOperator_, *outputAnalysisDistributionNormalizer_, *galacticFilter_, covarianceModel, covarianceBinomialBinsPerDecade, covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum"/>

    ! Count bins.
    self%binCount     =size(binCenter,kind=c_size_t)
    self%binCountTotal=self%binCount+2*bufferCount
    ! Determine bin minima and maxima. Allocate arrays for bins that include buffer regions.
    call allocateArray(self%binMinimum,dimensions=[self%binCountTotal],lowerBounds=[-bufferCount+1])
    call allocateArray(self%binMaximum,dimensions=[self%binCountTotal],lowerBounds=[-bufferCount+1])
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
    ! Allocate and initialize function values.
    call allocateArray(self%functionValue     ,[self%binCount              ])
    call allocateArray(self%functionCovariance,[self%binCount,self%binCount])
    self%functionValue     =0.0d0
    self%functionCovariance=0.0d0
    ! Allocate and initialize binomial covariance model halo arrays if necessary.
    if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       self%covarianceModelBinomialBinCount                  =int(log10(self%covarianceBinomialMassHaloMaximum/self%covarianceBinomialMassHaloMinimum)*dble(self%covarianceBinomialBinsPerDecade)+0.5d0)
       self%covarianceModelHaloMassMinimumLogarithmic        =log10(self%covarianceBinomialMassHaloMinimum)
       self%covarianceModelHaloMassIntervalLogarithmicInverse=dble(self%covarianceModelBinomialBinCount)/log10(self%covarianceBinomialMassHaloMaximum/self%covarianceBinomialMassHaloMinimum)
       call allocateArray(self%weightMainBranch       ,[self%binCount,self%covarianceModelBinomialBinCount])
       call allocateArray(self%weightSquaredMainBranch,[self%binCount,self%covarianceModelBinomialBinCount])
       self%weightMainBranch       =0.0d0
       self%weightSquaredMainBranch=0.0d0
    end if
    ! Initialize finalization status.
    self%finalized=.false.
    ! Initialize OpenMP accumulation lock.
    !$ call OMP_Init_Lock(self%accumulateLock)
   return
  end function volumeFunction1DConstructorInternal

  subroutine volumeFunction1DDestructor(self)
    !% Destructor for  the ``volumeFunction1D'' output analysis class.
    type(outputAnalysisVolumeFunction1D), intent(inout) :: self
    
    !# <objectDestructor name="self%outputAnalysisPropertyExtractor_"     />
    !# <objectDestructor name="self%outputAnalysisPropertyOperator_"      />
    !# <objectDestructor name="self%outputAnalysisPropertyUnoperator_"    />
    !# <objectDestructor name="self%outputAnalysisWeightOperator_"        />
    !# <objectDestructor name="self%outputAnalysisDistributionOperator_"  />
    !# <objectDestructor name="self%outputAnalysisDistributionNormalizer_"/>
    !# <objectDestructor name="self%galacticFilter_"                      />
    ! Destroy OpenMP lock.
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine volumeFunction1DDestructor
  
  subroutine volumeFunction1DAnalyze(self,node,iOutput)
    !% Implement a volumeFunction1D output analysis.
    implicit none
    class           (outputAnalysisVolumeFunction1D), intent(inout)                 :: self
    type            (treeNode                      ), intent(inout)                 :: node
    integer         (c_size_t                      ), intent(in   )                 :: iOutput
    double precision                                , allocatable  , dimension(:  ) :: distribution
    double precision                                , allocatable  , dimension(:,:) :: covariance
    class           (nodeComponentBasic            ), pointer                       :: basic
    double precision                                                                :: propertyValue         , weightValue  , &
         &                                                                             propertyValueIntrinsic
    integer                                                                         :: propertyType          , propertyQuantity
    integer         (c_size_t                      )                                :: j                     , k            , &
         &                                                                             indexHaloMass

    ! If weights for this output are all zero, we can skip analysis.
    if (all(self%outputWeight(:,iOutput) == 0.0d0)) return
    ! Filter this node.
    if (.not.self%galacticFilter_%passes(node)) return
    ! Allocate work arrays.
    allocate(distribution(-self%bufferCount+1:self%binCount+self%bufferCount              ))
    allocate(covariance  (                    self%binCount                 ,self%binCount))
    ! Extract the property from the node.
    propertyType          =self%outputAnalysisPropertyExtractor_%type    (    )
    propertyQuantity      =self%outputAnalysisPropertyExtractor_%quantity(    )
    propertyValue         =self%outputAnalysisPropertyExtractor_%extract (node)
    propertyValueIntrinsic=propertyValue
    ! Apply property operators.
    propertyValue=self%outputAnalysisPropertyOperator_%operate(propertyValue,node,propertyType,iOutput)
    ! Apply distribution operators.
    distribution=self%outputAnalysisDistributionOperator_%operateScalar(propertyValue,propertyType,self%binMinimum,self%binMaximum,iOutput,node)
    ! Compute the weight.
    weightValue=node%hostTree%volumeWeight
    ! Apply weight operators.
    weightValue=self%outputAnalysisWeightOperator_%operate(weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,iOutput)
    ! Apply weights.
    distribution(1:self%binCount)=+distribution              (1:self%binCount        ) &
         &                        *self        %outputWeight ( :             ,iOutput) &
         &                        *weightValue
    ! Accumulate the property, including weights from both the host tree and the output. Note that we accumulate only the
    ! non-buffer bins of the distribution.
    !$ call OMP_Set_Lock(self%accumulateLock)
    self%functionValue=+self%functionValue+distribution(1:self%binCount)
    !$ call OMP_Unset_Lock(self%accumulateLock)
    ! Accumulate covariance. If using the binomial model for main branch galaxies, handle them separately.
    if (node%isOnMainBranch() .and. self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       ! Find the bin to which this halo mass belongs.
       basic         => node%basic()
       indexHaloMass =  floor((log10(basic%mass())-self%covarianceModelHaloMassMinimumLogarithmic)*self%covarianceModelHaloMassIntervalLogarithmicInverse)+1
       ! Accumulate weights to halo mass arrays.
       if (indexHaloMass >= 1 .and. indexHaloMass <= self%covarianceModelBinomialBinCount) then
          self                 %weightMainBranch       ( :             ,indexHaloMass)=   &
               &  +self        %weightMainBranch       ( :             ,indexHaloMass)    &
               &  +distribution                        (1:self%binCount)
          self                 %weightSquaredMainBranch( :             ,indexHaloMass)=   &
               &  +self        %weightSquaredMainBranch( :             ,indexHaloMass)    &
               &  +distribution                        (1:self%binCount              )**2
       end if
    else
       ! Construct contribution to the covariance matrix assuming Poisson statistics.
       forall(j=1:self%binCount)
          forall(k=j:self%binCount)
             covariance(j,k)=+distribution(j  ) &
                  &          *distribution(  k)
             covariance(k,j)=+covariance  (j,k)
          end forall
       end forall
       ! Accumulate covariance.
       !$ call OMP_Set_Lock(self%accumulateLock)
       self        %functionCovariance= &
            & +self%functionCovariance  &
            & +             covariance
       !$ call OMP_Unset_Lock(self%accumulateLock)
    end if
    ! Deallocate workspace.
    deallocate(distribution)
    deallocate(covariance  )
    return
  end subroutine volumeFunction1DAnalyze

  subroutine volumeFunction1DFinalize(self)
    !% Implement a volumeFunction1D output analysis finalization.
    use IO_HDF5
    use Galacticus_HDF5
    implicit none
    class(outputAnalysisVolumeFunction1D), intent(inout) :: self
    type (hdf5Object                    )                :: analysesGroup        , analysisGroup, &
         &                                                  dataset

    ! Finalize analysis.
    call volumeFunction1DFinalizeAnalysis(self)
    ! Output.
    !$omp critical(HDF5_Access)
    analysesGroup=galacticusOutputFile%openGroup('analyses'                         )
    analysisGroup=analysesGroup       %openGroup(char(self%label),char(self%comment))
    call analysisGroup%writeDataset  (self%binCenter    (1:self%binCount                     ),char(self%    propertyLabel)              ,char(self%    propertyComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%    propertyUnits    )      ,'units'                                                                                                                        )
    call dataset      %writeAttribute(          self%    propertyUnitsInSI       ,'unitsInSI'                                                                                                                    )
    call dataset      %close()
    call analysisGroup%writeDataset  (self%functionValue(1:self%binCount                     ),char(self%distributionLabel)              ,char(self%distributionComment)                 ,datasetReturned=dataset)
    call dataset      %writeAttribute(     char(self%distributionUnits    )      ,'units'                                                                                                                        )
    call dataset      %writeAttribute(          self%distributionUnitsInSI       ,'unitsInSI'                                                                                                                    )
    call dataset      %close()
    call analysisGroup%writeDataset  (self%functionCovariance(1:self%binCount,1:self%binCount),char(self%distributionLabel)//"Covariance",char(self%distributionComment)//" [covariance]",datasetReturned=dataset)
    call dataset      %writeAttribute("["//char(self%distributionUnits    )//"]²",'units'                                                                                                                        )
    call dataset      %writeAttribute(          self%distributionUnitsInSI   **2 ,'unitsInSI'                                                                                                                    )
    call dataset      %close()
    call analysisGroup%close()
    call analysesGroup%close()
    !$omp end critical(HDF5_Access)    
    return
  end subroutine volumeFunction1DFinalize

  subroutine volumeFunction1DFinalizeAnalysis(self)
    !% Compute final covariances and normalize.
    implicit none
    class           (outputAnalysisVolumeFunction1D), intent(inout) :: self
    double precision                                , parameter     :: weightFractionMaximum=0.99d0
    integer         (c_size_t                      )                :: i                           , j, &
         &                                                             m
    double precision                                                :: weightMainBranchTotal

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
    ! If using the binomial covariance model, add on the contribution to covariance from main branch galaxies.
    if (self%covarianceModel == outputAnalysisCovarianceModelBinomial) then
       do m=1,self%covarianceModelBinomialBinCount
          weightMainBranchTotal=sum(self%weightMainBranch(:,m))
          if (weightMainBranchTotal > 0.0d0) then
             do i=1,self%binCount
                if (self%weightMainBranch(i,m) < weightFractionMaximum*weightMainBranchTotal) then
                   ! General case - multiple halo mass bins contributed to this bin of the volume function.
                   self               %functionCovariance     (i,i)=                       &
                        &         self%functionCovariance     (i,i)                        &
                        & +(1.0d0-self%weightMainBranch       (i,m)/weightMainBranchTotal) &
                        & *       self%weightSquaredMainBranch(i,m)
                else
                   ! Special case - only one halo mass bin contributed to this bin of the volume function. We revert to Poisson
                   ! statistics in this case to avoid having zero covariance.
                   self               %functionCovariance     (i,i)=                       &
                        &         self%functionCovariance     (i,i)                        &
                        & +       self%weightSquaredMainBranch(i,m)
                end if
                do j=1,self%binCount
                   if (i == j) cycle
                   self               %functionCovariance     (i,j)= &
                        &  +      self%functionCovariance     (i,j)  &
                        &  -      self%weightMainBranch       (i,m)  &
                        &  *      self%weightMainBranch       (j,m)  &
                        &  *sqrt(                                    &
                        &         self%weightSquaredMainBranch(i,m)  &
                        &        *self%weightSquaredMainBranch(j,m)  &
                        &       )                                    &
                        &  /weightMainBranchTotal
                end do
             end do
          end if
       end do
    end if
    ! Apply final distribution operators - pass only the non-buffer bin values here.
    call self%outputAnalysisDistributionNormalizer_%normalize(self%functionValue,self%functionCovariance,self%binMinimum(1:self%binCount),self%binMaximum(1:self%binCount))
    ! Apply any "unoperator" to output bin values. This can be used to reverse transformations (e.g. if masses were converted to
    ! log10 for analysis, that can be undone here).
    do i=1,self%binCount
       self%binCenter(i)=self%outputAnalysisPropertyUnoperator_%operate(self%binCenter(i))
    end do
    return
  end subroutine volumeFunction1DFinalizeAnalysis
    
  subroutine volumeFunction1DResults(self,binCenter,functionValue,functionCovariance)
    !% Implement a volumeFunction1D output analysis finalization.
    use Memory_Management
    implicit none
    class           (outputAnalysisVolumeFunction1D)                             , intent(inout)           :: self
    double precision                                , allocatable, dimension(:  ), intent(inout), optional :: binCenter         , functionValue
    double precision                                , allocatable, dimension(:,:), intent(inout), optional :: functionCovariance

    ! Finalize analysis.
    call volumeFunction1DFinalizeAnalysis(self)
    ! Return results.
    if (present(binCenter         )) then
       if (allocated(binCenter         )) call deallocateArray(binCenter         )
       call allocateArray(binCenter         ,shape(self%binCenter         ))
       binCenter         =self%binCenter
    end if
    if (present(functionValue     )) then
       if (allocated(functionValue     )) call deallocateArray(functionValue     )
       call allocateArray(functionValue     ,shape(self%functionValue     ))
       functionValue     =self%functionValue
    end if
    if (present(functionCovariance)) then
       if (allocated(functionCovariance)) call deallocateArray(functionCovariance)
       call allocateArray(functionCovariance,shape(self%functionCovariance))
       functionCovariance=self%functionCovariance
    end if
    return
  end subroutine volumeFunction1DResults
