!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements a spin parameter distribution output analysis class.

  use Cosmology_Functions
  use Dark_Matter_Halo_Scales
  use Dark_Matter_Profiles
  use Halo_Mass_Functions
  use Statistics_NBody_Halo_Mass_Errors
  use Dark_Matter_Profile_Scales       , only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass

  !# <outputAnalysis name="outputAnalysisSpinDistribution" defaultThreadPrivate="yes">
  !#  <description>A stellar mass function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisSpinDistribution
     !% A spinDistribution output analysis class.
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_
   contains
     final :: spinDistributionDestructor
  end type outputAnalysisSpinDistribution

  interface outputAnalysisSpinDistribution
     !% Constructors for the ``spinDistribution'' output analysis class.
     module procedure spinDistributionConstructorParameters
     module procedure spinDistributionConstructorInternal
  end interface outputAnalysisSpinDistribution

contains

  function spinDistributionConstructorParameters(parameters) result (self)
    !% Constructor for the ``spinDistribution'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisSpinDistribution   )                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass                 ), pointer       :: outputTimes_
    class           (nbodyHaloMassErrorClass          ), pointer       :: nbodyHaloMassError_ 
    class           (haloMassFunctionClass            ), pointer       :: haloMassFunction_            
    class           (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_         
    class           (darkMatterProfileClass           ), pointer       :: darkMatterProfile_
    class           (darkMatterProfileScaleRadiusClass), pointer       :: darkMatterProfileScaleRadius_
    double precision                                                   :: timeRecent
    
    !# <inputParameter>
    !#   <name>timeRecent</name>
    !#   <source>parameters</source>
    !#   <description>Halos which experienced a major node merger within a time $\Delta t=${\normalfont \ttfamily [timeRecent]} of the analysis time will be excluded from the analysis.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    !# <objectBuilder class="outputTimes"                  name="outputTimes_"                  source="parameters"/>
    !# <objectBuilder class="nbodyHaloMassError"           name="nbodyHaloMassError_"           source="parameters"/>
    !# <objectBuilder class="haloMassFunction"             name="haloMassFunction_"             source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    !# <objectBuilder class="darkMatterProfile"            name="darkMatterProfile_"            source="parameters"/>
    !# <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    self=outputAnalysisSpinDistribution(timeRecent,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfile_,darkMatterProfileScaleRadius_,outputTimes_)
    !# <inputParametersValidate source="parameters"/>
    nullify(nbodyHaloMassError_          )
    nullify(haloMassFunction_            )
    nullify(darkMatterHaloScale_         )
    nullify(darkMatterProfile_           )
    nullify(darkMatterProfileScaleRadius_)
    return
  end function spinDistributionConstructorParameters

  function spinDistributionConstructorInternal(timeRecent,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfile_,darkMatterProfileScaleRadius_,outputTimes_) result(self)
    !% Internal constructor for the ``spinDistribution'' output analysis class.
    use ISO_Varying_String
    use Output_Times
    use Output_Analyses_Options
    use Galacticus_Error
    use Galacticus_Paths
    use IO_HDF5
    use Memory_Management
    use Numerical_Comparison
    use Halo_Spin_Distributions
    implicit none
    type            (outputAnalysisSpinDistribution                   )                              :: self
    double precision                                                                , intent(in   )  :: timeRecent
    class           (cosmologyFunctionsClass                          ), target     , intent(in   )  :: cosmologyFunctions_
    class           (outputTimesClass                                 ), target     , intent(inout)  :: outputTimes_
    class           (nbodyHaloMassErrorClass                          ), target     , intent(in   )  :: nbodyHaloMassError_ 
    class           (haloMassFunctionClass                            ), target     , intent(in   )  :: haloMassFunction_            
    class           (darkMatterHaloScaleClass                         ), target     , intent(in   )  :: darkMatterHaloScale_         
    class           (darkMatterProfileClass                           ), target     , intent(in   )  :: darkMatterProfile_           
    class           (darkMatterProfileScaleRadiusClass                ), target     , intent(in   )  :: darkMatterProfileScaleRadius_           
    type            (outputAnalysisPropertyExtractorSpin              ), pointer                     :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyOperatorLog10              ), pointer                     :: outputAnalysisPropertyOperator_
    type            (outputAnalysisWeightOperatorIdentity             ), pointer                     :: outputAnalysisWeightOperator_
    type            (outputAnalysisDistributionNormalizerSequence     ), pointer                     :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisPropertyOperatorAntiLog10          ), pointer                     :: outputAnalysisPropertyOperatorAntiLog10_
    type            (outputAnalysisDistributionOperatorSpinNBodyErrors), pointer                     :: outputAnalysisDistributionOperator_
    type            (normalizerList                                   ), pointer                     :: normalizerSequence                             , normalizer_
    type            (propertyOperatorList                             ), pointer                     :: propertyOperatorSequence
    type            (galacticFilterHaloIsolated                       ), pointer                     :: galacticFilterHaloIsolated_
    type            (galacticFilterNodeMajorMergerRecent              ), pointer                     :: galacticFilterNodeMajorMergerRecent_
    type            (galacticFilterBasicMass                          ), pointer                     :: galacticFilterBasicMass_
    type            (galacticFilterNot                                ), pointer                     :: galacticFilterNot_
    type            (galacticFilterAll                                ), pointer                     :: galacticFilterAll_
    type            (filterList                                       ), pointer                     :: filters_                                       , filter_
    type            (haloSpinDistributionNbodyErrors                  ), pointer                     :: haloSpinDistribution_
    type            (haloSpinDistributionDeltaFunction                ), pointer                     :: haloSpinDistributionDeltaFunction_
    double precision                                                   , allocatable, dimension(:  ) :: spins
    double precision                                                   , allocatable, dimension(:,:) :: outputWeight
    double precision                                                   , parameter                   :: massParticleMillennium                  =1.178d9
    double precision                                                   , parameter                   :: countParticleMininum                    =3.000d2
    integer         (c_size_t                                         ), parameter                   :: bufferCountMinimum                      =5
    double precision                                                   , parameter                   :: bufferWidthLogarithmic                  =0.5d0
    integer                                                            , parameter                   :: covarianceBinomialBinsPerDecade         =2
    double precision                                                   , parameter                   :: covarianceBinomialMassHaloMinimum       =3.0d11, covarianceBinomialMassHaloMaximum=1.0d15
    integer         (c_size_t                                         )                              :: iOutput                                        , bufferCount
    type            (hdf5Object                                       )                              :: dataFile
    !# <constructorAssign variables="*cosmologyFunctions_"/>
    
    ! Construct spins matched to those used by Bett et al. (2007).
    !$ call hdf5Access%set()
    call dataFile%openFile   (char(galacticusPath(pathTypeDataStatic)//'darkMatter/bett2007HaloSpinDistribution.hdf5'),readOnly=.true.)
    call dataFile%readDataset(                              'spinParameter'                                     ,         spins )
    call dataFile%close      (                                                                                                  )
    !$ call hdf5Access%unset()
    self%binCount=size(spins)
    ! Compute weights that apply to each output redshift.
    call allocateArray(outputWeight,[self%binCount,outputTimes_%count()])
    outputWeight=0.0d0
    do iOutput=1,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(iOutput),0.0d0,absTol=1.0d-10)) outputWeight(:,iOutput)=1.0d0
    end do
    if (any(sum(outputWeight,dim=2) /= 1.0d0)) call Galacticus_Error_Report('zero redshift output is required'//{introspection:location})
    ! Build an N-body halo spin distribution class.
    allocate(haloSpinDistributionDeltaFunction_)
    allocate(haloSpinDistribution_             )
    haloSpinDistributionDeltaFunction_=  haloSpinDistributionDeltaFunction(                                                                                    &
         &                                                                 spin                              =                                       0.0d0     &
         &                                                                )
    haloSpinDistribution_             =  haloSpinDistributionNbodyErrors  (                                                                                    &
         &                                                                                                    haloSpinDistributionDeltaFunction_             , &
         &                                                                 massParticle                      =massParticleMillennium                         , &
         &                                                                 particleCountMinimum              =                                     300       , &
         &                                                                 energyEstimateParticleCountMaximum=                                    1000.000d0 , &
         &                                                                 time                              =self%cosmologyFunctions_%cosmicTime(   1.000d0), &
         &                                                                 nbodyHaloMassError_               =nbodyHaloMassError_                            , &
         &                                                                 haloMassFunction_                 =haloMassFunction_                              , &
         &                                                                 darkMatterHaloScale_              =darkMatterHaloScale_                           , &
         &                                                                 darkMatterProfile_                =darkMatterProfile_                             , &
         &                                                                 darkMatterProfileScaleRadius_     =darkMatterProfileScaleRadius_                    &
         &                                                                )
    ! Create a spin parameter property extractor.
    allocate(outputAnalysisPropertyExtractor_        )
    outputAnalysisPropertyExtractor_        =outputAnalysisPropertyExtractorSpin              (                                           )
    ! Create a log10 property operator.
    allocate(outputAnalysisPropertyOperator_         )
    outputAnalysisPropertyOperator_         =outputAnalysisPropertyOperatorLog10              (                                           )
    ! Create an identity weight operator.
    allocate(outputAnalysisWeightOperator_           )
    outputAnalysisWeightOperator_           =outputAnalysisWeightOperatorIdentity             (                                           )
    ! Create an N-body spin error distribution operator.
    allocate(outputAnalysisDistributionOperator_     )
    outputAnalysisDistributionOperator_     =outputAnalysisDistributionOperatorSpinNBodyErrors(haloSpinDistribution_                      )
    ! Create anit-log10 operator.
    allocate(outputAnalysisPropertyOperatorAntiLog10_)
    outputAnalysisPropertyOperatorAntiLog10_=outputAnalysisPropertyOperatorAntiLog10          (                                           )
    ! Create a filter to select isolated halos with no recent major merger.
    allocate(galacticFilterHaloIsolated_             )
    galacticFilterHaloIsolated_             =galacticFilterHaloIsolated                       (                                           )
    allocate(galacticFilterBasicMass_                )
    galacticFilterBasicMass_                =galacticFilterBasicMass                          (countParticleMininum*massParticleMillennium)
    allocate(galacticFilterNodeMajorMergerRecent_    )
    galacticFilterNodeMajorMergerRecent_    =galacticFilterNodeMajorMergerRecent              (timeRecent                                 )
    allocate(galacticFilterNot_                      )
    galacticFilterNot_                      =galacticFilterNot                                (galacticFilterNodeMajorMergerRecent_       )
    allocate(filters_                                )
    filter_ => filters_
    filter_%filter_ => galacticFilterHaloIsolated_
    allocate(filter_%next)
    filter_ => filter_%next
    filter_%filter_ => galacticFilterBasicMass_
    allocate(filter_%next)
    filter_ => filter_%next
    filter_%filter_ => galacticFilterNot_
    allocate(galacticFilterAll_                      )
    galacticFilterAll_                      =galacticFilterAll                                (filters_                            )
    ! Create a distribution normalizer which normalizes to unit integral, and then to bin width.
    allocate(normalizerSequence)
    normalizer_ => normalizerSequence
    allocate(outputAnalysisDistributionNormalizerUnitarity   :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerUnitarity )
       normalizer_=outputAnalysisDistributionNormalizerUnitarity ()
    end select
    allocate(normalizer_%next)
    normalizer_ => normalizer_%next
    allocate(outputAnalysisDistributionNormalizerBinWidth   :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerBinWidth  )
       normalizer_=outputAnalysisDistributionNormalizerBinWidth  ()
    end select
    allocate(normalizer_%next)
    normalizer_ => normalizer_%next
    allocate(outputAnalysisDistributionNormalizerLog10ToLog :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerLog10ToLog)
       normalizer_=outputAnalysisDistributionNormalizerLog10ToLog()
    end select
    allocate(outputAnalysisDistributionNormalizer_)
    outputAnalysisDistributionNormalizer_=outputAnalysisDistributionNormalizerSequence(normalizerSequence)
    ! Compute the number of buffer bins to add to either side of the mass function - these are needed to ensure that, e.g.,
    ! convolution operations on the distribution function are unaffected by edge effects.
    bufferCount=max(int(bufferWidthLogarithmic/log10(spins(2)/spins(1)))+1,bufferCountMinimum)
    ! Construct the object. We convert spins to log10(spins) here.
    self%outputAnalysisVolumeFunction1D=                                                       &
         & outputAnalysisVolumeFunction1D(                                                     &
         &                                var_str('spinDistribution')                        , &
         &                                var_str('Distribution of halo spin parameters'    ), &
         &                                var_str('spin'                                    ), &
         &                                var_str('Spin at the bin center'                  ), &
         &                                var_str('dimensionless'                           ), &
         &                                0.0d0                                              , &
         &                                var_str('spinDistributionFunction'                ), &
         &                                var_str('Spin distribution averaged over each bin'), &
         &                                var_str('dimensionless'                           ), &
         &                                0.0d0                                              , &
         &                                log10(spins)                                       , &
         &                                bufferCount                                        , &
         &                                outputWeight                                       , &
         &                                outputAnalysisPropertyExtractor_                   , &
         &                                outputAnalysisPropertyOperator_                    , &
         &                                outputAnalysisPropertyOperatorAntiLog10_           , &
         &                                outputAnalysisWeightOperator_                      , &
         &                                outputAnalysisDistributionOperator_                , &
         &                                outputAnalysisDistributionNormalizer_              , &
         &                                galacticFilterAll_                                 , &
         &                                outputTimes_                                       , &
         &                                outputAnalysisCovarianceModelPoisson               , &
         &                                covarianceBinomialBinsPerDecade                    , &
         &                                covarianceBinomialMassHaloMinimum                  , &
         &                                covarianceBinomialMassHaloMaximum                    &
         &                               )
    ! Clean up.
    nullify(outputAnalysisPropertyExtractor_        )
    nullify(outputAnalysisPropertyOperator_         )
    nullify(outputAnalysisPropertyOperatorAntiLog10_)
    nullify(outputAnalysisDistributionNormalizer_   )
    nullify(outputAnalysisWeightOperator_           )
    nullify(propertyOperatorSequence                )
    nullify(normalizerSequence                      )
    nullify(haloSpinDistribution_                   )
    nullify(haloSpinDistributionDeltaFunction_      )
    nullify(filters_                                )
    nullify(galacticFilterAll_                      )
    nullify(galacticFilterNot_                      )
    nullify(galacticFilterHaloIsolated_             )
    nullify(galacticFilterNodeMajorMergerRecent_    )
    nullify(galacticFilterBasicMass_                )
    return
  end function spinDistributionConstructorInternal

  subroutine spinDistributionDestructor(self)
    !% Destructor for the ``spinDistribution'' output analysis class.
    implicit none
    type(outputAnalysisSpinDistribution), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologyFunctions_" />    
    return
  end subroutine spinDistributionDestructor
