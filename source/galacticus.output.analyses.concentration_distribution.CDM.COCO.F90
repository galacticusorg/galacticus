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

  !% Contains a module which implements a concentration distribution output analysis class for COCO CDM data.

  !# <outputAnalysis name="outputAnalysisConcentrationDistributionCDMCOCO" defaultThreadPrivate="yes">
  !#  <description>A concentration distribution function output analysis class for COCO CDM data.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisConcentrationDistributionCDMCOCO
     !% A concentration distribution output analysis class for COCO CDM data.
     private
     integer :: distributionNumber
   contains
  end type outputAnalysisConcentrationDistributionCDMCOCO

  interface outputAnalysisConcentrationDistributionCDMCOCO
     !% Constructors for the ``concentrationDistributionCDMCOCO'' output analysis class.
     module procedure concentrationDistributionCDMCOCOConstructorParameters
     module procedure concentrationDistributionCDMCOCOConstructorInternal
  end interface outputAnalysisConcentrationDistributionCDMCOCO

contains

  function concentrationDistributionCDMCOCOConstructorParameters(parameters) result (self)
    !% Constructor for the ``concentrationDistributionCDMCOCO'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    use Statistics_NBody_Halo_Mass_Errors
    use Cosmology_Functions
    implicit none
    type            (outputAnalysisConcentrationDistributionCDMCOCO)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctions_
    class           (nbodyHaloMassErrorClass                       ), pointer       :: nbodyHaloMassError_ 
    integer                                                                         :: distributionNumber
    double precision                                                                :: formationTimeRecent

    !# <inputParameter>
    !#   <name>distributionNumber</name>
    !#   <source>parameters</source>
    !#   <description>The number (1-7) of the distribution to compute.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>formationTimeRecent</name>
    !#   <source>parameters</source>
    !#   <description>Halos which ``formed'' more recently than this time in the past will be excluded from the analysis.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !# <objectBuilder class="nbodyHaloMassError" name="nbodyHaloMassError_" source="parameters"/>
    self=outputAnalysisConcentrationDistributionCDMCOCO(distributionNumber,formationTimeRecent,cosmologyFunctions_,nbodyHaloMassError_)
    !# <inputParametersValidate source="parameters"/>
    nullify(cosmologyFunctions_)
    nullify(nbodyHaloMassError_)
    return
  end function concentrationDistributionCDMCOCOConstructorParameters

  function concentrationDistributionCDMCOCOConstructorInternal(distributionNumber,formationTimeRecent,cosmologyFunctions_,nbodyHaloMassError_) result(self)
    !% Internal constructor for the ``concentrationDistributionCDMCOCO'' output analysis class.
    use ISO_Varying_String
    use Galacticus_Output_Times
    use Output_Analyses_Options
    use Output_Analysis_Utilities
    use Galacticus_Error
    use Statistics_NBody_Halo_Mass_Errors
    use Cosmology_Functions
    use Galacticus_Input_Paths
    use IO_HDF5
    use Memory_Management
    use Numerical_Comparison
    use Virial_Density_Contrast
    implicit none
    type            (outputAnalysisConcentrationDistributionCDMCOCO  )                              :: self
    class           (cosmologyFunctionsClass                         ), target     , intent(in   )  :: cosmologyFunctions_
    class           (nbodyHaloMassErrorClass                         ), target     , intent(in   )  :: nbodyHaloMassError_ 
    integer                                                                        , intent(in   )  :: distributionNumber
    double precision                                                               , intent(in   )  :: formationTimeRecent
    type            (outputAnalysisPropertyExtractorConcentration    ), pointer                     :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyOperatorLog10             ), pointer                     :: outputAnalysisPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10         ), pointer                     :: outputAnalysisPropertyUnoperator_
    type            (outputAnalysisWeightOperatorNbodyMass           ), pointer                     :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorIdentity          ), pointer                     :: outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyExtractorMassHalo         ), pointer                     :: outputAnalysisWeightPropertyExtractor_
    type            (outputAnalysisDistributionNormalizerSequence    ), pointer                     :: outputAnalysisDistributionNormalizer_
    type            (outputAnalysisDistributionOperatorRndmErrNbdyCnc), pointer                     :: outputAnalysisDistributionOperator_
    type            (normalizerList                                  ), pointer                     :: normalizerSequence                                  , normalizer_
    type            (galacticFilterHaloIsolated                      ), pointer                     :: galacticFilterHaloIsolated_
    type            (galacticFilterFormationTime                     ), pointer                     :: galacticFilterFormationTime_
    type            (galacticFilterAll                               ), pointer                     :: galacticFilter_
    type            (filterList                                      ), pointer                     :: filters_
    type            (virialDensityContrastFixed                      ), pointer                     :: virialDensityContrast_
    double precision                                                  , allocatable, dimension(:  ) :: concentrations                                      , masses
    double precision                                                  , allocatable, dimension(:,:) :: outputWeight
    double precision                                                  , parameter                   :: massParticle                            = +1.612d+05
    double precision                                                  , parameter                   :: haloDensityContrast                     = +2.000d+02
    double precision                                                  , parameter   , dimension(3)  :: concentrationFitA                       =[                                                    &
         &                                                                                                                                       -0.270d+00,                                         &
         &                                                                                                                                       +1.780d+00,                                         &
         &                                                                                                                                       -0.490d+00                                          &
         &                                                                                                                                      ]
    double precision                                                  , parameter                   :: concentrationFitB                       = -0.550d+00
    integer                                                           , parameter                   :: covarianceBinomialBinsPerDecade         =  2
    double precision                                                  , parameter                   :: covarianceBinomialMassHaloMinimum       = +3.000d+11, covarianceBinomialMassHaloMaximum=1.0d15
    integer         (c_size_t                                        )                              :: iOutput                                             , bufferCount               
    type            (hdf5Object                                      )                              :: dataFile
    double precision                                                                                :: massMinimum                                         , massMaximum
    character       (len=16                                          )                              :: distributionName
    !# <constructorAssign variables="distributionNumber"/>

    ! Validate input.
    if (distributionNumber < 1 .or. distributionNumber > 10) call Galacticus_Error_Report('distributionNumber ∈ [1..7] is required'//{introspection:location})
    !$omp critical(HDF5_Access)
    call dataFile%openFile   (char(Galacticus_Input_Path()//'data/darkMatter/concentrationDistributionCocoCDM.hdf5'),readOnly=.true.        )
    call dataFile%readDataset(                              'concentration'                                         ,         concentrations)
    call dataFile%readDataset(                              'mass'                                                  ,         masses        )
    call dataFile%close      (                                                                                                              )
    !$omp end critical(HDF5_Access)
    self%binCount=size(concentrations)
    ! Determine minimum and maximum halo masses for this distribution.
    massMinimum=masses(distributionNumber)/sqrt(masses(2)/masses(1))
    massMaximum=masses(distributionNumber)*sqrt(masses(2)/masses(1))
    ! Compute weights that apply to each output redshift.
    call allocateArray(outputWeight,[self%binCount,Galacticus_Output_Time_Count()])
    outputWeight=0.0d0
    do iOutput=1,Galacticus_Output_Time_Count()
       if (Values_Agree(Galacticus_Output_Redshift(iOutput),0.0d0,absTol=1.0d-10)) outputWeight(:,iOutput)=1.0d0
    end do
    ! Build a filter which selects isolated halos, and rejects halos which formed too recently.
    allocate(galacticFilter_                  )
    allocate(galacticFilterHaloIsolated_      )
    allocate(galacticFilterFormationTime_     )
    allocate(filters_                         )
    allocate(filters_                    %next)
    filters_                         %filter_ => galacticFilterHaloIsolated_
    filters_                    %next%filter_ => galacticFilterFormationTime_
    galacticFilterHaloIsolated_               =  galacticFilterHaloIsolated  (                   )
    galacticFilterFormationTime_              =  galacticFilterFormationTime (formationTimeRecent)
    galacticFilter_                           =  galacticFilterAll           (filters_           )
    ! Create a distribution normalizer which normalizes to bin width and unitarity.
    allocate(normalizerSequence)
    normalizer_ => normalizerSequence
    allocate(outputAnalysisDistributionNormalizerUnitarity  :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerUnitarity)
       normalizer_=outputAnalysisDistributionNormalizerUnitarity ()
    end select
    allocate(normalizer_%next )
    normalizer_ => normalizer_%next
    allocate(outputAnalysisDistributionNormalizerBinWidth   :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerBinWidth )
       normalizer_=outputAnalysisDistributionNormalizerBinWidth  ()
    end select
    allocate(normalizer_%next )
    normalizer_ => normalizer_%next
    allocate(outputAnalysisDistributionNormalizerLog10ToLog :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerLog10ToLog )
       normalizer_=outputAnalysisDistributionNormalizerLog10ToLog()
    end select
    allocate(outputAnalysisDistributionNormalizer_ )
    outputAnalysisDistributionNormalizer_ =outputAnalysisDistributionNormalizerSequence    (                                        &
         &                                                                                  normalizerSequence                      &
         &                                                                                 )
    ! Build log10() property operator.
    allocate(outputAnalysisPropertyOperator_       )
    outputAnalysisPropertyOperator_       =outputAnalysisPropertyOperatorLog10             (                                        &
         &                                                                                 )
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_     )
    outputAnalysisPropertyUnoperator_     =outputAnalysisPropertyOperatorAntiLog10         (                                        &
         &                                                                                 )
    ! Create a virial density contrast object matched to the defintion used by Ludlow et al. (2016).
    allocate(virialDensityContrast_                )
    virialDensityContrast_                =virialDensityContrastFixed                      (                                        &
         &                                                                                  haloDensityContrast                   , &
         &                                                                                  fixedDensityTypeCritical              , &
         &                                                                                  cosmologyFunctions_                     &
         &                                                                                 )
    ! Create a concentration property extractor.
    allocate(outputAnalysisPropertyExtractor_      )
    outputAnalysisPropertyExtractor_      =outputAnalysisPropertyExtractorConcentration    (                                        &
         &                                                                                  virialDensityContrast_                  &
         &                                                                                 )
    ! Create a halo mass property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_)
    outputAnalysisWeightPropertyExtractor_=outputAnalysisPropertyExtractorMassHalo         (                                        &
         &                                                                                  virialDensityContrast_                  &
         &                                                                                 )
    ! Create an identity property operator.
    allocate(outputAnalysisWeightPropertyOperator_ )
    outputAnalysisWeightPropertyOperator_ =outputAnalysisPropertyOperatorIdentity          (                                        &
         &                                                                                 )
    ! Build error distribution operator.
    allocate(outputAnalysisDistributionOperator_   )
    outputAnalysisDistributionOperator_   =outputAnalysisDistributionOperatorRndmErrNbdyCnc(                                        &
         &                                                                                  concentrationFitA                     , &
         &                                                                                  concentrationFitB                     , &
         &                                                                                  massParticle                          , &
         &                                                                                  outputAnalysisWeightPropertyExtractor_  &
         &                                                                                 )
    ! Build N-body mass distribution weight operator.
    allocate(outputAnalysisWeightOperator_         )
    outputAnalysisWeightOperator_          =outputAnalysisWeightOperatorNbodyMass          (                                        &
         &                                                                                  massMinimum                           , &
         &                                                                                  massMaximum                           , &
         &                                                                                  outputAnalysisWeightPropertyExtractor_, &
         &                                                                                  outputAnalysisWeightPropertyOperator_ , &
         &                                                                                  nbodyHaloMassError_                     &
         &                                                                                 )
    ! Determine number of buffer bins.
    bufferCount=0
    ! Construct the object.
    write (distributionName,'(i2.2)') distributionNumber
    self%outputAnalysisVolumeFunction1D=                                                                       &
         & outputAnalysisVolumeFunction1D(                                                                     &
         &                                var_str('concentrationDistributionCDMCOCO')//trim(distributionName), &
         &                                var_str('Distribution of concentration, c₂₀₀{crit}'               ), &
         &                                var_str('concentration'                                           ), &
         &                                var_str('Concentration at the bin center'                         ), &
         &                                var_str('dimensionless'                                           ), &
         &                                0.0d0                                                              , &
         &                                var_str('concentrationFunction'                                   ), &
         &                                var_str('Concentration distribution averaged over each bin'       ), &
         &                                var_str('dimensionless'                                           ), &
         &                                0.0d0                                                              , &
         &                                log10(concentrations)                                              , &
         &                                bufferCount                                                        , &
         &                                outputWeight                                                       , &
         &                                outputAnalysisPropertyExtractor_                                   , &
         &                                outputAnalysisPropertyOperator_                                    , &
         &                                outputAnalysisPropertyUnoperator_                                  , &
         &                                outputAnalysisWeightOperator_                                      , &
         &                                outputAnalysisDistributionOperator_                                , &
         &                                outputAnalysisDistributionNormalizer_                              , &
         &                                galacticFilter_                                                    , &
         &                                outputAnalysisCovarianceModelPoisson                               , &
         &                                covarianceBinomialBinsPerDecade                                    , &
         &                                covarianceBinomialMassHaloMinimum                                  , &
         &                                covarianceBinomialMassHaloMaximum                                    &
         &                               )
    ! Clean up.
    nullify(outputAnalysisPropertyExtractor_      )
    nullify(outputAnalysisPropertyOperator_       )
    nullify(outputAnalysisPropertyUnoperator_     )
    nullify(outputAnalysisWeightOperator_         )
    nullify(outputAnalysisWeightPropertyOperator_ )
    nullify(outputAnalysisWeightPropertyExtractor_)
    nullify(outputAnalysisDistributionOperator_   )
    nullify(outputAnalysisDistributionNormalizer_ )
    nullify(normalizerSequence                    )
    nullify(galacticFilterHaloIsolated_           )
    nullify(galacticFilterFormationTime_          )
    nullify(galacticFilter_                       )
    nullify(filters_                              )
    nullify(virialDensityContrast_                )
  return
  end function concentrationDistributionCDMCOCOConstructorInternal
