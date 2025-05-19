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
  Implements a concentration vs. halo mass analysis class matched to the
  \cite{ludlow_mass-concentration-redshift_2016} CDM sample.
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <outputAnalysis name="outputAnalysisConcentrationVsHaloMassCDMLudlow2016">
   <description>A concentration vs. halo mass analysis class matched to the \cite{ludlow_mass-concentration-redshift_2016} CDM sample.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMeanFunction1D) :: outputAnalysisConcentrationVsHaloMassCDMLudlow2016
     !!{
     A concentration vs. halo mass analysis class matched to the \cite{ludlow_mass-concentration-redshift_2016} CDM sample.
     !!}
     private
    class(cosmologyParametersClass  ), pointer :: cosmologyParameters_   => null()
    class(cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
    class(virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
    class(nbodyHaloMassErrorClass   ), pointer :: nbodyHaloMassError_    => null()
    class(darkMatterProfileDMOClass ), pointer :: darkMatterProfileDMO_  => null()
  contains
    final :: concentrationVsHaloMassCDMLudlow2016Destructor
 end type outputAnalysisConcentrationVsHaloMassCDMLudlow2016

  interface outputAnalysisConcentrationVsHaloMassCDMLudlow2016
     !!{
     Constructors for the {\normalfont \ttfamily concentrationVsHaloMassCDMLudlow2016} output analysis class.
     !!}
     module procedure concentrationVsHaloMassCDMLudlow2016ConstructorParameters
     module procedure concentrationVsHaloMassCDMLudlow2016ConstructorInternal
  end interface outputAnalysisConcentrationVsHaloMassCDMLudlow2016

contains

  function concentrationVsHaloMassCDMLudlow2016ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily concentrationVsHaloMassCDMLudlow2016} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions , only : cosmologyFunctions , cosmologyFunctionsClass
    use :: Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass
    use :: Input_Parameters    , only : inputParameter     , inputParameters
    implicit none
    type (outputAnalysisConcentrationVsHaloMassCDMLudlow2016)                :: self
    type (inputParameters                                   ), intent(inout) :: parameters
    class(cosmologyParametersClass                          ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass                           ), pointer       :: cosmologyFunctions_
    class(virialDensityContrastClass                        ), pointer       :: virialDensityContrast_
    class(outputTimesClass                                  ), pointer       :: outputTimes_
    class(nbodyHaloMassErrorClass                           ), pointer       :: nbodyHaloMassError_
    class(darkMatterProfileDMOClass                         ), pointer       :: darkMatterProfileDMO_

    !![
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="outputTimes"           name="outputTimes_"           source="parameters"/>
    <objectBuilder class="nbodyHaloMassError"    name="nbodyHaloMassError_"    source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=outputAnalysisConcentrationVsHaloMassCDMLudlow2016(darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,nbodyHaloMassError_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    nullify(cosmologyParameters_)
    nullify(cosmologyFunctions_ )
    nullify(nbodyHaloMassError_ )
    nullify(outputTimes_        )
    !![
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="outputTimes_"          />
    <objectDestructor name="nbodyHaloMassError_"   />
    <objectDestructor name="darkMatterProfileDMO_" />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function concentrationVsHaloMassCDMLudlow2016ConstructorParameters

  function concentrationVsHaloMassCDMLudlow2016ConstructorInternal(darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,nbodyHaloMassError_,outputTimes_) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily concentrationVsHaloMassCDMLudlow2016} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters                  , only : cosmologyParametersClass
    use :: Galactic_Filters                      , only : filterList                                        , galacticFilterAll                  , galacticFilterBasicMass, galacticFilterHaloIsolated
    use :: Error                                 , only : Error_Report
    use :: Input_Paths                           , only : inputPath                                         , pathTypeDataStatic
    use :: HDF5_Access                           , only : hdf5Access
    use :: IO_HDF5                               , only : hdf5Object
    use :: ISO_Varying_String                    , only : var_str
    use :: Node_Property_Extractors              , only : nodePropertyExtractorConcentration                , nodePropertyExtractorMassHalo
    use :: Numerical_Comparison                  , only : Values_Agree
    use :: Numerical_Constants_Astronomical      , only : massSolar
    use :: Output_Analyses_Options               , only : outputAnalysisCovarianceModelBinomial
    use :: Output_Analysis_Distribution_Operators, only : outputAnalysisDistributionOperatorRndmErrNbodyMass
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorAntiLog10           , outputAnalysisPropertyOperatorLog10
    use :: Output_Analysis_Weight_Operators      , only : outputAnalysisWeightOperatorIdentity
    use :: Output_Times                          , only : outputTimesClass
    use :: Statistics_NBody_Halo_Mass_Errors     , only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast               , only : fixedDensityTypeCritical                          , virialDensityContrastFixed
    implicit none
    type            (outputAnalysisConcentrationVsHaloMassCDMLudlow2016)                                :: self
    class           (cosmologyParametersClass                          ), target       , intent(in   )  :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), target       , intent(in   )  :: cosmologyFunctions_
    class           (virialDensityContrastClass                        ), target       , intent(in   )  :: virialDensityContrast_
    class           (nbodyHaloMassErrorClass                           ), target       , intent(in   )  :: nbodyHaloMassError_
    class           (outputTimesClass                                  ), target       , intent(inout)  :: outputTimes_
    class           (darkMatterProfileDMOClass                         ), target       , intent(inout)  :: darkMatterProfileDMO_
    integer         (c_size_t                                          ), parameter                     :: massHaloCount                         =26
    double precision                                                    , parameter                     :: massHaloMinimum                       = 1.0d10, massHaloMaximum                    =1.0d15
    integer                                                             , parameter                     :: covarianceBinomialBinsPerDecade       =10
    double precision                                                    , parameter                     :: covarianceBinomialMassHaloMinimum     = 1.0d08, covarianceBinomialMassHaloMaximum  =1.0d16
    double precision                                                    , allocatable  , dimension(:  ) :: massHaloLogarithmic
    double precision                                                    , allocatable  , dimension(:,:) :: outputWeight
    type            (galacticFilterBasicMass                           ), pointer                       :: galacticFilterBasicMass_
    type            (galacticFilterHaloIsolated                        ), pointer                       :: galacticFilterHaloIsolated_
    type            (galacticFilterAll                                 ), pointer                       :: galacticFilterAll_
    type            (filterList                                        ), pointer                       :: filters_
    type            (outputAnalysisDistributionOperatorRndmErrNbodyMass), pointer                       :: outputAnalysisDistributionOperator_
    type            (outputAnalysisWeightOperatorIdentity              ), pointer                       :: outputAnalysisWeightOperator_
    type            (outputAnalysisPropertyOperatorLog10               ), pointer                       :: outputAnalysisPropertyOperator_              , outputAnalysisWeightPropertyOperator_
    type            (outputAnalysisPropertyOperatorAntiLog10           ), pointer                       :: outputAnalysisPropertyUnoperator_
    type            (nodePropertyExtractorMassHalo                     ), pointer                       :: nodePropertyExtractor_
    type            (nodePropertyExtractorConcentration                ), pointer                       :: outputAnalysisWeightPropertyExtractor_
    type            (virialDensityContrastFixed                        ), pointer                       :: virialDensityContrastDefinition_
    integer         (c_size_t                                          )                                :: iOutput
    type            (hdf5Object                                        )                                :: dataFile
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *virialDensityContrast_, *nbodyHaloMassError_, *outputTimes_"/>
    !!]
    
    ! Construct mass bins matched to those used by Ludlow et al. (2016).
    !$ call hdf5Access%set()
    dataFile=hdf5Object(char(inputPath(pathTypeDataStatic)//'darkMatter/concentrationMassRelationCDMLudlow2016.hdf5'),readOnly=.true.)
    call dataFile%readDataset('massHalo',massHaloLogarithmic)
    !$ call hdf5Access%unset()
    massHaloLogarithmic=log10(massHaloLogarithmic)
    ! Compute weights that apply to each output redshift.
    allocate(outputWeight(size(massHaloLogarithmic),outputTimes_%count()))
    outputWeight=0.0d0
    do iOutput=1,outputTimes_%count()
       if (Values_Agree(outputTimes_%redshift(iOutput),0.0d0,absTol=1.0d-10)) outputWeight(:,iOutput)=1.0d0
    end do
    if (any(sum(outputWeight,dim=2) /= 1.0d0)) call Error_Report('zero redshift output is required'//{introspection:location})
    ! Build a filter which select isolated halos with M200c mass above some coarse lower limit suitable for this sample.
    allocate(galacticFilterBasicMass_        )
    allocate(galacticFilterHaloIsolated_     )
    allocate(galacticFilterAll_              )
    allocate(filters_                        )
    allocate(filters_                   %next)
    filters_                        %filter_ => galacticFilterHaloIsolated_
    filters_                   %next%filter_ => galacticFilterBasicMass_
    galacticFilterBasicMass_                 =  galacticFilterBasicMass    (massThreshold=1.0d7   )
    galacticFilterHaloIsolated_              =  galacticFilterHaloIsolated (                      )
    galacticFilterAll_                       =  galacticFilterAll          (              filters_)
    ! Build N-body mass error distribution operator.
    allocate(outputAnalysisDistributionOperator_    )
    outputAnalysisDistributionOperator_    =  outputAnalysisDistributionOperatorRndmErrNbodyMass(nbodyHaloMassError_                                                                                                           )
    ! Build identity weight operator.
    allocate(outputAnalysisWeightOperator_          )
    outputAnalysisWeightOperator_          =  outputAnalysisWeightOperatorIdentity              (                                                                                                                              )
    ! Build log10() property operator.
    allocate(outputAnalysisPropertyOperator_        )
    outputAnalysisPropertyOperator_        =  outputAnalysisPropertyOperatorLog10               (                                                                                                                              )
    ! Build a log10 weight property operators.
    allocate(outputAnalysisWeightPropertyOperator_  )
    outputAnalysisWeightPropertyOperator_  =  outputAnalysisPropertyOperatorLog10               (                                                                                                                              )
    ! Build anti-log10() property operator.
    allocate(outputAnalysisPropertyUnoperator_      )
    outputAnalysisPropertyUnoperator_      =  outputAnalysisPropertyOperatorAntiLog10           (                                                                                                                              )
    ! Create a virial density contrast object matched to the definition used by Ludlow et al. (2016).
    allocate(virialDensityContrastDefinition_       )
    virialDensityContrastDefinition_       =  virialDensityContrastFixed                        (200.0d0,fixedDensityTypeCritical,2.0d0,cosmologyParameters_,cosmologyFunctions_                                               )
    ! Create a concentration weight property extractor.
    allocate(outputAnalysisWeightPropertyExtractor_ )
    outputAnalysisWeightPropertyExtractor_ =  nodePropertyExtractorConcentration                (.false.,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_)
    ! Create a halo mass property extractor.
    allocate(nodePropertyExtractor_       )
    nodePropertyExtractor_                 =  nodePropertyExtractorMassHalo                     (.false.,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_)
    ! Build the object.
    self%outputAnalysisMeanFunction1D=outputAnalysisMeanFunction1D(                                                       &
         &                                                         var_str('concentrationHaloMassRelationCDMLudlow2016'), &
         &                                                         var_str('Concentration vs. halo mass relation'      ), &
         &                                                         var_str('massHalo'                                  ), &
         &                                                         var_str('Halo mass [M₂₀₀c]'                         ), &
         &                                                         var_str('M☉'                                        ), &
         &                                                         massSolar                                            , &
         &                                                         var_str('log10Concentration'                        ), &
         &                                                         var_str('log₁₀ of Concentration [log₁₀(r₂₀₀c/rs)]'  ), &
         &                                                         var_str(' '                                         ), &
         &                                                         0.0d0                                                , &
         &                                                         massHaloLogarithmic                                  , &
         &                                                         5_c_size_t                                           , &
         &                                                         outputWeight                                         , &
         &                                                         nodePropertyExtractor_                               , &
         &                                                         outputAnalysisWeightPropertyExtractor_               , &
         &                                                         outputAnalysisPropertyOperator_                      , &
         &                                                         outputAnalysisWeightPropertyOperator_                , &
         &                                                         outputAnalysisPropertyUnoperator_                    , &
         &                                                         outputAnalysisWeightOperator_                        , &
         &                                                         outputAnalysisDistributionOperator_                  , &
         &                                                         galacticFilterAll_                                   , &
         &                                                         outputTimes_                                         , &
         &                                                         outputAnalysisCovarianceModelBinomial                , &
         &                                                         covarianceBinomialBinsPerDecade                      , &
         &                                                         covarianceBinomialMassHaloMinimum                    , &
         &                                                         covarianceBinomialMassHaloMaximum                      &
         &                                                        )
    ! Clean up.
    nullify(galacticFilterAll_                    )
    nullify(galacticFilterBasicMass_              )
    nullify(galacticFilterHaloIsolated_           )
    nullify(filters_                              )
    nullify(outputAnalysisDistributionOperator_   )
    nullify(outputAnalysisWeightOperator_         )
    nullify(outputAnalysisPropertyOperator_       )
    nullify(outputAnalysisPropertyUnoperator_     )
    nullify(outputAnalysisWeightPropertyOperator_ )
    nullify(outputAnalysisWeightPropertyExtractor_)
    nullify(nodePropertyExtractor_                )
    nullify(virialDensityContrastDefinition_      )
    return
  end function concentrationVsHaloMassCDMLudlow2016ConstructorInternal

  subroutine concentrationVsHaloMassCDMLudlow2016Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily concentrationVsHaloMassCDMLudlow2016} output analysis class.
    !!}
    implicit none
    type(outputAnalysisConcentrationVsHaloMassCDMLudlow2016), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%outputTimes_"          />
    <objectDestructor name="self%darkMatterProfileDMO_" />
    <objectDestructor name="self%nbodyHaloMassError_"   />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine concentrationVsHaloMassCDMLudlow2016Destructor
