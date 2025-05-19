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
  Implements a concentration distribution output analysis class for COCO CDM data.
  !!}

  !![
  <outputAnalysis name="outputAnalysisConcentrationDistributionCDMCOCO">
    <description>A concentration distribution function output analysis class for COCO CDM data.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisConcentrationDistribution) :: outputAnalysisConcentrationDistributionCDMCOCO
     !!{
     A concentration distribution output analysis class for COCO CDM data.
     !!}
     private
     integer :: distributionNumber
   contains
  end type outputAnalysisConcentrationDistributionCDMCOCO

  interface outputAnalysisConcentrationDistributionCDMCOCO
     !!{
     Constructors for the {\normalfont \ttfamily concentrationDistributionCDMCOCO} output analysis class.
     !!}
     module procedure concentrationDistributionCDMCOCOConstructorParameters
     module procedure concentrationDistributionCDMCOCOConstructorInternal
  end interface outputAnalysisConcentrationDistributionCDMCOCO

contains

  function concentrationDistributionCDMCOCOConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily concentrationDistributionCDMCOCO} output analysis class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctions , cosmologyFunctionsClass
    use :: Cosmology_Parameters             , only : cosmologyParameters, cosmologyParametersClass
    use :: Input_Parameters                 , only : inputParameter     , inputParameters
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassError , nbodyHaloMassErrorClass
    implicit none
    type            (outputAnalysisConcentrationDistributionCDMCOCO)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (cosmologyParametersClass                      ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass                              ), pointer       :: outputTimes_
    class           (nbodyHaloMassErrorClass                       ), pointer       :: nbodyHaloMassError_
    class           (virialDensityContrastClass                    ), pointer       :: virialDensityContrast_
    class           (darkMatterProfileDMOClass                     ), pointer       :: darkMatterProfileDMO_
    integer                                                                         :: distributionNumber
    double precision                                                                :: rootVarianceFractionalMinimum
    
    !![
    <inputParameter>
      <name>rootVarianceFractionalMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum fractional root variance (relative to the target dataset).</description>
    </inputParameter>
    <inputParameter>
      <name>distributionNumber</name>
      <source>parameters</source>
      <description>The number (1-7) of the distribution to compute.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="outputTimes"           name="outputTimes_"           source="parameters"/>
    <objectBuilder class="nbodyHaloMassError"    name="nbodyHaloMassError_"    source="parameters"/>
    !!]
    self=outputAnalysisConcentrationDistributionCDMCOCO(distributionNumber,rootVarianceFractionalMinimum,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,outputTimes_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters" />
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="outputTimes_"          />
    <objectDestructor name="nbodyHaloMassError_"   />
    <objectDestructor name="virialDensityContrast_"/>
    <objectDestructor name="darkMatterProfileDMO_" />
    !!]
    return
  end function concentrationDistributionCDMCOCOConstructorParameters

  function concentrationDistributionCDMCOCOConstructorInternal(distributionNumber,rootVarianceFractionalMinimum,darkMatterProfileDMO_,cosmologyParameters_,cosmologyFunctions_,nbodyHaloMassError_,outputTimes_,virialDensityContrast_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily concentrationDistributionCDMCOCO} output analysis class.
    !!}
    use :: Error                            , only : Error_Report
    use :: Cosmology_Functions              , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters             , only : cosmologyParametersClass
    use :: Input_Paths                      , only : inputPath                , pathTypeDataStatic
    use :: Output_Times                     , only : outputTimesClass
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass
    use :: String_Handling                  , only : operator(//)
    use :: Virial_Density_Contrast          , only : fixedDensityTypeCritical , virialDensityContrastClass, virialDensityContrastFixed
    implicit none
    type            (outputAnalysisConcentrationDistributionCDMCOCO)                           :: self
    class           (cosmologyParametersClass                      ), target   , intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass                       ), target   , intent(inout) :: cosmologyFunctions_
    class           (virialDensityContrastClass                    ), target   , intent(in   ) :: virialDensityContrast_
    class           (outputTimesClass                              ), target   , intent(inout) :: outputTimes_
    class           (nbodyHaloMassErrorClass                       ), target   , intent(in   ) :: nbodyHaloMassError_
    class           (darkMatterProfileDMOClass                     ), target   , intent(in   ) :: darkMatterProfileDMO_
    integer                                                                    , intent(in   ) :: distributionNumber
    double precision                                                           , intent(in   ) :: rootVarianceFractionalMinimum
    type            (virialDensityContrastFixed                    ), pointer                  :: virialDensityContrastDefinition_
    double precision                                                , parameter                :: haloDensityContrast             =200.0d0
    !![
    <constructorAssign variables="distributionNumber"/>
    !!]
    
    ! Validate input.
    if (distributionNumber < 1 .or. distributionNumber > 7) call Error_Report('distributionNumber ∈ [1..7] is required'//{introspection:location})
    ! Create a virial density contrast object matched to the definition used by Ludlow et al. (2016).
    allocate(virialDensityContrastDefinition_                )
    !![
    <referenceConstruct object="virialDensityContrastDefinition_">
     <constructor>
      virialDensityContrastFixed(                          &amp;
         &amp;                   haloDensityContrast     , &amp;
         &amp;                   fixedDensityTypeCritical, &amp;
         &amp;                   2.0d0                   , &amp;
         &amp;                   cosmologyParameters_    , &amp;
         &amp;                   cosmologyFunctions_       &amp;
         &amp;                  )
     </constructor>
    </referenceConstruct>
    !!]
    self%outputAnalysisConcentrationDistribution=outputAnalysisConcentrationDistribution(                                                                                                                    &
         &                                                                               char   (inputPath(pathTypeDataStatic)//'darkMatter/concentrationDistributionCocoCDM'//distributionNumber//'.hdf5'), &
         &                                                                               var_str(                                    'concentrationDistributionCDMCOCO'                                   ), &
         &                                                                               var_str(                                    'Distribution of halo concentrations'                                ), &
         &                                                                               rootVarianceFractionalMinimum                                                                                     , &
         &                                                                               darkMatterProfileDMO_                                                                                             , &
         &                                                                               cosmologyParameters_                                                                                              , &
         &                                                                               cosmologyFunctions_                                                                                               , &
         &                                                                               nbodyHaloMassError_                                                                                               , &
         &                                                                               virialDensityContrast_                                                                                            , &
         &                                                                               virialDensityContrastDefinition_                                                                                  , &
         &                                                                               outputTimes_                                                                                                        &
         &                                                                              )
    !![
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function concentrationDistributionCDMCOCOConstructorInternal
