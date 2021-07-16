!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a spin parameter distribution output analysis class.
!!}

  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass
  
  !![
  <outputAnalysis name="outputAnalysisSpinDistributionBett2007">
    <description>A stellar mass function output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisSpinDistribution) :: outputAnalysisSpinDistributionBett2007
     !!{
     A spinDistributionBett2007 output analysis class.
     !!}
     private
  end type outputAnalysisSpinDistributionBett2007

  interface outputAnalysisSpinDistributionBett2007
     !!{
     Constructors for the ``spinDistributionBett2007'' output analysis class.
     !!}
     module procedure spinDistributionBett2007ConstructorParameters
     module procedure spinDistributionBett2007ConstructorInternal
  end interface outputAnalysisSpinDistributionBett2007

contains

  function spinDistributionBett2007ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``spinDistributionBett2007'' output analysis class which takes a parameter set as input.
    !!}
    use :: Functions_Global, only : Virial_Density_Contrast_Percolation_Objects_Constructor_
    use :: Input_Parameters, only : inputParameter                                          , inputParameters
    implicit none
    type            (outputAnalysisSpinDistributionBett2007)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass                      ), pointer       :: outputTimes_
    class           (nbodyHaloMassErrorClass               ), pointer       :: nbodyHaloMassError_
    class           (haloMassFunctionClass                 ), pointer       :: haloMassFunction_
    class           (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass             ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterProfileScaleRadiusClass     ), pointer       :: darkMatterProfileScaleRadius_
    class           (*                                     ), pointer       :: percolationObjects_
    double precision                                                        :: logNormalRange
    logical                                                                 :: errorTolerant

    !![
    <inputParameter>
      <name>errorTolerant</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>Error tolerance for the N-body spin distribution operator.</description>
    </inputParameter>
    <inputParameter>
      <name>logNormalRange</name>
      <source>parameters</source>
      <defaultValue>100.0d0</defaultValue>
      <defaultSource>Approximately the range expected for the \cite{bett_spin_2007} ``QE'' cut.</defaultSource>
      <description>The multiplicative range of the log-normal distribution used to model the distribution of the mass and energy terms in the spin parameter. Specifically, the lognormal distribution is truncated outside the range $(\lambda_\mathrm{m}/R,\lambda_\mathrm{m} R$, where $\lambda_\mathrm{m}$ is the measured spin, and $R=${\normalfont \ttfamily [logNormalRange]}</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="outputTimes"                  name="outputTimes_"                  source="parameters"/>
    <objectBuilder class="nbodyHaloMassError"           name="nbodyHaloMassError_"           source="parameters"/>
    <objectBuilder class="haloMassFunction"             name="haloMassFunction_"             source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    percolationObjects_ => Virial_Density_Contrast_Percolation_Objects_Constructor_(parameters)
    self                =  outputAnalysisSpinDistributionBett2007(logNormalRange,errorTolerant,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfileDMO_,darkMatterProfileScaleRadius_,outputTimes_,percolationObjects_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="outputTimes_"                 />
    <objectDestructor name="nbodyHaloMassError_"          />
    <objectDestructor name="haloMassFunction_"            />
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function spinDistributionBett2007ConstructorParameters

  function spinDistributionBett2007ConstructorInternal(logNormalRange,errorTolerant,cosmologyFunctions_,nbodyHaloMassError_,haloMassFunction_,darkMatterHaloScale_,darkMatterProfileDMO_,darkMatterProfileScaleRadius_,outputTimes_,percolationObjects_) result(self)
    !!{
    Internal constructor for the ``spinDistributionBett2007'' output analysis class.
    !!}
    use :: Cosmology_Functions              , only : cosmologyFunctionsClass
    use :: Dark_Matter_Halo_Scales          , only : darkMatterHaloScaleClass
    use :: Galacticus_Paths                 , only : galacticusPath                  , pathTypeDataStatic
    use :: Halo_Mass_Functions              , only : haloMassFunctionClass
    use :: Output_Times                     , only : outputTimesClass
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass
    use :: Virial_Density_Contrast          , only : virialDensityContrastPercolation
    implicit none
    type            (outputAnalysisSpinDistributionBett2007)                              :: self
    double precision                                                     , intent(in   )  :: logNormalRange
    logical                                                              , intent(in   )  :: errorTolerant
    class           (cosmologyFunctionsClass               ), target     , intent(inout)  :: cosmologyFunctions_
    class           (outputTimesClass                      ), target     , intent(inout)  :: outputTimes_
    class           (nbodyHaloMassErrorClass               ), target     , intent(in   )  :: nbodyHaloMassError_
    class           (haloMassFunctionClass                 ), target     , intent(in   )  :: haloMassFunction_
    class           (darkMatterHaloScaleClass              ), target     , intent(in   )  :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass             ), target     , intent(in   )  :: darkMatterProfileDMO_
    class           (darkMatterProfileScaleRadiusClass     ), target     , intent(in   )  :: darkMatterProfileScaleRadius_
    class           (*                                     ), target     , intent(in   )  :: percolationObjects_
    type            (virialDensityContrastPercolation      ), pointer                     :: virialDensityContrast_

    ! Create a virial density contrast class to match the friends-of-friends halo definition used by Bett et al. (2007).
    allocate(virialDensityContrast_)
    !![
    <referenceConstruct object="virialDensityContrast_" constructor="virialDensityContrastPercolation(0.2d0,cosmologyFunctions_,percolationObjects_)"/>
    !!]
    self%outputAnalysisSpinDistribution=outputAnalysisSpinDistribution(                                                                                             &
         &                                                             char   (galacticusPath(pathTypeDataStatic)//'darkMatter/bett2007HaloSpinDistribution.hdf5'), &
         &                                                             var_str(                                    'spinDistributionBett2007'                    ), &
         &                                                             var_str(                                    'Distribution of halo spin parameters'        ), &
         &                                                             logNormalRange                                                                             , &
         &                                                             errorTolerant                                                                              , &
         &                                                             cosmologyFunctions_                                                                        , &
         &                                                             nbodyHaloMassError_                                                                        , &
         &                                                             haloMassFunction_                                                                          , &
         &                                                             darkMatterHaloScale_                                                                       , &
         &                                                             darkMatterProfileDMO_                                                                      , &
         &                                                             darkMatterProfileScaleRadius_                                                              , &
         &                                                             outputTimes_                                                                               , &
         &                                                             virialDensityContrast_                                                                       &
         &                                                            )
    !![
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function spinDistributionBett2007ConstructorInternal
