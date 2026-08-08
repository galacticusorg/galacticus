!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{RST
  An implementation of dark matter halo splashback radii using the model of :cite:t:`diemer_splashback_2020`.
  !!}

  !+ Implemented by Andrew Benson with assistance from Claude.

  !![
  <darkMatterHaloSplashbackRadius name="darkMatterHaloSplashbackRadiusDiemer2020" docformat="rst">
   <description>
   A splashback radius class implementing the model of :cite:t:`diemer_splashback_2020`. The splashback radius and mass are given by

   .. math::

      {X_\mathrm{sp} \over X_\mathrm{200m}} = A + B \exp\left( - {\Gamma \over C} \right),

   where :math:`X` is either radius or mass, and the coefficients :math:`A`, :math:`B`, and :math:`C` depend on the matter
   density parameter, :math:`\Omega_\mathrm{M}`, the peak height, :math:`\nu`, computed from :math:`M_\mathrm{200m}`,
   and on which statistic of the distribution of particle apocenters is being predicted (set by ``[definition]``).
   Here :math:`\Gamma` is the mass accretion rate measured over one dynamical time. This model supersedes
   that of :cite:t:`diemer_splashback_2017`, having been recalibrated against a larger set of simulations.

   The model was calibrated over the ranges :math:`0 \le \nu \le 5` and :math:`0 \le z \le 8`; it is evaluated here
   without restriction, so results outside of those ranges should be treated as extrapolations.
   </description>
  </darkMatterHaloSplashbackRadius>
  !!]
  type, extends(darkMatterHaloSplashbackRadiusDiemer) :: darkMatterHaloSplashbackRadiusDiemer2020
     !!{RST
     A dark matter halo splashback radius class implementing the model of :cite:t:`diemer_splashback_2020`.
     !!}
     private
   contains
     final :: diemer2020Destructor
  end type darkMatterHaloSplashbackRadiusDiemer2020

  interface darkMatterHaloSplashbackRadiusDiemer2020
     !!{RST
     Constructors for the :galacticus-class:`darkMatterHaloSplashbackRadiusDiemer2020` splashback radius class.
     !!}
     module procedure diemer2020ConstructorParameters
     module procedure diemer2020ConstructorInternal
  end interface darkMatterHaloSplashbackRadiusDiemer2020

contains

  function diemer2020ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterHaloSplashbackRadiusDiemer2020` splashback radius class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterHaloSplashbackRadiusDiemer2020)                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass                ), pointer       :: cosmologyParameters_
    class(criticalOverdensityClass                ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass           ), pointer       :: cosmologicalMassVariance_
    class(virialDensityContrastClass              ), pointer       :: virialDensityContrast_
    class(darkMatterHaloMassAccretionHistoryClass ), pointer       :: darkMatterHaloMassAccretionHistory_
    type (varying_string                          )                :: definition

    !![
    <inputParameter docformat="rst">
      <name>definition</name>
      <defaultValue>var_str('mean')</defaultValue>
      <description>
      The statistic of the distribution of particle apocenters which defines the splashback radius. Allowed values are
      ``mean`` (the mean of the distribution), and ``percentile50``, ``percentile75``, ``percentile87``, and
      ``percentile90`` (the corresponding percentiles of the distribution).
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"/>
    <objectBuilder class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"/>
    <objectBuilder class="criticalOverdensity"                name="criticalOverdensity_"                source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"           name="cosmologicalMassVariance_"           source="parameters"/>
    <objectBuilder class="virialDensityContrast"              name="virialDensityContrast_"              source="parameters"/>
    <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"/>
    !!]
    self=darkMatterHaloSplashbackRadiusDiemer2020(enumerationSplashbackDefinitionEncode(char(definition),includesPrefix=.false.),cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_,darkMatterHaloMassAccretionHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="cosmologyParameters_"               />
    <objectDestructor name="criticalOverdensity_"               />
    <objectDestructor name="cosmologicalMassVariance_"          />
    <objectDestructor name="virialDensityContrast_"             />
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    !!]
    return
  end function diemer2020ConstructorParameters

  function diemer2020ConstructorInternal(definition,cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_,darkMatterHaloMassAccretionHistory_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`darkMatterHaloSplashbackRadiusDiemer2020` splashback radius class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (darkMatterHaloSplashbackRadiusDiemer2020)                        :: self
    type (enumerationSplashbackDefinitionType     ), intent(in   )         :: definition
    class(cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologyParametersClass                ), intent(in   ), target :: cosmologyParameters_
    class(criticalOverdensityClass                ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass           ), intent(in   ), target :: cosmologicalMassVariance_
    class(virialDensityContrastClass              ), intent(in   ), target :: virialDensityContrast_
    class(darkMatterHaloMassAccretionHistoryClass ), intent(in   ), target :: darkMatterHaloMassAccretionHistory_
    !![
    <constructorAssign variables="definition, *cosmologyFunctions_, *cosmologyParameters_, *criticalOverdensity_, *cosmologicalMassVariance_, *virialDensityContrast_, *darkMatterHaloMassAccretionHistory_"/>
    !!]

    if (.not.enumerationSplashbackDefinitionIsValid(definition)) call Error_Report('invalid splashback definition'//{introspection:location})
    ! Select the coefficients of the fitting functions. Separate calibrations are provided for the mean of the distribution
    ! of particle apocenters, and for percentiles of that distribution.
    if (definition == splashbackDefinitionMean) then
       self%parametersRadius       =diemer2020ParametersRadiusMean
       self%parametersMass         =diemer2020ParametersMassMean
       self%parametersRadiusScatter=diemer2020ScatterRadiusMean
       self%parametersMassScatter  =diemer2020ScatterMassMean
    else
       self%parametersRadius       =diemer2020ParametersRadiusPercentile
       self%parametersMass         =diemer2020ParametersMassPercentile
       self%parametersRadiusScatter=diemer2020ScatterRadiusPercentile
       self%parametersMassScatter  =diemer2020ScatterMassPercentile
    end if
    return
  end function diemer2020ConstructorInternal

  subroutine diemer2020Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`darkMatterHaloSplashbackRadiusDiemer2020` splashback radius class.
    !!}
    implicit none
    type(darkMatterHaloSplashbackRadiusDiemer2020), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%cosmologyParameters_"               />
    <objectDestructor name="self%criticalOverdensity_"               />
    <objectDestructor name="self%cosmologicalMassVariance_"          />
    <objectDestructor name="self%virialDensityContrast_"             />
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    !!]
    return
  end subroutine diemer2020Destructor
