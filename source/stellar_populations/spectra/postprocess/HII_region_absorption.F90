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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a stellar population spectra postprocessor which removes ionizing light absorbed within HII regions.
  !!}

  use :: HII_Region_Escape_Fraction, only : hiiRegionEscapeFractionClass

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorHIIRegionAbsorption" docformat="rst">
   <description>
   A stellar population spectra postprocessor which removes the hydrogen-ionizing light absorbed within HII regions.
   Shortward of the Lyman limit the spectrum is multiplied by the escape fraction :math:`f_\mathrm{esc}(t)` of ionizing
   photons from HII regions of age :math:`t`, given by a :galacticus-class:`hiiRegionEscapeFractionClass` object;
   longward of it the spectrum is unchanged.

   The ionizing light which does not escape, a fraction :math:`1-f_\mathrm{esc}(t)`, is absorbed by the gas of the HII
   region and re-emitted as nebular emission. :galacticus-class:`nodePropertyExtractorLuminosityEmissionLine` weights
   its line luminosities by exactly that fraction. Applying this postprocessor to the stellar spectrum, with the same
   escape fraction, therefore avoids counting that energy twice---once in the stellar continuum and again in the
   lines---and in particular keeps it out of the light available to heat dust. Defining ``[hiiRegionEscapeFraction]``
   at the top level of the parameter file ensures that both use the same object.

   It is combined with other postprocessors, such as absorption by the intergalactic medium, through
   :galacticus-class:`stellarPopulationSpectraPostprocessorSequence`. It should not be combined with
   :galacticus-class:`stellarPopulationSpectraPostprocessorLycSuppress`, which removes all ionizing light whether or not
   it escapes.

   Its multiplier depends on age, but only as a factor applied at every age: no age is suppressed entirely. It
   therefore reports an unbounded age range and a sharp age window, which keeps it usable in the chains whose
   luminosities are differenced to isolate a range of ages---as dust attenuation by birth clouds requires---and that
   differencing remains exact, **provided that it appears in every chain being differenced** (both ``default`` and
   ``recent``, for example). If it is included in only one of them, light shortward of the Lyman limit is
   mis-apportioned between the age bins.

   Absorption of non-ionizing light by grains within HII regions is not included.
   </description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorHIIRegionAbsorption
     !!{RST
     A stellar population spectra postprocessor which removes ionizing light absorbed within HII regions.
     !!}
     private
     class(hiiRegionEscapeFractionClass), pointer :: hiiRegionEscapeFraction_ => null()
   contains
     final     ::                        hiiRegionAbsorptionDestructor
     procedure :: multiplier          => hiiRegionAbsorptionMultiplier
     procedure :: isRedshiftDependent => hiiRegionAbsorptionIsRedshiftDependent
     procedure :: ageWindowIsSharp    => hiiRegionAbsorptionAgeWindowIsSharp
  end type stellarPopulationSpectraPostprocessorHIIRegionAbsorption

  interface stellarPopulationSpectraPostprocessorHIIRegionAbsorption
     !!{RST
     Constructors for the :galacticus-class:`stellarPopulationSpectraPostprocessorHIIRegionAbsorption` stellar population
     spectra postprocessor class.
     !!}
     module procedure hiiRegionAbsorptionConstructorParameters
     module procedure hiiRegionAbsorptionConstructorInternal
  end interface stellarPopulationSpectraPostprocessorHIIRegionAbsorption

contains

  function hiiRegionAbsorptionConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`stellarPopulationSpectraPostprocessorHIIRegionAbsorption` stellar population
    spectra postprocessor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (stellarPopulationSpectraPostprocessorHIIRegionAbsorption)                :: self
    type (inputParameters                                         ), intent(inout) :: parameters
    class(hiiRegionEscapeFractionClass                            ), pointer       :: hiiRegionEscapeFraction_

    !![
    <objectBuilder class="hiiRegionEscapeFraction" name="hiiRegionEscapeFraction_" source="parameters"/>
    !!]
    self=stellarPopulationSpectraPostprocessorHIIRegionAbsorption(hiiRegionEscapeFraction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hiiRegionEscapeFraction_"/>
    !!]
    return
  end function hiiRegionAbsorptionConstructorParameters

  function hiiRegionAbsorptionConstructorInternal(hiiRegionEscapeFraction_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`stellarPopulationSpectraPostprocessorHIIRegionAbsorption` stellar
    population spectra postprocessor class.
    !!}
    implicit none
    type (stellarPopulationSpectraPostprocessorHIIRegionAbsorption)                        :: self
    class(hiiRegionEscapeFractionClass                            ), intent(in   ), target :: hiiRegionEscapeFraction_
    !![
    <constructorAssign variables="*hiiRegionEscapeFraction_"/>
    !!]

    return
  end function hiiRegionAbsorptionConstructorInternal

  subroutine hiiRegionAbsorptionDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`stellarPopulationSpectraPostprocessorHIIRegionAbsorption` stellar population
    spectra postprocessor class.
    !!}
    implicit none
    type(stellarPopulationSpectraPostprocessorHIIRegionAbsorption), intent(inout) :: self

    !![
    <objectDestructor name="self%hiiRegionEscapeFraction_"/>
    !!]
    return
  end subroutine hiiRegionAbsorptionDestructor

  double precision function hiiRegionAbsorptionMultiplier(self,wavelength,age,redshift) result(multiplier)
    !!{RST
    Return the fraction of the light at the given ``wavelength`` (in Å) from a population of the given ``age`` (in Gyr)
    which escapes its HII region: the escape fraction shortward of the Lyman limit, and unity longward of it.
    !!}
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen_atomic
    implicit none
    class           (stellarPopulationSpectraPostprocessorHIIRegionAbsorption), intent(inout) :: self
    double precision                                                          , intent(in   ) :: wavelength, age, &
         &                                                                                       redshift
    !$GLC attributes unused :: redshift

    if (wavelength < lymanSeriesLimitWavelengthHydrogen_atomic) then
       multiplier=self%hiiRegionEscapeFraction_%escapeFraction(age)
    else
       multiplier=1.0d0
    end if
    return
  end function hiiRegionAbsorptionMultiplier

  logical function hiiRegionAbsorptionIsRedshiftDependent(self) result(isRedshiftDependent)
    !!{RST
    Return false: absorption within HII regions does not depend on the redshift of the source.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorHIIRegionAbsorption), intent(inout) :: self
    !$GLC attributes unused :: self

    isRedshiftDependent=.false.
    return
  end function hiiRegionAbsorptionIsRedshiftDependent

  logical function hiiRegionAbsorptionAgeWindowIsSharp(self) result(ageWindowIsSharp)
    !!{RST
    Return true. The multiplier depends on age, but as a factor applied at every age rather than as a window which
    tapers: it suppresses no age entirely, and so differencing two chains which both include it isolates a range of
    ages exactly. Differencing is not exact if only one of the chains includes it---see the description of this class.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorHIIRegionAbsorption), intent(inout) :: self
    !$GLC attributes unused :: self

    ageWindowIsSharp=.true.
    return
  end function hiiRegionAbsorptionAgeWindowIsSharp
