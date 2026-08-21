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
  Implements a power-law dust extinction curve.
  !!}

  !![
  <dustExtinctionCurve name="dustExtinctionCurvePowerLaw" docformat="rst">
   <description>
   A power-law dust extinction curve, :math:`k(\lambda)/k_\mathrm{V} = (\lambda/\lambda_\mathrm{V})^{-\alpha}`, where
   :math:`\lambda_\mathrm{V}` is the :math:`V`-band effective wavelength. The default exponent, :math:`\alpha=0.7`, is
   that adopted by :cite:t:`charlot_simple_2000` for both the diffuse interstellar medium and stellar birth clouds.

   The curve is defined at all wavelengths, and diverges as :math:`\lambda \rightarrow 0`; it should not be used far
   into the ultraviolet without consideration of whether a power law remains appropriate there.
   </description>
  </dustExtinctionCurve>
  !!]
  type, extends(dustExtinctionCurveClass) :: dustExtinctionCurvePowerLaw
     !!{RST
     A power-law dust extinction curve.
     !!}
     private
     double precision :: exponent_, wavelengthReference
   contains
     procedure :: attenuationRelative => powerLawAttenuationRelative
  end type dustExtinctionCurvePowerLaw

  interface dustExtinctionCurvePowerLaw
     !!{RST
     Constructors for the :galacticus-class:`dustExtinctionCurvePowerLaw` dust extinction curve class.
     !!}
     module procedure powerLawConstructorParameters
     module procedure powerLawConstructorInternal
  end interface dustExtinctionCurvePowerLaw

contains

  function powerLawConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustExtinctionCurvePowerLaw` dust extinction curve class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustExtinctionCurvePowerLaw)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: exponent_ , wavelengthReference

    !![
    <inputParameter docformat="rst">
      <name>exponent</name>
      <variable>exponent_</variable>
      <defaultValue>0.7d0</defaultValue>
      <defaultSource>:cite:t:`charlot_simple_2000`</defaultSource>
      <description>
      The exponent :math:`\alpha` of the power-law dust extinction curve,
      :math:`k(\lambda)/k_\mathrm{V}=(\lambda/\lambda_\mathrm{V})^{-\alpha}`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavelengthReference</name>
      <defaultValue>wavelengthVBand</defaultValue>
      <description>
      The wavelength, in Å, at which the power law is normalized to unity. The default is the effective wavelength of
      the Buser :math:`V` filter, which is the definition used throughout this framework. Set it to
      :math:`5500\,`Å to reproduce results from the ``lmnstyStllrCF2000`` property extractor, which adopted that
      round value instead.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustExtinctionCurvePowerLaw(exponent_,wavelengthReference)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawConstructorParameters

  function powerLawConstructorInternal(exponent_,wavelengthReference) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustExtinctionCurvePowerLaw` dust extinction curve class.
    !!}
    implicit none
    type            (dustExtinctionCurvePowerLaw)                :: self
    double precision                             , intent(in   ) :: exponent_, wavelengthReference
    !![
    <constructorAssign variables="exponent_, wavelengthReference"/>
    !!]

    return
  end function powerLawConstructorInternal

  double precision function powerLawAttenuationRelative(self,wavelength) result(attenuationRelative)
    !!{RST
    Return the relative dust opacity for a power-law extinction curve.
    !!}
    implicit none
    class           (dustExtinctionCurvePowerLaw), intent(inout) :: self
    double precision                             , intent(in   ) :: wavelength

    attenuationRelative=+(                          &
         &                +     wavelength          &
         &                /self%wavelengthReference &
         &               )**(-self%exponent_)
    return
  end function powerLawAttenuationRelative
