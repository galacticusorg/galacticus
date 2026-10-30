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
  An identity power spectrum normalizer.
  !!}

  !![
  <powerSpectrumNormalizer name="powerSpectrumNormalizerIdentity" docformat="rst">
   <description>
   An identity power spectrum normalizer, which leaves the power spectrum unchanged. This is appropriate for power spectra which
   are already normalized, for example those read from file (see :galacticus-class:`powerSpectrumPrimordialTransferredFile`).
   </description>
  </powerSpectrumNormalizer>
  !!]
  type, extends(powerSpectrumNormalizerClass) :: powerSpectrumNormalizerIdentity
     !!{RST
     An identity power spectrum normalizer.
     !!}
     private
   contains
     procedure :: normalization => identityNormalization
  end type powerSpectrumNormalizerIdentity

  interface powerSpectrumNormalizerIdentity
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumNormalizerIdentity` power spectrum normalizer class.
     !!}
     module procedure identityConstructorParameters
  end interface powerSpectrumNormalizerIdentity

contains

  function identityConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumNormalizerIdentity` power spectrum normalizer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(powerSpectrumNormalizerIdentity)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=powerSpectrumNormalizerIdentity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identityConstructorParameters

  double precision function identityNormalization(self,powerUnnormalized) result(normalization)
    !!{RST
    Return a normalization factor of unity.
    !!}
    implicit none
    class    (powerSpectrumNormalizerIdentity         ), intent(inout), target :: self
    procedure(powerSpectrumNormalizerPowerUnnormalized)                        :: powerUnnormalized
    !$GLC attributes unused :: self, powerUnnormalized

    normalization=1.0d0
    return
  end function identityNormalization
