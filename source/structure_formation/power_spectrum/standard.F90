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
Implements a linear theory power spectrum class in which the power spectrum is just the transferred primordial power spectrum correctly normalized to :math:`z=0`.
!!}

  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredClass

  !![
  <powerSpectrum name="powerSpectrumStandard" docformat="rst">
   <description>
   Provides a linear theory power spectrum class in which the power spectrum is just the transferred primordial power spectrum (see :galacticus-class:`powerSpectrumPrimordialTransferred`), which is normalized by its :galacticus-class:`powerSpectrumNormalizer`.
   </description>
  </powerSpectrum>
  !!]
  type, extends(powerSpectrumClass) :: powerSpectrumStandard
     !!{RST
     A linear theory power spectrum class in which the power spectrum is just the transferred primordial power spectrum correctly normalized to :math:`z=0`.
     !!}
     private
     class(powerSpectrumPrimordialTransferredClass), pointer :: powerSpectrumPrimordialTransferred_ => null()
   contains
     final     ::                               standardDestructor
     procedure :: power                      => standardPower
     procedure :: powerLogarithmicDerivative => standardPowerLogarithmicDerivative
     procedure :: powerDimensionless         => standardPowerDimensionless
  end type powerSpectrumStandard

  interface powerSpectrumStandard
     !!{RST
     Constructors for the standard power spectrum class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface powerSpectrumStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the standard nonstandard power spectrum class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (powerSpectrumStandard                  )                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(powerSpectrumPrimordialTransferredClass), pointer       :: powerSpectrumPrimordialTransferred_

    !![
    <objectBuilder class="powerSpectrumPrimordialTransferred" name="powerSpectrumPrimordialTransferred_" source="parameters"/>
    !!]
    self=powerSpectrumStandard(powerSpectrumPrimordialTransferred_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="powerSpectrumPrimordialTransferred_"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(powerSpectrumPrimordialTransferred_) result(self)
    !!{RST
    Internal constructor for the standard power spectrum class.
    !!}
    implicit none
    type (powerSpectrumStandard)                                          :: self
    class(powerSpectrumPrimordialTransferredClass), intent(in   ), target :: powerSpectrumPrimordialTransferred_
    !![
    <constructorAssign variables="*powerSpectrumPrimordialTransferred_"/>
    !!]

    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{RST
    Destructor for the standard power spectrum class.
    !!}
    implicit none
    type(powerSpectrumStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%powerSpectrumPrimordialTransferred_"/>
    !!]
    return
  end subroutine standardDestructor

  double precision function standardPower(self,wavenumber,time)
    !!{RST
    Return the cosmological power spectrum for :math:`k=`\ ``wavenumber`` [Mpc\ :math:`^{-1}`].
    !!}
    implicit none
    class           (powerSpectrumStandard), intent(inout) :: self
    double precision                       , intent(in   ) :: wavenumber, time

    standardPower=self%powerSpectrumPrimordialTransferred_%power(wavenumber,time)
    return
  end function standardPower

  double precision function standardPowerLogarithmicDerivative(self,wavenumber,time)
    !!{RST
    Return the logarithmic derivative of the power spectrum, :math:`\mathrm{d}\ln P(k)/\mathrm{d}\ln k`, for :math:`k=`\ ``wavenumber`` [Mpc\ :math:`^{-1}`].
    !!}
    implicit none
    class           (powerSpectrumStandard), intent(inout) :: self
    double precision                       , intent(in   ) :: wavenumber, time

    standardPowerLogarithmicDerivative=self%powerSpectrumPrimordialTransferred_%logarithmicDerivative(wavenumber,time)
    return
  end function standardPowerLogarithmicDerivative

  double precision function standardPowerDimensionless(self,wavenumber,time)
    !!{RST
    Return the dimensionless power spectrum, :math:`\Delta^2(k)`, for :math:`k=`\ ``wavenumber`` [Mpc\ :math:`^{-1}`].
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumStandard), intent(inout) :: self
    double precision                       , intent(in   ) :: wavenumber, time

    standardPowerDimensionless=+4.0d0                          &
         &                     *Pi                             &
         &                     *           wavenumber      **3 &
         &                     *self%power(wavenumber,time)    &
         &                     /(                              &
         &                       +2.0d0                        &
         &                       *Pi                           &
         &                      )**3
    return
  end function standardPowerDimensionless
