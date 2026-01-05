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

!!{
Implements a linear theory power spectrum class in which the power spectrum is just the transferred primordial power spectrum
correctly normalized to $z=0$.
!!}

  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceClass
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredClass

  !![
  <powerSpectrum name="powerSpectrumStandard">
   <description>Provides a linear theory power spectrum class in which the power spectrum is just the transferred primordial power spectrum correctly normalized to $z=0$.</description>
  </powerSpectrum>
  !!]
  type, extends(powerSpectrumClass) :: powerSpectrumStandard
     !!{
     A linear theory power spectrum class in which the power spectrum is just the transferred primordial power spectrum
     correctly normalized to $z=0$.
     !!}
     private
     class(cosmologicalMassVarianceClass          ), pointer :: cosmologicalMassVariance_           => null()
     class(powerSpectrumPrimordialTransferredClass), pointer :: powerSpectrumPrimordialTransferred_ => null()
   contains
     final     ::                               standardDestructor
     procedure :: power                      => standardPower
     procedure :: powerLogarithmicDerivative => standardPowerLogarithmicDerivative
     procedure :: powerDimensionless         => standardPowerDimensionless
     procedure :: descriptor                 => standardDescriptor
  end type powerSpectrumStandard

  interface powerSpectrumStandard
     !!{
     Constructors for the standard power spectrum class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface powerSpectrumStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the standard nonstandard power spectrum class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (powerSpectrumStandard                  )                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(cosmologicalMassVarianceClass          ), pointer       :: cosmologicalMassVariance_
    class(powerSpectrumPrimordialTransferredClass), pointer       :: powerSpectrumPrimordialTransferred_

    !![
    <objectBuilder class="cosmologicalMassVariance"           name="cosmologicalMassVariance_"           source="parameters"/>
    <objectBuilder class="powerSpectrumPrimordialTransferred" name="powerSpectrumPrimordialTransferred_" source="parameters"/>
    !!]
    self=powerSpectrumStandard(cosmologicalMassVariance_,powerSpectrumPrimordialTransferred_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologicalMassVariance_"          />
    <objectDestructor name="powerSpectrumPrimordialTransferred_"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(cosmologicalMassVariance_,powerSpectrumPrimordialTransferred_) result(self)
    !!{
    Internal constructor for the standard power spectrum class.
    !!}
    implicit none
    type (powerSpectrumStandard)                                          :: self
    class(cosmologicalMassVarianceClass          ), intent(in   ), target :: cosmologicalMassVariance_
    class(powerSpectrumPrimordialTransferredClass), intent(in   ), target :: powerSpectrumPrimordialTransferred_
    !![
    <constructorAssign variables="*cosmologicalMassVariance_, *powerSpectrumPrimordialTransferred_"/>
    !!]

    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the standard power spectrum class.
    !!}
    implicit none
    type(powerSpectrumStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologicalMassVariance_"          />
    <objectDestructor name="self%powerSpectrumPrimordialTransferred_"/>
    !!]
    return
  end subroutine standardDestructor

  double precision function standardPower(self,wavenumber,time)
    !!{
    Return the cosmological power spectrum for $k=${\normalfont \ttfamily wavenumber} [Mpc$^{-1}$].
    !!}
    implicit none
    class           (powerSpectrumStandard), intent(inout) :: self
    double precision                       , intent(in   ) :: wavenumber, time

    ! Compute the power spectrum.
    standardPower=+self%powerSpectrumPrimordialTransferred_%power             (wavenumber,time) &
         &        *self%cosmologicalMassVariance_          %powerNormalization(               )
    return
  end function standardPower

  double precision function standardPowerLogarithmicDerivative(self,wavenumber,time)
    !!{
    Return the logarithmic derivative of the power spectrum, $\mathrm{d}\ln P(k)/\mathrm{d}\ln k$, for $k=${\normalfont
    \ttfamily wavenumber} [Mpc$^{-1}$].
    !!}
    implicit none
    class           (powerSpectrumStandard), intent(inout) :: self
    double precision                       , intent(in   ) :: wavenumber, time

    standardPowerLogarithmicDerivative=self%powerSpectrumPrimordialTransferred_%logarithmicDerivative(wavenumber,time)
    return
  end function standardPowerLogarithmicDerivative

  double precision function standardPowerDimensionless(self,wavenumber,time)
    !!{
    Return the dimensionless power spectrum, $\Delta^2(k)$, for $k=${\normalfont \ttfamily wavenumber} [Mpc$^{-1}$].
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

  subroutine standardDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
      !!{
      Generate a descriptor for the standard power spectrum class.
      !!}
      use Input_Parameters, only : inputParameters
      implicit none
      class  (powerSpectrumStandard), intent(inout)           :: self
      type   (inputParameters      ), intent(inout)           :: descriptor
      logical                       , intent(in   ), optional :: includeClass, includeFileModificationTimes
      type   (inputParameters      )                          :: parameters
      !![
      <optionalArgument name="includeClass" defaultsTo=".true." />
      !!]
      
      if (includeClass_) call descriptor%addParameter('powerSpectrum','standard')
      parameters=descriptor%subparameters('powerSpectrum')
      if (associated(self%cosmologicalMassVariance_          )) call self%cosmologicalMassVariance_          %descriptorNormalizationOnly(parameters,includeClass,includeFileModificationTimes)
      if (associated(self%powerSpectrumPrimordialTransferred_)) call self%powerSpectrumPrimordialTransferred_%descriptor                 (parameters,includeClass,includeFileModificationTimes)
      return      
    end subroutine standardDescriptor
