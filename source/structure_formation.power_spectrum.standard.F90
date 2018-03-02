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

!% Contains a module which implements a linear theory power spectrum class in which the power spectrum is just the transferred
!% primordial power spectrum correctly normalized to $z=0$.

  use Cosmological_Density_Field
  use Power_Spectra_Primordial_Transferred
  
  !# <powerSpectrum name="powerSpectrumStandard">
  !#  <description>Provides a linear theory power spectrum class in which the power spectrum is just the transferred primordial power spectrum correctly normalized to $z=0$.</description>
  !# </powerSpectrum>
  type, extends(powerSpectrumClass) :: powerSpectrumStandard
     !% A linear theory power spectrum class in which the power spectrum is just the transferred primordial power spectrum
     !% correctly normalized to $z=0$.
     private
     class(cosmologicalMassVarianceClass          ), pointer :: cosmologicalMassVariance_
     class(powerSpectrumPrimordialTransferredClass), pointer :: powerSpectrumPrimordialTransferred_
   contains
     final     ::                               standardDestructor
     procedure :: power                      => standardPower
     procedure :: powerLogarithmicDerivative => standardPowerLogarithmicDerivative
     procedure :: powerDimensionless         => standardPowerDimensionless
  end type powerSpectrumStandard

  interface powerSpectrumStandard
     !% Constructors for the standard power spectrum class.
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface powerSpectrumStandard

contains

  function standardConstructorParameters(parameters)
    !% Constructor for the standard nonstandard power spectrum class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(powerSpectrumStandard)                :: standardConstructorParameters
    type(inputParameters      ), intent(inout) :: parameters

    ! Check and read parameters.
    ! Build objects.
    !# <objectBuilder class="cosmologicalMassVariance"           name="standardConstructorParameters%cosmologicalMassVariance_"           source="parameters"/>
    !# <objectBuilder class="powerSpectrumPrimordialTransferred" name="standardConstructorParameters%powerSpectrumPrimordialTransferred_" source="parameters"/>
    !# <inputParametersValidate source="parameters"/>
    return
  end function standardConstructorParameters

  function standardConstructorInternal(cosmologicalMassVariance_,powerSpectrumPrimordialTransferred_)
    !% Internal constructor for the standard nonstandard power spectrum class.
    implicit none
    type (powerSpectrumStandard)                                          :: standardConstructorInternal
    class(cosmologicalMassVarianceClass          ), intent(in   ), target :: cosmologicalMassVariance_
    class(powerSpectrumPrimordialTransferredClass), intent(in   ), target :: powerSpectrumPrimordialTransferred_

    standardConstructorInternal%cosmologicalMassVariance_           => cosmologicalMassVariance_
    standardConstructorInternal%powerSpectrumPrimordialTransferred_ => powerSpectrumPrimordialTransferred_
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !% Destructor for the standard nonstandard power spectrum class.
    implicit none
    type(powerSpectrumStandard), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologicalMassVariance_"          />
    !# <objectDestructor name="self%powerSpectrumPrimordialTransferred_"/>
    return
  end subroutine standardDestructor

  double precision function standardPower(self,wavenumber)
    !% Return the cosmological power spectrum for $k=${\normalfont \ttfamily wavenumber} [Mpc$^{-1}$].
    implicit none
    class           (powerSpectrumStandard), intent(inout) :: self
    double precision                       , intent(in   ) :: wavenumber

    ! Compute the power spectrum.    
    standardPower=+self%powerSpectrumPrimordialTransferred_%power             (wavenumber) &
         &        *self%cosmologicalMassVariance_          %powerNormalization(          )
    return
  end function standardPower
  
  double precision function standardPowerLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the power spectrum, ${\mathrm d}\ln P(k)/{\mathrm d}\ln k$, for $k=${\normalfont
    !% \ttfamily wavenumber} [Mpc$^{-1}$].
    implicit none
    class           (powerSpectrumStandard), intent(inout) :: self
    double precision                       , intent(in   ) :: wavenumber

    standardPowerLogarithmicDerivative=self%powerSpectrumPrimordialTransferred_%logarithmicDerivative(wavenumber)
    return
  end function standardPowerLogarithmicDerivative

  double precision function standardPowerDimensionless(self,wavenumber)
    !% Return the dimensionless power spectrum, $\Delta^2(k)$, for $k=${\normalfont \ttfamily wavenumber} [Mpc$^{-1}$].
    use Numerical_Constants_Math
    implicit none
    class           (powerSpectrumStandard), intent(inout) :: self
    double precision                       , intent(in   ) :: wavenumber
    
    standardPowerDimensionless=+4.0d0                       &
         &                       *Pi                        &
         &                       *           wavenumber **3 &
         &                       *self%power(wavenumber)    &
         &                       /(                         &
         &                         +2.0d0                   &
         &                         *Pi                      &
         &                        )**3
    return
  end function standardPowerDimensionless
