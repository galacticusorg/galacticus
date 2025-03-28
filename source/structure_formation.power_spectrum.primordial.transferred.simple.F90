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
  A simple transferred primordial power spectrum class.
  !!}

  use :: Linear_Growth           , only : linearGrowthClass
  use :: Power_Spectra_Primordial, only : powerSpectrumPrimordialClass
  use :: Transfer_Functions      , only : transferFunctionClass

  !![
  <powerSpectrumPrimordialTransferred name="powerSpectrumPrimordialTransferredSimple">
   <description>Implements a simple transferred primordial power spectrum.</description>
  </powerSpectrumPrimordialTransferred>
  !!]
  type, extends(powerSpectrumPrimordialTransferredClass) :: powerSpectrumPrimordialTransferredSimple
     !!{
     A simple transferred primordial power spectrum class.
     !!}
     private
     class           (transferFunctionClass       ), pointer :: transferFunction_        => null()
     class           (linearGrowthClass           ), pointer :: linearGrowth_            => null()
     class           (powerSpectrumPrimordialClass), pointer :: powerSpectrumPrimordial_ => null()
     double precision                                        :: timeTransferFunction
   contains
     final     ::                                simpleDestructor
     procedure :: power                       => simplePower
     procedure :: logarithmicDerivative       => simpleLogarithmicDerivative
     procedure :: growthIsWavenumberDependent => simpleGrowthIsWavenumberDependent
  end type powerSpectrumPrimordialTransferredSimple

  interface powerSpectrumPrimordialTransferredSimple
     !!{
     Constructors for the {\normalfont \ttfamily simple} transferred primordial power spectrum class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface powerSpectrumPrimordialTransferredSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily simple} transferred primordial power spectrum class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (powerSpectrumPrimordialTransferredSimple)                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(transferFunctionClass                   ), pointer       :: transferFunction_
    class(powerSpectrumPrimordialClass            ), pointer       :: powerSpectrumPrimordial_
    class(linearGrowthClass                       ), pointer       :: linearGrowth_

    !![
    <objectBuilder class="powerSpectrumPrimordial" name="powerSpectrumPrimordial_" source="parameters"/>
    <objectBuilder class="transferFunction"        name="transferFunction_"        source="parameters"/>
    <objectBuilder class="linearGrowth"            name="linearGrowth_"            source="parameters"/>
    !!]
    self=powerSpectrumPrimordialTransferredSimple(powerSpectrumPrimordial_,transferFunction_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="powerSpectrumPrimordial_"/>
    <objectDestructor name="transferFunction_"       />
    <objectDestructor name="linearGrowth_"           />
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(powerSpectrumPrimordial_,transferFunction_,linearGrowth_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily simple} transferred primordial power spectrum class.
    !!}
    implicit none
    type (powerSpectrumPrimordialTransferredSimple)                        :: self
    class(powerSpectrumPrimordialClass            ), intent(in   ), target :: powerSpectrumPrimordial_
    class(transferFunctionClass                   ), intent(in   ), target :: transferFunction_
    class(linearGrowthClass                       ), intent(in   ), target :: linearGrowth_
    !![
    <constructorAssign variables="*powerSpectrumPrimordial_, *transferFunction_, *linearGrowth_"/>
    !!]

    self%timeTransferFunction=self%transferFunction_%epochTime()
    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily simple} transferred primordial power spectrum class.
    !!}
    implicit none
    type(powerSpectrumPrimordialTransferredSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%transferFunction_"       />
    <objectDestructor name="self%powerSpectrumPrimordial_"/>
    <objectDestructor name="self%linearGrowth_"           />
    !!]
   return
  end subroutine simpleDestructor

  double precision function simplePower(self,wavenumber,time)
    !!{
    Return the transferred primordial power spectrum at the given {\normalfont \ttfamily
    wavenumber}.
    !!}
    implicit none
    class           (powerSpectrumPrimordialTransferredSimple), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber, time

    simplePower=+(                                                                                           &
         &        +self%transferFunction_       %value(wavenumber=wavenumber)                                &
         &        *self%linearGrowth_           %value(wavenumber=wavenumber,time=     time                ) &
         &        /self%linearGrowth_           %value(wavenumber=wavenumber,time=self%timeTransferFunction) &
         &       )**2                                                                                        &
         &      *  self%powerSpectrumPrimordial_%power(wavenumber=wavenumber)
    return
  end function simplePower

  double precision function simpleLogarithmicDerivative(self,wavenumber,time)
    !!{
    Return the logarithmic derivative of the transferred primordial power spectrum at the
    given {\normalfont \ttfamily wavenumber}.
    !!}
    implicit none
    class           (powerSpectrumPrimordialTransferredSimple), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber, time

    simpleLogarithmicDerivative=+2.0d0*self%transferFunction_       %logarithmicDerivative          (wavenumber=wavenumber                               ) &
         &                      +2.0d0*self%linearGrowth_           %logarithmicDerivativeWavenumber(wavenumber=wavenumber,time=     time                ) &
         &                      -2.0d0*self%linearGrowth_           %logarithmicDerivativeWavenumber(wavenumber=wavenumber,time=self%timeTransferFunction) &
         &                      +      self%powerSpectrumPrimordial_%logarithmicDerivative          (wavenumber=wavenumber                               )
    return
  end function simpleLogarithmicDerivative

  logical function simpleGrowthIsWavenumberDependent(self)
    !!{
    Return true if the growth of the power spectrum is wavenumber-dependent.
    !!}
    implicit none
    class(powerSpectrumPrimordialTransferredSimple), intent(inout) :: self

    simpleGrowthIsWavenumberDependent=self%linearGrowth_%isWavenumberDependent()
    return
  end function simpleGrowthIsWavenumberDependent
