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
Implements a identity modifier for power spectra in the halo model of clustering.
!!}

  !![
  <haloModelPowerSpectrumModifier name="haloModelPowerSpectrumModifierIdentity">
   <description>
    A halo model power spectrum modifier class which applies an identity modifier.
   </description>
  </haloModelPowerSpectrumModifier>
  !!]
  type, extends(haloModelPowerSpectrumModifierClass) :: haloModelPowerSpectrumModifierIdentity
     private
   contains
     procedure :: modify => identityModify
  end type haloModelPowerSpectrumModifierIdentity

  interface haloModelPowerSpectrumModifierIdentity
     !!{
     Constructors for the \refClass{haloModelPowerSpectrumModifierIdentity} halo model power spectrum modifier class.
     !!}
     module procedure identityConstructorParameters
  end interface haloModelPowerSpectrumModifierIdentity

contains

  function identityConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily identity} hot halo outflow reincorporation class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(haloModelPowerSpectrumModifierIdentity)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self=haloModelPowerSpectrumModifierIdentity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identityConstructorParameters

  subroutine identityModify(self,wavenumber,term,powerSpectrum,powerSpectrumCovariance,mass)
    !!{
    Applies a identity modification to a halo model power spectrum.
    !!}
    implicit none
    class           (haloModelPowerSpectrumModifierIdentity), intent(inout)                           :: self
    double precision                                        , intent(in   ), dimension(:  )           :: wavenumber
    type            (enumerationHaloModelTermType          ), intent(in   )                           :: term
    double precision                                        , intent(inout), dimension(:  )           :: powerSpectrum
    double precision                                        , intent(inout), dimension(:,:), optional :: powerSpectrumCovariance
    double precision                                        , intent(in   )                , optional :: mass
    !$GLC attributes unused :: self, wavenumber, term, powerSpectrum, mass

    ! Do nothing, except to set covariance to zero.
    if (present(powerSpectrumCovariance)) powerSpectrumCovariance=0.0d0
    return
  end subroutine identityModify
