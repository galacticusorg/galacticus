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
  Implements calculations of attenuation of stellar spectra using the model of \cite{charlot_simple_2000}.
  !!}

  !![
  <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationCharlotFall2000">
   <description>Returns the dust attenuation of stellar spectra according to the model of \cite{charlot_simple_2000}.</description>
  </stellarSpectraDustAttenuation>
  !!]
  type, extends(stellarSpectraDustAttenuationClass) :: stellarSpectraDustAttenuationCharlotFall2000
     !!{
     A class implementing calculations of attenuation of stellar spectra using the model of \cite{charlot_simple_2000}.
     !!}
     private
     double precision :: opacityExponent, birthCloudLifetime, opticalDepthISM, opticalDepthBirthClouds
   contains
     procedure :: attenuation    => charlotFall2000Attenuation
     procedure :: isAgeDependent => charlotFall2000IsAgeDependent
  end type stellarSpectraDustAttenuationCharlotFall2000

  interface stellarSpectraDustAttenuationCharlotFall2000
     !!{
     Constructors for the \refClass{stellarSpectraDustAttenuationCharlotFall2000} stellar spectra dust attenuation class.
     !!}
     module procedure charlotFall2000ConstructorParameters
     module procedure charlotFall2000ConstructorInternal
  end interface stellarSpectraDustAttenuationCharlotFall2000

contains

  function charlotFall2000ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the \mono{charlotFall2000} stellar spectra dust attenuation class.
    !!}
    implicit none
    type            (stellarSpectraDustAttenuationCharlotFall2000)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: opacityExponent, birthCloudLifetime     , &
         &                                                                           opticalDepthISM, opticalDepthBirthClouds

    !![
    <inputParameter>
      <name>opacityExponent</name>
      <defaultValue>0.7d0</defaultValue>
      <description>The power-law exponent of wavelength in the opacity $\tau \propto \lambda^{-n}$ in the \cite{charlot_simple_2000} two-component dust attenuation model; a value of 0.7 corresponds to the standard parameterization for interstellar and birth cloud attenuation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>birthCloudLifetime</name>
      <defaultValue>1.0d-2</defaultValue>
      <description>The lifetime of stellar birth clouds (in Gyr) in the \cite{charlot_simple_2000} dust model; stars younger than this age are attenuated by both the ISM component and the birth cloud component, while older stars experience only ISM attenuation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>opticalDepthISM</name>
      <defaultValue>0.5d0</defaultValue>
      <description>The effective V-band optical depth $\hat{\tau}_\mathrm{V}^\mathrm{ISM}$ of the diffuse interstellar medium in the \cite{charlot_simple_2000} dust attenuation model, applied to all stellar populations regardless of age.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>opticalDepthBirthClouds</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The effective V-band optical depth $\hat{\tau}_\mathrm{V}^\mathrm{BC}$ of the stellar birth clouds in the \cite{charlot_simple_2000} dust model, applied only to stellar populations younger than \mono{birthCloudLifetime}; must exceed \mono{opticalDepthISM}.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=stellarSpectraDustAttenuationCharlotFall2000(opacityExponent,birthCloudLifetime,opticalDepthISM,opticalDepthBirthClouds)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function charlotFall2000ConstructorParameters

  function charlotFall2000ConstructorInternal(opacityExponent,birthCloudLifetime,opticalDepthISM,opticalDepthBirthClouds) result(self)
    !!{
    Constructor for the \refClass{stellarSpectraDustAttenuationCharlotFall2000} stellar spectra dust attenuation class.
    !!}
    implicit none
    type            (stellarSpectraDustAttenuationCharlotFall2000)                :: self
    double precision                                              , intent(in   ) :: opacityExponent, birthCloudLifetime     , &
         &                                                                           opticalDepthISM, opticalDepthBirthClouds
    !![
    <constructorAssign variables="opacityExponent,birthCloudLifetime,opticalDepthISM,opticalDepthBirthClouds"/>
    !!]

    return
  end function charlotFall2000ConstructorInternal

  double precision function charlotFall2000Attenuation(self,wavelength,age,vBandAttenuation)
    !!{
    Return attenuation of stellar spectra according to the model of \cite{charlot_simple_2000}. Note that the V-band
    attenuation is taken to be that due to the ISM alone (i.e. not including birth clouds).
    !!}
    implicit none
    class           (stellarSpectraDustAttenuationCharlotFall2000), intent(inout) :: self
    double precision                                              , intent(in   ) :: wavelength                         , age, &
         &                                                                           vBandAttenuation
    ! Effective wavelength of Buser V-band filter.
    double precision                                              , parameter     :: vBandWavelength =5504.61227375652d0
    double precision                                                              :: cloudFactor

    if (age <= self%birthCloudLifetime) then
       cloudFactor=1.0d0+self%opticalDepthBirthClouds/self%opticalDepthISM
    else
       cloudFactor=1.0d0
    end if
    charlotFall2000Attenuation=vBandAttenuation*cloudFactor/(wavelength/vBandWavelength)**self%opacityExponent
    return
  end function charlotFall2000Attenuation

  logical function charlotFall2000IsAgeDependent(self)
    !!{
    Return true since attenuation is age-dependent in the \cite{charlot_simple_2000} dust attenuation model.
    !!}
    implicit none
    class(stellarSpectraDustAttenuationCharlotFall2000), intent(inout) :: self
    !$GLC attributes unused :: self

    charlotFall2000IsAgeDependent=.true.
    return
  end function charlotFall2000IsAgeDependent
