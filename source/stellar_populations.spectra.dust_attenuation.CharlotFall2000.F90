!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Implements calculations of attenuation of stellar spectra using the model of \cite{charlot_simple_2000}.
  
  !# <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationCharlotFall2000">
  !#  <description>Returns the dust attenuation of stellar spectra according to the model of \cite{charlot_simple_2000}.</description>
  !# </stellarSpectraDustAttenuation>

  type, extends(stellarSpectraDustAttenuationClass) :: stellarSpectraDustAttenuationCharlotFall2000
     !% A class implementing calculations of attenuation of stellar spectra using the model of \cite{charlot_simple_2000}.
     private
     double precision :: opacityExponent, birthCloudLifetime, opticalDepthISM, opticalDepthBirthClouds
   contains
     procedure :: attenuation => charlotFall2000Attenuation
  end type stellarSpectraDustAttenuationCharlotFall2000

  interface stellarSpectraDustAttenuationCharlotFall2000
     !% Constructors for the ``charlotFall2000'' stellar spectra dust attenuation class.
     module procedure charlotFall2000DefaultConstructor
     module procedure charlotFall2000Constructor
  end interface stellarSpectraDustAttenuationCharlotFall2000

contains

  function charlotFall2000DefaultConstructor()
    !% Default constructor for the ``charlotFall2000'' stellar spectra dust attenuatio class.
    implicit none
    type            (stellarSpectraDustAttenuationCharlotFall2000)            :: charlotFall2000DefaultConstructor
    ! Exponent of the opacity wavelength dependence.
    double precision                                              , parameter :: opacityExponent        =    0.7d0
    ! Duration spent in birthclouds (Gyr).
    double precision                                              , parameter :: birthCloudLifetime     =    1.0d-2
    ! Effective optical depth of the ISM and birthclouds.
    double precision                                              , parameter :: opticalDepthISM        =    0.5d0 , &
         &                                                                       opticalDepthBirthClouds=    1.0d0

    charlotFall2000DefaultConstructor=charlotFall2000Constructor(opacityExponent,birthCloudLifetime,opticalDepthISM,opticalDepthBirthClouds)
    return
  end function charlotFall2000DefaultConstructor

  function charlotFall2000Constructor(opacityExponent,birthCloudLifetime,opticalDepthISM,opticalDepthBirthClouds)
    !% Constructor for the ``charlotFall2000'' stellar spectra dust attenuatio class.
    implicit none
    type            (stellarSpectraDustAttenuationCharlotFall2000)                :: charlotFall2000Constructor
    double precision                                              , intent(in   ) :: opacityExponent           , birthCloudLifetime     , &
         &                                                                           opticalDepthISM           , opticalDepthBirthClouds

    charlotFall2000Constructor%opacityExponent        =opacityExponent
    charlotFall2000Constructor%birthCloudLifetime     =birthCloudLifetime
    charlotFall2000Constructor%opticalDepthISM        =opticalDepthISM
    charlotFall2000Constructor%opticalDepthBirthClouds=opticalDepthBirthClouds
    return
  end function charlotFall2000Constructor

  double precision function charlotFall2000Attenuation(self,wavelength,age,vBandAttenuation)
    !% Return attenuation of stellar spectra according to the model of \cite{charlot_simple_2000}. Note that the V-band
    !% attenuation is taken to be that due to the ISM alone (i.e. not including birthclouds).
    implicit none
    class           (stellarSpectraDustAttenuationCharlotFall2000), intent(inout) :: self
    double precision                                              , intent(in   ) :: wavelength                                , age, &
         &                                                                           vBandAttenuation
    ! Effective wavelength of Buser V-band filter.
    double precision                                              , parameter     :: vBandWavelength        =5504.61227375652d0
    double precision                                                              :: cloudFactor

    if (age <= self%birthCloudLifetime) then
       cloudFactor=1.0d0+self%opticalDepthBirthClouds/self%opticalDepthISM
    else
       cloudFactor=1.0d0
    end if
    charlotFall2000Attenuation=vBandAttenuation*cloudFactor/(wavelength/vBandWavelength)**self%opacityExponent
    return
  end function charlotFall2000Attenuation
