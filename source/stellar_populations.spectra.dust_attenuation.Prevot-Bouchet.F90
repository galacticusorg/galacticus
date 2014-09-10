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

  !% Implements calculations of attenuation of stellar spectra using the model of \cite{prevot_typical_1984} and \cite{bouchet_visible_1985}.
  
  !# <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationPrevotBouchet">
  !#  <description>Returns the dust attenuation of stellar spectra according to the model of \cite{prevot_typical_1984} and \cite{bouchet_visible_1985}.</description>
  !# </stellarSpectraDustAttenuation>

  use Tables

  type, extends(stellarSpectraDustAttenuationTabulated) :: stellarSpectraDustAttenuationPrevotBouchet
     !% A class implementing calculations of attenuation of stellar spectra using the model of \cite{prevot_typical_1984} and \cite{bouchet_visible_1985}.
     private
   contains
  end type stellarSpectraDustAttenuationPrevotBouchet

  interface stellarSpectraDustAttenuationPrevotBouchet
     !% Constructors for the ``prevotBouchet'' stellar spectra dust attenuation class.
     module procedure prevotBouchetConstructor
  end interface stellarSpectraDustAttenuationPrevotBouchet

contains

  function prevotBouchetConstructor()
    !% Constructor for the ``prevotBouchet'' stellar spectra dust attenuation class. Data read directly from Table~3 of \cite{bouchet_visible_1985}.
    use Galacticus_Error
    implicit none
    type            (stellarSpectraDustAttenuationPrevotBouchet)                           :: prevotBouchetConstructor
    double precision                                            , parameter                :: Rv=2.7d0
    double precision                                            , parameter, dimension(27) :: ElambdaVOverEBV=[-2.56d0,-2.40d0,-2.11d0, 0.00d0, 1.00d0, 1.59d0, 2.28d0, 2.61d0, 2.96d0, 3.17d0, 3.49d0, 3.91d0, 4.28d0, 4.60d0, 5.31d0, 5.83d0, 6.40d0, 6.79d0, 6.89d0, 7.16d0, 7.74d0, 8.02d0, 8.53d0, 9.15d0, 9.36d0, 9.90d0,10.80d0]

    call prevotBouchetConstructor%attenuationTable%create([0.44d0,0.60d0,0.79d0,1.89d0,2.32d0,2.68d0,3.19d0,3.31d0,3.41d0,3.55d0,3.72d0,3.89d0,4.07d0,4.24d0,4.46d0,4.68d0,4.93d0,5.20d0,5.31d0,5.45d0,5.63d0,5.83d0,6.02d0,6.22d0,6.44d0,6.66d0,6.93d0],tableCount=1,extrapolationType=extrapolationTypeExtrapolate)
    call prevotBouchetConstructor%attenuationTable%populate(ElambdaVOverEBV/Rv+1.0d0)
    return
  end function prevotBouchetConstructor
