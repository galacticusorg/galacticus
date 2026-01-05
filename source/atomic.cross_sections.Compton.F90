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
Contains a module which implements calculation of the Compton cross-section.
!!}

module Atomic_Cross_Sections_Compton
  !!{
  Implements calculation of the Compton cross-section.
  !!}
  implicit none
  private
  public :: Atomic_Cross_Section_Compton

contains

  function Atomic_Cross_Section_Compton(photonEnergy)
    !!{
    Returns the Compton cross section (in cm$^2$) for the specified {\normalfont \ttfamily photonEnergy} (in keV) from \cite{klein_uber_1929}.
    !!}
    use :: Numerical_Constants_Physical, only : electronMass, speedLight, thomsonCrossSection
    use :: Numerical_Constants_Prefixes, only : hecto       , kilo
    use :: Numerical_Constants_Units   , only : electronVolt
    implicit none
    double precision, dimension(:                 ), intent(in   ) :: photonEnergy
    double precision, dimension(size(photonEnergy))                :: Atomic_Cross_Section_Compton                                             , eta
    double precision, parameter                                    :: electronEnergy              =electronMass*speedLight**2/electronVolt/kilo

    eta=photonEnergy/electronEnergy
    Atomic_Cross_Section_Compton=3.0d0*thomsonCrossSection*(hecto**2)*((1.0d0-(2.0d0*eta+2.0d0)/eta**2)*log(2.0d0*eta+1.0d0)&
         &+0.5d0+4.0d0/eta-0.5d0/(2.0d0*eta+1.0d0)**2)/8.0d0/eta
    return
  end function Atomic_Cross_Section_Compton

end module Atomic_Cross_Sections_Compton
