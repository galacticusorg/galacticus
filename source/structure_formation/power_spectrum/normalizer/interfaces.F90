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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module which provides interfaces used by ``powerSpectrumNormalizer`` classes.
!!}

module Power_Spectrum_Normalizer_Interfaces
  !!{RST
  Provides interfaces used by ``powerSpectrumNormalizer`` classes. These are kept separate from the class module so that they
  can be used by low-level modules (such as ``Functions_Global``) without creating a dependency on the class module.
  !!}
  private
  public :: powerSpectrumNormalizerPowerUnnormalized

  abstract interface
     double precision function powerSpectrumNormalizerPowerUnnormalized(wavenumber,time)
       !!{RST
       Interface for functions returning the unnormalized power spectrum to a :galacticus-class:`powerSpectrumNormalizer`.
       !!}
       double precision, intent(in   ) :: wavenumber, time
     end function powerSpectrumNormalizerPowerUnnormalized
  end interface

end module Power_Spectrum_Normalizer_Interfaces
