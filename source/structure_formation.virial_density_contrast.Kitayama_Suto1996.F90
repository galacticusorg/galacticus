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

  !% An implementation of \cite{kitayama_semianalytic_1996} dark matter halo virial density contrasts.

  !# <virialDensityContrast name="virialDensityContrastKitayamaSuto1996">
  !#  <description>\cite{kitayama_semianalytic_1996} dark matter halo virial density contrasts.</description>
  !# </virialDensityContrast>

  type, extends(virialDensityContrastClass) :: virialDensityContrastKitayamaSuto1996
     !% A dark matter halo virial density contrast class using the fitting functions of \cite{kitayama_semianalytic_1996}.
     private
   contains
     procedure :: densityContrast             => kitayamaSuto1996DensityContrast
     procedure :: densityContrastRateOfChange => kitayamaSuto1996DensityContrastRateOfChange
  end type virialDensityContrastKitayamaSuto1996

  interface virialDensityContrastKitayamaSuto1996
     !% Constructors for the {\tt kitayamaSuto1996} dark matter halo virial density contrast class.
     module procedure kitayamaSuto1996DefaultConstructor
  end interface virialDensityContrastKitayamaSuto1996

contains

  function kitayamaSuto1996DefaultConstructor()
    !% Default constructor for the {\tt kitayamaSuto1996} dark matter halo virial density contrast class.
    implicit none
    type (virialDensityContrastKitayamaSuto1996), target  :: kitayamaSuto1996DefaultConstructor

    return
  end function kitayamaSuto1996DefaultConstructor

  double precision function kitayamaSuto1996DensityContrast(self,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, assuming the fitting function of \cite{kitayama_semianalytic_1996}.
    use Cosmology_Functions
    use Numerical_Constants_Math
    implicit none
    class           (virialDensityContrastKitayamaSuto1996), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: time               , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass              ), pointer                 :: cosmologyFunctions_
    double precision                                                                 :: omegaf

    cosmologyFunctions_ => cosmologyFunctions()
    omegaf=max(0.0d0,1.0d0/cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)-1.0d0)
    kitayamaSuto1996DensityContrast=18.0d0*Pi**2*(1.0d0+0.4093d0*omegaf**0.9052d0)
    return
  end function kitayamaSuto1996DensityContrast

  double precision function kitayamaSuto1996DensityContrastRateOfChange(self,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, assuming the fitting function of \cite{kitayama_semianalytic_1996}.
    use Cosmology_Functions
    use Numerical_Constants_Math
    implicit none
    class           (virialDensityContrastKitayamaSuto1996), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: time      , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass              ), pointer                 :: cosmologyFunctions_
    double precision                                                                 :: omegaf

    cosmologyFunctions_ => cosmologyFunctions()
    omegaf=max(0.0d0,1.0d0/cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)-1.0d0)
    kitayamaSuto1996DensityContrastRateOfChange=                                            &
         & -18.0d0*Pi**2                                                                    &
         & *(0.9052d0*0.4093d0*omegaf**(0.9052d0-1.0d0))                                    &
         & *cosmologyFunctions_%omegaMatterRateOfChange(time,expansionFactor,collapsing)    &
         & /cosmologyFunctions_%omegaMatterEpochal     (time,expansionFactor,collapsing)**2
   return
  end function kitayamaSuto1996DensityContrastRateOfChange
