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
  An implementation of \cite{kitayama_semianalytic_1996} dark matter halo virial density contrasts.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <virialDensityContrast name="virialDensityContrastKitayamaSuto1996">
   <description>
  !!]

  !![
    A virial density contrast class using the fitting formula of \cite{kitayama_semianalytic_1996}, and so is valid only in
    flat cosmological models (an error will be reported in non-flat models). Specifically,
    \begin{equation}
     \Delta_\mathrm{virial}(t) = 18 \pi^2 [1+0.4093 \left\{{1\over \Omega_\mathrm{matter}(t)}-1\right\}(t)^{0.9052}].
    \end{equation}
   </description>
  </virialDensityContrast>
  !!]
  type, extends(virialDensityContrastClass) :: virialDensityContrastKitayamaSuto1996
     !!{
     A dark matter halo virial density contrast class using the fitting functions of \cite{kitayama_semianalytic_1996}.
     !!}
     private
     class(cosmologyFunctionsClass ), pointer :: cosmologyFunctions_ => null()
   contains
     final     ::                                kitayamaSuto1996Destructor
     procedure :: densityContrast             => kitayamaSuto1996DensityContrast
     procedure :: densityContrastRateOfChange => kitayamaSuto1996DensityContrastRateOfChange
  end type virialDensityContrastKitayamaSuto1996

  interface virialDensityContrastKitayamaSuto1996
     !!{
     Constructors for the \refClass{virialDensityContrastKitayamaSuto1996} dark matter halo virial density contrast class.
     !!}
     module procedure kitayamaSuto1996ConstructorParameters
     module procedure kitayamaSuto1996ConstructorInternal
  end interface virialDensityContrastKitayamaSuto1996

contains

  function kitayamaSuto1996ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialDensityContrastKitayamaSuto1996} dark matter halo virial density contrast class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (virialDensityContrastKitayamaSuto1996)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=virialDensityContrastKitayamaSuto1996(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function kitayamaSuto1996ConstructorParameters

  function kitayamaSuto1996ConstructorInternal(cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{virialDensityContrastKitayamaSuto1996} dark matter halo virial density contrast class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (virialDensityContrastKitayamaSuto1996)                        :: self
    class(cosmologyFunctionsClass              ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    return
  end function kitayamaSuto1996ConstructorInternal

  subroutine kitayamaSuto1996Destructor(self)
    !!{
    Destructor for the \refClass{virialDensityContrastKitayamaSuto1996} virial density contrast class.
    !!}
    implicit none
    type(virialDensityContrastKitayamaSuto1996), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine kitayamaSuto1996Destructor

  double precision function kitayamaSuto1996DensityContrast(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, assuming the fitting function of \cite{kitayama_semianalytic_1996}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (virialDensityContrastKitayamaSuto1996), intent(inout)           :: self
    double precision                                       , intent(in   )           :: mass
    double precision                                       , intent(in   ), optional :: time               , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    double precision                                                                 :: omegaf
    !$GLC attributes unused :: self, mass

    omegaf=max(0.0d0,1.0d0/self%cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)-1.0d0)
    kitayamaSuto1996DensityContrast=18.0d0*Pi**2*(1.0d0+0.4093d0*omegaf**0.9052d0)
    return
  end function kitayamaSuto1996DensityContrast

  double precision function kitayamaSuto1996DensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, assuming the fitting function of \cite{kitayama_semianalytic_1996}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (virialDensityContrastKitayamaSuto1996), intent(inout)           :: self
    double precision                                       , intent(in   )           :: mass
    double precision                                       , intent(in   ), optional :: time      , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    double precision                                                                 :: omegaf
    !$GLC attributes unused :: self, mass

    omegaf=max(0.0d0,1.0d0/self%cosmologyFunctions_%omegaMatterEpochal(time,expansionFactor,collapsing)-1.0d0)
    kitayamaSuto1996DensityContrastRateOfChange=                                                 &
         & -18.0d0*Pi**2                                                                         &
         & *(0.9052d0*0.4093d0*omegaf**(0.9052d0-1.0d0))                                         &
         & *self%cosmologyFunctions_%omegaMatterRateOfChange(time,expansionFactor,collapsing)    &
         & /self%cosmologyFunctions_%omegaMatterEpochal     (time,expansionFactor,collapsing)**2
   return
  end function kitayamaSuto1996DensityContrastRateOfChange
