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

  !+    Contributions to this file made by: Omid Sameie.

  !!{RST
  Provides a power spectrum window function class that implements the smooth-:math:`k` space filter of :cite:t:`leo_new_2018`.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionSmoothKSpace" docformat="rst">
   <description>
   A smooth window function for filtering of power spectra in wavenumber space, defined as :math:`W(kR) = 1/[1+(kR)^\beta]`, providing a tunable, sharp-but-smooth transition between large and small scales. The shape exponent :math:`\beta` and normalization relative to a top-hat filter are set by ``[beta]`` and ``[normalization]``.
   </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionSmoothKSpace
     !!{RST
     A smooth power spectrum window function class.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: beta                          , normalization
    contains
     final     ::                      smoothKSpaceDestructor
     procedure :: value             => smoothKSpaceValue
     procedure :: wavenumberMaximum => smoothKSpaceWavenumberMaximum
  end type powerSpectrumWindowFunctionSmoothKSpace

  interface powerSpectrumWindowFunctionSmoothKSpace
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumWindowFunctionSmoothKSpace` power spectrum window function class.
     !!}
     module procedure smoothKSpaceConstructorParameters
     module procedure smoothKSpaceConstructorInternal
  end interface powerSpectrumWindowFunctionSmoothKSpace

contains

  function smoothKSpaceConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumWindowFunctionSmoothKSpace` power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionSmoothKSpace)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    double precision                                                         :: beta                , normalization

    !![
    <inputParameter docformat="rst">
      <name>normalization</name>
      <source>parameters</source>
      <defaultValue>3.5d0</defaultValue>
      <description>
      The parameter "normalization" is equivalent to the normalization for a sharp-:math:`k` filter. It serves as the ratio of mass scales of the object to the one in the spherical model: :math:`R/R_\mathrm{topHat}`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>3.7d0</defaultValue>
      <description>
      The parameter "beta" is defined as the exponent of ":math:`kR`" in the denominator of the window function: :math:`W(kR)= 1/[1+(kR)^\beta]`.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self = powerSpectrumWindowFunctionSmoothKSpace(cosmologyParameters_,beta,normalization)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function smoothKSpaceConstructorParameters

  function smoothKSpaceConstructorInternal(cosmologyParameters_,beta,normalization) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumWindowFunctionSmoothKSpace` power spectrum window function class.
    !!}
    implicit none
    type            (powerSpectrumWindowFunctionSmoothKSpace)                        :: self
    class           (cosmologyParametersClass               ), intent(in   ), target :: cosmologyParameters_
    double precision                                         , intent(in   )         :: beta                , normalization
    !![
    <constructorAssign variables="beta, normalization, *cosmologyParameters_"/>
    !!]

    return
  end function smoothKSpaceConstructorInternal

  subroutine smoothKSpaceDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`powerSpectrumWindowFunctionSmoothKSpace` power spectrum window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionSmoothKSpace), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine smoothKSpaceDestructor

  double precision function smoothKSpaceValue(self,wavenumber,smoothingMass,time)
    !!{RST
    Smooth-:math:`k` space power spectrum window function proposed in :cite:t:`leo_new_2018`. spectrum.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumWindowFunctionSmoothKSpace), intent(inout) :: self
    double precision                                         , intent(in   ) :: smoothingMass, wavenumber, &
         &                                                                      time
    double precision                                                         :: smoothRadius , x
    !$GLC attributes unused :: time

    smoothRadius=+(                                             &
         &         +3.0d0                                       &
         &         /4.0d0                                       &
         &         /Pi                                          &
         &         *smoothingMass                               &
         &         /self%cosmologyParameters_%OmegaMatter    () &
         &         /self%cosmologyParameters_%densityCritical() &
         &        )**(1.0d0/3.0d0)                              &
         &        /self%normalization
    x           =+wavenumber                                    &
         &       *smoothRadius
    if (x <= 0.0d0) then
       smoothKSpaceValue=+0.0d0
    else
       ! For larger x, use the full expression.
       smoothKSpaceValue=1.0d0/(1.0d0+x**self%beta)
    end if
    return
  end function smoothKSpaceValue

  double precision function smoothKSpaceWavenumberMaximum(self,smoothingMass)
    !!{RST
    Maximum wavenumber for a top hat in real space window function Fourier transformed into :math:`k`-space used in computing the variance of the power spectrum.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionSmoothKSpace), intent(inout) :: self
    double precision                                         , intent(in   ) :: smoothingMass
    double precision                                         , parameter     :: wavenumberLarge=huge(1.0d0) ! Effective infinity.
    !$GLC attributes unused :: self, smoothingMass

    smoothKSpaceWavenumberMaximum=wavenumberLarge
    return
  end function smoothKSpaceWavenumberMaximum
