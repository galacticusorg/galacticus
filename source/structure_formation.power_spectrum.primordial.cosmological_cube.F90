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
  A primordial power spectrum class which truncates a primordial power spectrum into a cosmological simulation cube.
  !!}

  !![
  <powerSpectrumPrimordial name="powerSpectrumPrimordialCosmologicalCube">
   <description>
     A primordial power spectrum class which truncates a primordial power spectrum into a cosmological simulation
     cube. Specifically, it assumes that wave-vectors with all components ($k_x,k_y,k_z$) smaller than $\Delta k = 2 \pi f / L$,
     where $L=${\normalfont \ttfamily [lengthCube]} is the simulation cube length, and $f=${\normalfont \ttfamily
     wavenumberMinimumFactor} specifies the minimum wavenumber (in units of $2 \pi / L$) to which the power spectrum is
     integrated, are missing from the power spectrum. For wavenumbers of magnitude $k$ the fraction of power missed due to these
     missing wavenumbers is computed. The total power at that wavenumber is then reduced by that amount.

     The fractional suppression in power, $f(x)$, as a function of $x=k/\Delta k$ is given by:
     \begin{equation}
     f(x) = \left\{ \begin{array}{ll} 0 &amp; \hbox{ for } x \le 1, \\ 3(1-x^{-1}) &amp; \hbox{ for } 1 &lt; x \le \sqrt{2}, \\ 3 x{-1} [ 1 + 2 \pi^{-1} x \sin^{-1}(\{x^2-1\}^{-1}) - 4 \pi^{-1} x \sin^{-1}(\{1-x^2\}^{-1/2})  ] &amp; \hbox{ for } \sqrt{2} &lt; x \le \sqrt{3}, \\ 1 &amp; \hbox{ for } \sqrt{3} &lt; x. \end{array} \right. 
     \end{equation}

     For the $1 &lt; x \le \sqrt{2}$ the solution is found by considering the 6 spherical caps of the sphere which protrude from the
     faces of the cube. The solution for $\sqrt{2} &lt; x \le \sqrt{3}$ is more complicated and follows the solution given by
     \cite{achille2013}.
   </description>
  </powerSpectrumPrimordial>
  !!]
  type, extends(powerSpectrumPrimordialClass) :: powerSpectrumPrimordialCosmologicalCube
     !!{
     A primordial power spectrum class which truncates a primordial power spectrum into a cosmological simulation cube.
     !!}
     private
     class           (powerSpectrumPrimordialClass), pointer :: powerSpectrumPrimordial_ => null()
     double precision                                        :: lengthCube                        , wavenumberMinimum, &
          &                                                     wavenumberMinimumFactor
   contains
     final :: cosmologicalCubeDestructor
     procedure :: power                 => cosmologicalCubePower
     procedure :: logarithmicDerivative => cosmologicalCubeLogarithmicDerivative
  end type powerSpectrumPrimordialCosmologicalCube

  interface powerSpectrumPrimordialCosmologicalCube
     !!{
     Constructors for the \refClass{powerSpectrumPrimordialCosmologicalCube} primordial power spectrum class.
     !!}
     module procedure cosmologicalCubeConstructorParameters
     module procedure cosmologicalCubeConstructorInternal
  end interface powerSpectrumPrimordialCosmologicalCube

contains

  function cosmologicalCubeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{powerSpectrumPrimordialCosmologicalCube} primordial power spectrum class which takes a parameter set as input.
    !!}
    implicit none
    type            (powerSpectrumPrimordialCosmologicalCube)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (powerSpectrumPrimordialClass           ), pointer       :: powerSpectrumPrimordial_
    double precision                                                         :: lengthCube              , wavenumberMinimumFactor
    
    !![
    <inputParameter>
      <name>lengthCube</name>
      <source>parameters</source>
      <description>The length of the cosmological cube.</description>
    </inputParameter>
    <inputParameter>
      <name>wavenumberMinimumFactor</name>      
      <source>parameters</source>
      <defaultValue>0.5d0</defaultValue>
      <defaultSource>The default value of $1/2$ assumes that the power spectrum is integrated over cubic regions $\pm \pi/L$ centered on each grid point in the Fourier transform of the density field.</defaultSource>
      <description>The minimum wavenumber (in units of $2 \pi / L$) to which the power spectrum is integrated.</description>
    </inputParameter>
    <objectBuilder class="powerSpectrumPrimordial" name="powerSpectrumPrimordial_" source="parameters"/>
    !!]
    self=powerSpectrumPrimordialCosmologicalCube(lengthCube,wavenumberMinimumFactor,powerSpectrumPrimordial_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="powerSpectrumPrimordial_"/>
    !!]
    return
  end function cosmologicalCubeConstructorParameters

  function cosmologicalCubeConstructorInternal(lengthCube,wavenumberMinimumFactor,powerSpectrumPrimordial_) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumPrimordialCosmologicalCube} primordial power spectrum class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (powerSpectrumPrimordialCosmologicalCube)                        :: self
    double precision                                         , intent(in   )         :: lengthCube              , wavenumberMinimumFactor
    class           (powerSpectrumPrimordialClass           ), intent(in   ), target :: powerSpectrumPrimordial_
    !![
    <constructorAssign variables="lengthCube, wavenumberMinimumFactor, *powerSpectrumPrimordial_"/>
    !!]
    
    self%wavenumberMinimum=+2.0d0                        &
         &                 *Pi                           &
         &                 *self%wavenumberMinimumFactor &
         &                 /self%lengthCube
    return
  end function cosmologicalCubeConstructorInternal

  subroutine cosmologicalCubeDestructor(self)
    !!{
    Destructor for the \refClass{powerSpectrumPrimordialCosmologicalCube} primordial power spectrum class. 
    !!}
    implicit none
    type(powerSpectrumPrimordialCosmologicalCube), intent(inout) :: self

    !![
    <objectDestructor name="self%powerSpectrumPrimordial_"/>
    !!]
    return
  end subroutine cosmologicalCubeDestructor

  double precision function cosmologicalCubePower(self,wavenumber)
    !!{
    Return the primordial power spectrum at the given {\normalfont \ttfamily wavenumber}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumPrimordialCosmologicalCube), intent(inout) :: self
    double precision                                         , intent(in   ) :: wavenumber
    double precision                                                         :: wavenumberFractional

    wavenumberFractional =+     wavenumber        &
         &                /self%wavenumberMinimum
    if (wavenumberFractional < 1.0d0) then
       ! All modes are within the cubic exclusion region, so power is fully suppressed.
       cosmologicalCubePower=+0.0d0
    else
       cosmologicalCubePower=+self%powerSpectrumPrimordial_%power(wavenumber)
       if      (wavenumberFractional < sqrt(2.0d0)) then
          cosmologicalCubePower=+cosmologicalCubePower  &
               &                *3.0d0                  &
               &                *(                      &
               &                  +1.0d0                &
               &                  -1.0d0                &
               &                  /wavenumberFractional &
               &                 )
       else if (wavenumberFractional < sqrt(3.0d0)) then
          cosmologicalCubePower=+cosmologicalCubePower                                                           &
               &                *3.0d0                                                                           &
               &                /wavenumberFractional                                                            &
               &                *(                                                                               &
               &                  +1.0d0                                                                         &
               &                  +2.0d0*wavenumberFractional/Pi*asin(1.0d0/    (wavenumberFractional**2-1.0d0)) &
               &                  -4.0d0                     /Pi*asin(1.0d0/sqrt(wavenumberFractional**2-1.0d0)) &
               &                 )
       end if
    end if
    return
  end function cosmologicalCubePower

  double precision function cosmologicalCubeLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the primordial power spectrum at the given {\normalfont \ttfamily wavenumber}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumPrimordialCosmologicalCube), intent(inout) :: self
    double precision                                         , intent(in   ) :: wavenumber
    double precision                                                         :: wavenumberFractional

    wavenumberFractional =+     wavenumber        &
         &                /self%wavenumberMinimum
    if (wavenumberFractional < 1.0d0) then
       ! All modes are within the cubic exclusion region, so power is fully suppressed.
       cosmologicalCubeLogarithmicDerivative=+0.0d0
    else
       cosmologicalCubeLogarithmicDerivative=+self%powerSpectrumPrimordial_%logarithmicDerivative(wavenumber)
       if      (wavenumberFractional < sqrt(2.0d0)) then
          cosmologicalCubeLogarithmicDerivative=+cosmologicalCubeLogarithmicDerivative &
               &                                +1.0d0                                 &
               &                                /(                                     &
               &                                 +wavenumberFractional                 &
               &                                 -1.0d0                                &
               &                                 )
       else if (wavenumberFractional < sqrt(3.0d0)) then
          cosmologicalCubeLogarithmicDerivative=+cosmologicalCubeLogarithmicDerivative                                        &
               &                                -  1.0d0                                                                      &
               &                                -  2.0d0*wavenumberFractional*asin(1.0d0/    (wavenumberFractional**2-1.0d0)) &
               &                                /(                                                                            &
               &                                  +Pi                                                                         &
               &                                  +2.0d0*wavenumberFractional*asin(1.0d0/    (wavenumberFractional**2-1.0d0)) &
               &                                  -4.0d0                     *asin(1.0d0/sqrt(wavenumberFractional**2-1.0d0)) &
               &                                 )
       end if
    end if
    return
  end function cosmologicalCubeLogarithmicDerivative
