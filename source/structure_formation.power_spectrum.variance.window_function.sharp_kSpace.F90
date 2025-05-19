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
Implements a sharp $k$-space power spectrum window function class.
!!}
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionSharpKSpace">
   <description>
    A sharp $k$-space window function for filtering of power spectra. The window function is given by:
    \begin{equation}
     W(k) = \left\{ \begin{array}{ll} 1 &amp; \hbox{if } k &lt; k_\mathrm{s} \\ 0 &amp; \hbox{if } k &gt; k_\mathrm{s}, \end{array} \right.
    \end{equation}
    where if {\normalfont \ttfamily [normalization]}$=${\normalfont \ttfamily natural} then $k_\mathrm{s} = (6 \Pi^2 \bar{\rho}
    / M)^{1/3}$ for a smoothing scale $M$ and mean matter density $\bar{\rho}$. Otherwise, {\normalfont \ttfamily
    [normalization]} must be set to a numerical value, $\alpha$, in which case $k_\mathrm{s} = \alpha / R_\mathrm{th}$ with
    $R_\mathrm{th}=3M/4\pi\bar{\rho}$ for a smoothing scale $M$ and mean matter density $\bar{\rho}$.
   </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionSharpKSpace
     !!{
     A sharp $k$-space power spectrum window function class.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: cutOffNormalization
     type            (varying_string          )          :: normalization
   contains
     final     ::                               sharpKSpaceDestructor
     procedure :: value                      => sharpKSpaceValue
     procedure :: wavenumberMaximum          => sharpKSpaceWavenumberMaximum
     procedure :: amplitudeIsMassIndependent => sharpKSpaceAmplitudeIsMassIndependent
  end type powerSpectrumWindowFunctionSharpKSpace

  interface powerSpectrumWindowFunctionSharpKSpace
     !!{
     Constructors for the {\normalfont \ttfamily sharpKSpace} power spectrum window function class.
     !!}
     module procedure sharpKSpaceConstructorParameters
     module procedure sharpKSpaceConstructorInternal
  end interface powerSpectrumWindowFunctionSharpKSpace

contains

  function sharpKSpaceConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sharpKSpace} power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionSharpKSpace)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyParametersClass              ), pointer       :: cosmologyParameters_
    type            (varying_string                        )                :: normalization
    character       (len=32                                )                :: normalizationChar
    double precision                                                        :: normalizationValue

    ! Check parameters.
    !![
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <variable>normalization</variable>
      <defaultValue>var_str('natural')</defaultValue>
      <description>
        The parameter $a$ in the relation $k_\mathrm{s} = a/r_\mathrm{s}$, where $k_\mathrm{s}$ is the cut-off wavenumber for
        the sharp $k$-space window function and $r_\mathrm{s}$ is the radius of a sphere (in real-space) enclosing the
        requested smoothing mass. Alternatively, a value of {\normalfont \ttfamily natural} will be supplied in which case the normalization
        is chosen such that, in real-space, $W(r=0)=1$. This results in a contained mass of $M=6 \pi^2 \bar{\rho} k_\mathrm{s}^{-3}$.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    if (normalization == 'natural') then
       normalizationValue=0.0d0
    else
       normalizationChar=char(normalization)
       read (normalizationChar,*) normalizationValue
    end if
    self=sharpKSpaceConstructorInternal(cosmologyParameters_,normalizationValue)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function sharpKSpaceConstructorParameters

  function sharpKSpaceConstructorInternal(cosmologyParameters_,normalization) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily sharpKSpace} power spectrum window function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (powerSpectrumWindowFunctionSharpKSpace)                        :: self
    class           (cosmologyParametersClass              ), target, intent(in   ) :: cosmologyParameters_
    double precision                                                , intent(in   ) :: normalization
    character       (len=18                                )                        :: normalizationText
    !![
    <constructorAssign variables="*cosmologyParameters_"/>
    !!]

    ! Compute normalization.
    if (normalization <= 0.0d0) then
       ! Compute the "natural" normalization.
       self%cutOffNormalization=                                &
            & +(                                                &
            &   +6.0d0                                          &
            &   *Pi                                         **2 &
            &   *self%cosmologyParameters_%OmegaMatter    ()    &
            &   *self%cosmologyParameters_%densityCritical()    &
            &  )**(1.0d0/3.0d0)
       self%normalization='natural'
    else
       ! Use provided normalization.
       self%cutOffNormalization=                                &
            & +normalization                                    &
            & *(                                                &
            &   +4.0d0                                          &
            &   *Pi                                             &
            &   *self%cosmologyParameters_%OmegaMatter    ()    &
            &   *self%cosmologyParameters_%densityCritical()    &
            &   /3.0d0                                          &
            &  )**(1.0d0/3.0d0)
       write (normalizationText,'(e17.10)') normalization
       self%normalization=trim(normalizationText)
    end if
    return
  end function sharpKSpaceConstructorInternal

  subroutine sharpKSpaceDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily sharpKSpace} power spectrum window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionSharpKSpace), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine sharpKSpaceDestructor

  double precision function sharpKSpaceValue(self,wavenumber,smoothingMass)
    !!{
    Sharp $k$-space window function used in computing the variance of the power spectrum.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionSharpKSpace), intent(inout) :: self
    double precision                                        , intent(in   ) :: smoothingMass   , wavenumber
    double precision                                                        :: wavenumberCutOff

    wavenumberCutOff=self%wavenumberMaximum(smoothingMass)
    if      (wavenumber <=            0.0d0) then
       sharpKSpaceValue=0.0d0
    else if (wavenumber <= wavenumberCutOff) then
       sharpKSpaceValue=1.0d0
    else
       sharpKSpaceValue=0.0d0
    end if
    return
  end function sharpKSpaceValue

  double precision function sharpKSpaceWavenumberMaximum(self,smoothingMass)
    !!{
    Sharp $k$-space window function used in computing the variance of the power spectrum.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionSharpKSpace), intent(inout) :: self
    double precision                                        , intent(in   ) :: smoothingMass

    sharpKSpaceWavenumberMaximum=self%cutOffNormalization/smoothingMass**(1.0d0/3.0d0)
    return
  end function sharpKSpaceWavenumberMaximum

  logical function sharpKSpaceAmplitudeIsMassIndependent(self)
    !!{
    Indicate the the sharp $k$-space power spectrum window function has constant amplitude below the maximum wavenumber.
    !!}
    implicit none
    class(powerSpectrumWindowFunctionSharpKSpace), intent(inout) :: self
    !$GLC attributes unused :: self

    sharpKSpaceAmplitudeIsMassIndependent=.true.
    return
  end function sharpKSpaceAmplitudeIsMassIndependent

