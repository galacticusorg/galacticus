!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a sharp $k$-space power spectrum window function class.
  use Cosmology_Parameters

  !# <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionSharpKSpace">
  !#  <description>A sharp $k$-space window function for filtering of power spectra.</description>
  !# </powerSpectrumWindowFunction>
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionSharpKSpace
     !% A sharp $k$-space power spectrum window function class.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: cutOffNormalization
   contains
     final     ::                      sharpKSpaceDestructor
     procedure :: value             => sharpKSpaceValue
     procedure :: wavenumberMaximum => sharpKSpaceWavenumberMaximum
  end type powerSpectrumWindowFunctionSharpKSpace

  interface powerSpectrumWindowFunctionSharpKSpace
     !% Constructors for the {\normalfont \ttfamily sharpKSpace} power spectrum window function class.
     module procedure sharpKSpaceConstructorParameters
     module procedure sharpKSpaceConstructorInternal
  end interface powerSpectrumWindowFunctionSharpKSpace

contains

  function sharpKSpaceConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily sharpKSpace} power spectrum window function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (powerSpectrumWindowFunctionSharpKSpace)                :: sharpKSpaceConstructorParameters
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyParametersClass              ), pointer       :: cosmologyParameters_
    type            (varying_string                        )                :: normalizationText
    character       (len=32                                )                :: normalizationChar
    double precision                                                        :: normalization
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>normalization</name>
    !#   <source>parameters</source>
    !#   <variable>normalizationText</variable>
    !#   <defaultValue>var_str('natural')</defaultValue>
    !#   <description>
    !#     The parameter $a$ in the relation $k_{\mathrm s} = a/r_{\mathrm s}$, where $k_{\mathrm s}$ is the cut-off wavenumber for
    !#     the sharp $k$-space window function and $r_{\mathrm s}$ is the radius of a sphere (in real-space) enclosing the
    !#     requested smoothing mass. Alternatively, a value of {\normalfont \ttfamily natural} will be supplied in which case the normalization
    !#     is chosen such that, in real-space, $W(r=0)=1$. This results in a contained mass of $M=6 \pi^2 \bar{\rho} k_{\mathrm s}^{-3}$.
    !#   </description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    if (normalizationText == 'natural') then
       normalization=0.0d0
    else
       normalizationChar=normalizationText
       read (normalizationChar,*) normalization
    end if
    sharpKSpaceConstructorParameters=sharpKSpaceConstructorInternal(cosmologyParameters_,normalization)
    return
  end function sharpKSpaceConstructorParameters

  function sharpKSpaceConstructorInternal(cosmologyParameters_,normalization)
    !% Internal constructor for the {\normalfont \ttfamily sharpKSpace} power spectrum window function class.
    use Numerical_Constants_Math
    implicit none
    type            (powerSpectrumWindowFunctionSharpKSpace)                        :: sharpKSpaceConstructorInternal
    class           (cosmologyParametersClass              ), target, intent(in   ) :: cosmologyParameters_    
    double precision                                                                :: normalization
    
    sharpKSpaceConstructorInternal%cosmologyParameters_ => cosmologyParameters_
    ! Compute normalization.
    if (normalization <= 0.0d0) then
       ! Compute the "natural" normalization.
       sharpKSpaceConstructorInternal%cutOffNormalization=                                &
            & +(                                                                          &
            &   +6.0d0                                                                    &
            &   *Pi                                                                   **2 &
            &   *sharpKSpaceConstructorInternal%cosmologyParameters_%OmegaMatter    ()    &
            &   *sharpKSpaceConstructorInternal%cosmologyParameters_%densityCritical()    &
            &  )**(1.0d0/3.0d0)
    else
       ! Use provided normalization.
       sharpKSpaceConstructorInternal%cutOffNormalization=                                &
            & +normalization                                                              &
            & *(                                                                          &
            &   +4.0d0                                                                    &
            &   *Pi                                                                       &
            &   *sharpKSpaceConstructorInternal%cosmologyParameters_%OmegaMatter    ()    &
            &   *sharpKSpaceConstructorInternal%cosmologyParameters_%densityCritical()    &
            &   /3.0d0                                                                    &
            &  )**(1.0d0/3.0d0)
    end if
    return
  end function sharpKSpaceConstructorInternal

  subroutine sharpKSpaceDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sharpKSpace} power spectrum window function class.
    implicit none
    type(powerSpectrumWindowFunctionSharpKSpace), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    return
  end subroutine sharpKSpaceDestructor

  double precision function sharpKSpaceValue(self,wavenumber,smoothingMass)
    !% Sharp $k$-space window function used in computing the variance of the power spectrum.
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
    !% Sharp $k$-space window function used in computing the variance of the power spectrum.
    implicit none
    class           (powerSpectrumWindowFunctionSharpKSpace), intent(inout) :: self
    double precision                                        , intent(in   ) :: smoothingMass

    sharpKSpaceWavenumberMaximum=self%cutOffNormalization/smoothingMass**(1.0d0/3.0d0)
    return
  end function sharpKSpaceWavenumberMaximum
