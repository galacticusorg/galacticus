!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a hybrid top-hat/sharp $k$-space power spectrum window function class.
  use Cosmology_Parameters

  !# <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionTopHatSharpKHybrid">
  !#  <description>A hybrid top-hat/sharp $k$-space window function for filtering of power spectra.</description>
  !# </powerSpectrumWindowFunction>
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionTopHatSharpKHybrid
     !% A hybrid top-hat/sharp $k$-space power spectrum window function class.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: cutOffNormalization           , radiiRatio
   contains
     final     ::                      topHatSharpKHybridDestructor
     procedure :: value             => topHatSharpKHybridValue
     procedure :: wavenumberMaximum => topHatSharpKHybridWavenumberMaximum
  end type powerSpectrumWindowFunctionTopHatSharpKHybrid

  interface powerSpectrumWindowFunctionTopHatSharpKHybrid
     !% Constructors for the {\normalfont \ttfamily topHatSharpKHybrid} power spectrum window function class.
     module procedure topHatSharpKHybridConstructorParameters
     module procedure topHatSharpKHybridConstructorInternal
  end interface powerSpectrumWindowFunctionTopHatSharpKHybrid

contains

  function topHatSharpKHybridConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily topHatSharpKHybrid} power spectrum window function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (powerSpectrumWindowFunctionTopHatSharpKHybrid)                :: topHatSharpKHybridConstructorParameters
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (cosmologyParametersClass                     ), pointer       :: cosmologyParameters_
    type            (varying_string                               )                :: normalizationText
    character       (len=32                                       )                :: normalizationChar
    double precision                                                               :: normalization                          , radiiRatio
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
    !# <inputParameter>
    !#   <name>radiiRatio</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>
    !#     The parameter $\beta$ in the relation $r_{\mathrm s}=\beta r_{\mathrm th}$ between $k$-space sharp and top-hat window
    !#     function radii in the hybrid window function used for computing the variance in the power spectrum.
    !#   </description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    if (normalizationText == 'natural') then
       normalization=0.0d0
    else
       normalizationChar=normalizationText
       read (normalizationChar,*) normalization
    end if
    topHatSharpKHybridConstructorParameters=topHatSharpKHybridConstructorInternal(cosmologyParameters_,normalization,radiiRatio)
    return
  end function topHatSharpKHybridConstructorParameters

  function topHatSharpKHybridConstructorInternal(cosmologyParameters_,normalization,radiiRatio)
    !% Internal constructor for the {\normalfont \ttfamily topHatSharpKHybrid} power spectrum window function class.
    use Numerical_Constants_Math
    implicit none
    type            (powerSpectrumWindowFunctionTopHatSharpKHybrid)                        :: topHatSharpKHybridConstructorInternal
    class           (cosmologyParametersClass                     ), target, intent(in   ) :: cosmologyParameters_    
    double precision                                                                       :: normalization                        , radiiRatio
    
    topHatSharpKHybridConstructorInternal%cosmologyParameters_ => cosmologyParameters_
    topHatSharpKHybridConstructorInternal%radiiRatio           =  radiiRatio
    ! Compute normalization.
    if (normalization <= 0.0d0) then
       ! Compute the "natural" normalization.
       topHatSharpKHybridConstructorInternal%cutOffNormalization=                                &
            & +(                                                                                 &
            &   +6.0d0                                                                           &
            &   *Pi                                                                          **2 &
            &   *topHatSharpKHybridConstructorInternal%cosmologyParameters_%OmegaMatter    ()    &
            &   *topHatSharpKHybridConstructorInternal%cosmologyParameters_%densityCritical()    &
            &  )**(1.0d0/3.0d0)
    else
       ! Use provided normalization.
       topHatSharpKHybridConstructorInternal%cutOffNormalization=                                &
            & +normalization                                                                     &
            & *(                                                                                 &
            &   +4.0d0                                                                           &
            &   *Pi                                                                              &
            &   *topHatSharpKHybridConstructorInternal%cosmologyParameters_%OmegaMatter    ()    &
            &   *topHatSharpKHybridConstructorInternal%cosmologyParameters_%densityCritical()    &
            &   /3.0d0                                                                           &
            &  )**(1.0d0/3.0d0)
    end if
    return
  end function topHatSharpKHybridConstructorInternal

  subroutine topHatSharpKHybridDestructor(self)
    !% Destructor for the {\normalfont \ttfamily topHatSharpKHybrid} power spectrum window function class.
    implicit none
    type(powerSpectrumWindowFunctionTopHatSharpKHybrid), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologyParameters_"/>
    return
  end subroutine topHatSharpKHybridDestructor

  double precision function topHatSharpKHybridValue(self,wavenumber,smoothingMass)
    !% Computes a window function for calculations of the variance in the power spectrum. Specifically, uses a convolution of
    !% top-hat real-space and sharp $k$-space window functions. The top-hat radius is $r_{\mathrm th}$, while the $k$-space cut-off
    !% wavenumber is $k_{\mathrm s}=a/r_{\mathrm s}$, where $a=${\normalfont \ttfamily [normalization]}. The two radii
    !% are chosen such that $r_{\mathrm th}^2 + r_{\mathrm s}^2 = (3 M / 4 \pi \bar{rho})^{1/3}$ and $r_{\mathrm s}=\beta r_{\mathrm th}$ where
    !% $\beta=${\normalfont \ttfamily [pradiiRatio]}.
    use Numerical_Constants_Math
    implicit none
    class           (powerSpectrumWindowFunctionTopHatSharpKHybrid), intent(inout) :: self
    double precision                                               , intent(in   ) :: smoothingMass              , wavenumber
    double precision                                               , parameter     :: xSeriesMaximum      =1.0d-3
    double precision                                                               :: kSpaceSharpRadius          , topHatRadius, &
         &                                                                            topHatWindowFunction       , totalRadius , &
         &                                                                            wavenumberCutOff           , x           , &
         &                                                                            xSquared

    ! Find the radius enclosing this mass.
    totalRadius      =+(                                             &
         &              +3.0d0                                       &
         &              /4.0d0                                       &
         &              /Pi                                          &
         &              *smoothingMass                               &
         &              /self%cosmologyParameters_%OmegaMatter    () &
         &              /self%cosmologyParameters_%densityCritical() &
         &              )**(1.0d0/3.0d0)
    ! Find the top-hat and sharp k-space radii, and the k-space wavenumber.
    topHatRadius     =+totalRadius                                   &
         &            /sqrt(                                         &
         &                  +1.0d0                                   &
         &                  +self%radiiRatio**2                      &
         &                 )
    kSpaceSharpRadius=+self%radiiRatio                               &
         &            *topHatRadius
    wavenumberCutOff =+self%cutOffNormalization                      &
         &            /kSpaceSharpRadius
    ! Compute the top-hat window function.
    x                =+wavenumber*topHatRadius
    if      (x <= 0.0d0) then
       topHatWindowFunction=0.0d0
    else if (x <= xSeriesMaximum) then
       ! Use a series expansion of the window function for small x.
       xSquared            =+x**2
       topHatWindowFunction=+1.0d0                        &
            &               +xSquared*(  -1.0d0/   10.0d0 &
            &               +xSquared* ( +1.0d0/  280.0d0 &
            &               +xSquared*  (-1.0d0/15120.0d0 &
            &                           )                 &
            &                          )                  &
            &                         )
    else
       ! For larger x, use the full expression.
       topHatWindowFunction=3.0d0*(sin(x)-x*cos(x))/(x**3)
    end if
    ! Compute k-space sharp window function.
    if      (wavenumber <=            0.0d0) then
       wavenumberCutOff=0.0d0
    else if (wavenumber <= wavenumberCutOff) then
       wavenumberCutOff=1.0d0
    else
       wavenumberCutOff=0.0d0
    end if
    ! Compute the convolution (which is just the multiplication in k-space).
    topHatSharpKHybridValue=wavenumberCutOff*topHatWindowFunction
    return
  end function topHatSharpKHybridValue

  double precision function topHatSharpKHybridWavenumberMaximum(self,smoothingMass)
    !% Computes the maximum wavenumber at which the window function for calculations of the variance in the power spectrum is
    !% non-zero. Specifically, uses a convolution of top-hat real-space and sharp $k$-space window functions. The top-hat radius
    !% is $r_{\mathrm th}$, while the $k$-space cut-off wavenumber is $k_{\mathrm s}=a/r_{\mathrm s}$, where $a=${\tt
    !% [normalization]}. The two radii are chosen such that $r_{\mathrm th}^2 + r_{\mathrm s}^2 = (3 M / 4 \pi \bar{rho})^{1/3}$
    !% and $r_{\mathrm s}=\beta r_{\mathrm th}$ where $\beta=${\tt [radiiRatio]}.
    implicit none
    class           (powerSpectrumWindowFunctionTopHatSharpKHybrid), intent(inout) :: self
    double precision                                               , intent(in   ) :: smoothingMass
    double precision                                               , parameter     :: wavenumberLarge=1.0d30 !    Effective infinity.

    topHatSharpKHybridWavenumberMaximum=wavenumberLarge
    return
  end function topHatSharpKHybridWavenumberMaximum

