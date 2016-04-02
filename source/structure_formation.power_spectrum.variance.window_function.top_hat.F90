!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements a top-hat power spectrum window function class.
  use Cosmology_Parameters

  !# <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionTopHat">
  !#  <description>A top-hat in real space window function for filtering of power spectra.</description>
  !# </powerSpectrumWindowFunction>
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionTopHat
     !% A top-hat power spectrum window function class.
     private
     class(cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
    contains
     final     ::                      topHatDestructor
     procedure :: descriptor        => topHatDescriptor
     procedure :: value             => topHatValue
     procedure :: wavenumberMaximum => topHatWavenumberMaximum
  end type powerSpectrumWindowFunctionTopHat

  interface powerSpectrumWindowFunctionTopHat
     !% Constructors for the {\normalfont \ttfamily topHat} power spectrum window function class.
     module procedure topHatConstructorParameters
     module procedure topHatConstructorInternal
  end interface powerSpectrumWindowFunctionTopHat

contains

  function topHatConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily topHat} power spectrum window function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(powerSpectrumWindowFunctionTopHat)                :: topHatConstructorParameters
    type(inputParameters                  ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />    
    
    ! Check parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <objectBuilder class="cosmologyParameters" name="topHatConstructorParameters%cosmologyParameters_" source="parameters"/>
    return
  end function topHatConstructorParameters

  function topHatConstructorInternal(cosmologyParameters_)
    !% Internal constructor for the {\normalfont \ttfamily topHat} power spectrum window function class.
    implicit none
    type (powerSpectrumWindowFunctionTopHat)                        :: topHatConstructorInternal
    class(cosmologyParametersClass         ), target, intent(in   ) :: cosmologyParameters_    

    topHatConstructorInternal%cosmologyParameters_ => cosmologyParameters_
    return
  end function topHatConstructorInternal

  subroutine topHatDestructor(self)
    !% Destructor for the {\normalfont \ttfamily topHat} power spectrum window function class.
    implicit none
    type(powerSpectrumWindowFunctionTopHat), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    return
  end subroutine topHatDestructor

  double precision function topHatValue(self,wavenumber,smoothingMass)
    !% Top hat in real space window function Fourier transformed into $k$-space used in computing the variance of the power
    !% spectrum.
    use Numerical_Constants_Math
    use Cosmology_Parameters
    implicit none
    class           (powerSpectrumWindowFunctionTopHat), intent(inout) :: self
    double precision                                   , intent(in   ) :: smoothingMass        , wavenumber
    double precision                                   , parameter     :: xSeriesMaximum=1.0d-3
    double precision                                                   :: topHatRadius         , x         , &
         &                                                                xSquared

    topHatRadius=+(                                             &
         &         +3.0d0                                       &
         &         /4.0d0                                       &
         &         /Pi                                          &
         &         *smoothingMass                               &
         &         /self%cosmologyParameters_%OmegaMatter    () &
         &         /self%cosmologyParameters_%densityCritical() &
         &        )**(1.0d0/3.0d0)
    x           =+wavenumber                                    &
         &       *topHatRadius
    if      (x <= 0.0d0         ) then
       topHatValue=+0.0d0
    else if (x <= xSeriesMaximum) then
       ! Use a series expansion of the window function for small x.
       xSquared   =+x**2
       topHatValue=+1.0d0                        &
            &      +xSquared*(  -1.0d0/   10.0d0 &
            &      +xSquared* ( +1.0d0/  280.0d0 &
            &      +xSquared*  (-1.0d0/15120.0d0 &
            &                  )                 &
            &                 )                  &
            &                )
    else
       ! For larger x, use the full expression.
       topHatValue=3.0d0*(sin(x)-x*cos(x))/(x**3)
    end if
    return
  end function topHatValue

  double precision function topHatWavenumberMaximum(self,smoothingMass)
    !% Maximum wavenumber for a top hat in real space window function Fourier transformed into $k$-space used in computing the
    !% variance of the power spectrum.
    implicit none
    class           (powerSpectrumWindowFunctionTopHat), intent(inout) :: self
    double precision                                   , intent(in   ) :: smoothingMass
    double precision                                   , parameter     :: wavenumberLarge=1.0d30 !   Effective infinity.

    topHatWavenumberMaximum=wavenumberLarge
    return
  end function topHatWavenumberMaximum

  subroutine topHatDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class(powerSpectrumWindowFunctionTopHat), intent(inout) :: self
    type (inputParameters                  ), intent(inout) :: descriptor
    type (inputParameters                  )                :: subParameters

    call descriptor%addParameter("powerSpectrumWindowFunctionMethod","topHat")
    subParameters=descriptor%subparameters("powerSpectrumWindowFunctionMethod")
    call self%cosmologyParameters_%descriptor(subParameters)
    return
  end subroutine topHatDescriptor
