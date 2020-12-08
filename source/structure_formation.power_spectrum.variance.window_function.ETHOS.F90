!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements the ETHOS power spectrum window function class.
  
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !# <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionETHOS">
  !#  <description>
  !#   ETHOS window function for filtering of power spectra. This window function was chosen to give good matches to N-body halo
  !#   mass functions derived from the ETHOS transfer functions. The functional form was provided by Francis-Yan Cyr-Racine
  !#   ({\normalfont \ttfamily \textless fycr@unm.edu\textgreater}; private communication):
  !#   \begin{equation}
  !#    W(kR) = (\frac{1}{1+\left(\frac{kR}{c}\right)^\beta})
  !#   \end{equation}
  !#   with $c = 3.78062835$, $\beta = 3.4638743$, where $R$ is related to $M$ via the standard relation $M =
  !#   \frac{4\pi}{3}\bar\rho_m R^3$.
  !#  </description>
  !# </powerSpectrumWindowFunction>
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionETHOS
     !% ETHOS power spectrum window function class.
     private
     class(cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
   contains
     final     ::                      ETHOSDestructor
     procedure :: value             => ETHOSValue
     procedure :: wavenumberMaximum => ETHOSWavenumberMaximum
  end type powerSpectrumWindowFunctionETHOS

  interface powerSpectrumWindowFunctionETHOS
     !% Constructors for the ETHOS power spectrum window function class.
     module procedure ETHOSConstructorParameters
     module procedure ETHOSConstructorInternal
  end interface powerSpectrumWindowFunctionETHOS

contains

  function ETHOSConstructorParameters(parameters) result(self)
    !% Constructor for the ETHOS  power spectrum window function class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (powerSpectrumWindowFunctionETHOS)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    self=ETHOSConstructorInternal(cosmologyParameters_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    return
  end function ETHOSConstructorParameters

  function ETHOSConstructorInternal(cosmologyParameters_) result(self)
    !% Internal constructor for the ETHOS power spectrum window function class.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type (powerSpectrumWindowFunctionETHOS)                        :: self
    class(cosmologyParametersClass        ), target, intent(in   ) :: cosmologyParameters_
    !# <constructorAssign variables="*cosmologyParameters_"/>

    return
  end function ETHOSConstructorInternal

  subroutine ETHOSDestructor(self)
    !% Destructor for the ETHOS power spectrum window function class.
    implicit none
    type(powerSpectrumWindowFunctionETHOS), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    return
  end subroutine ETHOSDestructor

  double precision function ETHOSValue(self,wavenumber,smoothingMass)
    !% ETHOS window function used in computing the variance of the power spectrum.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumWindowFunctionETHOS), intent(inout) :: self
    double precision                                  , intent(in   ) :: smoothingMass               , wavenumber
    double precision                                  , parameter     :: b               =3.4638743d0, c         =3.78062835d0 
    double precision                                                  :: radius
    
    radius =+(                                             &
         &    +3.0d0                                       &
         &    /4.0d0                                       &
         &    /Pi                                          &
         &    *smoothingMass                               &
         &    /self%cosmologyParameters_%OmegaMatter    () &
         &    /self%cosmologyParameters_%densityCritical() &
         &   )**(1.0d0/3.0d0)
    if (wavenumber <= 0.0d0) then
       ETHOSValue=+0.0d0
    else
       ETHOSValue=+1.0d0          &
            &     /(              &
            &       +1.0d0        &
            &       +(            &
            &         +wavenumber &
            &         *radius     &
            &         /c          &
            &        )**b         &
            &      )
    end if
    return
  end function ETHOSValue

  double precision function ETHOSWavenumberMaximum(self,smoothingMass)
    !% Maximum wavenumber for a top hat in real space window function Fourier transformed into $k$-space used in computing the
    !% variance of the power spectrum.
    implicit none
    class           (powerSpectrumWindowFunctionETHOS), intent(inout) :: self
    double precision                                  , intent(in   ) :: smoothingMass
    double precision                                  , parameter     :: wavenumberLarge=huge(1.0d0) ! Effective infinity.
    !$GLC attributes unused :: self, smoothingMass

    ETHOSWavenumberMaximum=wavenumberLarge
    return
  end function ETHOSWavenumberMaximum


