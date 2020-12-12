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

!% Contains a module which implements the ETHOS power spectrum window function class.
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !# <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionETHOS">
  !#  <description>ETHOS window function for filtering of power spectra.</description>
  !# </powerSpectrumWindowFunction>
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionETHOS
     !% ETHOS power spectrum window function class.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: cutOffNormalization
     type            (varying_string          )          :: normalization
   contains
     final     ::                               ETHOSDestructor
     procedure :: value                      => ETHOSValue
     procedure :: wavenumberMaximum          => ETHOSWavenumberMaximum
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
    type            (powerSpectrumWindowFunctionETHOS)                      :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyParametersClass              ), pointer       :: cosmologyParameters_
    type            (varying_string                        )                :: normalization
    character       (len=32                                )                :: normalizationChar
    double precision                                                        :: normalizationValue

    ! Check parameters.
    !# <inputParameter>
    !#   <name>normalization</name>
    !#   <source>parameters</source>
    !#   <variable>normalization</variable>
    !#   <defaultValue>var_str('natural')</defaultValue>
    !#   <description>
    !#     The parameter $a$ in the relation $k_\mathrm{s} = a/r_\mathrm{s}$, where $k_\mathrm{s}$ is the cut-off wavenumber for
    !#     the sharp $k$-space window function and $r_\mathrm{s}$ is the radius of a sphere (in real-space) enclosing the
    !#     requested smoothing mass. Alternatively, a value of {\normalfont \ttfamily natural} will be supplied in which case the normalization
    !#     is chosen such that, in real-space, $W(r=0)=1$. This results in a contained mass of $M=6 \pi^2 \bar{\rho} k_\mathrm{s}^{-3}$.
    !#   </description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    if (normalization == 'natural') then
       normalizationValue=0.0d0
    else
       normalizationChar=char(normalization)
       read (normalizationChar,*) normalizationValue
    end if
    self=ETHOSConstructorInternal(cosmologyParameters_,normalizationValue)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    return
  end function ETHOSConstructorParameters

  function ETHOSConstructorInternal(cosmologyParameters_,normalization) result(self)
    !% Internal constructor for the ETHOS power spectrum window function class.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (powerSpectrumWindowFunctionETHOS)                              :: self
    class           (cosmologyParametersClass              ), target, intent(in   ) :: cosmologyParameters_
    double precision                                                , intent(in   ) :: normalization
    character       (len=18                                )                        :: normalizationText
    !# <constructorAssign variables="*cosmologyParameters_"/>

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
    implicit none
    class           (powerSpectrumWindowFunctionETHOS), intent(inout) :: self
    double precision                                        , intent(in   ) :: smoothingMass   , wavenumber
    double precision                                                        :: wavenumberCutOff, b, c, ETHOS_Radius

    b = 3.4638743d0
    c = 3.78062835d0 

    ETHOS_Radius = +( +3.0d0                                    &
         &         /4.0d0                                       &
         &         /3.14159265d0                                &
         &         *smoothingMass                               &
         &         /self%cosmologyParameters_%OmegaMatter    () &
         &         /self%cosmologyParameters_%densityCritical() &
         &        )**(1.0d0/3.0d0)

    wavenumberCutOff=self%wavenumberMaximum(smoothingMass)
    if      (wavenumber <=            0.0d0) then
       ETHOSValue=0.0d0
    else
       ETHOSValue= 1.0d0/(1.0d0+(wavenumber*ETHOS_Radius/c)**b)
    end if
    return
  end function ETHOSValue

   double precision function ETHOSWavenumberMaximum(self,smoothingMass)
    !% Sets maximum wavenumber to effectively infinity (really large number)
    implicit none
    class           (powerSpectrumWindowFunctionETHOS), intent(inout) :: self
    double precision                                   , intent(in   ) :: smoothingMass
    double precision                                   , parameter     :: wavenumberLarge=huge(1.0d0) ! Effective infinity.
    !$GLC attributes unused :: self, smoothingMass

    ETHOSWavenumberMaximum=wavenumberLarge
    return
  end function ETHOSWavenumberMaximum


