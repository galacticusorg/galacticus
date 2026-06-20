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

  !!{RST
  Implements the ETHOS power spectrum window function class from :cite:t:`bohr_halo_2021`.
  !!}
  
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionETHOS" docformat="rst">
   <description>
   ETHOS window function for filtering of power spectra from :cite:t:`bohr_halo_2021`. This window function was chosen to give good matches to N-body halo mass functions derived from the ETHOS transfer functions. Specifically the window function is given by:

   .. math::

      W(kR) = \frac{1}{1+\left(\frac{kR}{c_\mathrm{W}}\right)^\beta}

   with defaults of :math:`c_\mathrm{W} = 3.78062835`, :math:`\beta = 3.4638743`, where :math:`R` is related to :math:`M` via the standard relation :math:`M = \frac{4\pi}{3}\bar\rho_m R^3`.
   </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionETHOS
     !!{RST
     ETHOS power spectrum window function class.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: cW_                           , beta_
   contains
     !![
     <methods docformat="rst">
       <method method="cW"                      description="Compute the parameter :math:`c_\mathrm{W}` in the ETHOS window function."  />
       <method method="beta"                    description="Compute the parameter :math:`\beta` in the ETHOS window function."         />
       <method method="wavenumberScaledMinimum" description="Compute the parameter :math:`x_\mathrm{min}` in the ETHOS window function."/>
     </methods>
     !!]
     final     ::                            ETHOSDestructor
     procedure :: value                   => ETHOSValue
     procedure :: wavenumberMaximum       => ETHOSWavenumberMaximum
     procedure :: cW                      => ETHOSCW
     procedure :: beta                    => ETHOSBeta
     procedure :: wavenumberScaledMinimum => ETHOSWavenumberScaledMinimum
  end type powerSpectrumWindowFunctionETHOS

  interface powerSpectrumWindowFunctionETHOS
     !!{RST
     Constructors for the ETHOS power spectrum window function class.
     !!}
     module procedure ETHOSConstructorParameters
     module procedure ETHOSConstructorInternal
  end interface powerSpectrumWindowFunctionETHOS

contains

  function ETHOSConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ETHOS  power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionETHOS)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    double precision                                                  :: cW                  , beta

    !![
    <inputParameter docformat="rst">
      <name>cW</name>
      <source>parameters</source>
      <defaultValue>3.78062835d0</defaultValue>
      <description>
      The parameter :math:`c_\mathrm{W}` in the :cite:t:`bohr_halo_2021` power spectrum window function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>3.4638743d0</defaultValue>
      <description>
      The parameter :math:`\beta` in the :cite:t:`bohr_halo_2021` power spectrum window function.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=powerSpectrumWindowFunctionETHOS(cW,beta,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function ETHOSConstructorParameters

  function ETHOSConstructorInternal(cW_,beta_,cosmologyParameters_) result(self)
    !!{RST
    Internal constructor for the ETHOS power spectrum window function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (powerSpectrumWindowFunctionETHOS)                        :: self
    double precision                                          , intent(in   ) :: cW_                 , beta_
    class           (cosmologyParametersClass        ), target, intent(in   ) :: cosmologyParameters_
    !![
    <constructorAssign variables="cW_, beta_, *cosmologyParameters_"/>
    !!]

    return
  end function ETHOSConstructorInternal

  subroutine ETHOSDestructor(self)
    !!{RST
    Destructor for the ETHOS power spectrum window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionETHOS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine ETHOSDestructor

  double precision function ETHOSValue(self,wavenumber,smoothingMass,time)
    !!{RST
    ETHOS window function used in computing the variance of the power spectrum. Best fit values for parameters are from :cite:t:`bohr_halo_2021`.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumWindowFunctionETHOS), intent(inout) :: self
    double precision                                  , intent(in   ) :: smoothingMass          , wavenumber      , &
         &                                                               time
    double precision                                                  :: radius                 , wavenumberScaled, &
         &                                                               cW                     , beta            , &
         &                                                               wavenumberScaledMinimum
    
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
       cW  =self%cW  (wavenumber,time)
       beta=self%beta(wavenumber,time)
       if (cW > 0.0d0 .and. exponent(wavenumber)+exponent(radius)-exponent(cw) < maxExponent(0.0d0)) then
          wavenumberScaled       =+wavenumber &
               &                  *radius     &
               &                  /cW
          wavenumberScaledMinimum=+self%wavenumberScaledMinimum(wavenumber,time)
          if (wavenumberScaled <= wavenumberScaledMinimum) then
             ETHOSValue=1.0d0
          else
             if (beta*exponent(wavenumberScaled) < maxExponent(0.0d0)) then
                ETHOSValue=+1.0d0                                              &
                     &     /(                                                  &
                     &       +1.0d0                                            &
                     &       +(wavenumberScaled-wavenumberScaledMinimum)**beta &
                     &      )
             else
                ETHOSvalue=0.0d0
             end if
          end if
       else
          ETHOSvalue=0.0d0
       end if
    end if
    return
  end function ETHOSValue

  double precision function ETHOSWavenumberMaximum(self,smoothingMass)
    !!{RST
    Sets maximum wavenumber to effectively infinity (really large number).
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOS), intent(inout) :: self
    double precision                                  , intent(in   ) :: smoothingMass
    double precision                                  , parameter     :: wavenumberLarge=huge(1.0d0) ! Effective infinity.
    !$GLC attributes unused :: self, smoothingMass

    ETHOSWavenumberMaximum=wavenumberLarge
    return
  end function ETHOSWavenumberMaximum

  double precision function ETHOSCW(self,wavenumber,time) result(cW)
    !!{RST
    Compute the :math:`c_\mathrm{W}` parameter for the ETHOS window function.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOS), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber, time
    !$GLC attributes unused :: wavenumber, time

    cW=self%cW_
    return
  end function ETHOSCW
  
  double precision function ETHOSBeta(self,wavenumber,time) result(beta)
    !!{RST
    Compute the :math:`\beta` parameter for the ETHOS window function.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOS), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber, time
    !$GLC attributes unused :: wavenumber, time

    beta=self%beta_
    return
  end function ETHOSBeta
  
  double precision function ETHOSWavenumberScaledMinimum(self,wavenumber,time) result(wavenumberScaledMinimum)
    !!{RST
    Compute the :math:`\beta` parameter for the ETHOS window function.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOS), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber, time
    !$GLC attributes unused :: self, wavenumber, time

    wavenumberScaledMinimum=0.0d0
    return
  end function ETHOSWavenumberScaledMinimum
  

