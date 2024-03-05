!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implements a generalization of the ETHOS power spectrum window function class from \cite{bohr_halo_2021}.
  !!}
  
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredClass

  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionETHOSExtended">
   <description>
     A generalization of the ETHOS window function for filtering of power spectra from \cite{bohr_halo_2021}. The window function
     has the same functional form as in \cite{bohr_halo_2021}, but the parameters $c_\mathrm{W}$ and $\beta$ are now scale
     dependent following     
     \begin{equation}
     x = x_0 x_1^{n-n_0}
     \end{equation}     
     where $x$ refers to either $c_\mathrm{W}$ or $\beta$, $n = \mathrm{d}\log P / \mathrm{d} \log k$ is the logarithmic
     derivative of the linear theory power spectrum, and $n_0 = -2.6$ is a convenient zero-point.
   </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionETHOS) :: powerSpectrumWindowFunctionETHOSExtended
     !!{
     A generalization of the ETHOS power spectrum window function class.
     !!}
     private
     class           (powerSpectrumPrimordialTransferredClass), pointer :: powerSpectrumPrimordialTransferred_ => null()
     double precision                                                   :: cW0                                          , beta0, &
          &                                                                cW1                                          , beta1
   contains
     final     ::         ETHOSExtendedDestructor
     procedure :: cW   => ETHOSExtendedCW
     procedure :: beta => ETHOSExtendedBeta
  end type powerSpectrumWindowFunctionETHOSExtended

  interface powerSpectrumWindowFunctionETHOSExtended
     !!{
     Constructors for the ETHOS power spectrum window function class.
     !!}
     module procedure ETHOSExtendedConstructorParameters
     module procedure ETHOSExtendedConstructorInternal
  end interface powerSpectrumWindowFunctionETHOSExtended

  ! Zero-point in logarithmic slope of the power spectrum. Arbitrary, but convenient.
  double precision, parameter :: logarithmicDerivativeReference=-2.6d0

  ! Maximum allowed value for parameters.
  double precision, parameter :: parameterValueMaximum         =huge(0.0d0)/1.0d30
  double precision, parameter :: exponentPowerMaximum          =1.0d2

contains

  function ETHOSExtendedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ETHOS  power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionETHOSExtended)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (powerSpectrumPrimordialTransferredClass ), pointer       :: powerSpectrumPrimordialTransferred_
    class           (cosmologyParametersClass                ), pointer       :: cosmologyParameters_
    double precision                                                          :: cW0                                , beta0, &
         &                                                                       cW1                                , beta1

    !![
    <inputParameter>
      <name>cW0</name>
      <source>parameters</source>
      <defaultValue>3.78062835d0</defaultValue>
      <description>The parameter $c_\mathrm{W,0}$ in the generalized ETHOS power spectrum window function.</description>
    </inputParameter>
    <inputParameter>
      <name>beta0</name>
      <source>parameters</source>
      <defaultValue>3.4638743d0</defaultValue>
      <description>The parameter $\beta_0$ in the generalized ETHOS power spectrum window function.</description>
    </inputParameter>
    <inputParameter>
      <name>cW1</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $c_\mathrm{W,1}$ in the generalized ETHOS power spectrum window function.</description>
    </inputParameter>
    <inputParameter>
      <name>beta1</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $\beta_1$ in the generalized ETHOS power spectrum window function.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"/>
    <objectBuilder class="powerSpectrumPrimordialTransferred" name="powerSpectrumPrimordialTransferred_" source="parameters"/>
    !!]
    self=powerSpectrumWindowFunctionETHOSExtended(cW0,cW1,beta0,beta1,cosmologyParameters_,powerSpectrumPrimordialTransferred_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"               />
    <objectDestructor name="powerSpectrumPrimordialTransferred_"/>
    !!]
    return
  end function ETHOSExtendedConstructorParameters

  function ETHOSExtendedConstructorInternal(cW0,cW1,beta0,beta1,cosmologyParameters_,powerSpectrumPrimordialTransferred_) result(self)
    !!{
    Internal constructor for the ETHOS power spectrum window function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (powerSpectrumWindowFunctionETHOSExtended)                        :: self
    double precision                                          , intent(in   )         :: cW0                                , beta0, &
         &                                                                               cW1                                , beta1
    class           (cosmologyParametersClass                ), intent(in   ), target :: cosmologyParameters_
    class           (powerSpectrumPrimordialTransferredClass ), intent(in   ), target :: powerSpectrumPrimordialTransferred_
    !![
    <constructorAssign variables="cW0, cW1, beta0, beta1, *cosmologyParameters_, *powerSpectrumPrimordialTransferred_"/>
    !!]

    return
  end function ETHOSExtendedConstructorInternal

  subroutine ETHOSExtendedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily ETHOS} window function class.
    !!}
    implicit none
    type   (powerSpectrumWindowFunctionETHOSExtended), intent(inout) :: self

    !![
    <objectDestructor name="self%powerSpectrumPrimordialTransferred_"/>
    !!]
    return
  end subroutine ETHOSExtendedDestructor
  
  double precision function ETHOSExtendedCW(self,wavenumber,time) result(cW)
    !!{
    Compute the $c_\mathrm{W}$ parameter for the extended ETHOS window function.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOSExtended), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber   , time
    double precision                                                          :: exponentPower

    exponentPower=min(                                                                                                                     &
         &            max(                                                                                                                 &
         &                +self%powerSpectrumPrimordialTransferred_%logarithmicDerivative(wavenumber,time)-logarithmicDerivativeReference, &
         &                -exponentPowerMaximum                                                                                            &
         &               )                                                                                                               , &
         &                +exponentPowerMaximum                                                                                            &
         &           )
    if (exponent(self%cW0)+exponent(self%cW1)*exponentPower < exponent(parameterValueMaximum)) then
       cW     =+self%cW0                &
            &  *self%cW1**exponentPower
    else
       cW     =+parameterValueMaximum
    end if
    return
  end function ETHOSExtendedCW
  
  double precision function ETHOSExtendedBeta(self,wavenumber,time) result(beta)
    !!{
    Compute the $\beta$ parameter for the extended ETHOS window function.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOSExtended), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber   , time
    double precision                                                          :: exponentPower

    exponentPower=min(                                                                                                                     &
         &            max(                                                                                                                 &
         &                +self%powerSpectrumPrimordialTransferred_%logarithmicDerivative(wavenumber,time)-logarithmicDerivativeReference, &
         &                -exponentPowerMaximum                                                                                            &
         &               )                                                                                                               , &
         &                +exponentPowerMaximum                                                                                            &
         &           )
    if (exponent(self%beta0)+exponent(self%beta1)*exponentPower < maxExponent(parameterValueMaximum)) then
       beta   =+self%beta0                &
            &  *self%beta1**exponentPower
    else
       beta   =+parameterValueMaximum
    end if
    return
  end function ETHOSExtendedBeta
  

