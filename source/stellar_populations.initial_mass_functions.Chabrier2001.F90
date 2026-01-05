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

  !!{
  Implements a stellar initial mass function class based on \cite{chabrier_galactic_2001}.
  !!}

  !![
  <initialMassFunction name="initialMassFunctionChabrier2001">
   <description>
    A stellar initial mass function class based on \cite{chabrier_galactic_2001}:
    \begin{equation}
     \phi(M) \propto \left\{ \begin{array}{ll}
     M^{-1} \exp(-[\log_{10}(M/M_\mathrm{c})/\sigma_\mathrm{c}]^2/2) &amp; \hbox{ for } M_\mathrm{l} &lt; M &lt; M_\mathrm{t}  \\
     M^\alpha &amp; \hbox{ for } M_\mathrm{t} &lt; M &lt; M_\mathrm{u} \\
     0 &amp; \hbox {otherwise,} \end{array} \right.
    \end{equation}
    where $\sigma_\mathrm{c}=${\normalfont \ttfamily [sigma]}, $M_\mathrm{c}=${\normalfont \ttfamily
    [massCharacteristic]}$M_\odot$, $\alpha=${\normalfont \ttfamily [exponent]}, $M_\mathrm{t}=${\normalfont \ttfamily
    [massTransition]}$M_\odot$, $M_\mathrm{l}=${\normalfont \ttfamily [massLower]}$M_\odot$, and $M_\mathrm{u}=${\normalfont
    \ttfamily [massUpper]}$M_\odot$.
   </description>
  </initialMassFunction>
  !!]
  type, extends(initialMassFunctionClass) :: initialMassFunctionChabrier2001
     !!{
     A stellar initial mass function class based on \cite{chabrier_galactic_2001}.
     !!}
     private
     double precision :: massLower               , massTransition        , &
          &              massUpper               , exponent              , &
          &              massCharacteristic      , sigma                 , &
          &              normalizationExponential, normalizationLogNormal
   contains
     procedure :: massMinimum      => chabrier2001MassMinimum
     procedure :: massMaximum      => chabrier2001MassMaximum
     procedure :: phi              => chabrier2001Phi
     procedure :: numberCumulative => chabrier2001NumberCumulative
     procedure :: tabulate         => chabrier2001Tabulate
     procedure :: label            => chabrier2001Label
  end type initialMassFunctionChabrier2001

  interface initialMassFunctionChabrier2001
     !!{
     Constructors for the \refClass{initialMassFunctionChabrier2001} initial mass function class.
     !!}
     module procedure chabrier2001ConstructorParameters
     module procedure chabrier2001ConstructorInternal
  end interface initialMassFunctionChabrier2001

contains

  function chabrier2001ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{initialMassFunctionChabrier2001} initial mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (initialMassFunctionChabrier2001)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: massLower         , massTransition, &
          &                                                             massUpper         , exponent      , &
          &                                                             massCharacteristic, sigma

    !![
    <inputParameter>
      <name>massUpper</name>
      <defaultValue>125.0d0</defaultValue>
      <description>The upper mass limit for the \cite{chabrier_galactic_2001} \gls{imf}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massLower</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The lower mass limit for the \cite{chabrier_galactic_2001} \gls{imf}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massTransition</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The transition limit for the \cite{chabrier_galactic_2001} \gls{imf}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sigma</name>
      <defaultValue>0.69d0</defaultValue>
      <description>The width of the lognormal part of the \cite{chabrier_galactic_2001} \gls{imf}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <defaultValue>-2.3d0</defaultValue>
      <description>The exponent of the power law part of the \cite{chabrier_galactic_2001} \gls{imf}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massCharacteristic</name>
      <defaultValue>0.08d0</defaultValue>
      <description>Characteristic mass of the lognormal part of the \cite{chabrier_galactic_2001} \gls{imf}.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=initialMassFunctionChabrier2001(massLower,massTransition,massUpper,exponent,massCharacteristic,sigma)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function chabrier2001ConstructorParameters

  function chabrier2001ConstructorInternal(massLower,massTransition,massUpper,exponent,massCharacteristic,sigma) result(self)
    !!{
    Internal constructor for the \refClass{initialMassFunctionChabrier2001} initial mass function.
    !!}
    use :: Error_Functions         , only : Error_Function
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (initialMassFunctionChabrier2001)                :: self
    double precision                                 , intent(in   ) :: massLower         , massTransition, &
         &                                                              massUpper         , exponent      , &
         &                                                              massCharacteristic, sigma
    double precision                                                 :: normalization
    !![
    <constructorAssign variables="massLower,massTransition,massUpper,exponent,massCharacteristic,sigma"/>
    !!]

    self%normalizationLogNormal  =+sqrt(Pi/2.0d0)                                     &
         &                        *self%sigma                                         &
         &                        *self%massCharacteristic                            &
         &                        *log(10.0d0)                                        &
         &                        *exp(                                               &
         &                             +0.5d0                                         &
         &                             *self%sigma**2                                 &
         &                             *log(10.0d0)     **2                           &
         &                            )                                               &
         &                        *(                                                  &
         &                          -Error_Function(                                  &
         &                                          +self%sigma                       &
         &                                          *log (10.0d0)                     &
         &                                          /sqrt( 2.0d0)                     &
         &                                          -log10(                           &
         &                                                 +self%massTransition       &
         &                                                 /self%massCharacteristic   &
         &                                                )                           &
         &                                          /sqrt( 2.0d0)                     &
         &                                          /self%sigma                       &
         &                                         )                                  &
         &                          +Error_Function(                                  &
         &                                          +self%sigma                       &
         &                                          *log (10.0d0)                     &
         &                                          /sqrt( 2.0d0)                     &
         &                                          -log10(                           &
         &                                                 +self%massLower            &
         &                                                 /self%massCharacteristic   &
         &                                                )                           &
         &                                          /sqrt( 2.0d0)                     &
         &                                          /self%sigma                       &
         &                                         )                                  &
         &                         )
    self%normalizationExponential=+exp(                                               &
         &                              -0.50d0                                       &
         &                              *log10(                                       &
         &                                     +self%massTransition                   &
         &                                     /self%massCharacteristic               &
         &                                    )                        ** 2               &
         &                              /self%sigma                    ** 2               &
         &                             )                                                  &
         &                         *self%massTransition                                   &
         &                         /                                     (2.0d0+exponent) &
         &                         *(                                                     &
         &                           +(                                                   &
         &                             +self%massUpper                                    &
         &                             /self%massTransition                               &
         &                            )                                **(2.0d0+exponent) &
         &                           -1.0d0                                               &
         &                          )
    normalization                =+self%normalizationLogNormal   &
         &                        +self%normalizationExponential
    self%normalizationLogNormal  =1.0d0/normalization
    self%normalizationExponential=+exp(                                    &
         &                              -0.50d0                            &
         &                              *log10(                            &
         &                                     +self%massTransition        &
         &                                     /self%massCharacteristic    &
         &                                    )                        **2 &
         &                              /self%sigma                    **2 &
         &                             )                                   &
         &                         /self%massTransition                    &
         &                         /normalization
    return
  end function chabrier2001ConstructorInternal

  double precision function chabrier2001MassMinimum(self)
    !!{
    Return the minimum mass of stars in the \cite{chabrier_galactic_2001} \gls{imf}.
    !!}
    implicit none
    class(initialMassFunctionChabrier2001), intent(inout) :: self

    chabrier2001MassMinimum=self%massLower
    return
  end function chabrier2001MassMinimum

  double precision function chabrier2001MassMaximum(self)
    !!{
    Return the maximum mass of stars in the \cite{chabrier_galactic_2001} \gls{imf}.
    !!}
    implicit none
    class(initialMassFunctionChabrier2001), intent(inout) :: self

    chabrier2001MassMaximum=self%massUpper
    return
  end function chabrier2001MassMaximum

  double precision function chabrier2001Phi(self,massInitial)
    !!{
    Evaluate the \cite{chabrier_galactic_2001} stellar initial mass function.
    !!}
    implicit none
    class           (initialMassFunctionChabrier2001), intent(inout) :: self
    double precision                                 , intent(in   ) :: massInitial

    if      (                                    &
         &    massInitial >= self%massLower      &
         &   .and.                               &
         &    massInitial <  self%massTransition &
         &  ) then
       chabrier2001Phi=+self%normalizationLogNormal           &
            &          *exp(                                  &
            &               -0.5d0                            &
            &               *(                                &
            &                 +log10(                         &
            &                        +     massInitial        &
            &                        /self%massCharacteristic &
            &                       )                         &
            &                 /self%sigma                     &
            &                )**2                             &
            &              )                                  &
            &          /massInitial
    else if (                                    &
         &    massInitial >= self%massTransition &
         &   .and.                               &
         &    massInitial <  self%massUpper      &
         &  ) then
       chabrier2001Phi=+self%normalizationExponential &
            &          *massInitial**self%exponent
    else
       chabrier2001Phi=0.0d0
    end if
    return
  end function chabrier2001Phi

  double precision function chabrier2001NumberCumulative(self,massLower,massUpper) result(number)
    !!{
    Evaluate a piecewise power-law stellar initial mass function.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (initialMassFunctionChabrier2001), intent(inout) :: self
    double precision                                 , intent(in   ) :: massLower , massUpper
    double precision                                                 :: massLower_, massUpper_

    number=0.0d0
    ! Gaussian piece.
    massLower_=max(massLower,self%massLower     )
    massUpper_=min(massUpper,self%massTransition)
    if (massUpper_ > massLower_)                                                           &
         & number=+number                                                                  &
         &        +sqrt(Pi/ 2.0d0)                                                         &
         &        *log (   10.0d0)                                                         &
         &        *self%sigma                                                              &
         &        *self%normalizationLogNormal                                             &
         &        *(                                                                       &
         &          +erf(log10(massUpper_/self%massCharacteristic)/sqrt(2.0d0)/self%sigma) &
         &          -erf(log10(massLower_/self%massCharacteristic)/sqrt(2.0d0)/self%sigma) &
         &        )
    ! Power-law piece.
    massLower_=max(massLower,self%massTransition)
    massUpper_=min(massUpper,self%massUpper     )
    if (massUpper_ > massLower_)                       &
         & number=+     number                         &
         &        +self%normalizationExponential       &
         &        /              (1.0d0+self%exponent) &
         &        *(                                   &
         &          +massUpper_**(1.0d0+self%exponent) &
         &          +massLower_**(1.0d0+self%exponent) &
         &         )
    return
  end function chabrier2001NumberCumulative

  subroutine chabrier2001Tabulate(self,imfTable)
    !!{
    Construct and return a tabulation of the \cite{chabrier_galactic_2001} \gls{imf}.
    !!}
    use :: Tables, only : table1DLogarithmicLinear
    implicit none
    class  (initialMassFunctionChabrier2001)             , intent(inout) :: self
    class  (table1D                        ), allocatable, intent(inout) :: imfTable
    integer                                 , parameter                  :: countTable=100
    integer                                                              :: i

    allocate(table1DLogarithmicLinear :: imfTable)
    select type (imfTable)
    type is (table1DLogarithmicLinear)
       call imfTable%create(                 &
            &               self%massLower , &
            &               self%massUpper , &
            &                    countTable  &
            &              )
       do i=1,countTable
          call imfTable%populate(self%phi(imfTable%x(i)),i)
       end do
    end select
    return
  end subroutine chabrier2001Tabulate

  function chabrier2001Label(self)
    !!{
    Return a label for this \gls{imf}.
    !!}
    implicit none
    class(initialMassFunctionChabrier2001), intent(inout) :: self
    type (varying_string                    )             :: chabrier2001Label
    !$GLC attributes unused :: self

    chabrier2001Label="Chabrier2001"
    return
  end function chabrier2001Label
