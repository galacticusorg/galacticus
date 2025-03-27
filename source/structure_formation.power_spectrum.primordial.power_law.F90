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
  A primordial power spectrum class which provides a power-law power spectrum.
  !!}

  !![
  <powerSpectrumPrimordial name="powerSpectrumPrimordialPowerLaw">
   <description>
    Implements a power-law primordial power spectrum, possibly with a running index. The primordial power spectrum has the
    form:
    \begin{equation}
     P(k) \propto k^{n_\mathrm{eff}(k)},
    \end{equation}
    where
    \begin{equation}
     n_\mathrm{eff}(k) = n_\mathrm{s} + {1\over 2}{\d n \over \d \ln k} \ln \left( {k \over k_\mathrm{ref}} \right) + {1\over 6}{\d^2 n \over \d \ln k^2} \left[ \ln \left( {k \over k_\mathrm{ref}} \right) \right]^2,
    \end{equation}
    where $n_\mathrm{s}=${\normalfont \ttfamily [index]} is the power spectrum index at wavenumber
    $k_\mathrm{ref}=${\normalfont \ttfamily [wavenumberReference]}, $\d n / \d \ln k=${\normalfont \ttfamily [running]}, and $\d^2 n / \d \ln k^2=${\normalfont \ttfamily [runningRunning]}
    describes the running of this index with wavenumber.
   </description>
  </powerSpectrumPrimordial>
  !!]
  type, extends(powerSpectrumPrimordialClass) :: powerSpectrumPrimordialPowerLaw
     !!{
     A power-law primordial power spectrum class.
     !!}
     private
     double precision :: index_                , running            , &
          &              runningRunning        , wavenumberReference
     logical          :: runningSmallScalesOnly
   contains
     procedure :: power                 => powerLawPower
     procedure :: logarithmicDerivative => powerLawLogarithmicDerivative
  end type powerSpectrumPrimordialPowerLaw

  interface powerSpectrumPrimordialPowerLaw
     !!{
     Constructors for the {\normalfont \ttfamily powerLaw} primordial power spectrum class.
     !!}
     module procedure powerLawConstructorParameters
     module procedure powerLawConstructorInternal
  end interface powerSpectrumPrimordialPowerLaw

contains

  function powerLawConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily powerLaw} primordial power spectrum class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumPrimordialPowerLaw)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: index_                , wavenumberReference, &
         &                                                              running               , runningRunning
    logical                                                          :: runningSmallScalesOnly

    !![
    <inputParameter>
      <name>index</name>
      <variable>index_</variable>
      <source>parameters</source>
      <defaultValue>0.9649d0</defaultValue>
      <defaultSource>(\citealt{planck_collaboration_planck_2018}; TT,TE,EE$+$lowE$+$lensing)</defaultSource>
      <description>The index of the power-law primordial power spectrum.</description>
    </inputParameter>
    <inputParameter>
      <name>running</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The running, $\d n_\mathrm{s} / \d \ln k$, of the power spectrum index.</description>
    </inputParameter>
    <inputParameter>
      <name>runningRunning</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The running-of-the-running, $\d^2 n_\mathrm{s} / \d \ln k^2$, of the power spectrum index.</description>
    </inputParameter>
    <inputParameter>
      <name>wavenumberReference</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>When a running power spectrum index is used, this is the wavenumber, $k_\mathrm{ref}$, at which the index is equal to {\normalfont \ttfamily [index]}.</description>
    </inputParameter>
    <inputParameter>
      <name>runningSmallScalesOnly</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If {\normalfont \ttfamily true} then the index runs only for $k > k_\mathrm{ref}$, for smaller $k$ the index is constant.</description>
    </inputParameter>
    !!]
    self=powerSpectrumPrimordialPowerLaw(index_,running,runningRunning,wavenumberReference,runningSmallScalesOnly)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawConstructorParameters

  function powerLawConstructorInternal(index_,running,runningRunning,wavenumberReference,runningSmallScalesOnly) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily powerLaw} primordial power spectrum class.
    !!}
    use :: Error  , only : Warn
    use :: Display, only : displayBlue, displayYellow, displayGreen, displayReset
    implicit none
    type            (powerSpectrumPrimordialPowerLaw)                :: self
    double precision                                 , intent(in   ) :: index_                , wavenumberReference, &
         &                                                              running               , runningRunning
    logical                                          , intent(in   ) :: runningSmallScalesOnly
    !![
    <constructorAssign variables="index_, wavenumberReference, running, runningRunning, runningSmallScalesOnly"/>
    !!]

    if     (                                                                                                                                                                                               &
         &   (                                                                                                                                                                                             &
         &     running        > 0.0d0                                                                                                                                                                      &
         &    .or.                                                                                                                                                                                         &
         &     runningRunning > 0.0d0                                                                                                                                                                      &
         &   )                                                                                                                                                                                             &
         &  .and.                                                                                                                                                                                          &
         &   .not.runningSmallScalesOnly                                                                                                                                                                   &
         & )                                                                                                                                                                                               &
         & call Warn(                                                                                                                                                                                      &
         &           'primordial power spectra with positive running can lead to divergent σ(M) integrals - if this happens consider setting:'                                                //char(10)// &
         &           '    <'//displayBlue()//'powerSpectrumPrimordial'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"powerLaw"'//displayReset()//'>' //char(10)// &
         &           '      <'//displayBlue()//'runningSmallScalesOnly'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"true"'//displayReset()//'/>'   //char(10)// &
         &           '    </'//displayBlue()//'powerSpectrumPrimordial'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"powerLaw"'//displayReset()//'>'             &
         &          )
    return
  end function powerLawConstructorInternal

  double precision function powerLawPower(self,wavenumber)
    !!{
    Return the primordial power spectrum at the given {\normalfont \ttfamily wavenumber}.
    !!}
    implicit none
    class           (powerSpectrumPrimordialPowerLaw), intent(inout) :: self
    double precision                                 , intent(in   ) :: wavenumber
    double precision                                                 :: indexLocal

    if (self%runningSmallScalesOnly .and. wavenumber < self%wavenumberReference) then
       indexLocal=+self%index_
    else
       indexLocal=+self%index_                   &
            &     +1.0d0/2.0d0                   &
            &     *self%running                  &
            &     *log(                          &
            &          +     wavenumber          &
            &          /self%wavenumberReference &
            &         )                          &
            &     +1.0d0/6.0d0                   &
            &     *self%runningRunning           &
            &     *log(                          &
            &          +     wavenumber          &
            &          /self%wavenumberReference &
            &         )**2
    end if
    powerLawPower=+(                             &
         &          +     wavenumber             &
         &          /self%wavenumberReference    &
         &         )**indexLocal
    return
  end function powerLawPower

  double precision function powerLawLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the primordial power spectrum at the given {\normalfont \ttfamily wavenumber}.
    !!}
    implicit none
    class           (powerSpectrumPrimordialPowerLaw), intent(inout) :: self
    double precision                                 , intent(in   ) :: wavenumber

    if (self%runningSmallScalesOnly .and. wavenumber < self%wavenumberReference) then
       powerLawLogarithmicDerivative=+self%index_
    else
       powerLawLogarithmicDerivative=+self%index_                   &
            &                        +self%running                  &
            &                        *log(                          &
            &                             +     wavenumber          &
            &                             /self%wavenumberReference &
            &                            )                          &
            &                        +1.0d0/2.0d0                   &
            &                        *self%runningRunning           &
            &                        *log(                          &
            &                             +     wavenumber          &
            &                             /self%wavenumberReference &
            &                            )**2
    end if
    return
  end function powerLawLogarithmicDerivative
