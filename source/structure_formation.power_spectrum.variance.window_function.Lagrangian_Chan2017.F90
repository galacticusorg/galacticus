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
  Provides a power spectrum window function class that implements the Lagrangian filter of \cite{chan_effective_2017}.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionLagrangianChan2017">
   <description>A power spectrum window function class that implements the Lagrangian filter of \cite{chan_effective_2017}.</description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionLagrangianChan2017
     !!{
     The Lagrangian smooth power spectrum window function of \cite{chan_effective_2017} class.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: f                             , phi
    contains
     final     ::                      lagrangianChan2017Destructor
     procedure :: value             => lagrangianChan2017Value
     procedure :: wavenumberMaximum => lagrangianChan2017WavenumberMaximum
  end type powerSpectrumWindowFunctionLagrangianChan2017

  interface powerSpectrumWindowFunctionLagrangianChan2017
     !!{
     Constructors for the \refClass{powerSpectrumWindowFunctionLagrangianChan2017} power spectrum window function class.
     !!}
     module procedure lagrangianChan2017ConstructorParameters
     module procedure lagrangianChan2017ConstructorInternal
  end interface powerSpectrumWindowFunctionLagrangianChan2017

contains

  function lagrangianChan2017ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{powerSpectrumWindowFunctionLagrangianChan2017} power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionLagrangianChan2017)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (cosmologyParametersClass                     ), pointer       :: cosmologyParameters_
    double precision                                                               :: f

    !![
    <inputParameter>
      <name>f</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The parameter ``$f$'' which defines the scale of the Gaussian window.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self = powerSpectrumWindowFunctionLagrangianChan2017(cosmologyParameters_,f)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function lagrangianChan2017ConstructorParameters

  function lagrangianChan2017ConstructorInternal(cosmologyParameters_,f) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumWindowFunctionLagrangianChan2017} power spectrum window function class.
    !!}
    implicit none
    type            (powerSpectrumWindowFunctionLagrangianChan2017)                        :: self
    class           (cosmologyParametersClass                     ), intent(in   ), target :: cosmologyParameters_
    double precision                                               , intent(in   )         :: f
    !![
    <constructorAssign variables="f, *cosmologyParameters_"/>
    !!]

    self%phi=5.0d0/sqrt(f)
    return
  end function lagrangianChan2017ConstructorInternal

  subroutine lagrangianChan2017Destructor(self)
    !!{
    Destructor for the \refClass{powerSpectrumWindowFunctionLagrangianChan2017} power spectrum window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionLagrangianChan2017), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine lagrangianChan2017Destructor

  double precision function lagrangianChan2017Value(self,wavenumber,smoothingMass,time)
    !!{
    Smooth-$k$ space power spectrum window function proposed in \cite{leo_new_2018}.
    spectrum.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumWindowFunctionLagrangianChan2017), intent(inout) :: self
    double precision                                               , intent(in   ) :: smoothingMass          , wavenumber, &
         &                                                                            time
    double precision                                               , parameter     :: xSeriesMaximum  =1.0d-3
    double precision                                                               :: radiusLagrangian       , x          , &
         &                                                                            xSquared               , valueTopHat, &
         &                                                                            valueGaussian
    !$GLC attributes unused :: time

    radiusLagrangian=+(                                             &
         &             +3.0d0                                       &
         &             /4.0d0                                       &
         &             /Pi                                          &
         &             *smoothingMass                               &
         &             /self%cosmologyParameters_%OmegaMatter    () &
         &             /self%cosmologyParameters_%densityCritical() &
         &            )**(1.0d0/3.0d0)
    x               =+wavenumber                                    &
         &           *radiusLagrangian
    if (x <= 0.0d0) then
       lagrangianChan2017Value=+0.0d0
    else
       if (x <= xSeriesMaximum) then
          ! Use a series expansion of the window function for small x.
          xSquared   =+x**2
          valueTopHat=+1.0d0                        &
               &      +xSquared*(  -1.0d0/   10.0d0 &
               &      +xSquared* ( +1.0d0/  280.0d0 &
               &      +xSquared*  (-1.0d0/15120.0d0 &
               &                  )                 &
               &                 )                  &
               &                )
       else
          valueTopHat=+3.0d0       &
               &      *(           &
               &        +sin(x)    &
               &        -    x     &
               &        *cos(x)    &
               &       )           &
               &      /      x **3
       end if
       valueGaussian          =+exp(            &
            &                       -0.5d0      &
            &                       *(          &
            &                         +x        &
            &                         /self%phi &
            &                        )**2       &
            &                      )
       lagrangianChan2017Value=+valueTopHat     &
            &                  *valueGaussian
    end if
    return
  end function lagrangianChan2017Value

  double precision function lagrangianChan2017WavenumberMaximum(self,smoothingMass)
    !!{
    Maximum wavenumber for a top hat in real space window function Fourier transformed into $k$-space used in computing the
    variance of the power spectrum.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionLagrangianChan2017), intent(inout) :: self
    double precision                                               , intent(in   ) :: smoothingMass
    double precision                                               , parameter     :: wavenumberLarge=huge(1.0d0) ! Effective infinity.
    !$GLC attributes unused :: self, smoothingMass

    lagrangianChan2017WavenumberMaximum=wavenumberLarge
    return
  end function lagrangianChan2017WavenumberMaximum
