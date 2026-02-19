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
  Implements a transfer function class based on the non-cold dark matter fitting function of \cite{murgia_non-cold_2017},
  restricted to $\gamma=-5$ and with properties specified via the half-mode mass, $M_{1/2}$ and the slope of the transfer function
  at the corresponding half-mode wavenumber, $\mathrm{d}\log T / \mathrm{d} \log k$.
  !!}

  !![
  <transferFunction name="transferFunctionHalfModeSlope">
    <description>
      A transfer function class based on the non-cold dark matter fitting function of \cite{murgia_non-cold_2017}, restricted to
      $\gamma=-5$ and with properties specified via the half-mode mass, $M_{1/2}$ and the slope of the transfer function at the
      corresponding half-mode wavenumber, $\mathrm{d}\log T / \mathrm{d} \log k$. Specifically,
      \begin{equation}
        k_{1/2} = \left(\frac{3 M_{1/2}}{4 \pi \bar{\rho}}\right)^{1/3},
      \end{equation}
      where $\bar{\rho}$ is the mean density of the universe. In terms of the \cite{murgia_non-cold_2017} parameters we have:
      \begin{equation}
        \beta = \frac{1}{\gamma} \frac{2^{-1/\gamma}}{2^{-1/\gamma}-1} \frac{\mathrm{d}\log T}{\mathrm{d} \log k},
      \end{equation}
      and
      \begin{equation}
        \alpha = \frac{(2^{-1/\gamma}-1)^{1/\beta}}{k_\mathrm{1/2}}.
      \end{equation}
    </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionMurgia2017) :: transferFunctionHalfModeSlope
     !!{
     A transfer function class based on the non-cold dark matter fitting function of \cite{murgia_non-cold_2017}.
     !!}
     private
     double precision :: massHalfMode, slopeHalfMode
   contains
  end type transferFunctionHalfModeSlope

  interface transferFunctionHalfModeSlope
     !!{
     Constructors for the \refClass{transferFunctionHalfModeSlope} transfer function class.
     !!}
     module procedure halfModeSlopeConstructorParameters
     module procedure halfModeSlopeConstructorInternal
  end interface transferFunctionHalfModeSlope

contains
  
  function halfModeSlopeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{transferFunctionHalfModeSlope} transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionHalfModeSlope)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (transferFunctionClass        ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    double precision                                               :: slopeHalfMode       , massHalfMode, &
         &                                                            redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunction')) call Error_Report("an explicit 'transferFunction' must be given"//{introspection:location})
    ! Read parameters.
    !![
    <inputParameter>
      <name>massHalfMode</name>
      <source>parameters</source>
      <description>The half-mode mass of the transfer function.</description>
    </inputParameter>
    <inputParameter>
      <name>slopeHalfMode</name>
      <source>parameters</source>
      <description>The logarithmic slope of the transfer function at the half-mode mass.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%equalityEpochMatterRadiation(requestTypeExpansionFactor))</defaultValue>
      <description>The redshift of the epoch at which the transfer function is defined.</description>
    </inputParameter>
    !!]
    self=transferFunctionhalfModeSlope(transferFunctionCDM,massHalfMode,slopeHalfMode,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="transferFunctionCDM" />
    !!]
    return
  end function halfModeSlopeConstructorParameters

  function halfModeSlopeConstructorInternal(transferFunctionCDM,massHalfMode,slopeHalfMode,time,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{transferFunctionHalfModeSlope} transfer function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (transferFunctionHalfModeSlope)                        :: self
    class           (transferFunctionClass        ), target, intent(in   ) :: transferFunctionCDM
    double precision                                       , intent(in   ) :: massHalfMode        , slopeHalfMode, &
         &                                                                    time
    class           (cosmologyParametersClass     ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), target, intent(in   ) :: cosmologyFunctions_
    double precision                                                       :: wavenumberHalfMode
    !![
    <constructorAssign variables="*transferFunctionCDM, massHalfMode, slopeHalfMode, time, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    ! Compute the corresponding redshift.
    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    ! Compute the parameters for the underlying Murgia et al. (2017) transfer function.
    wavenumberHalfMode=+Pi                                                    &
         &             /(                                                     &
         &               +3.0d0                                               &
         &               *massHalfMode                                        &
         &               /4.0d0                                               &
         &               /Pi                                                  &
         &               /self%cosmologyParameters_%OmegaMatter    ()         &
         &               /self%cosmologyParameters_%densityCritical()         &
         &              )**(1.0d0/3.0d0)
    self%gamma        =-5.0d0
    self%beta         =+slopeHalfMode                                         &
         &             * 2.0d0**(-1.0d0/self%gamma)                           &
         &             /(2.0d0**(-1.0d0/self%gamma)-1.0d0)                    &
         &             /self%gamma
    self%alpha        =+(2.0d0**(-1.0d0/self%gamma)-1.0d0)**(1.0d0/self%beta) &
         &             /wavenumberHalfMode
    return
  end function halfModeSlopeConstructorInternal
