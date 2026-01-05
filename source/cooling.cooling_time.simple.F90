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
  Implementation of a simple cooling time class.
  !!}

  use :: Chemical_States  , only : chemicalState  , chemicalStateClass
  use :: Cooling_Functions, only : coolingFunction, coolingFunctionClass

  !![
  <coolingTime name="coolingTimeSimple">
   <description>
  !!]

  !![
    A cooling time class in which the cooling time is simply
    \begin{equation}
     t_\mathrm{cool} = {N \over 2} {\mathrm{k}_\mathrm{B} T n_\mathrm{tot} \over \Lambda},
    \end{equation}
    where $N=${\normalfont \ttfamily [degreesOfFreedom]} is the number of degrees of freedom in the cooling gas which has
    temperature $T$ and total particle number density (including electrons) $n_\mathrm{tot}$ and $\Lambda$ is the cooling
    function.
   </description>
  </coolingTime>
  !!]
  type, extends(coolingTimeClass) :: coolingTimeSimple
     !!{
     Implementation of cooling time calculation (based on the ratio of the thermal energy density to the volume cooling rate).
     !!}
     private
     class           (coolingFunctionClass), pointer :: coolingFunction_ => null()
     class           (chemicalStateClass  ), pointer :: chemicalState_   => null()
     double precision                                :: degreesOfFreedom
   contains
     final     ::                                   simpleDestructor
     procedure :: time                           => simpleTime
     procedure :: gradientDensityLogarithmic     => simpleGradientDensityLogarithmic
     procedure :: gradientTemperatureLogarithmic => simpleGradientTemperatureLogarithmic
  end type coolingTimeSimple

  interface coolingTimeSimple
     !!{
     Constructors for the simple cooling time class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface coolingTimeSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the simple cooling time class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coolingTimeSimple   )                :: self
    type            (inputParameters     ), intent(inout) :: parameters
    class           (coolingFunctionClass), pointer       :: coolingFunction_
    class           (chemicalStateClass  ), pointer       :: chemicalState_
    double precision                                      :: degreesOfFreedom

    !![
    <inputParameter>
      <name>degreesOfFreedom</name>
      <source>parameters</source>
      <defaultValue>3.0d0</defaultValue>
      <description>Number of degrees of freedom to assume when computing the energy density of cooling gas in the ``simple'' cooling time class.</description>
    </inputParameter>
    <objectBuilder class="coolingFunction" name="coolingFunction_" source="parameters"/>
    <objectBuilder class="chemicalState"   name="chemicalState_"   source="parameters"/>
    !!]
    self=coolingTimeSimple(degreesOfFreedom,coolingFunction_,chemicalState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingFunction_"/>
    <objectDestructor name="chemicalState_"  />
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(degreesOfFreedom,coolingFunction_,chemicalState_) result(self)
    !!{
    Internal constructor for the simple cooling time class.
    !!}
    implicit none
    type            (coolingTimeSimple   )                        :: self
    double precision                      , intent(in   )         :: degreesOfFreedom
    class           (coolingFunctionClass), intent(in   ), target :: coolingFunction_
    class           (chemicalStateClass  ), intent(in   ), target :: chemicalState_
    !![
    <constructorAssign variables="degreesOfFreedom, *coolingFunction_, *chemicalState_"/>
    !!]

    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{
    Destructor for the simple cooling time class.
    !!}
    implicit none
    type(coolingTimeSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingFunction_"/>
    <objectDestructor name="self%chemicalState_"  />
    !!]
    return
  end subroutine simpleDestructor

  double precision function simpleTime(self,node,temperature,density,gasAbundances,chemicalDensities,radiation)
    !!{
    Compute the cooling time (in Gyr) for gas at the given {\normalfont \ttfamily temperature} (in Kelvin), {\normalfont \ttfamily density} (in $M_\odot$
    Mpc$^{-3}$), composition specified by {\normalfont \ttfamily gasAbundances} and experiencing a radiation field as described by {\normalfont \ttfamily radiation}.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear          , massSolar, megaParsec
    use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : hecto
    use :: Numerical_Constants_Units       , only : ergs
    implicit none
    class           (coolingTimeSimple   ), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node
    double precision                      , intent(in   ) :: density                       , temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances  ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass ), intent(inout) :: radiation
    ! Effectively infinite time (for arbitrarily long cooling times).
    double precision                      , parameter     :: timeLarge              =1.0d10
    double precision                                      :: coolingFunctionValue          , energyDensityThermal , &
         &                                                   numberDensityAllSpecies       , numberDensityHydrogen
    !$GLC attributes unused :: node

    ! Compute number density of hydrogen (in cm⁻³).
    numberDensityHydrogen  =  +density                                    &
         &                    *gasAbundances   %hydrogenMassFraction()    &
         &                    *massSolar                                  &
         &                    /massHydrogenAtom                           &
         &                    /hecto                                  **3 &
         &                    /megaParsec                             **3
    ! Get the number density of all species, including electrons.
    numberDensityAllSpecies=+                                                numberDensityHydrogen                                                        &
         &                  /     gasAbundances %hydrogenNumberFraction(                                                                                ) &
         &                  +self%chemicalState_%electronDensity       (     numberDensityHydrogen,temperature,gasAbundances                  ,radiation)
    ! Get the cooling function (in ergs cm⁻³ s⁻¹).
    coolingFunctionValue   =+self%coolingFunction_%coolingFunction     (node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    ! Compute the cooling time.
    if (coolingFunctionValue > 0.0d0) then
       ! Determine the thermal energy density of the gas (in ergs cm⁻³).
       energyDensityThermal=+self%degreesOfFreedom   &
            &               /2.0d0                   &
            &               *boltzmannsConstant      &
            &               *temperature             &
            &               *numberDensityAllSpecies &
            &               /ergs
       simpleTime          =+energyDensityThermal    &
            &               /coolingFunctionValue    &
            &               /gigaYear
    else
       simpleTime          =+timeLarge
    end if
    return
  end function simpleTime

  double precision function simpleGradientDensityLogarithmic(self,node,temperature,density,gasAbundances,chemicalDensities,radiation)
    !!{
    Return $\d\ln t_\mathrm{cool}/\d\ln \rho$ for gas at the given {\normalfont \ttfamily temperature} (in Kelvin), {\normalfont \ttfamily density} (in $M_\odot$
    Mpc$^{-3}$), composition specified by {\normalfont \ttfamily gasAbundances} and experiencing a radiation field as described by {\normalfont \ttfamily radiation}.
    !!}
   implicit none
    class           (coolingTimeSimple   ), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node
    double precision                      , intent(in   ) :: density          , temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances  ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass ), intent(inout) :: radiation
    !$GLC attributes unused :: node

    simpleGradientDensityLogarithmic=+1.0d0                                                                                                                    &
         &                           -self%coolingFunction_%coolingFunctionDensityLogSlope(node,density,temperature,gasAbundances,chemicalDensities,radiation)
    return
  end function simpleGradientDensityLogarithmic

  double precision function simpleGradientTemperatureLogarithmic(self,node,temperature,density,gasAbundances,chemicalDensities,radiation)
    !!{
    Return $\d\ln t_\mathrm{cool}/\d\ln T$ for gas at the given {\normalfont \ttfamily temperature} (in Kelvin), {\normalfont \ttfamily density} (in $M_\odot$
    Mpc$^{-3}$), composition specified by {\normalfont \ttfamily gasAbundances} and experiencing a radiation field as described by {\normalfont \ttfamily radiation}.
    !!}
    implicit none
    class           (coolingTimeSimple   ), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node
    double precision                      , intent(in   ) :: density          , temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances  ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass ), intent(inout) :: radiation
    !$GLC attributes unused :: node

    simpleGradientTemperatureLogarithmic=-self%coolingFunction_%coolingFunctionTemperatureLogSlope(node,density,temperature,gasAbundances,chemicalDensities,radiation)
    return
  end function simpleGradientTemperatureLogarithmic
