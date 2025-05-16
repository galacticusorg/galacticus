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
  An implementation of atomic recombination cooling rates which are a fixed multiple of the recombination rate.
  !!}

  use :: Atomic_Rates_Recombination_Radiative, only : atomicRecombinationRateRadiativeClass
  
  !![
  <atomicRecombinationRateRadiativeCooling name="atomicRecombinationRateRadiativeCoolingFixed">
   <description>Atomic radiative cooling rateswhich are a fixed multiple of the recombination rate,  $\beta = \gamma \alpha$ where $\alpha$ is the corresponding radiative recombination coefficient and $\gamma$ is a parameter.</description>
  </atomicRecombinationRateRadiativeCooling>
  !!]
  type, extends(atomicRecombinationRateRadiativeCoolingClass) :: atomicRecombinationRateRadiativeCoolingFixed
     !!{
     A recombination cooling rate class assuming a thermal electron distribution.
     !!}
     private
     class           (atomicRecombinationRateRadiativeClass), pointer :: atomicRecombinationRateRadiative_ => null()
     double precision                                                 :: gamma
   contains
     final     ::         fixedDestructor
     procedure :: rate => fixedRate
  end type atomicRecombinationRateRadiativeCoolingFixed

  interface atomicRecombinationRateRadiativeCoolingFixed
     !!{
     Constructors for the \refClass{atomicRecombinationRateRadiativeCoolingFixed} atomic radiative recombination class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface atomicRecombinationRateRadiativeCoolingFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{atomicRecombinationRateRadiativeCoolingFixed} atomic radiative recombination class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (atomicRecombinationRateRadiativeCoolingFixed)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (atomicRecombinationRateRadiativeClass       ), pointer       :: atomicRecombinationRateRadiative_
    double precision                                                              :: gamma

    !![
    <inputParameter>
      <name>gamma</name>
      <description>The multiplicative factor, $\gamma$, used to compute the cooling coefficient.</description>
      <source>parameters</source>
      <defaultValue>0.67d0</defaultValue>
    </inputParameter>
    <objectBuilder class="atomicRecombinationRateRadiative" name="atomicRecombinationRateRadiative_" source="parameters"/>
    !!]
    self=atomicRecombinationRateRadiativeCoolingFixed(gamma,atomicRecombinationRateRadiative_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="atomicRecombinationRateRadiative_"/>
    !!]
    return
  end function fixedConstructorParameters
  
  function fixedConstructorInternal(gamma,atomicRecombinationRateRadiative_) result(self)
    !!{
    Internal constructor for the \refClass{atomicRecombinationRateRadiativeCoolingFixed} atomic radiative recombination class.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Table_Labels    , only : extrapolationTypeExtrapolate
    implicit none
    type            (atomicRecombinationRateRadiativeCoolingFixed)                        :: self
    class           (atomicRecombinationRateRadiativeClass       ), intent(in   ), target :: atomicRecombinationRateRadiative_
    double precision                                              , intent(in   )         :: gamma
    !![
    <constructorAssign variables="gamma, *atomicRecombinationRateRadiative_"/>
    !!]

    return
  end function fixedConstructorInternal

  subroutine fixedDestructor(self)
    !!{
    Destructor for the \refClass{atomicRecombinationRateRadiativeCoolingFixed} recombination cooling class.
    !!}
    implicit none
    type(atomicRecombinationRateRadiativeCoolingFixed), intent(inout) :: self

    !![
    <objectDestructor name="self%atomicRecombinationRateRadiative_"/>
    !!]
    return
  end subroutine fixedDestructor

  double precision function fixedRate(self,atomicNumber,ionizationState,temperature,level)
    !!{
    Returns the cooling rate coefficient.
    !!}
    implicit none
    class           (atomicRecombinationRateRadiativeCoolingFixed), intent(inout)           :: self
    integer                                                       , intent(in   )           :: atomicNumber, ionizationState
    double precision                                              , intent(in   )           :: temperature
    type            (enumerationRecombinationCaseType            ), intent(in   ), optional :: level

    fixedRate=+self%gamma                                                                                   &
         &     *self%atomicRecombinationRateRadiative_%rate(atomicNumber,ionizationState,temperature,level)
    return
  end function fixedRate
