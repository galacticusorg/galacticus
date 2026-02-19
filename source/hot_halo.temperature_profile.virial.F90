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
An implementation of the hot halo temperature class which uses an isothermal virial temperature.
!!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <hotHaloTemperatureProfile name="hotHaloTemperatureProfileVirial">
   <description>
    A hot halo temperature profile class which assumes an isothermal halo with a temperature equal to the virial temperature of
    the halo.
   </description>
  </hotHaloTemperatureProfile>
  !!]
  type, extends(hotHaloTemperatureProfileClass) :: hotHaloTemperatureProfileVirial
     !!{
     An implementation of the hot halo temperature profile class which uses an isothermal virial temperature.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::        virialDestructor
     procedure :: get => virialGet
  end type hotHaloTemperatureProfileVirial

  interface hotHaloTemperatureProfileVirial
     !!{
     Constructors for the virial hot halo temperature profile class.
     !!}
     module procedure virialConstructorParameters
     module procedure virialConstructorInternal
  end interface hotHaloTemperatureProfileVirial

contains

  function virialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the virial cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloTemperatureProfileVirial)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=hotHaloTemperatureProfileVirial(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function virialConstructorParameters

  function virialConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the virial cooling rate class.
    !!}
    implicit none
    type (hotHaloTemperatureProfileVirial)                        :: self
    class(darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function virialConstructorInternal

  subroutine virialDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloTemperatureProfileVirial} hot halo temperature profile class.
    !!}
    implicit none
    type(hotHaloTemperatureProfileVirial), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine virialDestructor

  function virialGet(self,node) result(kinematicsDistribution_)
    !!{
    Return the virial hot halo temperature distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Mass_Distributions              , only : kinematicsDistributionIsothermal
    use :: Numerical_Constants_Astronomical, only : meanAtomicMassPrimordial
    implicit none
    class(kinematicsDistributionClass    ), pointer       :: kinematicsDistribution_
    class(hotHaloTemperatureProfileVirial), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    ! Create an isothermal kinematics distribution.
    allocate(kinematicsDistributionIsothermal :: kinematicsDistribution_)
    select type(kinematicsDistribution_)
    type is (kinematicsDistributionIsothermal)
       !![
       <referenceConstruct object="kinematicsDistribution_">
	 <constructor>
           kinematicsDistributionIsothermal(                                                                         &amp;
             &amp;                          temperature_  =self%darkMatterHaloScale_%temperatureVirial       (node), &amp;
             &amp;                          massAtomicMean=                          meanAtomicMassPrimordial        &amp;
             &amp;                         )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function virialGet
