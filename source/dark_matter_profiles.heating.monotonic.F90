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

  !+    Contributions to this file made by: Xiaolong Du.

  !!{
  A dark matter halo profile heating class which takes another heating source and enforces monotonic heating energy perturbation.
  !!}

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingMonotonic">
    <description>
      A dark matter profile heating model builds \refClass{massDistributionHeatingMonotonic} objects to enforce monotonic heating
      energy perturbations.
    </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingMonotonic
     !!{
     A dark matter profile heating class which takes another heating source and enforces monotonic heating energy perturbation.
     !!}
     private
     class(darkMatterProfileHeatingClass), pointer :: darkMatterProfileHeating_ => null()
   contains
     final     ::        monotonicDestructor
     procedure :: get => monotonicGet
  end type darkMatterProfileHeatingMonotonic

  interface darkMatterProfileHeatingMonotonic
     !!{
     Constructors for the \refClass{darkMatterProfileHeatingMonotonic} dark matter profile heating class.
     !!}
     module procedure monotonicConstructorParameters
     module procedure monotonicConstructorInternal
  end interface darkMatterProfileHeatingMonotonic

contains

  function monotonicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileHeatingMonotonic} dark matter profile heating class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileHeatingMonotonic), target        :: self
    type (inputParameters                  ), intent(inout) :: parameters
    class(darkMatterProfileHeatingClass    ), pointer       :: darkMatterProfileHeating_

    !![
    <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    !!]
    self=darkMatterProfileHeatingMonotonic(darkMatterProfileHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileHeating_"/>
    !!]
    return
  end function monotonicConstructorParameters

  function monotonicConstructorInternal(darkMatterProfileHeating_) result(self)
    !!{
    Internal constructor for the ``monotonic'' dark matter profile heating class.
    !!}
    implicit none
    type (darkMatterProfileHeatingMonotonic)                        :: self
    class(darkMatterProfileHeatingClass    ), target, intent(in   ) :: darkMatterProfileHeating_
    !![
    <constructorAssign variables="*darkMatterProfileHeating_"/>
    !!]

    return
  end function monotonicConstructorInternal

  subroutine monotonicDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileHeatingMonotonic} dark matter profile heating class.
    !!}
    implicit none
    type(darkMatterProfileHeatingMonotonic), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileHeating_"/>
    !!]
    return
  end subroutine monotonicDestructor

  function monotonicGet(self,node) result(massDistributionHeating_)
    !!{
    Return the dark matter mass distribution heating for the given {\normalfont \ttfamily node}.
    !!}
    use :: Mass_Distributions, only : massDistributionHeatingMonotonic
    implicit none
    class(massDistributionHeatingClass     ), pointer       :: massDistributionHeating_
    class(darkMatterProfileHeatingMonotonic), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    class(massDistributionHeatingClass     ), pointer       :: massDistributionHeatingDecorated
 
    ! Create the mass distribution.
    allocate(massDistributionHeatingMonotonic :: massDistributionHeating_)
    select type(massDistributionHeating_)
    type is (massDistributionHeatingMonotonic)
       massDistributionHeatingDecorated => self%darkMatterProfileHeating_%get(node)
       !![
       <referenceConstruct object="massDistributionHeating_">
	 <constructor>
           massDistributionHeatingMonotonic(                                                          &amp;
           &amp;                            massDistributionHeating_=massDistributionHeatingDecorated &amp;
           &amp;                           )
	 </constructor>
       </referenceConstruct>
       <objectDestructor name="massDistributionHeatingDecorated"/>
       !!]
    end select
    return
  end function monotonicGet
