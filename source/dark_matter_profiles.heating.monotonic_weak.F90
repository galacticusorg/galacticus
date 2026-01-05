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
  A dark matter halo profile heating class which takes another heating source and enforces monotonic heating energy perturbation.
  !!}

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingMonotonicWeak">
    <description>
      A dark matter profile heating model builds \refClass{massDistributionHeatingMonotonicWeak} objects to enforce monotonic heating
      energy perturbations. This classes enforces a weaker condition (compared to \refClass{darkMatterProfileHeatingMonotonic}).
    </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingMonotonicWeak
     !!{
     A dark matter profile heating class which takes another heating source and enforces monotonic heating energy perturbation.
     !!}
     private
     class           (darkMatterProfileHeatingClass), pointer :: darkMatterProfileHeating_ => null()
     double precision                                         :: toleranceShellCrossing
  contains
     final     ::        monotonicWeakDestructor
     procedure :: get => monotonicWeakGet
  end type darkMatterProfileHeatingMonotonicWeak

  interface darkMatterProfileHeatingMonotonicWeak
     !!{
     Constructors for the \refClass{darkMatterProfileHeatingMonotonicWeak} dark matter profile heating class.
     !!}
     module procedure monotonicWeakConstructorParameters
     module procedure monotonicWeakConstructorInternal
  end interface darkMatterProfileHeatingMonotonicWeak

contains

  function monotonicWeakConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileHeatingMonotonicWeak} dark matter profile heating class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileHeatingMonotonicWeak), target        :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (darkMatterProfileHeatingClass        ), pointer       :: darkMatterProfileHeating_
    double precision                                                       :: toleranceShellCrossing
    !![
    <inputParameter>
      <name>toleranceShellCrossing</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The tolerance adopted in determining if the no-shell-crossing assumption is valid.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    !!]
    self=darkMatterProfileHeatingMonotonicWeak(toleranceShellCrossing,darkMatterProfileHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileHeating_"/>
    !!]
    return
  end function monotonicWeakConstructorParameters

  function monotonicWeakConstructorInternal(toleranceShellCrossing, darkMatterProfileHeating_) result(self)
    !!{
    Internal constructor for the ``monotonicWeak'' dark matter profile heating class.
    !!}
    implicit none
    type            (darkMatterProfileHeatingMonotonicWeak)                        :: self
    class           (darkMatterProfileHeatingClass        ), target, intent(in   ) :: darkMatterProfileHeating_
    double precision                                               , intent(in   ) :: toleranceShellCrossing
    !![
    <constructorAssign variables="toleranceShellCrossing, *darkMatterProfileHeating_"/>
    !!]

    return
  end function monotonicWeakConstructorInternal

  subroutine monotonicWeakDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileHeatingMonotonicWeak} dark matter profile heating class.
    !!}
    implicit none
    type(darkMatterProfileHeatingMonotonicWeak), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileHeating_"/>
    !!]
    return
  end subroutine monotonicWeakDestructor

  function monotonicWeakGet(self,node) result(massDistributionHeating_)
    !!{
    Return the dark matter mass distribution heating for the given {\normalfont \ttfamily node}.
    !!}
    use :: Mass_Distributions, only : massDistributionHeatingMonotonicWeak
    implicit none
    class(massDistributionHeatingClass         ), pointer       :: massDistributionHeating_
    class(darkMatterProfileHeatingMonotonicWeak), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    class(massDistributionHeatingClass         ), pointer       :: massDistributionHeatingDecorated
 
    ! Create the mass distribution.
    allocate(massDistributionHeatingMonotonicWeak :: massDistributionHeating_)
    select type(massDistributionHeating_)
    type is (massDistributionHeatingMonotonicWeak)
       massDistributionHeatingDecorated => self%darkMatterProfileHeating_%get(node)
       !![
       <referenceConstruct object="massDistributionHeating_">
	 <constructor>
           massDistributionHeatingMonotonicWeak(                                                            &amp;
           &amp;                            toleranceShellCrossing  =self%toleranceShellCrossing          , &amp;
           &amp;                            massDistributionHeating_=     massDistributionHeatingDecorated  &amp;
           &amp;                           )
	 </constructor>
       </referenceConstruct>
       <objectDestructor name="massDistributionHeatingDecorated"/>
       !!]
    end select
    return
  end function monotonicWeakGet
