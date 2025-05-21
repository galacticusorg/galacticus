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

  type, public :: stellarFeedbackOutflowsList
     class(stellarFeedbackOutflowsClass), pointer :: stellarFeedbackOutflows => null()
     type (stellarFeedbackOutflowsList ), pointer :: next                    => null()
  end type stellarFeedbackOutflowsList

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsSummation">
   <description>A photon source class for summation sources.</description>
   <linkedList type="stellarFeedbackOutflowsList" variable="stellarFeedbackOutflowss" next="next" object="stellarFeedbackOutflows" objectType="stellarFeedbackOutflowsClass"/>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsSummation
     !!{
     Implementation of a summation stellar feedback class.
     !!}
     private
     type   (stellarFeedbackOutflowsList), pointer :: stellarFeedbackOutflowss => null()
   contains
     final     ::                summationDestructor
     procedure :: outflowRate => summationOutflowRate
  end type stellarFeedbackOutflowsSummation

  interface stellarFeedbackOutflowsSummation
     !!{
     Constructors for the \refClass{stellarFeedbackOutflowsSummation} stellar feedback class.
     !!}
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface stellarFeedbackOutflowsSummation
  
contains

  function summationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarFeedbackOutflowsSummation} stellar feedback class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (stellarFeedbackOutflowsSummation), target        :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    type   (stellarFeedbackOutflowsList     ), pointer       :: stellarFeedbackOutflows_
    integer                                                  :: i

    stellarFeedbackOutflows_ => null()
    do i=1,parameters%copiesCount('stellarFeedbackOutflows',zeroIfNotPresent=.true.)
       if (associated(stellarFeedbackOutflows_)) then
          allocate(stellarFeedbackOutflows_%next)
          stellarFeedbackOutflows_ => stellarFeedbackOutflows_%next
       else
          allocate(self%stellarFeedbackOutflowss)
          stellarFeedbackOutflows_ => self%stellarFeedbackOutflowss
       end if
       !![
       <objectBuilder class="stellarFeedbackOutflows" name="stellarFeedbackOutflows_%stellarFeedbackOutflows" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="stellarFeedbackOutflows"/>
    !!]
    return
  end function summationConstructorParameters

  function summationConstructorInternal(stellarFeedbackOutflowss) result(self)
    !!{
    Internal constructor for the \refClass{stellarFeedbackOutflowsSummation} stellar feedback class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type(stellarFeedbackOutflowsSummation)                         :: self
    type(stellarFeedbackOutflowsList     ), target , intent(in   ) :: stellarFeedbackOutflowss
    type(stellarFeedbackOutflowsList     ), pointer                :: stellarFeedbackOutflows_

    self            %stellarFeedbackOutflowss => stellarFeedbackOutflowss
    stellarFeedbackOutflows_                  => stellarFeedbackOutflowss
    do while (associated(stellarFeedbackOutflows_))
       !![
       <referenceCountIncrement owner="stellarFeedbackOutflows_" object="stellarFeedbackOutflows"/>
       !!]
       stellarFeedbackOutflows_ => stellarFeedbackOutflows_%next
    end do
    return
  end function summationConstructorInternal

  subroutine summationDestructor(self)
    !!{
    Destructor for the \refClass{stellarFeedbackOutflowsSummation} stellar feedback class.
    !!}
    implicit none
    type(stellarFeedbackOutflowsSummation), intent(inout) :: self
    type(stellarFeedbackOutflowsList     ), pointer       :: stellarFeedbackOutflows, stellarFeedbackOutflowsNext

    stellarFeedbackOutflows => self%stellarFeedbackOutflowss
    do while (associated(stellarFeedbackOutflows))
       stellarFeedbackOutflowsNext => stellarFeedbackOutflows%next
       !![
       <objectDestructor name="stellarFeedbackOutflows%stellarFeedbackOutflows"/>
       !!]
       deallocate(stellarFeedbackOutflows)
       stellarFeedbackOutflows => stellarFeedbackOutflowsNext
    end do
    return
  end subroutine summationDestructor
  
  subroutine summationOutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Initialize the photon packet.
    !!}
    implicit none
    class           (stellarFeedbackOutflowsSummation), intent(inout) :: self
    type            (stellarFeedbackOutflowsList     ), pointer       :: stellarFeedbackOutflows
    class           (nodeComponent                   ), intent(inout) :: component
    double precision                                  , intent(in   ) :: rateEnergyInput        , rateStarFormation
    double precision                                  , intent(  out) :: rateOutflowEjective    , rateOutflowExpulsive
    double precision                                                  :: rateOutflowEjective_   , rateOutflowExpulsive_

    rateOutflowEjective  =  +0.0d0
    rateOutflowExpulsive =  +0.0d0
    stellarFeedbackOutflows      =>  self%stellarFeedbackOutflowss
    do while (associated(stellarFeedbackOutflows))
       call stellarFeedbackOutflows%stellarFeedbackOutflows%outflowRate(component,rateStarFormation,rateEnergyInput,rateOutflowEjective_,rateOutflowExpulsive_)
       rateOutflowEjective  =  +rateOutflowEjective +rateOutflowEjective_
       rateOutflowExpulsive =  +rateOutflowExpulsive+rateOutflowExpulsive_
       stellarFeedbackOutflows      =>  stellarFeedbackOutflows%next
    end do 
    return
  end subroutine summationOutflowRate
