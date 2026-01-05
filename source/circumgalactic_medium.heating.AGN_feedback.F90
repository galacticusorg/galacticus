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
  Implements a \gls{cgm} heating class for AGN feedback.
  !!}

  use :: Black_Hole_CGM_Heating, only : blackHoleCGMHeatingClass

  !![
  <circumgalacticMediumHeating name="circumgalacticMediumHeatingAGNFeedback">
   <description>
    A \gls{cgm} heating class for AGN feedback.
   </description>
  </circumgalacticMediumHeating>
  !!]
  type, extends(circumgalacticMediumHeatingClass) :: circumgalacticMediumHeatingAGNFeedback
     !!{
     A \gls{cgm} heating class for AGN feedback.
     !!}
     private
     class(blackHoleCGMHeatingClass), pointer :: blackHoleCGMHeating_ => null()
   contains
     final     ::                agnFeedbackDestructor
     procedure :: heatingRate => agnFeedbackHeatingRate
  end type circumgalacticMediumHeatingAGNFeedback
  
  interface circumgalacticMediumHeatingAGNFeedback
     !!{
     Constructors for the \refClass{circumgalacticMediumHeatingAGNFeedback} class.
     !!}
     module procedure agnFeedbackConstructorParameters
     module procedure agnFeedbackConstructorInternal
  end interface circumgalacticMediumHeatingAGNFeedback

contains

  function agnFeedbackConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{circumgalacticMediumHeatingAGNFeedback} class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (circumgalacticMediumHeatingAGNFeedback)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(blackHoleCGMHeatingClass              ), pointer       :: blackHoleCGMHeating_

    
    !![
    <objectBuilder class="blackHoleCGMHeating" name="blackHoleCGMHeating_" source="parameters"/>
    !!]
    self=circumgalacticMediumHeatingAGNFeedback(blackHoleCGMHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleCGMHeating_"/>
    !!]
    return
  end function agnFeedbackConstructorParameters

  function agnFeedbackConstructorInternal(blackHoleCGMHeating_) result(self)
    !!{
    Internal constructor for the \refClass{circumgalacticMediumHeatingAGNFeedback} node operator class.
    !!}
    implicit none
    type (circumgalacticMediumHeatingAGNFeedback)                        :: self
    class(blackHoleCGMHeatingClass              ), intent(in   ), target :: blackHoleCGMHeating_
    !![
    <constructorAssign variables="*blackHoleCGMHeating_"/>
    !!]
    
    return
  end function agnFeedbackConstructorInternal

  subroutine agnFeedbackDestructor(self)
    !!{
    Destructor for the \refClass{circumgalacticMediumHeatingAGNFeedback} node operator class.
    !!}
    implicit none
    type(circumgalacticMediumHeatingAGNFeedback), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleCGMHeating_"/>
    !!]
    return
  end subroutine agnFeedbackDestructor

  double precision function agnFeedbackHeatingRate(self,node) result(rateHeating)
    !!{
    Compute the heating rate of the \gls{cgm} due to feedback from AGN.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class  (circumgalacticMediumHeatingAGNFeedback), intent(inout) :: self
    type   (treeNode                              ), intent(inout) :: node
    class  (nodeComponentBlackHole                ), pointer       :: blackHole
    integer                                                        :: countBlackHole, indexBlackHole

    rateHeating   =0.0d0
    countBlackHole=node%blackHoleCount()
    if (countBlackHole < 1) return
    ! Iterate over all black holes in the node.
    do indexBlackHole=1,countBlackHole
       ! Find the wind power from this black hole.
       blackHole   =>  node                     %blackHole  (instance=indexBlackHole)
       rateHeating =  +                          rateHeating                          &
            &         +self%blackHoleCGMHeating_%heatingRate(              blackHole)
    end do
    return
  end function agnFeedbackHeatingRate
  
