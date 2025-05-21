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
  Implements a radiative transfer outputter class which combines multiple other outputters.
  !!}

  type, public :: multiOutputterList
     class(radiativeTransferOutputterClass), pointer :: outputter_ => null()
     type (multiOutputterList             ), pointer :: next       => null()
  end type multiOutputterList

  !![
  <radiativeTransferOutputter name="radiativeTransferOutputterMulti">
   <description>A radiative transfer outputter class which combines multiple other outputters.</description>
   <linkedList type="multiOutputterList" variable="outputters" next="next" object="outputter_" objectType="radiativeTransferOutputterClass"/>
  </radiativeTransferOutputter>
  !!]
  type, extends(radiativeTransferOutputterClass) :: radiativeTransferOutputterMulti
     !!{
     Implementation of a radiative transfer outputter class which combines multiple other outputters.
     !!}
     private
     type(multiOutputterList), pointer :: outputters => null()
   contains
     final     ::                        multiDestructor
     procedure :: reset               => multiReset
     procedure :: sourceProperties    => multiSourceProperties
     procedure :: photonPacketEscapes => multiPhotonPacketEscapes
     procedure :: finalize            => multiFinalize
     procedure :: output              => multiOutput
  end type radiativeTransferOutputterMulti

  interface radiativeTransferOutputterMulti
     !!{
     Constructors for the \refClass{radiativeTransferOutputterMulti} radiative transfer outputter class.
     !!}
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface radiativeTransferOutputterMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferOutputterMulti} radiative transfer outputter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (radiativeTransferOutputterMulti)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    type   (multiOutputterList             ), pointer       :: outputter_
    integer                                                 :: i

    self      %outputters => null()
    outputter_            => null()
    do i=1,parameters%copiesCount('radiativeTransferOutputter',zeroIfNotPresent=.true.)
       if (associated(outputter_)) then
          allocate(outputter_%next)
          outputter_ => outputter_%next
       else
          allocate(self%outputters)
          outputter_ => self%outputters
       end if
       !![
       <objectBuilder class="radiativeTransferOutputter" name="outputter_%outputter_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="radiativeTransferOutputter"/>
    !!]
    return
  end function multiConstructorParameters

  function multiConstructorInternal(outputters) result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferOutputterMulti} analysis class.
    !!}
    implicit none
    type(radiativeTransferOutputterMulti)                        :: self
    type(multiOutputterList             ), target, intent(in   ) :: outputters
    type(multiOutputterList             ), pointer               :: outputter_

    self      %outputters => outputters
    outputter_            => outputters
    do while (associated(outputter_))
       !![
       <referenceCountIncrement owner="outputter_" object="outputter_"/>
       !!]
       outputter_ => outputter_%next
    end do
    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !!{
    Destructor for the \refClass{radiativeTransferOutputterMulti} analysis class.
    !!}
    implicit none
    type(radiativeTransferOutputterMulti), intent(inout) :: self
    type(multiOutputterList             ), pointer       :: outputter_, analysisNext

    if (associated(self%outputters)) then
       outputter_ => self%outputters
       do while (associated(outputter_))
          analysisNext => outputter_%next
          !![
          <objectDestructor name="outputter_%outputter_"/>
          !!]
          deallocate(outputter_)
          outputter_ => analysisNext
       end do
    end if
    return
  end subroutine multiDestructor

  subroutine multiReset(self)
    !!{
    Reset all outputters.
    !!}
    implicit none
    class(radiativeTransferOutputterMulti), intent(inout) :: self
    type (multiOutputterList             ), pointer       :: outputter_

    outputter_ => self%outputters
    do while (associated(outputter_))
       call outputter_%outputter_%reset()
       outputter_ => outputter_%next
    end do
    return
  end subroutine multiReset

  subroutine multiSourceProperties(self,radiativeTransferSource_,outputGroup)
    !!{
    Compute and output all source properties.
    !!}
    implicit none
    class(radiativeTransferOutputterMulti), intent(inout) :: self
    class(radiativeTransferSourceClass   ), intent(inout) :: radiativeTransferSource_
    type (hdf5Object                     ), intent(inout) :: outputGroup
    type (multiOutputterList             ), pointer       :: outputter_

    outputter_ => self%outputters
    do while (associated(outputter_))
       call outputter_%outputter_%sourceProperties(radiativeTransferSource_,outputGroup)
       outputter_ => outputter_%next
    end do
    return
  end subroutine multiSourceProperties

  subroutine multiPhotonPacketEscapes(self,photonPacket)
    !!{
    Process an escaping photon packet.
    !!}
    implicit none
    class(radiativeTransferOutputterMulti   ), intent(inout) :: self
    class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket
    type (multiOutputterList                ), pointer       :: outputter_

    outputter_ => self%outputters
    do while (associated(outputter_))
       call outputter_%outputter_%photonPacketEscapes(photonPacket)
       outputter_ => outputter_%next
    end do
    return
  end subroutine multiPhotonPacketEscapes

  subroutine multiFinalize(self)
    !!{
    Finalize the results.
    !!}
    implicit none
    class(radiativeTransferOutputterMulti), intent(inout) :: self
    type (multiOutputterList             ), pointer       :: outputter_

    outputter_ => self%outputters
    do while (associated(outputter_))
       call outputter_%outputter_%finalize()
       outputter_ => outputter_%next
    end do
    return
  end subroutine multiFinalize
  
  subroutine multiOutput(self,outputGroup)
    !!{
    Output the results.
    !!}
    implicit none
    class(radiativeTransferOutputterMulti), intent(inout) :: self
    type (hdf5Object                     ), intent(inout) :: outputGroup
    type (multiOutputterList             ), pointer       :: outputter_

    outputter_ => self%outputters
    do while (associated(outputter_))
       call outputter_%outputter_%output(outputGroup)
       outputter_ => outputter_%next
    end do
    return
  end subroutine multiOutput
  
