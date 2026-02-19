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
Implements a stellar population spectra postprocessor class which applies a sequence of other postprocessors.
!!}

  type, public :: postprocessorList
     class(stellarPopulationSpectraPostprocessorClass), pointer :: postprocessor_ => null()
     type (postprocessorList                         ), pointer :: next           => null()
  end type postprocessorList

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorSequence">
   <description>A sequence stellar population spectra postprocessor class.</description>
   <linkedList type="postprocessorList" variable="postprocessors" next="next" object="postprocessor_" objectType="stellarPopulationSpectraPostprocessorClass"/>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorSequence
     !!{
     A sequence stellar population spectra postprocessor class.
     !!}
     private
     type(postprocessorList), pointer :: postprocessors => null()
   contains
     final     ::                        sequenceDestructor
     procedure :: multiplier          => sequenceMultiplier
     procedure :: isRedshiftDependent => sequenceIsRedshiftDependent
  end type stellarPopulationSpectraPostprocessorSequence

  interface stellarPopulationSpectraPostprocessorSequence
     !!{
     Constructors for the \refClass{stellarPopulationSpectraPostprocessorSequence} stellar population spectra postprocessor class.
     !!}
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface stellarPopulationSpectraPostprocessorSequence

contains

  function sequenceConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{stellarPopulationSpectraPostprocessorSequence} stellar population spectra postprocessor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (stellarPopulationSpectraPostprocessorSequence)                :: self
    type   (inputParameters                              ), intent(inout) :: parameters
    type   (postprocessorList                            ), pointer       :: postprocessor_
    integer                                                               :: i

    self     %postprocessors => null()
    postprocessor_           => null()
    do i=1,parameters%copiesCount('stellarPopulationSpectraPostprocessor',zeroIfNotPresent=.true.)
       if (associated(postprocessor_)) then
          allocate(postprocessor_%next)
          postprocessor_ => postprocessor_%next
       else
          allocate(self%postprocessors)
          postprocessor_ => self%postprocessors
       end if
       !![
       <objectBuilder class="stellarPopulationSpectraPostprocessor" name="postprocessor_%postprocessor_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="stellarPopulationSpectraPostprocessor"/>
    !!]
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(postprocessors) result (self)
    !!{
    Internal constructor for the sequence stellar population spectra postprocessor class.
    !!}
    implicit none
    type(stellarPopulationSpectraPostprocessorSequence)                        :: self
    type(postprocessorList                            ), target, intent(in   ) :: postprocessors
    type(postprocessorList                            ), pointer               :: postprocessor_

    self          %postprocessors => postprocessors
    postprocessor_                => postprocessors
    do while (associated(postprocessor_))
       !![
       <referenceCountIncrement owner="postprocessor_" object="postprocessor_"/>
       !!]
       postprocessor_ => postprocessor_%next
    end do
    return
  end function sequenceConstructorInternal

  subroutine sequenceDestructor(self)
    !!{
    Destructor for the sequence stellar population spectra postprocessor class.
    !!}
    implicit none
    type(stellarPopulationSpectraPostprocessorSequence), intent(inout) :: self
    type(postprocessorList                            ), pointer       :: postprocessor_, postprocessorNext

    if (associated(self%postprocessors)) then
       postprocessor_ => self%postprocessors
       do while (associated(postprocessor_))
          postprocessorNext => postprocessor_%next
          !![
          <objectDestructor name="postprocessor_%postprocessor_"/>
          !!]
          deallocate(postprocessor_)
          postprocessor_ => postprocessorNext
       end do
    end if
    return
  end subroutine sequenceDestructor

  double precision function sequenceMultiplier(self,wavelength,age,redshift)
    !!{
    Implement an sequence stellar population spectra postprocessor.
    !!}
    implicit none
    class           (stellarPopulationSpectraPostprocessorSequence), intent(inout) :: self
    double precision                                               , intent(in   ) :: age           , redshift, &
         &                                                                            wavelength
    type            (postprocessorList                            ), pointer       :: postprocessor_

    sequenceMultiplier =  1.0d0
    postprocessor_     => self%postprocessors
    do while (associated(postprocessor_))
       sequenceMultiplier =  +sequenceMultiplier                                                &
            &                *postprocessor_%postprocessor_%multiplier(wavelength,age,redshift)
       postprocessor_     =>  postprocessor_%next
    end do
    return
  end function sequenceMultiplier

  logical function sequenceIsRedshiftDependent(self) result(isRedshiftDependent)
    !!{
    Return true if the postprocessor is redshift dependent.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorSequence), intent(inout) :: self
    type(postprocessorList                             ), pointer       :: postprocessor_

    isRedshiftDependent =  .false.
    postprocessor_      => self%postprocessors
    do while (associated(postprocessor_))
       if (postprocessor_%postprocessor_%isRedshiftDependent()) then
          isRedshiftDependent=.true.
          return
       end if
       postprocessor_ => postprocessor_%next
    end do
    return
  end function sequenceIsRedshiftDependent
