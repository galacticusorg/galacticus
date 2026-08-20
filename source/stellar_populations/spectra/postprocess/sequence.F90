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

!!{RST
Implements a stellar population spectra postprocessor class which applies a sequence of other postprocessors.
!!}

  type, public :: postprocessorList
     class(stellarPopulationSpectraPostprocessorClass), pointer :: postprocessor_ => null()
     type (postprocessorList                         ), pointer :: next           => null()
  end type postprocessorList

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorSequence" docformat="rst">
   <description>
   A stellar population spectra postprocessor class that applies a sequence of other postprocessors in order, computing the total multiplicative correction as the product of all individual postprocessor multipliers, allowing multiple physical effects (e.g., IGM absorption and dust) to be combined.
   </description>
   <linkedList type="postprocessorList" variable="postprocessors" next="next" object="postprocessor_" objectType="stellarPopulationSpectraPostprocessorClass"/>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorSequence
     !!{RST
     A sequence stellar population spectra postprocessor class.
     !!}
     private
     type(postprocessorList), pointer :: postprocessors => null()
   contains
     final     ::                        sequenceDestructor
     procedure :: multiplier          => sequenceMultiplier
     procedure :: ageRange            => sequenceAgeRange
     procedure :: ageWindowIsSharp    => sequenceAgeWindowIsSharp
     procedure :: isRedshiftDependent => sequenceIsRedshiftDependent
  end type stellarPopulationSpectraPostprocessorSequence

  interface stellarPopulationSpectraPostprocessorSequence
     !!{RST
     Constructors for the :galacticus-class:`stellarPopulationSpectraPostprocessorSequence` stellar population spectra postprocessor class.
     !!}
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface stellarPopulationSpectraPostprocessorSequence

contains

  function sequenceConstructorParameters(parameters) result (self)
    !!{RST
    Constructor for the :galacticus-class:`stellarPopulationSpectraPostprocessorSequence` stellar population spectra postprocessor class which takes a parameter set as input.
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
    !!{RST
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
    !!{RST
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
    !!{RST
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
    !!{RST
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

  subroutine sequenceAgeRange(self,ageMinimum,ageMaximum)
    !!{RST
    Return the range of ages retained by a sequence of postprocessors, which is the *intersection* of the ranges
    retained by its members---the multiplier of a sequence being the product of those of its members, emission
    survives only if it survives every one.

    Because each member's range is a single interval, and the intersection of intervals is an interval, the result can
    never contain a gap. Members whose ranges do not overlap at all intersect to nothing, which is reported as an
    empty range: **on return, ``ageMaximum`` :math:`\le` ``ageMinimum`` means that no age survives the sequence, and
    consumers must test for this** rather than treating the result as a narrow but usable window. Such a sequence
    suppresses all emission at every age, and is almost always a misconfiguration.
    !!}
    implicit none
    class           (stellarPopulationSpectraPostprocessorSequence), intent(inout) :: self
    double precision                                               , intent(  out) :: ageMinimum      , ageMaximum
    type            (postprocessorList                            ), pointer       :: postprocessor_
    double precision                                                               :: ageMinimumMember, ageMaximumMember

    ageMinimum     =  0.0d0
    ageMaximum     =  huge(0.0d0)
    postprocessor_ => self%postprocessors
    do while (associated(postprocessor_))
       call postprocessor_%postprocessor_%ageRange(ageMinimumMember,ageMaximumMember)
       ageMinimum     =  max(ageMinimum,ageMinimumMember)
       ageMaximum     =  min(ageMaximum,ageMaximumMember)
       postprocessor_ => postprocessor_%next
    end do
    ! An empty intersection retains nothing. Report it as a zero-width range rather than an inverted one, so that the
    ! test `ageMaximum <= ageMinimum` identifies the empty case whichever way it arose.
    if (ageMaximum < ageMinimum) ageMaximum=ageMinimum
    return
  end subroutine sequenceAgeRange

  logical function sequenceAgeWindowIsSharp(self) result(ageWindowIsSharp)
    !!{RST
    Return true only if every member of the sequence has a sharp age window.

    Sharpness composes: writing a sharp member's multiplier as :math:`c_i(\lambda,z)` within its window and zero
    outside, the product over members is :math:`\prod_i c_i(\lambda,z)` throughout the intersection---still
    independent of age---and zero outside it, since at least one factor vanishes there. Note that "sharp" means
    independent of age within the window, *not* equal to unity, so members with differing multipliers still compose to
    a sharp window; that is precisely why a wavelength-dependent but age-independent member such as
    :galacticus-class:`stellarPopulationSpectraPostprocessorInoue2014` cancels when two chains differing only in an
    age window are divided.

    The test is sufficient but not necessary, and so is conservative. A tapering member whose taper lies entirely
    outside the intersection contributes only its flat part, and the product is then sharp despite this method
    reporting otherwise---for example
    :galacticus-class:`stellarPopulationSpectraPostprocessorBirthCloudsLacey2016` combined with a
    :galacticus-class:`stellarPopulationSpectraPostprocessorRecent` whose limit is below one birth cloud lifetime.
    Erring this way causes a consumer to decline to split a luminosity it could in fact have split, rather than to
    split one it should not have.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorSequence), intent(inout) :: self
    type (postprocessorList                            ), pointer       :: postprocessor_

    ageWindowIsSharp =  .true.
    postprocessor_   => self%postprocessors
    do while (associated(postprocessor_))
       if (.not.postprocessor_%postprocessor_%ageWindowIsSharp()) then
          ageWindowIsSharp=.false.
          return
       end if
       postprocessor_ => postprocessor_%next
    end do
    return
  end function sequenceAgeWindowIsSharp
