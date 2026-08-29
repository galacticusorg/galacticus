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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a dust attenuation class which applies a sequence of other attenuators.
  !!}

  type, public :: dustAttenuationList
     class(dustAttenuationClass), pointer :: dustAttenuation_ => null()
     type (dustAttenuationList ), pointer :: next             => null()
  end type dustAttenuationList

  !![
  <dustAttenuation name="dustAttenuationSequence" docformat="rst">
   <description>
   A dust attenuation class which applies a sequence of other attenuators, returning the product of their
   transmissions. Because transmission is multiplicative, this is exactly the result of the light passing through each
   in turn, and---for pure absorption---the order in which they are listed does not matter.

   The canonical use is to build a two-component model: a :galacticus-class:`dustAttenuationBirthCloud` attenuating
   young populations, followed by a diffuse interstellar medium screen attenuating everything. See
   :galacticus-class:`dustAttenuationCharlotFall2000` for that combination pre-assembled.

   The decomposition requested by a sequence is the union of those requested by its members: the age bin boundaries
   are merged, and an axis is resolved if any member depends upon it. A component is supported only if every member
   supports it.
   </description>
   <linkedList type="dustAttenuationList" variable="dustAttenuations" next="next" object="dustAttenuation_" objectType="dustAttenuationClass"/>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationSequence
     !!{RST
     A dust attenuation class which applies a sequence of other attenuators.
     !!}
     private
     type(dustAttenuationList), pointer :: dustAttenuations => null()
   contains
     final     ::                      sequenceDestructor
     procedure :: transmission      => sequenceTransmission
     procedure :: request           => sequenceRequest
     procedure :: supportsComponent => sequenceSupportsComponent
  end type dustAttenuationSequence

  interface dustAttenuationSequence
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationSequence` dust attenuation class.
     !!}
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface dustAttenuationSequence

contains

  function sequenceConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationSequence` dust attenuation class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (dustAttenuationSequence)                :: self
    type   (inputParameters        ), intent(inout) :: parameters
    type   (dustAttenuationList    ), pointer       :: dustAttenuation_
    integer                                         :: i

    self            %dustAttenuations => null()
    dustAttenuation_                  => null()
    do i=1,parameters%copiesCount('dustAttenuation',zeroIfNotPresent=.true.)
       if (associated(dustAttenuation_)) then
          allocate(dustAttenuation_%next)
          dustAttenuation_ => dustAttenuation_%next
       else
          allocate(self%dustAttenuations)
          dustAttenuation_ => self%dustAttenuations
       end if
       !![
       <objectBuilder class="dustAttenuation" name="dustAttenuation_%dustAttenuation_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="dustAttenuation"/>
    !!]
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(dustAttenuations) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationSequence` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationSequence)                        :: self
    type(dustAttenuationList    ), intent(in   ), target :: dustAttenuations
    type(dustAttenuationList    ), pointer               :: dustAttenuation_

    self            %dustAttenuations => dustAttenuations
    dustAttenuation_                  => dustAttenuations
    do while (associated(dustAttenuation_))
       !![
       <referenceCountIncrement owner="dustAttenuation_" object="dustAttenuation_"/>
       !!]
       dustAttenuation_ => dustAttenuation_%next
    end do
    return
  end function sequenceConstructorInternal

  subroutine sequenceDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationSequence` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationSequence), intent(inout) :: self
    type(dustAttenuationList    ), pointer       :: dustAttenuation_, dustAttenuationNext

    if (associated(self%dustAttenuations)) then
       dustAttenuation_ => self%dustAttenuations
       do while (associated(dustAttenuation_))
          dustAttenuationNext => dustAttenuation_%next
          !![
          <objectDestructor name="dustAttenuation_%dustAttenuation_"/>
          !!]
          deallocate(dustAttenuation_)
          dustAttenuation_ => dustAttenuationNext
       end do
    end if
    return
  end subroutine sequenceDestructor

  function sequenceTransmission(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the product of the transmissions of each attenuator in the sequence.

    Any ``inclination`` is passed on to each member, so that wrapping a sequence in
    :galacticus-class:`dustAttenuationInclinationAveraged` averages the sequence as a whole---the product evaluated
    at each orientation---rather than the product of separately averaged members, which is not the same thing.
    !!}
    implicit none
    class           (dustAttenuationSequence), intent(inout)                               :: self
    type            (treeNode               ), intent(inout), target                       :: node
    type            (emissionDescriptor     ), intent(in   ), dimension(:)                 :: descriptors
    double precision                         , intent(in   ), optional                     :: inclination
    double precision                                        , dimension(size(descriptors)) :: transmission
    type            (dustAttenuationList    ), pointer                                     :: dustAttenuation_

    transmission     =  1.0d0
    dustAttenuation_ => self%dustAttenuations
    do while (associated(dustAttenuation_))
       transmission     =  +transmission                                                                 &
            &              *dustAttenuation_%dustAttenuation_%transmission(node,descriptors,inclination)
       dustAttenuation_ =>  dustAttenuation_%next
    end do
    return
  end function sequenceTransmission

  function sequenceRequest(self) result(request)
    !!{RST
    Return the union of the decomposition requests of each attenuator in the sequence: age bin boundaries are merged
    (sorted, with duplicates removed), and each resolution flag is set if any member sets it.
    !!}
    implicit none
    type            (decompositionRequest   )                            :: request
    class           (dustAttenuationSequence), intent(inout)             :: self
    type            (dustAttenuationList    ), pointer                   :: dustAttenuation_
    type            (decompositionRequest   )                            :: requestMember
    double precision                         , allocatable, dimension(:) :: boundaries      , boundariesMerged
    integer                                                              :: countBoundaries , i               , &
         &                                                                  j
    double precision                                                     :: boundary

    request%resolveComponents =.false.
    request%resolveMetallicity=.false.
    request%resolveRadius     =.false.
    allocate(boundaries(0))
    dustAttenuation_ => self%dustAttenuations
    do while (associated(dustAttenuation_))
       requestMember=dustAttenuation_%dustAttenuation_%request()
       request%resolveComponents =request%resolveComponents  .or. requestMember%resolveComponents
       request%resolveMetallicity=request%resolveMetallicity .or. requestMember%resolveMetallicity
       request%resolveRadius     =request%resolveRadius      .or. requestMember%resolveRadius
       if (allocated(requestMember%ageBoundaries)) then
          boundariesMerged=[boundaries,requestMember%ageBoundaries]
          call move_alloc(boundariesMerged,boundaries)
       end if
       dustAttenuation_ => dustAttenuation_%next
    end do
    ! Sort the merged boundaries, then remove duplicates. The arrays here contain at most a handful of entries, so an
    ! insertion sort is entirely adequate and avoids pulling in a general sorting routine.
    do i=2,size(boundaries)
       boundary=boundaries(i)
       j       =i-1
       do while (j >= 1)
          if (boundaries(j) <= boundary) exit
          boundaries(j+1)=boundaries(j)
          j              =j-1
       end do
       boundaries(j+1)=boundary
    end do
    countBoundaries=0
    do i=1,size(boundaries)
       if (countBoundaries == 0) then
          countBoundaries            =1
          boundaries(countBoundaries)=boundaries(i)
       else if (boundaries(i) > boundaries(countBoundaries)) then
          countBoundaries            =countBoundaries+1
          boundaries(countBoundaries)=boundaries(i)
       end if
    end do
    request%ageBoundaries=boundaries(1:countBoundaries)
    return
  end function sequenceRequest

  logical function sequenceSupportsComponent(self,componentType) result(supportsComponent)
    !!{RST
    Return true only if every attenuator in the sequence supports the given component.
    !!}
    implicit none
    class(dustAttenuationSequence     ), intent(inout) :: self
    type (enumerationComponentTypeType), intent(in   ) :: componentType
    type (dustAttenuationList         ), pointer       :: dustAttenuation_

    supportsComponent=  .true.
    dustAttenuation_  => self%dustAttenuations
    do while (associated(dustAttenuation_))
       if (.not.dustAttenuation_%dustAttenuation_%supportsComponent(componentType)) then
          supportsComponent=.false.
          return
       end if
       dustAttenuation_ => dustAttenuation_%next
    end do
    return
  end function sequenceSupportsComponent
