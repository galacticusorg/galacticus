!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% Contains a module which implements a sequence of normalizers on on-the-fly outputs.

  !# <outputAnalysisDistributionNormalizer name="outputAnalysisDistributionNormalizerSequence" defaultThreadPrivate="yes">
  !#  <description>Provides a sequence of normalizers on on-the-fly outputs.</description>
  !# </outputAnalysisDistributionNormalizer>

  type, public :: normalizerList
     class(outputAnalysisDistributionNormalizerClass), pointer :: normalizer_
     type (normalizerList                           ), pointer :: next        => null()
  end type normalizerList

  type, extends(outputAnalysisDistributionNormalizerClass) :: outputAnalysisDistributionNormalizerSequence
     !% A sequence on-the-fly-output normalizer class.
     private
     type(normalizerList), pointer :: normalizers
  contains
     final     ::               sequenceDestructor
     procedure :: normalize  => sequenceNormalize
  end type outputAnalysisDistributionNormalizerSequence

  interface outputAnalysisDistributionNormalizerSequence
     !% Constructors for the sequence on-the-fly output distribution normalizer class.
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface outputAnalysisDistributionNormalizerSequence

contains

  function sequenceConstructorParameters(parameters)
    !% Constructor for the sequence on-the-fly output normalizer class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type   (outputAnalysisDistributionNormalizerSequence)                :: sequenceConstructorParameters
    type   (inputParameters                             ), intent(inout) :: parameters
    type   (normalizerList                              ), pointer       :: normalizer_
    integer                                                              :: i

    sequenceConstructorParameters%normalizers => null()
    normalizer_                               => null()
    do i=1,parameters%copiesCount('outputAnalysisDistributionNormalizerMethod',zeroIfNotPresent=.true.)
       if (associated(normalizer_)) then
          allocate(normalizer_%next)
          normalizer_ => normalizer_%next
       else
          allocate(sequenceConstructorParameters%normalizers)
          normalizer_ => sequenceConstructorParameters%normalizers
       end if
       normalizer_%normalizer_ => outputAnalysisDistributionNormalizer(parameters,i)
    end do
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(normalizers)
    !% Internal constructor for the sequence merger tree normalizer class.
    implicit none
    type(outputAnalysisDistributionNormalizerSequence)                        :: sequenceConstructorInternal
    type(normalizerList                              ), target, intent(in   ) :: normalizers

    sequenceConstructorInternal%normalizers => normalizers
    return
  end function sequenceConstructorInternal

  elemental subroutine sequenceDestructor(self)
    !% Destructor for the merger tree normalizer function class.
    implicit none
    type(outputAnalysisDistributionNormalizerSequence), intent(inout) :: self
    type(normalizerList                              ), pointer       :: normalizer_, normalizerNext

    if (associated(self%normalizers)) then
       normalizer_ => self%normalizers
       do while (associated(normalizer_))
          normalizerNext => normalizer_%next
          deallocate(normalizer_%normalizer_)
          deallocate(normalizer_            )
          normalizer_ => normalizerNext
       end do
    end if
    return
  end subroutine sequenceDestructor

  function sequenceNormalize(self,distribution,propertyValueMinimum,propertyValueMaximum)
    !% Perform a sequence normalization on an on-the-fly output distribution.
    implicit none
    class           (outputAnalysisDistributionNormalizerSequence), intent(inout)                                        :: self
    double precision                                              , intent(in   ), dimension(:)                          :: distribution
    double precision                                              , intent(in   ), dimension(:)                          :: propertyValueMinimum, propertyValueMaximum
    double precision                                                             , dimension(size(propertyValueMinimum)) :: sequenceNormalize
    type            (normalizerList                              ), pointer                                              :: normalizer_

    sequenceNormalize =  distribution
    normalizer_       => self%normalizers
    do while (associated(normalizer_))
       sequenceNormalize=normalizer_%normalizer_%normalize(sequenceNormalize,propertyValueMinimum,propertyValueMaximum)
       normalizer_ => normalizer_%next
    end do
    return
  end function sequenceNormalize
