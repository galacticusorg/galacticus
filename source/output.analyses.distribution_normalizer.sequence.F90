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
  Implements a sequence of normalizers on on-the-fly outputs.
  !!}

  !![
  <outputAnalysisDistributionNormalizer name="outputAnalysisDistributionNormalizerSequence">
   <description>Provides a sequence of normalizers on on-the-fly outputs.</description>
   <linkedList type="normalizerList" variable="normalizers" next="next" object="normalizer_" objectType="outputAnalysisDistributionNormalizerClass"/>
  </outputAnalysisDistributionNormalizer>
  !!]

  type, public :: normalizerList
     class(outputAnalysisDistributionNormalizerClass), pointer :: normalizer_ => null()
     type (normalizerList                           ), pointer :: next        => null()
  end type normalizerList

  type, extends(outputAnalysisDistributionNormalizerClass) :: outputAnalysisDistributionNormalizerSequence
     !!{
     A sequence on-the-fly-output normalizer class.
     !!}
     private
     type(normalizerList), pointer :: normalizers => null()
  contains
     final     ::              sequenceDestructor
     procedure :: normalize => sequenceNormalize
  end type outputAnalysisDistributionNormalizerSequence

  interface outputAnalysisDistributionNormalizerSequence
     !!{
     Constructors for the sequence on-the-fly output distribution normalizer class.
     !!}
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface outputAnalysisDistributionNormalizerSequence

contains

  function sequenceConstructorParameters(parameters) result(self)
    !!{
    Constructor for the sequence on-the-fly output normalizer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (outputAnalysisDistributionNormalizerSequence)                :: self
    type   (inputParameters                             ), intent(inout) :: parameters
    type   (normalizerList                              ), pointer       :: normalizer_
    integer                                                              :: i

    self       %normalizers => null()
    normalizer_             => null()
    do i=1,parameters%copiesCount('outputAnalysisDistributionNormalizer',zeroIfNotPresent=.true.)
       if (associated(normalizer_)) then
          allocate(normalizer_%next)
          normalizer_ => normalizer_%next
       else
          allocate(self%normalizers)
          normalizer_ => self%normalizers
       end if
       !![
       <objectBuilder class="outputAnalysisDistributionNormalizer" name="normalizer_%normalizer_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="outputAnalysisDistributionNormalizer"/>
    !!]
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(normalizers) result(self)
    !!{
    Internal constructor for the sequence merger tree normalizer class.
    !!}
    implicit none
    type(outputAnalysisDistributionNormalizerSequence)                        :: self
    type(normalizerList                              ), target, intent(in   ) :: normalizers
    type(normalizerList                              ), pointer               :: normalizer_

    self       %normalizers => normalizers
    normalizer_             => normalizers
    do while (associated(normalizer_))
       !![
       <referenceCountIncrement owner="normalizer_" object="normalizer_"/>
       !!]
       normalizer_ => normalizer_%next
    end do
    return
  end function sequenceConstructorInternal

  subroutine sequenceDestructor(self)
    !!{
    Destructor for the merger tree normalizer function class.
    !!}
    implicit none
    type(outputAnalysisDistributionNormalizerSequence), intent(inout) :: self
    type(normalizerList                              ), pointer       :: normalizer_, normalizerNext

    if (associated(self%normalizers)) then
       normalizer_ => self%normalizers
       do while (associated(normalizer_))
          normalizerNext => normalizer_%next
          !![
          <objectDestructor name="normalizer_%normalizer_"/>
          !!]
          deallocate(normalizer_)
          normalizer_ => normalizerNext
       end do
    end if
    return
  end subroutine sequenceDestructor

  subroutine sequenceNormalize(self,distribution,covariance,propertyValueMinimum,propertyValueMaximum)
    !!{
    Perform a sequence normalization on an on-the-fly output distribution.
    !!}
    implicit none
    class           (outputAnalysisDistributionNormalizerSequence), intent(inout)                           :: self
    double precision                                              , intent(inout), dimension(:  ), optional :: distribution
    double precision                                              , intent(inout), dimension(:,:), optional :: covariance
    double precision                                              , intent(in   ), dimension(:  )           :: propertyValueMinimum, propertyValueMaximum
    type            (normalizerList                              ), pointer                                 :: normalizer_

    normalizer_ => self%normalizers
    do while (associated(normalizer_))
       call normalizer_%normalizer_%normalize(distribution,covariance,propertyValueMinimum,propertyValueMaximum)
       normalizer_ => normalizer_%next
    end do
    return
  end subroutine sequenceNormalize
