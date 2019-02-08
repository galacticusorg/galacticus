!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !# <outputAnalysisDistributionNormalizer name="outputAnalysisDistributionNormalizerSequence">
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
     procedure :: deepCopy   => sequenceDeepCopy
  end type outputAnalysisDistributionNormalizerSequence

  interface outputAnalysisDistributionNormalizerSequence
     !% Constructors for the sequence on-the-fly output distribution normalizer class.
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface outputAnalysisDistributionNormalizerSequence

contains

  function sequenceConstructorParameters(parameters) result(self)
    !% Constructor for the sequence on-the-fly output normalizer class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type   (outputAnalysisDistributionNormalizerSequence)                :: self
    type   (inputParameters                             ), intent(inout) :: parameters
    type   (normalizerList                              ), pointer       :: normalizer_
    integer                                                              :: i

    self       %normalizers => null()
    normalizer_             => null()
    do i=1,parameters%copiesCount('outputAnalysisDistributionNormalizerMethod',zeroIfNotPresent=.true.)
       if (associated(normalizer_)) then
          allocate(normalizer_%next)
          normalizer_ => normalizer_%next
       else
          allocate(self%normalizers)
          normalizer_ => self%normalizers
       end if
       !# <objectBuilder class="outputAnalysisDistributionNormalizer" name="normalizer_%normalizer_" source="parameters" copy="i" />
    end do
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(normalizers) result(self)
    !% Internal constructor for the sequence merger tree normalizer class.
    implicit none
    type(outputAnalysisDistributionNormalizerSequence)                        :: self
    type(normalizerList                              ), target, intent(in   ) :: normalizers
    type(normalizerList                              ), pointer               :: normalizer_

    self       %normalizers => normalizers
    normalizer_             => normalizers
    do while (associated(normalizer_))
       !# <referenceCountIncrement owner="normalizer_" object="normalizer_"/>
       normalizer_ => normalizer_%next
    end do
    return
  end function sequenceConstructorInternal

  subroutine sequenceDestructor(self)
    !% Destructor for the merger tree normalizer function class.
    implicit none
    type(outputAnalysisDistributionNormalizerSequence), intent(inout) :: self
    type(normalizerList                              ), pointer       :: normalizer_, normalizerNext

    if (associated(self%normalizers)) then
       normalizer_ => self%normalizers
       do while (associated(normalizer_))
          normalizerNext => normalizer_%next
          !# <objectDestructor name="normalizer_%normalizer_"/>
          deallocate(normalizer_)
          normalizer_ => normalizerNext
       end do
    end if
    return
  end subroutine sequenceDestructor

  subroutine sequenceNormalize(self,distribution,covariance,propertyValueMinimum,propertyValueMaximum)
    !% Perform a sequence normalization on an on-the-fly output distribution.
    implicit none
    class           (outputAnalysisDistributionNormalizerSequence), intent(inout)                 :: self
    double precision                                              , intent(inout), dimension(:  ) :: distribution
    double precision                                              , intent(inout), dimension(:,:) :: covariance
    double precision                                              , intent(in   ), dimension(:  ) :: propertyValueMinimum, propertyValueMaximum
    type            (normalizerList                              ), pointer                       :: normalizer_

    normalizer_ => self%normalizers
    do while (associated(normalizer_))
       call normalizer_%normalizer_%normalize(distribution,covariance,propertyValueMinimum,propertyValueMaximum)
       normalizer_ => normalizer_%next
    end do
    return
  end subroutine sequenceNormalize

  subroutine sequenceDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily sequence} output analysis distribution normalizer operator class.
    use Galacticus_Error
    implicit none
    class(outputAnalysisDistributionNormalizerSequence), intent(inout) :: self
    class(outputAnalysisDistributionNormalizerClass   ), intent(inout) :: destination
    type (normalizerList                              ), pointer       :: normalizer_   , normalizerDestination_, &
         &                                                                normalizerNew_

    call self%outputAnalysisDistributionNormalizerClass%deepCopy(destination)
    select type (destination)
    type is (outputAnalysisDistributionNormalizerSequence)
       destination%normalizers => null          ()
       normalizerDestination_  => null          ()
       normalizer_             => self%normalizers
       do while (associated(normalizer_))
          allocate(normalizerNew_)
          if (associated(normalizerDestination_)) then
             normalizerDestination_%next       => normalizerNew_
             normalizerDestination_            => normalizerNew_             
          else
             destination          %normalizers => normalizerNew_
             normalizerDestination_            => normalizerNew_
          end if
          allocate(normalizerNew_%normalizer_,mold=normalizer_%normalizer_)
          !# <deepCopy source="normalizer_%normalizer_" destination="normalizerNew_%normalizer_"/>
          normalizer_ => normalizer_%next
       end do       
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine sequenceDeepCopy
