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
Implements a merger tree random number seed in which the seed is chosen at random (without repetition) from the available range.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_long
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <mergerTreeSeeds name="mergerTreeSeedsRandom">
   <description>A merger tree random number seed in which the seed is chosen at random (without repetition) from the available range.</description>
  </mergerTreeSeeds>
  !!]
  type, extends(mergerTreeSeedsClass) :: mergerTreeSeedsRandom
     !!{
     A merger tree random number seed in which the seed is chosen at random (without repetition) from the available range.
     !!}
     private
     class  (randomNumberGeneratorClass), pointer                   :: randomNumberGenerator_ => null()
     integer(c_long                    ), allocatable, dimension(:) :: seeds
     integer(c_long                    )                            :: rangeMinimum                    , rangeMaximum
   contains
     final     ::        randomDestructor
     procedure :: set => randomSet
  end type mergerTreeSeedsRandom

  interface mergerTreeSeedsRandom
     !!{
     Constructors for the {\normalfont \ttfamily random} merger tree seed class.
     !!}
     module procedure randomConstructorParameters
     module procedure randomConstructorInternal
  end interface mergerTreeSeedsRandom

  integer(c_long), parameter :: seedMaximum=2_c_long**31
  
contains

  function randomConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily random} merger tree seed class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (mergerTreeSeedsRandom     )                :: self
    type (inputParameters           ), intent(inout) :: parameters
    class(randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    
    !![
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=mergerTreeSeedsRandom(randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function randomConstructorParameters

  function randomConstructorInternal(randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily random} merger tree seed class.
    !!}
    implicit none
    type (mergerTreeSeedsRandom     )                        :: self
    class(randomNumberGeneratorClass), intent(inout), target :: randomNumberGenerator_

    ! Make a copy of the random number generator. We want to ensure that we have a copy to ourself here, so that no other class
    ! can use it (and potentially de-synchronize the random number sequence between threads/processes).
    allocate(self%randomNumberGenerator_,mold=randomNumberGenerator_)
    !![
    <deepCopyReset    variables="randomNumberGenerator_"/>
    <deepCopy         source="randomNumberGenerator_" destination="self%randomNumberGenerator_"/>
    <deepCopyFinalize variables="self%randomNumberGenerator_"/>
    !!]    
    ! Check that we have non-independent random numbers across threads and processes.
    if (self%randomNumberGenerator_%mpiIndependent   ()) call Error_Report('random number generator produces different sequences on each MPI process - uniqueness of seeds can not be guaranteed'  //{introspection:location})
    if (self%randomNumberGenerator_%openMPIndependent()) call Error_Report('random number generator produces different sequences on each OpenMP thread - uniqueness of seeds can not be guaranteed'//{introspection:location})
    ! Check that the random number generator has sufficient range.
    self%rangeMinimum=self%randomNumberGenerator_%rangeMinimum()
    self%rangeMaximum=self%randomNumberGenerator_%rangeMaximum()
    if (self%rangeMaximum-self%rangeMinimum < seedMaximum) call Error_Report('random number generator has insufficient range to produce seeds'//{introspection:location})
    return
  end function randomConstructorInternal

  subroutine randomDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily random} merger tree seed class.
    !!}
    implicit none
    type(mergerTreeSeedsRandom), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine randomDestructor

  subroutine randomSet(self,tree)
    !!{
    Set the random number seed in the given {\normalfont \ttfamily tree}.
    !!}
    !$ use :: omp_lib, only : omp_get_thread_num
    implicit none
    class  (mergerTreeSeedsRandom), intent(inout)               :: self
    type   (mergerTree           ), intent(inout)               :: tree
    integer(c_long               ), allocatable  , dimension(:) :: seedsTmp
    integer(c_size_t             )                              :: countSeedsPrevious, i

    if (allocated(self%seeds)) then
       countSeedsPrevious=size(self%seeds)
       if (tree%index > size(self%seeds)) then
          call move_alloc(self%seeds,seedsTmp)
          allocate(self%seeds(tree%index))
          self%seeds(1:countSeedsPrevious)=seedsTmp
       end if
    else
       countSeedsPrevious=0_c_size_t
       allocate(self%seeds(tree%index))
    end if
    do i=countSeedsPrevious+1,tree%index
       if (i == 1) then
          ! This is the first seed - choose any number.
          self%seeds   (i)=self%randomNumberGenerator_%sample(seedMaximum)+1_c_long
       else
          ! For subsequent seeds check that we do not duplicate any prior seed.
          self%seeds(i)=self%seeds(i-1)
          do while (any(self%seeds(i) == self%seeds(1:i-1)))
             self%seeds(i)=self%randomNumberGenerator_%sample(seedMaximum)+1_c_long
          end do
       end if
    end do
    call tree%randomNumberGenerator_%seedSet(seed=self%seeds(tree%index),offset=.false.)
    return
  end subroutine randomSet
