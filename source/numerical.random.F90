!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements pseudo-random numbers.

module Pseudo_Random
  !% Implements pseudo-random numbers.
  use FGSL
  implicit none
  private
  public :: Pseudo_Random_Get, Pseudo_Random_Free, Pseudo_Random_Store, Pseudo_Random_Retrieve

  type, public :: pseudoRandom
     !% Wrapper class for pseudo random sequence.
     private
     type   (fgsl_rng) :: pseudoSequence
     logical           :: pseudoSequenceReset=.true.
   contains
     !# <workaround type="gfortran" PR="58471 58470" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58471 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58470">
     !# final     ::           pseudoRandomDestructor
     !# </workaround>
     !@ <objectMethods>
     !@   <object>pseudoRandom</object>
     !@   <objectMethod>
     !@     <method>sample</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Return a pseudo-random number.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: sample => pseudoRandomSample
  end type pseudoRandom
  
  logical                 :: Seed_Is_Set=.false.
  integer                 :: randomSeed
  integer(kind=fgsl_long) :: randomSeedC=-1
  !$omp threadprivate(randomSeedC)

contains

  double precision function Pseudo_Random_Get(pseudoSequenceObject,reset,ompThreadOffset,mpiRankOffset,incrementSeed)
    !% Returns a scalar giving a pseudo-random number.
    use Input_Parameters
    use MPI
    !$ use OMP_Lib
    implicit none
    type   (fgsl_rng), intent(inout)           :: pseudoSequenceObject
    logical          , intent(inout), optional :: reset,ompThreadOffset
    logical          , intent(in   ), optional :: ompThreadOffset     , mpiRankOffset
    integer          , intent(in   ), optional :: incrementSeed
    logical                                    :: resetActual
    integer                                    :: mpiRank             , iError

    ! Determine if we need to reset.
    if (present(reset)) then
       resetActual=reset
       reset      =.false.
    else
       resetActual=.false.
    end if

    ! Read in the random number seed if necessary.
    !$omp critical (Pseudo_Random_Get)
    if (.not.Seed_Is_Set) then
       !@ <inputParameter>
       !@   <name>randomSeed</name>
       !@   <defaultValue>219</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     A seed value for the random number generator.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('randomSeed',randomSeed,defaultValue=219)
       Seed_Is_Set=.true.
    end if
    !$omp end critical (Pseudo_Random_Get)

    if (resetActual.or..not.FGSL_Well_Defined(pseudoSequenceObject)) then
       pseudoSequenceObject=FGSL_RNG_Alloc(FGSL_RNG_Default)
       if (randomSeedC < 0 .or. present(incrementSeed)) then
          randomSeedC=randomSeed
          !$ if (present(ompThreadOffset)) then
          !$    if (ompThreadOffset) randomSeedC=randomSeedC+omp_get_thread_num()
          !$ end if
          if (present(mpiRankOffset).and.mpiRankOffset) then
             call MPI_Comm_Rank(MPI_Comm_World,mpiRank,iError)
             randomSeedC=randomSeedC+mpiRank
          end if
          if (present(incrementSeed)) randomSeedC=randomSeedC+incrementSeed
       end if
       randomSeedC=randomSeedC+1
       call FGSL_RNG_Set(pseudoSequenceObject,randomSeedC)
    end if
    Pseudo_Random_Get=FGSL_RNG_Uniform(pseudoSequenceObject)
    return
  end function Pseudo_Random_Get

  subroutine Pseudo_Random_Free(pseudoSequenceObject)
    !% Frees a pseudo-random sequence object.
    implicit none
    type(fgsl_rng), intent(inout) :: pseudoSequenceObject

    call FGSL_RNG_Free(pseudoSequenceObject)
    return
  end subroutine Pseudo_Random_Free

  subroutine Pseudo_Random_Store(pseudoSequenceObject,fgslFile)
    !% Stores a pseudo-random sequence object to file.
    implicit none
    type   (fgsl_rng ), intent(in   ) :: pseudoSequenceObject
    type   (fgsl_file), intent(in   ) :: fgslFile
    integer                           :: iError

    iError=FGSL_Rng_FWrite(fgslFile,pseudoSequenceObject)
    return
  end subroutine Pseudo_Random_Store

  subroutine Pseudo_Random_Retrieve(pseudoSequenceObject,fgslFile)
    !% Stores a pseudo-random sequence object to file.
    implicit none
    type   (fgsl_rng ), intent(inout) :: pseudoSequenceObject
    type   (fgsl_file), intent(in   ) :: fgslFile
    integer                           :: iError

    if (.not.FGSL_Well_Defined(pseudoSequenceObject)) pseudoSequenceObject=FGSL_RNG_Alloc(FGSL_RNG_Default)
    iError=FGSL_Rng_FRead(fgslFile,pseudoSequenceObject)
    return
  end subroutine Pseudo_Random_Retrieve

  subroutine pseudoRandomDestructor(self)
    !% Destroy the pseudo-sequence wrapper classs.
    implicit none
    type(pseudoRandom), intent(inout) :: self
    
    if (.not.self%pseudoSequenceReset) call Pseudo_Random_Free(self%pseudoSequence)
    return
  end subroutine pseudoRandomDestructor

  double precision function pseudoRandomSample(self,ompThreadOffset,mpiRankOffset,incrementSeed)
    !% Sample from a pseudo-random sequence.
    implicit none
    class  (pseudoRandom), intent(inout)           :: self
    logical              , intent(in   ), optional :: ompThreadOffset     , mpiRankOffset
    integer              , intent(in   ), optional :: incrementSeed

    pseudoRandomSample=Pseudo_Random_Get(                          &
         &                               self%pseudoSequence     , &
         &                               self%pseudoSequenceReset, &
         &                               ompThreadOffset         , &
         &                               mpiRankOffset           , &
         &                               incrementSeed             &
         &                              )
    return
  end function pseudoRandomSample

end module Pseudo_Random
