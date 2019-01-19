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

!% Contains a module which implements pseudo-random numbers.

module Pseudo_Random
  !% Implements pseudo-random numbers.
  use FGSL, only : fgsl_rng        , fgsl_long       , FGSL_Well_Defined, FGSL_RNG_Alloc, &
       &           FGSL_RNG_Default, FGSL_RNG_Set    , FGSL_RNG_Free    , FGSL_Obj_C_Ptr, &
       &           FGSL_RNG_Uniform, FGSL_Ran_Poisson, FGSL_Ran_Gaussian, FGSL_Rng_FRead, &
       &           FGSL_Rng_Clone  , fgsl_file       , FGSL_Rng_FWrite
  implicit none
  private

  type, public :: pseudoRandom
     !% Wrapper class for pseudo random sequence.
     private
     type   (fgsl_rng) :: pseudoSequence
     logical           :: pseudoSequenceReset=.true.
   contains
     final     ::           pseudoRandomDestructor
     !@ <objectMethods>
     !@   <object>pseudoRandom</object>
     !@   <objectMethod>
     !@     <method>uniformSample</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Return a pseudo-random number drawn from a uniform distribution on the interval 0 to 1.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>poissonSample</method>
     !@     <type>\intzero</type>
     !@     <arguments>\doublezero\ mean</arguments>
     !@     <description>Return a pseudo-random number drawn from a Poisson distribution of the given mean.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>normalSample</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ mean</arguments>
     !@     <description>Return a pseudo-random number drawn from a standard normal distribution.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>clone</method>
     !@     <type>\textless type(pseudoRandom)\textgreater</type>
     !@     <arguments></arguments>
     !@     <description>Clone a pseudo-random number generator.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>store</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ stateFile\argin, \textless type(fgsl\_file)\textgreater\ fgslStateFile\argin</arguments>
     !@     <description>Store a pseudo-random number generator state to file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>restore</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ stateFile\argin, \textless type(fgsl\_file)\textgreater\ fgslStateFile\argin</arguments>
     !@     <description>Restore a pseudo-random number generator state from file.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: uniformSample => pseudoRandomUniformSample
     procedure :: poissonSample => pseudoRandomPoissonSample
     procedure :: normalSample  => pseudoRandomNormalSample
     procedure :: clone         => pseudoRandomClone
     procedure :: store         => pseudoRandomStore
     procedure :: restore       => pseudoRandomRestore
  end type pseudoRandom

  interface pseudoRandom
     !% Constructors for the pseudoRandom class.
     module procedure pseudoRandomConstructor
  end interface pseudoRandom
  
  logical                 :: Seed_Is_Set=.false.
  integer                 :: randomSeed
  integer(kind=fgsl_long) :: randomSeedC=-1
  !$omp threadprivate(randomSeedC)

contains

  subroutine pseudoRandomInitialize()
    !% Initialize pseudo-random number sequence generators.
    use Input_Parameters
    implicit none

    ! Read in the random number seed if necessary.
    if (.not.Seed_Is_Set) then
       !$omp critical (pseudoRandomInitialize)
       if (.not.Seed_Is_Set) then
          !# <inputParameter>
          !#   <name>randomSeed</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>219</defaultValue>
          !#   <description>A seed value for the random number generator.</description>
          !#   <source>globalParameters</source>
          !#   <type>integer</type>
          !# </inputParameter>
          Seed_Is_Set=.true.
       end if
       !$omp end critical (pseudoRandomInitialize)
    end if
    return
  end subroutine pseudoRandomInitialize

  subroutine pseudoRandomReset(pseudoSequenceObject,reset,ompThreadOffset,mpiRankOffset,incrementSeed)
    !$ use OMP_Lib
#ifdef USEMPI
    use MPI
#endif
    implicit none
#ifdef USEMPI
    integer                                    :: mpiRank             , iError
#endif
    type   (fgsl_rng), intent(inout)           :: pseudoSequenceObject
    logical          , intent(inout), optional :: reset
    logical          , intent(in   ), optional :: ompThreadOffset     , mpiRankOffset
    integer          , intent(in   ), optional :: incrementSeed
    !# <optionalArgument name="reset" defaultsTo=".false." />

    ! Determine if we need to reset.
    if (present(reset)) reset=.false.
    if (reset_.or..not.FGSL_Well_Defined(pseudoSequenceObject)) then
       pseudoSequenceObject=FGSL_RNG_Alloc(FGSL_RNG_Default)
       if (randomSeedC < 0 .or. present(incrementSeed)) then
          randomSeedC=randomSeed
          !$ if (present(ompThreadOffset)) then
          !$    if (ompThreadOffset) randomSeedC=randomSeedC+OMP_Get_Thread_Num()
          !$ end if
          if (present(mpiRankOffset).and.mpiRankOffset) then
#ifdef USEMPI
             call MPI_Comm_Rank(MPI_Comm_World,mpiRank,iError)
             randomSeedC=randomSeedC+mpiRank
#endif
          end if
          if (present(incrementSeed)) randomSeedC=randomSeedC+incrementSeed
       end if
       randomSeedC=randomSeedC+1
       call FGSL_RNG_Set(pseudoSequenceObject,randomSeedC)
    end if
    return
  end subroutine pseudoRandomReset

  subroutine pseudoRandomFree(pseudoSequenceObject)
    !% Frees a pseudo-random sequence object.
    implicit none
    type(fgsl_rng), intent(inout) :: pseudoSequenceObject

    call FGSL_RNG_Free(pseudoSequenceObject)
    return
  end subroutine pseudoRandomFree

  function pseudoRandomConstructor() result(self)
    !% Construct a pseudo-random sequence object.
    use, intrinsic :: ISO_C_Binding
    implicit none
    type(pseudoRandom) :: self

    self%pseudoSequenceReset=.true.
    call FGSL_Obj_C_Ptr(self%pseudoSequence,C_Null_Ptr)
    return
  end function pseudoRandomConstructor

  subroutine pseudoRandomDestructor(self)
    !% Destroy the pseudo-sequence wrapper classs.
    implicit none
    type(pseudoRandom), intent(inout) :: self
    
    if (.not.self%pseudoSequenceReset) call pseudoRandomFree(self%pseudoSequence)
    return
  end subroutine pseudoRandomDestructor

  double precision function pseudoRandomUniformSample(self,ompThreadOffset,mpiRankOffset,incrementSeed)
    !% Sample from a uniform distribution on the interval 0 to 1.
    implicit none
    class  (pseudoRandom), intent(inout)           :: self
    logical              , intent(in   ), optional :: ompThreadOffset, mpiRankOffset
    integer              , intent(in   ), optional :: incrementSeed

    call pseudoRandomInitialize(                                                                                        )
    call pseudoRandomReset     (self%pseudoSequence,self%pseudoSequenceReset,ompThreadOffset,mpiRankOffset,incrementSeed)
    pseudoRandomUniformSample=FGSL_RNG_Uniform(self%pseudoSequence)

    return
  end function pseudoRandomUniformSample

  integer function pseudoRandomPoissonSample(self,mean,ompThreadOffset,mpiRankOffset,incrementSeed)
    !% Sample from a Poisson distribution with the given mean.
    implicit none
    class           (pseudoRandom), intent(inout)           :: self
    double precision              , intent(in   )           :: mean
    logical                       , intent(in   ), optional :: ompThreadOffset, mpiRankOffset
    integer                       , intent(in   ), optional :: incrementSeed

    call pseudoRandomInitialize(                                                                                        )
    call pseudoRandomReset     (self%pseudoSequence,self%pseudoSequenceReset,ompThreadOffset,mpiRankOffset,incrementSeed)
    pseudoRandomPoissonSample=FGSL_Ran_Poisson(self%pseudoSequence,mean)
    return
  end function pseudoRandomPoissonSample

  double precision function pseudoRandomNormalSample(self,ompThreadOffset,mpiRankOffset,incrementSeed)
    !% Sample from a standard normal distribution.
    implicit none
    class           (pseudoRandom), intent(inout)           :: self
    logical                       , intent(in   ), optional :: ompThreadOffset, mpiRankOffset
    integer                       , intent(in   ), optional :: incrementSeed

    call pseudoRandomInitialize(                                                                                        )
    call pseudoRandomReset     (self%pseudoSequence,self%pseudoSequenceReset,ompThreadOffset,mpiRankOffset,incrementSeed)
    pseudoRandomNormalSample=FGSL_Ran_Gaussian(self%pseudoSequence,sigma=1.0d0)
    return
  end function pseudoRandomNormalSample

  function pseudoRandomClone(self)
    !% Clone a pseudo-random sequence object.
    implicit none
    type (pseudoRandom)                :: pseudoRandomClone
    class(pseudoRandom), intent(inout) :: self

    if (.not.self%pseudoSequenceReset) then
       if (FGSL_Well_Defined(pseudoRandomClone%pseudoSequence)) call pseudoRandomFree(pseudoRandomClone%pseudoSequence)
       pseudoRandomClone%pseudoSequence=FGSL_Rng_Clone(self%pseudoSequence)
    end if
    pseudoRandomClone%pseudoSequenceReset=self%pseudoSequenceReset
    return
  end function pseudoRandomClone

  subroutine pseudoRandomStore(self,stateFile,fgslStateFile)
    !% Store a pseudo-random sequence object state to file.
    implicit none
    class  (pseudoRandom), intent(inout) :: self
    integer              , intent(in   ) :: stateFile
    type   (fgsl_file   ), intent(in   ) :: fgslStateFile
    integer                              :: iStatus

    write (stateFile) self%pseudoSequenceReset
    if (.not.self%pseudoSequenceReset) iStatus=FGSL_Rng_FWrite(fgslStateFile,self%pseudoSequence)
    return
  end subroutine pseudoRandomStore

  subroutine pseudoRandomRestore(self,stateFile,fgslStateFile)
    !% Store a pseudo-random sequence object state to file.
    implicit none
    class  (pseudoRandom), intent(inout) :: self
    integer              , intent(in   ) :: stateFile
    type   (fgsl_file   ), intent(in   ) :: fgslStateFile
    integer                              :: iStatus

    read (stateFile) self%pseudoSequenceReset
    if (.not.self%pseudoSequenceReset) then
       if (.not.FGSL_Well_Defined(self%pseudoSequence)) self%pseudoSequence=FGSL_RNG_Alloc(FGSL_RNG_Default)
       iStatus=FGSL_Rng_FRead(fgslStateFile,self%pseudoSequence)
    end if
    return
  end subroutine pseudoRandomRestore

end module Pseudo_Random
