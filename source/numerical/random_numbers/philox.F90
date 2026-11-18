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
  Implements a counter-based random number generator class using the Philox4x64-10 algorithm.
  !!}

  ! Specify an explicit dependence on the philox_block.o object file.
  !: $(BUILDPATH)/numerical/random_numbers/philox_block.o

  use, intrinsic :: ISO_C_Binding, only : c_long, c_int64_t

  !![
  <randomNumberGenerator name="randomNumberGeneratorPhilox" docformat="rst">
   <description>
   A counter-based random number generator class using the Philox4x64-10 algorithm of :cite:`salmon_parallel_2011`. The
   generator has no internal state other than a 128-bit key and a 256-bit counter: each block of four 64-bit random integers is
   a fixed (bijective) function of the key and the counter value. The first word of the key is the seed, and the second holds
   the OpenMP thread number and MPI process rank if offsets by these are requested, so that the sequences of different threads
   and processes are guaranteed not to overlap. The output is identical to that of the Random123 library and of
   ``numpy.random.Philox``.

   As well as producing a conventional sequence of random numbers, the generator can be positioned at an arbitrary point by
   setting the upper three words of the counter. Each such setting gives an independent stream of :math:`2^{64}` blocks. A
   random number that is a function of some discrete labels (for example, of a node, a time interval and a vector component) can
   therefore be regenerated whenever it is needed, rather than stored.

   Uniform deviates use the upper 53 bits of each 64-bit integer. Gaussian deviates are drawn using the
   :cite:t:`box_note_1958` transform, and Poisson deviates by inversion (for means below 10) or by the transformed rejection
   method of :cite:t:`hormann_transformed_1993` (for larger means).
   </description>
  </randomNumberGenerator>
  !!]
  type, extends(randomNumberGeneratorClass) :: randomNumberGeneratorPhilox
     !!{RST
     A counter-based random number generator class using the Philox4x64-10 algorithm.
     !!}
     private
     integer         (c_int64_t), dimension(2) :: key
     integer         (c_int64_t), dimension(4) :: counter         , buffer
     integer                                   :: bufferPosition
     integer         (c_long   )               :: seed_           , seed__
     logical                                   :: ompThreadOffset , mpiRankOffset, &
          &                                       normalAvailable
     double precision                          :: normalStored
   contains
     !![
     <methods docformat="rst">
       <method method="counterSet" description="Position the generator at the start of the stream labeled by the upper three words of the counter."/>
     </methods>
     !!]
     procedure :: openMPIndependent    => philoxOpenMPIndependent
     procedure :: mpiIndependent       => philoxMPIIndependent
     procedure :: rangeMinimum         => philoxRangeMinimum
     procedure :: rangeMaximum         => philoxRangeMaximum
     procedure :: sample               => philoxSample
     procedure :: uniformSample        => philoxUniformSample
     procedure :: poissonSample        => philoxPoissonSample
     procedure :: standardNormalSample => philoxStandardNormalSample
     procedure :: seed                 => philoxSeed
     procedure :: seedSet              => philoxSeedSet
     procedure :: counterSet           => philoxCounterSet
  end type randomNumberGeneratorPhilox

  interface randomNumberGeneratorPhilox
     !!{RST
     Constructors for the :galacticus-class:`randomNumberGeneratorPhilox` random number generator class.
     !!}
     module procedure philoxConstructorParameters
     module procedure philoxConstructorInternal
  end interface randomNumberGeneratorPhilox

  interface
     subroutine philox4x64Next(counter,key,output) bind(c,name='philox4x64Next')
       !!{RST
       Template for the C function which increments the counter and returns the Philox4x64-10 block for the new counter value.
       !!}
       import c_int64_t
       integer(c_int64_t), dimension(4), intent(inout) :: counter
       integer(c_int64_t), dimension(2), intent(in   ) :: key
       integer(c_int64_t), dimension(4), intent(  out) :: output
     end subroutine philox4x64Next
  end interface

contains

  function philoxConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`randomNumberGeneratorPhilox` random number generator class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (randomNumberGeneratorPhilox)                :: self
    type   (inputParameters            ), intent(inout) :: parameters
    integer(c_long                     )                :: seed
    logical                                             :: ompThreadOffset, mpiRankOffset

    !![
    <inputParameter docformat="rst">
      <name>seed</name>
      <defaultValue>219_c_long</defaultValue>
      <description>
      A seed value for the random number generator.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>ompThreadOffset</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, offset the key by the OpenMP thread number.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>mpiRankOffset</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, offset the key by the MPI process rank.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=randomNumberGeneratorPhilox(seed,ompThreadOffset,mpiRankOffset)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function philoxConstructorParameters

  function philoxConstructorInternal(seed_,ompThreadOffset,mpiRankOffset) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`randomNumberGeneratorPhilox` random number generator class.
    !!}
#ifdef USEMPI
    use    :: MPI_F08, only : MPI_Comm_Rank     , MPI_Comm_World
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num
    implicit none
    type   (randomNumberGeneratorPhilox)                          :: self
#ifdef USEMPI
    integer                                                       :: mpiRank        , iError
#endif
    integer(c_long                     ), intent(in   )           :: seed_
    logical                             , intent(in   ), optional :: ompThreadOffset, mpiRankOffset
    !![
    <constructorAssign variables="seed_"/>
    <optionalArgument name="ompThreadOffset" defaultsTo=".false."/>
    <optionalArgument name="mpiRankOffset"   defaultsTo=".false."/>
    !!]

    self%ompThreadOffset=ompThreadOffset_
    self%mpiRankOffset  =mpiRankOffset_
    ! The second word of the key holds the OpenMP thread number (in its lower 32 bits) and the MPI process rank (in its upper 32
    ! bits). Different keys give unrelated sequences, so offsetting the key, rather than the seed, guarantees that the sequences
    ! of different threads and processes do not overlap.
    self%key(2)=0_c_int64_t
    !$ if (ompThreadOffset_) self%key(2)=self%key(2)+int(OMP_Get_Thread_Num(),kind=c_int64_t)
    if (mpiRankOffset_) then
#ifdef USEMPI
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,iError)
       self%key(2)=self%key(2)+shiftl(int(mpiRank,kind=c_int64_t),32)
#endif
    end if
    call self%seedSet(seed_,offset=.false.)
    return
  end function philoxConstructorInternal

  logical function philoxMPIIndependent(self)
    !!{RST
    Return true if this random number generator produces independent sequences per MPI process.
    !!}
    implicit none
    class(randomNumberGeneratorPhilox), intent(inout) :: self

#ifdef USEMPI
    philoxMPIIndependent=self%mpiRankOffset
#else
    philoxMPIIndependent=.false.
#endif
    return
  end function philoxMPIIndependent

  logical function philoxOpenMPIndependent(self)
    !!{RST
    Return true if this random number generator produces independent sequences per OpenMP thread.
    !!}
    implicit none
    class(randomNumberGeneratorPhilox), intent(inout) :: self

    philoxOpenMPIndependent=self%ompThreadOffset
    return
  end function philoxOpenMPIndependent

  function philoxRangeMinimum(self)
    !!{RST
    Return the minimum of the range.
    !!}
    implicit none
    integer(c_long                     )                :: philoxRangeMinimum
    class  (randomNumberGeneratorPhilox), intent(inout) :: self
    !$GLC attributes unused :: self

    philoxRangeMinimum=0_c_long
    return
  end function philoxRangeMinimum

  function philoxRangeMaximum(self)
    !!{RST
    Return the maximum of the range.
    !!}
    implicit none
    integer(c_long                     )                :: philoxRangeMaximum
    class  (randomNumberGeneratorPhilox), intent(inout) :: self
    !$GLC attributes unused :: self

    philoxRangeMaximum=huge(0_c_long)
    return
  end function philoxRangeMaximum

  function philoxBits(self) result(bits)
    !!{RST
    Return the next 64 random bits from the generator, generating a new block when the current one is used up.
    !!}
    implicit none
    integer(c_int64_t                  )                :: bits
    class  (randomNumberGeneratorPhilox), intent(inout) :: self

    if (self%bufferPosition > 4) then
       call philox4x64Next(self%counter,self%key,self%buffer)
       self%bufferPosition=1
    end if
    bits               =self%buffer(self%bufferPosition)
    self%bufferPosition=self%bufferPosition+1
    return
  end function philoxBits

  function philoxSample(self,n)
    !!{RST
    Sample a random integer. If ``n`` is present the integer is uniformly distributed on :math:`[0,n-1]`. Otherwise it is
    uniformly distributed on :math:`[1,2^{63}-1]`, matching the strictly positive integers returned by other generators.
    !!}
    use :: Error, only : Error_Report
    implicit none
    integer(c_long                     )                          :: philoxSample
    class  (randomNumberGeneratorPhilox), intent(inout)           :: self
    integer(c_long                     ), intent(in   ), optional :: n
    integer(c_long                     )                          :: limit

    if (present(n)) then
       if (n <= 0_c_long) call Error_Report('n must be positive'//{introspection:location})
       ! Use the upper 63 bits, which are uniform on [0,2⁶³-1], and reject values at or above the largest multiple of n not
       ! exceeding 2⁶³, so that every residue modulo n is equally likely. Here, 2⁶³ mod n is computed as ((2⁶³-1) mod n + 1) mod n.
       limit=huge(0_c_long)-mod(mod(huge(0_c_long),n)+1_c_long,n)
       philoxSample=limit+1_c_long
       do while (philoxSample > limit)
          philoxSample=shiftr(philoxBits(self),1)
       end do
       philoxSample=mod(philoxSample,n)
    else
       philoxSample=0_c_long
       do while (philoxSample == 0_c_long)
          philoxSample=shiftr(philoxBits(self),1)
       end do
    end if
    return
  end function philoxSample

  double precision function philoxUniformSample(self)
    !!{RST
    Sample from a uniform distribution on the interval :math:`[0,1)`, using the upper 53 bits of a 64-bit random integer.
    !!}
    implicit none
    class(randomNumberGeneratorPhilox), intent(inout) :: self

    philoxUniformSample=dble(shiftr(philoxBits(self),11))*2.0d0**(-53)
    return
  end function philoxUniformSample

  double precision function philoxUniformSamplePositive(self)
    !!{RST
    Sample from a uniform distribution on the interval :math:`(0,1]`, using the upper 53 bits of a 64-bit random integer.
    !!}
    implicit none
    class(randomNumberGeneratorPhilox), intent(inout) :: self

    philoxUniformSamplePositive=dble(shiftr(philoxBits(self),11)+1_c_int64_t)*2.0d0**(-53)
    return
  end function philoxUniformSamplePositive

  double precision function philoxStandardNormalSample(self)
    !!{RST
    Sample from a standard normal distribution, using the :cite:t:`box_note_1958` transform. Each transform gives two
    independent deviates, the second of which is returned by the next call.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (randomNumberGeneratorPhilox), intent(inout) :: self
    double precision                                             :: radius, angle

    if (self%normalAvailable) then
       philoxStandardNormalSample=self%normalStored
       self%normalAvailable      =.false.
    else
       ! Use a uniform deviate on (0,1] for the radius, so that its logarithm is finite.
       radius                    =sqrt(-2.0d0*log(philoxUniformSamplePositive(self)))
       angle                     =2.0d0*Pi*philoxUniformSample(self)
       philoxStandardNormalSample=radius*cos(angle)
       self%normalStored         =radius*sin(angle)
       self%normalAvailable      =.true.
    end if
    return
  end function philoxStandardNormalSample

  integer function philoxPoissonSample(self,mean)
    !!{RST
    Sample from a Poisson distribution with the given mean. For means below 10, the product of uniform deviates is compared with
    :math:`\exp(-\lambda)`. For larger means the transformed rejection method with squeeze (PTRS) of
    :cite:t:`hormann_transformed_1993` is used.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (randomNumberGeneratorPhilox), intent(inout) :: self
    double precision                             , intent(in   ) :: mean
    ! The mean below which inversion by multiplication of uniform deviates is used (Hörmann 1993 recommends PTRS for λ ≥ 10).
    double precision                             , parameter     :: meanTransformedRejection=10.0d0
    double precision                                             :: product                       , limit       , &
         &                                                          u                             , v           , &
         &                                                          uShifted                      , a           , &
         &                                                          b                             , alphaInverse, &
         &                                                          vR                            , meanLog     , &
         &                                                          meanRoot
    integer                                                      :: k

    if (mean < 0.0d0) call Error_Report('mean must be non-negative'//{introspection:location})
    if (mean < meanTransformedRejection) then
       limit  =exp(-mean)
       product=philoxUniformSample(self)
       k      =0
       do while (product > limit)
          k      =k      +1
          product=product*philoxUniformSample(self)
       end do
       philoxPoissonSample=k
    else
       ! Constants of the PTRS algorithm (Hörmann 1993, section 3, step 0).
       meanRoot    =sqrt(mean)
       meanLog     =log (mean)
       b           =+0.931d0+2.53d0*meanRoot
       a           =-0.059d0+0.02483d0*b
       alphaInverse=+1.1239d0+1.1328d0/(b-3.4d0)
       vR          =+0.9277d0-3.6224d0/(b-2.0d0)
       do while (.true.)
          u       =philoxUniformSample(self)-0.5d0
          v       =philoxUniformSample(self)
          uShifted=0.5d0-abs(u)
          k       =floor((2.0d0*a/uShifted+b)*u+mean+0.43d0)
          ! Accept immediately within the squeeze region.
          if (uShifted >= 0.07d0 .and. v <= vR) exit
          ! Reject outside the support, or in the tails of the hat function.
          if (k < 0 .or. (uShifted < 0.013d0 .and. v > uShifted)) cycle
          ! Otherwise accept by comparison with the Poisson probability.
          if     (                                                &
               &   log(v)+log(alphaInverse)-log(a/uShifted**2+b)  &
               &  <=                                              &
               &   -mean+dble(k)*meanLog-log_gamma(dble(k)+1.0d0) &
               & ) exit
       end do
       philoxPoissonSample=k
    end if
    return
  end function philoxPoissonSample

  function philoxSeed(self) result(seed)
    !!{RST
    Return the seed for this random number generator.
    !!}
    implicit none
    integer(c_long                     )                :: seed
    class  (randomNumberGeneratorPhilox), intent(inout) :: self

    seed=self%seed__
    return
  end function philoxSeed

  subroutine philoxSeedSet(self,seed,offset)
    !!{RST
    Reset the seed for this random number generator. The seed becomes the first word of the key, and the counter is reset to
    zero. If ``offset`` is true the seed is added to the seed given at construction.
    !!}
    implicit none
    class  (randomNumberGeneratorPhilox), intent(inout) :: self
    integer(c_long                     ), intent(in   ) :: seed
    logical                             , intent(in   ) :: offset
    integer(c_long                     )                :: seed_

    seed_=seed
    if (offset) seed_=seed_+self%seed_
    self%seed__=seed_
    self%key(1)=int(seed_,kind=c_int64_t)
    call self%counterSet([0_c_int64_t,0_c_int64_t,0_c_int64_t])
    return
  end subroutine philoxSeedSet

  subroutine philoxCounterSet(self,counter)
    !!{RST
    Position the generator at the start of the stream labeled by ``counter``, which gives the upper three words of the 256-bit
    counter. The lowest word, which is incremented as random numbers are drawn, is reset to zero, and any buffered random bits
    and stored Gaussian deviate are discarded. The sequence that follows depends only on the key and on ``counter``.
    !!}
    implicit none
    class  (randomNumberGeneratorPhilox), intent(inout)               :: self
    integer(c_int64_t                  ), intent(in   ), dimension(3) :: counter

    self%counter        (1  )=0_c_int64_t
    self%counter        (2:4)=counter
    self%bufferPosition      =5
    self%normalAvailable      =.false.
    self%normalStored         =0.0d0
    self%buffer               =0_c_int64_t
    return
  end subroutine philoxCounterSet
