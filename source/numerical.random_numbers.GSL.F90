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
  Implements a random number generator class which utilizes the \gls{gsl} random number generators.
  !!}

  ! Add dependency on GSL library.
  !; gsl

  use, intrinsic :: ISO_C_Binding, only : c_long, c_ptr, c_null_ptr
  
  !![
  <randomNumberGenerator name="randomNumberGeneratorGSL">
   <description>A random number generator class which utilizes the \gls{gsl} random number generators.</description>
  </randomNumberGenerator>
  !!]
  type, extends(randomNumberGeneratorClass) :: randomNumberGeneratorGSL
     !!{
     A random number generator class which utilizes the \gls{gsl} random number generators.
     !!}
     private
     integer(c_long) :: seed_                              , seed__
     logical         :: ompThreadOffset                    , mpiRankOffset
     type   (c_ptr ) :: gslRandomNumberGenerator=c_null_ptr
   contains
     final     ::                         gslDestructor
     procedure :: openMPIndependent    => gslOpenMPIndependent
     procedure :: mpiIndependent       => gslMPIIndependent
     procedure :: rangeMinimum         => gslRangeMinimum
     procedure :: rangeMaximum         => gslRangeMaximum
     procedure :: sample               => gslSample
     procedure :: uniformSample        => gslUniformSample
     procedure :: poissonSample        => gslPoissonSample
     procedure :: standardNormalSample => gslStandardNormalSample
     procedure :: seed                 => gslSeed
     procedure :: seedSet              => gslSeedSet
     procedure :: deepCopy             => gslDeepCopy
     procedure :: stateStore           => gslStateStore
     procedure :: stateRestore         => gslStateRestore
  end type randomNumberGeneratorGSL

  interface randomNumberGeneratorGSL
     !!{
     Constructors for the \refClass{randomNumberGeneratorGSL} merger tree evolve profiler class.
     !!}
     module procedure gslConstructorParameters
     module procedure gslConstructorInternal
  end interface randomNumberGeneratorGSL

  ! GSL RNG initialization status.
  logical                                         :: rngInitialized=.false.
  
  ! RNG types.
  type   (c_ptr), bind(C, name="gsl_rng_default") :: gsl_rng_default

  interface
     function gsl_rng_alloc(T) bind(c,name='gsl_rng_alloc')
       !!{
       Template for the GSL random number generator allocate function.
       !!}
       import c_ptr
       type(c_ptr)        :: gsl_rng_alloc
       type(c_ptr), value :: T
     end function gsl_rng_alloc
     
     subroutine gsl_rng_set(r,s) bind(c,name='gsl_rng_set')
       !!{
       Template for the GSL random number generator set function.
       !!}
       import c_ptr, c_long
       type   (c_ptr ), value :: r
       integer(c_long), value :: s
     end subroutine gsl_rng_set

     subroutine gsl_rng_free(r) bind(c,name='gsl_rng_free')
       !!{
       Template for the GSL random number generator free function.
       !!}
       import c_ptr
       type(c_ptr), value :: r
     end subroutine gsl_rng_free

     function gsl_rng_min(r) bind(c,name='gsl_rng_min')
       !!{
       Template for the GSL random number generator minimum function.
       !!}
       import c_ptr, c_long
       integer(c_long)        :: gsl_rng_min
       type   (c_ptr ), value :: r
     end function gsl_rng_min

     function gsl_rng_max(r) bind(c,name='gsl_rng_max')
       !!{
       Template for the GSL random number generator maximum function.
       !!}
       import c_ptr, c_long
       integer(c_long)        :: gsl_rng_max
       type   (c_ptr ), value :: r
     end function gsl_rng_max

     function gsl_rng_uniform_int(r,n) bind(c,name='gsl_rng_uniform_int')
       !!{
       Template for the GSL random number generator sampling function with a given upper bound.
       !!}
       import c_ptr, c_long
       integer(c_long)        :: gsl_rng_uniform_int
       type   (c_ptr ), value :: r
       integer(c_long), value :: n
     end function gsl_rng_uniform_int

     function gsl_rng_get(r) bind(c,name='gsl_rng_get')
       !!{
       Template for the GSL random number generator sampling function.
       !!}
       import c_ptr, c_long
       integer(c_long)        :: gsl_rng_get
       type   (c_ptr ), value :: r
     end function gsl_rng_get

     function gsl_rng_uniform(r) bind(c,name='gsl_rng_uniform')
       !!{
       Template for the GSL random number generator uniform sampling function.
       !!}
       import c_ptr, c_double
       real(c_double)        :: gsl_rng_uniform
       type(c_ptr   ), value :: r
     end function gsl_rng_uniform

     function gsl_ran_poisson(r,mu) bind(c,name='gsl_ran_poisson')
       !!{
       Template for the GSL random number generator Poisson sampling function.
       !!}
       import c_ptr, c_int, c_double
       integer(c_int   )        :: gsl_ran_poisson
       type   (c_ptr   ), value :: r
       real   (c_double), value :: mu
     end function gsl_ran_poisson

     function gsl_ran_gaussian(r,sigma) bind(c,name='gsl_ran_gaussian')
       !!{
       Template for the GSL random number generator Gaussian sampling function.
       !!}
       import c_ptr, c_double
       real(c_double)        :: gsl_ran_gaussian
       type(c_ptr   ), value :: r
       real(c_double), value :: sigma
     end function gsl_ran_gaussian

     function gsl_rng_clone(r) bind(c,name='gsl_rng_clone')
       !!{
       Template for the GSL random number generator state clone function.
       !!}
       import c_ptr
       type(c_ptr)        :: gsl_rng_clone
       type(c_ptr), value :: r
     end function gsl_rng_clone

     function gsl_rng_fwrite(stream,r) bind(c,name='gsl_rng_fwrite')
       !!{
       Template for the GSL random number generator state write function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_rng_fwrite
       type   (c_ptr), value :: stream        , r
     end function gsl_rng_fwrite

     function gsl_rng_fread(stream,r) bind(c,name='gsl_rng_fread')
       !!{
       Template for the GSL random number generator state read function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_rng_fread
       type   (c_ptr), value :: stream       , r
     end function gsl_rng_fread

     subroutine gsl_rng_env_setup() bind(c,name='gsl_rng_env_setup')
       !!{
       Template for function used to set up default GSL random number generator.
       !!}
     end subroutine gsl_rng_env_setup
  end interface
  
contains
  
  function gslConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{randomNumberGeneratorGSL} random number generator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (randomNumberGeneratorGSL)                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    integer(c_long                  )                :: seed
    logical                                          :: ompThreadOffset, mpiRankOffset

    !![
    <inputParameter>
      <name>seed</name>
      <defaultValue>219_c_long</defaultValue>
      <description>A seed value for the random number generator.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>ompThreadOffset</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, offset the seed by the OpenMP thread number.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mpiRankOffset</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, offset the seed by the MPI process rank.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=randomNumberGeneratorGSL(seed,ompThreadOffset,mpiRankOffset)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function gslConstructorParameters

  function gslConstructorInternal(seed_,ompThreadOffset,mpiRankOffset) result(self)
    !!{
    Internal constructor for the \refClass{randomNumberGeneratorGSL} random number generator class.
    !!}
#ifdef USEMPI
    use    :: MPI_F08, only : MPI_Comm_Rank     , MPI_Comm_World
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num
    implicit none
    type   (randomNumberGeneratorGSL)                          :: self
#ifdef USEMPI
    integer                                                    :: mpiRank        , iError
#endif
    integer(c_long                  ), intent(in   )           :: seed_
    logical                          , intent(in   ), optional :: ompThreadOffset, mpiRankOffset
    integer(c_long                  )                          :: seed__
    !![
    <constructorAssign variables="seed_"/>
    <optionalArgument name="ompThreadOffset" defaultsTo=".false."/>
    <optionalArgument name="mpiRankOffset"   defaultsTo=".false."/>
    !!]

    self%ompThreadOffset         =ompThreadOffset_
    self%mpiRankOffset           =mpiRankOffset_
    if (.not.rngInitialized) then
       !$omp critical(GSL_RNG_Initialize)
       if (.not.rngInitialized) then
          call gsl_rng_env_setup()
          rngInitialized=.true.
       end if
       !$omp end critical(GSL_RNG_Initialize)
    end if
    self%gslRandomNumberGenerator=GSL_RNG_Alloc(GSL_Rng_Default)
    seed__                         =seed_
    !$ if (ompThreadOffset_) seed__=seed__+OMP_Get_Thread_Num()
    if (mpiRankOffset_) then
#ifdef USEMPI
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,iError)
       seed__=seed__+mpiRank
#endif
    end if
    call GSL_RNG_Set(self%gslRandomNumberGenerator,seed__)
    self%seed_=seed__
    return
  end function gslConstructorInternal
  
  subroutine gslDestructor(self)
    !!{
    Destructor for the \refClass{randomNumberGeneratorGSL} random number generator class.
    !!}
    implicit none
    type(randomNumberGeneratorGSL), intent(inout) :: self

    if (c_associated(self%gslRandomNumberGenerator))        &
         & call GSL_RNG_Free(self%gslRandomNumberGenerator)
    return
  end subroutine gslDestructor

  logical function gslMPIIndependent(self)
    !!{
    Return true if this random number generator produces independent sequences per MPI process.
    !!}
    implicit none    
    class(randomNumberGeneratorGSL), intent(inout) :: self

#ifdef USEMPI
    gslMPIIndependent=self%mpiRankOffset
#else
    gslMPIIndependent=.false.
#endif
    return
  end function gslMPIIndependent

  logical function gslOpenMPIndependent(self)
    !!{
    Return true if this random number generator produces independent sequences per OpenMP thread.
    !!}
    implicit none    
    class(randomNumberGeneratorGSL), intent(inout) :: self
    
    gslOpenMPIndependent=self%ompThreadOffset
    return
  end function gslOpenMPIndependent

  function gslRangeMinimum(self)
    !!{
    Return the minimum of the range.
    !!}
    implicit none
    integer(c_long                  )                :: gslRangeMinimum
    class  (randomNumberGeneratorGSL), intent(inout) :: self

    gslRangeMinimum=GSL_RNG_Min(self%gslRandomNumberGenerator)
    return
  end function gslRangeMinimum

  function gslRangeMaximum(self)
    !!{
    Return the maximum of the range.
    !!}
    implicit none
    integer(c_long                  )                :: gslRangeMaximum
    class  (randomNumberGeneratorGSL), intent(inout) :: self

    gslRangeMaximum=GSL_RNG_Max(self%gslRandomNumberGenerator)
    return
  end function gslRangeMaximum

  function gslSample(self,n)
    !!{
    Sample a random integer.
    !!}
    implicit none
    integer(c_long                  )                          :: gslSample
    class  (randomNumberGeneratorGSL), intent(inout)           :: self
    integer(c_long                  ), intent(in   ), optional :: n

    if (present(n)) then
       gslSample=GSL_RNG_Uniform_Int(self%gslRandomNumberGenerator,n)
    else
       gslSample=-1_c_long
       do while (gslSample <= 0_c_long)
          gslSample=GSL_RNG_Get(self%gslRandomNumberGenerator)
       end do
    end if
    return
  end function gslSample

  double precision function gslUniformSample(self)
    !!{
    Sample from a uniform distribution on the interval [0,1).
    !!}
    implicit none
    class(randomNumberGeneratorGSL), intent(inout) :: self
    
    gslUniformSample=GSL_RNG_Uniform(self%gslRandomNumberGenerator)
    return
  end function gslUniformSample

  integer function gslPoissonSample(self,mean)
    !!{
    Sample from a Poisson distribution with the given mean.
    !!}
    implicit none
    class           (randomNumberGeneratorGSL), intent(inout) :: self
    double precision                          , intent(in   ) :: mean

    gslPoissonSample=GSL_Ran_Poisson(self%gslRandomNumberGenerator,mean)
    return
  end function gslPoissonSample

  double precision function gslStandardNormalSample(self)
    !!{
    Sample from a standard normal distribution.
    !!}
    implicit none
    class(randomNumberGeneratorGSL), intent(inout) :: self

    gslStandardNormalSample=GSL_Ran_Gaussian(self%gslRandomNumberGenerator,sigma=1.0d0)
    return
  end function gslStandardNormalSample

  function gslSeed(self) result(seed)
    !!{
    Return the seed for this random number generator.
    !!}
    implicit none
    integer(c_long                  )                :: seed
    class  (randomNumberGeneratorGSL), intent(inout) :: self

    seed=self%seed__
    return
  end function gslSeed

  subroutine gslSeedSet(self,seed,offset)
    !!{
    Reset the seed for this random number generator.
    !!}
    implicit none
    class  (randomNumberGeneratorGSL), intent(inout) :: self
    integer(c_long                  ), intent(in   ) :: seed
    logical                          , intent(in   ) :: offset
    integer(c_long                  )                :: seed_

    seed_=seed
    if (offset) seed_=seed_+self%seed_
    call GSL_RNG_Set(self%gslRandomNumberGenerator,seed_)
    self%seed__=seed_
    return
  end subroutine gslSeedSet

  subroutine gslDeepCopy(self,destination)
    !!{
    Perform a deep copy for the {\normalfont \ttfamily GSL} random number generator class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(randomNumberGeneratorGSL  ), intent(inout), target :: self
    class(randomNumberGeneratorClass), intent(inout)         :: destination

    call self%randomNumberGeneratorClass%deepCopy(destination)
    select type (destination)
    type is (randomNumberGeneratorGSL)
       destination%seed_                   =              self%seed_
       destination%seed__                  =              self%seed__
       destination%ompThreadOffset         =              self%ompThreadOffset
       destination%mpiRankOffset           =              self%mpiRankOffset
       destination%gslRandomNumberGenerator=GSL_Rng_Clone(self%gslRandomNumberGenerator)
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine gslDeepCopy
  
  subroutine gslStateStore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Store the state of this object to file.
    !!}
    use            :: Display           , only : displayIndent          , displayMessage, displayUnindent, displayVerbosity, &
          &                                      verbosityLevelWorking
    use            :: Error             , only : Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : var_str
    use            :: Interface_GSL     , only : gsl_success
    use            :: String_Handling   , only : operator(//)
    implicit none
    class    (randomNumberGeneratorGSL), intent(inout) :: self
    integer                            , intent(in   ) :: stateFile
    type     (c_ptr                   ), intent(in   ) :: gslStateFile
    integer  (c_size_t                ), intent(in   ) :: stateOperationID
    character(len=16                  )                :: label
    integer                                            :: status

    call displayIndent(var_str('storing state for "randomNumberGenerator" [position: ')//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
    if (self%stateOperationID == stateOperationID) then
       call displayUnindent('skipping - already stored',verbosity=verbosityLevelWorking)
       return
    end if
    self%stateOperationID=stateOperationID
    call displayMessage('object type "randomNumberGeneratorGSL"',verbosity=verbosityLevelWorking)
    if (displayVerbosity() >= verbosityLevelWorking) then
       write (label,'(i16)') sizeof(self%seed_)
       call displayMessage('storing "seed_" with size '//trim(adjustl(label))//' bytes')
    end if
    if (displayVerbosity() >= verbosityLevelWorking) then
       write (label,'(i16)') sizeof(self%seed__)
       call displayMessage('storing "seed__" with size '//trim(adjustl(label))//' bytes')
    end if
    if (displayVerbosity() >= verbosityLevelWorking) then
       write (label,'(i16)') sizeof(self%ompthreadoffset)
       call displayMessage('storing "ompthreadoffset" with size '//trim(adjustl(label))//' bytes')
    end if
    if (displayVerbosity() >= verbosityLevelWorking) then
       write (label,'(i16)') sizeof(self%mpirankoffset)
       call displayMessage('storing "mpirankoffset" with size '//trim(adjustl(label))//' bytes')
    end if
    write (stateFile) self%seed_,self%seed__,self%ompThreadOffset,self%mpiRankOffset
    status=GSL_Rng_FWrite(gslStateFile,self%gslRandomNumberGenerator)
    if (status /= GSL_Success) call Error_Report('failed to store GSL random number generator state'//{introspection:location})
    call displayUnindent('done',verbosity=verbosityLevelWorking)
    return
  end subroutine gslStateStore

  subroutine gslStateRestore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Restore the state of this object from file.
    !!}
    use            :: Display           , only : displayIndent          , displayMessage, displayUnindent, displayVerbosity, &
          &                                      verbosityLevelWorking
    use            :: Error             , only : Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : var_str
    use            :: Interface_GSL     , only : gsl_success
    use            :: String_Handling   , only : operator(//)
    implicit none
    class  (randomNumberGeneratorGSL), intent(inout) :: self
    integer                          , intent(in   ) :: stateFile
    type   (c_ptr                   ), intent(in   ) :: gslStateFile
    integer(c_size_t                ), intent(in   ) :: stateOperationID
    integer                                          :: status

    call displayIndent(var_str('restoring state for "randomNumberGenerator" [position: ')//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
    if (self%stateOperationID == stateOperationID) then
       call displayUnindent('skipping - already restored',verbosity=verbosityLevelWorking)
       return
    end if
    self%stateOperationID=stateOperationID
    call displayMessage('object type "randomNumberGeneratorGSL"',verbosity=verbosityLevelWorking)
    call displayMessage('restoring "seed_"'                     ,verbosity=verbosityLevelWorking)
    call displayMessage('restoring "seed__"'                    ,verbosity=verbosityLevelWorking)
    call displayMessage('restoring "ompthreadoffset"'           ,verbosity=verbosityLevelWorking)
    call displayMessage('restoring "mpirankoffset"'             ,verbosity=verbosityLevelWorking)
    read (stateFile) self%seed_,self%seed__,self%ompThreadOffset,self%mpiRankOffset
    status=GSL_Rng_FRead(gslStateFile,self%gslRandomNumberGenerator)
    if (status /= GSL_Success) call Error_Report('failed to store GSL random number generator state'//{introspection:location})
    call displayUnindent('done',verbosity=verbosityLevelWorking)
    return
  end subroutine gslStateRestore
