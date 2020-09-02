!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements a random number generator class which utilizes the \gls{gsl} random number generators.

  ! Add dependency on GSL library.
  !; gsl

  use, intrinsic :: ISO_C_Binding, only : c_long, c_ptr
  
  !# <randomNumberGenerator name="randomNumberGeneratorGSL">
  !#  <description>A random number generator class which utilizes the \gls{gsl} random number generators.</description>
  !# </randomNumberGenerator>
  type, extends(randomNumberGeneratorClass) :: randomNumberGeneratorGSL
     !% A random number generator class which utilizes the \gls{gsl} random number generators.
     private
     integer(c_long) :: seed
     logical         :: ompThreadOffset         , mpiRankOffset
     type   (c_ptr ) :: gslRandomNumberGenerator
   contains
     final     ::                         gslDestructor
     procedure :: mpiIndependent       => gslMPIIndependent
     procedure :: uniformSample        => gslUniformSample
     procedure :: poissonSample        => gslPoissonSample
     procedure :: standardNormalSample => gslStandardNormalSample
     procedure :: seedSet              => gslSeedSet
     procedure :: deepCopy             => gslDeepCopy
     procedure :: stateStore           => gslStateStore
     procedure :: stateRestore         => gslStateRestore
  end type randomNumberGeneratorGSL

  interface randomNumberGeneratorGSL
     !% Constructors for the {\normalfont \ttfamily gsl} merger tree evolve profiler class.
     module procedure gslConstructorParameters
     module procedure gslConstructorInternal
  end interface randomNumberGeneratorGSL

  ! GSL RNG initialization status.
  logical                                         :: gslRNGInitialized=.false.
  
  ! RNG types.
  type   (c_ptr), bind(C, name="gsl_rng_default") :: gsl_rng_default

  interface
     function gsl_rng_alloc(T) bind(c,name='gsl_rng_alloc')
       !% Template for the GSL random number generator allocate function.
       import c_ptr
       type(c_ptr)        :: gsl_rng_alloc
       type(c_ptr), value :: T
     end function gsl_rng_alloc
     
     subroutine gsl_rng_set(r,s) bind(c,name='gsl_rng_set')
       !% Template for the GSL random number generator set function.
       import c_ptr, c_long
       type   (c_ptr ), value :: r
       integer(c_long), value :: s
     end subroutine gsl_rng_set

     subroutine gsl_rng_free(r) bind(c,name='gsl_rng_free')
       !% Template for the GSL random number generator free function.
       import c_ptr
       type(c_ptr), value :: r
     end subroutine gsl_rng_free

     function gsl_rng_uniform(r) bind(c,name='gsl_rng_uniform')
       !% Template for the GSL random number generator uniform sampling function.
       import c_ptr, c_double
       real(c_double)        :: gsl_rng_uniform
       type(c_ptr   ), value :: r
     end function gsl_rng_uniform

     function gsl_ran_poisson(r,mu) bind(c,name='gsl_ran_poisson')
       !% Template for the GSL random number generator Poisson sampling function.
       import c_ptr, c_int, c_double
       integer(c_int   )        :: gsl_ran_poisson
       type   (c_ptr   ), value :: r
       real   (c_double), value :: mu
     end function gsl_ran_poisson

     function gsl_ran_gaussian(r,sigma) bind(c,name='gsl_ran_gaussian')
       !% Template for the GSL random number generator Gaussian sampling function.
       import c_ptr, c_double
       real(c_double)        :: gsl_ran_gaussian
       type(c_ptr   ), value :: r
       real(c_double), value :: sigma
     end function gsl_ran_gaussian

     function gsl_rng_clone(r) bind(c,name='gsl_rng_clone')
       !% Template for the GSL random number generator state clone function.
       import c_ptr
       type(c_ptr)        :: gsl_rng_clone
       type(c_ptr), value :: r
     end function gsl_rng_clone

     function gsl_rng_fwrite(stream,r) bind(c,name='gsl_rng_fwrite')
       !% Template for the GSL random number generator state write function.
       import c_ptr, c_int
       integer(c_int)        :: gsl_rng_fwrite
       type   (c_ptr), value :: stream        , r
     end function gsl_rng_fwrite

     function gsl_rng_fread(stream,r) bind(c,name='gsl_rng_fread')
       !% Template for the GSL random number generator state read function.
       import c_ptr, c_int
       integer(c_int)        :: gsl_rng_fread
       type   (c_ptr), value :: stream       , r
     end function gsl_rng_fread

     subroutine gsl_rng_env_setup() bind(c,name='gsl_rng_env_setup')
       !% Template for function used to set up default GSL random number generator.
     end subroutine gsl_rng_env_setup
  end interface
  
contains
  
  function gslConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily gsl} random number generator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (randomNumberGeneratorGSL)                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    integer(c_long                  )                :: seed
    logical                                          :: ompThreadOffset, mpiRankOffset

    !# <inputParameter>
    !#   <name>seed</name>
    !#   <defaultValue>219_c_long</defaultValue>
    !#   <description>A seed value for the random number generator.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>ompThreadOffset</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, offset the seed by the OpenMP thread number.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>mpiRankOffset</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, offset the seed by the MPI process rank.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=randomNumberGeneratorGSL(seed,ompThreadOffset,mpiRankOffset)
    !# <inputParametersValidate source="parameters"/>
    return
  end function gslConstructorParameters

  function gslConstructorInternal(seed,ompThreadOffset,mpiRankOffset) result(self)
    !% Internal constructor for the {\normalfont \ttfamily gsl} random number generator class.
#ifdef USEMPI
    use :: MPI, only : MPI_Comm_Rank, MPI_Comm_World
#endif   
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num
    implicit none
    type   (randomNumberGeneratorGSL)                          :: self
#ifdef USEMPI
    integer                                                    :: mpiRank        , iError
#endif
    integer(c_long                  ), intent(in   )           :: seed
    logical                          , intent(in   ), optional :: ompThreadOffset, mpiRankOffset
    integer(c_long                  )                          :: seed_
    !# <constructorAssign variables="seed"/>
    !# <optionalArgument name="ompThreadOffset" defaultsTo=".false."/>
    !# <optionalArgument name="mpiRankOffset"   defaultsTo=".false."/>

    self%ompThreadOffset         =ompThreadOffset_
    self%mpiRankOffset           =mpiRankOffset_
    if (.not.gslRNGInitialized) then
       !$omp critical(GSL_RNG_Initialize)
       if (.not.gslRNGInitialized) then
          call gsl_rng_env_setup()
          gslRNGInitialized=.true.
       end if
       !$omp end critical(GSL_RNG_Initialize)
    end if
    self%gslRandomNumberGenerator=GSL_RNG_Alloc(GSL_Rng_Default)
    seed_                        =seed
    !$ if (ompThreadOffset_) seed_=seed_+OMP_Get_Thread_Num()
    if (mpiRankOffset_) then
#ifdef USEMPI
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,iError)
       seed_=seed_+mpiRank
#endif
    end if
    call GSL_RNG_Set(self%gslRandomNumberGenerator,seed_)
    return
  end function gslConstructorInternal
  
  subroutine gslDestructor(self)
    !% Destructor for the {\normalfont \ttfamily gsl} random number generator class.
    implicit none
    type(randomNumberGeneratorGSL), intent(inout) :: self

    call GSL_RNG_Free(self%gslRandomNumberGenerator)
    return
  end subroutine gslDestructor

  logical function gslMPIIndependent(self)
    !% Return true if this random number generator produces independent sequences per MPI process.
    implicit none    
    class(randomNumberGeneratorGSL), intent(inout) :: self
    
    gslMPIIndependent=self%mpiRankOffset
    return
  end function gslMPIIndependent

  double precision function gslUniformSample(self)
    !% Sample from a uniform distribution on the interval [0,1).
    implicit none
    class(randomNumberGeneratorGSL), intent(inout) :: self
    
    gslUniformSample=GSL_RNG_Uniform(self%gslRandomNumberGenerator)
    return
  end function gslUniformSample

  integer function gslPoissonSample(self,mean)
    !% Sample from a Poisson distribution with the given mean.
    implicit none
    class           (randomNumberGeneratorGSL), intent(inout) :: self
    double precision                          , intent(in   ) :: mean

    gslPoissonSample=GSL_Ran_Poisson(self%gslRandomNumberGenerator,mean)
    return
  end function gslPoissonSample

  double precision function gslStandardNormalSample(self)
    !% Sample from a standard normal distribution.
    implicit none
    class(randomNumberGeneratorGSL), intent(inout) :: self

    gslStandardNormalSample=GSL_Ran_Gaussian(self%gslRandomNumberGenerator,sigma=1.0d0)
    return
  end function gslStandardNormalSample

  subroutine gslSeedSet(self,seed,offset)
    !% Reset the seed for this random number generator.
    implicit none
    class  (randomNumberGeneratorGSL), intent(inout) :: self
    integer(c_long                  ), intent(in   ) :: seed
    logical                          , intent(in   ) :: offset
    integer(c_long                  )                :: seed_

    seed_=seed
    if (offset) seed_=seed_+self%seed
    call GSL_RNG_Set(self%gslRandomNumberGenerator,seed_)
    return
  end subroutine gslSeedSet

  subroutine gslDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily GSL} random number generator class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(randomNumberGeneratorGSL  ), intent(inout) :: self
    class(randomNumberGeneratorClass), intent(inout) :: destination

    call self%randomNumberGeneratorClass%deepCopy(destination)
    select type (destination)
    type is (randomNumberGeneratorGSL)
       destination%seed                    =              self%seed
       destination%ompThreadOffset         =              self%ompThreadOffset
       destination%mpiRankOffset           =              self%mpiRankOffset
       destination%gslRandomNumberGenerator=GSL_Rng_Clone(self%gslRandomNumberGenerator)
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine gslDeepCopy
  
  subroutine gslStateStore(self,stateFile,gslStateFile,stateOperationID)
    !% Store the state of this object to file.
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: Interface_GSL     , only : gsl_success
    use            :: String_Handling   , only : operator(//)
    use            :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, Galacticus_Verbosity_Level, &
         &                                       verbosityWorking
    use            :: Galacticus_Error  , only : Galacticus_Error_Report
    use            :: ISO_Varying_String, only : var_str
    implicit none
    class    (randomNumberGeneratorGSL), intent(inout) :: self
    integer                            , intent(in   ) :: stateFile
    type     (c_ptr                   ), intent(in   ) :: gslStateFile
    integer  (c_size_t                ), intent(in   ) :: stateOperationID
    character(len=16                  )                :: label
    integer                                            :: status

    call Galacticus_Display_Indent(var_str('storing state for "randomNumberGenerator" [position: ')//FTell(stateFile)//']',verbosity=verbosityWorking)
    if (self%stateOperationID == stateOperationID) then
       call Galacticus_Display_Unindent('skipping - already stored',verbosity=verbosityWorking)
       return
    end if
    self%stateOperationID=stateOperationID
    call Galacticus_Display_Message('object type "randomNumberGeneratorGSL"',verbosity=verbosityWorking)
    if (Galacticus_Verbosity_Level() >= verbosityWorking) then
       write (label,'(i16)') sizeof(self%seed)
       call Galacticus_Display_Message('storing "seed" with size '//trim(adjustl(label))//' bytes')
    end if
    if (Galacticus_Verbosity_Level() >= verbosityWorking) then
       write (label,'(i16)') sizeof(self%ompthreadoffset)
       call Galacticus_Display_Message('storing "ompthreadoffset" with size '//trim(adjustl(label))//' bytes')
    end if
    if (Galacticus_Verbosity_Level() >= verbosityWorking) then
       write (label,'(i16)') sizeof(self%mpirankoffset)
       call Galacticus_Display_Message('storing "mpirankoffset" with size '//trim(adjustl(label))//' bytes')
    end if
    write (stateFile) self%seed,self%ompThreadOffset,self%mpiRankOffset
    status=GSL_Rng_FWrite(gslStateFile,self%gslRandomNumberGenerator)
    if (status /= GSL_Success) call Galacticus_Error_Report('failed to store GSL random number generator state'//{introspection:location})
    call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)
    return
  end subroutine gslStateStore

  subroutine gslStateRestore(self,stateFile,gslStateFile,stateOperationID)
    !% Restore the state of this object from file.
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: Interface_GSL     , only : gsl_success
    use            :: String_Handling   , only : operator(//)
    use            :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, Galacticus_Verbosity_Level, &
         &                                       verbosityWorking
    use            :: Galacticus_Error  , only : Galacticus_Error_Report
    use            :: ISO_Varying_String, only : var_str
    implicit none
    class  (randomNumberGeneratorGSL), intent(inout) :: self
    integer                          , intent(in   ) :: stateFile
    type   (c_ptr                   ), intent(in   ) :: gslStateFile
    integer(c_size_t                ), intent(in   ) :: stateOperationID
    integer                                          :: status

    call Galacticus_Display_Indent(var_str('restoring state for "randomNumberGenerator" [position: ')//FTell(stateFile)//']',verbosity=verbosityWorking)
    if (self%stateOperationID == stateOperationID) then
       call Galacticus_Display_Unindent('skipping - already restored',verbosity=verbosityWorking)
       return
    end if
    self%stateOperationID=stateOperationID
    call Galacticus_Display_Message('object type "randomNumberGeneratorGSL"',verbosity=verbosityWorking)
    call Galacticus_Display_Message('restoring "seed"'                      ,verbosity=verbosityWorking)
    call Galacticus_Display_Message('restoring "ompthreadoffset"'           ,verbosity=verbosityWorking)
    call Galacticus_Display_Message('restoring "mpirankoffset"'             ,verbosity=verbosityWorking)
    read (stateFile) self%seed,self%ompThreadOffset,self%mpiRankOffset
    status=GSL_Rng_FRead(gslStateFile,self%gslRandomNumberGenerator)
    if (status /= GSL_Success) call Galacticus_Error_Report('failed to store GSL random number generator state'//{introspection:location})
    call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)
    return
  end subroutine gslStateRestore
