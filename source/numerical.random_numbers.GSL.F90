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

  use :: FGSL, only : FGSL_RNG
  
  !# <randomNumberGenerator name="randomNumberGeneratorGSL">
  !#  <description>A random number generator class which utilizes the \gls{gsl} random number generators.</description>
  !# </randomNumberGenerator>
  type, extends(randomNumberGeneratorClass) :: randomNumberGeneratorGSL
     !% A random number generator class which utilizes the \gls{gsl} random number generators.
     private
     integer(FGSL_Long) :: seed
     logical            :: ompThreadOffset         , mpiRankOffset
     type   (FGSL_RNG ) :: gslRandomNumberGenerator
   contains
     final     ::                         gslDestructor
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

contains
  
  function gslConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily gsl} random number generator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (randomNumberGeneratorGSL)                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    integer(FGSL_Long               )                :: seed
    logical                                          :: ompThreadOffset, mpiRankOffset

    !# <inputParameter>
    !#   <name>seed</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>219_FGSL_Long</defaultValue>
    !#   <description>A seed value for the random number generator.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>ompThreadoffset</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, offset the seed by the OpenMP thread number.</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>mpiRankOffset</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, offset the seed by the MPI process rank.</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    self=randomNumberGeneratorGSL(seed,ompThreadOffset,mpiRankOffset)
    !# <inputParametersValidate source="parameters"/>
    return
  end function gslConstructorParameters

  function gslConstructorInternal(seed,ompThreadOffset,mpiRankOffset) result(self)
    !% Internal constructor for the {\normalfont \ttfamily gsl} random number generator class.
    use    :: FGSL   , only : FGSL_RNG_Alloc    , FGSL_RNG_Default, FGSL_RNG_Set
#ifdef USEMPI
    use    :: MPI    , only : MPI_Comm_Rank     , MPI_Comm_World
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num
    implicit none
    type   (randomNumberGeneratorGSL)                          :: self
#ifdef USEMPI
    integer                                                    :: mpiRank        , iError
#endif
    integer(FGSL_Long               ), intent(in   )           :: seed
    logical                          , intent(in   ), optional :: ompThreadOffset, mpiRankOffset
    integer(FGSL_Long               )                          :: seed_
    !# <constructorAssign variables="seed"/>
    !# <optionalArgument name="ompThreadOffset" defaultsTo=".false."/>
    !# <optionalArgument name="mpiRankOffset"   defaultsTo=".false."/>

    self%ompThreadOffset         =ompThreadOffset_
    self%mpiRankOffset           =mpiRankOffset_
    self%gslRandomNumberGenerator=FGSL_RNG_Alloc(FGSL_RNG_Default)
    seed_                        =seed
    !$ if (ompThreadOffset_) seed_=seed_+OMP_Get_Thread_Num()
    if (mpiRankOffset_) then
#ifdef USEMPI
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,iError)
       seed_=seed_+mpiRank
#endif
    end if
    call FGSL_RNG_Set(self%gslRandomNumberGenerator,seed_)
    return
  end function gslConstructorInternal
  
  subroutine gslDestructor(self)
    !% Destructor for the {\normalfont \ttfamily gsl} random number generator class.
    use :: FGSL, only : FGSL_RNG_Free
    implicit none
    type(randomNumberGeneratorGSL), intent(inout) :: self

    call FGSL_RNG_Free(self%gslRandomNumberGenerator)
    return
  end subroutine gslDestructor

  double precision function gslUniformSample(self)
    !% Sample from a uniform distribution on the interval [0,1).
    use :: FGSL, only : FGSL_RNG_Uniform
    implicit none
    class(randomNumberGeneratorGSL), intent(inout) :: self
    
    gslUniformSample=FGSL_RNG_Uniform(self%gslRandomNumberGenerator)
    return
  end function gslUniformSample

  integer function gslPoissonSample(self,mean)
    !% Sample from a Poisson distribution with the given mean.
    use :: FGSL, only : FGSL_Ran_Poisson
    implicit none
    class           (randomNumberGeneratorGSL), intent(inout) :: self
    double precision                          , intent(in   ) :: mean

    gslPoissonSample=FGSL_Ran_Poisson(self%gslRandomNumberGenerator,mean)
    return
  end function gslPoissonSample

  double precision function gslStandardNormalSample(self)
    !% Sample from a standard normal distribution.
    use :: FGSL, only : FGSL_Ran_Gaussian
    implicit none
    class(randomNumberGeneratorGSL), intent(inout) :: self

    gslStandardNormalSample=FGSL_Ran_Gaussian(self%gslRandomNumberGenerator,sigma=1.0d0)
    return
  end function gslStandardNormalSample

  subroutine gslSeedSet(self,seed,offset)
    !% Reset the seed for this random number generator.
    use :: FGSL, only : FGSL_RNG_Set
    implicit none
    class  (randomNumberGeneratorGSL), intent(inout) :: self
    integer(FGSL_Long               ), intent(in   ) :: seed
    logical                          , intent(in   ) :: offset
    integer(FGSL_Long               )                :: seed_

    seed_=seed
    if (offset) seed_=seed_+self%seed
    call FGSL_RNG_Set(self%gslRandomNumberGenerator,seed_)
    return
  end subroutine gslSeedSet

  subroutine gslDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily GSL} random number generator class.
    use :: FGSL            , only : FGSL_Rng_Clone
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(randomNumberGeneratorGSL  ), intent(inout) :: self
    class(randomNumberGeneratorClass), intent(inout) :: destination

    call self%randomNumberGeneratorClass%deepCopy(destination)
    select type (destination)
    type is (randomNumberGeneratorGSL)
       destination%seed                    =               self%seed
       destination%ompThreadOffset         =               self%ompThreadOffset
       destination%mpiRankOffset           =               self%mpiRankOffset
       destination%gslRandomNumberGenerator=FGSL_Rng_Clone(self%gslRandomNumberGenerator)
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine gslDeepCopy
  
  subroutine gslStateStore(self,stateFile,fgslStateFile,stateOperationID)
    !% Store the state of this object to file.
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: String_Handling   , only : operator(//)
    use            :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, Galacticus_Verbosity_Level, &
         &                                       verbosityWorking
    use            :: Galacticus_Error  , only : Galacticus_Error_Report
    use            :: FGSL              , only : FGSL_File                , FGSL_Rng_FWrite            , FGSL_Success
    use            :: ISO_Varying_String, only : var_str
    implicit none
    class    (randomNumberGeneratorGSL), intent(inout) :: self
    integer                            , intent(in   ) :: stateFile
    type     (FGSL_file               ), intent(in   ) :: fgslStateFile
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
    status=FGSL_Rng_FWrite(fgslStateFile,self%gslRandomNumberGenerator)
    if (status /= FGSL_Success) call Galacticus_Error_Report('failed to store GSL random number generator state'//{introspection:location})
    call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)
    return
  end subroutine gslStateStore
  
  subroutine gslStateRestore(self,stateFile,fgslStateFile,stateOperationID)
    !% Restore the state of this object from file.
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: String_Handling   , only : operator(//)
    use            :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, Galacticus_Verbosity_Level, &
         &                                       verbosityWorking
    use            :: Galacticus_Error  , only : Galacticus_Error_Report
    use            :: FGSL              , only : FGSL_File                , FGSL_Rng_FRead             , FGSL_Success
    use            :: ISO_Varying_String, only : var_str
    implicit none
    class  (randomNumberGeneratorGSL), intent(inout) :: self
    integer                          , intent(in   ) :: stateFile
    type   (fgsl_file               ), intent(in   ) :: fgslStateFile
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
    status=FGSL_Rng_FRead(fgslStateFile,self%gslRandomNumberGenerator)
    if (status /= FGSL_Success) call Galacticus_Error_Report('failed to store GSL random number generator state'//{introspection:location})
    call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)
    return
  end subroutine gslStateRestore
