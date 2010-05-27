!% Contains a module which implements pseudo-random numbers.

module Pseudo_Random
  !% Implements pseudo-random numbers.
  use FGSL
  use, intrinsic :: ISO_C_Binding
  private
  public :: Pseudo_Random_Get, Pseudo_Random_Free, Pseudo_Random_Store, Pseudo_Random_Retrieve
  
  logical           :: Seed_Is_Set=.false.
  integer(c_size_t) :: randomSeedC

contains

  double precision function Pseudo_Random_Get(pseudoSequenceObject,reset)
    !% Returns a scalar giving a pseudo-random number.
    use Input_Parameters
    implicit none
    type(fgsl_rng), intent(inout)           :: pseudoSequenceObject
    logical,        intent(inout), optional :: reset
    integer                                 :: randomSeed
    logical                                 :: resetActual

    ! Determine if we need to reset.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
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
       !@ </inputParameter>
       call Get_Input_Parameter('randomSeed',randomSeed,defaultValue=219)
       randomSeedC=randomSeed
       Seed_Is_Set=.true.
    end if
    !$omp end critical (Pseudo_Random_Get)
    
    if (resetActual.or..not.FGSL_Well_Defined(pseudoSequenceObject)) then
       pseudoSequenceObject=FGSL_RNG_Alloc(FGSL_RNG_Default)
       randomSeed=randomSeed+1
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
    type(fgsl_rng),  intent(in) :: pseudoSequenceObject
    type(fgsl_file), intent(in) :: fgslFile
    integer                     :: iError

    iError=FGSL_Rng_FWrite(fgslFile,pseudoSequenceObject)
    return
  end subroutine Pseudo_Random_Store

  subroutine Pseudo_Random_Retrieve(pseudoSequenceObject,fgslFile)
    !% Stores a pseudo-random sequence object to file.
    implicit none
    type(fgsl_rng),  intent(inout) :: pseudoSequenceObject
    type(fgsl_file), intent(in)    :: fgslFile
    integer                        :: iError

    if (.not.FGSL_Well_Defined(pseudoSequenceObject)) pseudoSequenceObject=FGSL_RNG_Alloc(FGSL_RNG_Default)
    iError=FGSL_Rng_FRead(fgslFile,pseudoSequenceObject)
    return
  end subroutine Pseudo_Random_Retrieve

end module Pseudo_Random
