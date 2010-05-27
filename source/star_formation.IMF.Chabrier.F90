!% Contains a module which implements the Chabrier stellar initial mass function \citep{}.

module Star_Formation_IMF_Chabrier
  !% Implements the Chabrier stellar initial mass function.
  private
  public :: Star_Formation_IMF_Register_Chabrier, Star_Formation_IMF_Register_Name_Chabrier,&
       & Star_Formation_IMF_Recycled_Instantaneous_Chabrier, Star_Formation_IMF_Yield_Instantaneous_Chabrier,&
       & Star_Formation_IMF_Tabulate_Chabrier, Star_Formation_IMF_Minimum_Mass_Chabrier, Star_Formation_IMF_Maximum_Mass_Chabrier&
       &, Star_Formation_IMF_Phi_Chabrier

  ! Index assigned to this IMF.
  integer :: imfIndex=-1

  ! Flag indicating if the module has been initialized.
  logical :: imfChabrierInitialized=.false.

  ! Parameters of the IMF.
  double precision :: imfChabrierRecycledInstantaneous, imfChabrierYieldInstantaneous

  ! Fixed parameters of the IMF.
  double precision, parameter :: chabrierMassLower=0.1d0, chabrierMassUpper=125.0d0, chabrierMassTransition=1.0d0
  double precision, parameter :: chabrierExponent=-2.3d0, chabrierSigma=0.69d0, chabrierMass=0.08d0
  double precision, parameter :: normalizationLogNormal=0.8326617874d0, normalizationExponential=0.2353362401d0

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_Chabrier</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_Chabrier(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex=imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_Chabrier

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_Chabrier</unitName>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_Chabrier(imfNames)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(inout) :: imfNames(:)

    imfNames(imfIndex)="Chabrier"
    return
  end subroutine Star_Formation_IMF_Register_Name_Chabrier

  subroutine Star_Formation_IMF_Initialize_Chabrier
    !% Initialize the Chabrier IMF module.
    use Input_Parameters
    use Star_Formation_IMF_PPL
    implicit none

    !$omp critical (IMF_Chabrier_Initialize)
    if (.not.imfChabrierInitialized) then
       !@ <inputParameter>
       !@   <name>imfChabrierRecycledInstantaneous</name>
       !@   <defaultValue>0.46 (internally computed)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The recycled fraction for the Chabrier \IMF\ in the instantaneous recycling approximation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('imfChabrierRecycledInstantaneous',imfChabrierRecycledInstantaneous,defaultValue=0.46d0)
       !@ <inputParameter>
       !@   <name>imfChabrierYieldInstantaneous</name>
       !@   <defaultValue>0.035 (internally computed)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The yield for the Chabrier \IMF\ in the instantaneous recycling approximation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('imfChabrierYieldInstantaneous'   ,imfChabrierYieldInstantaneous   ,defaultValue=0.035d0)

       imfChabrierInitialized=.true.
    end if
    !$omp end critical (IMF_Chabrier_Initialize)
    return
  end subroutine Star_Formation_IMF_Initialize_Chabrier

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_Chabrier</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_Chabrier(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       minimumMass=chabrierMassLower
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_Chabrier

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_Chabrier</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_Chabrier(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       maximumMass=chabrierMassUpper
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_Chabrier

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_Chabrier</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_Chabrier(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(in)    :: initialMass
    double precision, intent(out)   :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       imfPhi=Chabrier_Phi(initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_Chabrier

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_Chabrier</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_Chabrier(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       recycledFraction=imfChabrierRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_Chabrier

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_Chabrier</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_Chabrier(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       yield=imfChabrierYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_Chabrier

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_Chabrier</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_Chabrier(imfSelected,imfMatched,imfMass,imfPhi)
    !% Register the name of this IMF.
    use Memory_Management
    use Numerical_Ranges
    use Star_Formation_IMF_PPL
    implicit none
    integer,          intent(in)                               :: imfSelected
    logical,          intent(inout)                            :: imfMatched
    double precision, intent(inout), allocatable, dimension(:) :: imfMass,imfPhi
    integer,          parameter                                :: nPoints=100

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       call Alloc_Array(imfMass,nPoints,'imfMass')
       call Alloc_Array(imfPhi ,nPoints,'imfPhi' )
       imfMass=Make_Range(chabrierMassLower,chabrierMassUpper,nPoints,rangeType=rangeTypeLogarithmic)
       imfPhi =Chabrier_Phi(imfMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Tabulate_Chabrier

  elemental double precision function Chabrier_Phi(mass)
    !% Evaluates the Chabrier initial mass function.
    implicit none
    double precision, intent(in) :: mass

    if (mass >= chabrierMassLower .and. mass < chabrierMassTransition) then
       Chabrier_Phi=normalizationLogNormal*dexp(-0.5d0*(dlog10(mass/chabrierMass)/chabrierSigma)**2)/mass
    else if (mass >= chabrierMassTransition .and. mass < chabrierMassUpper) then
       Chabrier_Phi=normalizationExponential*(mass**chabrierExponent)
    else
       Chabrier_Phi=0.0d0
    end if
    return
  end function Chabrier_Phi

end module Star_Formation_IMF_Chabrier
