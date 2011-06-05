!% Contains a module which implements the virial overdensity for halos.

module Virial_Density_Contrast
  !% Implements the virial overdensity for halos.
  use ISO_Varying_String
  use FGSL
  private
  public :: Halo_Virial_Density_Contrast, Halo_Virial_Density_Contrast_Rate_of_Change

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: deltaVirialInitialized=.false.

  ! Variables to hold the tabulated virial overdensity data.
  integer                                        :: deltaVirialTableNumberPoints
  double precision,    allocatable, dimension(:) :: deltaVirialTableTime,deltaVirialTableDeltaVirial
  type(fgsl_interp)                              :: interpolationObject
  type(fgsl_interp_accel)                        :: interpolationAccelerator
  logical                                        :: resetInterpolation

  ! Name of virial overdensity method used.
  type(varying_string)                           :: virialDensityContrastMethod

  ! Pointer to the subroutine that tabulates the virial overdensity and template interface for that subroutine.
  procedure(Virial_Density_Contrast_Tabulate_Template), pointer :: Virial_Density_Contrast_Tabulate => null()
  abstract interface
     subroutine Virial_Density_Contrast_Tabulate_Template(time,deltaVirialTableNumberPoints,deltaVirialTime,deltaVirialDeltaVirial)
       double precision,                            intent(in)    :: time
       double precision, allocatable, dimension(:), intent(inout) :: deltaVirialTime,deltaVirialDeltaVirial
       integer,                                     intent(out)   :: deltaVirialTableNumberPoints
     end subroutine Virial_Density_Contrast_Tabulate_Template
  end interface
  
contains

  subroutine Virial_Density_Contrast_Initialize(time)
    !% Initializes the virial overdensity module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="virialDensityContrastMethod" type="moduleUse">
    include 'structure_formation.CDM.virial_overdensity.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in) :: time

    if (.not.deltaVirialInitialized) then
       ! Get the virial overdensity method parameter.
       !@ <inputParameter>
       !@   <name>virialDensityContrastMethod</name>
       !@   <defaultValue>spherical top hat</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Selects the method to be used for computing halo virial density contrasts.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('virialDensityContrastMethod',virialDensityContrastMethod,defaultValue='spherical top hat')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="virialDensityContrastMethod" type="code" action="subroutine">
       !#  <subroutineArgs>virialDensityContrastMethod,Virial_Density_Contrast_Tabulate</subroutineArgs>
       include 'structure_formation.CDM.virial_overdensity.inc'
       !# </include>
       if (.not.associated(Virial_Density_Contrast_Tabulate)) call Galacticus_Error_Report('Virial_Density_Contrast_Initialize','method ' &
            &//char(virialDensityContrastMethod)//' is unrecognized')
       ! Call routine to initialize the virial overdensity table.
       call Virial_Density_Contrast_Tabulate(time,deltaVirialTableNumberPoints,deltaVirialTableTime,deltaVirialTableDeltaVirial)
       ! Flag that the module is now initialized.
       deltaVirialInitialized=.true.
    end if
    return
  end subroutine Virial_Density_Contrast_Initialize
  
  subroutine Virial_Density_Contrast_Retabulate(time)
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: time
    logical                      :: remakeTable

    ! Check if we need to recompute our table.
    !$omp critical(Delta_Virial_Factor_Initialize)
    if (deltaVirialInitialized) then
       remakeTable=(time<deltaVirialTableTime(1).or.time>deltaVirialTableTime(deltaVirialTableNumberPoints))
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Virial_Density_Contrast_Initialize(time)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
    end if
    !$omp end critical(Delta_Virial_Factor_Initialize)
    return
  end subroutine Virial_Density_Contrast_Retabulate

  double precision function Halo_Virial_Density_Contrast(time,aExpansion,collapsing)
    !% Return the halo virial overdensity.
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: aExpansion,time
    logical,          intent(in), optional :: collapsing
    logical                                :: collapsingActual
    double precision                       :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Halo_Virial_Density_Contrast','only one argument can be specified')
       else
          timeActual=time
       end if
    else
       if (present(aExpansion)) then
          if (present(collapsing)) then
             collapsingActual=collapsing
          else
             collapsingActual=.false.
          end if
          timeActual=Cosmology_Age(aExpansion,collapsingActual)
       else
          call Galacticus_Error_Report('Halo_Virial_Density_Contrast','at least one argument must be given')
       end if
    end if

    ! Remake the table if necessary.
    call Virial_Density_Contrast_Retabulate(timeActual)

    ! Interpolate to get the expansion factor.
    Halo_Virial_Density_Contrast=Interpolate(deltaVirialTableNumberPoints,deltaVirialTableTime,deltaVirialTableDeltaVirial &
         &,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)
    return
  end function Halo_Virial_Density_Contrast

  double precision function Halo_Virial_Density_Contrast_Rate_of_Change(time,aExpansion,collapsing)
    !% Return the halo virial overdensity.
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: aExpansion,time
    logical,          intent(in), optional :: collapsing
    logical                                :: collapsingActual
    double precision                       :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Halo_Virial_Density_Contrast_Rate_of_Change','only one argument can be specified')
       else
          timeActual=time
       end if
    else
       if (present(aExpansion)) then
          if (present(collapsing)) then
             collapsingActual=collapsing
          else
             collapsingActual=.false.
          end if
          timeActual=Cosmology_Age(aExpansion,collapsingActual)
       else
          call Galacticus_Error_Report('Halo_Virial_Density_Contrast_Rate_of_Change','at least one argument must be given')
       end if
    end if

    ! Remake the table if necessary.
    call Virial_Density_Contrast_Retabulate(timeActual)

    ! Interpolate to get the expansion factor.
    Halo_Virial_Density_Contrast_Rate_of_Change=Interpolate_Derivative(deltaVirialTableNumberPoints,deltaVirialTableTime&
         &,deltaVirialTableDeltaVirial ,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)
    return
  end function Halo_Virial_Density_Contrast_Rate_of_Change

end module Virial_Density_Contrast
