!% Contains a module which implements linear growth factor calculations.

module Linear_Growth
  !% Implements the virial overdensity for halos.
  use ISO_Varying_String
  use FGSL
  private
  public :: Linear_Growth_Factor, Linear_Growth_Initialize, Linear_Growth_State_Retrieve

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: linearGrowthInitialized=.false., tablesInitialized=.false.

  ! Variables to hold the tabulated linear growth data.
  integer                                        :: linearGrowthTableNumberPoints
  double precision,    allocatable, dimension(:) :: linearGrowthTableTime,linearGrowthTableFactor
  integer,             parameter                 :: linearGrowthTableNPointsPerDecade=1000
  type(fgsl_interp)                              :: interpolationObject
  type(fgsl_interp_accel)                        :: interpolationAccelerator
  logical                                        :: resetInterpolation

  ! Name of virial overdensity method used.
  type(varying_string)                           :: linearGrowthMethod

  ! Normalization types.
  integer,             parameter,   public       :: normalizePresentDay=0, normalizeMatterDominated=1
  double precision                               :: normalizationMatterDominated

  ! Pointer to the subroutine that tabulates the linear growth factor and template interface for that subroutine.
  procedure(Linear_Growth_Tabulate_Template), pointer :: Linear_Growth_Tabulate => null()
  abstract interface 
     subroutine Linear_Growth_Tabulate_Template(time,linearGrowthTableNumberPoints,linearGrowthTime,linearGrowthFactor&
          &,normalizationMatterDominated)
       double precision,                            intent(in)    :: time
       double precision, allocatable, dimension(:), intent(inout) :: linearGrowthTime,linearGrowthFactor
       integer,                                     intent(out)   :: linearGrowthTableNumberPoints
       double precision,                            intent(out)   :: normalizationMatterDominated
     end subroutine Linear_Growth_Tabulate_Template
  end interface
  
contains

  subroutine Linear_Growth_Initialize(time)
    !% Initializes the growth factor module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="linearGrowthMethod" type="moduleUse">
    include 'structure_formation.CDM.linear_growth.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in) :: time

    if (.not.linearGrowthInitialized) then
       ! Get the virial overdensity method parameter.
       !@ <inputParameter>
       !@   <name>linearGrowthMethod</name>
       !@   <defaultValue>simple</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of the linear growth factor.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('linearGrowthMethod',linearGrowthMethod,defaultValue='simple')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="linearGrowthMethod" type="code" action="subroutine">
       !#  <subroutineArgs>linearGrowthMethod,Linear_Growth_Tabulate</subroutineArgs>
       include 'structure_formation.CDM.linear_growth.inc'
       !# </include>
       if (.not.associated(Linear_Growth_Tabulate)) call Galacticus_Error_Report('Linear_Growth_Initialize','method ' &
            &//char(linearGrowthMethod)//' is unrecognized')
    end if
    ! Call routine to initialize the virial overdensity table.
    call Linear_Growth_Tabulate(time,linearGrowthTableNumberPoints,linearGrowthTableTime,linearGrowthTableFactor&
         &,normalizationMatterDominated)
    tablesInitialized=.true.
    ! Flag that the module is now initialized.
    linearGrowthInitialized=.true.
    return
  end subroutine Linear_Growth_Initialize
  
  double precision function Linear_Growth_Factor(time,aExpansion,collapsing,normalize)
    !% Return the linear growth factor.
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: aExpansion,time
    logical,          intent(in), optional :: collapsing
    integer,          intent(in), optional :: normalize
    integer                                :: normalizeActual
    logical                                :: collapsingActual,remakeTable
    double precision                       :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Linear_Growth_Factor','only one argument can be specified')
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
          call Galacticus_Error_Report('Linear_Growth_Factor','at least one argument must be given')
       end if
    end if

    ! Validate the input.
    if (.not.Cosmic_Time_Is_Valid(timeActual)) call Galacticus_Error_Report('Linear_Growth_Factor','cosmic time is&
         & invalid')

    ! Determine normalization method.
    if (present(normalize)) then
       normalizeActual=normalize
    else
       normalizeActual=normalizePresentDay
    end if

    ! Check if we need to recompute our table.
    !$omp critical(Linear_Growth_Initialize)
    if (linearGrowthInitialized.and.tablesInitialized) then
       remakeTable=(time<linearGrowthTableTime(1).or.time>linearGrowthTableTime(linearGrowthTableNumberPoints))
    else
       call Linear_Growth_Initialize(timeActual)
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Linear_Growth_Tabulate(timeActual,linearGrowthTableNumberPoints,linearGrowthTableTime,linearGrowthTableFactor &
            &,normalizationMatterDominated)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
       linearGrowthInitialized=.true.
    end if
    !$omp end critical(Linear_Growth_Initialize)

    ! Interpolate to get the expansion factor.
    Linear_Growth_Factor=Interpolate(linearGrowthTableNumberPoints,linearGrowthTableTime,linearGrowthTableFactor &
         &,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)

    ! Normalize.
    select case (normalizeActual)
    case (normalizeMatterDominated)
       Linear_Growth_Factor=Linear_Growth_Factor*normalizationMatterDominated
    end select

    return
  end function Linear_Growth_Factor

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Linear_Growth_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Linear_Growth_State_Retrieve(stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use Memory_Management
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    linearGrowthTableNumberPoints=0
    if (allocated(linearGrowthTableTime  )) call Dealloc_Array(linearGrowthTableTime  )
    if (allocated(linearGrowthTableFactor)) call Dealloc_Array(linearGrowthTableFactor)
    tablesInitialized=.false.
    return
  end subroutine Linear_Growth_State_Retrieve
  
end module Linear_Growth
