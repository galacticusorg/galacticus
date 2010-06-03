!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements the critical linear theory overdensity for halo collapse.

module Critical_Overdensity
  !% Implements the critical linear theory overdensity for halo collapse.
  use ISO_Varying_String
  use FGSL
  private
  public :: Critical_Overdensity_for_Collapse, Time_of_Collapse, Critical_Overdensity_State_Retrieve

  ! Flag to indicate if this module and tables have been initialized.
  logical                                        :: deltaCriticalInitialized=.false., tablesInitialized=.false.

  ! Variables to hold the tabulated critical overdensity data.
  integer                                        :: deltaCritTableNumberPoints
  double precision,    allocatable, dimension(:) :: deltaCritTableTime,deltaCritTableDeltaCrit,deltaCritReverseTableTime&
       &,deltaCritReverseTableDeltaCrit
  type(fgsl_interp)                              :: interpolationObject,reverseInterpolationObject
  type(fgsl_interp_accel)                        :: interpolationAccelerator,reverseInterpolationAccelerator
  logical                                        :: resetInterpolation,reverseResetInterpolation

  ! Name of critical overdensity method used.
  type(varying_string)                           :: criticalOverdensityMethod

  ! Pointer to the subroutine that tabulates the critical overdensity and template interface for that subroutine.
  procedure(Critical_Overdensity_Tabulate_Template), pointer :: Critical_Overdensity_Tabulate => null()
  abstract interface
     subroutine Critical_Overdensity_Tabulate_Template(time,deltaCritTableNumberPoints,deltaCritTime,deltaCritDeltaCrit)
       double precision,                            intent(in)    :: time
       double precision, allocatable, dimension(:), intent(inout) :: deltaCritTime,deltaCritDeltaCrit
       integer,                                     intent(out)   :: deltaCritTableNumberPoints
     end subroutine Critical_Overdensity_Tabulate_Template
  end interface
  
contains

  subroutine Critical_Overdensity_Initialize(time)
    !% Initializes the critical overdensity module.
    use Galacticus_Error
    use Input_Parameters
    use Array_Utilities
    use Memory_Management
    !# <include directive="criticalOverdensityMethod" type="moduleUse">
    include 'structure_formation.CDM.critical_overdensity.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in) :: time

    if (.not.deltaCriticalInitialized) then
       ! Get the critical overdensity method parameter.
       !@ <inputParameter>
       !@   <name>criticalOverdensityMethod</name>
       !@   <defaultValue>spherical top hat</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for critical overdensities for halo collapse.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('criticalOverdensityMethod',criticalOverdensityMethod,defaultValue='spherical top hat')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="criticalOverdensityMethod" type="code" action="subroutine">
       !#  <subroutineArgs>criticalOverdensityMethod,Critical_Overdensity_Tabulate</subroutineArgs>
       include 'structure_formation.CDM.critical_overdensity.inc'
       !# </include>
       if (.not.associated(Critical_Overdensity_Tabulate)) call Galacticus_Error_Report('Critical_Overdensity_Initialize','method ' &
            &//char(criticalOverdensityMethod)//' is unrecognized')
    end if
    ! Call routine to initialize the critical overdensity table.
    call Critical_Overdensity_Tabulate(time,deltaCritTableNumberPoints,deltaCritTableTime,deltaCritTableDeltaCrit)
    ! Create the reversed arrays.
    if (allocated(deltaCritReverseTableTime     )) call Dealloc_Array(deltaCritReverseTableTime     )
    if (allocated(deltaCritReverseTableDeltaCrit)) call Dealloc_Array(deltaCritReverseTableDeltaCrit)
    call Alloc_Array(deltaCritReverseTableTime     ,deltaCritTableNumberPoints,'deltaCritReverseTableTime'     )
    call Alloc_Array(deltaCritReverseTableDeltaCrit,deltaCritTableNumberPoints,'deltaCritReverseTableDeltaCrit')
    deltaCritReverseTableTime     =Array_Reverse(deltaCritTableTime     )
    deltaCritReverseTableDeltaCrit=Array_Reverse(deltaCritTableDeltaCrit)
    ! Flag that the module and tables are now initialized.
    deltaCriticalInitialized=.true.
    tablesInitialized       =.true.
    return
  end subroutine Critical_Overdensity_Initialize
  
  double precision function Critical_Overdensity_for_Collapse(time,aExpansion,collapsing)
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision,        intent(in), optional :: aExpansion,time
    logical,                 intent(in), optional :: collapsing
    logical                                       :: collapsingActual,remakeTable
    double precision                              :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Critical_Overdensity_for_Collapse','only one argument can be specified')
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
          call Galacticus_Error_Report('Critical_Overdensity_for_Collapse','at least one argument must be given')
       end if
    end if

    ! Check if we need to recompute our table.
    !$omp critical(DeltaCrit_Factor_Initialize)
    if (deltaCriticalInitialized.and.tablesInitialized) then
       remakeTable=(time<deltaCritTableTime(1).or.time>deltaCritTableTime(deltaCritTableNumberPoints))
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Critical_Overdensity_Initialize(timeActual)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
       reverseResetInterpolation=.true.
       deltaCriticalInitialized=.true.
    end if
    !$omp end critical(DeltaCrit_Factor_Initialize)

    ! Interpolate to get the expansion factor.
    Critical_Overdensity_for_Collapse=Interpolate(deltaCritTableNumberPoints,deltaCritTableTime,deltaCritTableDeltaCrit&
         &,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)
    return
  end function Critical_Overdensity_for_Collapse

  double precision function Time_of_Collapse(criticalOverdensity)
    !% Returns the time of collapse for a perturbation of linear theory overdensity {\tt criticalOverdensity}.
    use Numerical_Interpolation
    use Cosmology_Functions
    implicit none
    double precision, intent(in), optional :: criticalOverdensity
    double precision                       :: time

    ! Check if we need to recompute our table.
    !$omp critical(DeltaCrit_Factor_Initialize)
    do while (.not.deltaCriticalInitialized.or.criticalOverdensity>deltaCritTableDeltaCrit(1).or.criticalOverdensity&
         &<deltaCritTableDeltaCrit(deltaCritTableNumberPoints))
       if (.not.(deltaCriticalInitialized.and.tablesInitialized)) then
          time=Cosmology_Age(1.0d0)
       else
          if (criticalOverdensity>deltaCritTableDeltaCrit(1)) then
             time=0.5d0*deltaCritTableTime(1)
          else
             time=2.0d0*deltaCritTableTime(deltaCritTableNumberPoints)
          end if
       end if
       call Critical_Overdensity_Initialize(time)
       call Interpolate_Done(reverseInterpolationObject,reverseInterpolationAccelerator,reverseResetInterpolation)
       reverseResetInterpolation=.true.
       resetInterpolation=.true.
       deltaCriticalInitialized=.true.
    end do
    !$omp end critical(DeltaCrit_Factor_Initialize)
    
    ! Interpolate to get the expansion factor.
    Time_of_Collapse=Interpolate(deltaCritTableNumberPoints,deltaCritReverseTableDeltaCrit,deltaCritReverseTableTime &
         &,reverseInterpolationObject,reverseInterpolationAccelerator,criticalOverdensity,reset=reverseResetInterpolation)

    return
  end function Time_of_Collapse

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Critical_Overdensity_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Critical_Overdensity_State_Retrieve(stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use Memory_Management
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    deltaCritTableNumberPoints=0
    if (allocated(deltaCritTableTime     )) call Dealloc_Array(deltaCritTableTime     )
    if (allocated(deltaCritTableDeltaCrit)) call Dealloc_Array(deltaCritTableDeltaCrit)
    tablesInitialized=.false.
    return
  end subroutine Critical_Overdensity_State_Retrieve
  
end module Critical_Overdensity
