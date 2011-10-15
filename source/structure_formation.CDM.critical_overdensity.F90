!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements the critical linear theory overdensity for halo collapse.

module Critical_Overdensity
  !% Implements the critical linear theory overdensity for halo collapse.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Critical_Overdensity_for_Collapse, Critical_Overdensity_for_Collapse_Time_Gradient, Time_of_Collapse,&
       & Critical_Overdensity_State_Retrieve, Critical_Overdensity_Collapsing_Mass

  ! Flag to indicate if this module and tables have been initialized.
  logical                                        :: deltaCriticalInitialized=.false., massScalingInitialized=.false., tablesInitialized=.false.

  ! Variables to hold the tabulated critical overdensity data.
  integer                                        :: deltaCritTableNumberPoints
  double precision,    allocatable, dimension(:) :: deltaCritTableTime,deltaCritTableDeltaCrit,deltaCritReverseTableTime&
       &,deltaCritReverseTableDeltaCrit
  type(fgsl_interp)                              :: interpolationObject,reverseInterpolationObject
  type(fgsl_interp_accel)                        :: interpolationAccelerator,reverseInterpolationAccelerator
  logical                                        :: resetInterpolation,reverseResetInterpolation

  ! Name of critical overdensity method used.
  type(varying_string)                           :: criticalOverdensityMethod,criticalOverdensityMassScalingMethod

  ! Global variable used in root finding.
  double precision                               :: collapseTime

  ! Pointer to the subroutine that tabulates the critical overdensity and template interface for that subroutine.
  procedure(Critical_Overdensity_Tabulate_Template), pointer :: Critical_Overdensity_Tabulate => null()
  abstract interface
     subroutine Critical_Overdensity_Tabulate_Template(time,deltaCritTableNumberPoints,deltaCritTime,deltaCritDeltaCrit)
       double precision,                            intent(in)    :: time
       double precision, allocatable, dimension(:), intent(inout) :: deltaCritTime,deltaCritDeltaCrit
       integer,                                     intent(out)   :: deltaCritTableNumberPoints
     end subroutine Critical_Overdensity_Tabulate_Template
  end interface

  ! Pointer to the mass scaling function.
  procedure(Critical_Overdensity_Mass_Scaling_Template), pointer :: Critical_Overdensity_Mass_Scaling_Get => null()
  abstract interface
     double precision function Critical_Overdensity_Mass_Scaling_Template(mass)
       double precision,                            intent(in)    :: mass
     end function Critical_Overdensity_Mass_Scaling_Template
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
       !@   <defaultValue>sphericalTopHat</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for critical overdensities for halo collapse.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('criticalOverdensityMethod',criticalOverdensityMethod,defaultValue='sphericalTopHat')
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
    if (.not.Array_Is_Monotonic(deltaCritTableDeltaCrit,direction=directionDecreasing)) call Galacticus_Error_Report('Critical_Overdensity_Initialize','critical overdensity must be monotonically decreasing with time')
    ! Create the reversed arrays.
    if (allocated(deltaCritReverseTableTime     )) call Dealloc_Array(deltaCritReverseTableTime     )
    if (allocated(deltaCritReverseTableDeltaCrit)) call Dealloc_Array(deltaCritReverseTableDeltaCrit)
    call Alloc_Array(deltaCritReverseTableTime     ,[deltaCritTableNumberPoints])
    call Alloc_Array(deltaCritReverseTableDeltaCrit,[deltaCritTableNumberPoints])
    deltaCritReverseTableTime     =Array_Reverse(deltaCritTableTime     )
    deltaCritReverseTableDeltaCrit=Array_Reverse(deltaCritTableDeltaCrit)
    ! Flag that the module and tables are now initialized.
    deltaCriticalInitialized=.true.
    tablesInitialized       =.true.
    return
  end subroutine Critical_Overdensity_Initialize

  double precision function Critical_Overdensity_Collapsing_Mass(time,aExpansion,collapsing)
    !% Return the mass scale just collapsing at the given cosmic time.
    use Root_Finder
    use FGSL
    use Cosmology_Functions
    implicit none
    double precision,        intent(in), optional :: aExpansion,time
    logical,                 intent(in), optional :: collapsing
    double precision,        parameter            :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6,massGuess=1.0d13
    type(fgsl_function),     save                 :: rootFunction
    type(fgsl_root_fsolver), save                 :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision                              :: massMinimum,massMaximum
    type(c_ptr)                                   :: parameterPointer

    ! Get the critical overdensity for collapse at this epoch.
    if (present(time)) then
       collapseTime=time
    else
       collapseTime=Cosmology_Age(aExpansion,collapsing)
    end if
    
    ! Find mass at which the root-variance (sigma) equals this critical overdensity.
    massMinimum=massGuess
    do while (Collapsing_Mass_Root(massMinimum,parameterPointer) <= 0.0d0)
       massMinimum=0.5d0*massMinimum
    end do
    massMaximum=massGuess
    do while (Collapsing_Mass_Root(massMaximum,parameterPointer) >= 0.0d0)
       massMaximum=2.0d0*massMaximum
    end do
    Critical_Overdensity_Collapsing_Mass=Root_Find(massMinimum,massMaximum,Collapsing_Mass_Root,parameterPointer,rootFunction&
         &,rootFunctionSolver,toleranceAbsolute,toleranceRelative)
    
    return
  end function Critical_Overdensity_Collapsing_Mass

  function Collapsing_Mass_Root(mass,parameterPointer) bind(c)
    !% Function used in finding the mass of halo just collapsing at a given cosmic epoch.
    use CDM_Power_Spectrum
    real(c_double)          :: Collapsing_Mass_Root
    real(c_double), value   :: mass
    type(c_ptr),    value   :: parameterPointer
    
    Collapsing_Mass_Root=sigma_CDM(mass)-Critical_Overdensity_for_Collapse(time=collapseTime,mass=mass)
    return
  end function Collapsing_Mass_Root

  double precision function Critical_Overdensity_for_Collapse(time,aExpansion,collapsing,mass)
    !% Return the linear theory critical overdensity for collapse at the given cosmic time.
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision,        intent(in), optional :: aExpansion,time,mass
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
    !$omp critical(Critical_Overdensity_for_Collapse_Interp)
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

    ! Interpolate to get the critical overdensity for collapse.
    Critical_Overdensity_for_Collapse=Interpolate(deltaCritTableNumberPoints,deltaCritTableTime,deltaCritTableDeltaCrit&
         &,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)
    !$omp end critical(Critical_Overdensity_for_Collapse_Interp)

    ! Scale by a mass dependent factor if necessary.
    if (present(mass)) Critical_Overdensity_for_Collapse=Critical_Overdensity_for_Collapse*Critical_Overdensity_Mass_Scaling(mass)
    return
  end function Critical_Overdensity_for_Collapse

  double precision function Critical_Overdensity_for_Collapse_Time_Gradient(time,aExpansion,collapsing,mass)
    !% Return the derivative with respect to time of the linear theory critical overdensity for collapse at the given cosmic time.
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision,        intent(in), optional :: aExpansion,time,mass
    logical,                 intent(in), optional :: collapsing
    logical                                       :: collapsingActual,remakeTable
    double precision                              :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Critical_Overdensity_for_Collapse_Time_Gradient','only one argument can be specified')
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
          call Galacticus_Error_Report('Critical_Overdensity_for_Collapse_Time_Gradient','at least one argument must be given')
       end if
    end if

    ! Check if we need to recompute our table.
    !$omp critical(Critical_Overdensity_for_Collapse_Interp)
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

    ! Interpolate to get the derivative.
    Critical_Overdensity_for_Collapse_Time_Gradient=Interpolate_Derivative(deltaCritTableNumberPoints,deltaCritTableTime&
         &,deltaCritTableDeltaCrit ,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)
    !$omp end critical(Critical_Overdensity_for_Collapse_Interp)

    ! Scale by a mass dependent factor if necessary.
    if (present(mass)) Critical_Overdensity_for_Collapse_Time_Gradient=Critical_Overdensity_for_Collapse_Time_Gradient*Critical_Overdensity_Mass_Scaling(mass)
    return
  end function Critical_Overdensity_for_Collapse_Time_Gradient

  double precision function Time_of_Collapse(criticalOverdensity,mass)
    !% Returns the time of collapse for a perturbation of linear theory overdensity {\tt criticalOverdensity}.
    use Numerical_Interpolation
    use Cosmology_Functions
    implicit none
    double precision, intent(in)           :: criticalOverdensity
    double precision, intent(in), optional :: mass
    double precision                       :: time,criticalOverdensityActual

    ! Scale by a mass dependent factor if necessary.
    if (present(mass)) then
       criticalOverdensityActual=criticalOverdensity/Critical_Overdensity_Mass_Scaling(mass)
    else
       criticalOverdensityActual=criticalOverdensity
    end if

    ! Check if we need to recompute our table.
    !$omp critical(Critical_Overdensity_for_Collapse_Interp)
    do while (.not.deltaCriticalInitialized.or.criticalOverdensityActual>deltaCritTableDeltaCrit(1).or.criticalOverdensityActual&
         &<deltaCritTableDeltaCrit(deltaCritTableNumberPoints))
       if (.not.(deltaCriticalInitialized.and.tablesInitialized)) then
          time=Cosmology_Age(1.0d0)
       else
          if (criticalOverdensityActual>deltaCritTableDeltaCrit(1)) then
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
    
    ! Interpolate to get the expansion factor.
    Time_of_Collapse=Interpolate(deltaCritTableNumberPoints,deltaCritReverseTableDeltaCrit,deltaCritReverseTableTime &
         &,reverseInterpolationObject,reverseInterpolationAccelerator,criticalOverdensityActual,reset=reverseResetInterpolation)
    !$omp end critical(Critical_Overdensity_for_Collapse_Interp)

    return
  end function Time_of_Collapse

  double precision function Critical_Overdensity_Mass_Scaling(mass)
    !% Return a multiplicative, mass-dependent factor by which the critical overdensity should be scaled.
    implicit none
    double precision, intent(in) :: mass

    ! Ensure the mass scaling method is initialized.
    call Critical_Overdensity_Mass_Scaling_Initialize    
  
    ! Perform the calculation.
    Critical_Overdensity_Mass_Scaling=Critical_Overdensity_Mass_Scaling_Get(mass)
    return
  end function Critical_Overdensity_Mass_Scaling

  subroutine Critical_Overdensity_Mass_Scaling_Initialize
    !% Initializes the critical overdensity mass scaling method.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="criticalOverdensityMassScalingMethod" type="moduleUse">
    include 'structure_formation.CDM.critical_overdensity.mass_scaling.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Critical_Overdensity_for_Collapse_Mass_Scaling_Initialize)
    if (.not.massScalingInitialized) then
       ! Get the critical overdensity mass scaling method parameter.
       !@ <inputParameter>
       !@   <name>criticalOverdensityMassScalingMethod</name>
       !@   <defaultValue>null</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for scaling critical overdensities for halo collapse with mass.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('criticalOverdensityMassScalingMethod',criticalOverdensityMassScalingMethod,defaultValue='null')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="criticalOverdensityMassScalingMethod" type="code" action="subroutine">
       !#  <subroutineArgs>criticalOverdensityMassScalingMethod,Critical_Overdensity_Mass_Scaling_Get</subroutineArgs>
       include 'structure_formation.CDM.critical_overdensity.mass_scaling.inc'
       !# </include>
       if (.not.associated(Critical_Overdensity_Mass_Scaling_Get)) call Galacticus_Error_Report('Critical_Overdensity_Initialize','method ' &
            &//char(criticalOverdensityMassScalingMethod)//' is unrecognized')
       ! Flag that mass scaling has been initialized.
       massScalingInitialized=.true.
    end if
    !$omp end critical(Critical_Overdensity_for_Collapse_Mass_Scaling_Initialize)
    return
  end subroutine Critical_Overdensity_Mass_Scaling_Initialize

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
