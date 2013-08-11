!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use Tables, only : table1D
  implicit none
  private
  public :: Critical_Overdensity_for_Collapse, Critical_Overdensity_for_Collapse_Time_Gradient, Time_of_Collapse,&
       & Critical_Overdensity_State_Retrieve, Critical_Overdensity_Collapsing_Mass, Critical_Overdensity_Mass_Scaling,&
       & Critical_Overdensity_Mass_Scaling_Gradient

  ! Flag to indicate if this module and tables have been initialized.
  logical                                                               :: deltaCriticalInitialized            =.false., massScalingInitialized   =.false., &
       &                                                                   tablesInitialized                   =.false.
  !$omp threadprivate(tablesInitialized)
  ! Variables to hold the tabulated critical overdensity data.
  class           (table1D                               ), allocatable :: deltaCritTable                              , deltaCritTableReversed
  !$omp threadprivate(deltaCritTable,deltaCritTableReversed)
  ! Name of critical overdensity method used.
  type            (varying_string                        )              :: criticalOverdensityMassScalingMethod        , criticalOverdensityMethod

  ! Global variable used in root finding.
  double precision                                                      :: collapseTime
  !$omp threadprivate(collapseTime)
  ! Pointer to the subroutine that tabulates the critical overdensity and template interface for that subroutine.
  procedure       (Critical_Overdensity_Tabulate_Template), pointer     :: Critical_Overdensity_Tabulate       =>null()
  abstract interface
     subroutine Critical_Overdensity_Tabulate_Template(time,deltaCritTable)
       import table1D
       double precision                      , intent(in   ) :: time
       class           (table1D), allocatable, intent(inout) :: deltaCritTable
     end subroutine Critical_Overdensity_Tabulate_Template
  end interface

  ! Pointer to the mass scaling function.
  procedure(Critical_Overdensity_Mass_Scaling_Template), pointer :: Critical_Overdensity_Mass_Scaling_Get         =>null()
  procedure(Critical_Overdensity_Mass_Scaling_Template), pointer :: Critical_Overdensity_Mass_Scaling_Gradient_Get=>null()
  abstract interface
     double precision function Critical_Overdensity_Mass_Scaling_Template(mass)
       double precision, intent(in   ) :: mass
     end function Critical_Overdensity_Mass_Scaling_Template
  end interface

contains

  subroutine Critical_Overdensity_Initialize(time)
    !% Initializes the critical overdensity module.
    use Galacticus_Error
    use Input_Parameters
    use Array_Utilities
    !# <include directive="criticalOverdensityMethod" type="moduleUse">
    include 'structure_formation.critical_overdensity.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in   ) :: time

    if (.not.deltaCriticalInitialized) then
       !$omp critical (Critical_Overdensity_Initialize)
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
          !# <include directive="criticalOverdensityMethod" type="functionCall" functionType="void">
          !#  <functionArgs>criticalOverdensityMethod,Critical_Overdensity_Tabulate</functionArgs>
          include 'structure_formation.critical_overdensity.inc'
          !# </include>
          if (.not.associated(Critical_Overdensity_Tabulate)) call Galacticus_Error_Report('Critical_Overdensity_Initialize','method ' &
               &//char(criticalOverdensityMethod)//' is unrecognized')
          deltaCriticalInitialized=.true.
       end if
       !$omp end critical (Critical_Overdensity_Initialize)
    end if
    ! Call routine to initialize the critical overdensity table.
    call Critical_Overdensity_Tabulate(time,deltaCritTable)
    if (.not.deltaCritTable%isMonotonic(direction=directionDecreasing)) call Galacticus_Error_Report('Critical_Overdensity_Initialize','critical overdensity must be monotonically decreasing with time')
    ! Create the reversed arrays.
    call deltaCritTable%reverse(deltaCritTableReversed)
    ! Flag that the module and tables are now initialized.
    tablesInitialized=.true.
    return
  end subroutine Critical_Overdensity_Initialize

  double precision function Critical_Overdensity_Collapsing_Mass(time,aExpansion,collapsing)
    !% Return the mass scale just collapsing at the given cosmic time.
    use Root_Finder
    use Cosmology_Functions
    implicit none
    double precision            , intent(in   ), optional :: aExpansion              , time
    logical                     , intent(in   ), optional :: collapsing
    double precision            , parameter               :: massGuess        =1.0d13, toleranceAbsolute=0.0d0, &
         &                                                   toleranceRelative=1.0d-6
    type            (rootFinder), save                    :: finder
    !$omp threadprivate(finder)
    ! Get the critical overdensity for collapse at this epoch.
    if (present(time)) then
       collapseTime=time
    else
       collapseTime=Cosmology_Age(aExpansion,collapsing)
    end if
    ! Find mass at which the root-variance (sigma) equals this critical overdensity.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(Collapsing_Mass_Root               )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
       call finder%rangeExpand (                                                             &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
    end if
    Critical_Overdensity_Collapsing_Mass=finder%find(rootGuess=massGuess)
    return
  end function Critical_Overdensity_Collapsing_Mass

  double precision function Collapsing_Mass_Root(mass)
    !% Function used in finding the mass of halo just collapsing at a given cosmic epoch.
    use Power_Spectra
    double precision, intent(in   ) :: mass

    Collapsing_Mass_Root=Cosmological_Mass_Root_Variance(mass)-Critical_Overdensity_for_Collapse(time=collapseTime,mass=mass)
    return
  end function Collapsing_Mass_Root

  double precision function Critical_Overdensity_for_Collapse(time,aExpansion,collapsing,mass)
    !% Return the linear theory critical overdensity for collapse at the given cosmic time.
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: aExpansion      , mass       , time
    logical         , intent(in   ), optional :: collapsing
    logical                                   :: collapsingActual, remakeTable
    double precision                          :: timeActual

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
    if (tablesInitialized) then
       remakeTable=(time<deltaCritTable%x(1).or.time>deltaCritTable%x(-1))
    else
       remakeTable=.true.
    end if
    if (remakeTable) call Critical_Overdensity_Initialize(timeActual)

    ! Interpolate to get the critical overdensity for collapse.
    Critical_Overdensity_for_Collapse=deltaCritTable%interpolate(timeActual)

    ! Scale by a mass dependent factor if necessary.
    if (present(mass)) Critical_Overdensity_for_Collapse=Critical_Overdensity_for_Collapse*Critical_Overdensity_Mass_Scaling(mass)
    return
  end function Critical_Overdensity_for_Collapse

  double precision function Critical_Overdensity_for_Collapse_Time_Gradient(time,aExpansion,collapsing,mass)
    !% Return the derivative with respect to time of the linear theory critical overdensity for collapse at the given cosmic time.
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: aExpansion      , mass       , time
    logical         , intent(in   ), optional :: collapsing
    logical                                   :: collapsingActual, remakeTable
    double precision                          :: timeActual

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
    if (tablesInitialized) then
       remakeTable=(time<deltaCritTable%x(1).or.time>deltaCritTable%x(-1))
    else
       remakeTable=.true.
    end if
    if (remakeTable) call Critical_Overdensity_Initialize(timeActual)

    ! Interpolate to get the derivative.
    Critical_Overdensity_for_Collapse_Time_Gradient=deltaCritTable%interpolateGradient(timeActual)

    ! Scale by a mass dependent factor if necessary.
    if (present(mass)) Critical_Overdensity_for_Collapse_Time_Gradient=Critical_Overdensity_for_Collapse_Time_Gradient*Critical_Overdensity_Mass_Scaling(mass)
    return
  end function Critical_Overdensity_for_Collapse_Time_Gradient

  double precision function Time_of_Collapse(criticalOverdensity,mass)
    !% Returns the time of collapse for a perturbation of linear theory overdensity {\tt criticalOverdensity}.
    use Cosmology_Functions
    implicit none
    double precision, intent(in   )           :: criticalOverdensity
    double precision, intent(in   ), optional :: mass
    double precision                          :: criticalOverdensityActual, time

    ! Scale by a mass dependent factor if necessary.
    if (present(mass)) then
       criticalOverdensityActual=criticalOverdensity/Critical_Overdensity_Mass_Scaling(mass)
    else
       criticalOverdensityActual=criticalOverdensity
    end if

    ! Check if we need to recompute our table.
    if (.not.tablesInitialized) call Critical_Overdensity_Initialize(Cosmology_Age(1.0d0))
    do while (criticalOverdensityActual<deltaCritTableReversed%x(1).or.criticalOverdensityActual&
         &>deltaCritTableReversed%x(-1))
       if (criticalOverdensityActual>deltaCritTableReversed%x(-1)) then
          time=0.5d0*deltaCritTable%x( 1)
       else
          time=2.0d0*deltaCritTable%x(-1)
       end if
       call Critical_Overdensity_Initialize(time)
    end do

    ! Interpolate to get the expansion factor.
    Time_of_Collapse=deltaCritTableReversed%interpolate(criticalOverdensityActual)
    return
  end function Time_of_Collapse

  double precision function Critical_Overdensity_Mass_Scaling(mass)
    !% Return a multiplicative, mass-dependent factor by which the critical overdensity should be scaled.
    implicit none
    double precision, intent(in   ) :: mass

    ! Ensure the mass scaling method is initialized.
    call Critical_Overdensity_Mass_Scaling_Initialize

    ! Perform the calculation.
    Critical_Overdensity_Mass_Scaling=Critical_Overdensity_Mass_Scaling_Get(mass)
    return
  end function Critical_Overdensity_Mass_Scaling

  double precision function Critical_Overdensity_Mass_Scaling_Gradient(mass)
    !% Return the gradient with mass of a multiplicative, mass-dependent factor by which the critical overdensity should be scaled.
    implicit none
    double precision, intent(in   ) :: mass

    ! Ensure the mass scaling method is initialized.
    call Critical_Overdensity_Mass_Scaling_Initialize

    ! Perform the calculation.
    Critical_Overdensity_Mass_Scaling_Gradient=Critical_Overdensity_Mass_Scaling_Gradient_Get(mass)
    return
  end function Critical_Overdensity_Mass_Scaling_Gradient

  subroutine Critical_Overdensity_Mass_Scaling_Initialize
    !% Initializes the critical overdensity mass scaling method.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="criticalOverdensityMassScalingMethod" type="moduleUse">
    include 'structure_formation.critical_overdensity.mass_scaling.modules.inc'
    !# </include>
    implicit none

    if (.not.massScalingInitialized) then
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
          !# <include directive="criticalOverdensityMassScalingMethod" type="functionCall" functionType="void">
          !#  <functionArgs>criticalOverdensityMassScalingMethod,Critical_Overdensity_Mass_Scaling_Get,Critical_Overdensity_Mass_Scaling_Gradient_Get</functionArgs>
          include 'structure_formation.critical_overdensity.mass_scaling.inc'
          !# </include>
          if (.not.(associated(Critical_Overdensity_Mass_Scaling_Get).and.associated(Critical_Overdensity_Mass_Scaling_Gradient_Get))) call Galacticus_Error_Report('Critical_Overdensity_Initialize','method ' &
               &//char(criticalOverdensityMassScalingMethod)//' is unrecognized')
          ! Flag that mass scaling has been initialized.
          massScalingInitialized=.true.
       end if
       !$omp end critical(Critical_Overdensity_for_Collapse_Mass_Scaling_Initialize)
    end if
    return
  end subroutine Critical_Overdensity_Mass_Scaling_Initialize

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Critical_Overdensity_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Critical_Overdensity_State_Retrieve(stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    tablesInitialized=.false.
    return
  end subroutine Critical_Overdensity_State_Retrieve

end module Critical_Overdensity
