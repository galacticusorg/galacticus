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

!% Contains a module which implements the virial overdensity for halos.

module Virial_Density_Contrast
  !% Implements the virial overdensity for halos.
  use ISO_Varying_String
  use Tables
  implicit none
  private
  public :: Halo_Virial_Density_Contrast, Halo_Virial_Density_Contrast_Rate_of_Change, Virial_Density_Contrast_State_Retrieve

  ! Flag to indicate if this module has been initialized.
  logical                                                           :: deltaVirialInitialized          =.false.

  ! Variables to hold the tabulated virial overdensity data.
  logical                                                           :: tablesInitialized               =.false.
  class    (table1D                                  ), allocatable :: deltaVirialTable
  !$omp threadprivate(deltaVirialTable,tablesInitialized)
  ! Name of virial overdensity method used.
  type     (varying_string                           )              :: virialDensityContrastMethod

  ! Pointer to the subroutine that tabulates the virial overdensity and template interface for that subroutine.
  procedure(Virial_Density_Contrast_Tabulate_Template), pointer     :: Virial_Density_Contrast_Tabulate=>null()
  abstract interface
     subroutine Virial_Density_Contrast_Tabulate_Template(time,deltaVirialTable)
       import table1D
       double precision         , intent(in   ) :: time
       class           (table1D), intent(inout) :: deltaVirialTable
     end subroutine Virial_Density_Contrast_Tabulate_Template
  end interface

contains

  subroutine Virial_Density_Contrast_Initialize(time)
    !% Initializes the virial overdensity module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="virialDensityContrastMethod" type="moduleUse">
    include 'structure_formation.virial_overdensity.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in   ) :: time

    ! Initialize the module if necessary.
    if (.not.deltaVirialInitialized) then
       !$omp critical (Virial_Density_Contrast_Initialize)
       if (.not.deltaVirialInitialized) then
          ! Get the virial overdensity method parameter.
          !@ <inputParameter>
          !@   <name>virialDensityContrastMethod</name>
          !@   <defaultValue>sphericalTopHat</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Selects the method to be used for computing halo virial density contrasts.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('virialDensityContrastMethod',virialDensityContrastMethod,defaultValue='sphericalTopHat')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="virialDensityContrastMethod" type="functionCall" functionType="void">
          !#  <functionArgs>virialDensityContrastMethod,Virial_Density_Contrast_Tabulate</functionArgs>
          include 'structure_formation.virial_overdensity.inc'
          !# </include>
          if (.not.associated(Virial_Density_Contrast_Tabulate)) call Galacticus_Error_Report('Virial_Density_Contrast_Initialize','method ' &
               &//char(virialDensityContrastMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          deltaVirialInitialized=.true.
       end if
       !$omp end critical (Virial_Density_Contrast_Initialize)
    end if
    return
  end subroutine Virial_Density_Contrast_Initialize

  subroutine Virial_Density_Contrast_Retabulate(time)
    !% Recompute the look-up tables for virial density contrast.
    implicit none
    double precision, intent(in   ) :: time
    logical                         :: remakeTable

    ! Ensure that the module is initialized.
    call Virial_Density_Contrast_Initialize(time)

    ! Check if we need to recompute our table.
    if (tablesInitialized) then
       remakeTable=(time<deltaVirialTable%x(1).or.time>deltaVirialTable%x(-1))
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Virial_Density_Contrast_Tabulate(time,deltaVirialTable)
       tablesInitialized=.true.
    end if
    return
  end subroutine Virial_Density_Contrast_Retabulate

  double precision function Halo_Virial_Density_Contrast(time,aExpansion,collapsing)
    !% Return the halo virial overdensity.
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: aExpansion      , time
    logical         , intent(in   ), optional :: collapsing
    logical                                   :: collapsingActual
    double precision                          :: timeActual

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
    Halo_Virial_Density_Contrast=deltaVirialTable%interpolate(timeActual)
    return
  end function Halo_Virial_Density_Contrast

  double precision function Halo_Virial_Density_Contrast_Rate_of_Change(time,aExpansion,collapsing)
    !% Return the halo virial overdensity.
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: aExpansion      , time
    logical         , intent(in   ), optional :: collapsing
    logical                                   :: collapsingActual
    double precision                          :: timeActual

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
    Halo_Virial_Density_Contrast_Rate_of_Change=deltaVirialTable%interpolateGradient(timeActual)
    return
  end function Halo_Virial_Density_Contrast_Rate_of_Change

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Virial_Density_Contrast_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Virial_Density_Contrast_State_Retrieve(stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    tablesInitialized=.false.
    return
  end subroutine Virial_Density_Contrast_State_Retrieve

end module Virial_Density_Contrast
