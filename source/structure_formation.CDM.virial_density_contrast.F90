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


!% Contains a module which implements the virial overdensity for halos.

module Virial_Density_Contrast
  !% Implements the virial overdensity for halos.
  use ISO_Varying_String
  use FGSL
  private
  public :: Halo_Virial_Density_Contrast, Halo_Virial_Density_Contrast_Rate_of_Change, Virial_Density_Contrast_State_Retrieve

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: deltaVirialInitialized=.false., tablesInitialized=.false.

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
       ! Flag that the module is now initialized.
       deltaVirialInitialized=.true.
       tablesInitialized     =.true.
    end if

    ! Call routine to initialize the virial overdensity table.
    call Virial_Density_Contrast_Tabulate(time,deltaVirialTableNumberPoints,deltaVirialTableTime,deltaVirialTableDeltaVirial)
    return
  end subroutine Virial_Density_Contrast_Initialize
  
  subroutine Virial_Density_Contrast_Retabulate(time)
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: time
    logical                      :: remakeTable

    ! Check if we need to recompute our table.
    !$omp critical(Delta_Virial_Factor_Initialize)
    if (deltaVirialInitialized.and.tablesInitialized) then
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
    !$omp critical(Halo_Virial_Density_Contrast_Interpolate)
    Halo_Virial_Density_Contrast=Interpolate(deltaVirialTableNumberPoints,deltaVirialTableTime,deltaVirialTableDeltaVirial &
         &,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)
    !$omp end critical(Halo_Virial_Density_Contrast_Interpolate)
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
    !$omp critical(Halo_Virial_Density_Contrast_Interpolate)
    Halo_Virial_Density_Contrast_Rate_of_Change=Interpolate_Derivative(deltaVirialTableNumberPoints,deltaVirialTableTime&
         &,deltaVirialTableDeltaVirial ,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)
    !$omp end critical(Halo_Virial_Density_Contrast_Interpolate)
    return
  end function Halo_Virial_Density_Contrast_Rate_of_Change

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Virial_Density_Contrast_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Virial_Density_Contrast_State_Retrieve(stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use Memory_Management
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    deltaVirialTableNumberPoints=0
    if (allocated(deltaVirialTableTime       )) call Dealloc_Array(deltaVirialTableTime       )
    if (allocated(deltaVirialTableDeltaVirial)) call Dealloc_Array(deltaVirialTableDeltaVirial)
    tablesInitialized=.false.
    return
  end subroutine Virial_Density_Contrast_State_Retrieve
  
end module Virial_Density_Contrast
