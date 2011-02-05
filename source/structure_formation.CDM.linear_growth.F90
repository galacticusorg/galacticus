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


!% Contains a module which implements linear growth factor calculations.

module Linear_Growth
  !% Implements the virial overdensity for halos.
  use ISO_Varying_String
  use FGSL
  private
  public :: Linear_Growth_Factor, Linear_Growth_Factor_Logarithmic_Derivative, Linear_Growth_State_Retrieve

  ! Flag to indicate if this module has been initialized.  
  logical                                                    :: linearGrowthInitialized=.false., tablesInitialized=.false.

  ! Variables to hold the tabulated linear growth data.
  integer                                                    :: linearGrowthTableNumberPoints
  double precision,    allocatable, dimension(:)             :: linearGrowthTableTime,linearGrowthTableWavenumber
  double precision,    allocatable, dimension(:,:,:), target :: linearGrowthTableFactor
  integer,             parameter                             :: linearGrowthTableNPointsPerDecade=1000
  type(fgsl_interp)                                          :: interpolationObject
  type(fgsl_interp_accel)                                    :: interpolationAccelerator,interpolationAcceleratorWavenumber
  logical                                                    :: resetInterpolation,resetInterpolationWavenumber

  ! Name of virial overdensity method used.
  type(varying_string)                                       :: linearGrowthMethod

  ! Normalization types.
  integer,             parameter,   public                   :: normalizePresentDay=0, normalizeMatterDominated=1
  double precision,                 dimension(3)             :: normalizationMatterDominated

  ! Component types.
  integer,             parameter,   public                   ::  linearGrowthComponentDarkMatter=1 &
       &                                                        ,linearGrowthComponentBaryons   =2 &
       &,linearGrowthComponentRadiation =3

  ! Pointer to the subroutine that tabulates the linear growth factor and template interface for that subroutine.
  procedure(Linear_Growth_Tabulate_Template), pointer :: Linear_Growth_Tabulate => null()
  abstract interface 
     subroutine Linear_Growth_Tabulate_Template(time,linearGrowthTableNumberPoints,linearGrowthTime,linearGrowthTableWavenumber&
          &,linearGrowthTableFactor ,normalizationMatterDominated)
       double precision,                                intent(in)    :: time
       double precision, allocatable, dimension(:),     intent(inout) :: linearGrowthTime,linearGrowthTableWavenumber
       double precision, allocatable, dimension(:,:,:), intent(inout) :: linearGrowthTableFactor
       integer,                                         intent(out)   :: linearGrowthTableNumberPoints
       double precision,              dimension(3),     intent(out)   :: normalizationMatterDominated
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
    call Linear_Growth_Tabulate(time,linearGrowthTableNumberPoints,linearGrowthTableTime,linearGrowthTableWavenumber&
         &,linearGrowthTableFactor ,normalizationMatterDominated)
    tablesInitialized=.true.
    ! Flag that the module is now initialized.
    linearGrowthInitialized=.true.
    return
  end subroutine Linear_Growth_Initialize
  
  double precision function Linear_Growth_Factor(time,aExpansion,collapsing,normalize,component,wavenumber)
    !% Return the linear growth factor.
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional       :: aExpansion,time
    logical,          intent(in), optional       :: collapsing
    integer,          intent(in), optional       :: normalize,component
    double precision, intent(in), optional       :: wavenumber
    double precision, pointer,    dimension(:)   :: thisLinearGrowthFactor
    integer,                      dimension(0:1) :: iWavenumber
    double precision,             dimension(0:1) :: hWavenumber
    integer                                      :: normalizeActual,componentActual,jWavenumber
    logical                                      :: collapsingActual,remakeTable
    double precision                             :: timeActual

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

    ! Determine the component type to use.
    if (present(component)) then
       componentActual=component
    else
       ! Assume dark matter growth factor is required by default.
       componentActual=linearGrowthComponentDarkMatter
    end if

    ! Check if we need to recompute our table.
    !$omp critical(Linear_Growth_Initialize)
    if (linearGrowthInitialized.and.tablesInitialized) then
       remakeTable=(timeActual<linearGrowthTableTime(1).or.timeActual>linearGrowthTableTime(linearGrowthTableNumberPoints))
    else
       call Linear_Growth_Initialize(timeActual)
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Linear_Growth_Tabulate(timeActual,linearGrowthTableNumberPoints,linearGrowthTableTime &
            &,linearGrowthTableWavenumber,linearGrowthTableFactor ,normalizationMatterDominated)
       call Interpolate_Done( interpolationObject     =interpolationObject                &
            &                ,interpolationAccelerator=interpolationAccelerator           &
            &                ,reset                   =resetInterpolation                 &
            &               )
       call Interpolate_Done( interpolationAccelerator=interpolationAcceleratorWavenumber &
            &                ,reset                   =resetInterpolationWavenumber       &
            &               )
       resetInterpolation          =.true.
       resetInterpolationWavenumber=.true.
       linearGrowthInitialized     =.true.
    end if
    !$omp end critical(Linear_Growth_Initialize)

    ! Determine which wavenumbers to use.
    call Interpolate_In_Wavenumber(iWavenumber,hWavenumber,wavenumber)

    ! Loop over wavenumbers and compute growth factor.
    Linear_Growth_Factor=0.0d0
    do jWavenumber=0,1
       
       ! Select the appropriate component.
       thisLinearGrowthFactor => linearGrowthTableFactor(componentActual,iWavenumber(jWavenumber),:)

       ! Interpolate to get the expansion factor.
       Linear_Growth_Factor=Linear_Growth_Factor+Interpolate(linearGrowthTableNumberPoints,linearGrowthTableTime &
            &,thisLinearGrowthFactor,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)&
            &*hWavenumber(jWavenumber)
    end do

    ! Normalize.
    select case (normalizeActual)
    case (normalizeMatterDominated)
       Linear_Growth_Factor=Linear_Growth_Factor*normalizationMatterDominated(componentActual)
    end select

    return
  end function Linear_Growth_Factor

  double precision function Linear_Growth_Factor_Logarithmic_Derivative(time,aExpansion,collapsing,component,wavenumber)
    !% Return the logarithmic derivative of the linear growth factor with respect to expansion factor., $\d \ln D / \d \ln a$.
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional       :: aExpansion,time
    logical,          intent(in), optional       :: collapsing
    integer,          intent(in), optional       :: component
    double precision, intent(in), optional       :: wavenumber
    double precision, pointer,    dimension(:)   :: thisLinearGrowthFactor
    integer,                      dimension(0:1) :: iWavenumber
    double precision,             dimension(0:1) :: hWavenumber
    integer                                      :: componentActual,jWavenumber
    logical                                      :: collapsingActual
    double precision                             :: timeActual,expansionFactor,linearGrowthFactor &
         &,linearGrowthFactorTimeDerivative

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Linear_Growth_Factor_Logarithmic_Derivative','only one argument can be specified')
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
          call Galacticus_Error_Report('Linear_Growth_Factor_Logarithmic_Derivative','at least one argument must be given')
       end if
    end if

    ! Validate the input.
    if (.not.Cosmic_Time_Is_Valid(timeActual)) call Galacticus_Error_Report('Linear_Growth_Factor_Logarithmic_Derivative','cosmic time is&
         & invalid')

    ! Determine the component type to use.
    if (present(component)) then
       componentActual=component
    else
       ! Assume dark matter growth factor is required by default.
       componentActual=linearGrowthComponentDarkMatter
    end if

    ! Get the linear growth factor (this will automatically build the tabulation if necessary).
    linearGrowthFactor=Linear_Growth_Factor(timeActual,component=componentActual,wavenumber=wavenumber)

    ! Determine which wavenumbers to use.
    call Interpolate_In_Wavenumber(iWavenumber,hWavenumber,wavenumber)

    ! Loop over wavenumbers and compute growth factor.
    linearGrowthFactorTimeDerivative=0.0d0
    do jWavenumber=0,1
       
       ! Select the appropriate component.
       thisLinearGrowthFactor => linearGrowthTableFactor(componentActual,iWavenumber(jWavenumber),:)

       ! Interpolate to get the expansion factor.
       linearGrowthFactorTimeDerivative=linearGrowthFactorTimeDerivative+Interpolate_Derivative(linearGrowthTableNumberPoints&
            &,linearGrowthTableTime ,thisLinearGrowthFactor,interpolationObject,interpolationAccelerator,timeActual,reset&
            &=resetInterpolation)*hWavenumber(jWavenumber)
    end do
    
    ! Get the expansion factor.
    expansionFactor=Expansion_Factor(timeActual)
    
    ! Construct the logarithmic derivative with respect to expansion factor.
    Linear_Growth_Factor_Logarithmic_Derivative=linearGrowthFactorTimeDerivative/linearGrowthFactor&
         &/Expansion_Rate(expansionFactor)

    return
  end function Linear_Growth_Factor_Logarithmic_Derivative

  subroutine Interpolate_In_Wavenumber(iWavenumber,hWavenumber,wavenumber)
    !% Find interpolating factors in the wavenumber dimesion for linear growth factor calculations.
    use Numerical_Interpolation
    implicit none
    integer,          intent(out), dimension(0:1) :: iWavenumber
    double precision, intent(out), dimension(0:1) :: hWavenumber
    double precision, intent(in),  optional       :: wavenumber

    if (present(wavenumber).and.size(linearGrowthTableWavenumber) > 1) then
       iWavenumber(0)=Interpolate_Locate              (                                    &
            &                                           size(linearGrowthTableWavenumber)  &
            &                                          ,     linearGrowthTableWavenumber   &
            &                                          ,interpolationAcceleratorWavenumber &
            &                                          ,wavenumber                         &
            &                                          ,reset=resetInterpolationWavenumber &
            &                                         )
       hWavenumber=Interpolate_Linear_Generate_Factors(                                    &
            &                                           size(linearGrowthTableWavenumber)  &
            &                                          ,     linearGrowthTableWavenumber   &
            &                                          ,iWavenumber(0)                     &
            &                                          ,wavenumber                         &
            &                                         )
       iWavenumber(1)=iWavenumber(0)+1
    else
       ! No wavenumber was specified, so use the largest scale (smallest wavenumber) tabulated.
       iWavenumber(0:1)=[1    ,1    ]
       hWavenumber(0:1)=[1.0d0,0.0d0]
    end if
    return
  end subroutine Interpolate_In_Wavenumber

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
    if (allocated(linearGrowthTableTime      )) call Dealloc_Array(linearGrowthTableTime      )
    if (allocated(linearGrowthTableWavenumber)) call Dealloc_Array(linearGrowthTableWavenumber)
    if (allocated(linearGrowthTableFactor    )) call Dealloc_Array(linearGrowthTableFactor    )
    tablesInitialized=.false.
    return
  end subroutine Linear_Growth_State_Retrieve
  
end module Linear_Growth
