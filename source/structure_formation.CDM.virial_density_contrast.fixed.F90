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


!% Contains a module which implements calculations of virial overdensity using a fixed value.

module Virial_Densities_Fixed
  !% Implements calculations of virial overdensity using a fixed value.
  implicit none
  private
  public :: Virial_Density_Fixed_Initialize

  ! The type of reference density to use.
  integer                     :: densityType
  integer,          parameter :: densityTypeCritical=0
  integer,          parameter :: densityTypeMean    =1

  ! The fixed overdensity to use.
  double precision            :: virialDensityConstrastFixed

  ! Variables to hold the tabulated critical overdensity data.
  double precision            :: deltaTableTimeMinimum=1.0d0, deltaTableTimeMaximum=20.0d0
  integer,          parameter :: deltaTableNPointsPerDecade=100

contains

  !# <virialDensityContrastMethod>
  !#  <unitName>Virial_Density_Fixed_Initialize</unitName>
  !# </virialDensityContrastMethod>
  subroutine Virial_Density_Fixed_Initialize(virialDensityContrastMethod,Virial_Density_Contrast_Tabulate)
    !% Initializes the $\Delta_{\rm vir}$ calculation for the fixed value implementation.
    use Input_Parameters
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type(varying_string),          intent(in)    :: virialDensityContrastMethod
    procedure(),          pointer, intent(inout) :: Virial_Density_Contrast_Tabulate
    type(varying_string)                         :: virialDensityConstrastFixedType

    if (virialDensityContrastMethod == 'fixed') then
       ! Return a pointer to our tabulation function.
       Virial_Density_Contrast_Tabulate => Virial_Density_Fixed
       ! Get the fixed value to use.
       !@ <inputParameter>
       !@   <name>virialDensityConstrastFixed</name>
       !@   <defaultValue>200</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The virial density contrast to use in the fixed value model.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("virialDensityConstrastFixed"    ,virialDensityConstrastFixed    ,defaultValue=200.0d0           )
       !@ <inputParameter>
       !@   <name>virialDensityConstrastFixedType</name>
       !@   <defaultValue>critical density</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The reference density to use in the fixed value virial density contrast model. Either of {\tt critical density} and {\tt mean density} are allowed.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("virialDensityConstrastFixedType",virialDensityConstrastFixedType,defaultValue='critical density')
       select case (char(virialDensityConstrastFixedType))
       case ("critical density")
          densityType=densityTypeCritical
       case ("mean density"    )
          densityType=densityTypeMean
       case default
          call Galacticus_Error_Report('Virial_Density_Fixed_Initialize','[virialDensityConstrastFixedType] must be either "critical density" or "mean density"')
       end select
    end if
    return
  end subroutine Virial_Density_Fixed_Initialize

  subroutine Virial_Density_Fixed(time,deltaTableNumberPoints,deltaTableTime,deltaTableDelta)
    !% Tabulate the virial density contrast assuming a fixed value.
    use Memory_Management
    use Cosmology_Functions
    use Numerical_Constants_Math
    use Numerical_Ranges
    implicit none
    double precision, intent(in)                               :: time
    integer,          intent(out)                              :: deltaTableNumberPoints
    double precision, intent(inout), allocatable, dimension(:) :: deltaTableTime,deltaTableDelta
    integer                                                    :: iTime

    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)
    
    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(dlog10(deltaTableTimeMaximum/deltaTableTimeMinimum)&
         &*dble(deltaTableNPointsPerDecade))
    
    ! Deallocate arrays if currently allocated.
    if (allocated(deltaTableTime )) call Dealloc_Array(deltaTableTime )
    if (allocated(deltaTableDelta)) call Dealloc_Array(deltaTableDelta)
    ! Allocate the arrays to current required size.
    call Alloc_Array(deltaTableTime ,[deltaTableNumberPoints])
    call Alloc_Array(deltaTableDelta,[deltaTableNumberPoints])
    
    ! Create the tabulation.
    deltaTableTime  =Make_Range(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints,rangeTypeLogarithmic)

    ! Set the fixed virial density contrast.
    deltaTableDelta=virialDensityConstrastFixed

    ! If the fixed value is defined with respect to the critical density, then translate it to be with respect to mean density.
    if (densityType == densityTypeCritical) then
       do iTime=1,deltaTableNumberPoints
          deltaTableDelta(iTime)= deltaTableDelta(iTime)/Omega_Matter_Total(deltaTableTime(iTime))
       end do
    end if
    return
  end subroutine Virial_Density_Fixed

end module Virial_Densities_Fixed
