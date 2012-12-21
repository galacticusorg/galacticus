!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  double precision            :: virialDensityContrastFixed

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
    type(varying_string)                         :: virialDensityContrastFixedType

    if (virialDensityContrastMethod == 'fixed') then
       ! Return a pointer to our tabulation function.
       Virial_Density_Contrast_Tabulate => Virial_Density_Fixed
       ! Get the fixed value to use.
       !@ <inputParameter>
       !@   <name>virialDensityContrastFixed</name>
       !@   <defaultValue>200</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The virial density contrast to use in the fixed value model.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("virialDensityContrastFixed"    ,virialDensityContrastFixed    ,defaultValue=200.0d0           )
       !@ <inputParameter>
       !@   <name>virialDensityContrastFixedType</name>
       !@   <defaultValue>criticalDensity</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The reference density to use in the fixed value virial density contrast model. Either of {\tt critical density} and {\tt mean density} are allowed.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("virialDensityContrastFixedType",virialDensityContrastFixedType,defaultValue='criticalDensity')
       select case (char(virialDensityContrastFixedType))
       case ("criticalDensity")
          densityType=densityTypeCritical
       case ("meanDensity"    )
          densityType=densityTypeMean
       case default
          call Galacticus_Error_Report('Virial_Density_Fixed_Initialize','[virialDensityContrastFixedType] must be either "critical density" or "mean density"')
       end select
    end if
    return
  end subroutine Virial_Density_Fixed_Initialize

  subroutine Virial_Density_Fixed(time,deltaVirialTable)
    !% Tabulate the virial density contrast assuming a fixed value.
    use Memory_Management
    use Cosmology_Functions
    use Numerical_Constants_Math
    use Numerical_Ranges
    use Tables
    implicit none
    double precision       , intent(in   )              :: time
    class           (table), intent(inout), allocatable :: deltaVirialTable
    integer                                             :: iTime,deltaTableNumberPoints
    double precision                                    :: densityContrast

    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)
    
    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(dlog10(deltaTableTimeMaximum/deltaTableTimeMinimum)&
         &*dble(deltaTableNPointsPerDecade))
    
    ! Deallocate table if currently allocated.
    if (allocated(deltaVirialTable)) then
       call deltaVirialTable%destroy()
       deallocate(deltaVirialTable)
    end if
    allocate(table1DLogarithmicLinear :: deltaVirialTable)
    select type (deltaVirialTable)
    type is (table1DLogarithmicLinear)
       ! Create the table.
       call deltaVirialTable%create(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints)
       ! Populate the table.
       do iTime=1,deltaTableNumberPoints
          densityContrast=virialDensityContrastFixed
          ! If the fixed value is defined with respect to the critical density, then translate it to be
          ! with respect to mean density.
          if (densityType == densityTypeCritical) densityContrast=densityContrast/Omega_Matter_Total(deltaVirialTable%x(iTime))
          call deltaVirialTable%populate(densityContrast,iTime)
       end do
    end select
    return
  end subroutine Virial_Density_Fixed

end module Virial_Densities_Fixed
