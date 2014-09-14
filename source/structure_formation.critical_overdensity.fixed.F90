!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of critical overdensity using a fixed value.

module Critical_Overdensities_Fixed
  implicit none
  private
  public :: Critical_Overdensity_Fixed_Initialize

  ! Value to use for critical overdenisty.
  double precision :: criticalOverdensityFixed

  ! Variables to hold the tabulated critical overdensity data.
  double precision            :: deltaTableTimeMaximum     =20.0d0, deltaTableTimeMinimum=1.0d0
  integer         , parameter :: deltaTableNPointsPerDecade=100

contains

  !# <criticalOverdensityMethod>
  !#  <unitName>Critical_Overdensity_Fixed_Initialize</unitName>
  !# </criticalOverdensityMethod>
  subroutine Critical_Overdensity_Fixed_Initialize(criticalOverdensityMethod,Critical_Overdensity_Contrast_Tabulate)
    !% Initializes the $\delta_{\rm c}$ calculation for the fixed module.
    use ISO_Varying_String
    use Input_Parameters
    use Numerical_Constants_Math
    implicit none
    type     (varying_string            ), intent(in   )          :: criticalOverdensityMethod
    procedure(Critical_Overdensity_Fixed), intent(inout), pointer :: Critical_Overdensity_Contrast_Tabulate

    if (criticalOverdensityMethod == 'fixed') then
       Critical_Overdensity_Contrast_Tabulate => Critical_Overdensity_Fixed     
       !@ <inputParameter>
       !@   <name>criticalOverdensityFixed</name>
       !@   <defaultValue>$(3/20)(12\pi)^{2/3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The value to use for the critical overdensity for collapse of dark matter halos when using a fixed value.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('criticalOverdensityFixed',criticalOverdensityFixed,defaultValue=(3.0d0/20.0d0)*(12.0d0*Pi)**(2.0d0/3.0d0))
    end if
    return
  end subroutine Critical_Overdensity_Fixed_Initialize

  subroutine Critical_Overdensity_Fixed(time,deltaTable)
    !% Tabulate the virial density contrast for the fixed module.
    use Tables
    use Linear_Growth
    implicit none
    double precision                      , intent(in   ) :: time
    class           (table1D), allocatable, intent(inout) :: deltaTable
    integer                                               :: deltaTableNumberPoints, iTime

    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)
    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(log10(deltaTableTimeMaximum/deltaTableTimeMinimum)*dble(deltaTableNPointsPerDecade))
    ! Create the table.
    if (allocated(deltaTable)) then
       call deltaTable%destroy()
       deallocate(deltaTable)
    end if
    allocate(table1DLogarithmicLinear :: deltaTable)
    select type (deltaTable)
    type is (table1DLogarithmicLinear)
       ! Create the table.
       call deltaTable%create(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints)
       ! Populate the table with our fixed value of the critical overdensity.
       do iTime=1,deltaTableNumberPoints
          call deltaTable%populate(criticalOverdensityFixed/Linear_Growth_Factor(deltaTable%x(iTime)),iTime)
       end do
    end select
    return
  end subroutine Critical_Overdensity_Fixed

end module Critical_Overdensities_Fixed
