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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements calculations of black hole binary separation growth rate.

module Black_Hole_Binary_Separations
  !% Implements calculations of black hole binary separation growth rate.
  use ISO_Varying_String
  implicit none
  private
  public :: Black_Hole_Binary_Separation_Growth_Rate

  ! Flag to indicate if this module has been initialized.  
  logical                                                      :: blackHoleBinarySeparationGrowthRateInitialized=.false.  
  
  ! Name of mass movement method used.                                                                                                                     
  type     (varying_string                          )          :: blackHoleBinarySeparationGrowthRateMethod               
  
  ! Pointer to the subroutine that returns descriptors for mass movement.                                                                                                                     
  procedure(Black_Hole_Binary_Separation_Growth_Rate), pointer :: Black_Hole_Binary_Separation_Growth_Rate_Get  =>null()  
                                                                                                                       
contains

  double precision function Black_Hole_Binary_Separation_Growth_Rate(thisBlackHoleComponent)
    !% Computes the separation growth rate of a black hole binary in units of Mpc/Gyr.
    use Galacticus_Error
    use Input_Parameters
    use Galacticus_Nodes
    !# <include directive="blackHoleBinarySeparationGrowthRateMethod" type="moduleUse">
    include 'black_holes.binary.separation_growth_rate.modules.inc'
    !# </include>
    implicit none
    class(nodeComponentBlackHole), intent(inout), pointer :: thisBlackHoleComponent  
                                                                                  
    if (.not.blackHoleBinarySeparationGrowthRateInitialized) then
       !$omp critical(blackHoleBinarySeparationGrowthRateInitialize)
       if (.not.blackHoleBinarySeparationGrowthRateInitialized) then
          ! Get the binary black hole separation growth rate method parameter.
          !@ <inputParameter>
          !@   <name>blackHoleBinarySeparationGrowthRateMethod</name>
          !@   <defaultValue>null</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing the separation growth rate of black hole binaries.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('blackHoleBinarySeparationGrowthRateMethod',blackHoleBinarySeparationGrowthRateMethod,defaultValue='null')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="blackHoleBinarySeparationGrowthRateMethod" type="functionCall" functionType="void">
          !#  <functionArgs>blackHoleBinarySeparationGrowthRateMethod,Black_Hole_Binary_Separation_Growth_Rate_Get</functionArgs>
          include 'black_holes.binaries.separation_growth_rate.inc'
          !# </include>
          if (.not.associated(Black_Hole_Binary_Separation_Growth_Rate_Get)) call Galacticus_Error_Report('Black_Hole_Binary_Separation_Growth_Rate','method ' &
               &//char(blackHoleBinarySeparationGrowthRateMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          blackHoleBinarySeparationGrowthRateInitialized=.true.
       end if
       !$omp end critical(blackHoleBinarySeparationGrowthRateInitialize)
    end if

    ! Call the routine to do the calculation.
    Black_Hole_Binary_Separation_Growth_Rate=Black_Hole_Binary_Separation_Growth_Rate_Get(thisBlackHoleComponent)

    return
  end function Black_Hole_Binary_Separation_Growth_Rate
  
end module Black_Hole_Binary_Separations
