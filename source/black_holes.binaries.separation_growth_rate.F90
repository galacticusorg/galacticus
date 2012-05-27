!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
  type(varying_string)                                         :: blackHoleBinarySeparationGrowthRateMethod

  ! Pointer to the subroutine that returns descriptors for mass movement.
  procedure(Black_Hole_Binary_Separation_Growth_Rate), pointer :: Black_Hole_Binary_Separation_Growth_Rate_Get => null()
  
contains

  double precision function Black_Hole_Binary_Separation_Growth_Rate(thisNode)
    !% Computes the separation growth rate of a black hole binary in units of Mpc/Gyr.
    use Galacticus_Error
    use Input_Parameters
    use Tree_Nodes
    !# <include directive="blackHoleBinarySeparationGrowthRateMethod" type="moduleUse">
    include 'black_holes.binary.separation_growth_rate.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

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
          !# <include directive="blackHoleBinarySeparationGrowthRateMethod" type="code" action="subroutine">
          !#  <subroutineArgs>blackHoleBinarySeparationGrowthRateMethod,Black_Hole_Binary_Separation_Growth_Rate_Get</subroutineArgs>
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
    Black_Hole_Binary_Separation_Growth_Rate=Black_Hole_Binary_Separation_Growth_Rate_Get(thisNode)

    return
  end function Black_Hole_Binary_Separation_Growth_Rate
  
end module Black_Hole_Binary_Separations
