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


!% Contains a module that implements calculations of ram pressure stripping of hot halos.

module Hot_Halo_Ram_Pressure_Stripping
  !% Implements calculations of ram pressure stripping of hot halos.
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Hot_Halo_Ram_Pressure_Stripping_Radius

  ! Flag to indicate if this module has been initialized.  
  logical              :: hotHaloRamPressureStrippingInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: hotHaloRamPressureStrippingMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Hot_Halo_Ram_Pressure_Stripping_Radius), pointer :: Hot_Halo_Ram_Pressure_Stripping_Radius_Get => null()

contains

  double precision function Hot_Halo_Ram_Pressure_Stripping_Radius(thisNode)
    !% Return the ram pressure stripping radius for the hot halo of {\tt thisNode} (in units of Mpc).
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="hotHaloRamPressureStrippingMethod" type="moduleUse">
    include 'hot_halo.ram_pressure_stripping.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize if necessary.
    if (.not.hotHaloRamPressureStrippingInitialized) then
       !$omp critical(Hot_Halo_Ram_Pressure_Stripping_Initialization) 
       if (.not.hotHaloRamPressureStrippingInitialized) then
          ! Get the cooling rate method parameter.
          !@ <inputParameter>
          !@   <name>hotHaloRamPressureStrippingMethod</name>
          !@   <defaultValue>Font2008</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used when computing ram pressure stripping of hot halos.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloRamPressureStrippingMethod',hotHaloRamPressureStrippingMethod,defaultValue='virialRadius')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="hotHaloRamPressureStrippingMethod" type="code" action="subroutine">
          !#  <subroutineArgs>hotHaloRamPressureStrippingMethod,Hot_Halo_Ram_Pressure_Stripping_Radius_Get</subroutineArgs>
          include 'hot_halo.ram_pressure_stripping.inc'
          !# </include>
          if (.not.associated(Hot_Halo_Ram_Pressure_Stripping_Radius_Get)) call Galacticus_Error_Report('Hot_Halo_Ram_Pressure_Stripping_Radius','method ' &
               &//char(hotHaloRamPressureStrippingMethod)//' is unrecognized')
          hotHaloRamPressureStrippingInitialized=.true.
       end if
       !$omp end critical(Hot_Halo_Ram_Pressure_Stripping_Initialization) 
    end if

    ! Get the cooling rate using the selected method.
    Hot_Halo_Ram_Pressure_Stripping_Radius=Hot_Halo_Ram_Pressure_Stripping_Radius_Get(thisNode)

    return
  end function Hot_Halo_Ram_Pressure_Stripping_Radius

end module Hot_Halo_Ram_Pressure_Stripping
