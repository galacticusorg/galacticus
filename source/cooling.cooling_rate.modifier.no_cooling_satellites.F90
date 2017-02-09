!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module which implements a cut off in the cooling rate in satellites.

module Cooling_Rates_Modifier_Satellite
  !% Implements a cut off in the cooling rate at given redshift and virial velocity.
  implicit none
  private
  public :: Cooling_Rate_Modifier_Satellite

  ! Record of whether module has been initialized.
  logical :: moduleInitialized                   =.false.

  ! Parameters controlling this modifier.
  logical :: coolingModifierNoCoolingInSatellites

contains

  !# <coolingRateModifierMethod>
  !#  <unitName>Cooling_Rate_Modifier_Satellite</unitName>
  !# </coolingRateModifierMethod>
  subroutine Cooling_Rate_Modifier_Satellite(thisNode,coolingRate)
    !% Modify cooling rates by truncating them to zero in satellites.
    use Input_Parameters
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout) :: thisNode
    double precision          , intent(inout) :: coolingRate

    if (.not.moduleInitialized) then
       !$omp critical (Cooling_Rate_Modifier_Satellite_Initialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>coolingModifierNoCoolingInSatellites</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Specifies whether to switch off cooling in satellite nodes.
          !@   </description>
          !@   <type>logical</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("coolingModifierNoCoolingInSatellites",coolingModifierNoCoolingInSatellites,defaultValue=.false.)
          ! Record that the module is now initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical (Cooling_Rate_Modifier_Satellite_Initialize)
    end if
    ! Apply the modifier.
    if (coolingModifierNoCoolingInSatellites.and.thisNode%isSatellite()) coolingRate=0.0d0
    return
  end subroutine Cooling_Rate_Modifier_Satellite

end module Cooling_Rates_Modifier_Satellite
