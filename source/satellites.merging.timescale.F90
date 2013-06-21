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

!% Contains a module that implements calculations of merging timescales for satellites.

module Satellite_Merging_Timescales
  !% Implements calculations of merging timescales for satellites.
  use ISO_Varying_String
  implicit none
  private
  public :: Satellite_Time_Until_Merging

  ! Flag to indicate if this module has been initialized.
  logical                                          :: satelliteMergeTimescaleInitialized=.false.

  ! Name of satellite merging timescale method used.
  type     (varying_string              )          :: satelliteMergingMethod

  ! Pointer for function that will be called to assign merging times to satellites.
  procedure(Satellite_Time_Until_Merging), pointer :: Satellite_Time_Until_Merging_Get  =>null()

contains

  subroutine Satellite_Merging_Timescales_Initialize
    !% Initialize the satellite merging timescale module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satelliteMergingMethod" type="moduleUse">
    include 'satellite.merging.timescale.moduleUse.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.satelliteMergeTimescaleInitialized) then
       !$omp critical(Satellite_Merging_Timescales_Initialization)
       if (.not.satelliteMergeTimescaleInitialized) then

          ! Get the satellite merging timescale method.
          !@ <inputParameter>
          !@   <name>satelliteMergingMethod</name>
          !@   <defaultValue>Jiang2008</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used to compute satellite merging timescales.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteMergingMethod',satelliteMergingMethod,defaultValue='Jiang2008')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="satelliteMergingMethod" type="functionCall" functionType="void">
          !#  <functionArgs>satelliteMergingMethod,Satellite_Time_Until_Merging_Get</functionArgs>
          include 'satellite.merging.timescale.inc'
          !# </include>
          if (.not.associated(Satellite_Time_Until_Merging_Get))                                           &
               & call Galacticus_Error_Report(                                                             &
               &                              'Tree_Node_Methods_Satellite_Orbit_Initialize'             , &
               &                              'method '//char(satelliteMergingMethod)//' is unrecognized'  &
               &                             )
          ! Record that this module is now initialized.
          satelliteMergeTimescaleInitialized=.true.
       end if
       !$omp end critical(Satellite_Merging_Timescales_Initialization)
    end if
    return
  end subroutine Satellite_Merging_Timescales_Initialize

  double precision function Satellite_Time_Until_Merging(thisNode,thisOrbit)
    !% Return the satellite merging timescale for {\tt thisNode} (in units of Gyr).
    use Galacticus_Nodes
    use Kepler_Orbits
    implicit none
    type(treeNode   ), intent(inout), pointer :: thisNode
    type(keplerOrbit), intent(inout)          :: thisOrbit

    ! Initialize the module.
    call Satellite_Merging_Timescales_Initialize

    ! Get the cooling radius using the selected method.
    Satellite_Time_Until_Merging=Satellite_Time_Until_Merging_Get(thisNode,thisOrbit)
    return
  end function Satellite_Time_Until_Merging

end module Satellite_Merging_Timescales
