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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!% Contains a module that implements calculations of the acceleration due to dynamical friction for satellites.

module Satellite_Dynamical_Friction
  !% Implements calculations of dynamical friction for satellites.
  use ISO_Varying_String
  implicit none
  private
  public :: Satellite_Dynamical_Friction_Acceleration

  ! Flag to indicate if this module has been initialized.  
  logical                                                       :: satelliteDynamicalFrictionAccelerationInitialized=.false.

  ! Name of satellite dynamical friction method used.
  type     (varying_string                           )          :: satelliteDynamicalFrictionMethod

  ! Pointer for function that will be called to assign dynamical friction
  ! accelerations to satellites.
  procedure(Satellite_Dynamical_Friction_Acceleration), pointer :: Satellite_Dynamical_Friction_Acceleration_Get => null()
  
contains

  subroutine Satellite_Dynamical_Friction_Acceleration_Initialize
    !% Initialize the satellite dynamical friction acceleration module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satelliteDynamicalFrictionMethod" type="moduleUse">
    include 'satellite.dynamical.friction.acceleration.moduleUse.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.satelliteDynamicalFrictionAccelerationInitialized) then
       !$omp critical(Satellite_Dynamical_Friction_Acceleration_Initialization) 
       if (.not.satelliteDynamicalFrictionAccelerationInitialized) then
          ! Get the satellite dynamical friction method.
          !@ <inputParameter>
          !@   <name>satelliteDynamicalFrictionMethod</name>
          !@   <defaultValue>Chandrasekhar1943</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used to compute satellite dynamical friction acceleration.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteDynamicalFrictionMethod',satelliteDynamicalFrictionMethod,defaultValue='Chandrasekhar1943')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="satelliteDynamicalFrictionMethod" type="functionCall" functionType="void">
          !#  <functionArgs>satelliteDynamicalFrictionMethod,Satellite_Dynamical_Friction_Acceleration_Get</functionArgs>
          include 'satellite.dynamical.friction.acceleration.inc'
          !# </include>
          if (.not.associated(Satellite_Dynamical_Friction_Acceleration_Get))                                        &
               & call Galacticus_Error_Report(                                                                       &
               &                              'Satellite_Dynamical_Friction_Acceleration_Initialize'               , &
               &                              'method '//char(satelliteDynamicalFrictionMethod)//' is unrecognized'  &
               &                             )
          ! Record that this module is now initialized.
          satelliteDynamicalFrictionAccelerationInitialized=.true.
       end if
       !$omp end critical(Satellite_Dynamical_Friction_Acceleration_Initialization) 
    end if
    return
  end subroutine Satellite_Dynamical_Friction_Acceleration_Initialize

  function Satellite_Dynamical_Friction_Acceleration(thisNode)
    !% Return the satellite acceleration due to dynamical friction for {\normalfont \ttfamily thisNode} (in units of km/s/Gyr).
    use Galacticus_Nodes
    implicit none
    double precision             , dimension(3)           :: Satellite_Dynamical_Friction_Acceleration
    type            (treeNode   ), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Satellite_Dynamical_Friction_Acceleration_Initialize

    ! Get the acceleration due to dynamical friction using the selected method.
    Satellite_Dynamical_Friction_Acceleration=Satellite_Dynamical_Friction_Acceleration_Get(thisNode)
    return
  end function Satellite_Dynamical_Friction_Acceleration
  
end module Satellite_Dynamical_Friction
