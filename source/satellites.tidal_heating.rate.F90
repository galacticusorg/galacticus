!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that implements calculations of the tidal heating rate for satellites.

module Satellite_Tidal_Heating
  !% Implements calculations of tidal heating for satellites.
  use ISO_Varying_String
  implicit none
  private
  public :: Satellite_Tidal_Heating_Rate

  ! Flag to indicate if this module has been initialized.  
  logical                                          :: satelliteTidalHeatingInitialized=.false.

  ! Name of satellite tidal heating method used.
  type     (varying_string              )          :: satelliteTidalHeatingMethod

  ! Pointer for function that will be called to assign tidal heating rates
  ! to satellites.
  procedure(Satellite_Tidal_Heating_Rate), pointer :: Satellite_Tidal_Heating_Rate_Get => null()
  
contains

  subroutine Satellite_Tidal_Heating_Rate_Initialize()
    !% Initialize the satellite tidal heating rate module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satelliteTidalHeatingMethod" type="moduleUse">
    include 'satellite.tidal.heating.rate.moduleUse.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.satelliteTidalHeatingInitialized) then
       !$omp critical(Satellite_Tidal_Heating_Initialization) 
       if (.not.satelliteTidalHeatingInitialized) then
          
          ! Get the satellite tidal heating method.
          !@ <inputParameter>
          !@   <name>satelliteTidalHeatingMethod</name>
          !@   <defaultValue>Gnedin1999</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used to compute satellite tidal heating rate.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteTidalHeatingMethod',satelliteTidalHeatingMethod,defaultValue='Gnedin1999')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="satelliteTidalHeatingMethod" type="functionCall" functionType="void">
          !#  <functionArgs>satelliteTidalHeatingMethod,Satellite_Tidal_Heating_Rate_Get</functionArgs>
          include 'satellite.tidal.heating.rate.inc'
          !# </include>
          if (.not.associated(Satellite_Tidal_Heating_Rate_Get))                                                &
               & call Galacticus_Error_Report(                                                                  &
               &                              'Satellite_Tidal_Heating_Rate_Initialize'                       , &
               &                              'method '//char(satelliteTidalHeatingMethod)//' is unrecognized'  &
               &                             )
          ! Record that this module is now initialized.
          satelliteTidalHeatingInitialized=.true.
       end if
       !$omp end critical(Satellite_Tidal_Heating_Initialization) 
    end if
    return
  end subroutine Satellite_Tidal_Heating_Rate_Initialize

  double precision function Satellite_Tidal_Heating_Rate(thisNode)
    !% Return the satellite tidal heating rate for {\normalfont \ttfamily thisNode} (in units of (km/s/Mpc)$^2$/Gyr).
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Satellite_Tidal_Heating_Rate_Initialize

    ! Get the tidal heating rate using the selected method.
    Satellite_Tidal_Heating_Rate=Satellite_Tidal_Heating_Rate_Get(thisNode)
    return
  end function Satellite_Tidal_Heating_Rate
  
end module Satellite_Tidal_Heating
