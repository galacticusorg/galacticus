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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!% Contains a module that implements calculations of the mass loss rate due to
!% tidal stripping for satellites.

module Satellite_Tidal_Stripping
  !% Implements calculations of tidal stripping for satellites.
  use ISO_Varying_String
  implicit none
  private
  public :: Satellite_Tidal_Stripping_Rate

  ! Flag to indicate if this module has been initialized.  
  logical                                            :: satelliteTidalStrippingRateInitialized=.false.

  ! Name of satellite tidal stripping method used.
  type     (varying_string                )          :: satelliteTidalStrippingMethod

  ! Pointer for function that will be called to assign tidal stripping rates to satellites.
  procedure(Satellite_Tidal_Stripping_Rate), pointer :: Satellite_Tidal_Stripping_Rate_Get => null()
  
contains

  subroutine Satellite_Tidal_Stripping_Initialize
    !% Initialize the satellite tidal stripping rate module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satelliteTidalStrippingMethod" type="moduleUse">
    include 'satellite.tidal.stripping.rate.moduleUse.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.satelliteTidalStrippingRateInitialized) then
       !$omp critical(Satellite_Tidal_Stripping_Rate_Initialization) 
       if (.not.satelliteTidalStrippingRateInitialized) then
          
          ! Get the satellite tidal stripping method.
          !@ <inputParameter>
          !@   <name>satelliteTidalStrippingMethod</name>
          !@   <defaultValue>Zentner2005</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used to compute satellite tidal stripping rate.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteTidalStrippingMethod',satelliteTidalStrippingMethod,defaultValue='Zentner2005')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="satelliteTidalStrippingMethod" type="functionCall" functionType="void">
          !#  <functionArgs>satelliteTidalStrippingMethod,Satellite_Tidal_Stripping_Rate_Get</functionArgs>
          include 'satellite.tidal.stripping.rate.inc'
          !# </include>
          if (.not.associated(Satellite_Tidal_Stripping_Rate_Get))                                                &
               & call Galacticus_Error_Report(                                                                    &
               &                              'Satellite_Tidal_Stripping_Rate_Initialize'                       , &
               &                              'method '//char(satelliteTidalStrippingMethod)//' is unrecognized'  &
               &                             )
          ! Record that this module is now initialized.
          satelliteTidalStrippingRateInitialized=.true.
       end if
       !$omp end critical(Satellite_Tidal_Stripping_Rate_Initialization) 
    end if
    return
  end subroutine Satellite_Tidal_Stripping_Initialize

  double precision function Satellite_Tidal_Stripping_Rate(thisNode)
    !% Return the satellite mass loss rate due to tidal stripping for {\normalfont \ttfamily thisNode} (in units of $M_\odot$/Gyr).
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Satellite_Tidal_Stripping_Initialize

    ! Get the mass loss rate due to tidal stripping using the selected method.
    Satellite_Tidal_Stripping_Rate=Satellite_Tidal_Stripping_Rate_Get(thisNode)
    return
  end function Satellite_Tidal_Stripping_Rate

end module Satellite_Tidal_Stripping
