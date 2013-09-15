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

!% Contains a module that implements calculations of tidal fields acting on satellites.

module Satellites_Tidal_Fields
  !% Implements calculations of tidal fields acting on satellites.
  use ISO_Varying_String
  implicit none
  private
  public :: Satellite_Tidal_Field

  ! Flag to indicate if this module has been initialized.  
  logical              :: satellitesTidalFieldInitialized=.false.

  ! Name of ram pressure field method used.
  type(varying_string) :: satellitesTidalFieldMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Satellite_Tidal_Field), pointer :: Satellites_Tidal_Field_Get => null()

contains

  double precision function Satellite_Tidal_Field(thisNode)
    !% Return the tidal field acting on a satellite {\tt thisNode}.
    use Galacticus_Nodes
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satellitesTidalFieldMethod" type="moduleUse">
    include 'satellites.tidal_fields.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize if necessary.
    if (.not.satellitesTidalFieldInitialized) then
       !$omp critical(Satellites_Tidal_Field_Initialization) 
       if (.not.satellitesTidalFieldInitialized) then
          ! Get the cooling rate method parameter.
          !@ <inputParameter>
          !@   <name>satellitesTidalFieldMethod</name>
          !@   <defaultValue>sphericalSymmetry</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used when computing the tidal field acting on a satellite.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('satellitesTidalFieldMethod',satellitesTidalFieldMethod,defaultValue='sphericalSymmetry')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="satellitesTidalFieldMethod" type="functionCall" functionType="void">
          !#  <functionArgs>satellitesTidalFieldMethod,Satellites_Tidal_Field_Get</functionArgs>
          include 'satellites.tidal_fields.inc'
          !# </include>
          if (.not.associated(Satellites_Tidal_Field_Get)) call Galacticus_Error_Report('Satellites_Tidal_Field','method ' &
               &//char(satellitesTidalFieldMethod)//' is unrecognized')
          satellitesTidalFieldInitialized=.true.
       end if
       !$omp end critical(Satellites_Tidal_Field_Initialization) 
    end if

    ! Get the cooling rate using the selected method.
    Satellite_Tidal_Field=Satellites_Tidal_Field_Get(thisNode)

    return
  end function Satellite_Tidal_Field

end module Satellites_Tidal_Fields
