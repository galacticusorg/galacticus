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

!% Contains a module which implements calculations for progenitor properties for merger remnant calculations.

module Satellite_Merging_Remnant_Sizes_Progenitors
  !% Implements calculations for progenitor properties for merger remnant calculations.
  use ISO_Varying_String
  implicit none
  private
  public :: Satellite_Merging_Remnant_Progenitor_Properties

  ! Flag to indicate if this module has been initialized.  
  logical                                            :: satelliteMergingRemnantProgenitorPropertiesInitialized=.false.

  ! Name of mass movement method used.
  type(varying_string)                               :: satelliteMergingRemnantProgenitorPropertiesMethod

  ! Pointer to the subroutine that returns properties of progenitors.
  procedure(Satellite_Merging_Remnant_Progenitor_Properties), pointer :: Satellite_Merging_Remnant_Progenitor_Properties_Get => null()
  
contains

  subroutine Satellite_Merging_Remnant_Progenitor_Properties(satelliteNode,hostNode,satelliteMass,hostMass,satelliteSpheroidMass&
       &,hostSpheroidMass,hostSpheroidMassPreMerger,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass&
       &,remnantSpheroidGasMass)
    !% Calculates progenitor properties for merger remnant calculations.
    use Galacticus_Nodes
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satelliteMergingRemnantProgenitorPropertiesMethod" type="moduleUse">
    include 'satellites.merging.remnant_sizes.progenitor_propeties.modules.inc'
    !# </include>
    implicit none
    type(treeNode),   intent(inout), pointer :: satelliteNode,hostNode
    double precision, intent(out)            :: satelliteMass,hostMass,satelliteSpheroidMass,hostSpheroidMass,hostSpheroidMassPreMerger&
         &,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass,remnantSpheroidGasMass
    
    if (.not.satelliteMergingRemnantProgenitorPropertiesInitialized) then
       !$omp critical(satelliteMergingRemnantProgenitorPropertiesInitialize)
       if (.not.satelliteMergingRemnantProgenitorPropertiesInitialized) then
          ! Get the progenitor properties method parameter.
          !@ <inputParameter>
          !@   <name>satelliteMergingRemnantProgenitorPropertiesMethod</name>
          !@   <defaultValue>standard</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing progenitor properties in merger remnant calculations.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteMergingRemnantProgenitorPropertiesMethod',satelliteMergingRemnantProgenitorPropertiesMethod,defaultValue='standard')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="satelliteMergingRemnantProgenitorPropertiesMethod" type="functionCall" functionType="void">
          !#  <functionArgs>satelliteMergingRemnantProgenitorPropertiesMethod,Satellite_Merging_Remnant_Progenitor_Properties_Get</functionArgs>
          include 'satellites.merging.remnant_sizes.progenitor_properties.inc'
          !# </include>
          if (.not.associated(Satellite_Merging_Remnant_Progenitor_Properties_Get)) call Galacticus_Error_Report('Satellite_Merging_Remnant_Progenitor_Properties','method ' &
               &//char(satelliteMergingRemnantProgenitorPropertiesMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          satelliteMergingRemnantProgenitorPropertiesInitialized=.true.
       end if
       !$omp end critical(satelliteMergingRemnantProgenitorPropertiesInitialize)
    end if

    ! Call the subroutine to perform the calculations.
    call Satellite_Merging_Remnant_Progenitor_Properties_Get(satelliteNode,hostNode,satelliteMass,hostMass,satelliteSpheroidMass&
         &,hostSpheroidMass,hostSpheroidMassPreMerger,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass&
         &,remnantSpheroidGasMass)
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties

end module Satellite_Merging_Remnant_Sizes_Progenitors
