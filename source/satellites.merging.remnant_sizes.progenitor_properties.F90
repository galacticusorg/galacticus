!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations for progenitor properties for merger remnant calculations.

module Satellite_Merging_Remnant_Sizes_Progenitors
  !% Implements calculations for progenitor properties for merger remnant calculations.
  use ISO_Varying_String
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
    use Tree_Nodes
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satelliteMergingRemnantProgenitorPropertiesMethod" type="moduleUse">
    include 'satellites.merging.remnant_sizes.progenitor_propeties.modules.inc'
    !# </include>
    implicit none
    type(treeNode),   intent(inout), pointer :: satelliteNode,hostNode
    double precision, intent(out)            :: satelliteMass,hostMass,satelliteSpheroidMass,hostSpheroidMass,hostSpheroidMassPreMerger&
         &,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass,remnantSpheroidGasMass
    
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
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteMergingRemnantProgenitorPropertiesMethod',satelliteMergingRemnantProgenitorPropertiesMethod,defaultValue='standard')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="satelliteMergingRemnantProgenitorPropertiesMethod" type="code" action="subroutine">
       !#  <subroutineArgs>satelliteMergingRemnantProgenitorPropertiesMethod,Satellite_Merging_Remnant_Progenitor_Properties_Get</subroutineArgs>
       include 'satellites.merging.remnant_sizes.progenitor_properties.inc'
       !# </include>
       if (.not.associated(Satellite_Merging_Remnant_Progenitor_Properties_Get)) call Galacticus_Error_Report('Satellite_Merging_Remnant_Progenitor_Properties','method ' &
            &//char(satelliteMergingRemnantProgenitorPropertiesMethod)//' is unrecognized')
       ! Flag that the module is now initialized.
       satelliteMergingRemnantProgenitorPropertiesInitialized=.true.
    end if
    !$omp end critical(satelliteMergingRemnantProgenitorPropertiesInitialize)

    ! Call the subroutine to perform the calculations.
    call Satellite_Merging_Remnant_Progenitor_Properties_Get(satelliteNode,hostNode,satelliteMass,hostMass,satelliteSpheroidMass&
         &,hostSpheroidMass,hostSpheroidMassPreMerger,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass&
         &,remnantSpheroidGasMass)
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties

end module Satellite_Merging_Remnant_Sizes_Progenitors
