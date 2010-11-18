!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which provides support for stellar population properties.

module Stellar_Population_Properties
  !% Provides support for stellar population properties.
  use ISO_Varying_String
  use Tree_Nodes
  use Abundances_Structure
  use Histories
  private
  public :: Stellar_Population_Properties_Rates, Stellar_Population_Properties_Scales,&
       & Stellar_Population_Properties_History_Count, Stellar_Population_Properties_History_Create

  ! Flag indicating whether this module has been initialized.
  logical              :: stellarPopulationPropertiesInitialized=.false.

  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationTimescaleDisksInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: stellarPopulationPropertiesMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Stellar_Population_Properties_Rates_Template), pointer :: Stellar_Population_Properties_Rates_Get => null()
  abstract interface
     subroutine Stellar_Population_Properties_Rates_Template(starFormationRate,fuelAbundances,thisNode,thisHistory &
          &,stellarMassRate,stellarAbundancesRates,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
       import treeNode, abundancesStructure, history
       double precision,          intent(out)                 :: stellarMassRate,fuelMassRate,energyInputRate
       type(abundancesStructure), intent(inout)               :: stellarAbundancesRates,fuelAbundancesRates
       double precision,          intent(out),   dimension(:) :: stellarLuminositiesRates
       double precision,          intent(in)                  :: starFormationRate
       type(abundancesStructure), intent(in)                  :: fuelAbundances
       type(treeNode),            intent(inout), pointer      :: thisNode
       type(history),             intent(inout)               :: thisHistory
     end subroutine Stellar_Population_Properties_Rates_Template
  end interface

  ! Pointer to the function that sets scale factors for error control of stellar population properties.
  procedure(Stellar_Population_Properties_Scales_Template), pointer :: Stellar_Population_Properties_Scales_Get => null()
  abstract interface
     subroutine Stellar_Population_Properties_Scales_Template(thisHistory,stellarMass,stellarAbundances)
       import abundancesStructure, history
       double precision,          intent(in)                  :: stellarMass
       type(abundancesStructure), intent(in)                  :: stellarAbundances
       type(history),             intent(inout)               :: thisHistory
     end subroutine Stellar_Population_Properties_Scales_Template
  end interface

  ! Pointer to the function that returns the size of any history required for stellar population properties.
  procedure(Stellar_Population_Properties_History_Count_Template), pointer :: Stellar_Population_Properties_History_Count_Get => null()
  abstract interface
     integer function Stellar_Population_Properties_History_Count_Template()
     end function Stellar_Population_Properties_History_Count_Template
  end interface

  ! Pointer to the subroutine that creates any history required for stellar population properties.
  procedure(Stellar_Population_Properties_History_Create_Template), pointer :: Stellar_Population_Properties_History_Create_Do => null()
  abstract interface
     subroutine Stellar_Population_Properties_History_Create_Template(thisNode,thisHistory)
       import treeNode, history
       type(treeNode), intent(inout), pointer :: thisNode
       type(history),  intent(inout)          :: thisHistory
     end subroutine Stellar_Population_Properties_History_Create_Template
  end interface

contains

  subroutine Stellar_Population_Properties_Rates_Initialize
    !% Initialize the disk star formation timecale module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="stellarPopulationPropertiesMethod" type="moduleUse">
    include 'stellar_populations.properties.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Stellar_Population_Properties_Rates_Initialization) 
    ! Initialize if necessary.
    if (.not.stellarPopulationPropertiesInitialized) then
       ! Get the halo spin distribution method parameter.
       !@ <inputParameter>
       !@   <name>stellarPopulationPropertiesMethod</name>
       !@   <defaultValue>instantaneous</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The method to use for computing properties of stellar populations.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationPropertiesMethod',stellarPopulationPropertiesMethod,defaultValue='instantaneous')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="stellarPopulationPropertiesMethod" type="code" action="subroutine">
       !#  <subroutineArgs>stellarPopulationPropertiesMethod,Stellar_Population_Properties_Rates_Get,Stellar_Population_Properties_Scales_Get,Stellar_Population_Properties_History_Count_Get,Stellar_Population_Properties_History_Create_Do</subroutineArgs>
       include 'stellar_populations.properties.inc'
       !# </include>
       if (.not.(associated(Stellar_Population_Properties_Rates_Get).and.associated(Stellar_Population_Properties_Scales_Get).and.associated(Stellar_Population_Properties_History_Count_Get).and.associated(Stellar_Population_Properties_History_Create_Do))) call Galacticus_Error_Report('Stellar_Population_Properties_Rates'&
            &,'method '//char(stellarPopulationPropertiesMethod)//' is unrecognized')
       stellarPopulationPropertiesInitialized=.true.
    end if
    !$omp end critical(Stellar_Population_Properties_Rates_Initialization) 

    return
  end subroutine Stellar_Population_Properties_Rates_Initialize

  subroutine Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,thisNode,thisHistory,stellarMassRate&
       &,stellarAbundancesRates ,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
    !% Return an array of stellar population property rates of change given a star formation rate and fuel abundances.
    implicit none
    double precision,          intent(out)                 :: stellarMassRate,fuelMassRate,energyInputRate
    type(abundancesStructure), intent(inout)               :: stellarAbundancesRates,fuelAbundancesRates
    double precision,          intent(out),   dimension(:) :: stellarLuminositiesRates
    double precision,          intent(in)                  :: starFormationRate
    type(abundancesStructure), intent(in)                  :: fuelAbundances
    type(treeNode),            intent(inout), pointer      :: thisNode
    type(history),             intent(inout)               :: thisHistory
    
    ! Ensure module is initialized.
    call Stellar_Population_Properties_Rates_Initialize

    ! Simply call the subroutine which does the actual work.
    call Stellar_Population_Properties_Rates_Get(starFormationRate,fuelAbundances,thisNode,thisHistory,stellarMassRate&
         &,stellarAbundancesRates,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
    return
  end subroutine Stellar_Population_Properties_Rates

  subroutine Stellar_Population_Properties_Scales(thisHistory,stellarMass,stellarAbundances)
    !% Set the scaling factors for error control on the absolute value of stellar population properties.
    implicit none
    double precision,          intent(in)    :: stellarMass
    type(abundancesStructure), intent(in)    :: stellarAbundances
    type(history),             intent(inout) :: thisHistory
    
    ! Ensure module is initialized.
    call Stellar_Population_Properties_Rates_Initialize

    ! Simply call the subroutine which does the actual work.
    call Stellar_Population_Properties_Scales_Get(thisHistory,stellarMass,stellarAbundances)
    return
  end subroutine Stellar_Population_Properties_Scales

  integer function Stellar_Population_Properties_History_Count()
    !% Return a count of the number of histories which must be stored for the selected stellar populations method.
    implicit none
  
    ! Ensure module is initialized.
    call Stellar_Population_Properties_Rates_Initialize

    ! Simply call the function which does the actual work.
    Stellar_Population_Properties_History_Count=Stellar_Population_Properties_History_Count_Get()

    return
  end function Stellar_Population_Properties_History_Count

  subroutine Stellar_Population_Properties_History_Create(thisNode,thisHistory)
    !% Create any history required for storing stellar population properties.
    use Histories
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(history),  intent(inout)          :: thisHistory
  
    ! Ensure module is initialized.
    call Stellar_Population_Properties_Rates_Initialize

    ! Simply call the function which does the actual work.
    call Stellar_Population_Properties_History_Create_Do(thisNode,thisHistory)

    return
  end subroutine Stellar_Population_Properties_History_Create

end module Stellar_Population_Properties
