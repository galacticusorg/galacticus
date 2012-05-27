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


!% Contains a module that implements calculations of chemical reaction rates.

module Chemical_Reaction_Rates
  !% Implements calculations of chemical reaction rates.
  use Chemical_Abundances_Structure
  use ISO_Varying_String 
  implicit none
  private
  public :: Chemical_Reaction_Rate

  ! Flag to indicate if this module has been initialized.  
  logical                                         :: chemicalReactionRateInitialized=.false.

  ! Name of chemical reaction rates methods used.
  type(varying_string), allocatable, dimension(:) :: chemicalReactionRateMethods

contains

  subroutine Chemical_Reaction_Rates_Initialize
    !% Initialize the chemical reaction rates module.
    use Galacticus_Error
    use Input_Parameters
    use Memory_Management
    !# <include directive="chemicalReactionRates" type="moduleUse">
    include 'chemical.reaction_rates.modules.inc'
    !# </include>
    implicit none
    integer :: chemicalReactionRatesCount
    
    ! Initialize if necessary.
    if (.not.chemicalReactionRateInitialized) then
       !$omp critical(Chemical_Reaction_Rates_Initialization) 
       if (.not.chemicalReactionRateInitialized) then
          ! Get the chemical reaction rates method parameter.
          !@ <inputParameter>
          !@   <name>chemicalReactionRateMethods</name>
          !@   <defaultValue>hydrogenNetwork</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The names of the methods to be used for computing chemical reaction rates.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1..*</cardinality>
          !@ </inputParameter>
          chemicalReactionRatesCount=max(1,Get_Input_Parameter_Array_Size('chemicalReactionRatesMethods'))
          allocate(chemicalReactionRateMethods(chemicalReactionRatesCount))
          call Memory_Usage_Record(sizeof(chemicalReactionRateMethods))
          call Get_Input_Parameter('chemicalReactionRateMethods',chemicalReactionRateMethods,defaultValue=['null'])
          
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="chemicalReactionRates" type="code" action="subroutine">
          !#  <subroutineArgs>chemicalReactionRateMethods</subroutineArgs>
          include 'chemical.reaction_rates.inc'
          !# </include>
          chemicalReactionRateInitialized=.true.
       end if
       !$omp end critical(Chemical_Reaction_Rates_Initialization) 
    end if
    return
  end subroutine Chemical_Reaction_Rates_Initialize

  subroutine Chemical_Reaction_Rate(chemicalRates,temperature,chemicalDensity,radiation)
    !% Return chemical reaction rates at the given temperature for the specified set of chemical densities (in cm$^{-3}$) and radiation
    !% field. Units of the returned rates are cm$^-3$ s$^{-1}$.
    use Abundances_Structure
    use Radiation_Structure
    !# <include directive="chemicalRatesCompute" type="moduleUse">
    include 'chemical.reaction_rates.compute.modules.inc'
    !# </include>
    implicit none
    type(chemicalAbundancesStructure), intent(inout) :: chemicalRates
    double precision,                  intent(in)    :: temperature
    type(chemicalAbundancesStructure), intent(in)    :: chemicalDensity
    type(radiationStructure),          intent(in)    :: radiation

    ! Initialize the module.
    call Chemical_Reaction_Rates_Initialize
  
    call chemicalRates%reset()
    !# <include directive="chemicalRatesCompute" type="code" action="subroutine">
    !#  <subroutineArgs>temperature,chemicalDensity,radiation,chemicalRates</subroutineArgs>
    include 'chemical.reaction_rates.compute.inc'
    !# </include>

    return
  end subroutine Chemical_Reaction_Rate
  
end module Chemical_Reaction_Rates
