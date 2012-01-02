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


!% Contains a module that implements calculations of atomic radiative recombination rates.

module Atomic_Rates_Recombination_Radiative
  !% Implements calculations of atomic radiative recombination rates.
  use ISO_Varying_String 
  implicit none
  private
  public :: Atomic_Rate_Recombination_Radiative

  ! Flag to indicate if this module has been initialized.  
  logical              :: recombinationRateInitialized=.false.

  ! Name of ionization state method used.
  type(varying_string) :: atomicRadiativeRecombinationMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Atomic_Rate_Recombination_Radiative_Template), pointer :: Atomic_Rate_Recombination_Radiative_Get => null()
  abstract interface
     double precision function Atomic_Rate_Recombination_Radiative_Template(atomicNumber,ionizationState,temperature)
       integer,          intent(in) :: atomicNumber,ionizationState
       double precision, intent(in) :: temperature
     end function Atomic_Rate_Recombination_Radiative_Template
  end interface
  
contains

  subroutine Atomic_Rate_Recombination_Radiative_Initialize
    !% Initialize the atomic radiative recombination rate module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="atomicRadiativeRecombinationMethod" type="moduleUse">
    include 'atomic.rates.recombination.radiative.modules.inc'
    !# </include>
    implicit none
    
    !$omp critical(Atomic_Rate_Recombination_Radiative_Initialization) 
    ! Initialize if necessary.
    if (.not.recombinationRateInitialized) then
       ! Get the ionization state method parameter.
       !@ <inputParameter>
       !@   <name>atomicRadiativeRecombinationMethod</name>
       !@   <defaultValue>Verner</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing atomic radiative recombination rates.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('atomicRadiativeRecombinationMethod',atomicRadiativeRecombinationMethod,defaultValue='Verner')

       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="atomicRadiativeRecombinationMethod" type="code" action="subroutine">
       !#  <subroutineArgs>atomicRadiativeRecombinationMethod,Atomic_Rate_Recombination_Radiative_Get</subroutineArgs>
       include 'atomic.rates.recombination.radiative.inc'
       !# </include>
       if (.not.associated(Atomic_Rate_Recombination_Radiative_Get)) call&
            & Galacticus_Error_Report('Atomic_Rate_Recombination_Radiative_Initialize','method '//char(atomicRadiativeRecombinationMethod)//' is unrecognized')
       recombinationRateInitialized=.true.
    end if
    !$omp end critical(Atomic_Rate_Recombination_Radiative_Initialization) 
    return
  end subroutine Atomic_Rate_Recombination_Radiative_Initialize

  double precision function Atomic_Rate_Recombination_Radiative(atomicNumber,ionizationState,temperature)
    implicit none
    integer,          intent(in) :: atomicNumber,ionizationState
    double precision, intent(in) :: temperature

    ! Initialize the module.
    call Atomic_Rate_Recombination_Radiative_Initialize

    ! Call the routine to do the calculation.
    Atomic_Rate_Recombination_Radiative=Atomic_Rate_Recombination_Radiative_Get(atomicNumber,ionizationState,temperature)
    
    return
  end function Atomic_Rate_Recombination_Radiative

end module Atomic_Rates_Recombination_Radiative
