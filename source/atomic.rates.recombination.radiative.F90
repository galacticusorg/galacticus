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
    
    ! Initialize if necessary.
    if (.not.recombinationRateInitialized) then
       !$omp critical(Atomic_Rate_Recombination_Radiative_Initialization) 
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
          !# <include directive="atomicRadiativeRecombinationMethod" type="functionCall" functionType="void">
          !#  <functionArgs>atomicRadiativeRecombinationMethod,Atomic_Rate_Recombination_Radiative_Get</functionArgs>
          include 'atomic.rates.recombination.radiative.inc'
          !# </include>
          if (.not.associated(Atomic_Rate_Recombination_Radiative_Get)) call&
               & Galacticus_Error_Report('Atomic_Rate_Recombination_Radiative_Initialize','method '//char(atomicRadiativeRecombinationMethod)//' is unrecognized')
          recombinationRateInitialized=.true.
       end if
       !$omp end critical(Atomic_Rate_Recombination_Radiative_Initialization) 
    end if
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
