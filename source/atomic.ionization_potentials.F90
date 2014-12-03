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

!+ Contributions to this file made by: Daniel McAndrew.

!% Contains a module that implements calculations of atomic ionization potentials.

module Atomic_Ionization_Potentials
  !% Implements calculations of atomic photo-ionization cross-sections.
  use ISO_Varying_String 
  implicit none
  private
  public :: Atomic_Ionization_Potential

  ! Flag to indicate if this module has been initialized.  
  logical                                         :: moduleInitialized               =  .false.

  ! Name of ionization state method used.
  type     (varying_string             )          :: atomicIonizationPotentialMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Atomic_Ionization_Potential), pointer :: Atomic_Ionization_Potential_Get => null()
  
contains

  subroutine Atomic_Ionization_Potential_Initialize()
    !% Initialize the atomic ionization potential module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="atomicIonizationPotentialMethod" type="moduleUse">
    include 'atomic.ionization_potentials.modules.inc'
    !# </include>
    implicit none
    
    ! Initialize if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Atomic_Ionization_Potential_Initialize)  
       if (.not.moduleInitialized) then
          ! Get the ionization potential method parameter.
          !@ <inputParameter>
          !@   <name>atomicIonizationPotentialMethod</name>
          !@   <defaultValue>Verner</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing atomic ionization potentials.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('atomicIonizationPotentialMethod',atomicIonizationPotentialMethod,defaultValue='Verner')
                    ! Include file that makes calls to all available method initialization routines.
          !# <include directive="atomicIonizationPotentialMethod" type="functionCall" functionType="void">
          !#  <functionArgs>atomicIonizationPotentialMethod,Atomic_Ionization_Potential_Get</functionArgs>
          include 'atomic.ionization_potentials.inc'
          !# </include>
          if (.not.associated(Atomic_Ionization_Potential_Get)) call&
               & Galacticus_Error_Report('Atomic_Ionization_Potential_Initialize','method '//char(atomicIonizationPotentialMethod)//' is unrecognized')
          moduleInitialized=.true.
       end if
       !$omp end critical(Atomic_Ionization_Potential_Initialize) 
    end if
    return
  end subroutine Atomic_Ionization_Potential_Initialize

  double precision function Atomic_Ionization_Potential(atomicNumber,electronNumber)
    !% Return the ionization potential (in units of eV) for a given atom in a given ionization state.
    implicit none
    integer, intent(in   ) :: atomicNumber,electronNumber

    ! Initialize the module.
    call Atomic_Ionization_Potential_Initialize()
    ! Call the routine to do the calculation.
    Atomic_Ionization_Potential=Atomic_Ionization_Potential_Get(atomicNumber,electronNumber)
    return
  end function Atomic_Ionization_Potential
  
end module Atomic_Ionization_Potentials
