!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that implements calculations of atomic collisional ionization rates.

module Atomic_Rates_Ionization_Collisional
  !% Implements calculations of atomic collisional ionization rates.
  use ISO_Varying_String
  implicit none
  private
  public :: Atomic_Rate_Ionization_Collisional

  ! Flag to indicate if this module has been initialized.
  logical                                                         :: ionizationRateInitialized             =.false.

  ! Name of ionization state method used.
  type     (varying_string                             )          :: atomicCollisionalIonizationMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Atomic_Rate_Ionization_Collisional_Template), pointer :: Atomic_Rate_Ionization_Collisional_Get=>null()
  abstract interface
     double precision function Atomic_Rate_Ionization_Collisional_Template(atomicNumber,ionizationState,temperature)
       integer         , intent(in   ) :: atomicNumber, ionizationState
       double precision, intent(in   ) :: temperature
     end function Atomic_Rate_Ionization_Collisional_Template
  end interface

contains

  subroutine Atomic_Rate_Ionization_Collisional_Initialize
    !% Initialize the atomic collisional ionization rate module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="atomicCollisionalIonizationMethod" type="moduleUse">
    include 'atomic.rates.ionization.collisional.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.ionizationRateInitialized) then
       !$omp critical(Atomic_Rate_Ionization_Collisional_Initialization)
       if (.not.ionizationRateInitialized) then
          ! Get the ionization state method parameter.
          !@ <inputParameter>
          !@   <name>atomicCollisionalIonizationMethod</name>
          !@   <defaultValue>Verner</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing atomic collisional ionization rates.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('atomicCollisionalIonizationMethod',atomicCollisionalIonizationMethod,defaultValue='Verner')

          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="atomicCollisionalIonizationMethod" type="functionCall" functionType="void">
          !#  <functionArgs>atomicCollisionalIonizationMethod,Atomic_Rate_Ionization_Collisional_Get</functionArgs>
          include 'atomic.rates.ionization.collisional.inc'
          !# </include>
          if (.not.associated(Atomic_Rate_Ionization_Collisional_Get)) call&
               & Galacticus_Error_Report('Atomic_Rate_Ionization_Collisional_Initialize','method '//char(atomicCollisionalIonizationMethod)//' is unrecognized')
          ionizationRateInitialized=.true.
       end if
       !$omp end critical(Atomic_Rate_Ionization_Collisional_Initialization)
    end if
    return
  end subroutine Atomic_Rate_Ionization_Collisional_Initialize

  double precision function Atomic_Rate_Ionization_Collisional(atomicNumber,ionizationState,temperature)
    implicit none
    integer         , intent(in   ) :: atomicNumber, ionizationState
    double precision, intent(in   ) :: temperature

    ! Initialize the module.
    call Atomic_Rate_Ionization_Collisional_Initialize

    ! Call the routine to do the calculation.
    Atomic_Rate_Ionization_Collisional=Atomic_Rate_Ionization_Collisional_Get(atomicNumber,ionizationState,temperature)

    return
  end function Atomic_Rate_Ionization_Collisional

end module Atomic_Rates_Ionization_Collisional
