!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module that implements calculations of collisional excitation rates.

module Atomic_Rates_Excitation_Collisional
  !% Implements calculations of collisional excitation rates.
  use ISO_Varying_String
  implicit none
  private
  public :: Collisional_Excitation_Cooling_Rate

  ! Flag to indicate if this module has been initialized.
  logical                 :: collisionalExcitationInitialized=.false.

  ! Name of collisional excitation method used.
  type   (varying_string) :: collisionalExcitationMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Collisional_Excitation_Rate_Template), pointer :: Collisional_Excitation_Rate_Get => null()
  abstract interface
     double precision function Collisional_Excitation_Rate_Template(atomicNumber,electronNumber,temperature)
       double precision,   intent(in) :: temperature                     
       integer,            intent(in) :: atomicNumber, electronNumber
     end function Collisional_Excitation_Rate_Template
   end interface

contains

  subroutine Collisional_Excitation_Rate_Initialize
    !% Initialize the collisional excitation rate module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="collisionalExcitationMethod" type="moduleUse">
    include 'atomic.rates.excitation.collisional.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.collisionalExcitationInitialized) then
       !$omp critical(Atomic_Rate_Excitation_Collisional_Initialization)
       if (.not.collisionalExcitationInitialized) then
          ! Get the collisional excitation method parameter.
          !@ <inputParameter>
          !@   <name>collisionalExcitationMethod</name>
          !@   <defaultValue>ScholzWalters91</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing collisional excitation rates.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('collisionalExcitationMethod',collisionalExcitationMethod,defaultValue='ScholzWalters91')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="collisionalExcitationMethod" type="functionCall" functionType="void">
          !#  <functionArgs>collisionalExcitationMethod,Collisional_Excitation_Rate_Get</functionArgs>
          include 'atomic.rates.excitation.collisional.inc'
          !# </include>
          if (.not.associated(Collisional_Excitation_Rate_Get)) call&
               & Galacticus_Error_Report('Collisional_Excitation_Rate_Initialize','method '//char(collisionalExcitationMethod)//' is unrecognized')
          collisionalExcitationInitialized = .true.
       end if
       !$omp end critical(Atomic_Rate_Excitation_Collisional_Initialization)
    end if
    return
  end subroutine Collisional_Excitation_Rate_Initialize

  double precision function Collisional_Excitation_Cooling_Rate(atomicNumber, electronNumber, temperature)
    !% Return the collisional excitation cooling rate , in units of J/m$^3$/s, for ion of given {\normalfont \ttfamily atomicNumber} and {\tt
    !% electronNumber} at temperature {\normalfont \ttfamily T} (in Kelvin).
    implicit none
    double precision, intent(in   ) :: temperature                     
    integer         , intent(in   ) :: atomicNumber, electronNumber

    ! Initialize the module.
    call Collisional_Excitation_Rate_Initialize()
    ! Call the routine to do the calculation.
    Collisional_Excitation_Cooling_Rate=Collisional_Excitation_Rate_Get(atomicNumber,electronNumber,temperature)
    return
  end function Collisional_Excitation_Cooling_Rate

end module Atomic_Rates_Excitation_Collisional
