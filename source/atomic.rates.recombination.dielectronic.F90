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

!% Contains a module that implements calculations of dielectronic recombination rates.

module Atomic_Rates_Recombination_Dielectronic
  !% Implements calculations of dielectronic recombination rates.
  use ISO_Varying_String
  implicit none
  private
  public :: Dielectronic_Recombination_Rate

  ! Flag to indicate if this module has been initialized.
  logical                 :: dielectronicRateInitialized    =.false.

  ! Name of ionization state method used.
  type   (varying_string) :: dielectronicRecombinationMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Dielectronic_Recombination_Ratebination_Rate_Template), pointer :: Dielectronic_Recombination_Ratebination_Rate_Get => null()
  abstract interface
     double precision function Dielectronic_Recombination_Ratebination_Rate_Template(atomicNumber, electronNumber, temperature)
       double precision, intent(in   ) :: temperature                     
       integer         , intent(in   ) :: atomicNumber, electronNumber
     end function Dielectronic_Recombination_Ratebination_Rate_Template
   end interface

contains

  subroutine Dielectronic_Recombination_Ratebination_Rate_Initialize
    !% Initialize the dielectronic recombination rate module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="dielectronicRecombinationMethod" type="moduleUse">
    include 'atomic.rates.recombination.dielectronic.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.dielectronicRateInitialized) then
       !$omp critical(Atomic_Rate_Recombination_Dielectronic_Initialization)
       if (.not.dielectronicRateInitialized) then
          ! Get the ionization state method parameter.
          !@ <inputParameter>
          !@   <name>dielectronicRecombinationMethod</name>
          !@   <defaultValue>Arnaud85</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing dielectronic recombination rates.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('dielectronicRecombinationMethod',dielectronicRecombinationMethod,defaultValue='Arnaud85')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="dielectronicRecombinationMethod" type="functionCall" functionType="void">
          !#  <functionArgs>dielectronicRecombinationMethod,Dielectronic_Recombination_Ratebination_Rate_Get</functionArgs>
          include 'atomic.rates.recombination.dielectronic.inc'
          !# </include>
          if (.not.associated(Dielectronic_Recombination_Ratebination_Rate_Get)) call&
               & Galacticus_Error_Report('Dielectronic_Recombination_Ratebination_Rate_Initialize','method '//char(dielectronicRecombinationMethod)//' is unrecognized')
          dielectronicRateInitialized = .true.
       end if
       !$omp end critical(Atomic_Rate_Recombination_Dielectronic_Initialization)
    end if
    return
  end subroutine Dielectronic_Recombination_Ratebination_Rate_Initialize

  double precision function Dielectronic_Recombination_Rate(atomicNumber, electronNumber, temperature)
    !% Return the dielectroninc recombination rate (in units of cm$^3$ s$^{-1}$) for the ion of given {\normalfont \ttfamily atomicNumber} and {\tt
    !% electronNumber} at the given {\normalfont \ttfamily temperature} (in Kelvin).
    implicit none
    double precision, intent(in   ) :: temperature                     
    integer         , intent(in   ) :: atomicNumber, electronNumber

    ! Initialize the module.
    call Dielectronic_Recombination_Ratebination_Rate_Initialize()
    ! Call the routine to do the calculation.
    Dielectronic_Recombination_Rate=Dielectronic_Recombination_Ratebination_Rate_Get(atomicNumber,electronNumber,temperature)
    return
  end function Dielectronic_Recombination_Rate

end module Atomic_Rates_Recombination_Dielectronic
