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

!% Contains a module that implements calculations of Gaunt factors for Bremsstrahlung emission from ions.

module Atomic_Radiation_Gaunt_Factors
  !% Implements calculations of Gaunt factors for Bremsstrahlung emission from ions.
  use ISO_Varying_String
  implicit none
  private
  public :: Gaunt_Factor

  ! Flag to indicate if this module has been initialized.
  logical                 :: gauntFactorInitialized=.false.

  ! Name of method used.
  type   (varying_string) :: gauntFactorMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Gaunt_Factor_Template), pointer :: Gaunt_Factor_Get => null()
  abstract interface
     double precision function Gaunt_Factor_Template(atomicNumber,electronNumber,temperature)
       double precision, intent(in  ) :: temperature                     
       integer         , intent(in  ) :: atomicNumber, electronNumber
     end function Gaunt_Factor_Template
  end interface

contains

  subroutine Gaunt_Factor_Initialize()
    !% Initialize the Gaunt factor module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="gauntFactorMethod" type="moduleUse">
    include 'gaunt_factor.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.gauntFactorInitialized) then
       !$omp critical(gauntFactorInitialization)
       if (.not.gauntFactorInitialized) then
          ! Get the method parameter.
          !@ <inputParameter>
          !@   <name>gauntFactorMethod</name>
          !@   <defaultValue>sutherland1998</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing Gaunt factors.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('gauntFactorMethod',gauntFactorMethod,defaultValue='sutherland1998')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="gauntFactorMethod" type="functionCall" functionType="void">
          !#  <functionArgs>gauntFactorMethod,Gaunt_Factor_Get</functionArgs>
          include 'gaunt_factor.inc'
          !# </include>
          if (.not.associated(Gaunt_Factor_Get)) call&
               & Galacticus_Error_Report('Gaunt_Factor_Initialize','method '//char(gauntFactorMethod)//' is unrecognized')
          gauntFactorInitialized = .true.
       end if
       !$omp end critical(gauntFactorInitialization)
    end if
    return
  end subroutine Gaunt_Factor_Initialize

  double precision function Gaunt_Factor(atomicNumber,electronNumber,temperature)
    !% Return the Gaunt factor for the given ion and temperature.
    implicit none
    double precision, intent(in   ) :: temperature                     
    integer         , intent(in   ) :: atomicNumber, electronNumber
    
    ! Initialize the module.
    call Gaunt_Factor_Initialize()
    ! Call the routine to do the calculation.
    Gaunt_Factor=Gaunt_Factor_Get(atomicNumber,electronNumber,temperature)
    return
  end function Gaunt_Factor

end module Atomic_Radiation_Gaunt_Factors
