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

!% Contains a module which provides calculations of stellar winds.

module Stellar_Astrophysics_Winds
  !% Provides calculations of stellar winds.
  use ISO_Varying_String
  implicit none
  private
  public :: Stellar_Winds_Mass_Loss_Rate, Stellar_Winds_Terminal_Velocity

  ! Flag indicating whether this module has been initialized.
  logical              :: stellarWindsInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: stellarWindsMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Stellar_Winds_Get_Template), pointer :: Stellar_Winds_Mass_Loss_Rate_Get    => null()
  procedure(Stellar_Winds_Get_Template), pointer :: Stellar_Winds_Terminal_Velocity_Get => null()
  abstract interface
    double precision function Stellar_Winds_Get_Template(initialMass,age,metallicity)
      double precision, intent(in) :: initialMass,age,metallicity
    end function Stellar_Winds_Get_Template
  end interface

contains

  subroutine Stellar_Winds_Initialize
    !% Initialize the stellar winds module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="stellarWindsMethod" type="moduleUse">
    include 'stellar_astrophysics.winds.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.stellarWindsInitialized) then
       !$omp critical(Stellar_Winds_Initialization) 
       if (.not.stellarWindsInitialized) then
          ! Get the halo spin distribution method parameter.
          !@ <inputParameter>
          !@   <name>stellarWindsMethod</name>
          !@   <defaultValue>standard</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The method to use for computing aspects of stellar winds.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('stellarWindsMethod',stellarWindsMethod,defaultValue='Leitherer1992')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="stellarWindsMethod" type="functionCall" functionType="void">
          !#  <functionArgs>stellarWindsMethod,Stellar_Winds_Mass_Loss_Rate_Get,Stellar_Winds_Terminal_Velocity_Get</functionArgs>
          include 'stellar_astrophysics.winds.inc'
          !# </include>
          if (.not.(associated(Stellar_Winds_Mass_Loss_Rate_Get).or.associated(Stellar_Winds_Terminal_Velocity_Get))) call&
               & Galacticus_Error_Report('Stellar_Winds_Initialize','method '//char(stellarWindsMethod)//' is unrecognized')
          stellarWindsInitialized=.true.
       end if
       !$omp end critical(Stellar_Winds_Initialization) 
    end if
    return
  end subroutine Stellar_Winds_Initialize

  double precision function Stellar_Winds_Mass_Loss_Rate(initialMass,age,metallicity)
    !% Return the mass loss rate (in $M_\odot$/Gyr) from stars of given {\tt initialMass}, {\tt age} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: initialMass,age,metallicity

    ! Ensure module is initialized.
    call Stellar_Winds_Initialize

    ! Simply call the function which does the actual work.
    Stellar_Winds_Mass_Loss_Rate=Stellar_Winds_Mass_Loss_Rate_Get(initialMass,age,metallicity)
    return
  end function Stellar_Winds_Mass_Loss_Rate

  double precision function Stellar_Winds_Terminal_Velocity(initialMass,age,metallicity)
    !% Return the terminal velocity (in km/s) of winds from stars of given {\tt initialMass}, {\tt age} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: initialMass,age,metallicity

    ! Ensure module is initialized.
    call Stellar_Winds_Initialize

    ! Simply call the function which does the actual work.
    Stellar_Winds_Terminal_Velocity=Stellar_Winds_Terminal_Velocity_Get(initialMass,age,metallicity)
    return
  end function Stellar_Winds_Terminal_Velocity

end module Stellar_Astrophysics_Winds
