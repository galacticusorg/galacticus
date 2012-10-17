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

!% Contains a module which provides calculations of Population III supernovae.

module Supernovae_Population_III
  !% Provides calculations of Population III supernovae.
  use ISO_Varying_String
  implicit none
  private
  public :: SNePopIII_Cumulative_Energy

  ! Flag indicating whether this module has been initialized.
  logical              :: supernovaePopIIIInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: supernovaePopIIIMethod

  ! Pointer to the function that actually does the calculation.
  procedure(SNePopIII_Cumulative_Template), pointer :: SNePopIII_Cumulative_Energy_Get => null()
  abstract interface
    double precision function SNePopIII_Cumulative_Template(initialMass,age,metallicity)
      double precision, intent(in) :: initialMass,age,metallicity
    end function SNePopIII_Cumulative_Template
  end interface

contains

  subroutine Supernovae_Population_III_Initialize
    !% Initialize the Population III supernovae module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="supernovaePopIIIMethod" type="moduleUse">
    include 'stellar_astrophysics.supernovae_type_PopIII.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.supernovaePopIIIInitialized) then
       !$omp critical(Supernovae_Population_III_Initialization) 
       if (.not.supernovaePopIIIInitialized) then
          ! Get the halo spin distribution method parameter.
          !@ <inputParameter>
          !@   <name>supernovaePopIIIMethod</name>
          !@   <defaultValue>Heger-Woosley2002</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The method to use for computing properties of Population III supernovae.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('supernovaePopIIIMethod',supernovaePopIIIMethod,defaultValue='Heger-Woosley2002')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="supernovaePopIIIMethod" type="code" action="subroutine">
          !#  <subroutineArgs>supernovaePopIIIMethod,SNePopIII_Cumulative_Energy_Get</subroutineArgs>
          include 'stellar_astrophysics.supernovae_type_PopIII.inc'
          !# </include>
          if (.not.associated(SNePopIII_Cumulative_Energy_Get)) call Galacticus_Error_Report('Supernovae_Population_III_Initialize'&
               &,'method '//char(supernovaePopIIIMethod)//' is unrecognized')
          supernovaePopIIIInitialized=.true.
       end if
       !$omp end critical(Supernovae_Population_III_Initialization) 
    end if
    return
  end subroutine Supernovae_Population_III_Initialize

  double precision function SNePopIII_Cumulative_Energy(initialMass,age,metallicity)
    !% Return the cumulative energy input from Population III supernovae from stars of given {\tt initialMass}, {\tt age} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: initialMass,age,metallicity

    ! Ensure module is initialized.
    call Supernovae_Population_III_Initialize

    ! Simply call the function which does the actual work.
    SNePopIII_Cumulative_Energy=SNePopIII_Cumulative_Energy_Get(initialMass,age,metallicity)
    return
  end function SNePopIII_Cumulative_Energy

end module Supernovae_Population_III
