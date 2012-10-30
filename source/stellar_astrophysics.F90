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

!% Contains a module which implements calculation of stellar astrophysics.

module Stellar_Astrophysics
  !% Implements calculation of stellar astrophysics.
  use ISO_Varying_String
  implicit none
  private
  public :: Star_Ejected_Mass, Star_Initial_Mass, Star_Metal_Yield_Mass, Star_Lifetime
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: stellarAstrophysicsInitialized=.false.

  ! Name of cosmology functions method used.
  type(varying_string) :: stellarAstrophysicsMethod

  ! Pointer to the functions that actually do the calculations.
  procedure(Stellar_Astrophysics_Template),       pointer :: Star_Ejected_Mass_Get     => null()
  procedure(Stellar_Astrophysics_Template),       pointer :: Star_Initial_Mass_Get     => null()
  procedure(Stellar_Astrophysics_Yield_Template), pointer :: Star_Metal_Yield_Mass_Get => null()
  procedure(Stellar_Astrophysics_Template),       pointer :: Star_Lifetime_Get         => null()

  abstract interface
     double precision function Stellar_Astrophysics_Template(inputParameter1,inputParameter2)
       double precision, intent(in) :: inputParameter1,inputParameter2
     end function Stellar_Astrophysics_Template
  end interface

  abstract interface
     double precision function Stellar_Astrophysics_Yield_Template(inputParameter1,inputParameter2,atomIndex)
       double precision, intent(in)           :: inputParameter1,inputParameter2
       integer,          intent(in), optional :: atomIndex
     end function Stellar_Astrophysics_Yield_Template
  end interface

contains

  subroutine Stellar_Astrophysics_Initialize
    !% Initialize the stellar astrophysics module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="stellarAstrophysicsMethod" type="moduleUse">
    include 'stellar_astrophysics.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.stellarAstrophysicsInitialized) then
       !$omp critical(Stellar_Astrophysics_Initialization) 
       if (.not.stellarAstrophysicsInitialized) then
          ! Get the stellar tracks method parameter.
          !@ <inputParameter>
          !@   <name>stellarAstrophysicsMethod</name>
          !@   <defaultValue>file</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for stellar astrophysics calculations.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('stellarAstrophysicsMethod',stellarAstrophysicsMethod,defaultValue='file')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="stellarAstrophysicsMethod" type="functionCall" functionType="void">
          !#  <functionArgs>stellarAstrophysicsMethod,Star_Ejected_Mass_Get,Star_Initial_Mass_Get,Star_Metal_Yield_Mass_Get,Star_Lifetime_Get</functionArgs>
          include 'stellar_astrophysics.inc'
          !# </include>
          if (.not.(associated(Star_Ejected_Mass_Get).and.associated(Star_Initial_Mass_Get).and.associated(Star_Metal_Yield_Mass_Get).and.associated(Star_Lifetime_Get))) &
               & call Galacticus_Error_Report('Stellar_Astrophysics','method '//char(stellarAstrophysicsMethod)//' is unrecognized')
          stellarAstrophysicsInitialized=.true.
       end if
       !$omp end critical(Stellar_Astrophysics_Initialization) 
    end if
    return
  end subroutine Stellar_Astrophysics_Initialize

  double precision function Star_Initial_Mass(lifetime,metallicity)
    !% Returns the initial mass of a star of given {\tt lifetime} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: lifetime,metallicity
    
    ! Initialize the module.
    call Stellar_Astrophysics_Initialize
    
    ! Get the answer using the selected method.
    Star_Initial_Mass=Star_Initial_Mass_Get(lifetime,metallicity)
    
    return
  end function Star_Initial_Mass
  
  double precision function Star_Ejected_Mass(initialMass,metallicity)
    !% Returns the mass ejected by a star of given {\tt initialMass} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: initialMass,metallicity
    
    ! Initialize the module.
    call Stellar_Astrophysics_Initialize
    
    ! Get the answer using the selected method.
    Star_Ejected_Mass=Star_Ejected_Mass_Get(initialMass,metallicity)
    
    return
  end function Star_Ejected_Mass
  
  double precision function Star_Metal_Yield_Mass(initialMass,metallicity,atomIndex)
    !% Returns the metal mass yielded by a star of given {\tt initialMass} and {\tt metallicity}.
    implicit none
    double precision, intent(in)           :: initialMass,metallicity
    integer,          intent(in), optional :: atomIndex
    
    ! Initialize the module.
    call Stellar_Astrophysics_Initialize
    
    ! Get the answer using the selected method.
    Star_Metal_Yield_Mass=Star_Metal_Yield_Mass_Get(initialMass,metallicity,atomIndex)
    
    return
  end function Star_Metal_Yield_Mass
  
  double precision function Star_Lifetime(initialMass,metallicity)
    !% Returns the lifetime of a star of given {\tt initialMass} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: initialMass,metallicity
    
    ! Initialize the module.
    call Stellar_Astrophysics_Initialize
    
    ! Get the answer using the selected method.
    Star_Lifetime=Star_Lifetime_Get(initialMass,metallicity)
    
    return
  end function Star_Lifetime
  
end module Stellar_Astrophysics
