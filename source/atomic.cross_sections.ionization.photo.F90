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

!% Contains a module that implements calculations of atomic photo-ionization cross-sections.

module Atomic_Cross_Sections_Ionization_Photo
  !% Implements calculations of atomic photo-ionization cross-sections.
  use ISO_Varying_String 
  implicit none
  private
  public :: Atomic_Cross_Section_Ionization_Photo

  ! Flag to indicate if this module has been initialized.  
  logical              :: ionizationCrossSectionInitialized=.false.

  ! Name of ionization state method used.
  type(varying_string) :: atomicPhotoIonizationMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Atomic_Cross_Section_Ionization_Photo_Template), pointer :: Atomic_Cross_Section_Ionization_Photo_Get => null()
  abstract interface
     double precision function Atomic_Cross_Section_Ionization_Photo_Template(atomicNumber,ionizationState,shellNumber,wavelength)
       integer,          intent(in) :: atomicNumber,ionizationState,shellNumber
       double precision, intent(in) :: wavelength
     end function Atomic_Cross_Section_Ionization_Photo_Template
  end interface
  
contains

  subroutine Atomic_Cross_Section_Ionization_Photo_Initialize
    !% Initialize the atomic photo ionization cross-section module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="atomicPhotoIonizationMethod" type="moduleUse">
    include 'atomic.cross_sections.ionization.photo.modules.inc'
    !# </include>
    implicit none
    
    ! Initialize if necessary.
    if (.not.ionizationCrossSectionInitialized) then
       !$omp critical(Atomic_Cross_Section_Ionization_Photo_Initialization)  
       if (.not.ionizationCrossSectionInitialized) then
          ! Get the ionization state method parameter.
          !@ <inputParameter>
          !@   <name>atomicPhotoIonizationMethod</name>
          !@   <defaultValue>Verner</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing atomic photo ionization rates.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('atomicPhotoIonizationMethod',atomicPhotoIonizationMethod,defaultValue='Verner')
          
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="atomicPhotoIonizationMethod" type="functionCall" functionType="void">
          !#  <functionArgs>atomicPhotoIonizationMethod,Atomic_Cross_Section_Ionization_Photo_Get</functionArgs>
          include 'atomic.cross_sections.ionization.photo.inc'
          !# </include>
          if (.not.associated(Atomic_Cross_Section_Ionization_Photo_Get)) call&
               & Galacticus_Error_Report('Atomic_Cross_Section_Ionization_Photo_Initialize','method '//char(atomicPhotoIonizationMethod)//' is unrecognized')
          ionizationCrossSectionInitialized=.true.
       end if
       !$omp end critical(Atomic_Cross_Section_Ionization_Photo_Initialization) 
    end if
    return
  end subroutine Atomic_Cross_Section_Ionization_Photo_Initialize

  double precision function Atomic_Cross_Section_Ionization_Photo(atomicNumber,ionizationState,shellNumber,wavelength)
    !% Return the cross-section (in units of cm$^2$) for a given atom in a given ionization state at the specified {\tt
    !% wavelength} (given in units of \AA).
    implicit none
    integer,          intent(in) :: atomicNumber,ionizationState,shellNumber
    double precision, intent(in) :: wavelength

    ! Initialize the module.
    call Atomic_Cross_Section_Ionization_Photo_Initialize

    ! Call the routine to do the calculation.
    Atomic_Cross_Section_Ionization_Photo=Atomic_Cross_Section_Ionization_Photo_Get(atomicNumber,ionizationState,shellNumber,wavelength)
    
    return
  end function Atomic_Cross_Section_Ionization_Photo

end module Atomic_Cross_Sections_Ionization_Photo
