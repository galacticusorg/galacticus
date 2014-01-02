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

!% Contains a module which provides calculations of Type Ia supernovae.

module Supernovae_Type_Ia
  !% Provides calculations of Type Ia supernovae.
  use ISO_Varying_String
  implicit none
  private
  public :: SNeIa_Cumulative_Number, SNeIa_Cumulative_Yield

  ! Flag indicating whether this module has been initialized.
  logical                                              :: supernovaeIaInitialized    =.false.

  ! Name of cooling rate available method used.
  type     (varying_string                  )          :: supernovaeIaMethod

  ! Pointer to the function that actually does the calculation.
  procedure(SNeIa_Cumulative_Number_Template), pointer :: SNeIa_Cumulative_Number_Get=>null()
  procedure(SNeIa_Cumulative_Yield_Template ), pointer :: SNeIa_Cumulative_Yield_Get =>null()
  abstract interface
    double precision function SNeIa_Cumulative_Number_Template(initialMass,age,metallicity)
      double precision, intent(in   ) :: age, initialMass, metallicity
    end function SNeIa_Cumulative_Number_Template
  end interface
  abstract interface
    double precision function SNeIa_Cumulative_Yield_Template(initialMass,age,metallicity,atomIndex)
      double precision, intent(in   )           :: age      , initialMass, metallicity
      integer         , intent(in   ), optional :: atomIndex
    end function SNeIa_Cumulative_Yield_Template
  end interface

contains

  subroutine Supernovae_Type_Ia_Initialize
    !% Initialize the Type Ia supernovae module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="supernovaeIaMethod" type="moduleUse">
    include 'stellar_astrophysics.supernovae_type_Ia.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.supernovaeIaInitialized) then
       !$omp critical(Supernovae_Type_Ia_Initialization)
       if (.not.supernovaeIaInitialized) then
          ! Get the halo spin distribution method parameter.
          !@ <inputParameter>
          !@   <name>supernovaeIaMethod</name>
          !@   <defaultValue>Nagashima</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The method to use for computing properties of Type Ia supernovae.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('supernovaeIaMethod',supernovaeIaMethod,defaultValue='Nagashima')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="supernovaeIaMethod" type="functionCall" functionType="void">
          !#  <functionArgs>supernovaeIaMethod,SNeIa_Cumulative_Number_Get,SNeIa_Cumulative_Yield_Get</functionArgs>
          include 'stellar_astrophysics.supernovae_type_Ia.inc'
          !# </include>
          if (.not.(associated(SNeIa_Cumulative_Number_Get).and.associated(SNeIa_Cumulative_Yield_Get))) call Galacticus_Error_Report('Supernovae_Type_Ia_Initialize'&
               &,'method '//char(supernovaeIaMethod)//' is unrecognized')
          supernovaeIaInitialized=.true.
       end if
       !$omp end critical(Supernovae_Type_Ia_Initialization)
    end if
    return
  end subroutine Supernovae_Type_Ia_Initialize

  double precision function SNeIa_Cumulative_Number(initialMass,age,metallicity)
    !% Return the cumulative number of Type Ia supernovae from stars of given {\tt initialMass}, {\tt age} and {\tt metallicity}.
    implicit none
    double precision, intent(in   ) :: age, initialMass, metallicity

    ! Ensure module is initialized.
    call Supernovae_Type_Ia_Initialize

    ! Simply call the function which does the actual work.
    SNeIa_Cumulative_Number=SNeIa_Cumulative_Number_Get(initialMass,age,metallicity)
    return
  end function SNeIa_Cumulative_Number

  double precision function SNeIa_Cumulative_Yield(initialMass,age,metallicity,atomIndex)
    !% Return the cumulative yield of Type Ia supernovae from stars of given {\tt initialMass}, {\tt age} and {\tt metallicity}.
    implicit none
    double precision, intent(in   )           :: age      , initialMass, metallicity
    integer         , intent(in   ), optional :: atomIndex

    ! Ensure module is initialized.
    call Supernovae_Type_Ia_Initialize

    ! Simply call the function which does the actual work.
    SNeIa_Cumulative_Yield=SNeIa_Cumulative_Yield_Get(initialMass,age,metallicity,atomIndex)
    return
  end function SNeIa_Cumulative_Yield

end module Supernovae_Type_Ia
