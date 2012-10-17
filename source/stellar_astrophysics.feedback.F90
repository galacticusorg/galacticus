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

!% Contains a module which provides calculations of stellar feedback.

module Stellar_Feedback
  !% Provides calculations of stellar feedback.
  use ISO_Varying_String
  implicit none
  private
  public :: Stellar_Feedback_Cumulative_Energy_Input

  ! Flag indicating whether this module has been initialized.
  logical              :: stellarFeedbackInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: stellarFeedbackMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Stellar_Feedback_Cumulative_Energy_Input_Get_Template), pointer :: Stellar_Feedback_Cumulative_Energy_Input_Get => null()
  abstract interface
    double precision function Stellar_Feedback_Cumulative_Energy_Input_Get_Template(initialMass,age,metallicity)
      double precision, intent(in) :: initialMass,age,metallicity
    end function Stellar_Feedback_Cumulative_Energy_Input_Get_Template
  end interface

  ! Canonical value of the total energy input from a single stellar population of $1 M_\odot$ after infinite time. All feedback
  ! calculations which don't specifically use the energy input should be scaled to this value if they want to have the correct
  ! time and IMF dependencies. Value was computed for a Salpeter IMF. Units are MSolar (km/s)^2.
  double precision, parameter, public :: feedbackEnergyInputAtInfinityCanonical=4.517d5

contains

  subroutine Stellar_Feedback_Initialize
    !% Initialize the stellar feedback module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="stellarFeedbackMethod" type="moduleUse">
    include 'stellar_astrophysics.feedback.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.stellarFeedbackInitialized) then
       !$omp critical(Stellar_Feedback_Initialization) 
       if (.not.stellarFeedbackInitialized) then
          ! Get the halo spin distribution method parameter.
          !@ <inputParameter>
          !@   <name>stellarFeedbackMethod</name>
          !@   <defaultValue>standard</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The method to use for computing aspects of stellar feedback.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('stellarFeedbackMethod',stellarFeedbackMethod,defaultValue='standard')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="stellarFeedbackMethod" type="code" action="subroutine">
          !#  <subroutineArgs>stellarFeedbackMethod,Stellar_Feedback_Cumulative_Energy_Input_Get</subroutineArgs>
          include 'stellar_astrophysics.feedback.inc'
          !# </include>
          if (.not.associated(Stellar_Feedback_Cumulative_Energy_Input_Get)) call Galacticus_Error_Report('Stellar_Feedback_Initialize'&
               &,'method '//char(stellarFeedbackMethod)//' is unrecognized')
          stellarFeedbackInitialized=.true.
       end if
       !$omp end critical(Stellar_Feedback_Initialization) 
    end if
    return
  end subroutine Stellar_Feedback_Initialize

  double precision function Stellar_Feedback_Cumulative_Energy_Input(initialMass,age,metallicity)
    !% Return the cumulative energy input per from stellar feedback from stars of given {\tt initialMass}, {\tt age} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: initialMass,age,metallicity

    ! Ensure module is initialized.
    call Stellar_Feedback_Initialize

    ! Simply call the function which does the actual work.
    Stellar_Feedback_Cumulative_Energy_Input=Stellar_Feedback_Cumulative_Energy_Input_Get(initialMass,age,metallicity)
    return
  end function Stellar_Feedback_Cumulative_Energy_Input

end module Stellar_Feedback
