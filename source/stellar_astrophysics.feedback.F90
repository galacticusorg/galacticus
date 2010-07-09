!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which provides calculations of stellar feedback.

module Stellar_Feedback
  !% Provides calculations of stellar feedback.
  use ISO_Varying_String
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

    !$omp critical(Stellar_Feedback_Initialization) 
    ! Initialize if necessary.
    if (.not.stellarFeedbackInitialized) then
       ! Get the halo spin distribution method parameter.
       !@ <inputParameter>
       !@   <name>stellarFeedbackMethod</name>
       !@   <defaultValue>standard</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The method to use for computing aspects of stellar feedback.
       !@   </description>
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
