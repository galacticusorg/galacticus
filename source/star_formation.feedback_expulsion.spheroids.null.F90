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

!% Contains a module which implements a null expulsive outflow rate in galactic spheroids.

module Star_Formation_Expulsive_Feedback_Spheroids_Null
  !% Implementss a null expulsive outflow rate in galactic spheroids.
  implicit none
  private
  public :: Star_Formation_Expulsive_Feedback_Spheroids_Null_Initialize

contains

  !# <starFormationExpulsiveFeedbackSpheroidsMethod>
  !#  <unitName>Star_Formation_Expulsive_Feedback_Spheroids_Null_Initialize</unitName>
  !# </starFormationExpulsiveFeedbackSpheroidsMethod>
  subroutine Star_Formation_Expulsive_Feedback_Spheroids_Null_Initialize(starFormationExpulsiveFeedbackSpheroidsMethod&
       &,Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Get)
    !% Initializes the ``null'' spheroid star formation expulsive feedback module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                              ), intent(in   )          :: starFormationExpulsiveFeedbackSpheroidsMethod
    procedure(Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Null), intent(inout), pointer :: Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Get

    if (starFormationExpulsiveFeedbackSpheroidsMethod == 'null') Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Get => Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Null

    return
  end subroutine Star_Formation_Expulsive_Feedback_Spheroids_Null_Initialize

  double precision function Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Null(thisNode,starFormationRate,energyInputRate)
    !% Implements a null expulsive outflow rate for spheroids.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: energyInputRate, starFormationRate

    ! Return a zero outflow rate.
    Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Null=0.0d0
    return
  end function Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Null

end module Star_Formation_Expulsive_Feedback_Spheroids_Null
