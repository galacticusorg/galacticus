!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements the first crossing distribution (assuming a linear barrier) for excursion set calculations
!% of dark matter halo formation.

module Excursion_Sets_First_Crossing_Linear_Barrier
  !% Implements the first crossing distribution (assuming a linear barrier) for excursion set calculations of dark matter halo
  !% formation.
  private
  public :: Excursion_Sets_First_Crossing_Linear_Barrier_Initialize

contains

  !# <excursionSetFirstCrossingMethod>
  !#  <unitName>Excursion_Sets_First_Crossing_Linear_Barrier_Initialize</unitName>
  !# </excursionSetFirstCrossingMethod>
  subroutine Excursion_Sets_First_Crossing_Linear_Barrier_Initialize(excursionSetFirstCrossingMethod&
       &,Excursion_Sets_First_Crossing_Probability_Get,Excursion_Sets_First_Crossing_Rate_Get&
       &,Excursion_Sets_Non_Crossing_Rate_Get)
    !% Initialize the linear barrier first crossing distribution for excursion sets module.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: excursionSetFirstCrossingMethod
    procedure(double precision), pointer, intent(inout) :: Excursion_Sets_First_Crossing_Probability_Get&
         &,Excursion_Sets_First_Crossing_Rate_Get,Excursion_Sets_Non_Crossing_Rate_Get

    if (excursionSetFirstCrossingMethod == 'linearBarrier') then
       Excursion_Sets_First_Crossing_Probability_Get => Excursion_Sets_First_Crossing_Probability_Linear
       Excursion_Sets_First_Crossing_Rate_Get        => Excursion_Sets_First_Crossing_Rate_Linear
       Excursion_Sets_Non_Crossing_Rate_Get          => Excursion_Sets_Non_Crossing_Rate_Linear
    end if
    return
  end subroutine Excursion_Sets_First_Crossing_Linear_Barrier_Initialize

  double precision function Excursion_Sets_First_Crossing_Probability_Linear(variance,time)
    !% Return the probability for excursion set first crossing assuming a linear barrier. Uses the analytic solution for this case
    !% \cite{sheth_excursion_1998,sheth_excursion_2002}.
    use Numerical_Constants_Math
    use Excursion_Sets_Barriers
    implicit none
    double precision, intent(in) :: variance,time

    Excursion_Sets_First_Crossing_Probability_Linear=Excursion_Sets_Barrier(0.0d0,time)*exp(-0.5d0&
         &*Excursion_Sets_Barrier(variance,time)**2/variance)/variance/sqrt(2.0d0*Pi*variance)
    return
  end function Excursion_Sets_First_Crossing_Probability_Linear
  
  double precision function Excursion_Sets_First_Crossing_Rate_Linear(variance,varianceProgenitor,time)
    !% Return the rate for excursion set first crossing assuming a linear barrier. Uses the analytic solution for this case
    !% \cite{sheth_excursion_1998,sheth_excursion_2002} with a simple offset in the starting coordinates. The rate of barrier
    !% crossing is computed by solving for the first crossing distribution at a slightly earlier time and then dividing through by
    !% that time interval.
    use Numerical_Constants_Math
    use Excursion_Sets_Barriers
    implicit none
    double precision, intent(in) :: variance,varianceProgenitor,time
    double precision, parameter  :: fractionalTimeChange=1.0d-3
    double precision             :: timeProgenitor

    ! Compute a slightly earlier time for the progenitor
    timeProgenitor=time*(1.0d0-fractionalTimeChange)
    if (variance >= varianceProgenitor) then
       Excursion_Sets_First_Crossing_Rate_Linear=0.0d0
    else
       Excursion_Sets_First_Crossing_Rate_Linear=Excursion_Sets_Barrier_Effective(variance,time,variance,timeProgenitor)*exp(&
            &-0.5d0*Excursion_Sets_Barrier_Effective(variance,time,varianceProgenitor,timeProgenitor)**2/(varianceProgenitor&
            &-variance))/(varianceProgenitor-variance)/sqrt(2.0d0*Pi*(varianceProgenitor-variance))/time/fractionalTimeChange
    end if
    return
  end function Excursion_Sets_First_Crossing_Rate_Linear
  
  double precision function Excursion_Sets_Non_Crossing_Rate_Linear(variance,time)
    !% Return the rate for excursion set non-crossing assuming a linear barrier. For a linear barrier the integral over the
    !% crossing probability (from zero to infinite variance) equals unity, so all trajectories cross. The non-crossing rate is
    !% therefore zero.
    use Numerical_Constants_Math
    use Excursion_Sets_Barriers
    implicit none
    double precision, intent(in) :: variance,time

    Excursion_Sets_Non_Crossing_Rate_Linear=0.0d0
    return
  end function Excursion_Sets_Non_Crossing_Rate_Linear
  
  double precision function Excursion_Sets_Barrier_Effective(variance0,time0,variance,time)
    !% The effective barrier for conditional excursion sets.
    use Excursion_Sets_Barriers
    implicit none
    double precision, intent(in) :: variance0,variance,time0,time
    
    Excursion_Sets_Barrier_Effective=                                         &
         &                            Excursion_Sets_Barrier(variance ,time ) &
         &                           -Excursion_Sets_Barrier(variance0,time0)
    return
  end function Excursion_Sets_Barrier_Effective

end module Excursion_Sets_First_Crossing_Linear_Barrier
