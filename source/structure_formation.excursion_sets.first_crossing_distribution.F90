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


!% Contains a module which implements calculations of first crossing distributions for excursion set calculations.

module Excursion_Sets_First_Crossings
  !% Implements calculations of first crossing distributions for excursion set calculations.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Excursion_Sets_First_Crossing_Probability,Excursion_Sets_First_Crossing_Rate,Excursion_Sets_Collapsed_Fraction&
       &,Excursion_Sets_Non_Crossing_Rate

  ! Flag to indicate if this module has been initialized.  
  logical              :: firstCrossingModuleInitalized=.false.

  ! Name of method to use for distribution function.
  type(varying_string) :: excursionSetFirstCrossingMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Excursion_Sets_First_Crossing_Probability), pointer :: Excursion_Sets_First_Crossing_Probability_Get => null()
  procedure(Excursion_Sets_First_Crossing_Rate       ), pointer :: Excursion_Sets_First_Crossing_Rate_Get        => null()
  procedure(Excursion_Sets_Non_Crossing_Rate         ), pointer :: Excursion_Sets_Non_Crossing_Rate_Get          => null()

contains

  subroutine Excursion_Sets_First_Crossings_Initialize
    !% Initialize the excursion sets first crossing distribution module.
    use Input_Parameters
    use Galacticus_Error
    !# <include directive="excursionSetFirstCrossingMethod" type="moduleUse">
    include 'structure_formation.excursion_sets.first_crossing_distribution.moduleUse.inc'
    !# </include>
    implicit none

    !$omp critical(Excursion_Sets_First_Crossing_Initialization) 
    ! Initialize if necessary.
    if (.not.firstCrossingModuleInitalized) then
       ! Get the barrier first crossing probability method parameter.
       !@ <inputParameter>
       !@   <name>excursionSetFirstCrossingMethod</name>
       !@   <defaultValue>linearBarrier</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of first crossing distributions for excursion sets.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('excursionSetFirstCrossingMethod',excursionSetFirstCrossingMethod,defaultValue='linearBarrier')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="excursionSetFirstCrossingMethod" type="code" action="subroutine">
       !#  <subroutineArgs>excursionSetFirstCrossingMethod,Excursion_Sets_First_Crossing_Probability_Get,Excursion_Sets_First_Crossing_Rate_Get,Excursion_Sets_Non_Crossing_Rate_Get</subroutineArgs>
       include 'structure_formation.excursion_sets.first_crossing_distribution.inc'
       !# </include>
       if     (                                                                     &
            &  .not.(                                                               &
            &             associated(Excursion_Sets_First_Crossing_Probability_Get) &
            &        .and.associated(Excursion_Sets_First_Crossing_Rate_Get       ) &
            &        .and.associated(Excursion_Sets_Non_Crossing_Rate_Get         ) &
            &       )                                                               &
            & ) call Galacticus_Error_Report('Excursion_Sets_Barrier_Initialize','method '//char(excursionSetFirstCrossingMethod)//' is unrecognized')
       firstCrossingModuleInitalized=.true.
    end if
    !$omp end critical(Excursion_Sets_First_Crossing_Initialization) 
    return
  end subroutine Excursion_Sets_First_Crossings_Initialize

  double precision function Excursion_Sets_First_Crossing_Probability(variance,time)
    !% Return the probability of first crossing for excursion sets at the given {\tt variance} and {\tt time}.
    implicit none
    double precision, intent(in) :: variance,time

    ! Initialize the module if necessary.
    call Excursion_Sets_First_Crossings_Initialize

    Excursion_Sets_First_Crossing_Probability=Excursion_Sets_First_Crossing_Probability_Get(variance,time)
    return
  end function Excursion_Sets_First_Crossing_Probability
  
  double precision function Excursion_Sets_First_Crossing_Rate(variance,varianceProgenitor,time)
    !% Return the rate of first crossing for excursion sets beginning at the given {\tt variance} and {\tt time} to transition to a first crossing at the given {\tt varianceProgenitor}.
    implicit none
    double precision, intent(in) :: variance,varianceProgenitor,time

    ! Initialize the module if necessary.
    call Excursion_Sets_First_Crossings_Initialize

    Excursion_Sets_First_Crossing_Rate=Excursion_Sets_First_Crossing_Rate_Get(variance,varianceProgenitor,time)
    return
  end function Excursion_Sets_First_Crossing_Rate

  double precision function Excursion_Sets_Non_Crossing_Rate(variance,time)
    !% Return the rate of non-crossing for excursion sets beginning at the given {\tt variance} and {\tt time}.
    implicit none
    double precision, intent(in) :: variance,time

    ! Initialize the module if necessary.
    call Excursion_Sets_First_Crossings_Initialize

    Excursion_Sets_Non_Crossing_Rate=Excursion_Sets_Non_Crossing_Rate_Get(variance,time)
    return
  end function Excursion_Sets_Non_Crossing_Rate

end module Excursion_Sets_First_Crossings

