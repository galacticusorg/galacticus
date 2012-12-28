!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of first crossing distributions for excursion set calculations.

module Excursion_Sets_First_Crossings
  !% Implements calculations of first crossing distributions for excursion set calculations.
  use ISO_Varying_String
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
       !# <include directive="excursionSetFirstCrossingMethod" type="functionCall" functionType="void">
       !#  <functionArgs>excursionSetFirstCrossingMethod,Excursion_Sets_First_Crossing_Probability_Get,Excursion_Sets_First_Crossing_Rate_Get,Excursion_Sets_Non_Crossing_Rate_Get</functionArgs>
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

