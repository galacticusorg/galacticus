!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module which implements a excursion set first crossing statistics class for linear barriers.

  use Excursion_Sets_Barriers
  
  !# <excursionSetFirstCrossing name="excursionSetFirstCrossingLinearBarrier">
  !#  <description>An excursion set first crossing statistics class for linear barriers.</description>
  !# </excursionSetFirstCrossing>
  type, extends(excursionSetFirstCrossingClass) :: excursionSetFirstCrossingLinearBarrier
     !% A linearBarrier excursion set barrier class.
     private
     class(excursionSetBarrierClass), pointer :: excursionSetBarrier_
   contains
     final     ::                    linearBarrierDestructor
     procedure :: probability     => linearBarrierProbability
     procedure :: rate            => linearBarrierRate
     procedure :: rateNonCrossing => linearBarrierRateNonCrossing
  end type excursionSetFirstCrossingLinearBarrier

  interface excursionSetFirstCrossingLinearBarrier
     !% Constructors for the linearBarrier excursion set barrier class.
     module procedure linearBarrierConstructorParameters
     module procedure linearBarrierConstructorInternal
  end interface excursionSetFirstCrossingLinearBarrier

contains

  function linearBarrierConstructorParameters(parameters) result(self)
    !% Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (excursionSetFirstCrossingLinearBarrier)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(excursionSetBarrierClass              ), pointer       :: excursionSetBarrier_

    !# <objectBuilder class="excursionSetBarrier" name="excursionSetBarrier_" source="parameters"/>
    self=excursionSetFirstCrossingLinearBarrier(excursionSetBarrier_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function linearBarrierConstructorParameters

  function linearBarrierConstructorInternal(excursionSetBarrier_) result(self)
    !% Constructor for the linear barrier excursion set class first crossing class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (excursionSetFirstCrossingLinearBarrier)                        :: self
    class(excursionSetBarrierClass              ), intent(in   ), target :: excursionSetBarrier_
    !# <constructorAssign variables="*excursionSetBarrier_"/>

    return
  end function linearBarrierConstructorInternal

  subroutine linearBarrierDestructor(self)
    !% Destructor for the critical overdensity excursion set barrier class.
    implicit none
    type(excursionSetFirstCrossingLinearBarrier), intent(inout) :: self

    !# <objectDestructor name="self%excursionSetBarrier_" />
    return
  end subroutine linearBarrierDestructor
  
  double precision function linearBarrierProbability(self,variance,time)
    !% Return the excursion set barrier at the given variance and time.
    use Numerical_Constants_Math
    implicit none
    class           (excursionSetFirstCrossingLinearBarrier), intent(inout) :: self
    double precision                                        , intent(in   ) :: variance, time

    linearBarrierProbability=+     self%excursionSetBarrier_%barrier(   0.0d0,time,rateCompute=.false.)    &
         &                   *exp(                                                                         &
         &                        -0.5d0                                                                   &
         &                        *self%excursionSetBarrier_%barrier(variance,time,rateCompute=.false.)**2 &
         &                        /variance                                                                &
         &                       )                                                                         &
         &                   /variance                                                                     &
         &                   /sqrt(                                                                        &
         &                         +2.0d0                                                                  &
         &                         *Pi                                                                     &
         &                         *variance                                                               &
         &                        )
    return
  end function linearBarrierProbability

  double precision function linearBarrierRate(self,variance,varianceProgenitor,time)
    !% Return the excursion set barrier at the given variance and time.
    use Numerical_Constants_Math
    implicit none
    class           (excursionSetFirstCrossingLinearBarrier), intent(inout) :: self
    double precision                                        , intent(in   ) :: variance                   , varianceProgenitor, &
         &                                                                     time
    double precision                                        , parameter     :: fractionalTimeChange=1.0d-3
    double precision                                                        :: timeProgenitor

    ! Compute a slightly earlier time for the progenitor
    timeProgenitor=time*(1.0d0-fractionalTimeChange)
    if (variance >= varianceProgenitor) then
       linearBarrierRate=0.0d0
    else
       linearBarrierRate=+     barrierEffective(variance,time,variance          ,timeProgenitor)    &
            &            *exp(                                                                      &
            &                 -0.5d0                                                                &
            &                 *barrierEffective(variance,time,varianceProgenitor,timeProgenitor)**2 &
            &                 / (+varianceProgenitor-variance)                                      &
            &                 )                                                                     &
            &            /      (+varianceProgenitor-variance)                                      &
            &            /sqrt(                                                                     &
            &                  +2.0d0                                                               &
            &                  *Pi                                                                  &
            &                  *(+varianceProgenitor-variance)                                      &
            &                 )                                                                     &
            &            /time                                                                      &
            &            /fractionalTimeChange
    end if
    return

  contains

    double precision function barrierEffective(variance0,time0,variance1,time1)
      !% The effective barrier for conditional excursion sets.
      implicit none
      double precision, intent(in   ) :: time1    , time0    , &
           &                             variance1, variance0

      barrierEffective=+self%excursionSetBarrier_%barrier(variance1,time1,rateCompute=.false.) &
           &                        -self%excursionSetBarrier_%barrier(variance0,time0,rateCompute=.false.)
      return
    end function barrierEffective

  end function linearBarrierRate

  double precision function linearBarrierRateNonCrossing(self,variance,time)
    !% Return the rate for excursion set non-crossing assuming a linearBarrier barrier. For a linearBarrier barrier the integral over the
    !% crossing probability (from zero to infinite variance) equals unity, so all trajectories cross. The non-crossing rate is
    !% therefore zero.
    implicit none
    class           (excursionSetFirstCrossingLinearBarrier), intent(inout) :: self
    double precision                                        , intent(in   ) :: time, variance
    !GCC$ attributes unused :: self, time, variance
    
    linearBarrierRateNonCrossing=0.0d0
    return
  end function linearBarrierRateNonCrossing
