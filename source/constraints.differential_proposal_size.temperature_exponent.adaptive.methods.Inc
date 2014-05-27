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

!% Defines methods for the ``adaptive'' differential evolution proposal size temperature exponent class.

function dePropSizeTempExpAdaptiveConstructor(exponent,exponentMinimum,exponentMaximum,exponentAdjustFactor,gradientMinimum,gradientMaximum,updateCount)
  !% Constructor for the ``adaptive'' differential evolution proposal size temperature exponent class.
  implicit none
  type            (dePropSizeTempExpAdaptive)                :: dePropSizeTempExpAdaptiveConstructor
  double precision                           , intent(in   ) :: exponent                            , exponentMinimum     , &
       &                                                        exponentMaximum                     , exponentAdjustFactor, &
       &                                                        gradientMinimum                     , gradientMaximum
  integer                                    , intent(in   ) :: updateCount

  dePropSizeTempExpAdaptiveConstructor%exponentCurrent     =exponent
  dePropSizeTempExpAdaptiveConstructor%exponentAdjustFactor=exponentAdjustFactor
  dePropSizeTempExpAdaptiveConstructor%exponentMinimum     =exponentMinimum
  dePropSizeTempExpAdaptiveConstructor%exponentMaximum     =exponentMaximum
  dePropSizeTempExpAdaptiveConstructor%gradientMinimum     =gradientMinimum
  dePropSizeTempExpAdaptiveConstructor%gradientMaximum     =gradientMaximum
  dePropSizeTempExpAdaptiveConstructor%updateCount         =updateCount
  dePropSizeTempExpAdaptiveConstructor%lastUpdateCount     =0
  return
end function dePropSizeTempExpAdaptiveConstructor

double precision function dePropSizeTempExpAdaptiveExponent(self,temperedStates,temperatures,simulationState,simulationConvergence)
  !% Return the adaptive differential evolution proposal size temperature exponent.
  use MPI_Utilities
  use ISO_Varying_String
  use Galacticus_Display
  use String_Handling
  implicit none
  class           (dePropSizeTempExpAdaptive), intent(inout)                                  :: self
  class           (state                    ), intent(in   )                                  :: simulationState
  class           (state                    ), intent(in   ), dimension(                   :) :: temperedStates
  double precision                           , intent(in   ), dimension(                   :) :: temperatures
  class           (convergence              ), intent(inout)                                  :: simulationConvergence
  double precision                                          , dimension(size(temperedStates)) :: acceptanceRates      , logTemperatures
  integer                                                                                     :: levelCount           , i
  double precision                                                                            :: gradient
  character       (len=8                    )                                                 :: label
  type            (varying_string           )                                                 :: message
 
  ! Should we consider updating the exponent?
  levelCount=size(temperedStates)
  if     (                                                                                               &
       &        temperedStates       (levelCount)%count      () >= self%lastUpdateCount+self%updateCount &
       &  .and.                                                                                          &
       &   .not.simulationConvergence            %isConverged()                                          &
       & ) then
     ! Reset the number of steps remaining.
     self%lastUpdateCount=temperedStates(levelCount)%count()
     ! Find the mean acceptance rates across all chains.
     do i=1,levelCount
        acceptanceRates(i)=temperedStates(i)%acceptanceRate()
     end do
     acceptanceRates=mpiSelf%average(acceptanceRates)
     ! Find the gradient of the acceptance rate with log temperature.
     logTemperatures=log(temperatures)
     gradient       =                                                 &
          &           (                                               &
          &            +sum (logTemperatures    *    acceptanceRates) &
          &            -sum (logTemperatures   )*sum(acceptanceRates) &
          &            /dble(levelCount        )                      &
          &           )                                               &
          &          /(                                               &
          &            +sum (logTemperatures**2)                      &
          &            -sum (logTemperatures   )**2                   &
          &            /dble(levelCount        )                      &
          &          )
     ! Report. 
     if (mpiSelf%rank() == 0 .and. Galacticus_Verbosity_Level() >= verbosityInfo) then
        message='Tempered acceptance rate report after '
        message=message//temperedStates(levelCount)%count()//' tempered steps:'
        call Galacticus_Display_Indent(message)
        call Galacticus_Display_Message('Temperature Acceptance Rate')
        call Galacticus_Display_Message('---------------------------')
        do i=1,levelCount
           write (label,'(f8.1)') temperatures   (i)
           message="   "//trim(label)
           write (label,'(f5.3)') acceptanceRates(i)
           message=message//"           "//trim(label)
           call Galacticus_Display_Message(message)
        end do
        call Galacticus_Display_Message('---------------------------')
        write (label,'(f8.3)') gradient
        message="Gradient [dR/dln(T)] = "//trim(label)
        call Galacticus_Display_Message(message)
        call Galacticus_Display_Unindent('done')
     end if
     ! If the gradient is out of range, adjust the exponent.
     if      (gradient > self%gradientMaximum .and. self%exponentCurrent < self%exponentMaximum) then
         self%exponentCurrent=min(self%exponentCurrent+self%exponentAdjustFactor,self%exponentMaximum)
         if (mpiSelf%rank() == 0 .and. Galacticus_Verbosity_Level() >= verbosityInfo) then
            write (label,'(f8.5)') self%exponentCurrent
            call Galacticus_Display_Message('Adjusting exponent up to '//label)
         end if
     else if (gradient < self%gradientMinimum .and. self%exponentCurrent > self%exponentMinimum) then
         self%exponentCurrent=max(self%exponentCurrent-self%exponentAdjustFactor,self%exponentMinimum)
         if (mpiSelf%rank() == 0 .and. Galacticus_Verbosity_Level() >= verbosityInfo) then
            write (label,'(f8.5)') self%exponentCurrent
            call Galacticus_Display_Message('Adjusting exponent down to '//label)
         end if
      end if
  end if
  ! Return the current exponent.
  dePropSizeTempExpAdaptiveExponent=self%exponentCurrent
  return
end function dePropSizeTempExpAdaptiveExponent
