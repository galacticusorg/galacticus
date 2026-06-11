!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{RST
  Implementation of a posterior sampling differential evolution proposal size temperature exponent class in which the exponent is adaptive.
  !!}

  !![
  <posteriorSampleDffrntlEvltnPrpslSzTmpExp name="posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive" docformat="rst">
   <description>
   This class adaptively changes :math:`\alpha` in an attempt to maintain the gradient of the acceptance rate with the logarithm of temperature, :math:`\mathrm{d} R/\mathrm{d}\ln T`, at an acceptable level. The algorithm is controlled by the following sub-parameters:

   ``[exponentInitial]``
      The initial value for :math:`\alpha`;

   ``[exponentFactor]``
      The additive factor by which :math:`\alpha` should be increased or decreased if the acceptance rate gradient is out of range;

   ``[exponentMinimum]``
      The smallest value allowed for :math:`\alpha`;

   ``[exponentMaximum]``
      The largest value allowed for :math:`\alpha`;

   ``[acceptanceRateMinimum]``
      The minimum acceptance rate gradient to accept before reducing :math:`\alpha`;

   ``[acceptanceRateMaximum]``
      The maximum acceptance rate gradient to accept before reducing :math:`\alpha`;

   ``[updateCount]``
      The number of steps between successive checks of the acceptance rate gradient.
   </description>
  </posteriorSampleDffrntlEvltnPrpslSzTmpExp>
  !!]
  type, extends(posteriorSampleDffrntlEvltnPrpslSzTmpExpClass) :: posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive
     !!{RST
     Implementation of a posterior sampling differential evolution proposal size class in which the exponent is adaptive.
     !!}
     private
     double precision :: exponentCurrent, exponentAdjustFactor, &
          &              exponentInitial
     double precision :: exponentMinimum, exponentMaximum
     double precision :: gradientMinimum, gradientMaximum
     integer          :: updateCount    , lastUpdateCount
   contains
     procedure :: exponent => adaptiveExponent
  end type posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive

  interface posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive
     !!{RST
     Constructors for the ``posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive`` posterior sampling differential evolution random jump class.
     !!}
     module procedure adaptiveConstructorParameters
     module procedure adaptiveConstructorInternal
  end interface posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive

contains

  function adaptiveConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive`` posterior sampling differential evolution random jump class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive)                 :: self
    type            (inputParameters                                 ), intent(inout)  :: parameters
    double precision                                                                   :: exponentInitial, exponentMinimum     , &
         &                                                                                exponentMaximum, exponentAdjustFactor, &
         &                                                                                gradientMinimum, gradientMaximum
    integer                                                                            :: updateCount

    !![
    <inputParameter docformat="rst">
      <name>exponentInitial</name>
      <description>
      The initial value of the temperature-scaling exponent :math:`\alpha` used before any adaptive adjustment based on the acceptance-rate gradient has been applied.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentMinimum</name>
      <description>
      The minimum value to which the temperature-scaling exponent :math:`\alpha` may be reduced during adaptive adjustment, preventing the temperature dependence from becoming negligibly weak.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentMaximum</name>
      <description>
      The maximum value to which the temperature-scaling exponent :math:`\alpha` may be increased during adaptive adjustment, preventing the proposal size from growing too steeply with temperature.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentAdjustFactor</name>
      <description>
      The additive increment by which the temperature-scaling exponent :math:`\alpha` is increased or decreased at each adaptation step when the acceptance-rate gradient falls outside the target range.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>gradientMinimum</name>
      <description>
      The minimum acceptable gradient of acceptance rate with log-temperature.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>gradientMaximum</name>
      <description>
      The maximum acceptable gradient of acceptance rate with log-temperature.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>updateCount</name>
      <description>
      The number of steps between potential updates of the temperature exponent.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive(exponentInitial,exponentMinimum,exponentMaximum,exponentAdjustFactor,gradientMinimum,gradientMaximum,updateCount)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function adaptiveConstructorParameters

  function adaptiveConstructorInternal(exponentInitial,exponentMinimum,exponentMaximum,exponentAdjustFactor,gradientMinimum,gradientMaximum,updateCount) result(self)
    !!{RST
    Constructor for the ``posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive`` posterior sampling differential evolution random jump class.
    !!}
    implicit none
    type            (posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive)                :: self
    double precision                                                  , intent(in   ) :: exponentInitial, exponentAdjustFactor, &
         &                                                                               exponentMinimum, exponentMaximum     , &
         &                                                                               gradientMinimum, gradientMaximum
    integer                                                           , intent(in   ) :: updateCount
    !![
    <constructorAssign variables="exponentInitial,exponentMinimum,exponentMaximum,exponentAdjustFactor,gradientMinimum,gradientMaximum,updateCount"/>
    !!]

    self%exponentCurrent=exponentInitial
    self%lastUpdateCount=0
    return
  end function adaptiveConstructorInternal

double precision function adaptiveExponent(self,temperedStates,temperatures,simulationState,simulationConvergence)
  !!{RST
  Return the adaptive differential evolution proposal size temperature exponent.
  !!}
  use :: Display           , only : displayIndent     , displayMessage, displayUnindent, displayVerbosity, &
          &                         verbosityLevelInfo
  use :: ISO_Varying_String, only : varying_string
  use :: MPI_Utilities     , only : mpiSelf
  use :: String_Handling   , only : operator(//)
  implicit none
  class           (posteriorSampleDffrntlEvltnPrpslSzTmpExpAdaptive), intent(inout)                                  :: self
  class           (posteriorSampleStateClass                       ), intent(inout)                                  :: simulationState
  class           (posteriorSampleStateClass                       ), intent(inout), dimension(                   :) :: temperedStates
  double precision                                                  , intent(in   ), dimension(                   :) :: temperatures
  class           (posteriorSampleConvergenceClass                 ), intent(inout)                                  :: simulationConvergence
  double precision                                                                 , dimension(size(temperedStates)) :: acceptanceRates      , logTemperatures
  integer                                                                                                            :: levelCount           , i
  double precision                                                                                                   :: gradient
  character       (len=8                                           )                                                 :: label
  type            (varying_string                                  )                                                 :: message
  !$GLC attributes unused :: simulationState

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
     if (mpiSelf%rank() == 0 .and. displayVerbosity() >= verbosityLevelInfo) then
        message='Tempered acceptance rate report after '
        message=message//temperedStates(levelCount)%count()//' tempered steps:'
        call displayIndent(message)
        call displayMessage('Temperature Acceptance Rate')
        call displayMessage('---------------------------')
        do i=1,levelCount
           write (label,'(f8.1)') temperatures   (i)
           message="   "//trim(label)
           write (label,'(f5.3)') acceptanceRates(i)
           message=message//"           "//trim(label)
           call displayMessage(message)
        end do
        call displayMessage('---------------------------')
        write (label,'(f8.3)') gradient
        message="Gradient [dR/dln(T)] = "//trim(label)
        call displayMessage(message)
        call displayUnindent('done')
     end if
     ! If the gradient is out of range, adjust the exponent.
     if      (gradient > self%gradientMaximum .and. self%exponentCurrent < self%exponentMaximum) then
         self%exponentCurrent=min(self%exponentCurrent+self%exponentAdjustFactor,self%exponentMaximum)
         if (mpiSelf%rank() == 0 .and. displayVerbosity() >= verbosityLevelInfo) then
            write (label,'(f8.5)') self%exponentCurrent
            call displayMessage('Adjusting exponent up to '//label)
         end if
     else if (gradient < self%gradientMinimum .and. self%exponentCurrent > self%exponentMinimum) then
         self%exponentCurrent=max(self%exponentCurrent-self%exponentAdjustFactor,self%exponentMinimum)
         if (mpiSelf%rank() == 0 .and. displayVerbosity() >= verbosityLevelInfo) then
            write (label,'(f8.5)') self%exponentCurrent
            call displayMessage('Adjusting exponent down to '//label)
         end if
      end if
  end if
  ! Return the current exponent.
  adaptiveExponent=self%exponentCurrent
  return
end function adaptiveExponent
