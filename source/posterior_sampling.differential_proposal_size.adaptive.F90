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
  Implementation of a posterior sampling differential evolution proposal size class in which the proposal size is adaptive.
  !!}

  use :: File_Utilities, only : file

  !![
  <posteriorSampleDffrntlEvltnProposalSize name="posteriorSampleDffrntlEvltnProposalSizeAdaptive" docformat="rst">
   <description>
   This class adaptively changes :math:`\gamma` in an attempt to maintain the acceptance rate at an acceptable level. The algorithm is controlled by the following sub-parameters:

   ``[gammaInitial]``
      The initial value for :math:`\gamma`.

   ``[gammaFactor]``
      The multiplicative factor by which :math:`\gamma` should be increased or decreased if the acceptance rate is out of range.

   ``[gammaMinimum]``
      The smallest value allowed for :math:`\gamma`.

   ``[gammaMaximum]``
      The largest value allowed for :math:`\gamma`.

   ``[acceptanceRateMinimum]``
      The minimum acceptance rate to accept before reducing :math:`\gamma`.

   ``[acceptanceRateMaximum]``
      The maximum acceptance rate to accept before reducing :math:`\gamma`.

   ``[updateCount]``
      The number of steps between successive checks of the acceptance rate.
   </description>
  </posteriorSampleDffrntlEvltnProposalSize>
  !!]
  type, extends(posteriorSampleDffrntlEvltnProposalSizeClass) :: posteriorSampleDffrntlEvltnProposalSizeAdaptive
     !!{RST
     Implementation of a posterior sampling differential evolution proposal size class in which the proposal size is adaptive.
     !!}
     private
     double precision                  :: gammaCurrent            , gammaAdjustFactor    , &
          &                               gammaInitial
     double precision                  :: gammaMinimum            , gammaMaximum
     double precision                  :: acceptanceRateMinimum   , acceptanceRateMaximum
     integer                           :: updateCount             , lastUpdateCount
     logical                           :: outliersInAcceptanceRate, appendLog            , &
          &                               restoreFromLog          , flushLog
     type            (varying_string ) :: logFileName
     type            (file           ) :: logFile
   contains
     procedure :: gamma => adaptiveGamma
  end type posteriorSampleDffrntlEvltnProposalSizeAdaptive

  interface posteriorSampleDffrntlEvltnProposalSizeAdaptive
     !!{RST
     Constructors for the ``posteriorSampleDffrntlEvltnProposalSizeAdaptive`` posterior sampling differential evolution random jump class.
     !!}
     module procedure adaptiveConstructorParameters
     module procedure adaptiveConstructorInternal
  end interface posteriorSampleDffrntlEvltnProposalSizeAdaptive

contains

  function adaptiveConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``posteriorSampleDffrntlEvltnProposalSizeAdaptive`` posterior sampling differential evolution random jump class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleDffrntlEvltnProposalSizeAdaptive)                 :: self
    type            (inputParameters                                ), intent(inout)  :: parameters
    double precision                                                                  :: gammaInitial            , gammaMinimum         , &
         &                                                                               gammaMaximum            , gammaAdjustFactor    , &
         &                                                                               acceptanceRateMinimum   , acceptanceRateMaximum
    integer                                                                           :: updateCount
    logical                                                                           :: outliersInAcceptanceRate, appendLog            , &
          &                                                                              restoreFromLog          , flushLog
    type            (varying_string                                 )                 :: logFileName

    !![
    <inputParameter docformat="rst">
      <name>logFileName</name>
      <description>
      The name of a file to which to log reports of adjustments to :math:`\gamma`. If empty, no reports are logged.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>gammaInitial</name>
      <description>
      The initial value of the proposal scaling parameter :math:`\gamma` used before the acceptance rate has been assessed and any adaptive adjustment has been made.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>gammaMinimum</name>
      <description>
      The minimum value to which the proposal scaling parameter :math:`\gamma` is permitted to be reduced during adaptive adjustment, preventing the step size from becoming vanishingly small.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>gammaMaximum</name>
      <description>
      The maximum value to which the proposal scaling parameter :math:`\gamma` is permitted to be increased during adaptive adjustment, preventing excessively large steps that would degrade acceptance rates.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>gammaAdjustFactor</name>
      <description>
      The multiplicative factor by which :math:`\gamma` is increased or decreased at each adaptation step when the current acceptance rate falls outside the target range.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>acceptanceRateMinimum</name>
      <description>
      The minimum acceptable chain acceptance rate; if the measured acceptance rate falls below this threshold :math:`\gamma` is reduced to produce smaller, more easily accepted proposals.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>acceptanceRateMaximum</name>
      <description>
      The maximum acceptable chain acceptance rate; if the measured acceptance rate exceeds this threshold :math:`\gamma` is increased to produce larger proposals that explore the posterior more efficiently.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>updateCount</name>
      <description>
      The number of steps between potential updates of the proposal size.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outliersInAcceptanceRate</name>
      <defaultValue>.true.</defaultValue>
      <description>
      The number of steps between potential updates of the proposal size.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>appendLog</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, append to the existing log file, otherwise overwrite.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>restoreFromLog</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, restore the value of :math:`\gamma` from the log file.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>flushLog</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, logs are flushed to file after every update.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleDffrntlEvltnProposalSizeAdaptive(logFileName,gammaInitial,gammaMinimum,gammaMaximum,gammaAdjustFactor,acceptanceRateMinimum,acceptanceRateMaximum,updateCount,outliersInAcceptanceRate,appendLog,restoreFromLog,flushLog)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function adaptiveConstructorParameters

  function adaptiveConstructorInternal(logFileName,gammaInitial,gammaMinimum,gammaMaximum,gammaAdjustFactor,acceptanceRateMinimum,acceptanceRateMaximum,updateCount,outliersInAcceptanceRate,appendLog,restoreFromLog,flushLog) result(self)
    !!{RST
    Constructor for the ``posteriorSampleDffrntlEvltnProposalSizeAdaptive`` posterior sampling differential evolution random jump class.
    !!}
    use :: MPI_Utilities, only : mpiSelf, mpiBarrier
    implicit none
    type            (posteriorSampleDffrntlEvltnProposalSizeAdaptive)                :: self
    type            (varying_string                                 ), intent(in   ) :: logFileName
    double precision                                                 , intent(in   ) :: gammaInitial            , gammaAdjustFactor    , &
         &                                                                              gammaMinimum            , gammaMaximum         , &
         &                                                                              acceptanceRateMinimum   , acceptanceRateMaximum
    integer                                                          , intent(in   ) :: updateCount
    logical                                                          , intent(in   ) :: outliersInAcceptanceRate, appendLog            , &
         &                                                                              restoreFromLog          , flushLog
    character       (len=32                                         )                :: line
    integer                                                                          :: ioStatus                , logFileUnit
    !![
    <constructorAssign variables="logFileName,gammaInitial,gammaMinimum,gammaMaximum,gammaAdjustFactor,acceptanceRateMinimum,acceptanceRateMaximum,updateCount,outliersInAcceptanceRate, appendLog, restoreFromLog, flushLog"/>
    !!]

    self%gammaCurrent   =gammaInitial
    if (self%restoreFromLog .and. logFileName /= '') then
       open(newunit=logFileUnit,file=char(logFileName),status='old',form='formatted',iostat=ioStatus)
       do while (ioStatus == 0)
          read (logFileUnit,'(a)',iostat=ioStatus) line
          if (ioStatus /= 0) exit
          if (index(line,"Adjusting") /= 0) read (line(index(trim(line)," ",back=.true.)+1:len_trim(line)),*) self%gammaCurrent
       end do
       close(logFileUnit)
       call mpiBarrier()
    end if
    self%lastUpdateCount=0
    if (mpiSelf%rank() == 0 .and. logFileName /= '') then
       if (self%appendLog) then
          self%logFile=file(logFileName,form='formatted',status='unknown',position='append')
       else
          self%logFile=file(logFileName,form='formatted',status='unknown'                  )
       end if
    end if
    return
  end function adaptiveConstructorInternal

  double precision function adaptiveGamma(self,simulationState,simulationConvergence)
    !!{RST
    Return the proposal size.
    !!}
    use :: Display           , only : displayMessage, displayVerbosity, verbosityLevelStandard
    use :: ISO_Varying_String, only : varying_string
    use :: MPI_Utilities     , only : mpiSelf
    use :: String_Handling   , only : operator(//)
    implicit none
    class           (posteriorSampleDffrntlEvltnProposalSizeAdaptive), intent(inout)              :: self
    class           (posteriorSampleStateClass                      ), intent(inout)              :: simulationState
    class           (posteriorSampleConvergenceClass                ), intent(inout)              :: simulationConvergence
    double precision                                                                              :: acceptanceRate
    character       (len=8                                          )                             :: label
    type            (varying_string                                 )                             :: message
    logical                                                          , dimension(:) , allocatable :: areOutliers
    
    ! Should we consider updating γ?
    if     (                                                                                   &
         &        simulationState      %count      () >= self%lastUpdateCount+self%updateCount &
         &  .and.                                                                              &
         &   .not.simulationConvergence%isConverged()                                          &
         & ) then
       ! Reset the number of steps remaining.
       self%lastUpdateCount=simulationState%count()
       ! Find the mean acceptance rate across all chains.

       if (self%outliersInAcceptanceRate) then
          acceptanceRate=mpiSelf%average(simulationState%acceptanceRate())
       else
          areOutliers=mpiSelf%gather(simulationConvergence%stateIsOutlier(mpiSelf%rank()))
          if (all(areOutliers)) then
             acceptanceRate=-huge(0.0d0)
          else
             acceptanceRate=mpiSelf%average(simulationState%acceptanceRate(),mask=.not.areOutliers)
          end if
       end if
       if (mpiSelf%rank() == 0) then
          if (acceptanceRate < 0.0d0) then
             label="unknown"
          else
             write (label,'(f5.3)') acceptanceRate
          end if
          message='After '
          message=message//simulationState%count()//' steps, acceptance rate is '//trim(label)
          if (displayVerbosity() >= verbosityLevelStandard) call displayMessage(message)
          if (associated(self%logFile%unit)) then
             write (self%logFile%unit,*) char(message)
             if (self%flushLog) call flush(self%logFile%unit)
          end if
       end if
       ! If the acceptance rate is out of range, adjust γ.
       if (acceptanceRate >= 0.0d0) then
          if      (acceptanceRate > self%acceptanceRateMaximum .and. self%gammaCurrent < self%gammaMaximum) then
             self%gammaCurrent=min(self%gammaCurrent*self%gammaAdjustFactor,self%gammaMaximum)
             if (mpiSelf%rank() == 0) then
                write (label,'(f8.5)') self%gammaCurrent
                if (displayVerbosity() >= verbosityLevelStandard) call displayMessage('Adjusting γ up to '//label)
                if (associated(self%logFile%unit)) then
                   write (self%logFile%unit,*) 'Adjusting γ up to '//label
                   if (self%flushLog) call flush(self%logFile%unit)
                end if
             end if
          else if (acceptanceRate < self%acceptanceRateMinimum .and. self%gammaCurrent > self%gammaMinimum) then
             self%gammaCurrent=max(self%gammaCurrent/self%gammaAdjustFactor,self%gammaMinimum)
             if (mpiSelf%rank() == 0) then
                write (label,'(f8.5)') self%gammaCurrent
                if (displayVerbosity() >= verbosityLevelStandard) call displayMessage('Adjusting γ down to '//label)
                if (associated(self%logFile%unit)) then
                   write (self%logFile%unit,*) 'Adjusting γ down to '//label
                   if (self%flushLog) call flush(self%logFile%unit)
                end if
             end if
          end if
       end if
    end if
    ! Return the current adaptive size.
    adaptiveGamma=self%gammaCurrent
    return
  end function adaptiveGamma
