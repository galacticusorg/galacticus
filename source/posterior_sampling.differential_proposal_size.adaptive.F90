!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  Implementation of a posterior sampling differential evolution proposal size class in which the proposal size is adaptive.
  !!}

  !![
  <posteriorSampleDffrntlEvltnProposalSize name="posteriorSampleDffrntlEvltnProposalSizeAdaptive">
   <description>
    This class adaptively changes $\gamma$ in an attempt to maintain the acceptance rate at an acceptable level. The algorithm is
    controlled by the following sub-parameters:
    \begin{description}
    \item[{\normalfont \ttfamily [gammaInitial]}] The initial value for $\gamma$.
    \item[{\normalfont \ttfamily [gammaFactor]}] The multiplicative factor by which $\gamma$ should be increased or decreased if the
      acceptance rate is out of range.
    \item[{\normalfont \ttfamily [gammaMinimum]}] The smallest value allowed for $\gamma$.
    \item[{\normalfont \ttfamily [gammaMaximum]}] The largest value allowed for $\gamma$.
    \item[{\normalfont \ttfamily [acceptanceRateMinimum]}] The minimum acceptance rate to accept before reducing $\gamma$.
    \item[{\normalfont \ttfamily [acceptanceRateMaximum]}] The maximum acceptance rate to accept before reducing $\gamma$.
    \item[{\normalfont \ttfamily [updateCount]}] The number of steps between successive checks of the acceptance rate.
    \end{description}
   </description>
  </posteriorSampleDffrntlEvltnProposalSize>
  !!]
  type, extends(posteriorSampleDffrntlEvltnProposalSizeClass) :: posteriorSampleDffrntlEvltnProposalSizeAdaptive
     !!{
     Implementation of a posterior sampling differential evolution proposal size class in which the proposal size is adaptive.
     !!}
     private
     double precision                 :: gammaCurrent            , gammaAdjustFactor    , &
          &                              gammaInitial
     double precision                 :: gammaMinimum            , gammaMaximum
     double precision                 :: acceptanceRateMinimum   , acceptanceRateMaximum
     integer                          :: updateCount             , lastUpdateCount
     logical                          :: outliersInAcceptanceRate, appendLog            , &
          &                              restoreFromLog          , flushLog
     type            (varying_string) :: logFileName
     integer                          :: logFileUnit
   contains
     final     ::          adaptiveDestructor
     procedure :: gamma => adaptiveGamma
  end type posteriorSampleDffrntlEvltnProposalSizeAdaptive

  interface posteriorSampleDffrntlEvltnProposalSizeAdaptive
     !!{
     Constructors for the {\normalfont \ttfamily adaptive} posterior sampling differential evolution random jump class.
     !!}
     module procedure adaptiveConstructorParameters
     module procedure adaptiveConstructorInternal
  end interface posteriorSampleDffrntlEvltnProposalSizeAdaptive

contains

  function adaptiveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily adaptive} posterior sampling differential evolution random jump class which
    builds the object from a parameter set.
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
    <inputParameter>
      <name>logFileName</name>
      <description>The name of a file to which to log reports of adjustments to $\gamma$. If empty, no reports are logged.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gammaInitial</name>
      <description>The initial proposal size, $\gamma$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gammaMinimum</name>
      <description>The minimum allowed proposal size, $\gamma$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gammaMaximum</name>
      <description>The maximum allowed proposal size, $\gamma$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gammaAdjustFactor</name>
      <description>The factor by which to adjust the proposal size, $\gamma$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>acceptanceRateMinimum</name>
      <description>The minimum acceptable acceptance rate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>acceptanceRateMaximum</name>
      <description>The maximum acceptable acceptance rate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>updateCount</name>
      <description>The number of steps between potential updates of the proposal size.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outliersInAcceptanceRate</name>
      <defaultValue>.true.</defaultValue>
      <description>The number of steps between potential updates of the proposal size.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>appendLog</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, append to the existing log file, otherwise overwrite.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>restoreFromLog</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, restore the value of $\gamma$ from the log file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>flushLog</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, logs are flushed to file after every update.</description>
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
    !!{
    Constructor for the {\normalfont \ttfamily adaptive} differential evolution proposal size class.
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
    integer                                                                          :: ioStatus
    !![
    <constructorAssign variables="logFileName,gammaInitial,gammaMinimum,gammaMaximum,gammaAdjustFactor,acceptanceRateMinimum,acceptanceRateMaximum,updateCount,outliersInAcceptanceRate, appendLog, restoreFromLog, flushLog"/>
    !!]

    self%gammaCurrent   =gammaInitial
    if (self%restoreFromLog .and. logFileName /= '') then
       open(newunit=self%logFileUnit,file=char(logFileName),status='old',form='formatted',iostat=ioStatus)
       do while (ioStatus == 0)
          read (self%logFileUnit,'(a)',iostat=ioStatus) line
          if (ioStatus /= 0) exit
          if (index(line,"Adjusting") /= 0) read (line(index(trim(line)," ",back=.true.)+1:len_trim(line)),*) self%gammaCurrent
       end do
       close(self%logFileUnit)
       call mpiBarrier()
    end if
    self%lastUpdateCount=0
    if (mpiSelf%rank() == 0 .and. logFileName /= '') then
       if (self%appendLog) then
          open(newunit=self%logFileUnit,file=char(logFileName),status='unknown',form='formatted',position='append')
       else
          open(newunit=self%logFileUnit,file=char(logFileName),status='unknown',form='formatted'                  )
       end if
    else
       self%logFileUnit=-huge(0)
    end if
    return
  end function adaptiveConstructorInternal

  subroutine adaptiveDestructor(self)
    !!{
    Destructor for the ``adaptive'' differential evolution proposal size class.
    !!}
    implicit none
    type(posteriorSampleDffrntlEvltnProposalSizeAdaptive), intent(inout) :: self
    
    if (self%logFileUnit /= -huge(0)) close(self%logFileUnit)
    return
  end subroutine adaptiveDestructor
  
  double precision function adaptiveGamma(self,simulationState,simulationConvergence)
    !!{
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
          if (self%logFileUnit /= -huge(0)) then
             write (self%logFileUnit,*) char(message)
             if (self%flushLog) call flush(self%logFileUnit)
          end if
       end if
       ! If the acceptance rate is out of range, adjust γ.
       if (acceptanceRate >= 0.0d0) then
          if      (acceptanceRate > self%acceptanceRateMaximum .and. self%gammaCurrent < self%gammaMaximum) then
             self%gammaCurrent=min(self%gammaCurrent*self%gammaAdjustFactor,self%gammaMaximum)
             if (mpiSelf%rank() == 0) then
                write (label,'(f8.5)') self%gammaCurrent
                if (displayVerbosity() >= verbosityLevelStandard) call displayMessage('Adjusting γ up to '//label)
                if (self%logFileUnit /= -huge(0)) then
                   write (self%logFileUnit,*) 'Adjusting γ up to '//label
                   if (self%flushLog) call flush(self%logFileUnit)
                end if
             end if
          else if (acceptanceRate < self%acceptanceRateMinimum .and. self%gammaCurrent > self%gammaMinimum) then
             self%gammaCurrent=max(self%gammaCurrent/self%gammaAdjustFactor,self%gammaMinimum)
             if (mpiSelf%rank() == 0) then
                write (label,'(f8.5)') self%gammaCurrent
                if (displayVerbosity() >= verbosityLevelStandard) call displayMessage('Adjusting γ down to '//label)
                if (self%logFileUnit /= -huge(0)) then
                   write (self%logFileUnit,*) 'Adjusting γ down to '//label
                   if (self%flushLog) call flush(self%logFileUnit)
                end if
             end if
          end if
       end if
    end if
    ! Return the current adaptive size.
    adaptiveGamma=self%gammaCurrent
    return
  end function adaptiveGamma
