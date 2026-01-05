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
     double precision :: gammaCurrent            , gammaAdjustFactor    , &
          &              gammaInitial
     double precision :: gammaMinimum            , gammaMaximum
     double precision :: acceptanceRateMinimum   , acceptanceRateMaximum
     integer          :: updateCount             , lastUpdateCount
     logical          :: outliersInAcceptanceRate
   contains
     procedure :: gamma => adaptiveGamma
  end type posteriorSampleDffrntlEvltnProposalSizeAdaptive

  interface posteriorSampleDffrntlEvltnProposalSizeAdaptive
     !!{
     Constructors for the \refClass{posteriorSampleDffrntlEvltnProposalSizeAdaptive} posterior sampling differential evolution random jump class.
     !!}
     module procedure adaptiveConstructorParameters
     module procedure adaptiveConstructorInternal
  end interface posteriorSampleDffrntlEvltnProposalSizeAdaptive

contains

  function adaptiveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleDffrntlEvltnProposalSizeAdaptive} posterior sampling differential evolution random jump class which
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
    logical                                                                           :: outliersInAcceptanceRate

    !![
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
    !!]
    self=posteriorSampleDffrntlEvltnProposalSizeAdaptive(gammaInitial,gammaMinimum,gammaMaximum,gammaAdjustFactor,acceptanceRateMinimum,acceptanceRateMaximum,updateCount,outliersInAcceptanceRate)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function adaptiveConstructorParameters

  function adaptiveConstructorInternal(gammaInitial,gammaMinimum,gammaMaximum,gammaAdjustFactor,acceptanceRateMinimum,acceptanceRateMaximum,updateCount,outliersInAcceptanceRate) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleDffrntlEvltnProposalSizeAdaptive} differential evolution proposal size class.
    !!}
    implicit none
    type            (posteriorSampleDffrntlEvltnProposalSizeAdaptive)                :: self
    double precision                                                 , intent(in   ) :: gammaInitial            , gammaAdjustFactor    , &
         &                                                                              gammaMinimum            , gammaMaximum         , &
         &                                                                              acceptanceRateMinimum   , acceptanceRateMaximum
    integer                                                          , intent(in   ) :: updateCount
    logical                                                          , intent(in   ) :: outliersInAcceptanceRate
    !![
    <constructorAssign variables="gammaInitial,gammaMinimum,gammaMaximum,gammaAdjustFactor,acceptanceRateMinimum,acceptanceRateMaximum,updateCount,outliersInAcceptanceRate"/>
    !!]

    self%gammaCurrent   =gammaInitial
    self%lastUpdateCount=0
    return
  end function adaptiveConstructorInternal

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
       if (mpiSelf%rank() == 0 .and. displayVerbosity() >= verbosityLevelStandard) then
          if (acceptanceRate < 0.0d0) then
             label="unknown"
          else
             write (label,'(f5.3)') acceptanceRate
          end if
          message='After '
          message=message//simulationState%count()//' steps, acceptance rate is '//trim(label)
          call displayMessage(message)
       end if
       ! If the acceptance rate is out of range, adjust γ.
       if (acceptanceRate >= 0.0d0) then
          if      (acceptanceRate > self%acceptanceRateMaximum .and. self%gammaCurrent < self%gammaMaximum) then
             self%gammaCurrent=min(self%gammaCurrent*self%gammaAdjustFactor,self%gammaMaximum)
             if (mpiSelf%rank() == 0 .and. displayVerbosity() >= verbosityLevelStandard) then
                write (label,'(f8.5)') self%gammaCurrent
                call displayMessage('Adjusting γ up to '//label)
             end if
          else if (acceptanceRate < self%acceptanceRateMinimum .and. self%gammaCurrent > self%gammaMinimum) then
             self%gammaCurrent=max(self%gammaCurrent/self%gammaAdjustFactor,self%gammaMinimum)
             if (mpiSelf%rank() == 0 .and. displayVerbosity() >= verbosityLevelStandard) then
                write (label,'(f8.5)') self%gammaCurrent
                call displayMessage('Adjusting γ down to '//label)
             end if
          end if
       end if
    end if
    ! Return the current adaptive size.
    adaptiveGamma=self%gammaCurrent
    return
  end function adaptiveGamma
