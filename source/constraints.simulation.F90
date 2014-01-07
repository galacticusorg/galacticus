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

!% Contains a module which implements simulation algorithms for use when constraining \glc.

module Constraints_Simulation
  !% Implements simulation algorithms for use when constraining \glc.
  use Statistics_Distributions
  use Constraints_Priors
  use Constraints_Likelihoods
  use Constraints_Convergence
  use Constraints_State
  use Constraints_Differential_Proposal_Size
  use ISO_Varying_String
  private
  public :: simulatorNew

  ! Define the basic simulator class.
  type, abstract, public :: simulator
   contains
     !@ <objectMethods>
     !@   <object>simulator</object>
     !@   <objectMethod>
     !@     <method>simulate</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Perform the simulation.</description>
     !@   </objectMethod>
     !@ </objectMethods>
    procedure(simulatorSimulate), deferred :: simulate
  end type simulator

  ! Interface for deferred functions.
  abstract interface
     subroutine simulatorSimulate(self)
       import :: simulator
       class(simulator), intent(inout) :: self
     end subroutine simulatorSimulate
  end interface

  ! Include all simualtion types.
  include 'constraints.simulation.differential_evolution.type.inc'
  include 'constraints.simulation.tempered_differential_evolution.type.inc'

contains

  function simulatorNew(definition,parameterPriors,randomDistributions,modelLikelihood,simulationConvergence ,simulationState&
       &,proposalSize) result (newSimulator)
    !% Create a new differential evolution proposal size from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    class           (simulator       )               , pointer                        :: newSimulator
    type            (node            ), intent(in   ), pointer                        :: definition
    type            (prior           ), intent(in   ), optional, target, dimension(:) :: parameterPriors
    type            (distributionList), intent(in   ), optional, target, dimension(:) :: randomDistributions
    class           (likelihood      ), intent(in   ), optional, target               :: modelLikelihood
    class           (convergence     ), intent(in   ), optional, target               :: simulationConvergence
    class           (state           ), intent(in   ), optional, target               :: simulationState
    class           (deProposalSize  ), intent(in   ), optional, target               :: proposalSize
    type            (node            ), pointer                                       :: simulatorStepsMaximumDefinition          , simulatorStepsPostConvergenceDefinition, &
         &                                                                               simulatorAcceptanceAverageCountDefinition, simulatorLogFileDefinition             , &
         &                                                                               simulatorTemperatureMaximumDefinition    , simulatorUntemperedStepCountDefinition , &
         &                                                                               simulatorTemperingLevelCountDefinition   , simulatorStepsPerLevelDefinition       , &
         &                                                                               simulatorExponentDefinition
    integer                                                                           :: simulatorStepsMaximum                    , simulatorStepsPostConvergence          , &
         &                                                                               simulatorAcceptanceAverageCount          , simulatorUntemperedStepCount           , &
         &                                                                               simulatorTemperingLevelCount             , simulatorStepsPerLevel
    double precision                                                                  :: simulatorTemperatureMaximum              , simulatorExponent
    type            (varying_string  )                                                :: simulatorLogFile

    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("differentialEvolution")
       allocate(simulatorDifferentialEvolution :: newSimulator)
       select type (newSimulator)
       type is (simulatorDifferentialEvolution)
          simulatorStepsMaximumDefinition           => XML_Get_First_Element_By_Tag_Name(definition,"stepsMaximum"          )
          simulatorStepsPostConvergenceDefinition   => XML_Get_First_Element_By_Tag_Name(definition,"stepsPostConvergence"  )
          simulatorAcceptanceAverageCountDefinition => XML_Get_First_Element_By_Tag_Name(definition,"acceptanceAverageCount")
          simulatorLogFileDefinition                => XML_Get_First_Element_By_Tag_Name(definition,"logFileRoot"           )
          call extractDataContent(simulatorStepsMaximumDefinition          ,simulatorStepsMaximum          )
          call extractDataContent(simulatorStepsPostConvergenceDefinition  ,simulatorStepsPostConvergence  )
          call extractDataContent(simulatorAcceptanceAverageCountDefinition,simulatorAcceptanceAverageCount)
          simulatorLogFile=XML_Extract_Text(simulatorLogFileDefinition)
          newSimulator=simulatorDifferentialEvolution(                                 &
               &                                      parameterPriors                , &
               &                                      randomDistributions            , &
               &                                      modelLikelihood                , &
               &                                      simulationConvergence          , &
               &                                      simulationState                , &
               &                                      proposalSize                   , &
               &                                      simulatorStepsMaximum          , &
               &                                      simulatorStepsPostConvergence  , &
               &                                      simulatorAcceptanceAverageCount, &
               &                                      char(simulatorLogFile)           &
               &                                     )
       end select
    case ("temperedDifferentialEvolution")
       allocate(simulatorTemperedDifferentialEvolution :: newSimulator)
       select type (newSimulator)
       type is (simulatorTemperedDifferentialEvolution)
          simulatorStepsMaximumDefinition           => XML_Get_First_Element_By_Tag_Name(definition,"stepsMaximum"            )
          simulatorStepsPostConvergenceDefinition   => XML_Get_First_Element_By_Tag_Name(definition,"stepsPostConvergence"    )
          simulatorAcceptanceAverageCountDefinition => XML_Get_First_Element_By_Tag_Name(definition,"acceptanceAverageCount"  )
          simulatorLogFileDefinition                => XML_Get_First_Element_By_Tag_Name(definition,"logFileRoot"             )
          simulatorTemperatureMaximumDefinition     => XML_Get_First_Element_By_Tag_Name(definition,"temperatureMaximum"      )
          simulatorUntemperedStepCountDefinition    => XML_Get_First_Element_By_Tag_Name(definition,"untemperedStepCount"     )
          simulatorTemperingLevelCountDefinition    => XML_Get_First_Element_By_Tag_Name(definition,"temperedLevels"          )
          simulatorStepsPerLevelDefinition          => XML_Get_First_Element_By_Tag_Name(definition,"stepsPerLevel"           )
          simulatorExponentDefinition               => XML_Get_First_Element_By_Tag_Name(definition,"gammaTemperatureExponent")
          call extractDataContent(simulatorStepsMaximumDefinition          ,simulatorStepsMaximum          )
          call extractDataContent(simulatorStepsPostConvergenceDefinition  ,simulatorStepsPostConvergence  )
          call extractDataContent(simulatorAcceptanceAverageCountDefinition,simulatorAcceptanceAverageCount)
          call extractDataContent(simulatorTemperatureMaximumDefinition    ,simulatorTemperatureMaximum    )
          call extractDataContent(simulatorUntemperedStepCountDefinition   ,simulatorUntemperedStepCount   )
          call extractDataContent(simulatorTemperingLevelCountDefinition   ,simulatorTemperingLevelCount   )
          call extractDataContent(simulatorStepsPerLevelDefinition         ,simulatorStepsPerLevel         )
          call extractDataContent(simulatorExponentDefinition              ,simulatorExponent              )
          simulatorLogFile=XML_Extract_Text(simulatorLogFileDefinition)
          newSimulator=simulatorTemperedDifferentialEvolution(                                 &
               &                                              parameterPriors                , &
               &                                              randomDistributions            , &
               &                                              modelLikelihood                , &
               &                                              simulationConvergence          , &
               &                                              simulationState                , &
               &                                              proposalSize                   , &
               &                                              simulatorStepsMaximum          , &
               &                                              simulatorStepsPostConvergence  , &
               &                                              simulatorAcceptanceAverageCount, &
               &                                              char(simulatorLogFile)         , &
               &                                              simulatorTemperatureMaximum    , &
               &                                              simulatorUntemperedStepCount   , &
               &                                              simulatorTemperingLevelCount   , &
               &                                              simulatorStepsPerLevel         , &
               &                                              simulatorExponent                &
               &                                             )
       end select
     case default
       call Galacticus_Error_Report('simulatorNew','simulator type is unrecognized')
    end select
    return
  end function simulatorNew

  ! Include all simulation methods.
  include 'constraints.simulation.differential_evolution.methods.inc'
  include 'constraints.simulation.tempered_differential_evolution.methods.inc'

end module Constraints_Simulation
