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

!% Contains a module which implements running a simulation to determine constraints.

module Constraints_Constrain
  !% Implements running a simulation to determine constraints.
  private
  public :: Constrain
  
contains
  
  subroutine Constrain(configFile)
    !% Run a simulation to determine constraints on a model.
    use ISO_Varying_String
    use FoX_DOM
    use IO_XML
    use Galacticus_Error
    use Statistics_Distributions
    use Constraints_Priors
    use Constraints_Likelihoods
    use Constraints_Convergence
    use Constraints_State
    use Constraints_Simulation
    use Constraints_Differential_Proposal_Size
    use Constraints_Differential_Random_Jump
    use System_Command
    use MPI_Utilities
    implicit none
    type   (varying_string  ), intent(in   )               :: configFile
    type   (prior           ), allocatable  , dimension(:) :: parameterPriors
    type   (distributionList), allocatable  , dimension(:) :: randomDistributions
    class  (likelihood      ), pointer                     :: modelLikelihood
    class  (convergence     ), pointer                     :: simulationConvergence
    class  (state           ), pointer                     :: simulationState
    class  (deProposalSize  ), pointer                     :: proposalSize
    class  (deRandomJump    ), pointer                     :: randomJump
    class  (simulator       ), pointer                     :: simulation
    type   (node            ), pointer                     :: configDoc            , priorDefinition       , &
         &                                                    parameterDefinition  , randomDefinition      , &
         &                                                    likelihoodDefinition , convergenceDefinition , &
         &                                                    stateDefinition      , proposalSizeDefinition, &
         &                                                    simulationDefinition , parametersElement     , &
         &                                                    randomJumpDefinition
    type   (nodeList        ), pointer                     :: parameterDefinitions , parametersList
    type   (varying_string  )                              :: filterCommand        , filteredFile
    integer                                                :: parameterCount       , ioError               , &
         &                                                    i                    , j                     , &
         &                                                    iParameter           , inactiveParameterCount

    ! Run the config file through an external XInclude filter to include any Xinclude'd files.
    filteredFile="/dev/shm/"//trim(configFile)//"_"//mpiSelf%rankLabel()
    filterCommand="scripts/aux/xmlInclude.pl "//trim(configFile)//" "//filteredFile
    call System_Command_Do(filterCommand)
    ! Parse the simulation config file.
    configDoc => parseFile(char(filteredFile),iostat=ioError)
    if (ioError /= 0) call Galacticus_Error_Report('Constrain','Unable to find or parse config file')
    call System_Command_Do("rm -f "//filteredFile)
    ! Determine the number of parameters.
    parameterCount         =  0
    inactiveParameterCount =  0
    parametersList         => getElementsByTagName(configDoc,"parameters")
    do i=0,getLength(parametersList)-1
       parametersElement    => item(parametersList,i)
       parameterDefinitions => getElementsByTagName(parametersElement,"parameter")
       do j=1,getLength(parameterDefinitions)
          parameterDefinition => item(parameterDefinitions,j-1)
          if (XML_Path_Exists(parameterDefinition,"prior")) then
             parameterCount        =parameterCount        +1
          else
             inactiveParameterCount=inactiveParameterCount+1
          end if
      end do
    end do
    if (parameterCount <= 0) call Galacticus_Error_Report('Constrain','at least one parameter must be specified in config file')
    if (mpiSelf%isMaster()) write (0,*) 'Found ',parameterCount,' active parameters (and ',inactiveParameterCount,' inactive parameters)'
    ! Initialize priors and random perturbers.
    allocate(parameterPriors    (parameterCount))
    allocate(randomDistributions(parameterCount))
    parametersList => getElementsByTagName(configDoc,"parameters")
    iParameter=0
    do i=0,getLength(parametersList)-1
       parametersElement    => item(parametersList,i)
       parameterDefinitions => getElementsByTagName(parametersElement,"parameter")
       do j=1,getLength(parameterDefinitions)
          parameterDefinition => item(parameterDefinitions,j-1)
          if (XML_Path_Exists(parameterDefinition,"prior")) then
             iParameter=iParameter+1
             priorDefinition                                  => XML_Get_First_Element_By_Tag_Name(parameterDefinition ,"prior"   )
             parameterPriors    (iParameter)                  =  prior                            (priorDefinition     ,iParameter)
             randomDefinition                                 => XML_Get_First_Element_By_Tag_Name(parameterDefinition ,"random"  )
             randomDistributions(iParameter)%thisDistribution => distributionNew                  (randomDefinition               )
          end if
       end do
    end do
    ! Initialize likelihood.
    likelihoodDefinition   => XML_Get_First_Element_By_Tag_Name(configDoc,"likelihood"  )
    modelLikelihood        =>     likelihoodNew(  likelihoodDefinition,configFile    )
    ! Initialize convergence.
    convergenceDefinition => XML_Get_First_Element_By_Tag_Name(configDoc,"convergence"  )
    simulationConvergence =>     convergenceNew( convergenceDefinition               )
    ! Initialize state.
    stateDefinition        => XML_Get_First_Element_By_Tag_Name(configDoc,"state"       )
    simulationState        =>          stateNew(       stateDefinition,parameterCount)
    ! Initialize proposal size.
    proposalSizeDefinition => XML_Get_First_Element_By_Tag_Name(configDoc,"proposalSize")
    proposalSize           => deProposalSizeNew(proposalSizeDefinition               )
    ! Initialize random jump.
    randomJumpDefinition   => XML_Get_First_Element_By_Tag_Name(configDoc,"randomJump"  )
    randomJump             => deRandomJumpNew  (randomJumpDefinition                 )
    ! Initialize simulation.
    simulationDefinition   => XML_Get_First_Element_By_Tag_Name(configDoc,"simulation"  )
    simulation             =>      simulatorNew(                       &
         &                                      simulationDefinition , &
         &                                      parameterPriors      , &
         &                                      randomDistributions  , &
         &                                      modelLikelihood      , &
         &                                      simulationConvergence, &
         &                                      simulationState      , &
         &                                      proposalSize         , &
         &                                      randomJump             &
         &                                     )
    ! Destroy the simulation config document.
    call destroy(configDoc)
    ! Perform the simulation.
    call simulation%simulate()
    ! Clean up.
    deallocate(parameterPriors      )
    deallocate(modelLikelihood      )
    deallocate(simulationConvergence)
    deallocate(simulationState      )
    deallocate(proposalSize         )
    deallocate(randomJump           )
    deallocate(simulation           )
    do i=1,parameterCount
       deallocate(randomDistributions(i)%thisDistribution)
    end do
    return
  end subroutine Constrain

end module Constraints_Constrain
