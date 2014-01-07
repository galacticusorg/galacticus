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
    implicit none
    type   (varying_string  ), intent(in   )               :: configFile
    type   (prior           ), allocatable  , dimension(:) :: parameterPriors
    type   (distributionList), allocatable  , dimension(:) :: randomDistributions
    class  (likelihood      ), pointer                     :: modelLikelihood
    class  (convergence     ), pointer                     :: simulationConvergence
    class  (state           ), pointer                     :: simulationState
    class  (deProposalSize  ), pointer                     :: proposalSize
    class  (simulator       ), pointer                     :: simulation
    type   (node            ), pointer                     :: configDoc            , priorDefinition       , &
         &                                                    parameterDefinition  , randomDefinition      , &
         &                                                    likelihoodDefinition , convergenceDefinition , &
         &                                                    stateDefinition      , proposalSizeDefinition, &
         &                                                    simulationDefinition , parametersElement
    type   (nodeList        ), pointer                     :: parameterDefinitions
    integer                                                :: parameterCount       , ioError               , &
         &                                                    i

    !! AJB: TODO
    ! Document.
    ! Tempered evolution:
    !   When tempering can we cheat: Instead of actually dividing the likelihood by the
    !   temperature, what if we just run a smaller number of merger tree realizations? This
    !   would make the model covariance matrix larger, giving the same effect as a higher
    !   temperature. We'd have to artificially increase the data covariance matrix by the same
    !   factor. Probably we would then multiply the resulting likelihood by the temperature so
    !   that once it was passed back to the simulation object it would get re-divided by
    !   temperature there. Could save a lot of run time when tempering.
    ! Adaptive random perturbers - e.g. Cauchy with width equal to some fraction of the current min/max range between chains in the parameters.
    ! Add likelihood emulator.

    ! Parse the simulation config file.
    configDoc => parseFile(char(configFile),iostat=ioError)
    if (ioError /= 0) call Galacticus_Error_Report('Constrain','Unable to find or parse config file')
    ! Determine the number of parameters.
    parametersElement => XML_Get_First_Element_By_Tag_Name(configDoc        ,"parameters")
    parameterCount    =  XML_Array_Length                 (parametersElement,"parameter" )
    if (parameterCount <= 0) call Galacticus_Error_Report('Constrain','at least one parameter must be specified in config file')
    ! Initialize priors.
    allocate(parameterPriors(parameterCount))
    parameterDefinitions => getElementsByTagName(parametersElement,"parameter")
    do i=1,parameterCount
       parameterDefinition => item                             (parameterDefinitions,i-1    )
       priorDefinition     => XML_Get_First_Element_By_Tag_Name(parameterDefinition ,"prior")
       parameterPriors(i)  =  prior                            (priorDefinition     ,i      )
    end do
    ! Initialize random perturbers.
    allocate(randomDistributions(parameterCount))
    do i=1,parameterCount
       parameterDefinition => item(parameterDefinitions   ,i-1    )
       randomDefinition    => XML_Get_First_Element_By_Tag_Name(parameterDefinition ,"random")
       randomDistributions(i)%thisDistribution => distributionNew(randomDefinition)
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
    ! Initialize simulation.
    simulationDefinition   => XML_Get_First_Element_By_Tag_Name(configDoc,"simulation"  )
    simulation             =>      simulatorNew(                       &
         &                                      simulationDefinition , &
         &                                      parameterPriors      , &
         &                                      randomDistributions  , &
         &                                      modelLikelihood      , &
         &                                      simulationConvergence, &
         &                                      simulationState      , &
         &                                      proposalSize           &
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
    deallocate(simulation           )
    do i=1,parameterCount
       deallocate(randomDistributions(i)%thisDistribution)
    end do
    return
  end subroutine Constrain

end module Constraints_Constrain
