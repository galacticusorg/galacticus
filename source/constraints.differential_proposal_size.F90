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

!% Contains a module which implements algorithms for the proposal size in differential evolution algorithms.

module Constraints_Differential_Proposal_Size
  !% Implements algorithms for the proposal size in differential evolution algorithms.
  use Constraints_State
  use Constraints_Convergence
  private
  public :: deProposalSizeNew

  ! Define the basic simulator class.
  type, abstract, public :: deProposalSize
   contains
     !@ <objectMethods>
     !@   <object>deProposalSize</object>
     !@   <objectMethod>
     !@     <method>gamma</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless class(state)\textgreater} simulationState\argin, \textcolor{red}{\textless class(convergence)\textgreater} simulationConvergence\arginout</arguments>
     !@     <description>Return the proposal step size for differential evolution simulationsx2.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(deProposalSizeGamma), deferred :: gamma
  end type deProposalSize

  ! Interface for deferred functions.
  abstract interface
     double precision function deProposalSizeGamma(self,simulationState,simulationConvergence)
       import :: deProposalSize, state, convergence
       class(deProposalSize), intent(inout) :: self
       class(state         ), intent(in   ) :: simulationState
       class(convergence   ), intent(inout) :: simulationConvergence
    end function deProposalSizeGamma
  end interface

  ! Include all proposal size types.
  include 'constraints.differential_proposal_size.fixed.type.inc'
  include 'constraints.differential_proposal_size.adaptive.type.inc'

contains

  function deProposalSizeNew(definition) result (newDeProposalSize)
    !% Create a new differential evolution proposal size from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    class  (deProposalSize), pointer                :: newDeProposalSize
    type   (node          ), pointer, intent(in   ) :: definition
    type   (node          ), pointer                :: deProposalSizeGammaInitialDefinition         , deProposalSizeGammaFactorDefinition          , &
         &                                             deProposalSizeAccpetanceRateMinimumDefinition, deProposalSizeAccpetanceRateMaximumDefinition, &
         &                                             deProposalSizeUpdateCountDefinition
    double precision                                :: deProposalSizeGammaInitial                   , deProposalSizeGammaFactor                    , &
         &                                             deProposalSizeAccpetanceRateMinimum          , deProposalSizeAccpetanceRateMaximum
    integer                                         :: deProposalSizeUpdateCount
   
    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("fixed")
       allocate(deProposalSizeFixed :: newDeProposalSize)
       select type (newDeProposalSize)
       type is (deProposalSizeFixed)
          deProposalSizeGammaInitialDefinition => XML_Get_First_Element_By_Tag_Name(definition,"gammaInitial")
          call extractDataContent(deProposalSizeGammaInitialDefinition,deProposalSizeGammaInitial)
          newDeProposalSize=deProposalSizeFixed(deProposalSizeGammaInitial)
       end select
    case ("adaptive")
       allocate(deProposalSizeAdaptive :: newDeProposalSize)
       select type (newDeProposalSize)
       type is (deProposalSizeAdaptive)
          deProposalSizeGammaInitialDefinition          => XML_Get_First_Element_By_Tag_Name(definition,"gammaInitial"         )
          deProposalSizeGammaFactorDefinition           => XML_Get_First_Element_By_Tag_Name(definition,"gammaFactor"          )
          deProposalSizeAccpetanceRateMinimumDefinition => XML_Get_First_Element_By_Tag_Name(definition,"acceptanceRateMinimum")
          deProposalSizeAccpetanceRateMaximumDefinition => XML_Get_First_Element_By_Tag_Name(definition,"acceptanceRateMaximum")
          deProposalSizeUpdateCountDefinition           => XML_Get_First_Element_By_Tag_Name(definition,"updateCount"          )
          call extractDataContent(deProposalSizeGammaInitialDefinition         ,deProposalSizeGammaInitial         )
          call extractDataContent(deProposalSizeGammaFactorDefinition          ,deProposalSizeGammaFactor          )
          call extractDataContent(deProposalSizeAccpetanceRateMinimumDefinition,deProposalSizeAccpetanceRateMinimum)
          call extractDataContent(deProposalSizeAccpetanceRateMaximumDefinition,deProposalSizeAccpetanceRateMaximum)
          call extractDataContent(deProposalSizeUpdateCountDefinition          ,deProposalSizeUpdateCount          )
          newDeProposalSize=deProposalSizeAdaptive(                                     &
               &                                   deProposalSizeGammaInitial         , &
               &                                   deProposalSizeGammaFactor          , &
               &                                   deProposalSizeAccpetanceRateMinimum, &
               &                                   deProposalSizeAccpetanceRateMaximum, &
               &                                   deProposalSizeUpdateCount            &
               &                                  )
       end select
     case default
       call Galacticus_Error_Report('deProposalSizeNew','deProposalSize type is unrecognized')
    end select
    return
  end function deProposalSizeNew

  ! Include all proposal methods.
  include 'constraints.differential_proposal_size.fixed.methods.inc'
  include 'constraints.differential_proposal_size.adaptive.methods.inc'

end module Constraints_Differential_Proposal_Size
