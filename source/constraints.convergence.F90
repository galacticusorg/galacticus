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

!% Contains a module which implements convergence criteria for use when constraining \glc.

module Constraints_Convergence
  !% Implements convergence criteria for use when constraining \glc.
  use Constraints_State
  private
  public :: convergenceNew

  ! Define the basic convergence class.
  type, abstract, public :: convergence
   contains
     !@ <objectMethods>
     !@   <object>convergence</object>
     !@   <objectMethod>
     !@     <method>isConverged</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if the simulation is converged.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Reset the convergence object.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(convergenceIsConverged), deferred :: isConverged
     procedure(convergenceReset      ), deferred :: reset
  end type convergence

  ! Interface for deferred functions.
  abstract interface
     logical function convergenceIsConverged(self,simulationState,logLikelihood)
       import :: convergence, state
       class           (convergence), intent(inout)           :: self
       class           (state      ), intent(in   ), optional :: simulationState
       double precision             , intent(in   ), optional :: logLikelihood
     end function convergenceIsConverged
  end interface
  abstract interface
     subroutine convergenceReset(self)
       import :: convergence
       class(convergence), intent(inout) :: self
     end subroutine convergenceReset
  end interface

  ! Include all convergence types.
  include 'constraints.convergence.never.type.inc'
  include 'constraints.convergence.Gelman-Rubin.type.inc'

contains

  function convergenceNew(definition) result (newConvergence)
    !% Create a new convergence from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    class           (convergence), pointer                :: newConvergence
    type            (node       ), pointer, intent(in   ) :: definition
    type            (node       ), pointer                :: convergenceRhatDefinition               , convergenceBurnCountDefinition    , &
         &                                                   convergenceTestCountDefinition          , convergenceOutlierCountDefinition , &
         &                                                   convergenceOutlierSignificanceDefinition, convergenceOutlierOffsetDefinition
    integer                                               :: convergenceBurnCount                    , convergenceTestCount              , &
         &                                                   convergenceOutlierCount
    double precision                                      :: convergenceRhat                         , convergenceOutlierSignificance    , &
         &                                                   convergenceOutlierOffset

    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("GelmanRubin")
       allocate(convergenceGelmanRubin :: newConvergence)
       select type (newConvergence)
       type is (convergenceGelmanRubin)
          convergenceRhatDefinition                => XML_Get_First_Element_By_Tag_Name(definition,"Rhat"                      )
          convergenceBurnCountDefinition           => XML_Get_First_Element_By_Tag_Name(definition,"burnCount"                 )
          convergenceTestCountDefinition           => XML_Get_First_Element_By_Tag_Name(definition,"testCount"                 )
          convergenceOutlierCountDefinition        => XML_Get_First_Element_By_Tag_Name(definition,"outlierCountMaximum"       )
          convergenceOutlierSignificanceDefinition => XML_Get_First_Element_By_Tag_Name(definition,"outlierSignificance"       )
          convergenceOutlierOffsetDefinition       => XML_Get_First_Element_By_Tag_Name(definition,"outlierLogLikelihoodOffset")
          call extractDataContent(convergenceRhatDefinition               ,convergenceRhat               )
          call extractDataContent(convergenceBurnCountDefinition          ,convergenceBurnCount          )
          call extractDataContent(convergenceTestCountDefinition          ,convergenceTestCount          )
          call extractDataContent(convergenceOutlierCountDefinition       ,convergenceOutlierCount       )
          call extractDataContent(convergenceOutlierSignificanceDefinition,convergenceOutlierSignificance)
          call extractDataContent(convergenceOutlierOffsetDefinition      ,convergenceOutlierOffset)
          newConvergence=convergenceGelmanRubin(convergenceRhat,convergenceBurnCount,convergenceTestCount,convergenceOutlierCount,convergenceOutlierSignificance,convergenceOutlierOffset)
       end select
    case ("never")
       allocate(convergenceNever :: newConvergence)
       select type (newConvergence)
       type is (convergenceNever)
          newConvergence=convergenceNever()
       end select
    case default
       call Galacticus_Error_Report('convergenceNew','convergence type is unrecognized')
    end select
    return
  end function convergenceNew

  ! Include all convergence methods.
  include 'constraints.convergence.never.methods.inc'
  include 'constraints.convergence.Gelman-Rubin.methods.inc'

end module Constraints_Convergence
