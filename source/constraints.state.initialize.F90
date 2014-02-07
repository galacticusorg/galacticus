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

!% Contains a module which implements state initialization for use when constraining \glc.

module Constraints_State_Initialize
  !% Implements state initialization for use when constraining \glc.
  use Constraints_Priors
  use Constraints_State
  use ISO_Varying_String
  use Pseudo_Random
  private
  public :: stateInitializorNew
  
  ! Define the basic state initialization class.
  type, abstract, public :: stateInitializor
   contains
     !@ <objectMethods>
     !@   <object>stateInitializor</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless class(state)\textgreater} simulationState\arginout, \textcolor{red}{\textless type(prior)\textgreater}[:] parameterPriors\arginout</arguments>
     !@     <description>Initialize the simulation state.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(stateInitializorInitialize), deferred :: initialize
  end type stateInitializor

  ! Interface for deferred functions.
  abstract interface
     subroutine stateInitializorInitialize(self,simulationState,parameterPriors)
       import :: stateInitializor, state, prior
       class(stateInitializor), intent(inout)               :: self
       class(state           ), intent(inout)               :: simulationState
       type (prior           ), intent(inout), dimension(:) :: parameterPriors
     end subroutine stateInitializorInitialize
  end interface

  ! Include all stateInitializor types.
  include 'constraints.state.initialize.prior_random.type.inc'
  include 'constraints.state.initialize.resume.type.inc'

contains

  function stateInitializorNew(definition,configFileName) result (newStateInitializor)
    !% Create a new {\tt stateInitializor} from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    use String_Handling
    implicit none
    class(stateInitializor), pointer                    :: newStateInitializor
    type (node            ), pointer    , intent(in   ) :: definition
    type (varying_string  ), optional   , intent(in   ) :: configFileName
    type (node            ), pointer                    :: stateInitializorFileDefinition
    type (varying_string  )                             :: logFileRoot

    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("priorRandom")
       allocate(stateInitializorPriorRandom :: newStateInitializor)
       select type (newStateInitializor)
       type is (stateInitializorPriorRandom)
          newStateInitializor=stateInitializorPriorRandom()
       end select
    case ("resume"     )
       allocate(stateInitializorResume :: newStateInitializor)
       select type (newStateInitializor)
       type is (stateInitializorResume)
          stateInitializorFileDefinition => XML_Get_First_Element_By_Tag_Name(definition,"logFileRoot")
          logFileRoot=getTextContent(stateInitializorFileDefinition)
          newStateInitializor=stateInitializorResume(char(logFileRoot))
       end select
    case default
       call Galacticus_Error_Report('stateInitializorNew','stateInitializor type is unrecognized')
    end select
    return
  end function stateInitializorNew

  ! Include all stateInitializor methods.
  include 'constraints.state.initialize.prior_random.methods.inc'
  include 'constraints.state.initialize.resume.methods.inc'

end module Constraints_State_Initialize
