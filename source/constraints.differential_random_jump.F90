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

!% Contains a module which implements algorithms for the random jump component in differential evolution algorithms.

module Constraints_Differential_Random_Jump
  !% Implements algorithms for the random jump component in differential evolution algorithms.
  use Statistics_Distributions
  use Constraints_State
  private
  public :: deRandomJumpNew

  ! Define the basic random jump class.
  type, abstract, public :: deRandomJump
   contains
     !@ <objectMethods>
     !@   <object>deRandomJump</object>
     !@   <objectMethod>
     !@     <method>sample</method>
     !@     <type>\doubleone</type>
     !@     <arguments></arguments>
     !@     <description>Return a random jump sampeld from the appropriate distribution.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(deRandomJumpSample), deferred :: sample
  end type deRandomJump

  ! Interface for deferred functions.
  abstract interface
     function deRandomJumpSample(self,randomDistributions,simulationState)
       import :: deRandomJump, distributionList, state
       class           (deRandomJump    )                                      , intent(inout) :: self
       type            (distributionList), dimension(:)                        , intent(in   ) :: randomDistributions
       class           (state           )                                      , intent(in   ) :: simulationState
      double precision                   , dimension(size(randomDistributions))                :: deRandomJumpSample
     end function deRandomJumpSample
  end interface
  
  ! Include all proposal size types.
  include 'constraints.differential_random_jump.simple.type.inc'
  include 'constraints.differential_random_jump.adaptive.type.inc'

contains

  function deRandomJumpNew(definition) result (newDeRandomJump)
    !% Create a new differential evolution random jump from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    class  (deRandomJump), pointer                :: newDeRandomJump
    type   (node        ), pointer, intent(in   ) :: definition
   
    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("simple")
       allocate(deRandomJumpSimple :: newDeRandomJump)
       select type (newDeRandomJump)
       type is (deRandomJumpSimple)
          newDeRandomJump=deRandomJumpSimple()
       end select
    case ("adaptive")
       allocate(deRandomJumpAdaptive :: newDeRandomJump)
       select type (newDeRandomJump)
       type is (deRandomJumpAdaptive)
          newDeRandomJump=deRandomJumpAdaptive()
       end select
     case default
       call Galacticus_Error_Report('deRandomJumpNew','deRandomJump type is unrecognized')
    end select
    return
  end function deRandomJumpNew

  ! Include all proposal methods.
  include 'constraints.differential_random_jump.simple.methods.inc'
  include 'constraints.differential_random_jump.adaptive.methods.inc'

end module Constraints_Differential_Random_Jump
