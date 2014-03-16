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

!% Contains a module which implements algorithms for the temperature exponent of proposal size in tempered differential evolution
!% algorithms.

module Constraints_Differential_Prop_Size_Temp_Exp
  !% Implements algorithms for the temperature exponent of proposal size in tempered differential evolution
  !% algorithms.
  use Constraints_State
  use Constraints_Convergence
  private
  public :: dePropSizeTempExpNew

  ! Define the basic proposal size temperature exponent class.
  type, abstract, public :: dePropSizeTempExp
   contains
     !@ <objectMethods>
     !@   <object>dePropSizeTempExp</object>
     !@   <objectMethod>
     !@     <method>exponent</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless class(state)[:]\textgreater} temperedStates\argin, \doubleone temperatures\argin, \textcolor{red}{\textless class(convergence)\textgreater} simulationConvergence\arginout</arguments>
     !@     <description>Return the temperature exponent of proposal step size for tempered differential evolution simulations.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(dePropSizeTempExpExponent), deferred :: exponent
  end type dePropSizeTempExp

  ! Interface for deferred functions.
  abstract interface
     double precision function dePropSizeTempExpExponent(self,temperedStates,temperatures,simulationState,simulationConvergence)
       import :: dePropSizeTempExp, state, convergence
       class           (dePropSizeTempExp), intent(inout)               :: self
       class           (state            ), intent(in   )               :: simulationState
       class           (state            ), intent(in   ), dimension(:) :: temperedStates
       double precision                   , intent(in   ), dimension(:) :: temperatures
       class           (convergence      ), intent(inout)               :: simulationConvergence
     end function dePropSizeTempExpExponent
  end interface

  ! Include all proposal size types.
  include 'constraints.differential_proposal_size.temperature_exponent.fixed.type.inc'
  include 'constraints.differential_proposal_size.temperature_exponent.adaptive.type.inc'

contains

  function dePropSizeTempExpNew(definition) result (newDePropSizeTempExp)
    !% Create a new differential evolution proposal size temperature exponent object from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    class          (dePropSizeTempExp), pointer                :: newDePropSizeTempExp
    type           (node             ), pointer, intent(in   ) :: definition
    type           (node             ), pointer                :: dePropSizeTempExpExponentInitialDefinition, dePropSizeTempExpExponentFactorDefinition , &
         &                                                        dePropSizeTempExpGradientMinimumDefinition, dePropSizeTempExpGradientMaximumDefinition, &
         &                                                        dePropSizeTempExpUpdateCountDefinition    , dePropSizeTempExpExponentMinimumDefinition, &
         &                                                        dePropSizeTempExpExponentMaximumDefinition
    double precision                                           :: dePropSizeTempExpExponentInitial          , dePropSizeTempExpExponentFactor           , &
         &                                                        dePropSizeTempExpGradientMinimum          , dePropSizeTempExpGradientMaximum          , &
         &                                                        dePropSizeTempExpExponentMinimum          , dePropSizeTempExpExponentMaximum
    integer                                                    :: dePropSizeTempExpUpdateCount
   
    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("fixed")
       allocate(dePropSizeTempExpFixed :: newDePropSizeTempExp)
       select type (newDePropSizeTempExp)
       type is (dePropSizeTempExpFixed)
          dePropSizeTempExpExponentInitialDefinition => XML_Get_First_Element_By_Tag_Name(definition,"exponentInitial")
          call extractDataContent(dePropSizeTempExpExponentInitialDefinition,dePropSizeTempExpExponentInitial)
          newDePropSizeTempExp=dePropSizeTempExpFixed(dePropSizeTempExpExponentInitial)
       end select
    case ("adaptive")
       allocate(dePropSizeTempExpAdaptive :: newDePropSizeTempExp)
       select type (newDePropSizeTempExp)
       type is (dePropSizeTempExpAdaptive)
          dePropSizeTempExpExponentInitialDefinition => XML_Get_First_Element_By_Tag_Name(definition,"exponentInitial")
          dePropSizeTempExpExponentMinimumDefinition => XML_Get_First_Element_By_Tag_Name(definition,"exponentMinimum")
          dePropSizeTempExpExponentMaximumDefinition => XML_Get_First_Element_By_Tag_Name(definition,"exponentMaximum")
          dePropSizeTempExpExponentFactorDefinition  => XML_Get_First_Element_By_Tag_Name(definition,"exponentFactor" )
          dePropSizeTempExpGradientMinimumDefinition => XML_Get_First_Element_By_Tag_Name(definition,"gradientMinimum")
          dePropSizeTempExpGradientMaximumDefinition => XML_Get_First_Element_By_Tag_Name(definition,"gradientMaximum")
          dePropSizeTempExpUpdateCountDefinition     => XML_Get_First_Element_By_Tag_Name(definition,"updateCount"    )
          call extractDataContent(dePropSizeTempExpExponentInitialDefinition,dePropSizeTempExpExponentInitial)
          call extractDataContent(dePropSizeTempExpExponentMinimumDefinition,dePropSizeTempExpExponentMinimum)
          call extractDataContent(dePropSizeTempExpExponentMaximumDefinition,dePropSizeTempExpExponentMaximum)
          call extractDataContent(dePropSizeTempExpExponentFactorDefinition ,dePropSizeTempExpExponentFactor )
          call extractDataContent(dePropSizeTempExpGradientMinimumDefinition,dePropSizeTempExpGradientMinimum)
          call extractDataContent(dePropSizeTempExpGradientMaximumDefinition,dePropSizeTempExpGradientMaximum)
          call extractDataContent(dePropSizeTempExpUpdateCountDefinition    ,dePropSizeTempExpUpdateCount    )
          newDePropSizeTempExp=dePropSizeTempExpAdaptive(                                  &
               &                                         dePropSizeTempExpExponentInitial, &
               &                                         dePropSizeTempExpExponentMinimum, &
               &                                         dePropSizeTempExpExponentMaximum, &
               &                                         dePropSizeTempExpExponentFactor , &
               &                                         dePropSizeTempExpGradientMinimum, &
               &                                         dePropSizeTempExpGradientMaximum, &
               &                                         dePropSizeTempExpUpdateCount      &
               &                                        )
       end select
     case default
       call Galacticus_Error_Report('dePropSizeTempExpNew','dePropSizeTempExp type is unrecognized')
    end select
    return
  end function dePropSizeTempExpNew

  ! Include all proposal size temperature exponent methods.
  include 'constraints.differential_proposal_size.temperature_exponent.fixed.methods.inc'
  include 'constraints.differential_proposal_size.temperature_exponent.adaptive.methods.inc'

end module Constraints_Differential_Prop_Size_Temp_Exp
