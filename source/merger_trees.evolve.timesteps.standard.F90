!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !# <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepStandard">
  !#  <description>A merger tree evolution timestepping class which limits the step to the minimum of that given by the {\normalfont \ttfamily simple} and {\normalfont \ttfamily satellite} timesteps.</description>
  !#  <deepCopy>
  !#   <functionClass variables="simple, satellite"/>
  !#  </deepCopy>
  !#  <stateStorable>
  !#   <functionClass variables="simple, satellite"/>
  !#  </stateStorable>
  !# </mergerTreeEvolveTimestep>
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepStandard
     !% Implementation of a merger tree evolution timestepping class which limits the step to the minimum of that given by the
     !% {\normalfont \ttfamily simple} and {\normalfont \ttfamily satellite} timesteps.
     private
     type(mergerTreeEvolveTimestepSimple   ), pointer :: simple    => null()
     type(mergerTreeEvolveTimestepSatellite), pointer :: satellite => null()
   contains
     final     ::                 standardDestructor
     procedure :: timeEvolveTo => standardTimeEvolveTo
  end type mergerTreeEvolveTimestepStandard

  interface mergerTreeEvolveTimestepStandard
     !% Constructors for the {\normalfont \ttfamily standard} merger tree evolution timestep class.
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface mergerTreeEvolveTimestepStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily standard} merger tree evolution timestep class which takes a parameter set as
    !% input.
    use :: Input_Parameters   , only : inputParameter         , inputParameters
    use :: Cosmology_Functions, only : cosmologyFunctionsClass
    use :: Nodes_Operators    , only : nodeOperatorClass
    implicit none
    type (mergerTreeEvolveTimestepStandard)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class(nodeOperatorClass               ), pointer       :: nodeOperator_

    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !# <objectBuilder class="nodeOperator"       name="nodeOperator_"       source="parameters"/>
    self=mergerTreeEvolveTimestepStandard(cosmologyFunctions_,nodeOperator_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    !# <objectDestructor name="nodeOperator_"/>
    return
  end function standardConstructorParameters

  function standardConstructorInternal(cosmologyFunctions_,nodeOperator_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily standard} merger tree evolution timestep class.
    implicit none
    type (mergerTreeEvolveTimestepStandard)                        :: self
    class(cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    class(nodeOperatorClass               ), intent(in   ), target :: nodeOperator_

    allocate(self%simple   )
    allocate(self%satellite)
    !# <referenceConstruct isResult="yes" owner="self" object="simple"    constructor="mergerTreeEvolveTimestepSimple   (timeStepAbsolute         =1.0d+0,timeStepRelative         =1.0d-1,cosmologyFunctions_=cosmologyFunctions_)"/>
    !# <referenceConstruct isResult="yes" owner="self" object="satellite" constructor="mergerTreeEvolveTimestepSatellite(timeOffsetMaximumAbsolute=1.0d-2,timeOffsetMaximumRelative=1.0d-3,nodeOperator_      =nodeOperator_      )"/>
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !% Destructor for the {\normalfont \ttfamily standard} merger tree evolution timestep class.
    implicit none
    type(mergerTreeEvolveTimestepStandard), intent(inout) :: self

    !# <objectDestructor name="self%simple"   />
    !# <objectDestructor name="self%satellite"/>
    return
  end subroutine standardDestructor

  double precision function standardTimeEvolveTo(self,node,task,taskSelf,report,lockNode,lockType)
    !% Determine a suitable timestep for {\normalfont \ttfamily node} by combining the {\normalfont \ttfamily simple} and
    !% {\normalfont \ttfamily satellite} timesteps.
    use :: ISO_Varying_String, only : varying_string
    implicit none
    class           (mergerTreeEvolveTimestepStandard), intent(inout), target            :: self
    type            (treeNode                        ), intent(inout), target            :: node
    procedure       (timestepTask                    ), intent(  out), pointer           :: task
    class           (*                               ), intent(  out), pointer           :: taskSelf
    logical                                           , intent(in   )                    :: report
    type            (treeNode                        ), intent(  out), pointer, optional :: lockNode
    type            (varying_string                  ), intent(  out)         , optional :: lockType
    double precision                                                                     :: timeEvolveToSimple, timeEvolveToSatellite
    procedure       (timestepTask                    )               , pointer           :: taskSimple        , taskSatellite
    type            (treeNode                        )               , pointer           :: lockNodeSimple    , lockNodeSatellite
    class           (*                               )               , pointer           :: taskSelfSimple    , taskSelfSatellite
    type            (varying_string                  ), save                             :: lockTypeSimple    , lockTypeSatellite
    !$omp threadprivate(lockTypeSimple,lockTypeSatellite)

    timeEvolveToSimple   =self%simple   %timeEvolveTo(node,taskSimple   ,taskSelfSimple   ,report,lockNode,lockType)
    if (present(lockNode)) lockNodeSimple    => lockNode
    if (present(lockType)) lockTypeSimple    =  lockType
    timeEvolveToSatellite=self%satellite%timeEvolveTo(node,taskSatellite,taskSelfSatellite,report,lockNode,lockType)
    if (present(lockNode)) lockNodeSatellite => lockNode
    if (present(lockType)) lockTypeSatellite =  lockType
    if (timeEvolveToSatellite <= timeEvolveToSimple) then
       standardTimeEvolveTo            =  timeEvolveToSatellite
       task                            => taskSatellite
       taskSelf                        => taskSelfSatellite
       if (present(lockNode)) lockNode => lockNodeSatellite
       if (present(lockType)) lockType =  lockTypeSatellite
    else
       standardTimeEvolveTo            =  timeEvolveToSimple
       task                            => taskSimple
       taskSelf                        => taskSelfSimple
       if (present(lockNode)) lockNode => lockNodeSimple
       if (present(lockType)) lockType =  lockTypeSimple
    end if
    return
  end function standardTimeEvolveTo
