!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  !![
  <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepSimple">
   <description>
    A merger tree evolution timestepping class enforces that
    \begin{eqnarray}
    \Delta t &amp;\le&amp; t_\mathrm{simple}, \\
    \Delta t &amp;\le&amp; \epsilon_\mathrm{simple} (a/\dot{a}),
    \end{eqnarray}
    where $t_\mathrm{simple}=${\normalfont \ttfamily [timestepSimpleAbsolute]}, $\epsilon_\mathrm{simple}=${\normalfont \ttfamily
    [timestepSimpleRelative]}, and $a$ is expansion factor. These criteria are intended to prevent any one node evolving over an
    excessively large time in one step. In general, these criteria are not necessary, as nodes should be free to evolve as far as
    possible unless prevented by some physical requirement. These criteria are therefore present to provide a simple example of how
    timestep criteria work.
   </description>
  </mergerTreeEvolveTimestep>
  !!]
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepSimple
     !!{
     Implementation of an output times class which reads a simple of output times from a parameter.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: timeStepAbsolute             , timeStepRelative
   contains
     final     ::                 simpleDestructor
     procedure :: timeEvolveTo => simpleTimeEvolveTo
  end type mergerTreeEvolveTimestepSimple

  interface mergerTreeEvolveTimestepSimple
     !!{
     Constructors for the \refClass{mergerTreeEvolveTimestepSimple} merger tree evolution timestep class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface mergerTreeEvolveTimestepSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepSimple} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeEvolveTimestepSimple)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    double precision                                                :: timeStepAbsolute   , timeStepRelative

    !![
    <inputParameter>
      <name>timeStepRelative</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The maximum allowed relative change in time for a single step in the evolution of a node.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeStepAbsolute</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The maximum allowed absolute change in time (in Gyr) for a single step in the evolution of a node.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=mergerTreeEvolveTimestepSimple(timeStepAbsolute,timeStepRelative,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(timeStepAbsolute,timeStepRelative,cosmologyFunctions_) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepSimple} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    implicit none
    type            (mergerTreeEvolveTimestepSimple)                        :: self
    double precision                                , intent(in   )         :: timeStepAbsolute   , timeStepRelative
    class           (cosmologyFunctionsClass       ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="timeStepAbsolute, timeStepRelative, *cosmologyFunctions_"/>
    !!]

    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolveTimestepSimple} merger tree evolution timestep class.
    !!}
    implicit none
    type(mergerTreeEvolveTimestepSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine simpleDestructor

  double precision function simpleTimeEvolveTo(self,timeEnd,node,task,taskSelf,report,lockNode,lockType)
    !!{
    Determine a suitable timestep for {\normalfont \ttfamily node} using the simple method. This simply selects the smaller of {\normalfont \ttfamily
    timeStepAbsolute} and {\normalfont \ttfamily timeStepRelative}$H^{-1}(t)$.
    !!}
    use :: Evolve_To_Time_Reports, only : Evolve_To_Time_Report
    use :: Galacticus_Nodes      , only : nodeComponentBasic   , treeNode
    use :: ISO_Varying_String    , only : varying_string
    implicit none
    class           (mergerTreeEvolveTimestepSimple), intent(inout), target            :: self
    double precision                                , intent(in   )                    :: timeEnd
    type            (treeNode                      ), intent(inout), target            :: node
    procedure       (timestepTask                  ), intent(  out), pointer           :: task
    class           (*                             ), intent(  out), pointer           :: taskSelf
    logical                                         , intent(in   )                    :: report
    type            (treeNode                      ), intent(  out), pointer, optional :: lockNode
    type            (varying_string                ), intent(  out)         , optional :: lockType
    class           (nodeComponentBasic            )               , pointer           :: basic
    double precision                                                                   :: expansionFactor, timescaleExpansion, &
         &                                                                                time
    !$GLC attributes unused :: timeEnd

    ! Find current expansion timescale.
    basic => node%basic()
    if (self%timeStepRelative > 0.0d0) then
       time               =        basic                    %time           (               )
       expansionFactor    =        self %cosmologyFunctions_%expansionFactor(           time)
       timescaleExpansion =  1.0d0/self %cosmologyFunctions_%expansionRate  (expansionFactor)
       simpleTimeEvolveTo =  min(self%timestepRelative*timescaleExpansion,self%timeStepAbsolute)+basic%time()
    else
       simpleTimeEvolveTo =                                               self%timeStepAbsolute +basic%time()
    end if
    task                            => null()
    taskSelf                        => null()
    if (present(lockNode)) lockNode => node
    if (present(lockType)) lockType =  "simple"
    if (        report   ) call Evolve_To_Time_Report("simple: ",simpleTimeEvolveTo)
    return
  end function simpleTimeEvolveTo
