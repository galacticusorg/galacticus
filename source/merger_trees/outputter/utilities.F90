!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{RST
Contains a module of globally-accessible functions supporting the :galacticus-class:`mergerTreeOutputterClass` class.
!!}

module Merger_Tree_Outputters_Utilities
  !!{RST
  Provides globally-accessible functions supporting the :galacticus-class:`mergerTreeOutputterClass` class.

  These exist so that a class may own a :galacticus-class:`mergerTreeOutputterClass` object without its module having to
  ``use`` the ``Merger_Tree_Outputters`` module. That matters because ``Merger_Tree_Outputters`` itself (transitively, via
  :galacticus-class:`mergerTreeOutputterFullState`) depends on ``Merger_Tree_Construction``, which in turn depends on
  ``Nodes_Operators``. A :galacticus-class:`nodeOperatorClass` which named the outputter class directly would therefore close a
  circular dependency between those modules; routing the accesses through these globally-callable functions avoids that.
  !!}
  private
  public :: mergerTreeOutputterConstruct       , mergerTreeOutputterDestruct        , &
       &    mergerTreeOutputterDeepCopy        , mergerTreeOutputterDeepCopyReset   , &
       &    mergerTreeOutputterDeepCopyFinalize, mergerTreeOutputterOutputTrajectory

  ! Module-scope pointer to our outputter object. This is used for reference counting so that debugging information is consistent
  ! between the increments and decrements.
  class(*), pointer :: mergerTreeOutputter__
  !$omp threadprivate(mergerTreeOutputter__)

contains

  !![
  <functionGlobal>
   <unitName>mergerTreeOutputterConstruct</unitName>
   <type>void</type>
   <module>Input_Parameters, only : inputParameters</module>
   <arguments>type (inputParameters), intent(inout), target  :: parameters          </arguments>
   <arguments>class(*              ), intent(  out), pointer :: mergerTreeOutputter_</arguments>
  </functionGlobal>
  !!]
  subroutine mergerTreeOutputterConstruct(parameters,mergerTreeOutputter_)
    !!{RST
    Build a ``mergerTreeOutputter`` object from a given parameter set. This is a globally-callable function to allow us to subvert
    the class/module hierarchy.

    Note that, unlike some similar functions, we do *not* search parent parameter sets for a ``mergerTreeOutputter`` definition.
    Doing so would cause a caller which failed to specify its own outputter to silently acquire the outputter used for regular
    outputs---and so to write its data into the regular ``Outputs`` group. Instead, if no outputter is specified, a default
    outputter writing to the ``nodeTrajectories`` group is constructed.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Parameters  , only : inputParameter          , inputParameters
    use :: ISO_Varying_String, only : varying_string          , var_str
    use :: Merger_Tree_Outputters, only : mergerTreeOutputterClass, mergerTreeOutputter
    implicit none
    type (inputParameters), intent(inout), target  :: parameters
    class(*              ), intent(  out), pointer :: mergerTreeOutputter_
    type (inputParameters)                         :: parametersDefault
    type (varying_string ), dimension(1)           :: allowedParameterName__
    type (varying_string )                         :: defaultXML__

    if (parameters%isPresent('mergerTreeOutputter')) then
       mergerTreeOutputter__ => mergerTreeOutputter(parameters       )
    else
       defaultXML__             =var_str('<parameters><mergerTreeOutputter value="standard"><outputsGroupName value="nodeTrajectories"/><outputReferences value="false"/><nodePropertyExtractor value="multi"><nodePropertyExtractor value="time"/><nodePropertyExtractor value="indicesTree"/><nodePropertyExtractor value="nodeIndices"/><nodePropertyExtractor value="trajectoryEvent"/></nodePropertyExtractor></mergerTreeOutputter></parameters>')
       allowedParameterName__(1)=var_str('mergerTreeOutputter')
       parametersDefault        =inputParameters(defaultXML__,allowedParameterNames=allowedParameterName__,noOutput=.true.)
       mergerTreeOutputter__ => mergerTreeOutputter(parametersDefault)
    end if
    select type (mergerTreeOutputter__)
    class is (mergerTreeOutputterClass)
       !![
       <referenceCountIncrement object="mergerTreeOutputter__"/>
       !!]
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    mergerTreeOutputter_ => mergerTreeOutputter__
    return
  end subroutine mergerTreeOutputterConstruct

  !![
  <functionGlobal>
   <unitName>mergerTreeOutputterOutputTrajectory</unitName>
   <type>void</type>
   <module>Galacticus_Nodes, only : treeNode</module>
   <module>ISO_C_Binding, only : c_size_t</module>
   <arguments>class  (*       ), intent(inout) :: mergerTreeOutputter_</arguments>
   <arguments>type   (treeNode), intent(inout) :: node                </arguments>
   <arguments>integer(c_size_t), intent(in   ) :: indexOutput         </arguments>
  </functionGlobal>
  !!]
  subroutine mergerTreeOutputterOutputTrajectory(mergerTreeOutputter_,node,indexOutput)
    !!{RST
    Output a single node as a trajectory record, using a ``mergerTreeOutputter`` object passed to us as an unlimited polymorphic
    object. The output type is fixed to ``trajectory`` here so that callers need not name the ``outputGroupType`` enumeration
    (which would reintroduce the module dependency that these functions exist to avoid).
    !!}
    use, intrinsic :: ISO_C_Binding         , only : c_size_t
    use            :: Error                 , only : Error_Report
    use            :: Galacticus_Nodes      , only : treeNode
    use            :: Merger_Tree_Outputters, only : mergerTreeOutputterClass, outputGroupTypeTrajectory
    implicit none
    class  (*       ), intent(inout) :: mergerTreeOutputter_
    type   (treeNode), intent(inout) :: node
    integer(c_size_t), intent(in   ) :: indexOutput

    select type (mergerTreeOutputter_)
    class is (mergerTreeOutputterClass)
       call mergerTreeOutputter_%outputNode(node,indexOutput,outputGroupTypeTrajectory)
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine mergerTreeOutputterOutputTrajectory

  !![
  <functionGlobal>
   <unitName>mergerTreeOutputterDestruct</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout), pointer :: mergerTreeOutputter_</arguments>
  </functionGlobal>
  !!]
  subroutine mergerTreeOutputterDestruct(mergerTreeOutputter_)
    !!{RST
    Destruct a ``mergerTreeOutputter`` object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error                 , only : Error_Report
    use :: Merger_Tree_Outputters, only : mergerTreeOutputterClass
    implicit none
    class(*), intent(inout), pointer :: mergerTreeOutputter_

    mergerTreeOutputter__ => mergerTreeOutputter_
    select type (mergerTreeOutputter__)
    class is (mergerTreeOutputterClass)
       !![
       <objectDestructor name="mergerTreeOutputter__"/>
       !!]
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine mergerTreeOutputterDestruct

  !![
  <functionGlobal>
   <unitName>mergerTreeOutputterDeepCopyReset</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine mergerTreeOutputterDeepCopyReset(self)
    !!{RST
    Perform a deep copy reset of a ``mergerTreeOutputter`` object.
    !!}
    use :: Error                 , only : Error_Report
    use :: Merger_Tree_Outputters, only : mergerTreeOutputterClass
    implicit none
    class(*), intent(inout) :: self

    select type (self)
    class is (mergerTreeOutputterClass)
       call self%deepCopyReset()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine mergerTreeOutputterDeepCopyReset

  !![
  <functionGlobal>
   <unitName>mergerTreeOutputterDeepCopyFinalize</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine mergerTreeOutputterDeepCopyFinalize(self)
    !!{RST
    Finalize a deep copy of a ``mergerTreeOutputter`` object.
    !!}
    use :: Error                 , only : Error_Report
    use :: Merger_Tree_Outputters, only : mergerTreeOutputterClass
    implicit none
    class(*), intent(inout) :: self

    select type (self)
    class is (mergerTreeOutputterClass)
       call self%deepCopyFinalize()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine mergerTreeOutputterDeepCopyFinalize

  !![
  <functionGlobal>
   <unitName>mergerTreeOutputterDeepCopy</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self, destination</arguments>
  </functionGlobal>
  !!]
  subroutine mergerTreeOutputterDeepCopy(self,destination)
    !!{RST
    Perform a deep copy of a ``mergerTreeOutputter`` object.
    !!}
    use :: Error                 , only : Error_Report
    use :: Merger_Tree_Outputters, only : mergerTreeOutputterClass
    implicit none
    class(*), intent(inout) :: self, destination

    select type (self)
    class is (mergerTreeOutputterClass)
       select type (destination)
       class is (mergerTreeOutputterClass)
          call self%deepCopy(destination)
       class default
          call Error_Report("unexpected class"//{introspection:location})
       end select
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine mergerTreeOutputterDeepCopy

end module Merger_Tree_Outputters_Utilities
