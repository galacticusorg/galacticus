!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

!!{
Contains a module of globally-accessible functions supporting the {\normalfont \ttfamily accretionFlows} class.
!!}

module Accretion_Flows_Utilities
  !!{
  Provides globally-accessible functions supporting the {\normalfont \ttfamily accretionFlows} class.
  !!}
  private
  public :: accretionFlowsConstruct    , accretionFlowsDensity         , accretionFlowsDestruct, accretionFlowsDeepCopy, &
       &    accretionFlowsDeepCopyReset, accretionFlowsDeepCopyFinalize

  ! Module-scope pointer to our accretionFlow object. This is used for reference counting so that debugging information is consistent
  ! between the increments and decrements.
  class(*), pointer :: accretionFlows__
  !$omp threadprivate(accretionFlows__)

contains

  !![
  <functionGlobal>
   <unitName>accretionFlowsConstruct</unitName>
   <type>void</type>
   <module>Input_Parameters, only : inputParameters</module>
   <arguments>type (inputParameters), intent(inout), target  :: parameters</arguments>
   <arguments>class(*              ), intent(  out), pointer :: accretionFlows_</arguments>
  </functionGlobal>
  !!]
  subroutine accretionFlowsConstruct(parameters,accretionFlows_)
    !!{
    Build a {\normalfont \ttfamily accretionFlowsClass} object from a given parameter set. This is a globally-callable function
    to allow us to subvert the class/module hierarchy.
    !!}
    use :: Error                             , only : Error_Report
    use :: Input_Parameters                  , only : inputParameter, inputParameters
    use :: Spherical_Collapse_Accretion_Flows, only : accretionFlows, accretionFlowsClass
    implicit none
    type (inputParameters), intent(inout), target  :: parameters
    class(*              ), intent(  out), pointer :: accretionFlows_
    type (inputParameters)               , pointer :: parametersCurrent

    parametersCurrent => parameters
    do while (.not.parametersCurrent%isPresent('accretionFlows').and.associated(parametersCurrent%parent))
       parametersCurrent => parametersCurrent%parent
    end do
    if (.not.parametersCurrent%isPresent('accretionFlows')) parametersCurrent => parameters
    accretionFlows__ => accretionFlows(parametersCurrent)
    select type (accretionFlows__)
    class is (accretionFlowsClass)
       !![
       <referenceCountIncrement object="accretionFlows__"/>
       !!]
       call accretionFlows__%autoHook()
    class default
       call Error_Report('accretionFlow must be of the "accretionFlowsClass" class'//{introspection:location})
    end select
    accretionFlows_ => accretionFlows__
    return
  end subroutine accretionFlowsConstruct

  !![
  <functionGlobal>
   <unitName>accretionFlowsDensity</unitName>
   <type>double precision</type>
   <module>Galacticus_Nodes, only : treeNode</module>
   <arguments>class           (*       ), intent(inout) :: accretionFlows_</arguments>
   <arguments>type            (treeNode), intent(inout) :: node</arguments>
   <arguments>double precision          , intent(in   ) :: radius</arguments>
  </functionGlobal>
  !!]
  double precision function accretionFlowsDensity(accretionFlows_,node,radius)
    !!{
    Perform the accretionFlow for a {\normalfont \ttfamily accretionFlowsClass} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error                             , only : Error_Report
    use :: Galacticus_Nodes                  , only : treeNode
    use :: Spherical_Collapse_Accretion_Flows, only : accretionFlowsClass
    implicit none
    class           (*       ), intent(inout) :: accretionFlows_
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: radius

    select type (accretionFlows_)
    class is (accretionFlowsClass)
       accretionFlowsDensity=accretionFlows_%density(node,radius)
    class default
       accretionFlowsDensity=0.0d0
       call Error_Report('accretionFlow must be of the "accretionFlowsClass" class'//{introspection:location})
    end select
    return
  end function accretionFlowsDensity

  !![
  <functionGlobal>
   <unitName>accretionFlowsDestruct</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout), pointer :: accretionFlows_</arguments>
  </functionGlobal>
  !!]
  subroutine accretionFlowsDestruct(accretionFlows_)
    !!{
    Destruct a {\normalfont \ttfamily accretionFlowsClass} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error                             , only : Error_Report
    use :: Spherical_Collapse_Accretion_Flows, only : accretionFlowsClass
    use :: Function_Classes
    use :: iso_varying_string
    implicit none
    class(*), intent(inout), pointer :: accretionFlows_

    accretionFlows__ => accretionFlows_
    select type (accretionFlows__)
    class is (accretionFlowsClass)
       !![
       <objectDestructor name="accretionFlows__"/>
       !!]
    class default
       call Error_Report('accretionFlow must be of the "accretionFlowsClass" class'//{introspection:location})
    end select
    return
  end subroutine accretionFlowsDestruct

  !![
  <functionGlobal>
   <unitName>accretionFlowsDeepCopyReset</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine accretionFlowsDeepCopyReset(self)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use :: Error                             , only : Error_Report
    use :: Spherical_Collapse_Accretion_Flows, only : accretionFlowsClass
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (accretionFlowsClass)
       call self%deepCopyReset()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine accretionFlowsDeepCopyReset
  
  !![
  <functionGlobal>
   <unitName>accretionFlowsDeepCopyFinalize</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine accretionFlowsDeepCopyFinalize(self)
    !!{
    Finalize a deep copy of galactic structure objects.
    !!}
    use :: Error                             , only : Error_Report
    use :: Spherical_Collapse_Accretion_Flows, only : accretionFlowsClass
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (accretionFlowsClass)
       call self%deepCopyFinalize()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine accretionFlowsDeepCopyFinalize
  
  !![
  <functionGlobal>
   <unitName>accretionFlowsDeepCopy</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self, destination</arguments>
  </functionGlobal>
  !!]
  subroutine accretionFlowsDeepCopy(self,destination)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use :: Error                             , only : Error_Report
    use :: Spherical_Collapse_Accretion_Flows, only : accretionFlowsClass
    implicit none
    class(*), intent(inout) :: self, destination

    select type (self)
    class is (accretionFlowsClass)
       select type (destination)
       class is (accretionFlowsClass)
          call self%deepCopy(destination)
       class default
          call Error_Report("unexpected class"//{introspection:location})
       end select
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine accretionFlowsDeepCopy

end module Accretion_Flows_Utilities
