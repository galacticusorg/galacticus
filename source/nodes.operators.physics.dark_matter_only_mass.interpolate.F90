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

  !!{
  Implements a node operator class that interpolates the ``\gls{dmou}'' mass of the halo linearly between child and parent nodes.
  !!}

  !![
  <nodeOperator name="nodeOperatorDMOInterpolate">
   <description>
    A node operator class that interpolates the ``\gls{dmou}'' mass of the halo linearly between child and parent nodes.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDMOInterpolate
     !!{
     A node operator class that causes dark matter profile scale radius to be interpolated linearly between child and parent nodes.
     !!}
     private
     integer :: massDMOTargetID
   contains
     procedure :: nodeInitialize                      => dmoInterpolateNodeInitialize
     procedure :: nodePromote                         => dmoInterpolateNodePromote
     procedure :: nodesMerge                          => dmoInterpolateNodesMerge
     procedure :: differentialEvolutionAnalytics      => dmoInterpolateDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionSolveAnalytics => dmoInterpolateDifferentialEvolutionSolveAnalytics
  end type nodeOperatorDMOInterpolate
  
  interface nodeOperatorDMOInterpolate
     !!{
     Constructors for the \refClass{nodeOperatorDMOInterpolate} node operator class.
     !!}
     module procedure dmoInterpolateConstructorParameters
     module procedure dmoInterpolateConstructorInternal
  end interface nodeOperatorDMOInterpolate
  
contains
  
  function dmoInterpolateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDMOInterpolate} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDMOInterpolate)                :: self
    type (inputParameters           ), intent(inout) :: parameters
 
    self=nodeOperatorDMOInterpolate()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function dmoInterpolateConstructorParameters

  function dmoInterpolateConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorDMOInterpolate} node operator class.
    !!}
    implicit none
    type(nodeOperatorDMOInterpolate) :: self

    !![
    <addMetaProperty component="basic" name="massDMOTarget" id="self%massDMOTargetID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function dmoInterpolateConstructorInternal

  recursive subroutine dmoInterpolateNodeInitialize(self,node)
    !!{
    Compute the rate of growth of the ``\gls{dmou}'' mass of a halo assuming a constant growth rate.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorDMOInterpolate), intent(inout), target  :: self
    type            (treeNode                  ), intent(inout), target  :: node
    type            (treeNode                  )               , pointer :: nodeChild          , nodeParent
    class           (nodeComponentBasic        )               , pointer :: basicChild         , basicParent     , &
         &                                                                  basic
    double precision                                                     :: timeDelta          , massUnresolved  , &
         &                                                                  massTotalProgenitor, timeLastIsolated
    !$GLC attributes unused :: self

    basic => node%basic()
    ! Set the last isolated time to the current time at the farthest point along the future of this branch.
    if (.not.associated(node%firstChild)) then
       nodeParent => node
       do while (associated(nodeParent%parent).and.nodeParent%isPrimaryProgenitor())
          nodeParent => nodeParent%parent
       end do
       basicParent      => nodeParent %basic()
       timeLastIsolated =  basicParent%time ()
       ! Set the last isolated time for all nodes along this branch - this avoids having to re-walk this branch for each node
       ! along it.
       nodeParent => node
       do while (associated(nodeParent%parent).and.nodeParent%isPrimaryProgenitor())
          basicParent => nodeParent%basic()
          call basicParent%timeLastIsolatedSet(timeLastIsolated)
          nodeParent => nodeParent%parent
       end do
    end if
    ! Determine node status.
    if (node%isSatellite()) then
       ! Node is a satellite - we assume no accretion.
       call basic%         accretionRateSet(                     0.0d0       )
       call basic%floatRank0MetaPropertySet(self%massDMOTargetID,basic%mass())
    else if (.not.associated(node%parent)) then
       ! For parent-less nodes (i.e. the root node of the tree), the rate is set equal to that of the
       ! progenitor, if it has one.
       nodeChild => node%firstChild
       if (associated(nodeChild)) then
          ! Get the basic component of the child node.
          basicChild => nodeChild%basic()
          ! Ensure the child has a mass growth rate computed.
          call self%nodeInitialize(nodeChild)
          ! Get the growth rate of the child.
          call basic%         accretionRateSet(                     basicChild%accretionRate())
          call basic%floatRank0MetaPropertySet(self%massDMOTargetID,basic     %mass         ())
       else
          ! Parentless node has no child - set a zero growth rate.
          call basic%         accretionRateSet(                     0.0d0       )
          call basic%floatRank0MetaPropertySet(self%massDMOTargetID,basic%mass())
       end if
    else
       ! Get the parent node.
       nodeParent => node%parent
       ! Get the basic component of the parent node.
       basicParent => nodeParent%basic()
       ! Compute the unresolved mass.
       massUnresolved=nodeMassUnresolved(nodeParent)
       if (massUnresolved > 0.0d0) then
          ! Positive mass growth - assume this occurs entirely in the main progenitor.
          if (node%isPrimaryProgenitor()) then
             ! Main progenitor - compute required growth rate.
             timeDelta=basicParent%time()-basic%time()
             if (timeDelta > 0.0d0) &
                  & call basic%         accretionRateSet(                     massUnresolved/timeDelta   )
             call        basic%floatRank0MetaPropertySet(self%massDMOTargetID,basic%mass()+massUnresolved)
          else
             ! Non-main progenitor - assume zero growth rate.
             call basic%         accretionRateSet(                     0.0d0       )
             call basic%floatRank0MetaPropertySet(self%massDMOTargetID,basic%mass())
          end if
       else
          ! Negative mass growth - assume all progenitors lose mass at proportionally equal rates.
          ! Compute the total mass in progenitors.
          massTotalProgenitor=basicParent%mass()-massUnresolved
          ! Compute the time available for accretion.
          timeDelta=basicParent%time()-basic%time()
          ! Compute mass growth rate.
          if (timeDelta > 0.0d0) &
               & call basic%         accretionRateSet(                     (massUnresolved/timeDelta)*(basic%mass()/massTotalProgenitor))
          call        basic%floatRank0MetaPropertySet(self%massDMOTargetID,basic%mass()+massUnresolved*basic%mass()/massTotalProgenitor )
       end if
    end if
    return

  contains
    
    double precision function nodeMassUnresolved(node)
      !!{
      Return the unresolved mass for {\normalfont \ttfamily node}.
      !!}
      use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
      implicit none
      type (treeNode          ), intent(inout), pointer :: node
      type (treeNode          )               , pointer :: nodeChild
      class(nodeComponentBasic)               , pointer :: basicChild, basic

      ! Get the basic component.
      basic => node%basic()
      ! Initialize the unresolved mass to the mass of the current node's basic component.
      nodeMassUnresolved=basic%mass()
      ! Remove the mass of all child nodes.
      nodeChild => node%firstChild
      do while (associated(nodeChild))
         basicChild         =>  nodeChild %basic             ()
         nodeMassUnresolved =  +           nodeMassUnresolved   &
              &                -basicChild%mass              ()
         nodeChild          =>  nodeChild %sibling
      end do
      return
    end function nodeMassUnresolved

  end subroutine dmoInterpolateNodeInitialize

  subroutine dmoInterpolateDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorDMOInterpolate), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    class(nodeComponentBasic        ), pointer       :: basic

    basic => node%basic()
    call basic%massAnalytic()
    return
  end subroutine dmoInterpolateDifferentialEvolutionAnalytics

  subroutine dmoInterpolateDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Evolve ``\gls{dmou}'' mass at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorDMOInterpolate), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: time
    class           (nodeComponentBasic        ), pointer       :: basic            , basicParent
    double precision                                            :: massTarget       , timeTarget , &
         &                                                         massRateAccretion
    
    basic             => node %basic        ()
    massRateAccretion =  basic%accretionRate()
    if (massRateAccretion == 0.0d0) return
    basicParent => node       %parent%basic                    (                    )
    massTarget  =  basic             %floatRank0MetaPropertyGet(self%massDMOTargetID)
    timeTarget  =  basicParent       %time                     (                    )
    ! The mass is assumed to grow linearly with time.
    call basic%massSet(massTarget+massRateAccretion*(time-timeTarget))
    return
  end subroutine dmoInterpolateDifferentialEvolutionSolveAnalytics

  subroutine dmoInterpolateNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent.
    !!}
    use :: Error             , only : Error_Report
    use :: Galacticus_Nodes  , only : nodeComponentBasic
    use :: ISO_Varying_String, only : var_str           , varying_string, operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    class    (nodeOperatorDMOInterpolate), intent(inout) :: self
    type     (treeNode                  ), intent(inout) :: node
    type     (treeNode                  ), pointer       :: nodeParent
    class    (nodeComponentBasic        ), pointer       :: basicParent, basic
    type     (varying_string            )                :: message
    character(len=12                    )                :: label
    !$GLC attributes unused :: self

    nodeParent  => node      %parent
    basic       => node      %basic ()
    basicParent => nodeParent%basic ()
    ! Ensure the two halos exist at the same time.
    if (basic%time() /= basicParent%time()) then
       message=var_str("node [")//node%index()//"] has not been evolved to its parent ["//nodeParent%index()//"]"//char(10)
       write (label,'(f12.6)') basic%time()
       message=message//"    node is at time: "//label//" Gyr"//char(10)
       write (label,'(f12.6)') basicParent%time()
       message=message//"  parent is at time: "//label//" Gyr"
       call Error_Report(message//{introspection:location})
    end if
    ! Adjust the mass, target, and accretion rate to that of the parent node.
    call basic%massSet                  (                            basicParent%mass                     (                           ))
    call basic%floatRank0MetaPropertySet(self%       massDMOTargetID,basicParent%floatRank0MetaPropertyGet(self%       massDMOTargetID))
    call basic%         accretionRateSet(                            basicParent%accretionRate            (                           ))
    return
  end subroutine dmoInterpolateNodePromote

  subroutine dmoInterpolateNodesMerge(self,node)
    !!{
    Act on a merger between nodes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorDMOInterpolate), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    class(nodeComponentBasic        ), pointer       :: basic

    ! Shut down mass accretion onto the halo now that it is a satellite.
    basic => node%basic()
    call basic%accretionRateSet(0.0d0)
    return
  end subroutine dmoInterpolateNodesMerge
