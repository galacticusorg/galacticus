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
  Implements a node operator class that sets dark matter profile scale radius.
  !!}

  use :: Galacticus_Nodes          , only : treeNodeLinkedList
  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadiusClass
  
  type :: branchStack
     !!{
     Type used to create a stack of branches being processed.
     !!}
     double precision                              :: mass              , radiusScale
     type            (branchStack       ), pointer :: next     => null()
     type            (treeNodeLinkedList), pointer :: node     => null()
     integer         (kind_int8         )          :: indexTip
  end type branchStack
  
  !![
  <nodeOperator name="nodeOperatorDarkMatterProfileScaleSet">
   <description>
    A node operator class that sets dark matter profile scale radius.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfileScaleSet
     !!{
     A node operator class that sets the dark matter profile scale radius in halos.
     !!}
     private
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     type            (branchStack                      ), pointer :: branch                        => null()
     double precision                                             :: factorReset
     logical                                                      :: forward
   contains
     final     ::                       darkMatterProfileScaleSetConstructorDestructor
     procedure :: nodeTreeInitialize => darkMatterProfileScaleSetNodeTreeInitialize
     procedure :: nodePromote        => darkMatterProfileScaleSetNodePromote
  end type nodeOperatorDarkMatterProfileScaleSet
  
  interface nodeOperatorDarkMatterProfileScaleSet
     !!{
     Constructors for the \refClass{nodeOperatorDarkMatterProfileScaleSet} node operator class.
     !!}
     module procedure darkMatterProfileScaleSetConstructorParameters
     module procedure darkMatterProfileScaleSetConstructorInternal
  end interface nodeOperatorDarkMatterProfileScaleSet

contains
  
  function darkMatterProfileScaleSetConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileScaleSet} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkMatterProfileScaleSet)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (darkMatterProfileScaleRadiusClass    ), pointer       :: darkMatterProfileScaleRadius_
    double precision                                                       :: factorReset
    logical                                                                :: forward

    !![
    <inputParameter>
      <name>factorReset</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The factor by which a node must increase in mass before its scale radius is reset.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>forward</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, updates to the scale radius are determined by walking forward along the branch until the mass exceeds by {\normalfont \ttfamily [factor]} that for which the scale radius was last computed. If false, updates are computed by walking backward along the branch until the mass is less than $1/${\normalfont \ttfamily [factor]} of that for which the scale radius was last computed.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    self=nodeOperatorDarkMatterProfileScaleSet(factorReset,forward,darkMatterProfileScaleRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function darkMatterProfileScaleSetConstructorParameters

  function darkMatterProfileScaleSetConstructorInternal(factorReset,forward,darkMatterProfileScaleRadius_) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileScaleSet} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkMatterProfileScaleSet)                        :: self
    class           (darkMatterProfileScaleRadiusClass    ), intent(in   ), target :: darkMatterProfileScaleRadius_
    double precision                                       , intent(in   )         :: factorReset
    logical                                                , intent(in   )         :: forward
    !![
    <constructorAssign variables="factorReset, forward, *darkMatterProfileScaleRadius_"/>
    !!]

    return
  end function darkMatterProfileScaleSetConstructorInternal

  subroutine darkMatterProfileScaleSetConstructorDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorDarkMatterProfileScaleSet} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfileScaleSet), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    !!]
    return
  end subroutine darkMatterProfileScaleSetConstructorDestructor

  subroutine darkMatterProfileScaleSetNodeTreeInitialize(self,node)
    !!{
    Initialize dark matter profile scale radii.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileScaleSet), intent(inout), target  :: self
    type            (treeNode                             ), intent(inout), target  :: node
    class           (nodeComponentBasic                   )               , pointer :: basic              , basicDescendant
    class           (nodeComponentDarkMatterProfile       )               , pointer :: darkMatterProfile
    type            (branchStack                          )               , pointer :: branchPrevious
    type            (treeNodeLinkedList                   )               , pointer :: nodeNext
    type            (treeNode                             )               , pointer :: nodeDescendant
    double precision                                                                :: radiusScalePrevious, radiusScaleUpdated
    logical                                                                         :: radiusScaleNew

    if (self%factorReset > 0.0d0) then
       basic => node%basic()
       if (.not.associated(node%firstChild)) then
          ! This is the tip of a branch. Push the branch onto the stack, and build a stack of nodes at which new scale radii must
          ! be computed.
          branchPrevious => self%branch
          nullify (self%branch)
          allocate(self%branch)
          self%branch%next     => branchPrevious
          self%branch%indexTip =  node          %uniqueID()
          if (self%forward) then
             self%branch%mass=basic%mass()
          else
             nodeDescendant => node
             do while (nodeDescendant%isPrimaryProgenitor())
                nodeDescendant => nodeDescendant%parent
             end do
             allocate(self%branch%node)
             basicDescendant                  => nodeDescendant %basic()
             self           %branch%node%node => nodeDescendant
             self           %branch%node%next => null()
             self           %branch%mass      =  basicDescendant%mass ()
             do while (associated(nodeDescendant))
                basicDescendant => nodeDescendant%basic()
                if (basicDescendant%mass() < self%branch%mass/self%factorReset) then
                   nodeNext => self%branch%node
                   nullify (self%branch%node)
                   allocate(self%branch%node)
                   self%branch%mass      =  basicDescendant%mass()
                   self%branch%node%node => nodeDescendant
                   self%branch%node%next => nodeNext
                end if
                nodeDescendant => nodeDescendant%firstChild
             end do
          end if
          ! In this case we must always compute a new scale radius.
          radiusScaleNew=.true.
       else
          ! This is not a branch tip.
          if (self%forward) then
             ! A new scale radius is computed only if the prior mass is exceeded by the required factor.
             radiusScaleNew=basic%mass() > self%factorReset*self%branch%mass
             ! If a new scale radius is to be computed, record the mass at which this occurred.
             if (radiusScaleNew) self%branch%mass=basic%mass()
          else
             ! A new scale radius is computed only if we have reached the prior node for which it was computed.
             radiusScaleNew=associated(self%branch%node%node,node)
             if (radiusScaleNew) then
                nodeNext => self%branch%node%next
                deallocate(self%branch%node)
                self%branch%node => nodeNext
             end if
             ! If there are no more nodes along this branch at which to compute the scale radius, mark that scale radius is not to
             ! be computed.
             if (.not.associated(self%branch%node)) radiusScaleNew=.false.
          end if
       end if
    else
       ! The reset factor is zero - we assign a new scale radius for every node.
       radiusScaleNew=.true.
    end if
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    if (radiusScaleNew) then
       ! Compute a new scale radius, and save it on our current branch.
       if (self%forward .or. self%factorReset <= 0.0d0) then
          nodeDescendant =>             node
       else
          nodeDescendant => self%branch%node%node
       end if
       radiusScaleUpdated=self%darkMatterProfileScaleRadius_%radius(nodeDescendant)
       ! If the current node was one for which we needed to compute a new scale radius, pop the node stack here. This can happen
       ! if the branch tip node is a node for which we must compute the scale radius.
       if (self%factorReset > 0.0d0 .and. .not.self%forward  .and. associated(self%branch%node%node,node)) then
          nodeNext => self%branch%node%next
          deallocate(self%branch%node)
          self%branch%node => nodeNext
       end if
       ! If using the "forward" approach, or this is the branch tip, we update the actual scale radius now, as we want to apply it
       ! to this node. Otherwise, use the previous value.
       if (self%factorReset <= 0.0d0 .or. self%forward .or. .not.associated(node%firstChild)) then
          radiusScalePrevious=            radiusScaleUpdated
       else
          radiusScalePrevious=self%branch%radiusScale
       end if
       if (associated(self%branch)) self%branch%radiusScale=radiusScaleUpdated
    else
       ! Use the previously computed scale radius for this branch.
       radiusScalePrevious=self%branch%radiusScale
    end if
    call darkMatterProfile%scaleSet(radiusScalePrevious)
    ! If the node is not the primary progenitor, then we have reached the end of this branch.
    if (self%factorReset > 0.0d0 .and. .not.node%isPrimaryProgenitor()) then
       ! Remove the branch from the stack and restore the previous branch.
       branchPrevious        => self          %branch%next
       deallocate(self%branch)
       self          %branch => branchPrevious
   end if
    return
  end subroutine darkMatterProfileScaleSetNodeTreeInitialize
    
  subroutine darkMatterProfileScaleSetNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the scale radius
    of {\normalfont \ttfamily node} to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileScaleSet), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile       ), pointer       :: darkMatterProfile, darkMatterProfileParent
    
    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    call darkMatterProfile%scaleSet(darkMatterProfileParent%scale())
    return
  end subroutine darkMatterProfileScaleSetNodePromote
