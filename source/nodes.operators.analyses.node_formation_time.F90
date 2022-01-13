!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Implements a node operator class that computes the formation time for each node.
  !!}

  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryClass

  !![
  <nodeOperator name="nodeOperatorNodeFormationTime">
   <description>A node operator class that computes the formation time for each node.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorNodeFormationTime
     !!{
     A node operator class that computes the formation time for each node.
     !!}
     private
     class           (darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
     integer                                                            :: nodeFormationTimeID
     logical                                                            :: assumeMonotonicGrowth
     double precision                                                   :: fractionMassFormation
   contains
     final     ::                       nodeFormationTimeDestructor
     procedure :: nodeTreeInitialize => nodeFormationTimeNodeTreeInitialize
     procedure :: nodePromote        => nodeFormationTimeNodePromote
  end type nodeOperatorNodeFormationTime
  
  interface nodeOperatorNodeFormationTime
     !!{
     Constructors for the {\normalfont \ttfamily nodeFormationTime} node operator class.
     !!}
     module procedure nodeFormationTimeConstructorParameters
     module procedure nodeFormationTimeConstructorInternal
  end interface nodeOperatorNodeFormationTime
  
contains

  function nodeFormationTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily nodeFormationTime} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorNodeFormationTime          )                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (darkMatterHaloMassAccretionHistoryClass), pointer       :: darkMatterHaloMassAccretionHistory_
    logical                                                                  :: assumeMonotonicGrowth
    double precision                                                         :: fractionMassFormation

    !![
    <inputParameter>
      <name>fractionMassFormation</name>
      <defaultValue>0.5d0</defaultValue>
      <description>The mass fraction in the main branch progenitor used to define the formation time of each halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>assumeMonotonicGrowth</name>
      <defaultValue>.false.</defaultValue>
      <description>If true assume that halo mass growth is monotonic along each branch when computing node formation times.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"/>
    !!]
    self=nodeOperatorNodeFormationTime(fractionMassFormation,assumeMonotonicGrowth,darkMatterHaloMassAccretionHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    !!]
    return
  end function nodeFormationTimeConstructorParameters

  function nodeFormationTimeConstructorInternal(fractionMassFormation,assumeMonotonicGrowth,darkMatterHaloMassAccretionHistory_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily nodeFormationTime} node operator class.
    !!}
    implicit none
    type            (nodeOperatorNodeFormationTime          )                        :: self
    class           (darkMatterHaloMassAccretionHistoryClass), intent(in   ), target :: darkMatterHaloMassAccretionHistory_
    logical                                                  , intent(in   )         :: assumeMonotonicGrowth
    double precision                                         , intent(in   )         :: fractionMassFormation
    !![
    <constructorAssign variables="fractionMassFormation, assumeMonotonicGrowth, *darkMatterHaloMassAccretionHistory_"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="nodeFormationTime" id="self%nodeFormationTimeID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function nodeFormationTimeConstructorInternal

  subroutine nodeFormationTimeDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily nodeFormationTime} node operator class.
    !!}
    implicit none
    type(nodeOperatorNodeFormationTime), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    !!]
    return
  end subroutine nodeFormationTimeDestructor

  subroutine nodeFormationTimeNodeTreeInitialize(self,node)
    !!{
    Initialize node formation times.
    !!}
    use :: Dark_Matter_Halo_Formation_Times, only : Dark_Matter_Halo_Formation_Time
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    use :: Merger_Tree_Walkers             , only : mergerTreeWalkerAllNodes
    implicit none
    class           (nodeOperatorNodeFormationTime), intent(inout), target  :: self
    type            (treeNode                     ), intent(inout), target  :: node
    type            (treeNode                     )               , pointer :: nodeFormation, nodeWork, &
         &                                                                     nodeParent
    class           (nodeComponentBasic           )               , pointer :: basic
    type            (mergerTreeWalkerAllNodes     )                         :: treeWalker
    double precision                                                        :: timeFormation

    
    treeWalker=mergerTreeWalkerAllNodes(node%hostTree,spanForest=.false.)
    do while (treeWalker%next(nodeWork))
       if (.not.self%assumeMonotonicGrowth) then
          basic         => nodeWork%basic()
          timeFormation =  Dark_Matter_Halo_Formation_Time(                                                                              &
               &                                           node                               =nodeWork                                , &
               &                                           formationMassFraction              =self%fractionMassFormation              , &
               &                                           darkMatterHaloMassAccretionHistory_=self%darkMatterHaloMassAccretionHistory_  &
               &                                          )
          call basic%metaPropertySet(self%nodeFormationTimeID,timeFormation)
       else if (.not.associated(nodeWork%firstChild)) then
          nodeParent    => nodeWork
          nodeFormation => nodeWork
          do while (associated(nodeParent))
             basic         => nodeParent%basic()
             timeFormation =  Dark_Matter_Halo_Formation_Time(                                                                              &
                  &                                           node                               =nodeParent                              , &
                  &                                           nodeFormation                      =nodeFormation                           , &
                  &                                           formationMassFraction              =self%fractionMassFormation              , &
                  &                                           darkMatterHaloMassAccretionHistory_=self%darkMatterHaloMassAccretionHistory_  &
                  &                                          )
             call basic%metaPropertySet(self%nodeFormationTimeID,timeFormation)
             if (nodeParent%isPrimaryProgenitor()) then
                nodeParent => nodeParent%parent
             else
                nodeParent => null()
             end if
          end do
       end if
    end do
    return
  end subroutine nodeFormationTimeNodeTreeInitialize
 
  subroutine nodeFormationTimeNodePromote(self,node)
    !!{
    Promote node major merger times.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorNodeFormationTime), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(nodeComponentBasic           ), pointer       :: basic, basicParent
    
    basic       => node       %basic()
    basicParent => node%parent%basic()
    call basic%metaPropertySet(self%nodeFormationTimeID,basicParent%metaPropertyGet(self%nodeFormationTimeID))
    return
  end subroutine nodeFormationTimeNodePromote
