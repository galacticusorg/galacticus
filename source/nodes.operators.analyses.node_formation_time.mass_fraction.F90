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

  !!{
  Implements a node operator class that computes the formation time for each node based on a mass fraction definition.
  !!}

  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryClass

  !![
  <nodeOperator name="nodeOperatorNodeFormationTimeMassFraction">
   <description>A node operator class that computes the formation time for each node based on a mass fraction definition.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorNodeFormationTimeMassFraction
     !!{
     A node operator class that computes the formation time for each node based on a mass fraction definition.
     !!}
     private
     class           (darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
     integer                                                            :: nodeFormationTimeID
     logical                                                            :: assumeMonotonicGrowth
     double precision                                                   :: fractionMassFormation
   contains
     final     ::                       nodeFormationTimeMassFractionDestructor
     procedure :: nodeTreeInitialize => nodeFormationTimeMassFractionNodeTreeInitialize
     procedure :: nodePromote        => nodeFormationTimeMassFractionNodePromote
  end type nodeOperatorNodeFormationTimeMassFraction
  
  interface nodeOperatorNodeFormationTimeMassFraction
     !!{
     Constructors for the \refClass{nodeOperatorNodeFormationTimeMassFraction} node operator class.
     !!}
     module procedure nodeFormationTimeMassFractionConstructorParameters
     module procedure nodeFormationTimeMassFractionConstructorInternal
  end interface nodeOperatorNodeFormationTimeMassFraction
  
contains

  function nodeFormationTimeMassFractionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorNodeFormationTimeMassFraction} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorNodeFormationTimeMassFraction)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (darkMatterHaloMassAccretionHistoryClass  ), pointer       :: darkMatterHaloMassAccretionHistory_
    logical                                                                    :: assumeMonotonicGrowth
    double precision                                                           :: fractionMassFormation

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
    self=nodeOperatorNodeFormationTimeMassFraction(fractionMassFormation,assumeMonotonicGrowth,darkMatterHaloMassAccretionHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    !!]
    return
  end function nodeFormationTimeMassFractionConstructorParameters

  function nodeFormationTimeMassFractionConstructorInternal(fractionMassFormation,assumeMonotonicGrowth,darkMatterHaloMassAccretionHistory_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorNodeFormationTimeMassFraction} node operator class.
    !!}
    implicit none
    type            (nodeOperatorNodeFormationTimeMassFraction)                        :: self
    class           (darkMatterHaloMassAccretionHistoryClass  ), intent(in   ), target :: darkMatterHaloMassAccretionHistory_
    logical                                                    , intent(in   )         :: assumeMonotonicGrowth
    double precision                                           , intent(in   )         :: fractionMassFormation
    !![
    <constructorAssign variables="fractionMassFormation, assumeMonotonicGrowth, *darkMatterHaloMassAccretionHistory_"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="nodeFormationTime" id="self%nodeFormationTimeID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function nodeFormationTimeMassFractionConstructorInternal

  subroutine nodeFormationTimeMassFractionDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorNodeFormationTimeMassFraction} node operator class.
    !!}
    implicit none
    type(nodeOperatorNodeFormationTimeMassFraction), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    !!]
    return
  end subroutine nodeFormationTimeMassFractionDestructor

  subroutine nodeFormationTimeMassFractionNodeTreeInitialize(self,node)
    !!{
    Initialize node formation times.
    !!}
    use :: Dark_Matter_Halo_Formation_Times, only : Dark_Matter_Halo_Formation_Time
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    implicit none
    class           (nodeOperatorNodeFormationTimeMassFraction), intent(inout), target  :: self
    type            (treeNode                                 ), intent(inout), target  :: node
    type            (treeNode                                 )               , pointer :: nodeFormation, nodeParent
    class           (nodeComponentBasic                       )               , pointer :: basic
    double precision                                                                    :: timeFormation
    
    if (.not.self%assumeMonotonicGrowth) then
       basic         => node%basic()
       timeFormation =  Dark_Matter_Halo_Formation_Time(                                                                              &
            &                                           node                               =node                                    , &
            &                                           formationMassFraction              =self%fractionMassFormation              , &
            &                                           darkMatterHaloMassAccretionHistory_=self%darkMatterHaloMassAccretionHistory_  &
            &                                          )
       call basic%floatRank0MetaPropertySet(self%nodeFormationTimeID,timeFormation)
    else if (.not.associated(node%firstChild)) then
       nodeParent    => node
       nodeFormation => node
       do while (associated(nodeParent))
          basic         => nodeParent%basic()
          timeFormation =  Dark_Matter_Halo_Formation_Time(                                                                              &
               &                                           node                               =nodeParent                              , &
               &                                           nodeFormation                      =nodeFormation                           , &
               &                                           formationMassFraction              =self%fractionMassFormation              , &
               &                                           darkMatterHaloMassAccretionHistory_=self%darkMatterHaloMassAccretionHistory_  &
               &                                          )
          call basic%floatRank0MetaPropertySet(self%nodeFormationTimeID,timeFormation)
          if (nodeParent%isPrimaryProgenitor()) then
             nodeParent => nodeParent%parent
          else
             nodeParent => null()
          end if
       end do
    end if
    return
  end subroutine nodeFormationTimeMassFractionNodeTreeInitialize
 
  subroutine nodeFormationTimeMassFractionNodePromote(self,node)
    !!{
    Promote node major merger times.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorNodeFormationTimeMassFraction), intent(inout) :: self
    type (treeNode                                 ), intent(inout) :: node
    class(nodeComponentBasic                       ), pointer       :: basic, basicParent
    
    basic       => node       %basic()
    basicParent => node%parent%basic()
    call basic%floatRank0MetaPropertySet(self%nodeFormationTimeID,basicParent%floatRank0MetaPropertyGet(self%nodeFormationTimeID))
    return
  end subroutine nodeFormationTimeMassFractionNodePromote
