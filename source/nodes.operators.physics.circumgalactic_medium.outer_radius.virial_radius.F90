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
  Implements a node operator class that fixes the :term:`CGM` outer radius to the virial radius.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <nodeOperator name="nodeOperatorCGMOuterRadiusVirialRadius" docformat="rst">
   <description>
   Sets the outer boundary of the circumgalactic medium (:term:`CGM`) equal to the halo virial radius as computed by a :galacticus-class:`darkMatterHaloScaleClass` object, analytically enforcing this constraint during differential evolution.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMOuterRadiusVirialRadius
     !!{RST
     A node operator class that fixes the :term:`CGM` outer radius to the virial radius.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                                        cgmOuterRadiusVirialRadiusDestructor
     procedure :: differentialEvolutionAnalytics      => cgmOuterRadiusVirialRadiusDifferentialEvolutionAnalytics
     procedure :: preDeterminedSolveAnalytics         => cgmOuterRadiusVirialRadiusSolveAnalytics
     procedure :: differentialEvolutionSolveAnalytics => cgmOuterRadiusVirialRadiusSolveAnalytics
     procedure :: nodesMerge                          => cgmOuterRadiusVirialRadiusNodesMerge
  end type nodeOperatorCGMOuterRadiusVirialRadius
  
  interface nodeOperatorCGMOuterRadiusVirialRadius
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorCGMOuterRadiusVirialRadius` node operator class.
     !!}
     module procedure cgmOuterRadiusVirialRadiusConstructorParameters
     module procedure cgmOuterRadiusVirialRadiusConstructorInternal
  end interface nodeOperatorCGMOuterRadiusVirialRadius
  
contains
  
  function cgmOuterRadiusVirialRadiusConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorCGMOuterRadiusVirialRadius` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorCGMOuterRadiusVirialRadius)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodeOperatorCGMOuterRadiusVirialRadius(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function cgmOuterRadiusVirialRadiusConstructorParameters

  function cgmOuterRadiusVirialRadiusConstructorInternal(darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorCGMOuterRadiusVirialRadius` node operator class.
    !!}
    implicit none
    type (nodeOperatorCGMOuterRadiusVirialRadius)                        :: self
    class(darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function cgmOuterRadiusVirialRadiusConstructorInternal

  subroutine cgmOuterRadiusVirialRadiusDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorCGMOuterRadiusVirialRadius` node operator class.
    !!}
    implicit none
    type(nodeOperatorCGMOuterRadiusVirialRadius), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine cgmOuterRadiusVirialRadiusDestructor
  
  subroutine cgmOuterRadiusVirialRadiusDifferentialEvolutionAnalytics(self,node)
    !!{RST
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class(nodeOperatorCGMOuterRadiusVirialRadius), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node
    class(nodeComponentHotHalo                  ), pointer       :: hotHalo

    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo exists - nothing to do.
    class default
       call hotHalo%outerRadiusAnalytic()
    end select
    return
  end subroutine cgmOuterRadiusVirialRadiusDifferentialEvolutionAnalytics

  subroutine cgmOuterRadiusVirialRadiusSolveAnalytics(self,node,time)
    !!{RST
    Set values of analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class           (nodeOperatorCGMOuterRadiusVirialRadius), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    double precision                                        , intent(in   ) :: time
    class           (nodeComponentHotHalo                  ), pointer       :: hotHalo
    !$GLC attributes unused :: time
    
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo exists - nothing to do.
    class default
       call hotHalo%outerRadiusSet(self%darkMatterHaloScale_%radiusVirial(node))
    end select
    return
  end subroutine cgmOuterRadiusVirialRadiusSolveAnalytics

  subroutine cgmOuterRadiusVirialRadiusNodesMerge(self,node)
    !!{RST
    Act on a merger between nodes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class(nodeOperatorCGMOuterRadiusVirialRadius), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node
    class(nodeComponentHotHalo                  ), pointer       :: hotHalo

    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo exists - nothing to do.
    class default
       call hotHalo%outerRadiusSet(self%darkMatterHaloScale_%radiusVirial(node))
    end select
    return
  end subroutine cgmOuterRadiusVirialRadiusNodesMerge
  
