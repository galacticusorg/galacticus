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
  Implementation of a stellar feedback model which scales with formation halo velocity.
  !!}

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsHaloScalingFormationHalo">
   <description>An stellar feedback model which scales with halo velocity.</description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsHaloScaling) :: stellarFeedbackOutflowsHaloScalingFormationHalo
     !!{
     Implementation of a stellar feedback model which scales with formation halo velocity.
     !!}
     private
   contains
     procedure :: node => haloScalingFormationHaloNode
  end type stellarFeedbackOutflowsHaloScalingFormationHalo

  interface stellarFeedbackOutflowsHaloScalingFormationHalo
     !!{
     Constructors for the halo scaling fraction stellar feedback class.
     !!}
     module procedure haloScalingFormationHaloConstructorParameters
     module procedure haloScalingFormationHaloConstructorInternal
  end interface stellarFeedbackOutflowsHaloScalingFormationHalo

contains

  function haloScalingFormationHaloConstructorParameters(parameters) result(self)
    !!{
    Constructor for the halo scaling (formation halo) fraction stellar feedback class which takes a parameter set as input.
    !!}
    implicit none
    type            (stellarFeedbackOutflowsHaloScalingFormationHalo)                :: self
    type            (inputParameters                                ), intent(inout) :: parameters
    double precision                                                                 :: fraction            , exponentRedshift, &
         &                                                                              exponentVelocity
    class           (cosmologyFunctionsClass                        ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                       ), pointer       :: darkMatterHaloScale_

    !![
    <inputParameter>
      <name>fraction</name>
      <source>parameters</source>
      <defaultValue>0.01d0</defaultValue>
      <description>The ratio of outflow rate to star formation rate in disks.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentVelocity</name>
      <source>parameters</source>
      <defaultValue>-2.0d0</defaultValue>
      <description>The exponent of virial velocity in the outflow rate in disks.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentRedshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent of redshift in the outflow rate in disks.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=stellarFeedbackOutflowsHaloScalingFormationHalo(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function haloScalingFormationHaloConstructorParameters

  function haloScalingFormationHaloConstructorInternal(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the halo scaling (formation halo) stellar feedback class.
    !!}
    implicit none
    type            (stellarFeedbackOutflowsHaloScalingFormationHalo)                        :: self
    double precision                                                 , intent(in   )         :: fraction            , exponentRedshift, &
         &                                                                                      exponentVelocity
    class           (cosmologyFunctionsClass                        ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                       ), intent(in   ), target :: darkMatterHaloScale_

    self%stellarFeedbackOutflowsHaloScaling=stellarFeedbackOutflowsHaloScaling(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_,darkMatterHaloScale_)
    return
  end function haloScalingFormationHaloConstructorInternal

  function haloScalingFormationHaloNode(self,component) result(node)
    !!{
    Returns a pointer to the node from which to extract halo properties.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type (treeNode                                       ), pointer       :: node
    class(stellarFeedbackOutflowsHaloScalingFormationHalo), intent(inout) :: self
    class(nodeComponent                                  ), intent(in   ) :: component
    !$GLC attributes unused :: self

    if (associated(component%hostNode%formationNode)) then
       node => component%hostNode%formationNode
    else
       node => null()
       call Error_Report('no formation node exists'//{introspection:location})
    end if
    return
  end function haloScalingFormationHaloNode
