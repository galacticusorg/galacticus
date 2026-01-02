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
  Implementation of a simple infall torque calculation in which a fixed fraction of angular momentum is lost during infall.
  !!}

  !![
  <coolingInfallTorque name="coolingInfallTorqueFixed">
   <description>
    A cooling torque class in which a fixed fraction of angular momentum is lost during infall.
   </description>
  </coolingInfallTorque>
  !!]
  type, extends(coolingInfallTorqueClass) :: coolingInfallTorqueFixed
     !!{
     Implementation of a simple infall torque calculation in which a fixed fraction of angular momentum is lost during infall.
     !!}
     private
     double precision :: fractionLossAngularMomentum
   contains
     procedure :: fractionAngularMomentumLoss => fixedFractionAngularMomentumLoss
  end type coolingInfallTorqueFixed

  interface coolingInfallTorqueFixed
     !!{
     Constructors for the cooling radius infall radii class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface coolingInfallTorqueFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{coolingInfallTorqueFixed} class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (coolingInfallTorqueFixed)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    double precision                                          :: fractionLossAngularMomentum
    
    !![
    <inputParameter>
      <name>fractionLossAngularMomentum</name>
      <defaultValue>0.3d0</defaultValue>
      <description>Specifies the fraction of angular momentum that is lost from cooling/infalling gas.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=coolingInfallTorqueFixed(fractionLossAngularMomentum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(fractionLossAngularMomentum) result(self)
    !!{
    Internal constructor for the \refClass{coolingInfallTorqueFixed} class.
    !!}
    implicit none
    type            (coolingInfallTorqueFixed)                :: self
    double precision                          , intent(in   ) :: fractionLossAngularMomentum
    !![
    <constructorAssign variables="fractionLossAngularMomentum"/>
    !!]

    return
  end function fixedConstructorInternal

  double precision function fixedFractionAngularMomentumLoss(self,node) result(fractionAngularMomentumLoss)
    !!{
    Return the fraction of angular momentum lost during infall.
    !!}
    implicit none
    class(coolingInfallTorqueFixed), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node

    fractionAngularMomentumLoss=self%fractionLossAngularMomentum
    return
  end function fixedFractionAngularMomentumLoss
