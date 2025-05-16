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
An implementation of the hot halo outflow reincorporation class in which reincorporation occurs on a multiple of the halo
dynamical timescale.
!!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <hotHaloOutflowReincorporation name="hotHaloOutflowReincorporationHaloDynamicalTime">
   <description>An implementation of the hot halo outflow reincorporation class in which reincorporation occurs on a multiple of the halo dynamical timescale.</description>
  </hotHaloOutflowReincorporation>
  !!]
  type, extends(hotHaloOutflowReincorporationClass) :: hotHaloOutflowReincorporationHaloDynamicalTime
     !!{
     An implementation of the hot halo outflow reincorporation class in which reincorporation occurs on a multiple of the halo dynamical timescale.
     !!}
     private
     double precision                                    :: multiplier
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::         haloDynamicalTimeDestructor
     procedure :: rate => haloDynamicalTimeRate
  end type hotHaloOutflowReincorporationHaloDynamicalTime

  interface hotHaloOutflowReincorporationHaloDynamicalTime
     !!{
     Constructors for the \refClass{hotHaloOutflowReincorporationHaloDynamicalTime} hot halo outflow reincorporation class.
     !!}
     module procedure haloDynamicalTimeConstructorParameters
     module procedure haloDynamicalTimeConstructorInternal
  end interface hotHaloOutflowReincorporationHaloDynamicalTime

contains

  function haloDynamicalTimeConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily haloDynamicalTime} hot halo outflow reincorporation class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hotHaloOutflowReincorporationHaloDynamicalTime)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                      ), pointer       :: darkMatterHaloScale_
    double precision                                                                :: multiplier

    !![
    <inputParameter>
      <name>multiplier</name>
      <defaultValue>5.0d0</defaultValue>
      <description>Specifies the rate at which reheated mass is returned to the hot phase in units of the inverse halo dynamical timed.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=hotHaloOutflowReincorporationHaloDynamicalTime(multiplier,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function haloDynamicalTimeConstructorParameters

  function haloDynamicalTimeConstructorInternal(multiplier,darkMatterHaloScale_) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily haloDynamicalTime} hot halo outflow reincorporation class.
    !!}
    implicit none
    type            (hotHaloOutflowReincorporationHaloDynamicalTime)                        :: self
    double precision                                                , intent(in   )         :: multiplier
    class           (darkMatterHaloScaleClass                      ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="multiplier, *darkMatterHaloScale_"/>
    !!]

    return
  end function haloDynamicalTimeConstructorInternal

  subroutine haloDynamicalTimeDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloOutflowReincorporationHaloDynamicalTime} hot halo outflow reincorporation class.
    !!}
    implicit none
    type(hotHaloOutflowReincorporationHaloDynamicalTime), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine haloDynamicalTimeDestructor

  double precision function haloDynamicalTimeRate(self,node)
    !!{
    Return the rate of mass reincorporation for outflowed gas in the hot halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class(hotHaloOutflowReincorporationHaloDynamicalTime), intent(inout) :: self
    type (treeNode                                      ), intent(inout) :: node
    class(nodeComponentHotHalo                          ), pointer       :: hotHalo

    hotHalo               =>  node   %hotHalo                                (    )
    haloDynamicalTimeRate =  +hotHalo%outflowedMass                          (    )  &
         &                   *self   %multiplier                                     &
         &                   /self   %darkMatterHaloScale_%timescaleDynamical(node)
    return
  end function haloDynamicalTimeRate
