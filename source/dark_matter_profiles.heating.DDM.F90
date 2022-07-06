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
  A dark matter halo profile heating class which accounts for heating from decays.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingDDM">
   <description>
    $E_{\mathrm{specific}} = 0.5 (1 - e^{-\Delta t/ \tau}) \epsilon^2 c^2$
   </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingDDM
     !!{
     A dark matter profile heating class which accounts for heating due to decays.
     !!}
     private
     double precision :: lifetime, massSplitting
  contains
     procedure :: specificEnergy                  => DDMSpecificEnergy
     procedure :: specificEnergyGradient          => DDMSpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero  => DDMSpecificEnergyIsEverywhereZero
  end type darkMatterProfileHeatingDDM

  interface darkMatterProfileHeatingDDM
     !!{
     Constructors for the {\normalfont \ttfamily DDM} dark matter profile heating class.
     !!}
     module procedure DDMConstructorParameters
     module procedure DDMConstructorInternal
  end interface darkMatterProfileHeatingDDM

contains

  function DDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily DDM} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileHeatingDDM), target        :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: lifetime  , massSplitting
         
    !![
    <inputParameter>
      <name>lifetime</name>
      <defaultValue>-1.0d0</defaultValue> !! remove default values
      <source>parameters</source>
      <description>The lifetime of the dark matter particles.</description>
    </inputParameter>
    <inputParameter>
      <name>massSplitting</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>It is equal to $(M - m) / M$ where $M$ is the mass of the original particle and $m$ is the mass of the daughter particle.</description>
    </inputParameter>
    !!]
    self=darkMatterProfileHeatingDDM(lifetime, massSplitting)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function DDMConstructorParameters

  function DDMConstructorInternal(lifetime, massSplitting) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily DDM} dark matter profile heating scales class.
    !!}
    implicit none
    type            (darkMatterProfileHeatingDDM)                :: self
    double precision                               , intent(in   ) :: lifetime, massSplitting
    !![
    <constructorAssign variables="lifetime, massSplitting"/>
    !!]
    return
  end function DDMConstructorInternal

  double precision function DDMSpecificEnergy(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class           (darkMatterProfileHeatingDDM), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    double precision                             , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass  ), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic         ), pointer       :: basic

    basic             => node%basic()
    DDMSpecificEnergy = +0.5d0*(                                   &
         &              +1.0d0 - exp(-basic%time() / self%lifetime)&
         &                     )                                   &
         &              *self%massSplitting**2                     &
         &              *(speedLight/kilo)**2
    return
  end function DDMSpecificEnergy

  double precision function DDMSpecificEnergyGradient(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class           (darkMatterProfileHeatingDDM), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    
    DDMSpecificEnergyGradient=+0.0d0
    return
  end function DDMSpecificEnergyGradient

  logical function DDMSpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !!{
    Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileHeatingDDM), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    !$GLC attributes unused :: darkMatterProfileDMO_

    DDMSpecificEnergyIsEverywhereZero=(self%lifetime <= 0.0d0) .or. (self%massSplitting <= 0.0d0)
    return
  end function DDMSpecificEnergyIsEverywhereZero
