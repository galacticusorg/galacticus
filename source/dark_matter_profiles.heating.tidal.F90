!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% A dark matter halo profile heating class which accounts for heating from tidal shocking.

  !# <darkMatterProfileHeating name="darkMatterProfileHeatingTidal">
  !#  <description>A dark matter profile heating model which accounts for heating due to tidal shocking.</description>
  !# </darkMatterProfileHeating>

  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingTidal
     !% A dark matter profile heating class which accounts for heating due to tidal shocking.
     private
   contains
     procedure :: specificEnergy                 => tidalSpecificEnergy
     procedure :: specificEnergyGradient         => tidalSpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero => tidalSpecificEnergyIsEverywhereZero
  end type darkMatterProfileHeatingTidal

  interface darkMatterProfileHeatingTidal
     !% Constructors for the {\normalfont \ttfamily tidal} dark matter profile heating class.
     module procedure tidalConstructorParameters
  end interface darkMatterProfileHeatingTidal

contains

  function tidalConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily tidal} dark matter profile heating scales class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(darkMatterProfileHeatingTidal), target        :: self
    type(inputParameters              ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=darkMatterProfileHeatingTidal()
    return
  end function tidalConstructorParameters

  double precision function tidalSpecificEnergy(self,node,darkMatterProfileDMO_,radius)
    !% Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    implicit none
    class           (darkMatterProfileHeatingTidal), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    double precision                               , intent(in   ) :: radius
    class           (nodeComponentSatellite       ), pointer       :: satellite
    double precision                                               :: specificEnergyOverRadiusSquared
    !$GLC attributes unused :: self, darkMatterProfileDMO_

    satellite                       =>      node     %satellite             ()
    specificEnergyOverRadiusSquared =  max(                                     &
         &                                 +0.0d0                             , &
         &                                 +satellite%tidalHeatingNormalized()  &
         &                                )
    tidalSpecificEnergy             =      +specificEnergyOverRadiusSquared     &
         &                                 *radius                      **2
    return
  end function tidalSpecificEnergy

  double precision function tidalSpecificEnergyGradient(self,node,darkMatterProfileDMO_,radius)
    !% Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    implicit none
    class           (darkMatterProfileHeatingTidal), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    double precision                               , intent(in   ) :: radius
    class           (nodeComponentSatellite       ), pointer       :: satellite
    double precision                                               :: specificEnergyOverRadiusSquared
    !$GLC attributes unused :: self, darkMatterProfileDMO_

    satellite                       =>      node     %satellite             ()
    specificEnergyOverRadiusSquared =  max(                                     &
         &                                 +0.0d0                             , &
         &                                 +satellite%tidalHeatingNormalized()  &
         &                                )
    tidalSpecificEnergyGradient     =      +2.0d0                               &
         &                                 *specificEnergyOverRadiusSquared     &
         &                                 *radius
    return
  end function tidalSpecificEnergyGradient

  logical function tidalSpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !% Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    implicit none
    class(darkMatterProfileHeatingTidal), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    class(nodeComponentSatellite       ), pointer       :: satellite
    !$GLC attributes unused :: self, darkMatterProfileDMO_

    satellite                           => node     %satellite             ()
    tidalSpecificEnergyIsEverywhereZero =  satellite%tidalHeatingNormalized() <= 0.0d0
    return
  end function tidalSpecificEnergyIsEverywhereZero
