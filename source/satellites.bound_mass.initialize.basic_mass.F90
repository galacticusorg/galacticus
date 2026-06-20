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
  Implementation of a satellite bound mass initializor class that sets the initial bound mass to the basic node mass.
  !!}

  !![
  <satelliteMassBoundInitializor name="satelliteMassBoundInitializorBasicMass" docformat="rst">
   <description>
   A satellite bound mass initializor class that sets the initial bound mass of the satellite halo to the basic node mass, i.e.\ the total mass of the halo at the time it becomes a satellite.
   </description>
  </satelliteMassBoundInitializor>
  !!]
  type, extends(satelliteMassBoundInitializorClass) :: satelliteMassBoundInitializorBasicMass
     !!{RST
     Implementation of a satellite bound mass initializor class that sets the initial bound mass to the basic node mass.
     !!}
     private
   contains
     procedure :: massBound => basicMassMassBound
  end type satelliteMassBoundInitializorBasicMass

  interface satelliteMassBoundInitializorBasicMass
     !!{RST
     Constructors for the :galacticus-class:`satelliteMassBoundInitializorBasicMass` satellite bound mass initializor class.
     !!}
     module procedure basicMassConstructorParameters
  end interface satelliteMassBoundInitializorBasicMass

contains

  function basicMassConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteMassBoundInitializorBasicMass` satellite bound mass initializor class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(satelliteMassBoundInitializorBasicMass)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self=satelliteMassBoundInitializorBasicMass()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function basicMassConstructorParameters

  double precision function basicMassMassBound(self,node) result(massBound)
    !!{RST
    Returns the initial bound mass of a satellite halo set equal to the node basic mass.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(satelliteMassBoundInitializorBasicMass), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node
    class(nodeComponentBasic                    ), pointer       :: basic
    !$GLC attributes unused :: self

    basic     => node %basic()
    massBound =  basic%mass ()
    return
  end function basicMassMassBound

