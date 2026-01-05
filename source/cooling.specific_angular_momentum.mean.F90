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
  Implementation of a specific angular momentum of cooling gas class in which all gas has the mean specific angular momentum of
  the hot gas halo.
  !!}

  !![
  <coolingSpecificAngularMomentum name="coolingSpecificAngularMomentumMean">
   <description>
    A cooling specific angular momentum class in which the specific angular momentum of cooling gas is given by
  \begin{equation}
     j_\mathrm{cool} = J_\mathrm{hot}/M_\mathrm{hot},
    \end{equation}
    where $J_\mathrm{hot}$ and $M_\mathrm{hot}$ are the total angular momentum and mass of the hot halo respectively.
   </description>
  </coolingSpecificAngularMomentum>
  !!]
  type, extends(coolingSpecificAngularMomentumClass) :: coolingSpecificAngularMomentumMean
     !!{
     Implementation of the specific angular momentum of cooling gas class in which all gas has the mean specific angular momentum of the hot gas halo.
     !!}
     private
   contains
     procedure :: angularMomentumSpecific => meanAngularMomentumSpecific
  end type coolingSpecificAngularMomentumMean

  interface coolingSpecificAngularMomentumMean
     !!{
     Constructors for the mean specific angular momentum of cooling gas class.
     !!}
     module procedure meanConstructorParameters
  end interface coolingSpecificAngularMomentumMean

contains

  function meanConstructorParameters(parameters) result(self)
    !!{
    Constructor for the mean freefall radius class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(coolingSpecificAngularMomentumMean)                :: self
    type(inputParameters                   ), intent(inout) :: parameters

    self=coolingSpecificAngularMomentumMean()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function meanConstructorParameters

  double precision function meanAngularMomentumSpecific(self,node,radius)
    !!{
    Return the specific angular momentum of cooling gas in the mean model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class           (coolingSpecificAngularMomentumMean ), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: radius
    class           (nodeComponentHotHalo               ), pointer       :: hotHalo
    !$GLC attributes unused :: self, radius

    ! Compute mean specific angular momentum from the hot halo component.
    hotHalo                     =>  node   %hotHalo        ()
    meanAngularMomentumSpecific =  +hotHalo%angularMomentum() &
         &                         /hotHalo%mass           ()
    return
  end function meanAngularMomentumSpecific
