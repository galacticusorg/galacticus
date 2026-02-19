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
Contains a module functions supporting the black hole physics in the \refClass{nodeOperatorClass} class.
!!}

module Nodes_Operators_Black_Holes_Utilities
  !!{
  Provides functions supporting the black hole physics in the \refClass{nodeOperatorClass} class.
  !!}
  private
  public :: blackHolesRecoilEscapes

contains

  logical function blackHolesRecoilEscapes(node,radius,radiusEscape,velocityRecoil,ignoreCentralBlackHole) result(escapes)
    !!{
    Return true if the given recoil velocity is sufficient to eject a black hole from the halo.
    !!}
    use :: Galacticus_Nodes          , only : treeNode
    use :: Coordinates               , only : coordinateSpherical   , assignment(=)
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Galactic_Structure_Options, only : componentTypeBlackHole
    implicit none
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: radius                , radiusEscape     , &
         &                                                                    velocityRecoil
    logical                                                , intent(in   ) :: ignoreCentralBlackHole
    class           (massDistributionClass                ), pointer       :: massDistribution_
    double precision                                                       :: potential             , potentialSelf
    type            (coordinateSpherical                  )                :: coordinates           , coordinatesVirial
    
    ! Compute relevant potentials.
    coordinates       =  [radius      ,0.0d0,0.0d0]
    coordinatesVirial =  [radiusEscape,0.0d0,0.0d0]
    massDistribution_ => node             %massDistribution   (                             )
    potential         =  massDistribution_%potentialDifference(coordinates,coordinatesVirial)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    if (ignoreCentralBlackHole) then
       ! Compute potential of central black hole to be subtracted off of total value.
       massDistribution_ => node             %massDistribution   (componentType=componentTypeBlackHole                  )
       potentialSelf     =  massDistribution_%potentialDifference(              coordinates           ,coordinatesVirial)
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    else
       ! No correction for central black hole as it is to be included.
       potentialSelf=0.0d0
    end if
    ! Evaluate the escape condition.
    escapes=                          &
         &   +0.5d0*velocityRecoil**2 &
         &   +      potential         &
         &  >                         &
         &   +      potentialSelf
    return
  end function blackHolesRecoilEscapes

end module Nodes_Operators_Black_Holes_Utilities
