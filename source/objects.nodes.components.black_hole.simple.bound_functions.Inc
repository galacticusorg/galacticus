!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains custom functions for the simple black hole component.
!!}

double precision function Node_Component_Black_Hole_Simple_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the mass within a given radius for a central black hole. Black hole is treated as a point mass.
  !!}
  use :: Galactic_Structure_Options, only : componentTypeAll, componentTypeBlackHole, massTypeAll, massTypeBlackHole, &
          &                                 weightByMass    , massTypeGalactic
  implicit none
  class           (nodeComponentBlackHoleSimple), intent(inout) :: self
  integer                                       , intent(in   ) :: componentType, massType, weightBy, weightIndex
  double precision                              , intent(in   ) :: radius
  !$GLC attributes unused :: weightIndex

  ! Set zero enclosed mass by default.
  Node_Component_Black_Hole_Simple_Enclosed_Mass=0.0d0
  ! Return the black hole mass only if massType and componentType are of black hole type.
  if (.not.(componentType == componentTypeAll                                   .or. componentType == componentTypeBlackHole)) return
  if (.not.(massType      == massTypeAll      .or. massType == massTypeGalactic .or. massType      == massTypeBlackHole     )) return
  if (.not.(weightBy      == weightByMass                                                                                   )) return
  if (      radius        <  0.0d0                                                                                           ) return
  ! Return the mass of the black hole.
  Node_Component_Black_Hole_Simple_Enclosed_Mass=self%mass()
  return
end function Node_Component_Black_Hole_Simple_Enclosed_Mass

 function Node_Component_Black_Hole_Simple_Acceleration(self,positionCartesian,componentType,massType)
  !!{
  Computes the acceleration due to a central black hole. Black hole is treated as a point mass.
  !!}
   use :: Galactic_Structure_Options     , only : componentTypeAll               , componentTypeBlackHole, massTypeAll, massTypeBlackHole, &
        &                                         massTypeGalactic
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus, gigaYear              , megaParsec
  use :: Numerical_Constants_Prefixes    , only : kilo
  implicit none
  double precision                                             , dimension(3) :: Node_Component_Black_Hole_Simple_Acceleration
  class           (nodeComponentBlackHoleSimple), intent(inout)               :: self
  integer                                       , intent(in   )               :: componentType                                , massType
  double precision                              , intent(in   ), dimension(3) :: positionCartesian
  double precision                                                            :: radius
  
  ! Set zero enclosed acceleration by default.
  Node_Component_Black_Hole_Simple_Acceleration=0.0d0
  ! Return the black hole mass only if massType and componentType are of black hole type.
  if (.not.(componentType == componentTypeAll                                   .or. componentType == componentTypeBlackHole)) return
  if (.not.(massType      == massTypeAll      .or. massType == massTypeGalactic .or. massType      == massTypeBlackHole     )) return
  radius=sqrt(sum(positionCartesian**2))
  if (      radius        <  0.0d0                                                                                           ) return
  ! Compute the acceleration.
  Node_Component_Black_Hole_Simple_Acceleration=-kilo                            &
         &                                      *gigaYear                        &
         &                                      /megaParsec                      &
         &                                      *gravitationalConstantGalacticus &
         &                                      *self%mass()                     &
         &                                      *positionCartesian               &
         &                                      /radius**3
  return
end function Node_Component_Black_Hole_Simple_Acceleration

function Node_Component_Black_Hole_Simple_Tidal_Tensor(self,positionCartesian,componentType,massType)
  !!{
  Computes the acceleration due to a central black hole. Black hole is treated as a point mass.
  !!}
  use :: Galactic_Structure_Options      , only : componentTypeAll               , componentTypeBlackHole, massTypeAll          , massTypeBlackHole, &
       &                                          massTypeGalactic
  use :: Numerical_Constants_Math        , only : Pi
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  use :: Tensors                         , only : tensorRank2Dimension3Symmetric , tensorNullR2D3Sym     , tensorIdentityR2D3Sym, assignment(=)    , &
       &                                         operator(*)  
  use :: Vectors                         , only : Vector_Outer_Product
  implicit none
  type            (tensorRank2Dimension3Symmetric)                              :: Node_Component_Black_Hole_Simple_Tidal_Tensor
  class           (nodeComponentBlackHoleSimple  ), intent(inout)               :: self
  integer                                         , intent(in   )               :: componentType                                , massType
  double precision                                , intent(in   ), dimension(3) :: positionCartesian
  double precision                                                              :: radius
  type            (tensorRank2Dimension3Symmetric)                              :: positionTensor

  ! Set zero enclosed acceleration by default.
  Node_Component_Black_Hole_Simple_Tidal_Tensor=tensorNullR2D3Sym
  ! Return the black hole mass only if massType and componentType are of black hole type.
  if (.not.(componentType == componentTypeAll                                   .or. componentType == componentTypeBlackHole)) return
  if (.not.(massType      == massTypeAll      .or. massType == massTypeGalactic .or. massType      == massTypeBlackHole     )) return
  radius=sqrt(sum(positionCartesian**2))
  if (      radius        <  0.0d0                                                                                           ) return
  ! Compute the tidal tensor.
  if (self%radialPosition() <= 0.0d0) then
     positionTensor                               = Vector_Outer_Product(positionCartesian,symmetrize=.true.)
     Node_Component_Black_Hole_Simple_Tidal_Tensor=+gravitationalConstantGalacticus                       &
          &                                        *(                                                     &
          &                                          -(self%mass()      /radius**3)*tensorIdentityR2D3Sym &
          &                                          +(self%mass()*3.0d0/radius**5)*positionTensor        &
          &                                         )
  end if
  return
end function Node_Component_Black_Hole_Simple_Tidal_Tensor
