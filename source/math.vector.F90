!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of vectors.

module Vectors
  !% Implements calculations of vectors.
  implicit none
  private
  public :: Vector_Magnitude, Vector_Product, Vector_Self_Outer_Product, Vector_Outer_Product

contains

  pure double precision function Vector_Magnitude(vector1)
    !% Computes the magnitude of {\normalfont \ttfamily vector1}.
    implicit none
    double precision, dimension(3), intent(in   ) :: vector1

    Vector_Magnitude=sqrt(sum(vector1**2))
    return
  end function Vector_Magnitude

  pure function Vector_Product(vector1,vector2) result(vector3)
    !% Computes the vector product of {\normalfont \ttfamily vector1} and {\normalfont \ttfamily vector2}.
    implicit none
    double precision, dimension(3)                :: vector3
    double precision, dimension(3), intent(in   ) :: vector1, vector2

    vector3(1)=vector1(2)*vector2(3)-vector1(3)*vector2(2)
    vector3(2)=vector1(3)*vector2(1)-vector1(1)*vector2(3)
    vector3(3)=vector1(1)*vector2(2)-vector1(2)*vector2(1)
    return
  end function Vector_Product

  pure function Vector_Outer_Product(vector1,vector2)
    !% Returns the outer product of two vectors.
    implicit none
    double precision, dimension(:                          ), intent(in   ) :: vector1, vector2
    double precision, dimension(size(vector1),size(vector2))                :: Vector_Outer_Product
    integer                                                                 :: i

    forall(i=1:size(vector1)) Vector_Outer_Product(i,:)=vector1(i)*vector2(:)
    return
  end function Vector_Outer_Product
  
  function Vector_Self_Outer_Product(vector)
    !% Find the outer product of a vector with itself, returning the result as a rank 2, 3D, symmetric tensor.
    use Tensors
    implicit none
    type            (tensorRank2Dimension3Symmetric)                :: Vector_Self_Outer_Product
    double precision                                , intent(in   ) :: vector(3)
    double precision                                                :: matrix(3,3)
    integer                                                         :: i,j

    forall(i=1:3) 
       forall(j=1:3)
          matrix(i,j)=vector(i)*vector(j)
       end forall
    end forall
    call Vector_Self_Outer_Product%fromMatrix(matrix)
    return
  end function Vector_Self_Outer_Product

end module Vectors
