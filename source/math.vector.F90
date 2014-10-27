!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  public :: Vector_Magnitude, Vector_Product, Vector_Outer_Product

  interface Vector_Outer_Product
     module procedure Vector_Outer_Product_Distinct
     module procedure Vector_Outer_Product_Self
  end interface Vector_Outer_Product

contains

  pure double precision function Vector_Magnitude(vector1)
    !% Computes the magnitude of {\tt vector1}.
    implicit none
    double precision, dimension(3), intent(in   ) :: vector1

    Vector_Magnitude=sqrt(sum(vector1**2))
    return
  end function Vector_Magnitude

  pure function Vector_Product(vector1,vector2) result(vector3)
    !% Computes the vector product of {\tt vector1} and {\tt vector2}.
    implicit none
    double precision, dimension(3)                :: vector3
    double precision, dimension(3), intent(in   ) :: vector1, vector2

    vector3(1)=vector1(2)*vector2(3)-vector1(3)*vector2(2)
    vector3(2)=vector1(3)*vector2(1)-vector1(1)*vector2(3)
    vector3(3)=vector1(1)*vector2(2)-vector1(2)*vector2(1)
    return
  end function Vector_Product

 function Vector_Outer_Product_Distinct(vector1,vector2)
    !% Returns the outer product of two vectors.
    implicit none
    double precision, dimension(:                          ), intent(in   ) :: vector1, vector2
    double precision, dimension(size(vector1),size(vector2))                :: Vector_Outer_Product_Distinct

    Vector_Outer_Product_Distinct=0.0d0
    ! Call the appropriate BLAS routine.
    call dger(size(vector1),size(vector2),1.0d0,vector1,1,vector2,1,Vector_Outer_Product_Distinct,size(vector1))
    return
  end function Vector_Outer_Product_Distinct

 function Vector_Outer_Product_Self(vector1,symmetrize)
    !% Returns the outer product of a vector with itself.
    implicit none
    double precision, dimension(:                          ), intent(in   ) :: vector1
    double precision, dimension(size(vector1),size(vector1))                :: Vector_Outer_Product_Self
    logical         , optional                              , intent(in   ) :: symmetrize
    integer                                                                 :: i                        , j

    Vector_Outer_Product_Self=0.0d0
    ! Call the appropriate BLAS routine.
    call dsyr("u",size(vector1),1.0d0,vector1,1,Vector_Outer_Product_Self,size(vector1))
    if (present(symmetrize).and.symmetrize) then
       forall(i=1:size(vector1)-1)
          forall(j=i+1:size(vector1))
             Vector_Outer_Product_Self(j,i)=Vector_Outer_Product_Self(i,j)
          end forall
       end forall
    end if
    return
  end function Vector_Outer_Product_Self

end module Vectors
