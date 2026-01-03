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
Contains a module which implements calculations of vectors.
!!}

module Vectors
  !!{
  Implements calculations of vectors.
  !!}
  implicit none
  private
  public :: Vector_Magnitude, Vector_Product, Vector_Outer_Product, Vector_Outer_Product_Accumulate, Matrix_Copy_Upper_To_Lower_Triangle, Vector_Matrix_Multiply

  interface Vector_Outer_Product
     module procedure Vector_Outer_Product_Distinct
     module procedure Vector_Outer_Product_Self
  end interface Vector_Outer_Product

  interface Vector_Outer_Product_Accumulate
     module procedure Vector_Outer_Product_Accumulate_Self
  end interface Vector_Outer_Product_Accumulate

contains

  pure double precision function Vector_Magnitude(vector1)
    !!{
    Computes the magnitude of {\normalfont \ttfamily vector1}.
    !!}
    implicit none
    double precision, dimension(3), intent(in   ) :: vector1

    Vector_Magnitude=sqrt(sum(vector1**2))
    return
  end function Vector_Magnitude

  pure function Vector_Product(vector1,vector2) result(vector3)
    !!{
    Computes the vector product of {\normalfont \ttfamily vector1} and {\normalfont \ttfamily vector2}.
    !!}
    implicit none
    double precision, dimension(3)                :: vector3
    double precision, dimension(3), intent(in   ) :: vector1, vector2

    vector3(1)=vector1(2)*vector2(3)-vector1(3)*vector2(2)
    vector3(2)=vector1(3)*vector2(1)-vector1(1)*vector2(3)
    vector3(3)=vector1(1)*vector2(2)-vector1(2)*vector2(1)
    return
  end function Vector_Product

 function Vector_Outer_Product_Distinct(vector1,vector2)
    !!{
    Returns the outer product of two vectors.
    !!}
    implicit none
    double precision, dimension(:                          ), intent(in   ) :: vector1, vector2
    double precision, dimension(size(vector1),size(vector2))                :: Vector_Outer_Product_Distinct

    Vector_Outer_Product_Distinct=0.0d0
    ! Call the appropriate BLAS routine.
    call dger(size(vector1),size(vector2),1.0d0,vector1,1,vector2,1,Vector_Outer_Product_Distinct,size(vector1))
    return
  end function Vector_Outer_Product_Distinct

 function Vector_Outer_Product_Self(vector1,symmetrize)
    !!{
    Returns the outer product of a vector with itself.
    !!}
    implicit none
    double precision, dimension(:                          ), intent(in   ) :: vector1
    double precision, dimension(size(vector1),size(vector1))                :: Vector_Outer_Product_Self
    logical         , optional                              , intent(in   ) :: symmetrize

    Vector_Outer_Product_Self=0.0d0
    ! Call the appropriate BLAS routine.
    call dsyr("u",size(vector1),1.0d0,vector1,1,Vector_Outer_Product_Self,size(vector1))
    if (present(symmetrize).and.symmetrize) Vector_Outer_Product_Self=Matrix_Copy_Upper_To_Lower_Triangle(Vector_Outer_Product_Self)
    return
  end function Vector_Outer_Product_Self

  subroutine Vector_Outer_Product_Accumulate_Self(vector1,matrix,symmetrize,sparse)
    !!{
    Compute the outer product of a vector with itself and accumulate it to the given
    matrix. Compute only the upper triangle unless the symmetrize option is set to true. If
    the sparse option is set to true, assume a sparse matrix and accumulate only non-zero
    terms.
    !!}
    implicit none
    double precision, dimension(:  ), intent(in   ) :: vector1
    double precision, dimension(:,:), intent(inout) :: matrix
    logical         , optional      , intent(in   ) :: symmetrize, sparse
    integer         , dimension(:  ), allocatable   :: nonZero
    integer                                         :: i         , j

    if (.not.present(sparse).or..not.sparse) then
       ! Non-sparse calculation. Call the appropriate BLAS routine.
       call dsyr("u",size(vector1),1.0d0,vector1,1,matrix,size(vector1))
       ! Symmetrize if necessary.
       if (present(symmetrize).and.symmetrize) matrix=Matrix_Copy_Upper_To_Lower_Triangle(matrix)
    else
       ! Sparse calculation.
       ! Find non-zero elements of the vector.
       nonZero=pack([(i,i=1,size(vector1))],vector1 /= 0.0d0)
       if (size(nonZero) > 0) then
          ! Upper triangle.
          do i=1,size(nonZero)
             do j=i,size(nonZero)
                matrix(nonZero(i),nonZero(j))=matrix(nonZero(i),nonZero(j))+vector1(nonZero(i))*vector1(nonZero(j))
             end do
          end do
          ! Lower triangle (if necessary).
          if (present(symmetrize).and.symmetrize) then
             do i=2,size(nonZero)
                do j=1,i-1
                   matrix(nonZero(i),nonZero(j))=matrix(nonZero(i),nonZero(j))+vector1(nonZero(i))*vector1(nonZero(j))
                end do
             end do
          end if
       end if
    end if
    return
  end subroutine Vector_Outer_Product_Accumulate_Self

  function Matrix_Copy_Upper_To_Lower_Triangle(matrix)
    !!{
    Copies the upper triangle of a square matrix to the lower triangle.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision, dimension(:                 ,:                 ), intent(in   ) :: matrix
    double precision, dimension(size(matrix,dim=1),size(matrix,dim=2))                :: Matrix_Copy_Upper_To_Lower_Triangle
    integer                                                                           :: i                                  , j

    if (size(matrix,dim=1) /= size(matrix,dim=2)) call Error_Report('matrix must be square'//{introspection:location})
    do i=1,size(matrix,dim=1)
       do j=i,size(matrix,dim=2)
          Matrix_Copy_Upper_To_Lower_Triangle(i,j)=matrix(i,j)
          Matrix_Copy_Upper_To_Lower_Triangle(j,i)=matrix(i,j)
       end do
    end do
    return
  end function Matrix_Copy_Upper_To_Lower_Triangle

 function Vector_Matrix_Multiply(vector,matrix)
    !!{
    Returns the product of a vector with a matrix.
    !!}
    implicit none
    double precision, dimension(:                   ), intent(in   ) :: vector
    double precision, dimension(:,:                 ), intent(in   ) :: matrix
    double precision, dimension(  size(matrix,dim=2))                :: Vector_Matrix_Multiply

    ! Call the appropriate BLAS routine.
    call dgemv("t",size(matrix,dim=1),size(matrix,dim=2),1.0d0,matrix,size(matrix,dim=1),vector,1,0.0d0,Vector_Matrix_Multiply,1)
    return
  end function Vector_Matrix_Multiply

end module Vectors
