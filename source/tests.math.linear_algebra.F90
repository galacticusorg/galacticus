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
Contains a program to test linear algebra functions.
!!}

program Test_Math_Linear_Algebra
  !!{
  Tests of linear algebra functions.
  !!}
  use            :: Display                 , only : displayVerbositySet, verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Linear_Algebra          , only : assignment(=)      , matrix                , matrixLU            , matrixRotation   , &
          &                                          operator(*)        , vector
  use            :: Numerical_Constants_Math, only : Pi
  use            :: Sorting                 , only : sortIndex
  use            :: Unit_Tests              , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (vector  ), allocatable    :: vector1          , vector2          , &
       &                                        vector3          , vectorE
  type            (matrix  ), allocatable    :: matrix_          , matrixI          , &
       &                                        matrixT          , matrixP          , &
       &                                        matrixC          , matrixCD         , &
       &                                        matrixE          , matrixR          , &
       &                                        matrix1
  type            (matrixLU), allocatable    :: matrixLU_
  double precision          , dimension(3  ) :: vectorComponents
  integer         (c_size_t), dimension(3  ) :: vectorOrder
  type            (vector  ), dimension(3  ) :: vectors          , vectorsRotated
  double precision          , dimension(3,3) :: matrixComponents , matrixComponentsT
  double precision          , dimension(3,4) :: matrixComponents1
  double precision                           :: angle

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Math: linear algebra")
  !! Build the vectors.
  allocate(vector1)
  allocate(vector2)
  allocate(vectorE)
  vector1=vector([+4.0d0,2.0d0, 3.0d0])
  vector2=vector([-1.0d0,5.0d0,11.0d0])
  !! Build the matrix.
  allocate(matrix_)
  allocate(matrixI)
  allocate(matrixT)
  allocate(matrixP)
  allocate(matrixC)
  allocate(matrixE)
  allocate(matrixR)
  allocate(matrixCD)
  allocate(matrix1)
  matrix_=matrix(reshape([6.0d0,1.0d0,1.0d0,4.0d0,8.0d0,5.0d0,2.0d0,8.0d0,7.0d0                  ],[3,3]))
  matrix1=matrix(reshape([1.0d0,2.0d0,2.0d0,3.0d0,1.0d0,1.0d0,5.0d0,6.0d0,7.0d0,2.0d0,3.0d0,4.0d0],[3,4]))
  matrixC=matrix(reshape([1.0d0,1.1d0,1.4d0,1.1d0,2.0d0,1.9d0,1.4d0,1.9d0,3.0d0                  ],[3,3]))
  ! Tests of vector operations.
  call Unit_Tests_Begin_Group("Vector operations")
  !! Test vector addition.
  allocate(vector3)
  vector3         =vector1+vector2
  vectorComponents=vector3
  deallocate(vector3)
  call Assert("vector addition",vectorComponents,[3.0d0,7.0d0,14.0d0],relTol=1.0d-6)
  !! Test vector subtraction.
  allocate(vector3)
  vector3         =vector1-vector2
  vectorComponents=vector3
  deallocate(vector3)
  call Assert("vector subtraction",vectorComponents,[5.0d0,-3.0d0,-8.0d0],relTol=1.0d-6)
  !! Test vector magnitude.
  call Assert("vector magnitude",vector1%magnitude(),29.0d0)
  !! Test vector dot-product.
  call Assert("vector dot product",vector1.dot.vector2,39.0d0)
  !! Test vector cross-product.
  allocate(vector3)
  vector3         =vector1.cross.vector2
  vectorComponents=vector3
  deallocate(vector3)
  call Assert("vector cross product",vectorComponents,[7.0d0,-47.0d0,22.0d0],relTol=1.0d-6)
  call Unit_Tests_End_Group()
  ! Tests of matrix operations.
  call Unit_Tests_Begin_Group("Matrix operations")
  !! Test matrix-vector multiplication.
  deallocate(vector2)
  allocate(vector2)
  vector2         =matrix_*vector1
  vectorComponents=vector2
  deallocate(vector2)
  call Assert("matrix-vector product",vectorComponents,[38.0d0,44.0d0,35.0d0])
  !! Test matrix-matrix multiplication.
  matrixP          =matrix_*matrix1
  matrixComponents1=matrixP
  call Assert("matrix-matrix product",matrixComponents1,reshape([18.0d0,33.0d0,25.0d0,24.0d0,19.0d0,15.0d0,68.0d0,109.0d0,84.0d0,32.0d0,58.0d0,45.0d0],[3,4]))  
  !! Test determinant functions. Determinant computed using Mathematica.
  call Assert("log determinant" ,matrix_%logarithmicDeterminant(),log(94.0d0),absTol=1.0d-6)
  call Assert("determinant"     ,matrix_%           determinant(),    94.0d0 ,absTol=1.0d-6)
  call Assert("determinant sign",matrix_%       signDeterminant(),     1                   )
  !! Test transpose function.
  matrixT          =matrix_%transpose()
  matrixComponents =matrixT
  matrixComponentsT=matrix_
  matrixComponentsT=transpose(matrixComponentsT)
  call Assert("matrix transpose",matrixComponents,matrixComponentsT)
  !! Test inverse function. Inverse computed using Mathematica.
  matrixI         =matrix_%inverse()
  matrixComponents=matrixI
  call Assert("matrix inverse",matrixComponents,reshape([8.0d0/47.0d0,1.0d0/94.0d0,-3.0d0/94.0d0,-9.0d0/47.0d0,20.0d0/47.0d0,-13.0d0/47.0d0,8.0d0/47.0d0,-23.0d0/47.0d0,22.0d0/47.0d0],[3,3]),relTol=1.0d-3)
  !! Test covarince product function. Inverse computed using Mathematica.
  call Assert("covariance product",matrixC%covarianceProduct(vector1),25.8815d0,relTol=1.0d-3)
  !! Test eigensystem. Values computed using Mathematica.
  call matrixC%eigenSystem(matrixE,vectorE)
  vectorComponents=vectorE
  matrixComponents=matrixE
  vectorOrder=sortIndex(vectorComponents)
  call Assert("matrix eigenvalues",[vectorComponents(vectorOrder(3)),vectorComponents(vectorOrder(2)),vectorComponents(vectorOrder(1))],[5.21645d0,0.5361d0,0.247449d0],relTol=1.0d-3)
  call Assert(                                                                                                                   &
       &      "matrix eigenvectors",                                                                                             &
       &      reshape(                                                                                                           &
       &              [                                                                                                          &
       &               matrixComponents(1,vectorOrder(3)),matrixComponents(2,vectorOrder(3)),matrixComponents(3,vectorOrder(3)), &
       &               matrixComponents(1,vectorOrder(2)),matrixComponents(2,vectorOrder(2)),matrixComponents(3,vectorOrder(2)), &
       &               matrixComponents(1,vectorOrder(1)),matrixComponents(2,vectorOrder(1)),matrixComponents(3,vectorOrder(1))  &
       &              ]                                                                                                        , &
       &              [3,3]                                                                                                      &
       &             )                                                                                                         , &
       &      reshape(                                                                                                           &
       &              [                                                                                                          &
       &               +0.3889940d0,+0.563531d0,+0.728778d0                                                                    , &
       &               -0.0478684d0,-0.777650d0,+0.626872d0                                                                    , &
       &               +0.9199960d0,-0.278735d0,-0.275526d0                                                                      &
       &              ]                                                                                                        , &
       &              [3,3]                                                                                                      &
       &             )                                                                                                        ,  &
       &      relTol=1.0d-3                                                                                                      &
       &     )
  !! Test Cholesky decomposition. Values computed using Mathematica.
  matrixCD=matrix(matrixC)
  call matrixCD%choleskyDecomposition()
  matrixComponents=matrixCD
  call Assert("matrix Cholesky decomposition",matrixComponents,reshape([1.0d0, 1.1d0, 1.4d0, 1.1d0, 0.888819d0, 0.405032d0, 1.4d0,0.405032d0, 0.935922d0],[3,3]),relTol=1.0d-3)
  !! Test linear system solve. Values computed using Mathematica.
  allocate(vector3)
  vector3         =matrixC%linearSystemSolve(vector1)
  vectorComponents=vector3
  deallocate(vector3)
  call Assert("linear system solve",vectorComponents,[8.88728d0,-2.25434d0,-1.71965d0],relTol=1.0d-3)
  !! Test linear system solve using an LU matrix. Values computed using Mathematica.
  allocate(vector3  )
  allocate(matrixLU_)
  matrixLU_       =matrixLU                   (matrixC)
  vector3         =matrixLU_%squareSystemSolve(vector1)
  vectorComponents=vector3
  deallocate(vector3  )
  deallocate(matrixLU_)
  call Assert("LU square system solve",vectorComponents,[8.88728d0,-2.25434d0,-1.71965d0],relTol=1.0d-3)
  !! Done.
  deallocate(vector1 )
  deallocate(vectorE )
  deallocate(matrix_ )
  deallocate(matrixI )
  deallocate(matrixT )
  deallocate(matrixP )
  deallocate(matrixC )
  deallocate(matrixE )
  deallocate(matrixCD)
  deallocate(matrix1 )
  call Unit_Tests_End_Group()
  ! Tests of geometrical transformations.
  call Unit_Tests_Begin_Group("Geometrical transformations")
  ! Test rotation matrix construction.
  vectors       (1)=vector([+1.0d0,+0.0d0,+0.0d0])
  vectors       (2)=vector([+0.0d0,+1.0d0,+0.0d0])
  vectors       (3)=vector([+0.0d0,+0.0d0,+1.0d0])
  vectorsRotated(1)=vector([+0.0d0,+1.0d0,+0.0d0])
  vectorsRotated(2)=vector([-1.0d0,+0.0d0,+0.0d0])
  vectorsRotated(3)=vector([+0.0d0,+0.0d0,+1.0d0])
  angle            =-Pi/2.0d0
  matrixR=matrixRotation(vectors,vectorsRotated)
  matrixComponents=matrixR
  call Assert("rotation matrix construct (-π/2, around z-axis)", &
       &      matrixComponents                                 , &
       &      reshape(                                           &
       &              [                                          &
       &               +cos(angle),-sin(angle),+0.0d0,           &
       &               +sin(angle),+cos(angle),+0.0d0,           &
       &               +    0.0d0 ,+    0.0d0 ,+1.0d0            &
       &              ]                                        , &
       &              [3,3]                                      &
       &             )                                         , &
       &      absTol=1.0d-9                                    , &
       &      relTol=1.0d-9                                      &
       &     )
  deallocate(matrixR)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Math_Linear_Algebra
