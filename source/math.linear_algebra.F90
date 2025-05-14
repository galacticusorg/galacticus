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
Contains a module which implements linear algebra calculations.
!!}

module Linear_Algebra
  !!{
  Implements linear algebra calculations.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_ptr, c_double, c_size_t, c_int
  implicit none
  private
  public :: vector        , matrix         , matrixRotation      , matrixRotationPlusTranslation, &
       &    matrixLU      , assignment(=)  , operator(*)         , gsl_vector_get               , &
       &    gsl_vector_set, gsl_vector_free, matrixRotationRandom

  type, public :: vector
     !!{
     Vector class.
     !!}
     private
     type   (c_ptr   ), allocatable :: vector_
     integer(c_size_t)              :: size_
   contains
     !![
     <methods>
       <method description="Compute the magnitude of a vector."                                                method="magnitude"        />
       <method description="Compute {\normalfont \ttfamily vector1} $\cdot$ {\normalfont \ttfamily vector2}."  method="operator(.dot.)"  />
       <method description="Compute {\normalfont \ttfamily vector1}-{\normalfont \ttfamily vector2}."          method="operator(-)"      />
       <method description="Compute {\normalfont \ttfamily vector1}+{\normalfont \ttfamily vector2}."          method="operator(+)"      />
       <method description="Compute {\normalfont \ttfamily vector1} $\times$ {\normalfont \ttfamily vector2}." method="operator(.cross.)"/>
       <method description="Return a C pointer to the GSL vector object."                                      method="gslObject"        />
     </methods>
     !!]
     final     ::                        vectorDestructorRank0, vectorDestructorRank1
     procedure :: magnitude           => vectorMagnitude
     procedure ::                        vectorDotProduct
     generic   :: operator  (.dot.  ) => vectorDotProduct
     procedure ::                        vectorSubtract
     generic   :: operator  (-      ) => vectorSubtract
     procedure ::                        vectorAdd
     generic   :: operator  (+      ) => vectorAdd
     procedure ::                        vectorCrossProduct
     generic   :: operator  (.cross.) => vectorCrossProduct
     procedure :: gslObject           => vectorGSLObject
  end type vector

  interface vector
     !!{
     Interface to vector constructors.
     !!}
     module procedure vectorConstructor
     module procedure vectorZeroConstructor
     module procedure vectorCopyConstructor
  end interface vector
  
  type, public :: matrix
     !!{
     Matrix class.
     !!}
     private
     type   (c_ptr   ), allocatable  :: matrix_
     integer(c_size_t), dimension(2) :: size_
     logical                         :: nonZeroRowColumnsChecked          , isSquare, &
          &                             hasZeroRowColumns
     type   (matrix  ), pointer      :: LUdecomposition          => null()
     type   (c_ptr   ), pointer      :: LUpermutation            => null()
   contains
     !![
     <methods>
       <method description="Compute the product of two matrices."                                                              method="operator(*)"           />
       <method description="Compute the sum of two matrices."                                                                  method="operator(+)"           />
       <method description="Compute and return the determinant of the matrix."                                                 method="determinant"           />
       <method description="Compute and return the logarithm of the determinant of the matrix."                                method="logarithmicDeterminant"/>
       <method description="Compute and return the sign of the determinant of the matrix."                                     method="signDeterminant"       />
       <method description="Compute and return the matrix inverse."                                                            method="inverse"               />
       <method description="Return the transpose of a matrix."                                                                 method="transpose"             />
       <method description="Compute $y C^{-1} y^\mathrm{T}$ as appears in likelihood functions utilizing covariance matrices." method="covarianceProduct"     />
       <method description="Solve the linear system $y = A \cdot x$ where $A$ is ourself and $y$ is the specified vector."     method="linearSystemSolve"     />
       <method description="Compute eigenvectors and eigenvalues of the matrix."                                               method="eigenSystem"           />
       <method description="Compute the Cholesky decomposition of the matrix in place."                                        method="choleskyDecomposition" />
     </methods>
     !!]
     final     ::                           matrixDestructorRank0       , matrixDestructorRank1
     procedure ::                           matrixMatrixProduct
     procedure ::                           matrixMatrixAdd
     generic   :: operator(*)            => matrixMatrixProduct
     generic   :: operator(+)            => matrixMatrixAdd
     procedure :: determinant            => matrixDeterminant
     procedure :: logarithmicDeterminant => matrixLogarithmicDeterminant
     procedure :: signDeterminant        => matrixSignDeterminant
     procedure :: transpose              => matrixTranspose
     procedure :: inverse                => matrixInverse
     procedure :: linearSystemSolve      => matrixLinearSystemSolve
     procedure :: covarianceProduct      => matrixCovarianceProduct
     procedure :: eigenSystem            => matrixEigenSystem
     procedure :: choleskyDecomposition  => matrixCholeskyDecomposition
  end type matrix

  interface matrix
     !!{
     Interface to matrix constructors.
     !!}
     module procedure matrixConstructor
     module procedure matrixZeroConstructor
     module procedure matrixCopyConstructor
  end interface matrix
  
  type, public, extends(matrix) :: matrixLU
     !!{
     Matrix class for LU matrices.
     !!}
     private
     type   (c_ptr), allocatable :: permutation
     integer(c_int)              :: decompositionSign
   contains
     !![
     <methods>
       <method description="Solve the linear system $y = A \cdot x$ where $A$ is ourself and $y$ is the specified vector." method="squareSystemSolve" />
     </methods>
     !!]
     final     ::                      matrixLUDestructorRank0  , matrixLUDestructorRank1
     procedure :: squareSystemSolve => matrixLUSquareSystemSolve
  end type matrixLU

  interface matrixLU
     !!{
     Interface to LU matrix constructors.
     !!}
     module procedure matrixLUConstructor
  end interface matrixLU
  
  ! Assignment interfaces.
  interface assignment(=)
     module procedure vectorAssignmentConstructor
     module procedure vectorUnassignment
     module procedure matrixAssignmentConstructor
     module procedure matrixUnassignment
     module procedure matrixLUUnassignment
  end interface assignment(=)

  ! Operator interfaces.
  interface operator(*)
     module procedure matrixVectorMultiply
     module procedure matrixScalarMultiply
  end interface operator(*)
  
  interface
     function gsl_vector_alloc(n) bind(c,name='gsl_vector_alloc')
       !!{
       Template for the GSL vector alloc function.
       !!}
       import c_ptr, c_size_t
       type   (c_ptr   )        :: gsl_vector_alloc
       integer(c_size_t), value :: n
     end function gsl_vector_alloc

     subroutine gsl_vector_free(v) bind(c,name='gsl_vector_free')
       !!{
       Template for the GSL vector free function.
       !!}
       import c_ptr, c_size_t
       type(c_ptr), value :: v
     end subroutine gsl_vector_free

     subroutine gsl_vector_set(v,i,x) bind(c,name='gsl_vector_set')
       !!{
       Template for the GSL vector set function.
       !!}
       import c_ptr, c_size_t, c_double
       type   (c_ptr   ), value :: v
       integer(c_size_t), value :: i
       real   (c_double), value :: x
     end subroutine gsl_vector_set

     subroutine gsl_vector_set_zero(v) bind(c,name='gsl_vector_set_zero')
       !!{
       Template for the GSL vector set zero function.
       !!}
       import c_ptr
       type(c_ptr), value :: v
     end subroutine gsl_vector_set_zero

     function gsl_vector_get(v,i) bind(c,name='gsl_vector_get')
       !!{
       Template for the GSL vector set function.
       !!}
       import c_ptr, c_size_t, c_double
       real   (c_double)        :: gsl_vector_get
       type   (c_ptr   ), value :: v
       integer(c_size_t), value :: i
     end function gsl_vector_get
     
     function gsl_vector_memcpy(dest,src) bind(c,name='gsl_vector_memcpy')
       !!{
       Template for the GSL vector copy function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_vector_memcpy
       type   (c_ptr), value :: dest, src
     end function gsl_vector_memcpy

     function gsl_vector_add(a,b) bind(c,name='gsl_vector_add')
       !!{
       Template for the GSL vector addition function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_vector_add
       type   (c_ptr), value :: a             , b
     end function gsl_vector_add

     function gsl_vector_sub(a,b) bind(c,name='gsl_vector_sub')
       !!{
       Template for the GSL vector subtract function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_vector_sub
       type   (c_ptr), value :: a             , b
     end function gsl_vector_sub

     function gsl_matrix_alloc(n1,n2) bind(c,name='gsl_matrix_alloc')
       !!{
       Template for the GSL matrix alloc function.
       !!}
       import c_ptr, c_size_t
       type   (c_ptr   )        :: gsl_matrix_alloc
       integer(c_size_t), value :: n1              , n2
     end function gsl_matrix_alloc

     subroutine gsl_matrix_free(m) bind(c,name='gsl_matrix_free')
       !!{
       Template for the GSL matrix free function.
       !!}
       import c_ptr, c_size_t
       type(c_ptr), value :: m
     end subroutine gsl_matrix_free

     subroutine gsl_matrix_set(m,i,j,x) bind(c,name='gsl_matrix_set')
       !!{
       Template for the GSL matrix set function.
       !!}
       import c_ptr, c_size_t, c_double
       type   (c_ptr   ), value :: m
       integer(c_size_t), value :: i, j
       real   (c_double), value :: x
     end subroutine gsl_matrix_set
     
     subroutine gsl_matrix_set_zero(m) bind(c,name='gsl_matrix_set_zero')
       !!{
       Template for the GSL matrix set zero function.
       !!}
       import c_ptr
       type(c_ptr), value :: m
     end subroutine gsl_matrix_set_zero

     function gsl_matrix_get(m,i,j) bind(c,name='gsl_matrix_get')
       !!{
       Template for the GSL matrix set function.
       !!}
       import c_ptr, c_size_t, c_double
       real   (c_double)        :: gsl_matrix_get
       type   (c_ptr   ), value :: m
       integer(c_size_t), value :: i             , j
     end function gsl_matrix_get
     
     function gsl_matrix_memcpy(dest,src) bind(c,name='gsl_matrix_memcpy')
       !!{
       Template for the GSL matrix copy function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_matrix_memcpy
       type   (c_ptr), value :: dest             , src
     end function gsl_matrix_memcpy

     function gsl_matrix_transpose_memcpy(dest,src) bind(c,name='gsl_matrix_transpose_memcpy')
       !!{
       Template for the GSL matrix transpose function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_matrix_transpose_memcpy
       type   (c_ptr), value :: dest                       , src
     end function gsl_matrix_transpose_memcpy

     function gsl_permutation_alloc(n) bind(c,name='gsl_permutation_alloc')
       !!{
       Template for the GSL permutation alloc function.
       !!}
       import c_ptr, c_size_t
       type   (c_ptr   )        :: gsl_permutation_alloc
       integer(c_size_t), value :: n
     end function gsl_permutation_alloc

     subroutine gsl_permutation_free(p) bind(c,name='gsl_permutation_free')
       !!{
       Template for the GSL permutation free function.
       !!}
       import c_ptr, c_size_t
       type(c_ptr), value :: p
     end subroutine gsl_permutation_free

     function gsl_matrix_scale(a,x) bind(c,name='gsl_matrix_scale')
       !!{
       Template for the GSL matrix scale function.
       !!}
       import c_ptr, c_int, c_double
       integer(c_int   )        :: gsl_matrix_scale
       type   (c_ptr   ), value :: a
       real   (c_double), value :: x
     end function gsl_matrix_scale

     function gsl_matrix_add(a,b) bind(c,name='gsl_matrix_add')
       !!{
       Template for the GSL matrix addition function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_matrix_add
       type   (c_ptr), value :: a             , b
     end function gsl_matrix_add

     function gsl_linalg_LU_decomp(A,p,signum) bind(c,name='gsl_linalg_LU_decomp')
       !!{
       Template for the GSL LU decomposition function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_linalg_LU_decomp
       type   (c_ptr), value :: A                   , p
       integer(c_int)        :: signum
     end function gsl_linalg_LU_decomp

     function gsl_linalg_LU_solve(LU,p,b,x) bind(c,name='gsl_linalg_LU_solve')
       !!{
       Template for the GSL LU decomposition function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_linalg_LU_solve
       type   (c_ptr), value :: LU                 , p, &
            &                   b                  , x
     end function gsl_linalg_LU_solve

     function gsl_linalg_LU_det(LU,signum) bind(c,name='gsl_linalg_LU_det')
       !!{
       Template for the GSL LU determinant function.
       !!}
       import c_ptr, c_double, c_int
       real   (c_double)        :: gsl_linalg_LU_det
       type   (c_ptr   ), value :: LU
       integer(c_int   ), value :: signum
     end function gsl_linalg_LU_det

     function gsl_linalg_LU_invert(LU,p,inverse) bind(c,name='gsl_linalg_LU_invert')
       !!{
       Template for the GSL LU invert function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_linalg_LU_invert
       type   (c_ptr), value :: LU                  , p, &
            &                   inverse
     end function gsl_linalg_LU_invert

     function gsl_linalg_LU_lndet(LU) bind(c,name='gsl_linalg_LU_lndet')
       !!{
       Template for the GSL LU logarithmic determinant function.
       !!}
       import c_ptr, c_double
       real(c_double)        :: gsl_linalg_LU_lndet
       type(c_ptr   ), value :: LU
     end function gsl_linalg_LU_lndet

     function gsl_linalg_LU_sgndet(LU,signum) bind(c,name='gsl_linalg_LU_sgndet')
       !!{
       Template for the GSL LU determinant sign function.
       !!}
       import c_ptr, c_int
       integer(c_int)           :: gsl_linalg_LU_sgndet
       type   (c_ptr   ), value :: LU
       integer(c_int   ), value :: signum
     end function gsl_linalg_LU_sgndet

     function gsl_linalg_cholesky_decomp(A) bind(c,name='gsl_linalg_cholesky_decomp')
       !!{
       Template for the GSL Cholesky decomposition function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_linalg_cholesky_decomp
       type   (c_ptr), value :: A
     end function gsl_linalg_cholesky_decomp

     function gsl_eigen_symmv_alloc(n) bind(c,name='gsl_eigen_symmv_alloc')
       !!{
       Template for the GSL symmetric eigenvalues workspace alloc function.
       !!}
       import c_ptr, c_size_t
       type   (c_ptr   )        :: gsl_eigen_symmv_alloc
       integer(c_size_t), value :: n
     end function gsl_eigen_symmv_alloc

     subroutine gsl_eigen_symmv_free(w) bind(c,name='gsl_eigen_symmv_free')
       !!{
       Template for the GSL symmetric eigenvalues workspace free function.
       !!}
       import c_ptr
       type(c_ptr), value :: w
     end subroutine gsl_eigen_symmv_free

     function gsl_eigen_symmv(A,eval,evec,w) bind(c,name='gsl_eigen_symmv')
       !!{
       Template for the GSL symmetric eigenvalue system function.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gsl_eigen_symmv
       type   (c_ptr), value :: A              , eval, &
            &                   evec           , w
     end function gsl_eigen_symmv

     function gsl_blas_ddot(x,y,res) bind(c,name='gsl_blas_ddot')
       !!{
       Template for the GSL BLAS vector dot product function.
       !!}
       import c_ptr, c_double ,c_int
       integer(c_int   )        :: gsl_blas_ddot
       type   (c_ptr   ), value :: x            , y
       real   (c_double)        :: res
     end function gsl_blas_ddot
     
     function gsl_blas_dgemv(TransA,alpha,A,x,beta,y) bind(c,name='gsl_blas_dgemv')
       !!{
       Template for the GSL BLAS matrix-vector multiply-add function.
       !!}
       import c_ptr, c_double, c_int
       integer(c_int             )        :: gsl_blas_dgemv
       integer(kind(CblasNoTrans)), value :: TransA
       real   (c_double          ), value :: alpha         , beta
       type   (c_ptr             ), value :: A             , x   , &
            &                                y
       !$GLC attributes interoperable :: TransA
     end function gsl_blas_dgemv

     function gsl_blas_dgemm(TransA,TransB,alpha,A,B,beta,C) bind(c,name='gsl_blas_dgemm')
       !!{
       Template for the GSL BLAS matrix-vector multiply-add function.
       !!}
       import c_ptr, c_double, c_int
       integer(c_int             )        :: gsl_blas_dgemm
       integer(kind(CblasNoTrans)), value :: TransA        , TransB
       real   (c_double          ), value :: alpha         , beta
       type   (c_ptr             ), value :: A             , B   , &
            &                                C
       !$GLC attributes interoperable :: TransA, TransB
     end function gsl_blas_dgemm
  end interface
  
  ! Extract enums needed for BLAS functions.
  !![
  <constant variable="CBLAS_Transpose" gslSymbol="CBLAS_TRANSPOSE" gslHeader="gsl_cblas" type="enum" members="CblasNoTrans, CblasTrans, CblasConjTrans" description="Enumeration of CBLAS transpose options." reference="Gnu Scientific Library" group="GSL"/>
  <constant variable="CBLAS_UpLo"      gslSymbol="CBLAS_UPLO"      gslHeader="gsl_cblas" type="enum" members="CblasUpper  , CblasLower"                 description="Enumeration of CBLAS matrix upper/lower options." reference="Gnu Scientific Library" group="GSL"/>
  <constant variable="CBLAS_Diag"      gslSymbol="CBLAS_DIAG"      gslHeader="gsl_cblas" type="enum" members="CblasNonUnit, CblasUnit"                  description="Enumeration of whether a CBLAS triangular matrix has unit diagonal or not." reference="Gnu Scientific Library" group="GSL"/>
  !!]
  
contains

  !! Vector functions.

  function vectorConstructor(array) result(self)
    !!{
    Constructor for {\normalfont \ttfamily vector} class which builds the vector from an array.
    !!}
    implicit none
    type            (vector  )                              :: self
    double precision          , intent(in   ), dimension(:) :: array
    integer         (c_size_t)                              :: i

    allocate(self%vector_)
    self%size_  =size(array,kind=c_size_t)
    self%vector_=gsl_vector_alloc(self%size_)
    do i=1,size(array,dim=1,kind=c_size_t)
       call gsl_vector_set(self%vector_,i-1_c_size_t,array(i))
    end do
    return
  end function vectorConstructor

  function vectorZeroConstructor(n) result(self)
    !!{
    Constructor for {\normalfont \ttfamily vector} class which builds the vector and initializes all elements to zero.
    !!}
    implicit none
    type   (vector  )                :: self
    integer(c_size_t), intent(in   ) :: n
    
    allocate(self%vector_)
    self%size_  =n
    self%vector_=gsl_vector_alloc(self%size_)
    call gsl_vector_set_zero(self%vector_)
    return
  end function vectorZeroConstructor

  function vectorCopyConstructor(source) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily vector} class which builds the vector by copying a {\normalfont \ttfamily source} vector.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (vector)             :: self
    type   (vector), intent(in) :: source
    integer(c_int )             :: status

    allocate(self%vector_)
    self%size_  =source%size_
    self%vector_=gsl_vector_alloc(self%size_)
    status      =gsl_vector_memcpy(self%vector_,source%vector_)
    if (status /= GSL_Success) call Error_Report('vector copy failed'//{introspection:location})
    return
  end function vectorCopyConstructor
  
  subroutine vectorAssignmentConstructor(self,array)
    !!{
    Constructor for {\normalfont \ttfamily vector} class which overloads the assignment operator.
    !!}
    implicit none
    type            (vector), intent(  out)               :: self
    double precision        , intent(in   ), dimension(:) :: array
    
    self=vector(array)
    return
  end subroutine vectorAssignmentConstructor
  
  subroutine vectorUnassignment(array,self)
    !!{
    Assign elements of a {\normalfont \ttfamily vector} class to an array
    !!}
    implicit none
    double precision          , intent(  out), dimension(:) :: array
    type            (vector  ), intent(in   )               :: self
    integer         (c_size_t)                              :: i

    do i=1_c_size_t,self%size_
       array(i)=gsl_vector_get(self%vector_,i-1_c_size_t)
    end do
    return
  end subroutine vectorUnassignment
  
  subroutine vectorDestructorRank0(self)
    !!{
    Rank-0 destructor for the {\normalfont \ttfamily vector} class.
    !!}
    implicit none
    type(vector), intent(inout) :: self
    
    if (allocated(self%vector_)) then
       call gsl_vector_free(self%vector_)
       deallocate(self%vector_)
    end if
    return
  end subroutine vectorDestructorRank0

  subroutine vectorDestructorRank1(self)
    !!{
    Rank-1 destructor for the {\normalfont \ttfamily vector} class.
    !!}
    implicit none
    type   (vector), intent(inout), dimension(:) :: self
    integer                                      :: i
    
    do i=1,size(self)
       if (allocated(self(i)%vector_)) then
          call gsl_vector_free(self(i)%vector_)
          deallocate(self(i)%vector_)
       end if
    end do
    return
  end subroutine vectorDestructorRank1

  function vectorAdd(vector1,vector2)
    !!{
    Add one vector to another.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (vector)                :: vectorAdd
    class  (vector), intent(in   ) :: vector1  , vector2
    integer(c_int )                :: status

    vectorAdd=vector(vector1)
    status   =gsl_vector_add(vectorAdd%vector_,vector2%vector_)
    if (status /= GSL_Success) call Error_Report('vector addition failed'//{introspection:location})
    return
  end function vectorAdd

  function vectorSubtract(vector1,vector2)
    !!{
    Subtract one vector from another.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (vector)                :: vectorSubtract
    class  (vector), intent(in   ) :: vector1       , vector2
    integer(c_int )                :: status

    vectorSubtract=vector(vector1)
    status        =gsl_vector_sub(vectorSubtract%vector_,vector2%vector_)
    if (status /= GSL_Success) call Error_Report('vector subtract failed'//{introspection:location})
    return
  end function vectorSubtract

  double precision function vectorMagnitude(vector_)
    !!{
    Compute the magnitude of a vector.
    !!}
    implicit none
    class(vector), intent(in   ) :: vector_

   vectorMagnitude=vector_.dot.vector_
    return
  end function vectorMagnitude
  
  double precision function vectorDotProduct(vector1,vector2)
    !!{
    Compute the dot product of two vectors.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class  (vector), intent(in   ) :: vector1, vector2
    integer(c_int )                :: status

    status=gsl_blas_ddot(vector1%vector_,vector2%vector_,vectorDotProduct)
    if (status /= GSL_Success) call Error_Report('vector dot product failed'//{introspection:location})
    return
  end function vectorDotProduct
  
  function vectorCrossProduct(vector1,vector2)
    !!{
    Compute the cross product of two vectors.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (vector)                :: vectorCrossProduct
    class(vector), intent(in   ) :: vector1           , vector2

    if (vector1%size_ /= 3_c_size_t .or. vector2%size_ /= 3_c_size_t) &
         & call Error_Report('vector cross product only defined for 3D vectors'//{introspection:location})
    vectorCrossProduct=vector(                                                                                         &
         &                    [                                                                                        &
         &                     +gsl_vector_get(vector1%vector_,1_c_size_t)*gsl_vector_get(vector2%vector_,2_c_size_t)  &
         &                     -gsl_vector_get(vector1%vector_,2_c_size_t)*gsl_vector_get(vector2%vector_,1_c_size_t), &
         &                     +gsl_vector_get(vector1%vector_,2_c_size_t)*gsl_vector_get(vector2%vector_,0_c_size_t)  &
         &                     -gsl_vector_get(vector1%vector_,0_c_size_t)*gsl_vector_get(vector2%vector_,2_c_size_t), &
         &                     +gsl_vector_get(vector1%vector_,0_c_size_t)*gsl_vector_get(vector2%vector_,1_c_size_t)  &
         &                     -gsl_vector_get(vector1%vector_,1_c_size_t)*gsl_vector_get(vector2%vector_,0_c_size_t)  &
         &                    ]                                                                                        &
         &                   )
    return
  end function vectorCrossProduct

  function vectorGSLObject(self)
    !!{
    Return a C pointer to the GSL vector object.
    !!}
    implicit none
    type (c_ptr )                :: vectorGSLObject
    class(vector), intent(in   ) :: self

    vectorGSLObject=self%vector_
    return
  end function vectorGSLObject
  
  !! Matrix functions.

  function matrixConstructor(array) result(self)
    !!{
    Constructor for {\normalfont \ttfamily matrix} class which builds the matrix from an array.
    !!}
    implicit none
    type            (matrix  )                                :: self
    double precision          , intent(in   ), dimension(:,:) :: array
    integer         (c_size_t)                                :: i    , j

    allocate(self%matrix_)
    self%size_                   =shape(array,kind=c_size_t)
    self%isSquare                =self%size_(1) == self%size_(2)
    self%matrix_                 =gsl_matrix_alloc(self%size_(1),self%size_(2))
    self%nonZeroRowColumnsChecked=.false.
    self%hasZeroRowColumns       =.false.
    do i=1,size(array,dim=1,kind=c_size_t)
       do j=1,size(array,dim=2,kind=c_size_t)
          call gsl_matrix_set(self%matrix_,i-1_c_size_t,j-1_c_size_t,array(i,j))
       end do
    end do
    return
  end function matrixConstructor

  function matrixZeroConstructor(n1,n2) result(self)
    !!{
    Constructor for {\normalfont \ttfamily matrix} class which builds the matrix and initializes all elements to zero.
    !!}
    implicit none
    type   (matrix  )                :: self
    integer(c_size_t), intent(in   ) :: n1  , n2

    allocate(self%matrix_)
    self%isSquare                =n1 == n2
    self%size_                   =[n1,n2]
    self%nonZeroRowColumnsChecked=.false.
    self%hasZeroRowColumns       =.false.
    self%matrix_                 =gsl_matrix_alloc(self%size_(1),self%size_(2))
    call gsl_matrix_set_zero(self%matrix_)
    return
  end function matrixZeroConstructor

  function matrixCopyConstructor(source) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily matrix} class which builds the matrix by copying a {\normalfont \ttfamily source} matrix.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (matrix)             :: self
    type   (matrix), intent(in) :: source
    integer(c_int )             :: status

    allocate(self%matrix_)
    self%size_                   =source%size_
    self%isSquare                =source%isSquare
    self%nonZeroRowColumnsChecked=source%nonZeroRowColumnsChecked
    self%hasZeroRowColumns       =source%hasZeroRowColumns
    self%matrix_                 =gsl_matrix_alloc(self%size_(1),self%size_(2))
    status                       =gsl_matrix_memcpy(self%matrix_,source%matrix_)
    if (status /= GSL_Success) call Error_Report('matrix copy failed'//{introspection:location})
    return
  end function matrixCopyConstructor

  subroutine matrixAssignmentConstructor(self,array)
    !!{
    Constructor for {\normalfont \ttfamily matrix} class which overloads the assignment operator.
    !!}
    implicit none
    type            (matrix), intent(  out)                 :: self
    double precision        , intent(in   ), dimension(:,:) :: array
    
    self=matrix(array)
    return
  end subroutine matrixAssignmentConstructor
  
  subroutine matrixUnassignment(array,self)
    !!{
    Assign elements of a {\normalfont \ttfamily matrix} class to an array
    !!}
    implicit none
    double precision          , intent(  out), dimension(:,:) :: array
    type            (matrix  ), intent(in   )                 :: self
    integer         (c_size_t)                                :: i    , j

    do i=1_c_size_t,self%size_(1)
       do j=1_c_size_t,self%size_(2)
          array(i,j)=gsl_matrix_get(self%matrix_,i-1_c_size_t,j-1_c_size_t)
       end do
    end do
    return
  end subroutine matrixUnassignment
  
  subroutine matrixDestructorRank0(self)
    !!{
    Destructor for the {\normalfont \ttfamily matrix} class
    !!}
    implicit none
    type(matrix), intent(inout) :: self

    if (allocated(self%matrix_)) then
       call gsl_matrix_free(self%matrix_)
       deallocate(self%matrix_)
    end if
    if (associated(self%LUdecomposition)) then
       call gsl_permutation_free(self%LUpermutation)
       deallocate(self%LUdecomposition)
       deallocate(self%LUpermutation  )
    end if
    return
  end subroutine matrixDestructorRank0

  subroutine matrixDestructorRank1(self)
    !!{
    Destructor for the {\normalfont \ttfamily matrix} class
    !!}
    implicit none
    type   (matrix), intent(inout), dimension(:) :: self
    integer                                      :: i

    do i=1,size(self)
       if (allocated(self(i)%matrix_)) then
          call gsl_matrix_free(self(i)%matrix_)
          deallocate(self(i)%matrix_)
       end if
    end do
    return
  end subroutine matrixDestructorRank1

  function matrixMatrixProduct(matrix1,matrix2)
    !!{
    Multiply a matrix by a matrix, returning a matrix.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (matrix)                :: matrixMatrixProduct
    class  (matrix), intent(in   ) :: matrix1             , matrix2
    integer(c_int )                :: status

    if (matrix1%size_(2) /= matrix2%size_(1)) call Error_Report('matrices can not be multiplied'//{introspection:location})
    matrixMatrixProduct=matrix(matrix1%size_(1),matrix2%size_(2))
    status             =gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0d0,matrix1%matrix_,matrix2%matrix_,0.0d0,matrixMatrixProduct%matrix_)
    if (status /= GSL_Success) call Error_Report('matrix-matrix multiply failed'//{introspection:location})
    return
  end function matrixMatrixProduct
  
  function matrixMatrixAdd(matrix1,matrix2)
    !!{
    Sum two matrices, returning a matrix.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (matrix)                :: matrixMatrixAdd
    class  (matrix), intent(in   ) :: matrix1        , matrix2
    integer(c_int )                :: status

    if (any(matrix1%size_ /= matrix2%size_)) call Error_Report('matrices can not be summed'//{introspection:location})
    matrixMatrixAdd=matrix(matrix1)
    status         =gsl_matrix_add(matrixMatrixAdd%matrix_,matrix2%matrix_)
    if (status /= GSL_Success) call Error_Report('matrix-matrix summation failed'//{introspection:location})
    return
  end function matrixMatrixAdd
  
  double precision function matrixDeterminant(self)
    !!{
    Compute the determinant of a matrix.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class  (matrix), intent(in   ) :: self
    type   (c_ptr )                :: permutation
    integer(c_int )                :: status     , decompositionSign
    type   (matrix)                :: LU

    if (.not.self%isSquare) call Error_Report('LU decomposition can only be performed on square matrices'//{introspection:location})
    LU         =matrix(self)
    permutation=GSL_Permutation_Alloc(self%size_  (1)                              )
    status     =GSL_LinAlg_LU_Decomp (LU  %matrix_   ,permutation,decompositionSign)
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    matrixDeterminant=GSL_LinAlg_LU_Det(LU%matrix_,decompositionSign)
    call gsl_permutation_free(permutation)
    return
  end function matrixDeterminant

  double precision function matrixLogarithmicDeterminant(self)
    !!{
    Compute the logarithmic determinant of a matrix.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class  (matrix), intent(in   ) :: self
    type   (c_ptr )                :: permutation
    integer(c_int )                :: status     , decompositionSign
    type   (matrix)                :: LU

    if (.not.self%isSquare) call Error_Report('LU decomposition can only be performed on square matrices'//{introspection:location})
    LU         =matrix(self)
    permutation=GSL_Permutation_Alloc(self%size_  (1)                              )
    status     =GSL_LinAlg_LU_Decomp (LU  %matrix_   ,permutation,decompositionSign)
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    matrixLogarithmicDeterminant=GSL_LinAlg_LU_lnDet(LU%matrix_)
    call gsl_permutation_free(permutation)
    return
  end function matrixLogarithmicDeterminant

  integer function matrixSignDeterminant(self)
    !!{
    Compute the sign of the determinant of a matrix.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    class  (matrix), intent(in   ) :: self
    type   (c_ptr )                :: permutation
    integer(c_int )                :: status     , decompositionSign
    type   (matrix)                :: LU

    if (.not.self%isSquare) call Error_Report('LU decomposition can only be performed on square matrices'//{introspection:location})
    LU         =matrix(self)
    permutation=GSL_Permutation_Alloc(self%size_  (1)                              )
    status     =GSL_LinAlg_LU_Decomp (LU  %matrix_   ,permutation,decompositionSign)
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    matrixSignDeterminant=GSL_LinAlg_LU_sgnDet(LU%matrix_,decompositionSign)
    call gsl_permutation_free(permutation)
    return
  end function matrixSignDeterminant

  function matrixTranspose(self)
    !!{
    Transpose a matrix.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (matrix)                :: matrixTranspose
    class  (matrix), intent(in   ) :: self
    integer(c_int )                :: status

    matrixTranspose=matrix(self%size_(2),self%size_(1))
    status=gsl_matrix_transpose_memcpy(matrixTranspose%matrix_,self%matrix_)
    if (status /= GSL_Success) call Error_Report('matrix transpose failed'//{introspection:location})
    return
  end function matrixTranspose

  function matrixInverse(self)
    !!{
    Compute the inverse of a matrix.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (matrix)                :: matrixInverse
    class  (matrix), intent(in   ) :: self
    type   (c_ptr )                :: permutation
    integer(c_int )                :: status       , decompositionSign
    type   (matrix)                :: LU

    if (.not.self%isSquare) call Error_Report('LU decomposition can only be performed on square matrices'//{introspection:location})
    LU           =matrix(self                       )
    matrixInverse=matrix(self%size_(1),self%size_(2))
    permutation  =GSL_Permutation_Alloc(self%size_  (1)                                      )
    status       =GSL_LinAlg_LU_Decomp (LU  %matrix_   ,permutation,decompositionSign        )
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    status       =GSL_LinAlg_LU_Invert (LU  %matrix_   ,permutation,matrixInverse    %matrix_)
    if (status /= GSL_Success) call Error_Report('LU invert failed'       //{introspection:location})
    call gsl_permutation_free(permutation)
    return
  end function matrixInverse

  !! Scalar-matrix functions.
  
  function matrixScalarMultiply(matrix_,scalar_)
    !!{
    Multiply a matrix by a scalar, returning a scalar.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type            (matrix)                :: matrixScalarMultiply
    class           (matrix), intent(in   ) :: matrix_
    double precision        , intent(in   ) :: scalar_
    integer         (c_int )                :: status

    matrixScalarMultiply=matrix(matrix_)
    status              =gsl_matrix_scale(matrixScalarMultiply%matrix_,scalar_)
    if (status /= GSL_Success) call Error_Report('matrix-scalar multiply failed'//{introspection:location})
    return
  end function matrixScalarMultiply

  !! Vector-matrix functions.
  
  function matrixVectorMultiply(matrix_,vector_)
    !!{
    Multiply a matrix by a vector, returning a vector.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (vector)                :: matrixVectorMultiply
    class  (matrix), intent(in   ) :: matrix_
    class  (vector), intent(in   ) :: vector_
    type   (vector)                :: vectorX
    integer(c_int )                :: status

    vectorX             =vector(vector_      )
    matrixVectorMultiply=vector(vector_%size_)
    status              =gsl_blas_dgemv(CblasNoTrans,1.0d0,matrix_%matrix_,vectorX%vector_,0.0d0,matrixVectorMultiply%vector_)
    if (status /= GSL_Success) call Error_Report('matrix-vector multiply failed'//{introspection:location})
    return
  end function matrixVectorMultiply

  function matrixLinearSystemSolve(self,y)
    !!{
    Solve the linear system $y = A \cdot x$.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (vector)                :: matrixLinearSystemSolve
    class  (matrix), intent(inout) :: self
    type   (vector), intent(in   ) :: y
    integer(c_int )                :: status                 , decompositionSign

    if (.not.associated(self%LUdecomposition)) then
       allocate(self%LUdecomposition)
       allocate(self%LUpermutation  )
       self%LUdecomposition=matrix               (self                                                                )
       self%LUpermutation  =GSL_Permutation_Alloc(self%size_                  (1)                                     )
       status              =GSL_LinAlg_LU_Decomp (self%LUdecomposition%matrix_   ,self%LUpermutation,decompositionSign)
       if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    end if    
    matrixLinearSystemSolve=vector             (                                                y%size_                                  )
    status                 =GSL_LinAlg_LU_Solve(self%LUdecomposition%matrix_,self%LUpermutation,y%vector_,matrixLinearSystemSolve%vector_)
    if (status /= GSL_Success) call Error_Report('LU solve failed'//{introspection:location})
    return
  end function matrixLinearSystemSolve

  double precision function matrixCovarianceProduct(self,y,status)
    !!{
    Compute the quantity $y C^{-1} y^\mathrm{T}$ as appears in likelihood functions utilizing covariance matrices. Instead of
    directly inverting the covariance matrix (which is computationally slow and can be inaccurate), we solve the linear system
    $y = C x$ for $x$, and then evaluate $y x$.
    !!}
    use :: Interface_GSL, only : GSL_Success , GSL_ESing
    use :: Error        , only : Error_Report
    implicit none
    class  (matrix  ), intent(inout)           :: self
    type   (vector  ), intent(in   )           :: y
    integer          , intent(  out), optional :: status
    type   (vector  )                          :: CyT
    integer(c_size_t)                          :: i      , j
    logical                                    :: allZero

    if (present(status)) then
       status=GSL_Success
       ! Check that the matrix has no zero rows/columns.
       if (.not.self%nonZeroRowColumnsChecked) then
          do i=1_c_size_t,self%size_(1)
             allZero=.true.
             do j=1_c_size_t,self%size_(2)
                if (gsl_matrix_get(self%matrix_,i-1_c_size_t,j-1_c_size_t) /= 0.0d0) then
                   allZero=.false.
                   exit
                endif
             end do
             if (allZero) self%hasZeroRowColumns=.true.
          end do
          do j=1_c_size_t,self%size_(2)
             allZero=.true.
             do i=1_c_size_t,self%size_(1)
                if (gsl_matrix_get(self%matrix_,i-1_c_size_t,j-1_c_size_t) /= 0.0d0) then
                   allZero=.false.
                   exit
                endif
             end do
             if (allZero) self%hasZeroRowColumns=.true.
          end do
          self%nonZeroRowColumnsChecked=.true.
       end if
       if (self%hasZeroRowColumns) then
          matrixCovarianceProduct=0.0d0
          status                 =GSL_ESing
          return
       end if
    end if
    CyT=self%linearSystemSolve(y)
    matrixCovarianceProduct=y.dot.CyT
    if (matrixCovarianceProduct < 0.0d0) then
       if (present(status)) then
          status=GSL_ESing
       else
          call Error_Report('matrix is not semi-positive definite'//{introspection:location})
       end if
    end if
    return
  end function matrixCovarianceProduct
  
  subroutine matrixEigensystem(self,eigenVectors,eigenValues)
    !!{
    Find eigenvectors and eigenvalues of a real symmetric matrix.
    !!}
    use :: Interface_GSL, only : GSL_Success
    use :: Error        , only : Error_Report
    implicit none
    class  (matrix), intent(inout) :: self
    type   (matrix), intent(  out) :: eigenVectors
    type   (vector), intent(  out) :: eigenValues
    type   (c_ptr )                :: workspace
    integer(c_int )                :: status
    type   (matrix)                :: matrix_

    matrix_     =matrix               (self                                                                    )
    eigenValues =vector               (self   %size_  (1)                                                      )
    eigenVectors=matrix               (self   %size_  (1),self       %size_  (1)                               )
    workspace   =GSL_Eigen_SymmV_Alloc(self   %size_  (1)                                                      )
    status      =GSL_Eigen_SymmV      (matrix_%matrix_   ,eigenValues%vector_   ,eigenVectors%matrix_,workspace)
    if (status /= GSL_Success) call Error_Report('eigensystem evaluation failed'//{introspection:location})
    call GSL_Eigen_Symmv_Free(workspace)
    return
  end subroutine matrixEigensystem

  subroutine matrixCholeskyDecomposition(self)
    !!{
    Find the Cholesky decomposition of a matrix.
    !!}
    use :: Interface_GSL, only : GSL_Success
    use :: Error        , only : Error_Report
    implicit none
    class  (matrix     ), intent(inout) :: self
    integer(c_int      )                :: status

    status=GSL_LinAlg_Cholesky_Decomp(self%matrix_)
    if (status /= GSL_Success) call Error_Report('Cholesky decomposition failed'//{introspection:location})
    return
  end subroutine matrixCholeskyDecomposition

  !! LU matrix functions.

  !! Matrix functions.
  
  function matrixLUConstructor(matrix_) result(self)
    !!{
    Constructor for {\normalfont \ttfamily matrixLU} class which builds the matrix from an array.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (matrixLU)                :: self
    type   (matrix  ), intent(in   ) :: matrix_
    integer(c_int   )                :: status

    if (.not.matrix_%isSquare) call Error_Report('can not find LU decomposition of a non-square matrix'//{introspection:location})
    self%matrix     =matrix(matrix_)
    self%permutation=GSL_Permutation_Alloc(self%size_  (1)                                        )
    status          =GSL_LinAlg_LU_Decomp (self%matrix_   ,self%permutation,self%decompositionSign)
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    return
  end function matrixLUConstructor

  subroutine matrixLUDestructorRank0(self)
    !!{
    Destructor for the {\normalfont \ttfamily matrixLU} class
    !!}
    implicit none
    type(matrixLU), intent(inout) :: self
    
    if (allocated(self%matrix_)) then
       call gsl_matrix_free(self%matrix_)
       deallocate(self%matrix_)
    end if
    if (allocated(self%permutation)) then
       call gsl_permutation_free(self%permutation)
       deallocate(self%permutation)
    end if
    return
  end subroutine matrixLUDestructorRank0

  subroutine matrixLUDestructorRank1(self)
    !!{
    Destructor for the {\normalfont \ttfamily matrixLU} class
    !!}
    implicit none
    type   (matrixLU), intent(inout), dimension(:) :: self
    integer                                      :: i

    do i=1,size(self)
       if (allocated(self(i)%matrix_)) then
          call gsl_matrix_free(self(i)%matrix_)
          deallocate(self(i)%matrix_)
       end if
       if (allocated(self(i)%permutation)) then
          call gsl_permutation_free(self(i)%permutation)
          deallocate(self(i)%permutation)
       end if
    end do
    return
  end subroutine matrixLUDestructorRank1

  function matrixLUSquareSystemSolve(self,y)
    !!{
    Solve the square linear system $y = A \cdot x$.
    !!}
    use :: Error        , only : Error_Report
    use :: Interface_GSL, only : GSL_Success
    implicit none
    type   (vector  )                :: matrixLUSquareSystemSolve
    class  (matrixLU), intent(inout) :: self
    type   (vector  ), intent(in   ) :: y
    integer(c_int   )                :: status

    matrixLUSquareSystemSolve=vector(y%size_)
    status=GSL_LinAlg_LU_Solve(self%matrix_,self%permutation,y%vector_,matrixLUSquareSystemSolve%vector_)
    if (status /= GSL_Success) call Error_Report('LU solve failed'//{introspection:location})
    return
  end function matrixLUSquareSystemSolve

  subroutine matrixLUUnassignment(array,self)
    !!{
    Assign elements of a {\normalfont \ttfamily matrixLU} class to an array
    !!}
    implicit none
    double precision          , intent(  out), dimension(:,:) :: array
    type            (matrixLU), intent(in   )                 :: self
    integer         (c_size_t)                                :: i    , j

    do i=1_c_size_t,self%size_(1)
       do j=1_c_size_t,self%size_(2)
          array(i,j)=gsl_matrix_get(self%matrix_,i-1_c_size_t,j-1_c_size_t)
       end do
    end do
    return
  end subroutine matrixLUUnassignment
  
  !! Geometrical transformations.

  function matrixRotationRandom(randomNumberGenerator_)
    !!{
    Generate a random 3-D rotation matrix. ``Random'' here means that the distribution is invariant when composed with an
    arbitrary rotation (see \href{https://en.wikipedia.org/wiki/Rotation_matrix\#Uniform_random_rotation_matrices}{here} for
    further details).
    !!}
    use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (matrix                    )                 :: matrixRotationRandom
    class           (randomNumberGeneratorClass), intent(inout)  :: randomNumberGenerator_
    double precision                            , dimension(3,3) :: matrixComponents
    double precision                                             :: theta                 , phi   , &
         &                                                          psi                   , x     , &
         &                                                          y                     , z     , &
         &                                                          cosPsi                , sinPsi, &
         &                                                          mcosPsi
    
    theta  =+acos(+2.0d0   *randomNumberGenerator_%uniformSample()-1.0d0)
    phi    =      +2.0d0*Pi*randomNumberGenerator_%uniformSample()
    psi    =      +2.0d0*Pi*randomNumberGenerator_%uniformSample()
    x      =+sin(theta)*cos(phi)
    y      =+sin(theta)*sin(phi)
    z      =+cos(theta)
    cosPsi =+cos(psi)
    sinPsi =+sin(psi)
    mcosPsi=+1.0d0  &
         &  -cosPsi
    matrixComponents=reshape(                                                                    &
         &                   [                                                                   &
         &                    +x*x*mcosPsi+  cosPsi,+x*y*mcosPsi-z*sinPsi,+x*z*mcosPsi+y*sinPsi, &
         &                    +y*x*mcosPsi+z*sinPsi,+y*y*mcosPsi+  cosPsi,+y*z*mcosPsi-x*sinPsi, &
         &                    +z*x*mcosPsi-y*sinPsi,+z*y*mcosPsi+x*sinPsi,+z*z*mcosPsi+  cosPsi  &
         &                   ]                                                                 , &
         &                   [3,3]                                                               &
         &                  )
    matrixRotationRandom=matrix(matrixComponents)
    return
  end function matrixRotationRandom
  
  function matrixRotation(points,pointsRotated)
    !!{
    Given a set of 3 points, and a corresponding set of points to which some rotation has been applied, construct the
    corresponding rotation matrix. The distances between the points must be the same---currently this is not checked. The
    method used is that of \cite{andrei2016}
    !!}
    implicit none
    type            (matrix)                                :: matrixRotation
    type            (vector), intent(in   ), dimension(3  ) :: points           , pointsRotated
    type            (matrix)                                :: P                , Q            , &
         &                                                     Pinverse
    double precision                       , dimension(3,3) :: matrixComponents

    matrixComponents(:,1)=points          (1)
    matrixComponents(:,2)=points          (2)
    matrixComponents(:,3)=points          (3)
    P                    =matrix(matrixComponents)
    matrixComponents(:,1)=pointsRotated   (1)
    matrixComponents(:,2)=pointsRotated   (2)
    matrixComponents(:,3)=pointsRotated   (3)
    Q                    =matrix(matrixComponents)
    Pinverse             =P%inverse()
    matrixRotation       =Q*Pinverse
    return
  end function matrixRotation

  function matrixRotationPlusTranslation(points,pointsRotated,translation) result(matrixRotation)
    !!{
    Given a set of 3 points, and a corresponding set of points to which some rotation has been applied, construct the
    corresponding rotation matrix (and, optionally, any translation between the points). The distances between the points must
    be the same---currently this is not checked. The method used is that of \cite[][their ``More information, easier
    computation'' solution]{robjohn2012}
    !!}
    implicit none
    type            (matrix)                                :: matrixRotation
    type            (vector), intent(in   ), dimension(3  ) :: points           , pointsRotated
    type            (vector), intent(  out), optional       :: translation
    type            (vector)                                :: point4           , pointRotated4 , &
         &                                                     point21          , point31       , &
         &                                                     pointRotated21   , pointRotated31
    type            (matrix)                                :: P                , Q
    double precision                       , dimension(3,3) :: matrixComponents

    point21       =points       (2)-points       (1)
    pointRotated21=pointsRotated(2)-pointsRotated(1)
    point31       =points       (3)-points       (1)
    pointRotated31=pointsRotated(3)-pointsRotated(1)
    point4        =points       (1)+point21       .cross.point31
    pointRotated4 =pointsRotated(1)+pointRotated21.cross.pointRotated31
    matrixComponents(:,1)=point21
    matrixComponents(:,2)=point31
    matrixComponents(:,3)=point4       -points       (1)
    P                    =matrix(matrixComponents)
    matrixComponents(:,1)=pointRotated21
    matrixComponents(:,2)=pointRotated31
    matrixComponents(:,3)=pointRotated4-pointsRotated(1)
    Q                    =matrix(matrixComponents)
    matrixRotation       =Q*P%inverse()
    if (present(translation)) &
         & translation   = pointsRotated (1) &
         &                -matrixRotation    &
         &                *points        (1)
    return
  end function matrixRotationPlusTranslation

end module Linear_Algebra
