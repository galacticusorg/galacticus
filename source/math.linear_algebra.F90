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
Contains a module which implements linear algebra calculations.
!!}

module Linear_Algebra
  !!{
  Implements linear algebra calculations.
  !!}
  use, intrinsic :: ISO_C_Binding   , only : c_ptr          , c_double, c_size_t, c_int, &
       &                                     c_null_ptr
  use            :: Resource_Manager, only : resourceManager
  implicit none
  private
  public :: vector        , matrix         , matrixRotation      , matrixRotationPlusTranslation, &
       &    matrixLU      , assignment(=)  , operator(*)         , gsl_vector_get               , &
       &    gsl_vector_set, gsl_vector_free, matrixRotationRandom

  type :: vectorWrapper
     !!{
     Wrapper class for managing GSL vectors.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: vectorWrapperDestructor
  end type vectorWrapper
  
  type, public :: vector
     !!{
     Vector class.
     !!}
     private
     type   (resourceManager)          :: vectorManager
     type   (vectorWrapper  ), pointer :: vector_       => null()
     integer(c_size_t       )          :: size_
   contains
     !![
     <methods>
       <method description="Compute the magnitude of a vector."                                                method="magnitude"        />
       <method description="Compute {\normalfont \ttfamily vector1} $\cdot$ {\normalfont \ttfamily vector2}."  method="operator(.dot.)"  />
       <method description="Compute {\normalfont \ttfamily vector1}-{\normalfont \ttfamily vector2}."          method="operator(-)"      />
       <method description="Compute {\normalfont \ttfamily vector1}+{\normalfont \ttfamily vector2}."          method="operator(+)"      />
       <method description="Compute {\normalfont \ttfamily vector1} $\times$ {\normalfont \ttfamily vector2}." method="operator(.cross.)"/>
       <method description="Return a C pointer to the GSL vector object."                                      method="gslObject"        />
       <method description="Assign vector objects."                                                            method="assignment(=)"    />
     </methods>
     !!]
     procedure ::                        vectorAssign
     generic   :: assignment(=)       => vectorAssign
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
  
  type :: matrixWrapper
     !!{
     Wrapper class for managing GSL matrices.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: matrixWrapperDestructor
  end type matrixWrapper
  
  type, public :: matrix
     !!{
     Matrix class.
     !!}
     private
     type   (resourceManager)               :: matrixManager
     type   (matrixWrapper  ), pointer      :: matrix_       => null()
     integer(c_size_t       ), dimension(2) :: size_
     logical                                :: isSquare
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
       <method description="Assign matrix objects."                                                                            method="assignment(=)"         />
     </methods>
     !!]
     procedure ::                           matrixAssign
     generic   :: assignment(=)          => matrixAssign
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
  
  type :: permutationWrapper
     !!{
     Wrapper class for managing GSL permutations.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: permutationWrapperDestructor
  end type permutationWrapper
  
  type, public, extends(matrix) :: matrixLU
     !!{
     Matrix class for LU matrices.
     !!}
     private
     type   (resourceManager   )          :: permutationManager
     type   (permutationWrapper), pointer :: permutation_       => null()
     integer(c_int             )          :: decompositionSign
   contains
     !![
     <methods>
       <method description="Solve the linear system $y = A \cdot x$ where $A$ is ourself and $y$ is the specified vector." method="squareSystemSolve" />
     </methods>
     !!]
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
    class           (*       ), pointer                     :: dummyPointer_
    
    allocate(self%vector_)
    self%size_      =size(array,kind=c_size_t)
    self%vector_%gsl=gsl_vector_alloc(self%size_)
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_      => self%vector_
    self%vectorManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    do i=1,size(array,dim=1,kind=c_size_t)
       call gsl_vector_set(self%vector_%gsl,i-1_c_size_t,array(i))
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
    class  (*       ), pointer       :: dummyPointer_

    allocate(self%vector_)
    self%size_      =n
    self%vector_%gsl=gsl_vector_alloc(self%size_)
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_      => self%vector_
    self%vectorManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    call gsl_vector_set_zero(self%vector_%gsl)
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
    class  (*     ), pointer    :: dummyPointer_

    allocate(self%vector_)
    self%size_      =source%size_
    self%vector_%gsl=gsl_vector_alloc(self%size_)
    status          =gsl_vector_memcpy(self%vector_%gsl,source%vector_%gsl)
    if (status /= GSL_Success) call Error_Report('vector copy failed'//{introspection:location})
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_      => self%vector_
    self%vectorManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
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
       array(i)=gsl_vector_get(self%vector_%gsl,i-1_c_size_t)
    end do
    return
  end subroutine vectorUnassignment

  subroutine vectorAssign(self,from)
    !!{
    Perform assignment of vectors.
    !!}
    implicit none
    class(vector), intent(inout) :: self
    class(vector), intent(in   ) :: from

    self%vectorManager =  from%vectorManager
    self%vector_       => from%vector_
    self%size_         =  from%size_
    return
  end subroutine vectorAssign

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
    status   =gsl_vector_add(vectorAdd%vector_%gsl,vector2%vector_%gsl)
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
    status        =gsl_vector_sub(vectorSubtract%vector_%gsl,vector2%vector_%gsl)
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

    status=gsl_blas_ddot(vector1%vector_%gsl,vector2%vector_%gsl,vectorDotProduct)
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
    vectorCrossProduct=vector(                                                                                                 &
         &                    [                                                                                                &
         &                     +gsl_vector_get(vector1%vector_%gsl,1_c_size_t)*gsl_vector_get(vector2%vector_%gsl,2_c_size_t)  &
         &                     -gsl_vector_get(vector1%vector_%gsl,2_c_size_t)*gsl_vector_get(vector2%vector_%gsl,1_c_size_t), &
         &                     +gsl_vector_get(vector1%vector_%gsl,2_c_size_t)*gsl_vector_get(vector2%vector_%gsl,0_c_size_t)  &
         &                     -gsl_vector_get(vector1%vector_%gsl,0_c_size_t)*gsl_vector_get(vector2%vector_%gsl,2_c_size_t), &
         &                     +gsl_vector_get(vector1%vector_%gsl,0_c_size_t)*gsl_vector_get(vector2%vector_%gsl,1_c_size_t)  &
         &                     -gsl_vector_get(vector1%vector_%gsl,1_c_size_t)*gsl_vector_get(vector2%vector_%gsl,0_c_size_t)  &
         &                    ]                                                                                                &
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

    vectorGSLObject=self%vector_%gsl
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
    integer         (c_size_t)                                :: i            , j
    class           (*       ), pointer                       :: dummyPointer_

    allocate(self%matrix_)
    self%size_       =shape(array,kind=c_size_t)
    self%isSquare    =self%size_(1) == self%size_(2)
    self%matrix_ %gsl=gsl_matrix_alloc(self%size_(1),self%size_(2))
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_      => self%matrix_
    self%matrixManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    do i=1,size(array,dim=1,kind=c_size_t)
       do j=1,size(array,dim=2,kind=c_size_t)
          call gsl_matrix_set(self%matrix_%gsl,i-1_c_size_t,j-1_c_size_t,array(i,j))
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
    integer(c_size_t), intent(in   ) :: n1           , n2
    class  (*       ), pointer       :: dummyPointer_

    allocate(self%matrix_)
    self%isSquare    =n1 == n2
    self%size_       =[n1,n2]
    self%matrix_ %gsl=gsl_matrix_alloc(self%size_(1),self%size_(2))
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_      => self%matrix_
    self%matrixManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    call gsl_matrix_set_zero(self%matrix_%gsl)
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
    class  (*     ), pointer    :: dummyPointer_

    allocate(self%matrix_)
    self%size_       =source%size_
    self%isSquare    =source%isSquare
    self%matrix_ %gsl=gsl_matrix_alloc (self%size_      (1),self  %size_      (2))
    status           =gsl_matrix_memcpy(self%matrix_%gsl   ,source%matrix_%gsl   )
    if (status /= GSL_Success) call Error_Report('matrix copy failed'//{introspection:location})
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_      => self%matrix_
    self%matrixManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
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
          array(i,j)=gsl_matrix_get(self%matrix_%gsl,i-1_c_size_t,j-1_c_size_t)
       end do
    end do
    return
  end subroutine matrixUnassignment

  subroutine matrixAssign(self,from)
    !!{
    Perform assignment of matrices.
    !!}
    implicit none
    class(matrix), intent(inout) :: self
    class(matrix), intent(in   ) :: from

    self%matrixManager =  from%matrixManager
    self%matrix_       => from%matrix_
    self%size_         =  from%size_
    self%isSquare      =  from%isSquare
    select type (self)
    class is (matrixLU)
       select type (from)
       class is (matrixLU)
          self%decompositionSign  =  from%decompositionSign
          self%permutationManager =  from%permutationManager
          self%permutation_       => from%permutation_
       end select
    end select
    return
  end subroutine matrixAssign

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
    status             =gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0d0,matrix1%matrix_%gsl,matrix2%matrix_%gsl,0.0d0,matrixMatrixProduct%matrix_%gsl)
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
    status         =gsl_matrix_add(matrixMatrixAdd%matrix_%gsl,matrix2%matrix_%gsl)
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
    permutation=GSL_Permutation_Alloc(self%size_      (1)                              )
    status     =GSL_LinAlg_LU_Decomp (LU  %matrix_%gsl   ,permutation,decompositionSign)
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    matrixDeterminant=GSL_LinAlg_LU_Det(LU%matrix_%gsl,decompositionSign)
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
    permutation=GSL_Permutation_Alloc(self%size_      (1)                              )
    status     =GSL_LinAlg_LU_Decomp (LU  %matrix_%gsl   ,permutation,decompositionSign)
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    matrixLogarithmicDeterminant=GSL_LinAlg_LU_lnDet(LU%matrix_%gsl)
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
    permutation=GSL_Permutation_Alloc(self%size_      (1)                              )
    status     =GSL_LinAlg_LU_Decomp (LU  %matrix_%gsl   ,permutation,decompositionSign)
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    matrixSignDeterminant=GSL_LinAlg_LU_sgnDet(LU%matrix_%gsl,decompositionSign)
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
    status=gsl_matrix_transpose_memcpy(matrixTranspose%matrix_%gsl,self%matrix_%gsl)
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
    permutation  =GSL_Permutation_Alloc(self%size_      (1)                                      )
    status       =GSL_LinAlg_LU_Decomp (LU  %matrix_%gsl   ,permutation,decompositionSign        )
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    status       =GSL_LinAlg_LU_Invert (LU  %matrix_%gsl   ,permutation,matrixInverse    %matrix_%gsl)
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
    status              =gsl_matrix_scale(matrixScalarMultiply%matrix_%gsl,scalar_)
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
    status              =gsl_blas_dgemv(CblasNoTrans,1.0d0,matrix_%matrix_%gsl,vectorX%vector_%gsl,0.0d0,matrixVectorMultiply%vector_%gsl)
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
    type   (c_ptr )                :: permutation
    integer(c_int )                :: status                 , decompositionSign
    type   (matrix)                :: LU

    matrixLinearSystemSolve=vector(y   %size_)
    LU                     =matrix(self      )
    permutation            =GSL_Permutation_Alloc(self%size_      (1)                                                                      )
    status                 =GSL_LinAlg_LU_Decomp (LU  %matrix_%gsl   ,permutation,decompositionSign                                        )
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    status                 =GSL_LinAlg_LU_Solve  (LU%matrix_%gsl     ,permutation,y%vector_%gsl        ,matrixLinearSystemSolve%vector_%gsl)
    if (status /= GSL_Success) call Error_Report('LU solve failed'        //{introspection:location})
    call gsl_permutation_free(permutation)
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
       do i=1_c_size_t,self%size_(1)
          allZero=.true.
          do j=1_c_size_t,self%size_(2)
             if (gsl_matrix_get(self%matrix_%gsl,i-1_c_size_t,j-1_c_size_t) /= 0.0d0) then
                allZero=.false.
                exit
             endif
          end do
          if (allZero) then
             matrixCovarianceProduct=0.0d0
             status                 =GSL_ESing
             return
          end if
       end do
       do j=1_c_size_t,self%size_(2)
          allZero=.true.
          do i=1_c_size_t,self%size_(1)
             if (gsl_matrix_get(self%matrix_%gsl,i-1_c_size_t,j-1_c_size_t) /= 0.0d0) then
                allZero=.false.
                exit
             endif
          end do
          if (allZero) then
             matrixCovarianceProduct=0.0d0
             status                 =GSL_ESing
             return
          end if
       end do
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

    matrix_     =matrix               (self                                                                                )
    eigenValues =vector               (self   %size_      (1)                                                              )
    eigenVectors=matrix               (self   %size_      (1),self       %size_      (1)                                   )
    workspace   =GSL_Eigen_SymmV_Alloc(self   %size_      (1)                                                              )
    status      =GSL_Eigen_SymmV      (matrix_%matrix_%gsl   ,eigenValues%vector_%gsl   ,eigenVectors%matrix_%gsl,workspace)
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

    status=GSL_LinAlg_Cholesky_Decomp(self%matrix_%gsl)
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
    class  (*       ), pointer       :: dummyPointer_
    
    if (.not.matrix_%isSquare) call Error_Report('can not find LU decomposition of a non-square matrix'//{introspection:location})
    allocate(self%permutation_)
    self%matrix          =matrix(matrix_)
    self%permutation_%gsl=GSL_Permutation_Alloc(self%size_      (1)                                             )
    status               =GSL_LinAlg_LU_Decomp (self%matrix_%gsl   ,self%permutation_%gsl,self%decompositionSign)
    if (status /= GSL_Success) call Error_Report('LU decomposition failed'//{introspection:location})
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_           => self%permutation_
    self%permutationManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    return
  end function matrixLUConstructor

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
    status=GSL_LinAlg_LU_Solve(self%matrix_%gsl,self%permutation_%gsl,y%vector_%gsl,matrixLUSquareSystemSolve%vector_%gsl)
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
          array(i,j)=gsl_matrix_get(self%matrix_%gsl,i-1_c_size_t,j-1_c_size_t)
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

  subroutine vectorWrapperDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily vectorWrapper} class.
    !!}
    implicit none
    type(vectorWrapper), intent(inout) :: self

    call gsl_vector_free(self%gsl)
    return
  end subroutine vectorWrapperDestructor
  
  subroutine matrixWrapperDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily matrixWrapper} class.
    !!}
    implicit none
    type(matrixWrapper), intent(inout) :: self

    call gsl_matrix_free(self%gsl)
    return
  end subroutine matrixWrapperDestructor
  
  subroutine permutationWrapperDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily permutationWrapper} class.
    !!}
    implicit none
    type(permutationWrapper), intent(inout) :: self

    call gsl_permutation_free(self%gsl)
    return
  end subroutine permutationWrapperDestructor
  
end module Linear_Algebra
