!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements linear algebra calculations.

module Linear_Algebra
  !% Implements linear algebra calculations.
  use :: FGSL, only : FGSL_Eigen_SymmV          , FGSL_Eigen_SymmV_Alloc, FGSL_Eigen_Symmv_Free, FGSL_LinAlg_Cholesky_Decomp, &
          &           FGSL_LinAlg_LU_Decomp     , FGSL_LinAlg_LU_Det    , FGSL_LinAlg_LU_Invert, FGSL_LinAlg_LU_Solve       , &
          &           FGSL_LinAlg_LU_lnDet      , FGSL_Matrix_Align     , FGSL_Matrix_Free     , FGSL_Permutation_Alloc     , &
          &           FGSL_Permutation_Free     , FGSL_Vector_Align     , FGSL_Vector_Free     , FGSL_Vector_Init           , &
          &           fgsl_eigen_symmv_workspace, fgsl_matrix           , fgsl_matrix_init     , fgsl_permutation           , &
          &           fgsl_vector
  implicit none
  private
  public :: assignment(=), operator(*), matrixRotation, matrixRotationPlusTranslation

  type, public :: vector
     !% Vector class.
     private
     double precision, dimension(:), allocatable :: elements
   contains
     !@ <objectMethods>
     !@   <object>vector</object>
     !@   <objectMethod>
     !@     <method>subtract</method>
     !@     <type>\textcolor{red}{\textless type(vector)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless class(vector)\textgreater} vector1\argin, \textcolor{red}{\textless class(vector)\textgreater} vector2\argin</arguments>
     !@     <description>Compute {\normalfont \ttfamily vector1}-{\normalfont \ttfamily vector2}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>add</method>
     !@     <type>\textcolor{red}{\textless type(vector)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless class(vector)\textgreater} vector1\argin, \textcolor{red}{\textless class(vector)\textgreater} vector2\argin</arguments>
     !@     <description>Compute {\normalfont \ttfamily vector1}+{\normalfont \ttfamily vector2}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>crossProduct</method>
     !@     <type>\textcolor{red}{\textless type(vector)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless class(vector)\textgreater} vector1\argin, \textcolor{red}{\textless class(vector)\textgreater} vector2\argin</arguments>
     !@     <description>Compute {\normalfont \ttfamily vector1} $\times$ {\normalfont \ttfamily vector2}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dotProduct</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless class(vector)\textgreater} vector1\argin, \textcolor{red}{\textless class(vector)\textgreater} vector2\argin</arguments>
     !@     <description>Compute {\normalfont \ttfamily vector1} $\cdot$ {\normalfont \ttfamily vector2}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>magnitude</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Compute the magnitude of a vector.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                      vectorDestroy
     procedure :: subtract          => vectorSubtract
     procedure :: add               => vectorAdd
     procedure :: crossProduct      => vectorCrossProduct
     procedure :: dotProduct        => vectorDotProduct
     procedure :: magnitude         => vectorMagnitude
     generic   :: operator(-)       => subtract
     generic   :: operator(+)       => add
     generic   :: operator(.cross.) => crossProduct
     generic   :: operator(.dot.  ) => dotProduct
  end type vector

  type, public :: matrix
     !% Matrix class.
     private
     double precision, dimension(:,:), allocatable :: elements
   contains
     !@ <objectMethods>
     !@   <object>matrix</object>
     !@   <objectMethod>
     !@     <method>invert</method>
     !@     <type>\textcolor{red}{\textless type(matrix) \textgreater}</type>
     !@     <arguments></arguments>
     !@     <description>Compute and return the matrix inverse.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>logarithmicDeterminant</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Compute and return the logarithm of the determinant of the matrix.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>determinant</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Compute and return the determinant of the matrix.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>linearSystemSolve</method>
     !@     <type>\textcolor{red}{\textless type(vector)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless type(vector) y\argin \textgreater}</arguments>
     !@     <description>Solve the linear system $y = A \cdot x$ where $A$ is ourself and $y$ is the specified vector.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>transpose</method>
     !@     <type>\textcolor{red}{\textless type(matrix)\textgreater}</type>
     !@     <arguments></arguments>
     !@     <description>Return the transpose of a matrix.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>makeSemiPositiveDefinite</method>
     !@     <type>\textcolor{red}{\textless type(matrix)\textgreater}</type>
     !@     <arguments></arguments>
     !@     <description>Make a matrix semi-positive definite by setting any negative eigenvalues to zero.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>eigenSystem</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(matrix)\textgreater} eigenVectors\argout, \textcolor{red}{\textless type(vector)\textgreater} eigenValues\argout</arguments>
     !@     <description>Compute eigenvectors and eigenvalues of the matrix.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>symmetrize</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Make the matrix symmetric using $A_{ij} \rightarrow (A_{ij}+A_{ji})/2$.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>choleskyDecompose</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Compute the Cholesky decomposition of the matrix in place.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>covarianceProduct</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(vector)\textgreater} y\argin</arguments>
     !@     <description>Compute $y C^{-1} y^\mathrm{T}$ as appears in likelihood functions utilizing covariance matrices.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                             matrixDestroy
     procedure :: invert                   => matrixInvert
     procedure :: logarithmicDeterminant   => matrixLogarithmicDeterminant
     procedure :: determinant              => matrixDeterminant
     procedure :: linearSystemSolve        => matrixLinearSystemSolve
     procedure :: transpose                => matrixTranspose
     procedure :: eigenSystem              => matrixEigensystem
     procedure :: symmetrize               => matrixSymmetrize
     procedure :: makeSemiPositiveDefinite => matrixMakeSemiPositiveDefinite
     procedure :: choleskyDecompose        => matrixCholeskyDecompose
     procedure :: covarianceProduct        => matrixCovarianceProduct
  end type matrix

  ! Assignment interfaces.
  interface assignment(=)
     module procedure arrayToVectorAssign
     module procedure vectorToArrayAssign
     module procedure vectorToVectorAssign
     module procedure arrayToMatrixAssign
     module procedure matrixToArrayAssign
  end interface assignment(=)

  ! Operator interfaces.
  interface operator(*)
     module procedure vectorVectorMultiply
     module procedure matrixVectorMultiply
     module procedure matrixMatrixMultiply
  end interface operator(*)

contains

  subroutine arrayToVectorAssign(self,array)
    !% Assign an array to a vector.
    implicit none
    type            (vector), intent(  out)               :: self
    double precision        , intent(in   ), dimension(:) :: array

    if (allocated(self%elements)) deallocate(self%elements)
    self%elements=array
    return
  end subroutine arrayToVectorAssign

  subroutine vectorToArrayAssign(array,vector1)
    !% Assign a vector to an array.
    implicit none
    type            (vector), intent(in   )               :: vector1
    double precision        , intent(  out), dimension(:) :: array

    array=vector1%elements
    return
  end subroutine vectorToArrayAssign

  subroutine vectorToVectorAssign(vector1,vector2)
    !% Assign a vector to an array.
    implicit none
    type            (vector), intent(  out)                 :: vector1
    type            (vector), intent(in   )                 :: vector2

    if (allocated(vector1%elements)) deallocate(vector1%elements)
    vector1%elements=vector2%elements
    return
  end subroutine vectorToVectorAssign

  subroutine vectorDestroy(self)
    !% Destroy a vector object.
    implicit none
    type(vector), intent(inout) :: self

    if (allocated(self%elements)) deallocate(self%elements)
    return
  end subroutine vectorDestroy

  double precision function vectorMagnitude(self)
    !% Return the magnitude of a vector.
    implicit none
    class(vector), intent(in   ) :: self

    vectorMagnitude=sqrt(sum(self%elements**2))
    return
  end function vectorMagnitude

  function vectorSubtract(vector1,vector2)
    !% Subtract one vector from another.
    implicit none
    type (vector)                :: vectorSubtract
    class(vector), intent(in   ) :: vector1       , vector2

    allocate(vectorSubtract%elements(size(vector1%elements,dim=1)))
    vectorSubtract%elements=vector1%elements-vector2%elements
    return
  end function vectorSubtract

  function vectorAdd(vector1,vector2)
    !% Add one vector to another.
    implicit none
    type (vector)                :: vectorAdd
    class(vector), intent(in   ) :: vector1  , vector2

    allocate(vectorAdd%elements(size(vector1%elements,dim=1)))
    vectorAdd%elements=vector1%elements+vector2%elements
    return
  end function vectorAdd

  double precision function vectorDotProduct(vector1,vector2)
    !% Compute the dot product of two vectors.
    implicit none
    class(vector), intent(in   ) :: vector1, vector2

    vectorDotProduct=sum(vector1%elements*vector2%elements)
    return
  end function vectorDotProduct

  function vectorCrossProduct(vector1,vector2)
    !% Compute the cross product of two vectors.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type (vector)                :: vectorCrossProduct
    class(vector), intent(in   ) :: vector1           , vector2

    if (size(vector1%elements) /= 3 .or. size(vector2%elements) /= 3) &
         & call Galacticus_Error_Report('vector cross product only defined for 3D vectors'//{introspection:location})
    allocate(vectorCrossProduct%elements(3))
    vectorCrossProduct%elements=[                                                                                 &
         &                       vector1%elements(2)*vector2%elements(3)-vector1%elements(3)*vector2%elements(2), &
         &                       vector1%elements(3)*vector2%elements(1)-vector1%elements(1)*vector2%elements(3), &
         &                       vector1%elements(1)*vector2%elements(2)-vector1%elements(2)*vector2%elements(1)  &
         &                      ]
    return
  end function vectorCrossProduct

  subroutine arrayToMatrixAssign(self,array)
    !% Assign an array to a matrix.
    implicit none
    type            (matrix), intent(  out)                 :: self
    double precision        , intent(in   ), dimension(:,:) :: array

    if (allocated(self%elements)) deallocate(self%elements)
    self%elements=array
    return
  end subroutine arrayToMatrixAssign

  subroutine matrixToArrayAssign(array,matrix1)
    !% Assign a matrix to an array.
    implicit none
    type            (matrix), intent(in   )                 :: matrix1
    double precision        , intent(  out), dimension(:,:) :: array

    array=matrix1%elements
    return
  end subroutine matrixToArrayAssign

  subroutine matrixDestroy(self)
    !% Destroy a matrix object.
    implicit none
    type(matrix), intent(inout) :: self

    if (allocated(self%elements)) deallocate(self%elements)
    return
  end subroutine matrixDestroy

  subroutine matrixMakeSemiPositiveDefinite(self)
    !% Make a matrix semi-positive definite by setting any negative
    !% eigenvalues to zero and reconstructing the matrix from its
    !% eigenvectors.
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (matrix  ), intent(inout)                 :: self
    double precision          , allocatable  , dimension(:,:) :: eigenValuesMatrixArray
    double precision          , allocatable  , dimension(:  ) :: eigenValueArray
    integer         (c_size_t)                                :: selfMatrixSize        , i
    type            (matrix  )                                :: eigenValuesMatrix     , eigenVectors, eigenVectorsInverse
    type            (vector  )                                :: eigenValues

    selfMatrixSize=size(self%elements,dim=1)
    call self%eigenSystem(eigenVectors,eigenValues)
    allocate(eigenValuesMatrixArray(selfMatrixSize,selfMatrixSize))
    allocate(eigenValueArray(selfMatrixSize))
    eigenValueArray=eigenValues
    eigenValuesMatrixArray=0.0d0
    do i=1,selfMatrixSize
       if (eigenValueArray(i) > 0.0d0) eigenValuesMatrixArray(i,i)=eigenValueArray(i)
    end do
    eigenValuesMatrix  =eigenValuesMatrixArray
    eigenVectors       =eigenVectors%transpose()
    eigenVectorsInverse=eigenVectors%invert()
    select type (self)
    type is (matrix)
       self            =eigenVectors*eigenValuesMatrix*eigenVectorsInverse
    end select
    return
  end subroutine matrixMakeSemiPositiveDefinite

  function matrixInvert(self)
    !% Invert a matrix.
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_int, c_double
    implicit none
    type            (matrix          )                                                                         :: matrixInvert
    class           (matrix          ), intent(inout)                                                          :: self
    type            (fgsl_matrix     )                                                                         :: selfMatrix       , selfInverse
    type            (fgsl_permutation)                                                                         :: permutations
    integer         (c_int           )                                                                         :: decompositionSign, status
    integer         (c_size_t        )                                                                         :: selfMatrixSize
    double precision                  , dimension(size(self%elements,dim=1),size(self%elements,dim=2)), target :: inverse
    double precision                  , dimension(size(self%elements,dim=1),size(self%elements,dim=2))         :: selfArray

    selfMatrixSize=size(self%elements,dim=1)
    selfMatrix    =FGSL_Matrix_Init(type=1.0_c_double)
    selfArray     =self%elements
    status        =FGSL_Matrix_Align(self%elements,selfMatrixSize,selfMatrixSize,selfMatrixSize,selfMatrix )
    selfInverse   =FGSL_Matrix_Init      (type=1.0_c_double)
    status        =FGSL_Matrix_Align     (inverse      ,selfMatrixSize,selfMatrixSize,selfMatrixSize,selfInverse)
    permutations  =FGSL_Permutation_Alloc(selfMatrixSize                               )
    status        =FGSL_LinAlg_LU_Decomp (selfMatrix    ,permutations,decompositionSign)
    status        =FGSL_LinAlg_LU_Invert (selfMatrix    ,permutations,selfInverse      )
    allocate(matrixInvert%elements(size(inverse,dim=1),size(inverse,dim=2)))
    matrixInvert%elements=inverse
    call FGSL_Matrix_Free     (selfMatrix  )
    call FGSL_Matrix_Free     (selfInverse )
    call FGSL_Permutation_Free(permutations)
    ! Restore the original matrix.
    self%elements=selfArray
    return
  end function matrixInvert

  double precision function matrixLogarithmicDeterminant(self)
    !% Return the logarithm of the determinant of a matrix.
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_int, c_double
    implicit none
    class           (matrix          ), intent(inout)                                                  :: self
    type            (fgsl_matrix     )                                                                 :: selfMatrix
    type            (fgsl_permutation)                                                                 :: permutations
    integer         (c_int           )                                                                 :: decompositionSign, status
    integer         (c_size_t        )                                                                 :: selfMatrixSize
    double precision                  , dimension(size(self%elements,dim=1),size(self%elements,dim=2)) :: selfArray

    selfMatrixSize              =size(self%elements,dim=1)
    selfMatrix                  =FGSL_Matrix_Init(type=1.0_c_double)
    permutations                =FGSL_Permutation_Alloc(selfMatrixSize)
    selfArray                   =self%elements
    status                      =FGSL_Matrix_Align(self%elements,selfMatrixSize,selfMatrixSize,selfMatrixSize,selfMatrix )
    status                      =FGSL_LinAlg_LU_Decomp(selfMatrix,permutations,decompositionSign)
    matrixLogarithmicDeterminant=FGSL_LinAlg_LU_lnDet(selfMatrix)
    call FGSL_Permutation_Free(permutations)
    call FGSL_Matrix_Free     (selfMatrix  )
    ! Restore the original matrix.
    self%elements               =selfArray
    return
  end function matrixLogarithmicDeterminant

  double precision function matrixDeterminant(self)
    !% Return the of a matrix.
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_int, c_double
    implicit none
    class           (matrix          ), intent(inout)                                                  :: self
    type            (fgsl_matrix     )                                                                 :: selfMatrix
    type            (fgsl_permutation)                                                                 :: permutations
    integer         (c_int           )                                                                 :: decompositionSign, status
    integer         (c_size_t        )                                                                 :: selfMatrixSize
    double precision                  , dimension(size(self%elements,dim=1),size(self%elements,dim=2)) :: selfArray

    selfMatrixSize   =size(self%elements,dim=1)
    selfMatrix       =FGSL_Matrix_Init(type=1.0_c_double)
    permutations     =FGSL_Permutation_Alloc(selfMatrixSize)
    selfArray        =self%elements
    status           =FGSL_Matrix_Align(self%elements,selfMatrixSize,selfMatrixSize,selfMatrixSize,selfMatrix )
    status           =FGSL_LinAlg_LU_Decomp(selfMatrix,permutations,decompositionSign)
    matrixDeterminant=FGSL_LinAlg_LU_Det(selfMatrix,decompositionSign)
    call FGSL_Permutation_Free(permutations)
    call FGSL_Matrix_Free     (selfMatrix  )
    ! Restore the original matrix.
    self%elements    =selfArray
    return
  end function matrixDeterminant

  function vectorVectorMultiply(vector1,vector2)
    !% Multiply a vector by a vector, returning a scalar.
    implicit none
    double precision                        :: vectorVectorMultiply
    class           (vector), intent(in   ) :: vector1             , vector2

    vectorVectorMultiply=sum(vector1%elements*vector2%elements)
    return
  end function vectorVectorMultiply

  function matrixVectorMultiply(matrix1,vector2)
    !% Multiply a matrix by a vector, returning a vector.
    implicit none
    type   (vector     )                :: matrixVectorMultiply
    class  (matrix     ), intent(in   ) :: matrix1
    class  (vector     ), intent(in   ) :: vector2
    integer                             :: i

    allocate(matrixVectorMultiply%elements(size(vector2%elements)))
    forall(i=1:size(vector2%elements))
       matrixVectorMultiply%elements(i)=sum(matrix1%elements(i,:)*vector2%elements)
    end forall
    return
  end function matrixVectorMultiply

  function matrixMatrixMultiply(matrix1,matrix2)
    !% Multiply a matrix by a matrix, returning a matrix.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type   (matrix)                :: matrixMatrixMultiply
    class  (matrix), intent(in   ) :: matrix1, matrix2

    if (size(matrix1%elements,dim=2) /= size(matrix2%elements,dim=1)) call Galacticus_Error_Report('dimension mismatch'//{introspection:location})
    allocate(matrixMatrixMultiply%elements(size(matrix1%elements,dim=1),size(matrix2%elements,dim=2)))
    matrixMatrixMultiply%elements=matmul(matrix1%elements,matrix2%elements)
    return
  end function matrixMatrixMultiply

  function matrixLinearSystemSolve(self,y)
    !% Solve the linear system $y = A \cdot x$.
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_int, c_double
    implicit none
    type            (vector          )                                                                 :: matrixLinearSystemSolve
    class           (matrix          ), intent(inout)                                                  :: self
    type            (vector          ), intent(in   )                                                  :: y
    type            (fgsl_matrix     )                                                                 :: selfMatrix
    type            (fgsl_vector     )                                                                 :: xVector          , yVector
    type            (fgsl_permutation)                                                                 :: permutations
    integer         (c_int           )                                                                 :: decompositionSign, status
    integer         (c_size_t        )                                                                 :: selfMatrixSize
    double precision                  , dimension(size(self%elements,dim=1),size(self%elements,dim=2)) :: selfArray

    selfMatrixSize              =size(self%elements,dim=1)
    selfMatrix                  =FGSL_Matrix_Init(type=1.0_c_double)
    xVector                     =FGSL_Vector_Init(type=1.0_c_double)
    yVector                     =FGSL_Vector_Init(type=1.0_c_double)
    permutations                =FGSL_Permutation_Alloc(selfMatrixSize)
    selfArray                   =self%elements
    matrixLinearSystemSolve     =y
    status                      =FGSL_Vector_Align(matrixLinearSystemSolve%elements,selfMatrixSize,xVector,selfMatrixSize,0_c_size_t,1_c_size_t)
    status                      =FGSL_Vector_Align(                      y%elements,selfMatrixSize,yVector,selfMatrixSize,0_c_size_t,1_c_size_t)
    status                      =FGSL_Matrix_Align(self%elements,selfMatrixSize,selfMatrixSize,selfMatrixSize,selfMatrix)
    status                      =FGSL_LinAlg_LU_Decomp(selfMatrix,permutations,decompositionSign)
    status                      =FGSL_LinAlg_LU_Solve(selfMatrix,permutations,yVector,xVector)
    call FGSL_Permutation_Free(permutations)
    call FGSL_Matrix_Free     (selfMatrix  )
    call FGSL_Vector_Free     (yVector     )
    call FGSL_Vector_Free     (xVector     )
    ! Restore the original matrix.
    self%elements               =selfArray
    return
  end function matrixLinearSystemSolve

  double precision function matrixCovarianceProduct(self,y,status)
    !% Compute the quantity $y C^{-1} y^\mathrm{T}$ as appears in likelihood functions utilizing covariance matrices. Instead of
    !% directly inverting the covariance matrix (which is computationally slow and can be inaccurate), we solve the linear system
    !% $y = C x$ for $x$, and then evaluate $y x$.
    use :: Interface_GSL, only : GSL_Success, GSL_ESing
    implicit none
    class  (matrix), intent(inout)           :: self
    type   (vector), intent(in   )           :: y
    integer        , intent(  out), optional :: status
    integer                                  :: i

    if (present(status)) then
       status=GSL_Success
       ! Check that the matrix has no zero rows/columns.
       do i=1,size(self%elements,dim=1)
          if (all(self%elements(i,:) == 0.0d0)) then
             matrixCovarianceProduct=0.0d0
             status                 =GSL_ESing
             return
          end if
       end do
       do i=1,size(self%elements,dim=2)
          if (all(self%elements(:,i) == 0.0d0)) then
             matrixCovarianceProduct=0.0d0
             status                 =GSL_ESing
             return
          end if
       end do
    end if
    matrixCovarianceProduct=y*self%linearSystemSolve(y)
    return
  end function matrixCovarianceProduct

  function matrixTranspose(self)
    !% Transpose a matrix.
    implicit none
    type (matrix)                :: matrixTranspose
    class(matrix), intent(inout) :: self

    allocate(matrixTranspose%elements(size(self%elements,dim=2),size(self%elements,dim=1)))
    matrixTranspose%elements=transpose(self%elements)
    return
  end function matrixTranspose

  subroutine matrixEigensystem(self,eigenVectors,eigenValues)
    !% Find eigenvectors and eigenvalues of a real symmetric matrix.
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_int, c_double
    implicit none
    class           (matrix                    ), intent(inout)                                                          :: self
    type            (matrix                    ), intent(  out)                                                          :: eigenVectors
    type            (vector                    ), intent(  out)                                                          :: eigenValues
    type            (fgsl_matrix               )                                                                         :: selfMatrix      , eigenVectorMatrix
    type            (fgsl_vector               )                                                                         :: eigenValueVector
    type            (fgsl_eigen_symmv_workspace)                                                                         :: workspace
    integer         (c_int                     )                                                                         :: status
    integer         (c_size_t                  )                                                                         :: selfMatrixSize
    double precision                            , dimension(size(self%elements,dim=1),size(self%elements,dim=2)), target :: eigenVectorArray
    double precision                            , dimension(size(self%elements,dim=1)                          ), target :: eigenValueArray
    double precision                            , dimension(size(self%elements,dim=1),size(self%elements,dim=2))         :: selfArray

    selfMatrixSize   =size(self%elements,dim=1)
    selfArray        =self%elements
    selfMatrix       =FGSL_Matrix_Init      (type=1.0_c_double)
    eigenVectorMatrix=FGSL_Matrix_Init      (type=1.0_c_double)
    eigenValueVector =FGSL_Vector_Init      (type=1.0_c_double)
    status           =FGSL_Matrix_Align     (self%elements   ,selfMatrixSize,selfMatrixSize  ,selfMatrixSize,selfMatrix                  )
    status           =FGSL_Matrix_Align     (eigenVectorArray,selfMatrixSize,selfMatrixSize  ,selfMatrixSize,eigenVectorMatrix           )
    status           =FGSL_Vector_Align     (eigenValueArray ,selfMatrixSize,eigenValueVector,selfMatrixSize,0_c_size_t       ,1_c_size_t)
    workspace        =FGSL_Eigen_SymmV_Alloc(                 selfMatrixSize                                                             )
    status           =FGSL_Eigen_SymmV      (selfMatrix                     ,eigenValueVector,eigenVectorMatrix,workspace                )
    eigenVectors     =eigenVectorArray
    eigenValues      =eigenValueArray
    call FGSL_Eigen_Symmv_Free(workspace        )
    call FGSL_Matrix_Free     (selfMatrix       )
    call FGSL_Matrix_Free     (eigenVectorMatrix)
    call FGSL_Vector_Free     (eigenValueVector )
    ! Restore the original matrix.
    self%elements               =selfArray
    return
  end subroutine matrixEigensystem

  subroutine matrixSymmetrize(self)
    !% Symmetrize a matrix.
    implicit none
    class(matrix), intent(inout) :: self

    self%elements=0.5d0*(self%elements+transpose(self%elements))
    return
  end subroutine matrixSymmetrize

  subroutine matrixCholeskyDecompose(self)
    !% Find the Cholesky decomposition of a matrix.
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_int, c_double
    implicit none
    class  (matrix     ), intent(inout) :: self
    type   (fgsl_matrix)                :: selfMatrix
    integer(c_int      )                :: status
    integer(c_size_t   )                :: selfMatrixSize

    selfMatrixSize=size(self%elements,dim=1)
    selfMatrix    =FGSL_Matrix_Init           (type=1.0_c_double)
    status        =FGSL_Matrix_Align          (self%elements,selfMatrixSize,selfMatrixSize,selfMatrixSize,selfMatrix)
    status        =FGSL_LinAlg_Cholesky_Decomp(selfMatrix)
    call FGSL_Matrix_Free(selfMatrix)
    return
  end subroutine matrixCholeskyDecompose

  function matrixRotation(points,pointsRotated)
    !% Given a set of 3 points, and a corresponding set of points to which some rotation has been applied, construct the
    !% corresponding rotation matrix. The distances between the points must be the same---currently this is not checked. The
    !% method used is that of \cite{andrei2016}
    implicit none
    type            (matrix)                                :: matrixRotation
    type            (vector), intent(in   ), dimension(3  ) :: points           , pointsRotated
    type            (matrix)                                :: P                , Q
    double precision                       , dimension(3,3) :: matrixComponents

    matrixComponents(:,1)=points          (1)
    matrixComponents(:,2)=points          (2)
    matrixComponents(:,3)=points          (3)
    P                    =matrixComponents
    matrixComponents(:,1)=pointsRotated   (1)
    matrixComponents(:,2)=pointsRotated   (2)
    matrixComponents(:,3)=pointsRotated   (3)
    Q                    =matrixComponents
    matrixRotation       =Q*P%invert()
    return
  end function matrixRotation

  function matrixRotationPlusTranslation(points,pointsRotated,translation) result(matrixRotation)
    !% Given a set of 3 points, and a corresponding set of points to which some rotation has been applied, construct the
    !% corresponding rotation matrix (and, optionally, any translation between the points). The distances between the points must
    !% be the same---currently this is not checked. The method used is that of \cite[][their ``More information, easier
    !% computation'' solution]{robjohn2012}
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
    P                    =matrixComponents
    matrixComponents(:,1)=pointRotated21
    matrixComponents(:,2)=pointRotated31
    matrixComponents(:,3)=pointRotated4-pointsRotated(1)
    Q                    =matrixComponents
    matrixRotation       =Q*P%invert()
    if (present(translation)) &
         & translation   = pointsRotated (1) &
         &                -matrixRotation    &
         &                *points        (1)
    return
  end function matrixRotationPlusTranslation

end module Linear_Algebra
