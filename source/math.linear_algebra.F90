!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use FGSL
  implicit none
  private
  public :: assignment(=), operator(*)

  type, public :: vector
     !% Vector class.
     private
     double precision, dimension(:), allocatable :: elements
   contains
     !# <workaround type="gfortran" PR="58471 58470" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58471 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58470">
     !# final     :: vectorDestroy
     !# </workaround>
     !@ <objectMethods>
     !@   <object>vector</object>
     !@   <objectMethod>
     !@     <method>subtract</method>
     !@     <type>\textcolor{red}{\textless type(vector)</type>
     !@     <arguments>\textcolor{red}{\textless class(vector)\textgreater} vector1\argin, \textcolor{red}{\textless class(vector)\textgreater} vector2\argin</arguments>
     !@     <description>Compute {\tt vector1}-{\tt vector2}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>add</method>
     !@     <type>\textcolor{red}{\textless type(vector)</type>
     !@     <arguments>\textcolor{red}{\textless class(vector)\textgreater} vector1\argin, \textcolor{red}{\textless class(vector)\textgreater} vector2\argin</arguments>
     !@     <description>Compute {\tt vector1}+{\tt vector2}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: subtract    => vectorSubtract
     procedure :: add         => vectorAdd
     generic   :: operator(-) => subtract
     generic   :: operator(+) => add
  end type vector
  
  type, public :: matrix
     !% Matrix class.
     private
     double precision, dimension(:,:), allocatable :: elements
   contains
     !# <workaround type="gfortran" PR="58471 58470" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58471 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58470">
     !# final     :: matrixDestroy
     !# </workaround>
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
     !@     <method>linearSystemSolve</method>
     !@     <type>\textcolor{red}{\textless type(vector)</type>
     !@     <arguments>\textcolor{red}{\textless type(vector) y\argin \textgreater}</arguments>
     !@     <description>Solve the linear system $y = A \cdot x$ where $A$ is ourself and $y$ is the specified vector.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: invert                 => matrixInvert
     procedure :: logarithmicDeterminant => matrixLogarithmicDeterminant
     procedure :: linearSystemSolve      => matrixLinearSystemSolve
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
  end interface operator(*)

contains

  subroutine arrayToVectorAssign(self,array)
    !% Assign an array to a vector.
    implicit none
    type            (vector), intent(  out)               :: self
    double precision        , intent(in   ), dimension(:) :: array

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

  function vectorSubtract(vector1,vector2)
    !% Subtract one vector from another.
    implicit none   
    type (vector)                :: vectorSubtract
    class(vector), intent(in   ) :: vector1        , vector2

    vectorSubtract%elements=vector1%elements-vector2%elements
    return
  end function vectorSubtract

  function vectorAdd(vector1,vector2)
    !% Add one vector to another.
    implicit none   
    type (vector)                :: vectorAdd
    class(vector), intent(in   ) :: vector1  , vector2

    vectorAdd%elements=vector1%elements+vector2%elements
    return
  end function vectorAdd

  subroutine arrayToMatrixAssign(self,array)
    !% Assign an array to a matrix.
    implicit none
    type            (matrix), intent(  out)                 :: self
    double precision        , intent(in   ), dimension(:,:) :: array

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
  
  function matrixInvert(self)
    !% Invert a matrix.
    implicit none
    type            (matrix          )                                                                         :: matrixInvert
    class           (matrix          ), intent(inout)                                                          :: self
    type            (fgsl_matrix     )                                                                         :: selfMatrix       , selfInverse
    type            (fgsl_permutation)                                                                         :: permutations
    integer         (kind=fgsl_int   )                                                                         :: decompositionSign, status
    integer         (kind=fgsl_size_t)                                                                         :: selfMatrixSize
    double precision                  , dimension(size(self%elements,dim=1),size(self%elements,dim=2)), target :: inverse
    double precision                  , dimension(size(self%elements,dim=1),size(self%elements,dim=2))         :: selfArray

    selfMatrixSize       =size(self%elements,dim=1)
    selfMatrix           =FGSL_Matrix_Init(type=1.0_fgsl_double)
    selfInverse          =FGSL_Matrix_Init(type=1.0_fgsl_double)
    permutations         =FGSL_Permutation_Alloc(selfMatrixSize)
    selfArray            =self%elements
    status               =FGSL_Matrix_Align(self%elements,selfMatrixSize,selfMatrixSize,selfMatrixSize,selfMatrix )
    status               =FGSL_Matrix_Align(inverse      ,selfMatrixSize,selfMatrixSize,selfMatrixSize,selfInverse)
    status               =FGSL_LinAlg_LU_Decomp(selfMatrix,permutations,decompositionSign)
    status               =FGSL_LinAlg_LU_Invert(selfMatrix,permutations,selfInverse      )
    matrixInvert%elements=inverse
    call FGSL_Permutation_Free(permutations)
    call FGSL_Matrix_Free     (selfMatrix)
    call FGSL_Matrix_Free     (selfInverse)
    ! Restore the original matrix.
    self%elements        =selfArray
    return
  end function matrixInvert

  double precision function matrixLogarithmicDeterminant(self)
    !% Return the logarithm of the determinant of a matrix.
    implicit none
    class           (matrix          ), intent(inout)                                                  :: self
    type            (fgsl_matrix     )                                                                 :: selfMatrix
    type            (fgsl_permutation)                                                                 :: permutations
    integer         (kind=fgsl_int   )                                                                 :: decompositionSign, status
    integer         (kind=fgsl_size_t)                                                                 :: selfMatrixSize
    double precision                  , dimension(size(self%elements,dim=1),size(self%elements,dim=2)) :: selfArray

    selfMatrixSize              =size(self%elements,dim=1)
    selfMatrix                  =FGSL_Matrix_Init(type=1.0_fgsl_double)
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

  function matrixLinearSystemSolve(self,y)
    !% Solve the linear system $y = A \cdot x$.
    implicit none
    type            (vector          )                                                                 :: matrixLinearSystemSolve
    class           (matrix          ), intent(inout)                                                  :: self
    type            (vector          ), intent(in   )                                                  :: y
    type            (fgsl_matrix     )                                                                 :: selfMatrix
    type            (fgsl_vector     )                                                                 :: xVector          , yVector
    type            (fgsl_permutation)                                                                 :: permutations
    integer         (kind=fgsl_int   )                                                                 :: decompositionSign, status
    integer         (kind=fgsl_size_t)                                                                 :: selfMatrixSize
    double precision                  , dimension(size(self%elements,dim=1),size(self%elements,dim=2)) :: selfArray

    selfMatrixSize              =size(self%elements,dim=1)
    selfMatrix                  =FGSL_Matrix_Init(type=1.0_fgsl_double)
    xVector                     =FGSL_Vector_Init(type=1.0_fgsl_double)
    yVector                     =FGSL_Vector_Init(type=1.0_fgsl_double)
    permutations                =FGSL_Permutation_Alloc(selfMatrixSize)
    selfArray                   =self%elements
    matrixLinearSystemSolve     =y
    status                      =FGSL_Vector_Align(matrixLinearSystemSolve%elements,selfMatrixSize,xVector,selfMatrixSize,0_FGSL_Size_T,1_FGSL_Size_T)
    status                      =FGSL_Vector_Align(                      y%elements,selfMatrixSize,yVector,selfMatrixSize,0_FGSL_Size_T,1_FGSL_Size_T)
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
  
end module Linear_Algebra
