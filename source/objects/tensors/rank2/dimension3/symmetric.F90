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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!!{RST
Contains a submodule which provides implementations of functions for rank-2, dimension-3, symmetric tensors.
!!}

submodule (Tensors) Tensor_R2_D3_Sym
  !!{RST
  Provides implementations of functions for rank-2, dimension-3, symmetric tensors.

  The six independent elements of a symmetric 3×3 tensor are stored, packed, in the
  ``c(6)`` component as the upper triangle in row-major order:

  .. math::

     c = \left( T_{00}, T_{01}, T_{02}, T_{11}, T_{12}, T_{22} \right).

  Almost every operation on a symmetric tensor is therefore just an array operation on
  ``c``, with three small mappings capturing the packing: ``packedIndex`` maps a pair of
  (0-based) tensor indices to a storage index, ``multiplicity`` gives how many times each
  packed element appears in the full 3×3 tensor (twice for off-diagonal elements), and
  ``diagonalIndices`` lists the storage indices of the diagonal elements.
  !!}
  implicit none

  ! Mapping from a pair of (0-based) tensor indices (i,j) to the packed storage index in c.
  integer         , dimension(0:2,0:2), parameter :: packedIndex     =reshape([1,2,3,2,4,5,3,5,6],[3,3])
  ! Names of the packed elements (used for XML input and for dumps).
  character(len=3), dimension(6      ), parameter :: elementNames    =['x00','x01','x02','x11','x12','x22']
  ! Multiplicity of each packed element in the full 3×3 tensor (off-diagonal elements
  ! appear twice); used to weight contractions and projections.
  double precision, dimension(6      ), parameter :: multiplicity    =[1.0d0,2.0d0,2.0d0,1.0d0,2.0d0,1.0d0]
  ! Storage indices of the three diagonal elements; used for the trace.
  integer         , dimension(3      ), parameter :: diagonalIndices =[1,4,6]

contains

  module procedure tensorRank2Dimension3SymmetricNull
    !!{RST
    Constructor for ``tensorRank2Dimension3Symmetric`` objects which sets all components to zero.
    !!}

    self%c=0.0d0
    return
  end procedure tensorRank2Dimension3SymmetricNull

  function tensorRank2Dimension3SymmetricInternal(x00,x01,x02,x11,x12,x22) result(self)
    !!{RST
    Constructor for ``tensorRank2Dimension3Symmetric`` objects.
    !!}
    implicit none
    type            (tensorRank2Dimension3Symmetric)                :: self
    double precision                                , intent(in   ) :: x00 , x01, &
         &                                                             x02 , x11, &
         &                                                             x12 , x22

    self%c=[x00,x01,x02,x11,x12,x22]
    return
  end function tensorRank2Dimension3SymmetricInternal

  module procedure Tensor_R2_D3_Sym_Destroy
    !!{RST
    Destroy a ``tensorRank2Dimension3Symmetric`` symmetric object.
    !!}
    implicit none
    !$GLC attributes unused :: self

    return
  end procedure Tensor_R2_D3_Sym_Destroy

  module procedure Tensor_R2_D3_Sym_Builder
    !!{RST
    Build a ``tensorRank2Dimension3Symmetric`` object from the given XML ``tensorDefinition``.
    !!}
    use :: FoX_DOM, only : node                        , extractDataContent
    use :: Error  , only : Error_Report
    use :: IO_XML , only : XML_Get_Elements_By_Tag_Name, xmlNodeList
    implicit none
    type   (node       )               , pointer     :: element
    type   (xmlNodeList), dimension(:) , allocatable :: elementList
    integer                                          :: i

    ! Get the elements.
    do i=1,6
       !$omp critical (FoX_DOM_Access)
       call XML_Get_Elements_By_Tag_Name(tensorDefinition,elementNames(i),elementList)
       !$omp end critical (FoX_DOM_Access)
       if (size(elementList) > 1) call Error_Report('multiple "'//elementNames(i)//'" values specified'//{introspection:location})
       if (size(elementList) < 1) call Error_Report('no "'     //elementNames(i)//'" value specified'  //{introspection:location})
       !$omp critical (FoX_DOM_Access)
       element => elementList(0)%element
       call extractDataContent(element,self%c(i))
       !$omp end critical (FoX_DOM_Access)
    end do
    return
  end procedure Tensor_R2_D3_Sym_Builder

  module procedure Tensor_R2_D3_Sym_Dump
    !!{RST
    Dump properties of a ``tensorRank2Dimension3Symmetric`` symmetric object.
    !!}
    use :: Display           , only : displayMessage
    use :: ISO_Varying_String, only : assignment(=) , varying_string
    implicit none
    character(len=22        ) :: label
    type     (varying_string) :: message
    integer                   :: i

    do i=1,6
       write (label,'(e22.16)') self%c(i)
       message=elementNames(i)//': '//label
       call displayMessage(message,verbosityLevel)
    end do
    return
  end procedure Tensor_R2_D3_Sym_Dump

  module procedure Tensor_R2_D3_Sym_Dump_Raw
    !!{RST
    Dump a ``tensorRank2Dimension3Symmetric`` object to binary.
    !!}
    implicit none

    ! Dump the content.
    write (fileHandle) self%c
    return
  end procedure Tensor_R2_D3_Sym_Dump_Raw

  module procedure Tensor_R2_D3_Sym_Read_Raw
    !!{RST
    Read a ``tensorRank2Dimension3Symmetric`` object from binary.
    !!}
    implicit none

    ! Read the content.
    read (fileHandle) self%c
    return
  end procedure Tensor_R2_D3_Sym_Read_Raw

  module procedure Tensor_R2_D3_Sym_Reset
    !!{RST
    Reset a ``tensorRank2Dimension3Symmetric`` object.
    !!}
    implicit none

    ! Zero all elements.
    self%c=0.0d0
    return
  end procedure Tensor_R2_D3_Sym_Reset

  module procedure Tensor_R2_D3_Sym_Set_To_Unity
    !!{RST
    Set a ``tensorRank2Dimension3Symmetric`` object to unity.
    !!}
    implicit none

    ! Set values to unity.
    self%c=1.0d0
    return
  end procedure Tensor_R2_D3_Sym_Set_To_Unity

  module procedure Tensor_R2_D3_Sym_Set_To_Identity
    !!{RST
    Set a ``tensorRank2Dimension3Symmetric`` object to the identity matrix.
    !!}
    implicit none

    ! Set diagonal elements to unity and off-diagonal elements to zero.
    self%c                 =0.0d0
    self%c(diagonalIndices)=1.0d0
    return
  end procedure Tensor_R2_D3_Sym_Set_To_Identity

  module procedure Tensor_R2_D3_Sym_Is_Zero
    !!{RST
    Test whether a ``tensorRank2Dimension3Symmetric`` object is zero.
    !!}
    implicit none

    ! Detect if all tensor elements are zero.
    Tensor_R2_D3_Sym_Is_Zero=all(self%c == 0.0d0)
    return
  end procedure Tensor_R2_D3_Sym_Is_Zero

  module procedure Tensor_R2_D3_Sym_Element
    !!{RST
    Return the enumeration element of a ``tensorRank2Dimension3Symmetric`` object.
    !!}
    use :: Error, only : Error_Report
    implicit none

    if (i < 0 .or. i > 2 .or. j < 0 .or. j > 2) then
       Tensor_R2_D3_Sym_Element=0.0d0
       call Error_Report('invalid indices'//{introspection:location})
    else
       Tensor_R2_D3_Sym_Element=self%c(packedIndex(i,j))
    end if
    return
  end procedure Tensor_R2_D3_Sym_Element

  module procedure Tensor_R2_D3_Sym_Add
    !!{RST
    Add two ``tensorRank2Dimension3Symmetric`` objects.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Add%c=tensor1%c
    if (present(tensor2)) Tensor_R2_D3_Sym_Add%c=Tensor_R2_D3_Sym_Add%c+tensor2%c
    return
  end procedure Tensor_R2_D3_Sym_Add

  module procedure Tensor_R2_D3_Sym_Increment
    !!{RST
    Increment a ``tensorRank2Dimension3Symmetric`` object.
    !!}
    implicit none

    self%c=self%c+increment%c
    return
  end procedure Tensor_R2_D3_Sym_Increment

  module procedure Tensor_R2_D3_Sym_Subtract
    !!{RST
    Subtract two ``tensorRank2Dimension3Symmetric`` objects.
    !!}
    implicit none

    if (present(tensor2)) then
       Tensor_R2_D3_Sym_Subtract%c=+tensor1%c-tensor2%c
    else
       Tensor_R2_D3_Sym_Subtract%c=-tensor1%c
    end if
    return
  end procedure Tensor_R2_D3_Sym_Subtract

  module procedure Tensor_R2_D3_Sym_Scalar_Multiply
    !!{RST
    Multiply a ``tensorRank2Dimension3Symmetric`` object by a scalar.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Scalar_Multiply%c=tensor1%c*multiplier
    return
  end procedure Tensor_R2_D3_Sym_Scalar_Multiply

  module procedure Tensor_R2_D3_Sym_Scalar_Multiply_Switched
    !!{RST
    Multiply a scalar by a ``tensorRank2Dimension3Symmetric`` object.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Scalar_Multiply_Switched=Tensor_R2_D3_Sym_Scalar_Multiply(tensor1,multiplier)
    return
  end procedure Tensor_R2_D3_Sym_Scalar_Multiply_Switched

  module procedure Tensor_R2_D3_Sym_Max
    !!{RST
    Return an element-by-element ``max()`` on two ``tensorRank2Dimension3Symmetric`` objects.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Max%c=max(tensor1%c,tensor2%c)
    return
  end procedure Tensor_R2_D3_Sym_Max

  module procedure Tensor_R2_D3_Sym_Scalar_Divide
    !!{RST
    Multiply a ``tensorRank2Dimension3Symmetric`` object by a scalar.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Scalar_Divide%c=tensor1%c/divisor
    return
  end procedure Tensor_R2_D3_Sym_Scalar_Divide

  module procedure Tensor_R2_D3_Sym_Vector_Project
    !!{RST
    Find the magnitude of the projection of a ``tensorRank2Dimension3Symmetric``/vector dot product onto the same vector, :math:`\mathbf{x} \cdot \mathbf{A} \cdot \mathbf{x}`.
    !!}

    Tensor_R2_D3_Sym_Vector_Project=            &
         & +      self%c(1)*vector(1)*vector(1) &
         & +      self%c(4)*vector(2)*vector(2) &
         & +      self%c(6)*vector(3)*vector(3) &
         & +2.0d0*self%c(2)*vector(1)*vector(2) &
         & +2.0d0*self%c(3)*vector(1)*vector(3) &
         & +2.0d0*self%c(5)*vector(2)*vector(3)
    return
  end procedure Tensor_R2_D3_Sym_Vector_Project

  module procedure Tensor_R2_D3_Sym_Double_Contract
    !!{RST
    Find the double contraction of two ``tensorRank2Dimension3Symmetric`` objects, :math:`\mathbf{A}:\mathbf{B}`.
    !!}
    implicit none

    ! Off-diagonal elements are counted twice (they appear both above and below the diagonal).
    Tensor_R2_D3_Sym_Double_Contract=sum(multiplicity*self%c*tensor1%c)
    return
  end procedure Tensor_R2_D3_Sym_Double_Contract

  module procedure Tensor_R2_D3_Sym_Contract
    !!{RST
    Return the contraction (trace) of a ``tensorRank2Dimension3Symmetric``.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Contract=sum(self%c(diagonalIndices))
    return
  end procedure Tensor_R2_D3_Sym_Contract

  module procedure Tensor_R2_D3_Sym_Property_Count
    !!{RST
    Return the number of properties required to track a rank 2, 3 dimensional, symmetric tensor. This is equal to 6.
    !!}
    implicit none
    integer, parameter :: propertyCount=6

    Tensor_R2_D3_Sym_Property_Count=propertyCount
    return
  end procedure Tensor_R2_D3_Sym_Property_Count

  module procedure Tensor_R2_D3_Sym_Deserialize
    !!{RST
    Pack an array into a ``tensorRank2Dimension3Symmetric`` symmetric structure.
    !!}
    implicit none

    self%c=tensorArray
    return
  end procedure Tensor_R2_D3_Sym_Deserialize

  module procedure Tensor_R2_D3_Sym_Serialize
    !!{RST
    Pack a ``tensorRank2Dimension3Symmetric`` into an array.
    !!}
    implicit none

    ! Place tensor into array.
    tensorArray=self%c
    return
  end procedure Tensor_R2_D3_Sym_Serialize

  module procedure Tensor_R2_D3_Sym_From_Matrix
    !!{RST
    Construct a ``tensorRank2Dimension3Symmetric`` object from a matrix.
    !!}
    use :: Error, only : Error_Report
    implicit none

    ! Test for symmetry
    if     (                            &
         &   matrix(1,2) /= matrix(2,1) &
         &  .or.                        &
         &   matrix(1,3) /= matrix(3,1) &
         &  .or.                        &
         &   matrix(2,3) /= matrix(3,2) &
         & ) call Error_Report('supplied matrix is not symmetric'//{introspection:location})
    ! Set the values of the tensor object from the upper triangle of the matrix.
    self%c=[matrix(1,1),matrix(1,2),matrix(1,3),matrix(2,2),matrix(2,3),matrix(3,3)]
    return
  end procedure Tensor_R2_D3_Sym_From_Matrix

  module procedure Tensor_R2_D3_Sym_To_Matrix
    !!{RST
    Construct a matrix from a ``tensorRank2Dimension3Symmetric``.
    !!}
    implicit none
    integer :: i, j

    do i=1,3
       do j=1,3
          Tensor_R2_D3_Sym_To_Matrix(i,j)=self%c(packedIndex(i-1,j-1))
       end do
    end do
    return
  end procedure Tensor_R2_D3_Sym_To_Matrix

  module procedure Tensor_R2_D3_Sym_Assign_To
    !!{RST
    Assign a matrix to a ``tensorRank2Dimension3Symmetric`` object.
    !!}
    implicit none

    call tensor%fromMatrix(matrix)
    return
  end procedure Tensor_R2_D3_Sym_Assign_To

  module procedure Tensor_R2_D3_Sym_Assign_From
    !!{RST
    Assign a ``tensorRank2Dimension3Symmetric`` to a matrix.
    !!}
    implicit none

    matrix=tensor%toMatrix()
    return
  end procedure Tensor_R2_D3_Sym_Assign_From

  module procedure Tensor_R2_D3_Sym_Matrix_Equality
    !!{RST
    Return true if the supplied tensor and matrix are equal.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Matrix_Equality=all(matrix == self%toMatrix())
    return
  end procedure Tensor_R2_D3_Sym_Matrix_Equality

  module procedure Tensor_R2_D3_Sym_Non_Static_Size_Of
    !!{RST
    Return the size of any non-static components of the object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    !$GLC attributes unused :: self

    Tensor_R2_D3_Sym_Non_Static_Size_Of=0_c_size_t
    return
  end procedure Tensor_R2_D3_Sym_Non_Static_Size_Of

end submodule Tensor_R2_D3_Sym
