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

!!{
Contains a submodule which provides implementations of functions for rank-2, dimension-3, symmetric tensors.
!!}

submodule (Tensors) Tensor_R2_D3_Sym
  !!{
  Provides implementations of functions for rank-2, dimension-3, symmetric tensors.
  !!}
  implicit none

contains

  module procedure tensorRank2Dimension3SymmetricNull
    !!{
    Constructor for {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects which sets all components to zero.
    !!}

    self%x00=0.0d0
    self%x01=0.0d0
    self%x02=0.0d0
    self%x11=0.0d0
    self%x12=0.0d0
    self%x22=0.0d0
    return
  end procedure tensorRank2Dimension3SymmetricNull
  
  function tensorRank2Dimension3SymmetricInternal(x00,x01,x02,x11,x12,x22) result(self)
    !!{
    Constructor for {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
    !!}
    implicit none
    type            (tensorRank2Dimension3Symmetric)                :: self
    double precision                                , intent(in   ) :: x00 , x01, &
         &                                                             x02 , x11, &
         &                                                             x12 , x22
    !![
    <constructorAssign variables="x00, x01, x02, x11, x12, x22"/>
    !!]
    
    return
  end function tensorRank2Dimension3SymmetricInternal
  
  module procedure Tensor_R2_D3_Sym_Destroy
    !!{
    Destroy a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} symmetric object.
    !!}
    implicit none
    !$GLC attributes unused :: self

    return
  end procedure Tensor_R2_D3_Sym_Destroy

  module procedure Tensor_R2_D3_Sym_Builder
    !!{
    Build a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object from the given XML {\normalfont \ttfamily tensorDefinition}.
    !!}
    use :: FoX_DOM, only : node                        , extractDataContent
    use :: Error  , only : Error_Report
    use :: IO_XML , only : XML_Get_Elements_By_Tag_Name, xmlNodeList
    implicit none
    type     (node       )               , pointer     :: element
    type     (xmlNodeList), dimension(:) , allocatable :: elementList
    character(len=3      ), dimension(6) , parameter   :: elementNames=['x00','x01','x02','x11','x12','x22']
    integer                                            :: i

    ! Get the elements
    do i=1,6
       !$omp critical (FoX_DOM_Access)
       call XML_Get_Elements_By_Tag_Name(tensorDefinition,elementNames(i),elementList)
       !$omp end critical (FoX_DOM_Access)
       if (size(elementList) > 1) call Error_Report('multiple "'//elementNames(i)//'" values specified'//{introspection:location})
       if (size(elementList) < 1) call Error_Report('no "'//elementNames(i)//'" value specified'       //{introspection:location})
       !$omp critical (FoX_DOM_Access)
       element => elementList(0)%element
       !$omp end critical (FoX_DOM_Access)
       select case (elementNames(i))
       case ( 'x00' )
          !$omp critical (FoX_DOM_Access)
          call extractDataContent(element,self%x00)
          !$omp end critical (FoX_DOM_Access)
       case ( 'x01' )
          !$omp critical (FoX_DOM_Access)
          call extractDataContent(element,self%x01)
          !$omp end critical (FoX_DOM_Access)
       case ( 'x02' )
          !$omp critical (FoX_DOM_Access)
          call extractDataContent(element,self%x02)
          !$omp end critical (FoX_DOM_Access)
       case ( 'x11' )
          !$omp critical (FoX_DOM_Access)
          call extractDataContent(element,self%x11)
          !$omp end critical (FoX_DOM_Access)
       case ( 'x12' )
          !$omp critical (FoX_DOM_Access)
          call extractDataContent(element,self%x12)
          !$omp end critical (FoX_DOM_Access)
       case ( 'x22' )
          !$omp critical (FoX_DOM_Access)
          call extractDataContent(element,self%x22)
          !$omp end critical (FoX_DOM_Access)
       end select
    end do
    return
  end procedure Tensor_R2_D3_Sym_Builder

  module procedure Tensor_R2_D3_Sym_Dump
    !!{
    Dump properties of a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} symmetric object.
    !!}
    use :: Display           , only : displayMessage
    use :: ISO_Varying_String, only : assignment(=) , varying_string
    implicit none
    character(len=22        ) :: label
    type     (varying_string):: message

    write (label,'(e22.16)') self%x00
    message='x00: '//label
    call displayMessage(message,verbosityLevel)
    write (label,'(e22.16)') self%x01
    message='x01: '//label
    call displayMessage(message,verbosityLevel)
    write (label,'(e22.16)') self%x02
    message='x02: '//label
    call displayMessage(message,verbosityLevel)
    write (label,'(e22.16)') self%x11
    message='x11: '//label
    call displayMessage(message,verbosityLevel)
    write (label,'(e22.16)') self%x12
    message='x12: '//label
    call displayMessage(message,verbosityLevel)
    write (label,'(e22.16)') self%x22
    message='x22: '//label
    call displayMessage(message,verbosityLevel)
    return
  end procedure Tensor_R2_D3_Sym_Dump

  module procedure Tensor_R2_D3_Sym_Dump_Raw
    !!{
    Dump a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object to binary.
    !!}
    implicit none

    ! Dump the content.
    write (fileHandle) self%x00
    write (fileHandle) self%x01
    write (fileHandle) self%x02
    write (fileHandle) self%x11
    write (fileHandle) self%x12
    write (fileHandle) self%x22
    return
  end procedure Tensor_R2_D3_Sym_Dump_Raw

  module procedure Tensor_R2_D3_Sym_Read_Raw
    !!{
    Read a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object from binary.
    !!}
    implicit none

    ! Read the content.
    read (fileHandle) self%x00
    read (fileHandle) self%x01
    read (fileHandle) self%x02
    read (fileHandle) self%x11
    read (fileHandle) self%x12
    read (fileHandle) self%x22
    return
  end procedure Tensor_R2_D3_Sym_Read_Raw

  module procedure Tensor_R2_D3_Sym_Reset
    !!{
    Reset a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
    !!}
    implicit none

    ! Zero all elements.
    self                       %x00=0.0d0
    self                       %x11=0.0d0
    self                       %x22=0.0d0
    self                       %x01=0.0d0
    self                       %x02=0.0d0
    self                       %x12=0.0d0
    return
  end procedure Tensor_R2_D3_Sym_Reset

  module procedure Tensor_R2_D3_Sym_Set_To_Unity
    !!{
    Set a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object to unity.
    !!}
    implicit none

    ! Set values to unity.
    self                       %x00=1.0d0
    self                       %x11=1.0d0
    self                       %x22=1.0d0
    self                       %x01=1.0d0
    self                       %x02=1.0d0
    self                       %x12=1.0d0
    return
  end procedure Tensor_R2_D3_Sym_Set_To_Unity

  module procedure Tensor_R2_D3_Sym_Set_To_Identity
    !!{
    Set a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object to the identity matrix.
    !!}
    implicit none

    ! Set values to unity.
    self                       %x00=1.0d0
    self                       %x11=1.0d0
    self                       %x22=1.0d0
    self                       %x01=0.0d0
    self                       %x02=0.0d0
    self                       %x12=0.0d0
    return
  end procedure Tensor_R2_D3_Sym_Set_To_Identity

  module procedure Tensor_R2_D3_Sym_Is_Zero
    !!{
    Test whether a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object is zero.
    !!}
    implicit none

    ! Detect if all tensor elements are zero.
    Tensor_R2_D3_Sym_Is_Zero=          &
         &         (                   &
         &           self%x00 == 0.0d0 &
         &          .and.              &
         &           self%x01 == 0.0d0 &
         &          .and.              &
         &           self%x02 == 0.0d0 &
         &          .and.              &
         &           self%x11 == 0.0d0 &
         &          .and.              &
         &           self%x12 == 0.0d0 &
         &          .and.              &
         &           self%x22 == 0.0d0 &
         &         )
    return
  end procedure Tensor_R2_D3_Sym_Is_Zero

  module procedure Tensor_R2_D3_Sym_Element
    !!{
    Return the enumeration element of a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
    !!}
    use :: Error, only : Error_Report
    implicit none

    select case (i)
    case (0)
       select case (j)
       case (0)
          Tensor_R2_D3_Sym_Element=self%x00
          return
       case (1)
          Tensor_R2_D3_Sym_Element=self%x01
          return
       case (2)
          Tensor_R2_D3_Sym_Element=self%x02
          return
       end select
    case (1)
       select case (j)
       case (0)
          Tensor_R2_D3_Sym_Element=self%x01
          return
       case (1)
          Tensor_R2_D3_Sym_Element=self%x11
          return
       case (2)
          Tensor_R2_D3_Sym_Element=self%x12
          return
       end select
    case (2)
       select case (j)
       case (0)
          Tensor_R2_D3_Sym_Element=self%x02
          return
       case (1)
          Tensor_R2_D3_Sym_Element=self%x12
          return
       case (2)
          Tensor_R2_D3_Sym_Element=self%x22
          return
       end select
    end select
    Tensor_R2_D3_Sym_Element=0.0d0
    call Error_Report('invalid indices'//{introspection:location})
    return
  end procedure Tensor_R2_D3_Sym_Element

  module procedure Tensor_R2_D3_Sym_Add
    !!{
    Add two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Add%x00=tensor1%x00
    Tensor_R2_D3_Sym_Add%x01=tensor1%x01
    Tensor_R2_D3_Sym_Add%x02=tensor1%x02
    Tensor_R2_D3_Sym_Add%x11=tensor1%x11
    Tensor_R2_D3_Sym_Add%x12=tensor1%x12
    Tensor_R2_D3_Sym_Add%x22=tensor1%x22
    if (present(tensor2)) then
       Tensor_R2_D3_Sym_Add%x00=Tensor_R2_D3_Sym_Add%x00+tensor2%x00
       Tensor_R2_D3_Sym_Add%x01=Tensor_R2_D3_Sym_Add%x01+tensor2%x01
       Tensor_R2_D3_Sym_Add%x02=Tensor_R2_D3_Sym_Add%x02+tensor2%x02
       Tensor_R2_D3_Sym_Add%x11=Tensor_R2_D3_Sym_Add%x11+tensor2%x11
       Tensor_R2_D3_Sym_Add%x12=Tensor_R2_D3_Sym_Add%x12+tensor2%x12
       Tensor_R2_D3_Sym_Add%x22=Tensor_R2_D3_Sym_Add%x22+tensor2%x22
    end if
    return
  end procedure Tensor_R2_D3_Sym_Add

  module procedure Tensor_R2_D3_Sym_Increment
    !!{
    Increment a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
    !!}
    implicit none

    self%x00=self%x00+increment%x00
    self%x01=self%x01+increment%x01
    self%x02=self%x02+increment%x02
    self%x11=self%x11+increment%x11
    self%x12=self%x12+increment%x12
    self%x22=self%x22+increment%x22
    return
  end procedure Tensor_R2_D3_Sym_Increment

  module procedure Tensor_R2_D3_Sym_Subtract
    !!{
    Subtract two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
    !!}
    implicit none

    if (present(tensor2)) then
       Tensor_R2_D3_Sym_Subtract%x00=+tensor1%x00-tensor2%x00
       Tensor_R2_D3_Sym_Subtract%x01=+tensor1%x01-tensor2%x01
       Tensor_R2_D3_Sym_Subtract%x02=+tensor1%x02-tensor2%x02
       Tensor_R2_D3_Sym_Subtract%x11=+tensor1%x11-tensor2%x11
       Tensor_R2_D3_Sym_Subtract%x12=+tensor1%x12-tensor2%x12
       Tensor_R2_D3_Sym_Subtract%x22=+tensor1%x22-tensor2%x22
    else
       Tensor_R2_D3_Sym_Subtract%x00=-tensor1%x00
       Tensor_R2_D3_Sym_Subtract%x01=-tensor1%x01
       Tensor_R2_D3_Sym_Subtract%x02=-tensor1%x02
       Tensor_R2_D3_Sym_Subtract%x11=-tensor1%x11
       Tensor_R2_D3_Sym_Subtract%x12=-tensor1%x12
       Tensor_R2_D3_Sym_Subtract%x22=-tensor1%x22
    end if
    return
  end procedure Tensor_R2_D3_Sym_Subtract

  module procedure Tensor_R2_D3_Sym_Scalar_Multiply
    !!{
    Multiply a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object by a scalar.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Scalar_Multiply%x00=tensor1%x00*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x01=tensor1%x01*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x02=tensor1%x02*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x11=tensor1%x11*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x12=tensor1%x12*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x22=tensor1%x22*multiplier
    return
  end procedure Tensor_R2_D3_Sym_Scalar_Multiply

  module procedure Tensor_R2_D3_Sym_Scalar_Multiply_Switched
    !!{
    Multiply a scalar by a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Scalar_Multiply_Switched=Tensor_R2_D3_Sym_Scalar_Multiply(tensor1,multiplier)
    return
  end procedure Tensor_R2_D3_Sym_Scalar_Multiply_Switched

  module procedure Tensor_R2_D3_Sym_Max
    !!{
    Return an element-by-element {\normalfont \ttfamily max()} on two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Max%x00=max(tensor1%x00,tensor2%x00)
    Tensor_R2_D3_Sym_Max%x01=max(tensor1%x01,tensor2%x01)
    Tensor_R2_D3_Sym_Max%x02=max(tensor1%x02,tensor2%x02)
    Tensor_R2_D3_Sym_Max%x11=max(tensor1%x11,tensor2%x11)
    Tensor_R2_D3_Sym_Max%x12=max(tensor1%x12,tensor2%x12)
    Tensor_R2_D3_Sym_Max%x22=max(tensor1%x22,tensor2%x22)
    return
  end procedure Tensor_R2_D3_Sym_Max

  module procedure Tensor_R2_D3_Sym_Scalar_Divide
    !!{
    Multiply a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object by a scalar.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Scalar_Divide%x00=tensor1%x00/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x01=tensor1%x01/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x02=tensor1%x02/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x11=tensor1%x11/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x12=tensor1%x12/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x22=tensor1%x22/divisor
    return
  end procedure Tensor_R2_D3_Sym_Scalar_Divide
  
  module procedure Tensor_R2_D3_Sym_Vector_Project
    !!{
     Find the magnitude of the projection of a {\normalfont \ttfamily tensorRank2Dimension3Symmetric}/vector dot product onto
    the same vector, $\mathbf{x} \cdot \mathbf{A} \cdot \mathbf{x}$.
    !!}

    Tensor_R2_D3_Sym_Vector_Project=           &
         & +      self%x00*vector(1)*vector(1) &
         & +      self%x11*vector(2)*vector(2) &
         & +      self%x22*vector(3)*vector(3) &
         & +2.0d0*self%x01*vector(1)*vector(2) &
         & +2.0d0*self%x02*vector(1)*vector(3) &
         & +2.0d0*self%x12*vector(2)*vector(3)
    return
  end procedure Tensor_R2_D3_Sym_Vector_Project

  module procedure Tensor_R2_D3_Sym_Double_Contract
    !!{
    Find the double contraction of two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects, $\mathbf{A}:\mathbf{B}$.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Double_Contract=  &
         & +      self%x00*tensor1%x00 &
         & +      self%x11*tensor1%x11 &
         & +      self%x22*tensor1%x22 &
         & +2.0d0*self%x01*tensor1%x01 &
         & +2.0d0*self%x02*tensor1%x02 &
         & +2.0d0*self%x12*tensor1%x12
    return
  end procedure Tensor_R2_D3_Sym_Double_Contract

 module procedure Tensor_R2_D3_Sym_Contract
    !!{
    Return the contraction (trace) of a {\normalfont \ttfamily tensorRank2Dimension3Symmetric}.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Contract=self%x00+self%x11+self%x22
    return
  end procedure Tensor_R2_D3_Sym_Contract

  module procedure Tensor_R2_D3_Sym_Property_Count
    !!{
    Return the number of properties required to track a rank 2, 3 dimensional, symmetric tensor. This is equal to 6.
    !!}
    implicit none
    integer, parameter :: propertyCount=6

    Tensor_R2_D3_Sym_Property_Count=propertyCount
    return
  end procedure Tensor_R2_D3_Sym_Property_Count

  module procedure Tensor_R2_D3_Sym_Deserialize
    !!{
    Pack an array into a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} symmetric structure.
    !!}
    implicit none

    self%x00=tensorArray(1)
    self%x01=tensorArray(2)
    self%x02=tensorArray(3)
    self%x11=tensorArray(4)
    self%x12=tensorArray(5)
    self%x22=tensorArray(6)
    return
  end procedure Tensor_R2_D3_Sym_Deserialize

  module procedure Tensor_R2_D3_Sym_Serialize
    !!{
    Pack a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} into an array.
    !!}
    implicit none

    ! Place tensor into array.
    tensorArray(1)=self%x00
    tensorArray(2)=self%x01
    tensorArray(3)=self%x02
    tensorArray(4)=self%x11
    tensorArray(5)=self%x12
    tensorArray(6)=self%x22
    return
  end procedure Tensor_R2_D3_Sym_Serialize

  module procedure Tensor_R2_D3_Sym_From_Matrix
    !!{
    Construct a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object from a matrix.
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
    ! Set the values of the tensor object.
    self%x00=matrix(1,1)
    self%x01=matrix(1,2)
    self%x02=matrix(1,3)
    self%x11=matrix(2,2)
    self%x12=matrix(2,3)
    self%x22=matrix(3,3)
    return
  end procedure Tensor_R2_D3_Sym_From_Matrix

  module procedure Tensor_R2_D3_Sym_To_Matrix
    !!{
    Construct a matrix from a {\normalfont \ttfamily tensorRank2Dimension3Symmetric}.
    !!}
    implicit none

    Tensor_R2_D3_Sym_To_Matrix(1,1)=self%x00
    Tensor_R2_D3_Sym_To_Matrix(1,2)=self%x01
    Tensor_R2_D3_Sym_To_Matrix(1,3)=self%x02
    Tensor_R2_D3_Sym_To_Matrix(2,1)=self%x01
    Tensor_R2_D3_Sym_To_Matrix(2,2)=self%x11
    Tensor_R2_D3_Sym_To_Matrix(2,3)=self%x12
    Tensor_R2_D3_Sym_To_Matrix(3,1)=self%x02
    Tensor_R2_D3_Sym_To_Matrix(3,2)=self%x12
    Tensor_R2_D3_Sym_To_Matrix(3,3)=self%x22
    return
  end procedure Tensor_R2_D3_Sym_To_Matrix

  module procedure Tensor_R2_D3_Sym_Assign_To
    !!{
    Assign a matrix to a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
    !!}
    implicit none

    call tensor%fromMatrix(matrix)
    return
  end procedure Tensor_R2_D3_Sym_Assign_To

  module procedure Tensor_R2_D3_Sym_Assign_From
    !!{
    Assign a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} to a matrix.
    !!}
    implicit none

    matrix(1,1)=tensor%x00
    matrix(1,2)=tensor%x01
    matrix(1,3)=tensor%x02
    matrix(2,1)=tensor%x01
    matrix(2,2)=tensor%x11
    matrix(2,3)=tensor%x12
    matrix(3,1)=tensor%x02
    matrix(3,2)=tensor%x12
    matrix(3,3)=tensor%x22
    return
  end procedure Tensor_R2_D3_Sym_Assign_From

  module procedure Tensor_R2_D3_Sym_Matrix_Equality
    !!{
    Return true if the supplied tensor and matrix are equal.
    !!}
    implicit none

    Tensor_R2_D3_Sym_Matrix_Equality= &
         &  matrix(1,1) == self%x00   &
         & .and.                      &
         &  matrix(1,2) == self%x01   &
         & .and.                      &
         &  matrix(1,3) == self%x02   &
         & .and.                      &
         &  matrix(2,1) == self%x01   &
         & .and.                      &
         &  matrix(2,2) == self%x11   &
         & .and.                      &
         &  matrix(2,3) == self%x12   &
         & .and.                      &
         &  matrix(3,1) == self%x02   &
         & .and.                      &
         &  matrix(3,2) == self%x12   &
         & .and.                      &
         &  matrix(3,3) == self%x22
    return
  end procedure Tensor_R2_D3_Sym_Matrix_Equality

  module procedure Tensor_R2_D3_Sym_Non_Static_Size_Of
    !!{
    Return the size of any non-static components of the object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    !$GLC attributes unused :: self

    Tensor_R2_D3_Sym_Non_Static_Size_Of=0_c_size_t
    return
  end procedure Tensor_R2_D3_Sym_Non_Static_Size_Of

end submodule Tensor_R2_D3_Sym
