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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!!{
Contains a module which defines the structure used for describing tensors.
!!}

module Tensors
  !!{
  Defines the structure used for describing tensors.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  use            :: FoX_DOM      , only : node
  use            :: Display      , only : enumerationVerbosityLevelType
  implicit none
  private
  public :: tensorRank2Dimension3Symmetric, assignment(=), operator(*), max

  ! Interface to assignment functions.
  interface assignment(=)
     module procedure Tensor_R2_D3_Sym_Assign_To
     module procedure Tensor_R2_D3_Sym_Assign_From
  end interface assignment(=)

  ! Interface to multiplication operators with tensor symmetric objects as their second argument.
  interface operator(*)
     module procedure Tensor_R2_D3_Sym_Scalar_Multiply_Switched
  end interface operator(*)

  ! Interface to max() function for tensor symmetric objects.
  interface max
     module procedure Tensor_R2_D3_Sym_Max
  end interface max

  type :: tensor
     !!{
     A generic tensor type.
     !!}
  end type tensor

  type, extends(tensor) :: tensorRank2Dimension3Symmetric
     !!{
     A rank 2, three dimensional, symmetric tensor.
     !!}
     private
     double precision :: x00,x01,x02,x11,x12,x22
   contains
     !![
     <methods>
       <method description="Return the enumerated element." method="element" />
       <method description="Return true if a tensor object is zero." method="isZero" />
       <method description="Destroy a tensor object." method="destroy" />
       <method description="Multiply a tensor by a scalar." method="operator(*)" />
       <method description="Divide a tensor by a scalar." method="operator(/)" />
       <method description="Add two tensors." method="operator(+)" />
       <method description="Subtract one tensor from another." method="operator(-)" />
       <method description="Return true if a tensor and a matrix are equal." method="operator(==)" />
       <method description="Reset elements to zero." method="reset" />
       <method description="Build tensor object from a provided XML description." method="builder" />
       <method description="Dump the tensor object." method="dump" />
       <method description="Dump the tensor object to binary." method="dumpRaw" />
       <method description="Read the tensor object from binary." method="readRaw" />
       <method description="Set all elements of the tensor object to unity." method="setToUnity" />
       <method description="Set a tensor object to the indentity (i.e. all diagonal elements 1, all other elements 0)." method="setToIdentity" />
       <method description="Return a count of the number of properties in a serialized tensor object." method="serializeCount" />
       <method description="Serialize the tensor object to an array." method="serialize" />
       <method description="Deserialize the tensor object from an array." method="deserialize" />
       <method description="Increment the tensor object." method="increment" />
       <method description="Construct a matrix from a tensor object." method="toMatrix" />
       <method description="Construct a matrix from a tensor object." method="fromMatrix" />
       <method description="Contract a tensor, returning $\mathbf{T}^i_i$." method="contract" />
       <method description="Return the double contraction of two tensors, $\mathbf{A}^i_j \mathbf{B}^j_i$." method="doubleContract" />
       <method description="Return the projection of a tensor/vector dot product onto the same vector, $\sum_{i=1}^3 \sum_{j=1}^3 \mathbf{A}^i_j x_i x_j$." method="vectorProject" />
       <method description="Returns the size of any non-static components of the type." method="nonStaticSizeOf" />
     </methods>
     !!]
     procedure         ::                    Tensor_R2_D3_Sym_Add
     procedure         ::                    Tensor_R2_D3_Sym_Subtract
     procedure         ::                    Tensor_R2_D3_Sym_Scalar_Multiply
     procedure         ::                    Tensor_R2_D3_Sym_Scalar_Divide
     procedure         ::                    Tensor_R2_D3_Sym_Matrix_Equality
     generic           :: operator(+)     => Tensor_R2_D3_Sym_Add
     generic           :: operator(-)     => Tensor_R2_D3_Sym_Subtract
     generic           :: operator(*)     => Tensor_R2_D3_Sym_Scalar_Multiply
     generic           :: operator(/)     => Tensor_R2_D3_Sym_Scalar_Divide
     generic           :: operator(==)    => Tensor_R2_D3_Sym_Matrix_Equality
     procedure         :: nonStaticSizeOf => Tensor_R2_D3_Sym_Non_Static_Size_Of
     procedure         :: element         => Tensor_R2_D3_Sym_Element
     procedure         :: isZero          => Tensor_R2_D3_Sym_Is_Zero
     procedure         :: destroy         => Tensor_R2_D3_Sym_Destroy
     procedure         :: setToIdentity   => Tensor_R2_D3_Sym_Set_To_Identity
     procedure         :: reset           => Tensor_R2_D3_Sym_Reset
     procedure         :: builder         => Tensor_R2_D3_Sym_Builder
     procedure         :: dump            => Tensor_R2_D3_Sym_Dump
     procedure         :: dumpRaw         => Tensor_R2_D3_Sym_Dump_Raw
     procedure         :: readRaw         => Tensor_R2_D3_Sym_Read_Raw
     procedure         :: setToUnity      => Tensor_R2_D3_Sym_Set_To_Unity
     procedure, nopass :: serializeCount  => Tensor_R2_D3_Sym_Property_Count
     procedure         :: serialize       => Tensor_R2_D3_Sym_Serialize
     procedure         :: deserialize     => Tensor_R2_D3_Sym_Deserialize
     procedure         :: increment       => Tensor_R2_D3_Sym_Increment
     procedure         :: fromMatrix      => Tensor_R2_D3_Sym_From_Matrix
     procedure         :: toMatrix        => Tensor_R2_D3_Sym_To_Matrix
     procedure         :: contract        => Tensor_R2_D3_Sym_Contract
     procedure         :: doubleContract  => Tensor_R2_D3_Sym_Double_Contract
     procedure         :: vectorProject   => Tensor_R2_D3_Sym_Vector_Project
  end type tensorRank2Dimension3Symmetric

  ! Constructors for {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
  interface tensorRank2Dimension3Symmetric
     module procedure tensorRank2Dimension3SymmetricNull
     module procedure tensorRank2Dimension3SymmetricInternal
  end interface tensorRank2Dimension3Symmetric
  
  ! Identity, unitary, and null tensors.
  type(tensorRank2Dimension3Symmetric), public :: tensorIdentityR2D3Sym=tensorRank2Dimension3Symmetric(1.0d0,0.0d0,0.0d0,1.0d0,0.0d0,1.0d0)
  type(tensorRank2Dimension3Symmetric), public :: tensorUnitR2D3Sym    =tensorRank2Dimension3Symmetric(1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0)
  type(tensorRank2Dimension3Symmetric), public :: tensorNullR2D3Sym    =tensorRank2Dimension3Symmetric(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0)

  ! Interfaces to type-bound functions.
  interface
     module function tensorRank2Dimension3SymmetricNull() result(self)
       !!{
       Constructor for {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects which sets all components to zero.
       !!}
       type(tensorRank2Dimension3Symmetric) :: self
     end function tensorRank2Dimension3SymmetricNull
     module function tensorRank2Dimension3SymmetricInternal(x00,x01,x02,x11,x12,x22) result(self)
       !!{
       Constructor for {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
       !!}
       type            (tensorRank2Dimension3Symmetric)                :: self
       double precision                                , intent(in   ) :: x00 , x01, &
            &                                                             x02 , x11, &
            &                                                             x12 , x22
     end function tensorRank2Dimension3SymmetricInternal
     module subroutine Tensor_R2_D3_Sym_Destroy(self)
       !!{
       Destroy a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} symmetric object.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(inout) :: self
     end subroutine Tensor_R2_D3_Sym_Destroy
     module subroutine Tensor_R2_D3_Sym_Builder(self,tensorDefinition)
       !!{
       Build a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object from the given XML {\normalfont \ttfamily tensorDefinition}.
       !!}
       implicit none
       class(tensorRank2Dimension3Symmetric), intent(inout)              :: self
       type (node                          ), intent(in   ), pointer     :: tensorDefinition
     end subroutine Tensor_R2_D3_Sym_Builder
     module subroutine Tensor_R2_D3_Sym_Dump(self,verbosityLevel)
       !!{
       Reset a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} symmetric object.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(in   ) :: self
       type (enumerationVerbosityLevelType ), intent(in   ) :: verbosityLevel
     end subroutine Tensor_R2_D3_Sym_Dump
     module subroutine Tensor_R2_D3_Sym_Dump_Raw(self,fileHandle)
       !!{
       Dump a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object to binary.
       !!}
       class  (tensorRank2Dimension3Symmetric), intent(in   ) :: self
       integer                                , intent(in   ) :: fileHandle
     end subroutine Tensor_R2_D3_Sym_Dump_Raw
     module subroutine Tensor_R2_D3_Sym_Read_Raw(self,fileHandle)
       !!{
       Read a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object from binary.
       !!}
       class  (tensorRank2Dimension3Symmetric), intent(inout) :: self
       integer                                , intent(in   ) :: fileHandle
     end subroutine Tensor_R2_D3_Sym_Read_Raw
     module subroutine Tensor_R2_D3_Sym_Reset(self)
       !!{
       Reset a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(inout) :: self
     end subroutine Tensor_R2_D3_Sym_Reset
     module subroutine Tensor_R2_D3_Sym_Set_To_Unity(self)
       !!{
       Set a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object to unity.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(inout) :: self
     end subroutine Tensor_R2_D3_Sym_Set_To_Unity
     module subroutine Tensor_R2_D3_Sym_Set_To_Identity(self)
       !!{
       Set a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object to the identity matrix.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(inout) :: self
     end subroutine Tensor_R2_D3_Sym_Set_To_Identity
     module double precision function Tensor_R2_D3_Sym_Element(self,i,j)
       !!{
       Return the enumeration element of a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
       !!}
       class  (tensorRank2Dimension3Symmetric), intent(in) :: self
       integer                                , intent(in) :: i   , j
     end function Tensor_R2_D3_Sym_Element
     module logical function Tensor_R2_D3_Sym_Is_Zero(self)
       !!{
       Test whether a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object is zero.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(in) :: self
     end function Tensor_R2_D3_Sym_Is_Zero
     module function Tensor_R2_D3_Sym_Add(tensor1,tensor2)
       !!{
       Add two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
       !!}
       type (tensorRank2Dimension3Symmetric)                       :: Tensor_R2_D3_Sym_Add
       class(tensorRank2Dimension3Symmetric), intent(in)           :: tensor1
       class(tensorRank2Dimension3Symmetric), intent(in), optional :: tensor2
     end function Tensor_R2_D3_Sym_Add
     module subroutine Tensor_R2_D3_Sym_Increment(self,increment)
       !!{
       Increment a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(inout) :: self
       class(tensorRank2Dimension3Symmetric), intent(in   ) :: increment
     end subroutine Tensor_R2_D3_Sym_Increment
     module function Tensor_R2_D3_Sym_Subtract(tensor1,tensor2)
       !!{
       Subtract two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
       !!}
       type (tensorRank2Dimension3Symmetric)                       :: Tensor_R2_D3_Sym_Subtract
       class(tensorRank2Dimension3Symmetric), intent(in)           :: tensor1
       class(tensorRank2Dimension3Symmetric), intent(in), optional :: tensor2
     end function Tensor_R2_D3_Sym_Subtract
     module function Tensor_R2_D3_Sym_Scalar_Multiply(tensor1,multiplier)
       !!{
       Multiply a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object by a scalar.
       !!}
       type            (tensorRank2Dimension3Symmetric)             :: Tensor_R2_D3_Sym_Scalar_Multiply
       class           (tensorRank2Dimension3Symmetric), intent(in) :: tensor1
       double precision                                , intent(in) :: multiplier
     end function Tensor_R2_D3_Sym_Scalar_Multiply
     module function Tensor_R2_D3_Sym_Scalar_Multiply_Switched(multiplier,tensor1)
       !!{
       Multiply a scalar by a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
       !!}
       type            (tensorRank2Dimension3Symmetric)                :: Tensor_R2_D3_Sym_Scalar_Multiply_Switched
       type            (tensorRank2Dimension3Symmetric), intent(in   ) :: tensor1
       double precision                                , intent(in   ) :: multiplier
     end function Tensor_R2_D3_Sym_Scalar_Multiply_Switched
     module function Tensor_R2_D3_Sym_Max(tensor1,tensor2)
       !!{
       Return an element-by-element {\normalfont \ttfamily max()} on two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
       !!}
       type(tensorRank2Dimension3Symmetric)                :: Tensor_R2_D3_Sym_Max
       type(tensorRank2Dimension3Symmetric), intent(in   ) :: tensor1,tensor2
     end function Tensor_R2_D3_Sym_Max
     module function Tensor_R2_D3_Sym_Scalar_Divide(tensor1,divisor)
       !!{
       Multiply a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object by a scalar.
       !!}
       type            (tensorRank2Dimension3Symmetric)                :: Tensor_R2_D3_Sym_Scalar_Divide
       class           (tensorRank2Dimension3Symmetric), intent(in   ) :: tensor1
       double precision                                , intent(in   ) :: divisor
     end function Tensor_R2_D3_Sym_Scalar_Divide
     module double precision function Tensor_R2_D3_Sym_Vector_Project(self,vector)
       !!{
       Find the magnitude of the projection of a {\normalfont \ttfamily tensorRank2Dimension3Symmetric}/vector dot product onto
       the same vector, $\mathbf{x} \cdot \mathbf{A} \cdot \mathbf{x}$.
       !!}
       class           (tensorRank2Dimension3Symmetric), intent(in   )               :: self
       double precision                                , intent(in   ), dimension(3) :: vector
     end function Tensor_R2_D3_Sym_Vector_Project
     module double precision function Tensor_R2_D3_Sym_Double_Contract(self,tensor1)
       !!{
       Find the double contraction of two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects, $\mathbf{A}:\mathbf{B}$.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(in   ) :: self,tensor1
     end function Tensor_R2_D3_Sym_Double_Contract
     module double precision function Tensor_R2_D3_Sym_Contract(self)
       !!{
       Return the contraction (trace) of a {\normalfont \ttfamily tensorRank2Dimension3Symmetric}.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(in   ) :: self
     end function Tensor_R2_D3_Sym_Contract
     module integer function Tensor_R2_D3_Sym_Property_Count()
       !!{
       Return the number of properties required to track a rank 2, 3 dimensional, symmetric tensor. This is equal to 6.
       !!}
     end function Tensor_R2_D3_Sym_Property_Count
     module subroutine Tensor_R2_D3_Sym_Deserialize(self,tensorArray)
       !!{
       Pack an array into a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} symmetric structure.
       !!}
       class           (tensorRank2Dimension3Symmetric), intent(inout)               :: self
       double precision                                , intent(in   ), dimension(6) :: tensorArray
     end subroutine Tensor_R2_D3_Sym_Deserialize
     module subroutine Tensor_R2_D3_Sym_Serialize(self,tensorArray)
       !!{
       Pack a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} into an array.
       !!}
       double precision                                , intent(  out), dimension(:) :: tensorArray(6)
       class           (tensorRank2Dimension3Symmetric), intent(in   )               :: self
     end subroutine Tensor_R2_D3_Sym_Serialize
     module subroutine Tensor_R2_D3_Sym_From_Matrix(self,matrix)
       !!{
       Construct a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object from a matrix.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(inout)                 :: self
       double precision                     , intent(in   ), dimension(3,3) :: matrix
     end subroutine Tensor_R2_D3_Sym_From_Matrix
     module function Tensor_R2_D3_Sym_To_Matrix(self)
       !!{
       Construct a matrix from a {\normalfont \ttfamily tensorRank2Dimension3Symmetric}.
       !!}
       double precision                     , dimension(3,3) :: Tensor_R2_D3_Sym_To_Matrix
       class(tensorRank2Dimension3Symmetric), intent(in   )  :: self
     end function Tensor_R2_D3_Sym_To_Matrix
     module subroutine Tensor_R2_D3_Sym_Assign_To(tensor,matrix)
       !!{
       Assign a matrix to a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(inout)                 :: tensor
       double precision                     , intent(in   ), dimension(3,3) :: matrix
     end subroutine Tensor_R2_D3_Sym_Assign_To
     module subroutine Tensor_R2_D3_Sym_Assign_From(matrix,tensor)
       !!{
       Assign a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} to a matrix.
       !!}
       double precision                     , intent(  out), dimension(3,3) :: matrix
       class(tensorRank2Dimension3Symmetric), intent(in   )                 :: tensor
     end subroutine Tensor_R2_D3_Sym_Assign_From
     module logical function Tensor_R2_D3_Sym_Matrix_Equality(self,matrix)
       !!{
       Return true if the supplied tensor and matrix are equal.
       !!}
       class(tensorRank2Dimension3Symmetric), intent(in   )                 :: self
       double precision                     , intent(in   ), dimension(3,3) :: matrix
     end function Tensor_R2_D3_Sym_Matrix_Equality
     module function Tensor_R2_D3_Sym_Non_Static_Size_Of(self)
       !!{
       Return the size of any non-static components of the object.
       !!}
       integer(c_size_t                      )                :: Tensor_R2_D3_Sym_Non_Static_Size_Of
       class  (tensorRank2Dimension3Symmetric), intent(in   ) :: self
     end function Tensor_R2_D3_Sym_Non_Static_Size_Of
  end interface
 
end module Tensors
