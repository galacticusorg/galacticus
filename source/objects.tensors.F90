!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which defines the tensor symmetric structure used for describing symmetric tensors.

module Tensors
  !% Defines the tensor symmetric structure used for describing symmetric tensors.
  use ISO_Varying_String
  use Numerical_Constants_Astronomical
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
     !% A generic tensor type.
  end type tensor

  type, extends(tensor) :: tensorRank2Dimension3Symmetric
     !% A rank 2, three dimensional, symmetric tensor.
     private
     double precision :: x00,x01,x02,x11,x12,x22
   contains
     !@ <objectMethods>
     !@   <object>tensorRank2Dimension3Symmetric</object>
     !@   <objectMethod>
     !@     <method>isZero</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if a tensor object is zero.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Destroy a tensor object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>multiply</method>
     !@     <type>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater}</type>
     !@     <arguments>\doublezero\ multiplier\argin</arguments>
     !@     <description>Multiply a tensor by a scalar.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>divide</method>
     !@     <type>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater}</type>
     !@     <arguments>\doublezero\ divisor\argin</arguments>
     !@     <description>Divide a tensor by a scalar.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>add</method>
     !@     <type>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater} tensorRank2Dimension3Symmetric2\argin</arguments>
     !@     <description>Add two tensors.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>subtract</method>
     !@     <type>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater} tensor1\argin, \textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater} tensor2\argin</arguments>
     !@     <description>Subtract one tensor from another.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>equality</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater} self\argin, \textcolor{red}{\textless double(3,3)\textgreater} matrix\argin</arguments>
     !@     <description>Return true if a tensor and a matrix are equal.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <description>Reset elements to zero.</description>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>builder</method>
     !@     <description>Build tensor object from a provided XML description.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless *type(node)\textgreater} tensorRank2Dimension3SymmetricDefinition\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dump</method>
     !@     <description>Dump the tensor object.</description>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dumpRaw</method>
     !@     <description>Dump the tensor object to binary.</description>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ fileHandle\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>readRaw</method>
     !@     <description>Read the tensor object from binary.</description>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ fileHandle\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setToUnity</method>
     !@     <description>Set all elements of the tensor object to unity.</description>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setToIdentity</method>
     !@     <description>Set a tensor object to the indentity (i.e. all diagonal elements 1, all other elements 0).</description>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serializeCount</method>
     !@     <description>Return a count of the number of properties in a serialized tensor object.</description>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serialize</method>
     !@     <description>Serialize the tensor object to an array.</description>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ tensorArray\argout</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>deserialize</method>
     !@     <description>Deserialize the tensor object from an array.</description>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ tensorArray\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>increment</method>
     !@     <description>Increment the tensor object.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater} addTensor\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>toMatrix</method>
     !@     <description>Construct a matrix from a tensor object.</description>
     !@     <type>\textcolor{red}{double(3,3)}</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>fromMatrix</method>
     !@     <description>Construct a matrix from a tensor object.</description>
     !@     <arguments>\textcolor{red}{double(3,3)} matrix</arguments>
     !@     <type>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater}</type>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>contract</method>
     !@     <description>Contract a tensor, returning ${\mathbf T}^i_i$.</description>
     !@     <arguments></arguments>
     !@     <type>\doublezero</type>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>doubleContract</method>
     !@     <description>Return the double contraction of two tensors, ${\mathbf A}^i_j {\mathbf B}^j_i$.</description>
     !@     <arguments>\textcolor{red}{\textless type(tensorRank2Dimension3Symmetric)\textgreater} tensor1</arguments>
     !@     <type>\doublezero</type>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure         :: add            => Tensor_R2_D3_Sym_Add
     procedure         :: subtract       => Tensor_R2_D3_Sym_Subtract
     procedure         :: multiply       => Tensor_R2_D3_Sym_Scalar_Multiply
     procedure         :: divide         => Tensor_R2_D3_Sym_Scalar_Divide
     procedure         :: equality       => Tensor_R2_D3_Sym_Matrix_Equality
     generic           :: operator(+)    => add
     generic           :: operator(-)    => subtract
     generic           :: operator(*)    => multiply
     generic           :: operator(/)    => divide
     generic           :: operator(==)   => equality
     procedure         :: isZero         => Tensor_R2_D3_Sym_Is_Zero
     procedure         :: destroy        => Tensor_R2_D3_Sym_Destroy
     procedure         :: setToIdentity  => Tensor_R2_D3_Sym_Set_To_Identity
     procedure         :: reset          => Tensor_R2_D3_Sym_Reset
     procedure         :: builder        => Tensor_R2_D3_Sym_Builder
     procedure         :: dump           => Tensor_R2_D3_Sym_Dump
     procedure         :: dumpRaw        => Tensor_R2_D3_Sym_Dump_Raw
     procedure         :: readRaw        => Tensor_R2_D3_Sym_Read_Raw
     procedure         :: setToUnity     => Tensor_R2_D3_Sym_Set_To_Unity
     procedure, nopass :: serializeCount => Tensor_R2_D3_Sym_Property_Count
     procedure         :: serialize      => Tensor_R2_D3_Sym_Serialize
     procedure         :: deserialize    => Tensor_R2_D3_Sym_Deserialize
     procedure         :: increment      => Tensor_R2_D3_Sym_Increment
     procedure         :: fromMatrix     => Tensor_R2_D3_Sym_From_Matrix
     procedure         :: toMatrix       => Tensor_R2_D3_Sym_To_Matrix
     procedure         :: contract       => Tensor_R2_D3_Sym_Contract
     procedure         :: doubleContract => Tensor_R2_D3_Sym_Double_Contract
  end type tensorRank2Dimension3Symmetric

  ! Identity, unitary, and null tensors.
  type(tensorRank2Dimension3Symmetric), public :: tensorIdentityR2D3Sym=tensorRank2Dimension3Symmetric(1.0d0,0.0d0,0.0d0,1.0d0,0.0d0,1.0d0)
  type(tensorRank2Dimension3Symmetric), public :: tensorUnitR2D3Sym    =tensorRank2Dimension3Symmetric(1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0)
  type(tensorRank2Dimension3Symmetric), public :: tensorNullR2D3Sym    =tensorRank2Dimension3Symmetric(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0)

contains

  subroutine Tensor_R2_D3_Sym_Destroy(self)
    !% Destroy a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} symmetric object.
    use Memory_Management
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(inout) :: self

    return
  end subroutine Tensor_R2_D3_Sym_Destroy

  subroutine Tensor_R2_D3_Sym_Builder(self,tensorDefinition)
    !% Build a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object from the given XML {\normalfont \ttfamily tensorDefinition}.
    use FoX_DOM
    use Galacticus_Error
    implicit none
    class    (tensorRank2Dimension3Symmetric), intent(inout)            :: self
    type     (node                          ), intent(in   ), pointer   :: tensorDefinition
    type     (node                          )               , pointer   :: element
    type     (nodeList                      )               , pointer   :: elementList
    character(len=3                         ), dimension(6) , parameter :: elementNames=['x00','x01','x02','x11','x12','x22']
    integer                                                             :: i

    ! Get the elements
    do i=1,6
       elementList => getElementsByTagName(tensorDefinition,elementNames(i))
       if (getLength(elementList) > 1) call Galacticus_Error_Report('Tensor_R2_D3_Sym_Builder','multiple "'//elementNames(i)//'" values specified')
       if (getLength(elementList) < 1) call Galacticus_Error_Report('Tensor_R2_D3_Sym_Builder','no "'//elementNames(i)//'" value specified'       )
       element => item(elementList,0)
       select case (elementNames(i))
       case ( 'x00' )
          call extractDataContent(element,self%x00)
       case ( 'x01' )
          call extractDataContent(element,self%x01)
       case ( 'x02' )
          call extractDataContent(element,self%x02)
       case ( 'x11' )
          call extractDataContent(element,self%x11)
       case ( 'x12' )
          call extractDataContent(element,self%x12)
       case ( 'x22' )
          call extractDataContent(element,self%x22)
       end select
    end do
    return
  end subroutine Tensor_R2_D3_Sym_Builder

  subroutine Tensor_R2_D3_Sym_Dump(self)
    !% Reset a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} symmetric object.
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class    (tensorRank2Dimension3Symmetric), intent(in   ) :: self
    character(len=12                        )                :: label
    type     (varying_string                )                :: message

    write (label,'(e12.6)') self%x00
    message='x00: '//label
    call Galacticus_Display_Message(message)
    write (label,'(e12.6)') self%x01
    message='x01: '//label
    call Galacticus_Display_Message(message)
    write (label,'(e12.6)') self%x02
    message='x02: '//label
    call Galacticus_Display_Message(message)
    write (label,'(e12.6)') self%x11
    message='x11: '//label
    call Galacticus_Display_Message(message)
    write (label,'(e12.6)') self%x12
    message='x12: '//label
    call Galacticus_Display_Message(message)
    write (label,'(e12.6)') self%x22
    message='x22: '//label
    call Galacticus_Display_Message(message)
    return
  end subroutine Tensor_R2_D3_Sym_Dump

  subroutine Tensor_R2_D3_Sym_Dump_Raw(self,fileHandle)
    !% Dump a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object to binary.
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class  (tensorRank2Dimension3Symmetric), intent(in   ) :: self
    integer                                , intent(in   ) :: fileHandle

    ! Dump the content.
    write (fileHandle) self%x00
    write (fileHandle) self%x01
    write (fileHandle) self%x02
    write (fileHandle) self%x11
    write (fileHandle) self%x12
    write (fileHandle) self%x22
    return
  end subroutine Tensor_R2_D3_Sym_Dump_Raw

  subroutine Tensor_R2_D3_Sym_Read_Raw(self,fileHandle)
    !% Read a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object from binary.
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class  (tensorRank2Dimension3Symmetric), intent(inout) :: self
    integer                                , intent(in   ) :: fileHandle

    ! Read the content.
    read (fileHandle) self%x00
    read (fileHandle) self%x01
    read (fileHandle) self%x02
    read (fileHandle) self%x11
    read (fileHandle) self%x12
    read (fileHandle) self%x22
    return
  end subroutine Tensor_R2_D3_Sym_Read_Raw

  subroutine Tensor_R2_D3_Sym_Reset(self)
    !% Reset a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(inout) :: self

    ! Zero all elements.
    self                       %x00=0.0d0
    self                       %x11=0.0d0
    self                       %x22=0.0d0
    self                       %x01=0.0d0
    self                       %x02=0.0d0
    self                       %x12=0.0d0
    return
  end subroutine Tensor_R2_D3_Sym_Reset

  subroutine Tensor_R2_D3_Sym_Set_To_Unity(self)
    !% Set a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object to unity.
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(inout) :: self

    ! Set values to unity.
    self                       %x00=1.0d0
    self                       %x11=1.0d0
    self                       %x22=1.0d0
    self                       %x01=1.0d0
    self                       %x02=1.0d0
    self                       %x12=1.0d0
    return
  end subroutine Tensor_R2_D3_Sym_Set_To_Unity

  subroutine Tensor_R2_D3_Sym_Set_To_Identity(self)
    !% Set a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object to the identity matrix.
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(inout) :: self

    ! Set values to unity.
    self                       %x00=1.0d0
    self                       %x11=1.0d0
    self                       %x22=1.0d0
    self                       %x01=0.0d0
    self                       %x02=0.0d0
    self                       %x12=0.0d0
    return
  end subroutine Tensor_R2_D3_Sym_Set_To_Identity

  logical function Tensor_R2_D3_Sym_Is_Zero(self)
    !% Test whether a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object is zero.
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(in) :: self

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
  end function Tensor_R2_D3_Sym_Is_Zero
  
  function Tensor_R2_D3_Sym_Add(tensor1,tensor2)
    !% Add two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
    implicit none
    type (tensorRank2Dimension3Symmetric)                       :: Tensor_R2_D3_Sym_Add
    class(tensorRank2Dimension3Symmetric), intent(in)           :: tensor1
    class(tensorRank2Dimension3Symmetric), intent(in), optional :: tensor2
    
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
  end function Tensor_R2_D3_Sym_Add

  subroutine Tensor_R2_D3_Sym_Increment(self,increment)
    !% Increment a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(inout) :: self
    class(tensorRank2Dimension3Symmetric), intent(in   ) :: increment

    self%x00=self%x00+increment%x00
    self%x01=self%x01+increment%x01
    self%x02=self%x02+increment%x02
    self%x11=self%x11+increment%x11
    self%x12=self%x12+increment%x12
    self%x22=self%x22+increment%x22
    return
  end subroutine Tensor_R2_D3_Sym_Increment

  function Tensor_R2_D3_Sym_Subtract(tensor1,tensor2)
    !% Subtract two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
    implicit none
    type (tensorRank2Dimension3Symmetric)                       :: Tensor_R2_D3_Sym_Subtract
    class(tensorRank2Dimension3Symmetric), intent(in)           :: tensor1
    class(tensorRank2Dimension3Symmetric), intent(in), optional :: tensor2

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
  end function Tensor_R2_D3_Sym_Subtract

  function Tensor_R2_D3_Sym_Scalar_Multiply(tensor1,multiplier)
    !% Multiply a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object by a scalar.
    implicit none
    type (tensorRank2Dimension3Symmetric)             :: Tensor_R2_D3_Sym_Scalar_Multiply
    class(tensorRank2Dimension3Symmetric), intent(in) :: tensor1
    double precision , intent(in) :: multiplier

    Tensor_R2_D3_Sym_Scalar_Multiply%x00=tensor1%x00*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x01=tensor1%x01*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x02=tensor1%x02*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x11=tensor1%x11*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x12=tensor1%x12*multiplier
    Tensor_R2_D3_Sym_Scalar_Multiply%x22=tensor1%x22*multiplier
    return
  end function Tensor_R2_D3_Sym_Scalar_Multiply

  function Tensor_R2_D3_Sym_Scalar_Multiply_Switched(multiplier,tensor1)
    !% Multiply a scalar by a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
    implicit none
    type            (tensorRank2Dimension3Symmetric)                :: Tensor_R2_D3_Sym_Scalar_Multiply_Switched
    type            (tensorRank2Dimension3Symmetric), intent(in   ) :: tensor1
    double precision                                , intent(in   ) :: multiplier

    Tensor_R2_D3_Sym_Scalar_Multiply_Switched=Tensor_R2_D3_Sym_Scalar_Multiply(tensor1,multiplier)
    return
  end function Tensor_R2_D3_Sym_Scalar_Multiply_Switched

  function Tensor_R2_D3_Sym_Max(tensor1,tensor2)
    !% Return an element-by-element {\normalfont \ttfamily max()} on two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects.
    implicit none
    type(tensorRank2Dimension3Symmetric)                :: Tensor_R2_D3_Sym_Max
    type(tensorRank2Dimension3Symmetric), intent(in   ) :: tensor1,tensor2

    Tensor_R2_D3_Sym_Max%x00=max(tensor1%x00,tensor2%x00)
    Tensor_R2_D3_Sym_Max%x01=max(tensor1%x01,tensor2%x01)
    Tensor_R2_D3_Sym_Max%x02=max(tensor1%x02,tensor2%x02)
    Tensor_R2_D3_Sym_Max%x11=max(tensor1%x11,tensor2%x11)
    Tensor_R2_D3_Sym_Max%x12=max(tensor1%x12,tensor2%x12)
    Tensor_R2_D3_Sym_Max%x22=max(tensor1%x22,tensor2%x22)
    return
  end function Tensor_R2_D3_Sym_Max

  function Tensor_R2_D3_Sym_Scalar_Divide(tensor1,divisor)
    !% Multiply a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object by a scalar.
    implicit none
    type            (tensorRank2Dimension3Symmetric)                :: Tensor_R2_D3_Sym_Scalar_Divide
    class           (tensorRank2Dimension3Symmetric), intent(in   ) :: tensor1
    double precision                                , intent(in   ) :: divisor
 
    Tensor_R2_D3_Sym_Scalar_Divide%x00=tensor1%x00/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x01=tensor1%x01/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x02=tensor1%x02/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x11=tensor1%x11/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x12=tensor1%x12/divisor
    Tensor_R2_D3_Sym_Scalar_Divide%x22=tensor1%x22/divisor
    return
  end function Tensor_R2_D3_Sym_Scalar_Divide

  double precision function Tensor_R2_D3_Sym_Double_Contract(self,tensor1)
    !% Find the double contraction of two {\normalfont \ttfamily tensorRank2Dimension3Symmetric} objects, ${\mathbf A}:{\mathbf B}$.
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(in   ) :: self,tensor1

    Tensor_R2_D3_Sym_Double_Contract=  &
         & +      self%x00*tensor1%x00 &
         & +      self%x11*tensor1%x11 &
         & +      self%x22*tensor1%x22 &
         & +2.0d0*self%x01*tensor1%x01 &
         & +2.0d0*self%x02*tensor1%x02 &
         & +2.0d0*self%x12*tensor1%x12
    return
  end function Tensor_R2_D3_Sym_Double_Contract

  double precision function Tensor_R2_D3_Sym_Contract(self)
    !% Return the contraction (trace) of a {\normalfont \ttfamily tensorRank2Dimension3Symmetric}.
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(in   ) :: self

    Tensor_R2_D3_Sym_Contract=self%x00+self%x11+self%x22
    return
  end function Tensor_R2_D3_Sym_Contract

  integer function Tensor_R2_D3_Sym_Property_Count()
    !% Return the number of properties required to track a rank 2, 3 dimensional, symmetric tensor. This is equal to 6.
    implicit none
    integer, parameter :: propertyCount=6

    Tensor_R2_D3_Sym_Property_Count=propertyCount
    return
  end function Tensor_R2_D3_Sym_Property_Count

  subroutine Tensor_R2_D3_Sym_Deserialize(self,tensorArray)
    !% Pack an array into a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} symmetric structure.
    implicit none
    class           (tensorRank2Dimension3Symmetric), intent(inout)               :: self
    double precision                                , intent(in   ), dimension(6) :: tensorArray

    self%x00=tensorArray(1)
    self%x01=tensorArray(2)
    self%x02=tensorArray(3)
    self%x11=tensorArray(4)
    self%x12=tensorArray(5)
    self%x22=tensorArray(6)
    return
  end subroutine Tensor_R2_D3_Sym_Deserialize

  subroutine Tensor_R2_D3_Sym_Serialize(self,tensorArray)
    !% Pack a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} into an array.
    implicit none
    double precision                                , intent(  out), dimension(:) :: tensorArray(6)
    class           (tensorRank2Dimension3Symmetric), intent(in   )               :: self

    ! Place tensor into array.
    tensorArray(1)=self%x00
    tensorArray(2)=self%x01
    tensorArray(3)=self%x02
    tensorArray(4)=self%x11
    tensorArray(5)=self%x12
    tensorArray(6)=self%x22
    return
  end subroutine Tensor_R2_D3_Sym_Serialize

  subroutine Tensor_R2_D3_Sym_From_Matrix(self,matrix)
    !% Construct a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object from a matrix.
    use Galacticus_Error
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(inout)                 :: self
    double precision                     , intent(in   ), dimension(3,3) :: matrix
    
    ! Test for symmetry
    if     (                            &
         &   matrix(1,2) /= matrix(2,1) &
         &  .or.                        &
         &   matrix(1,3) /= matrix(3,1) &
         &  .or.                        &
         &   matrix(2,3) /= matrix(3,2) &
         & ) call Galacticus_Error_Report('Tensor_R2_D3_Sym_From_Matrix','supplied matrix is not symmetric')
    ! Set the values of the tensor object.
    self%x00=matrix(1,1)
    self%x01=matrix(1,2)
    self%x02=matrix(1,3)
    self%x11=matrix(2,2)
    self%x12=matrix(2,3)
    self%x22=matrix(3,3)
    return
  end subroutine Tensor_R2_D3_Sym_From_Matrix

  function Tensor_R2_D3_Sym_To_Matrix(self)
    !% Construct a matrix from a {\normalfont \ttfamily tensorRank2Dimension3Symmetric}.
    implicit none
    double precision                     , dimension(3,3) :: Tensor_R2_D3_Sym_To_Matrix
    class(tensorRank2Dimension3Symmetric), intent(in   )  :: self

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
  end function Tensor_R2_D3_Sym_To_Matrix

  subroutine Tensor_R2_D3_Sym_Assign_To(tensor,matrix)
    !% Assign a matrix to a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} object.
    use Galacticus_Error
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(inout)                 :: tensor
    double precision                     , intent(in   ), dimension(3,3) :: matrix
    
    call tensor%fromMatrix(matrix)
    return
  end subroutine Tensor_R2_D3_Sym_Assign_To

  subroutine Tensor_R2_D3_Sym_Assign_From(matrix,tensor)
    !% Assign a {\normalfont \ttfamily tensorRank2Dimension3Symmetric} to a matrix.
    implicit none
    double precision                     , intent(  out), dimension(3,3) :: matrix
    class(tensorRank2Dimension3Symmetric), intent(in   )                 :: tensor

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
  end subroutine Tensor_R2_D3_Sym_Assign_From

  logical function Tensor_R2_D3_Sym_Matrix_Equality(self,matrix)
    !% Return true if the supplied tensor and matrix are equal.
    implicit none
    class(tensorRank2Dimension3Symmetric), intent(in   )                 :: self
    double precision                     , intent(in   ), dimension(3,3) :: matrix

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
  end function Tensor_R2_D3_Sym_Matrix_Equality

end module Tensors
