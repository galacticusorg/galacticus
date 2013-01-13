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

!% Contains a module which defines a {\tt table} class with optimized interpolation operators.

module Tables
  !% Defines a {\tt table} class with optimized interpolation operators.
  use FGSL
  private
  public :: table,table1D,table1DLinearLinear,table1DLogarithmicLinear,table1DLinearCSpline,table1DLogarithmicCSpline

  type, abstract :: table
     !% Basic table type.
   contains
     procedure(Table_Destroy), deferred :: destroy
  end type table

  interface
     subroutine Table_Destroy(self)
       !% Interface to {\tt table} destructor.
       import table
       implicit none
       class(table), intent(inout) :: self
     end subroutine Table_Destroy
  end interface

  type, abstract, extends(table) :: table1D
     !% Basic table type.
     integer                                       :: xCount
     double precision, allocatable, dimension(:  ) :: xv
     double precision, allocatable, dimension(:,:) :: yv
   contains
     procedure(Table1D_Interpolate ), deferred :: interpolate
     procedure(Table1D_Interpolate ), deferred :: interpolateGradient
     procedure                                 :: destroy             => Table_1D_Destroy
     procedure                                 :: reverse             => Table_1D_Reverse
     procedure                                 :: isMonotonic         => Table1D_Is_Monotonic
     procedure                                 :: x                   => Table1D_X
     procedure                                 :: y                   => Table1D_Y
     procedure                                 :: xs                  => Table1D_Xs
     procedure                                 :: ys                  => Table1D_Ys
  end type table1D
  
  interface
     double precision function Table1D_Interpolate(self,x,table)
       !% Interface to {\tt table} interpolator.
       import table1D
       implicit none
       class(table1D)  , intent(inout)           :: self
       double precision, intent(in   )           :: x
       integer         , intent(in   ), optional :: table
     end function Table1D_Interpolate
  end interface
  
  type, extends(table1D) :: table1DGeneric
     !% Table type supporting generic one dimensional tables.
     type   (fgsl_interp      ) :: interpolator
     type   (fgsl_interp_accel) :: accelerator
     logical                    :: reset
   contains
     procedure :: create              => Table_Generic_1D_Create
     procedure :: destroy             => Table_Generic_1D_Destroy
     procedure ::                        Table_Generic_1D_Populate
     procedure ::                        Table_Generic_1D_Populate_Single
     generic   :: populate            => Table_Generic_1D_Populate, Table_Generic_1D_Populate_Single
     procedure :: interpolate         => Table_Generic_1D_Interpolate
     procedure :: interpolateGradient => Table_Generic_1D_Interpolate_Gradient
  end type table1DGeneric

  type, extends(table1D) :: table1DLinearLinear
     !% Table type supporting one dimensional table with linear spacing in $x$.
     double precision :: inverseDeltaX,xPrevious,yPrevious,dxPrevious,dyPrevious
     integer          :: tablePrevious,dTablePrevious
   contains
     procedure :: create              => Table_Linear_1D_Create
     procedure ::                        Table_Linear_1D_Populate
     procedure ::                        Table_Linear_1D_Populate_Single
     generic   :: populate            => Table_Linear_1D_Populate, Table_Linear_1D_Populate_Single
     procedure :: interpolate         => Table_Linear_1D_Interpolate
     procedure :: interpolateGradient => Table_Linear_1D_Interpolate_Gradient
  end type table1DLinearLinear
  
  type, extends(table1DLinearLinear) :: table1DLogarithmicLinear
     !% Table type supporting one dimensional table with logarithmic spacing in $x$.
   contains
     procedure :: create              => Table_Logarithmic_1D_Create
     procedure :: interpolate         => Table_Logarithmic_1D_Interpolate
     procedure :: interpolateGradient => Table_Logarithmic_1D_Interpolate_Gradient
     procedure :: x                   => Table_Logarithmic_1D_X
     procedure :: xs                  => Table_Logarithmic_1D_Xs
  end type table1DLogarithmicLinear

  type, extends(table1D) :: table1DLinearCSpline
     !% Table type supporting one dimensional table with linear spacing in $x$ and cubic spline interpolation.
     double precision, allocatable, dimension(:,:) :: sv
     integer                                       :: iPrevious,tablePrevious,dTablePrevious
     double precision                              :: deltaX,inverseDeltaX
     double precision                              :: xPrevious,yPrevious,dxPrevious,dyPrevious,aPrevious,bPrevious,cPrevious,dPrevious
   contains
     procedure :: create              => Table_Linear_CSpline_1D_Create
     procedure ::                        Table_Linear_CSpline_1D_Populate
     procedure ::                        Table_Linear_CSpline_1D_Populate_Single
     generic   :: populate            => Table_Linear_CSpline_1D_Populate, Table_Linear_CSpline_1D_Populate_Single
     procedure :: interpolate         => Table_Linear_CSpline_1D_Interpolate
     procedure :: interpolateGradient => Table_Linear_CSpline_1D_Interpolate_Gradient
  end type table1DLinearCSpline
  
  type, extends(table1DLinearCSpline) :: table1DLogarithmicCSpline
     !% Table type supporting one dimensional table with logarithmic spacing in $x$ and cubic spline interpolation.
   contains
     procedure :: create              => Table_Logarithmic_CSpline_1D_Create
     procedure :: interpolate         => Table_Logarithmic_CSpline_1D_Interpolate
     procedure :: interpolateGradient => Table_Logarithmic_CSpline_1D_Interpolate_Gradient
     procedure :: x                   => Table_Logarithmic_CSpline_1D_X
     procedure :: xs                  => Table_Logarithmic_CSpline_1D_Xs
  end type table1DLogarithmicCSpline

contains

  subroutine Table_1D_Destroy(self)
    !% Destroy a 1-D table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class(table1D), intent(inout) :: self
    
    if (allocated(self%xv)) call Dealloc_Array(self%xv)
    if (allocated(self%yv)) call Dealloc_Array(self%yv)
    return
  end subroutine Table_1D_Destroy
  
  double precision function Table1D_X(self,i)
    !% Return the {\tt i}$^{\rm th}$ $x$-value for a 1D table.
    use Galacticus_Error
    implicit none
    class  (table1D), intent(inout) :: self
    integer         , intent(in   ) :: i
    integer                         :: ii
    
    ii=i
    if (ii < 0) ii=ii+size(self%xv)+1
    Table1D_X=self%xv(ii)
    return
  end function Table1D_X

  function Table1D_Xs(self)
    !% Return the $x$-values for a 1D table.
    use Galacticus_Error
    implicit none
    class           (table1D), intent(in   )             :: self
    double precision         , dimension(size(self%xv))  :: Table1D_Xs

    Table1D_Xs=self%xv
    return
  end function Table1D_Xs

  double precision function Table1D_Y(self,i,table)
    !% Return the {\tt i}$^{\rm th}$ $y$-value for a 1D table.
    use Galacticus_Error
    implicit none
    class  (table1D), intent(in   )           :: self
    integer         , intent(in   )           :: i
    integer         , intent(in   ), optional :: table
    integer                                   :: ii,tableActual
    
    tableActual=1
    if (present(table)) tableActual=table
    ii=i
    if (ii < 0) ii=ii+size(self%xv)+1
    Table1D_Y=self%yv(ii,tableActual)
    return
  end function Table1D_Y

  function Table1D_Ys(self)
    !% Return the $y$-values for a 1D table.
    use Galacticus_Error
    implicit none
    class           (table1D), intent(in   )                                      :: self
    double precision         , dimension(size(self%yv,dim=1),size(self%yv,dim=2)) :: Table1D_Ys

    Table1D_Ys=self%yv
    return
  end function Table1D_Ys

  subroutine Table_1D_Reverse(self,reversedSelf,table)
    !% Reverse a 1D table (i.e. swap $x$ and $y$ components). Optionally allows specification of
    !% which $y$ table to swap with.
    use Array_Utilities
    use Galacticus_Error
    implicit none
    class  (table1D), intent(in   )              :: self
    class  (table1D), intent(inout), allocatable :: reversedSelf
    integer         , intent(in   ), optional    :: table
    integer                                      :: tableActual,i

    tableActual=1
    if (present(table)) tableActual=table
    if (.not.Array_Is_Monotonic(self%yv(:,tableActual))) call Galacticus_Error_Report('Table_1D_Reverse','reversed table would not be monotonic')
    if (allocated(reversedSelf)) deallocate(reversedSelf)
    allocate(table1DGeneric :: reversedSelf)
    select type (reversedSelf)
    type is (table1DGeneric)
       reversedSelf%reset                =.true.
       reversedSelf%xCount               =self        %xCount
       reversedSelf%yv                   =self        %ys    (             )
       reversedSelf%xv                   =reversedSelf%yv    (:,tableActual)
       reversedSelf%yv    (:,tableActual)=self        %xs    (             )
       ! If the table is monotonically decreasing when reversed, we will also need to switch the order.
       if (Array_Is_Monotonic(reversedSelf%xv,direction=directionDecreasing)) then
          reversedSelf%xv=Array_Reverse(reversedSelf%xv)
          do i=1,size(reversedSelf%yv,dim=2)
             reversedSelf%yv(:,i)=Array_Reverse(reversedSelf%yv(:,i))
          end do
       end if
    end select
    return
  end subroutine Table_1D_Reverse

  logical function Table1D_Is_Monotonic(self,direction,allowEqual,table)
    !% Return true if a 1D table is monotonic. Optionally allows specification of the direction,
    !% and whether or not equal elements are allowed for monotonicity.
    use Array_Utilities
    implicit none
    class  (table1D), intent(in   )           :: self
    integer         , intent(in   ), optional :: direction,table
    logical         , intent(in   ), optional :: allowEqual
    integer                                   :: tableActual

    tableActual=1
    if (present(table)) tableActual=table
    Table1D_Is_Monotonic=Array_Is_Monotonic(self%yv(:,tableActual),direction,allowEqual)
    return
  end function Table1D_Is_Monotonic

  subroutine Table_Generic_1D_Create(self,x,tableCount)
    !% Create a 1-D generic table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class           (table1DGeneric), intent(inout)               :: self
    double precision                , intent(in   ), dimension(:) :: x
    integer                         , intent(in   ), optional     :: tableCount
    integer                                                       :: tableCountActual
    
    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    
    ! Allocate arrays and construct the x-range.
    self%xCount=size(x)
    call Alloc_Array(self%xv,[size(x)                 ])
    call Alloc_Array(self%yv,[size(x),tableCountActual])
    self%xv   =x
    self%reset=.true.
    return
  end subroutine Table_Generic_1D_Create

  subroutine Table_Generic_1D_Destroy(self)
    !% Destroy a generic 1-D table.
    use Numerical_Interpolation
    implicit none
    class(table1DGeneric), intent(inout) :: self
    
    call Table_1D_Destroy(self)
    call Interpolate_Done(self%interpolator,self%accelerator,self%reset)
    return
  end subroutine Table_Generic_1D_Destroy

  subroutine Table_Generic_1D_Populate(self,y,table)
    !% Populate a 1-D generic table.
    use Galacticus_Error
    implicit none
    class           (table1DGeneric), intent(inout)               :: self
    double precision                , intent(in   ), dimension(:) :: y
    integer                         , intent(in   ), optional     :: table
    integer                                                       :: tableActual
    
    ! Validate the input.
    if (.not.allocated(self%yv)       ) call Galacticus_Error_Report("Table_Generic_1D_Populate","create the table before populating it")
    if (size(self%yv,dim=1) /= size(y)) call Galacticus_Error_Report("Table_Generic_1D_Populate","provided y array is of wrong size"    )

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    
    ! Store the y values.
    self%yv(:,tableActual)=y
    return
  end subroutine Table_Generic_1D_Populate
  
  subroutine Table_Generic_1D_Populate_Single(self,y,i,table)
    !% Populate a single element of a 1-D generic table.
    use Galacticus_Error
    implicit none
    class           (table1DGeneric), intent(inout)           :: self
    double precision                , intent(in   )           :: y
    integer                         , intent(in   )           :: i
    integer                         , intent(in   ), optional :: table
    integer                                                   :: tableActual
    
    ! Validate the input.
    if (.not.allocated(self%yv)           ) call Galacticus_Error_Report("Table_Generic_1D_Populate_Single","create the table before populating it")
    if (i < 1 .or. i > size(self%yv,dim=1)) call Galacticus_Error_Report("Table_Generic_1D_Populate_Single","provided i value is out of bounds"    )
    
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    
    ! Store the y values.
    self%yv(i,tableActual)=y
    return
  end subroutine Table_Generic_1D_Populate_Single
  
  double precision function Table_Generic_1D_Interpolate(self,x,table)
    !% Perform generic interpolation in a generic 1D table.
    use Numerical_Interpolation
    implicit none
    class           (table1DGeneric), intent(inout)           :: self
    double precision                , intent(in   )           :: x
    integer                         , intent(in   ), optional :: table
    integer                                                   :: tableActual

    tableActual=1
    if (present(table)) tableActual=table
    Table_Generic_1D_Interpolate=Interpolate(size(self%xv),self%xv,self%yv(:,tableActual),self%interpolator,self%accelerator,x,reset=self%reset)
    return
  end function Table_Generic_1D_Interpolate
  
  double precision function Table_Generic_1D_Interpolate_Gradient(self,x,table)
    !% Perform generic interpolation in a generic 1D table.
    use Numerical_Interpolation
    implicit none
    class           (table1DGeneric), intent(inout)           :: self
    double precision                , intent(in   )           :: x
    integer                         , intent(in   ), optional :: table
    integer                                                   :: tableActual
    
    tableActual=1
    if (present(table)) tableActual=table
    Table_Generic_1D_Interpolate_Gradient=Interpolate_Derivative(size(self%xv),self%xv,self%yv(:,tableActual),self%interpolator,self%accelerator,x,reset=self%reset)
    return
  end function Table_Generic_1D_Interpolate_Gradient
  
  subroutine Table_Linear_1D_Create(self,xMinimum,xMaximum,xCount,tableCount)
    !% Create a 1-D linear table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class           (table1DLinearLinear), intent(inout)           :: self
    double precision                     , intent(in   )           :: xMinimum,xMaximum
    integer                              , intent(in   )           :: xCount
    integer                              , intent(in   ), optional :: tableCount
    integer                                                        :: tableCountActual
    
    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    
    ! Allocate arrays and construct the x-range.
    self%xCount=xCount
    call Alloc_Array(self%xv,[xCount                 ])
    call Alloc_Array(self%yv,[xCount,tableCountActual])
    self%xv           =Make_Range(xMinimum,xMaximum,xCount,rangeType=rangeTypeLinear)
    self%inverseDeltaX=1.0d0/(self%xv(2)-self%xv(1))
    return
  end subroutine Table_Linear_1D_Create

  subroutine Table_Linear_1D_Populate(self,y,table)
    !% Populate a 1-D linear table.
    use Galacticus_Error
    implicit none
    class           (table1DLinearLinear), intent(inout)               :: self
    double precision                     , intent(in   ), dimension(:) :: y
    integer                              , intent(in   ), optional     :: table
    integer                                                            :: tableActual
    
    ! Validate the input.
    if (.not.allocated(self%yv)       ) call Galacticus_Error_Report("Table_Linear_1D_Populate","create the table before populating it")
    if (size(self%yv,dim=1) /= size(y)) call Galacticus_Error_Report("Table_Linear_1D_Populate","provided y array is of wrong size"    )
    
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    
    ! Reset all previously stored values.
    self% tablePrevious=-1
    self%dTablePrevious=-1
    
    ! Store the y values.
    self%yv(:,tableActual)=y
    return
  end subroutine Table_Linear_1D_Populate
  
  subroutine Table_Linear_1D_Populate_Single(self,y,i,table)
    !% Populate a single element of a 1-D linear table.
    use Galacticus_Error
    implicit none
    class           (table1DLinearLinear), intent(inout)           :: self
    double precision                     , intent(in   )           :: y
    integer                              , intent(in   )           :: i
    integer                              , intent(in   ), optional :: table
    integer                                                        :: tableActual
    
    ! Validate the input.
    if (.not.allocated(self%yv)           ) call Galacticus_Error_Report("Table_Linear_1D_Populate_Single","create the table before populating it")
    if (i < 1 .or. i > size(self%yv,dim=1)) call Galacticus_Error_Report("Table_Linear_1D_Populate_Single","provided i value is out of bounds"    )
    
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    
    ! Reset all previously stored values.
    self% tablePrevious=-1
    self%dTablePrevious=-1
    
    ! Store the y values.
    self%yv(i,tableActual)=y
    return
  end subroutine Table_Linear_1D_Populate_Single
  
  double precision function Table_Linear_1D_Interpolate(self,x,table)
    !% Perform linear interpolation in a linear 1D table.
    implicit none
    class           (table1DLinearLinear), intent(inout)           :: self
    double precision                     , intent(in   )           :: x
    integer                              , intent(in   ), optional :: table
    integer                                                        :: i,tableActual
    double precision                                               :: h

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%xPrevious .or. tableActual /= self%tablePrevious) then    
       ! Determine the location in the table.
       if      (x <  self%xv(          1)) then
          i=1
       else if (x >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=max(min(int((x-self%xv(1))*self%inverseDeltaX)+1,self%xCount-1),1)
       end if       
       ! Compute the interpolating factor.
       h=(x-self%xv(i))*self%inverseDeltaX    
       ! Interpolate in the table.
       self%tablePrevious=tableActual
       self%    yPrevious=self%yv(i,tableActual)*(1.0d0-h)+self%yv(i+1,tableActual)*h
    end if
    ! Interpolate in the table.
    Table_Linear_1D_Interpolate=self%yPrevious
    return
  end function Table_Linear_1D_Interpolate
  
  double precision function Table_Linear_1D_Interpolate_Gradient(self,x,table)
    !% Perform linear interpolation in a linear 1D table.
    implicit none
    class           (table1DLinearLinear), intent(inout)           :: self
    double precision                     , intent(in   )           :: x
    integer                              , intent(in   ), optional :: table
    integer                                                        :: i,tableActual
    
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then    
       ! Determine the location in the table.
       if      (x <  self%xv(          1)) then
          i=1
       else if (x >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=int((x-self%xv(1))*self%inverseDeltaX)+1
       end if
       ! Interpolate in the table.
       self%dTablePrevious=tableActual
       self%    dyPrevious=(self%yv(i+1,tableActual)-self%yv(i,tableActual))*self%inverseDeltaX
    end if
    Table_Linear_1D_Interpolate_Gradient=self%dyPrevious
    return
  end function Table_Linear_1D_Interpolate_Gradient
  
  subroutine Table_Logarithmic_1D_Create(self,xMinimum,xMaximum,xCount,tableCount)
    !% Create a 1-D logarithmic table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class           (table1DLogarithmicLinear), intent(inout)           :: self
    double precision                          , intent(in   )           :: xMinimum,xMaximum
    integer                                   , intent(in   )           :: xCount
    integer                                   , intent(in   ), optional :: tableCount
    
    ! Call the creator for linear tables with the logarithms of the input x range.
    call self%table1DLinearLinear%create(log(xMinimum),log(xMaximum),xCount,tableCount)
    return
  end subroutine Table_Logarithmic_1D_Create
  
  double precision function Table_Logarithmic_1D_X(self,i)
    !% Return the {\tt i}$^{\rm th}$ $x$-value for a logarithmic 1D table.
    use Galacticus_Error
    implicit none
    class  (table1DLogarithmicLinear), intent(inout) :: self
    integer                          , intent(in   ) :: i
    
    Table_Logarithmic_1D_X=exp(self%table1DLinearLinear%x(i))
    return
  end function Table_Logarithmic_1D_X
  
  function Table_Logarithmic_1D_Xs(self)
    !% Return the $x$-values for a 1D table.
    use Galacticus_Error
    implicit none
    class           (table1DLogarithmicLinear), intent(in   )             :: self
    double precision                          , dimension(size(self%xv))  :: Table_Logarithmic_1D_Xs

    Table_Logarithmic_1D_Xs=exp(self%table1DLinearLinear%xs())
    return
  end function Table_Logarithmic_1D_Xs

  double precision function Table_Logarithmic_1D_Interpolate(self,x,table)
    !% Perform linear interpolation in a logarithmic 1D table.
    implicit none
    class           (table1DLogarithmicLinear), intent(inout)           :: self
    double precision                          , intent(in   )           :: x
    integer                                   , intent(in   ), optional :: table
    
    Table_Logarithmic_1D_Interpolate=self%table1DLinearLinear%interpolate(log(x),table)
    return
  end function Table_Logarithmic_1D_Interpolate
  
  double precision function Table_Logarithmic_1D_Interpolate_Gradient(self,x,table)
    !% Perform linear interpolation in a logarithmic 1D table.
    implicit none
    class           (table1DLogarithmicLinear), intent(inout)           :: self
    double precision                          , intent(in   )           :: x
    integer                                   , intent(in   ), optional :: table
    
    Table_Logarithmic_1D_Interpolate_Gradient=self%table1DLinearLinear%interpolateGradient(log(x),table)/x
    return
  end function Table_Logarithmic_1D_Interpolate_Gradient
  
  subroutine Table_Linear_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount)
    !% Create a 1-D linear table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class           (table1DLinearCSpline), intent(inout)           :: self
    double precision                      , intent(in   )           :: xMinimum,xMaximum
    integer                               , intent(in   )           :: xCount
    integer                               , intent(in   ), optional :: tableCount
    integer                                                        :: tableCountActual
    
    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    
    ! Allocate arrays and construct the x-range.
    self%xCount=xCount
    call Alloc_Array(self%xv,[xCount                 ])
    call Alloc_Array(self%yv,[xCount,tableCountActual])
    call Alloc_Array(self%sv,[xCount,tableCountActual])
    self%xv           =Make_Range(xMinimum,xMaximum,xCount,rangeType=rangeTypeLinear)
    self%       deltaX=self%xv(2)-self%xv(1)
    self%inverseDeltaX=1.0d0/self%deltaX
    return
  end subroutine Table_Linear_CSpline_1D_Create

  subroutine Table_Linear_CSpline_1D_Destroy(self)
    !% Destroy a linear cubic-sline 1-D table.
    use Memory_Management
    implicit none
    class(table1DLinearCSpline), intent(inout) :: self
    
    call Table_1D_Destroy(self)
    if (allocated(self%sv)) call Dealloc_Array(self%sv)
    return
  end subroutine Table_Linear_CSpline_1D_Destroy

  subroutine Table_Linear_CSpline_1D_Populate(self,y,table,computeSpline)
    !% Populate a 1-D linear table.
    use Galacticus_Error
    implicit none
    class           (table1DLinearCSpline), intent(inout)               :: self
    double precision                      , intent(in   ), dimension(:) :: y
    integer                               , intent(in   ), optional     :: table
    logical                               , intent(in   ), optional     :: computeSpline
    integer                                                             :: tableActual
    logical                                                             :: computeSplineActual

    ! Validate the input.
    if (.not.(allocated(self%yv).and.allocated(self%sv))) call Galacticus_Error_Report("Table_Linear_CSpline_1D_Populate","create the table before populating it")
    if (size(self%yv,dim=1) /= size(y)                  ) call Galacticus_Error_Report("Table_Linear_CSpline_1D_Populate","provided y array is of wrong size"    )
    
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    
    ! Store the y values.
    self%yv(:,tableActual)=y

    ! Compute the spline interpolation for this table.
    computeSplineActual=.true.
    if (present(computeSpline)) computeSplineActual=computeSpline
    if (computeSplineActual) call Table_Linear_CSpline_1D_Compute_Spline(self,tableActual)
    return
  end subroutine Table_Linear_CSpline_1D_Populate
  
  subroutine Table_Linear_CSpline_1D_Populate_Single(self,y,i,table,computeSpline)
    !% Populate a single element of a 1-D linear table.
    use Galacticus_Error
    implicit none
    class           (table1DLinearCSpline), intent(inout)           :: self
    double precision                      , intent(in   )           :: y
    integer                               , intent(in   )           :: i
    integer                               , intent(in   ), optional :: table
    logical                               , intent(in   ), optional :: computeSpline
    integer                                                         :: tableActual
    logical                                                         :: computeSplineActual
    
    ! Validate the input.
    if (.not.(allocated(self%yv).and.allocated(self%sv))) call Galacticus_Error_Report("Table_Linear_CSpline_1D_Populate_Single","create the table before populating it")
    if (i < 1 .or. i > size(self%yv,dim=1)              ) call Galacticus_Error_Report("Table_Linear_CSpline_1D_Populate_Single","provided i value is out of bounds"    )
    
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    
    ! Store the y values.
    self%yv(i,tableActual)=y

    ! Compute the spline interpolation for this table.
    computeSplineActual=.true.
    if (present(computeSpline)) computeSplineActual=computeSpline
    if (computeSplineActual) call Table_Linear_CSpline_1D_Compute_Spline(self,tableActual)
    return
  end subroutine Table_Linear_CSpline_1D_Populate_Single
  
  subroutine Table_Linear_CSpline_1D_Compute_Spline(self,table)
    !% Compute the interpolating spline factors for a 1-D linear spline.
    implicit none
    type            (table1DLinearCSpline), intent(inout)               :: self
    integer                               , intent(in   )               :: table
    double precision                      , allocatable  , dimension(:) :: b,u,v
    integer                                                             :: i
    double precision                                                    :: h

    ! Reset all previously stored values.
    self% tablePrevious=-1
    self%dTablePrevious=-1
    self%     iPrevious=-1
    ! Compute the spline interpolation factors.
    allocate(b(size(self%xv)))
    allocate(u(size(self%xv)))
    allocate(v(size(self%xv)))
    h=self%deltaX
    forall(i=1:size(self%xv)-1)
       b(i)=(self%yv(i+1,table)-self%yv(i,table))/h
    end forall
    u(2)=4.0d0*h
    v(2)=6.0d0*(b(2)-b(1))
    do i=3,size(self%xv)-1
       u(i)=4.0d0*h-h**2/u(i-1)
       v(i)=6.0d0*(b(i)-b(i-1))-h*v(i-1)/u(i-1)
    end do
    self%sv(size(self%xv),table)=0.0d0
    do i=size(self%xv)-1,2,-1
       self%sv(i,table)=(v(i)-h*self%sv(i+1,table))/u(i)
    end do
    self%sv(            1,table)=0.0d0
    deallocate(b,u,v)
    return
  end subroutine Table_Linear_CSpline_1D_Compute_Spline

  subroutine Table_Linear_CSpline_1D_Coefficients(self,table,x,i,a,b,c,d,dx)
    !% Compute coefficients for a spline interpolation.
    implicit none
    type            (table1DLinearCSpline), intent(inout) :: self
    double precision                      , intent(in   ) :: x
    integer                               , intent(in   ) :: table,i
    double precision                      , intent(  out) :: a,b,c,d,dx

    if (i /= self%iPrevious) then
       self%aPrevious=                      self%yv(i  ,table)
       self%bPrevious=-self%       deltaX*  self%sv(i+1,table)/6.0d0 &
            &         -self%       deltaX*  self%sv(i  ,table)/3.0d0 &
            &         +self%inverseDeltaX*( self%yv(i+1,table)       &
            &                              -self%yv(i  ,table)       &
            &                             )
       self%cPrevious=                      self%sv(i  ,table)/2.0d0
       self%dPrevious= self%inverseDeltaX*(                          &
            &                               self%sv(i+1,table)       &
            &                              -self%sv(i  ,table)       &
            &                             )                   /6.0d0
       self%iPrevious=i
    end if
    a =  self%aPrevious
    b =  self%bPrevious
    c =  self%cPrevious
    d =  self%dPrevious
    dx=x-self%xv       (i)
    return
  end subroutine Table_Linear_CSpline_1D_Coefficients

  double precision function Table_Linear_CSpline_1D_Interpolate(self,x,table)
    !% Perform linear interpolation in a linear 1D table.
    implicit none
    class           (table1DLinearCSpline), intent(inout)           :: self
    double precision                      , intent(in   )           :: x
    integer                               , intent(in   ), optional :: table
    integer                                                         :: i,tableActual
    double precision                                                :: a,b,c,d,dx
    
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table       
    ! Check for recall with same value as previous call.
    if (x /= self%xPrevious .or. tableActual /= self%tablePrevious) then
       ! Determine the location in the table.
       if      (x <  self%xv(          1)) then
          i=1
       else if (x >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=int((x-self%xv(1))*self%inverseDeltaX)+1
       end if      
       ! Compute polynomial coefficients.
       call Table_Linear_CSpline_1D_Coefficients(self,tableActual,x,i,a,b,c,d,dx)
       ! Interpolate in the table.
       self%tablePrevious=tableActual
       self%    yPrevious=a+dx*(b+dx*(c+dx*d))
    end if
    Table_Linear_CSpline_1D_Interpolate=self%yPrevious
    return
  end function Table_Linear_CSpline_1D_Interpolate
  
  double precision function Table_Linear_CSpline_1D_Interpolate_Gradient(self,x,table)
    !% Perform linear interpolation in a linear 1D table.
    implicit none
    class           (table1DLinearCSpline), intent(inout)           :: self
    double precision                      , intent(in   )           :: x
    integer                               , intent(in   ), optional :: table
    integer                                                         :: i,tableActual
    double precision                                                :: a,b,c,d,dx
        
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then
       ! Determine the location in the table.
       if      (x <  self%xv(          1)) then
          i=1
       else if (x >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=int((x-self%xv(1))*self%inverseDeltaX)+1
       end if
       ! Compute polynomial coefficients.
       call Table_Linear_CSpline_1D_Coefficients(self,tableActual,x,i,a,b,c,d,dx)       
       ! Interpolate in the table.
       self%dTablePrevious=tableActual
       self%    dyPrevious=b+dx*(2.0d0*c+dx*3.0d0*d)
    end if
    Table_Linear_CSpline_1D_Interpolate_Gradient=self%dyPrevious
    return
  end function Table_Linear_CSpline_1D_Interpolate_Gradient
  
  subroutine Table_Logarithmic_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount)
    !% Create a 1-D logarithmic table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class           (table1DLogarithmicCSpline), intent(inout)           :: self
    double precision                           , intent(in   )           :: xMinimum,xMaximum
    integer                                    , intent(in   )           :: xCount
    integer                                    , intent(in   ), optional :: tableCount
    
    ! Call the creator for linear tables with the logarithms of the input x range.
    call self%table1DLinearCSpline%create(log(xMinimum),log(xMaximum),xCount,tableCount)
    return
  end subroutine Table_Logarithmic_CSpline_1D_Create
  
  double precision function Table_Logarithmic_CSpline_1D_X(self,i)
    !% Return the {\tt i}$^{\rm th}$ $x$-value for a logarithmic 1D table.
    use Galacticus_Error
    implicit none
    class  (table1DLogarithmicCSpline), intent(inout) :: self
    integer                           , intent(in   ) :: i
    
    Table_Logarithmic_CSpline_1D_X=exp(self%table1DLinearCSpline%x(i))
    return
  end function Table_Logarithmic_CSpline_1D_X
  
  function Table_Logarithmic_CSpline_1D_Xs(self)
    !% Return the $x$-values for a 1D table.
    use Galacticus_Error
    implicit none
    class           (table1DLogarithmicCSpline), intent(in   )             :: self
    double precision                           , dimension(size(self%xv))  :: Table_Logarithmic_CSpline_1D_Xs

    Table_Logarithmic_CSpline_1D_Xs=exp(self%table1DLinearCSpline%xs())
    return
  end function Table_Logarithmic_CSpline_1D_Xs

  double precision function Table_Logarithmic_CSpline_1D_Interpolate(self,x,table)
    !% Perform linear interpolation in a logarithmic 1D table.
    implicit none
    class           (table1DLogarithmicCSpline), intent(inout)           :: self
    double precision                           , intent(in   )           :: x
    integer                                    , intent(in   ), optional :: table
    
    Table_Logarithmic_CSpline_1D_Interpolate=self%table1DLinearCSpline%interpolate(log(x),table)
    return
  end function Table_Logarithmic_CSpline_1D_Interpolate
  
  double precision function Table_Logarithmic_CSpline_1D_Interpolate_Gradient(self,x,table)
    !% Perform linear interpolation in a logarithmic 1D table.
    implicit none
    class           (table1DLogarithmicCSpline), intent(inout)           :: self
    double precision                           , intent(in   )           :: x
    integer                                    , intent(in   ), optional :: table
    
    Table_Logarithmic_CSpline_1D_Interpolate_Gradient=self%table1DLinearCSpline%interpolateGradient(log(x),table)/x
    return
  end function Table_Logarithmic_CSpline_1D_Interpolate_Gradient

end module Tables
