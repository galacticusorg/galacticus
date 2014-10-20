!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  public :: table,table1D,table1DGeneric,table1DLinearLinear,table1DLogarithmicLinear,table1DNonUniformLinearLogarithmic,table1DLinearCSpline,table1DLogarithmicCSpline

  !@ <enumeration>
  !@  <name>extrapolationType</name>
  !@  <description>Used to specify the type of extrapolation to use when interpolating in tables.</description>
  !@  <entry label="extrapolationTypeExtrapolate"/>
  !@  <entry label="extrapolationTypeFix"        />
  !@  <entry label="extrapolationTypeAbort"      />
  !@ </enumeration>
  integer, parameter, public :: extrapolationTypeExtrapolate=1
  integer, parameter, public :: extrapolationTypeFix        =2
  integer, parameter, public :: extrapolationTypeAbort      =3

  type, abstract :: table
     !% Basic table type.
   contains
     !@ <objectMethods>
     !@   <object>table</object>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <arguments></arguments>
     !@     <type>\void</type>
     !@     <description>Destroy the table.</description>
     !@   </objectMethod>
     !@ </objectMethods>
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
     integer                                       :: extrapolationType, xCount
     double precision, allocatable, dimension(:  ) :: xv
     double precision, allocatable, dimension(:,:) :: yv
   contains
     !@ <objectMethods>
     !@   <object>table1D</object>
     !@   <objectMethod>
     !@     <method>interpolate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ x,\intzero\ [table]</arguments>
     !@     <description>Interpolate to {\tt x} in the {\tt table}$^{\rm th}$ table.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>interpolateGradient</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ x,\intzero\ [table]</arguments>
     !@     <description>Interpolate the gradient to {\tt x} in the {\tt table}$^{\rm th}$ table.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reverse</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(table)\textgreater} reversedSelf,\intzero\ [table], \logicalzero\ [precise]</arguments>
     !@     <description>Reverse the table (i.e. swap $x$ and $y$ components) and return in {\tt reversedSelf}. If {\tt table} is specified then the {\tt table}$^{\rm th}$ table is used for the $y$-values, otherwise the first table is used. If the optional {\tt precise} argument is set to {\tt true} then the reversal must be precisely invertible---if this is not possible the method will abort.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isMonotonic</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\enum\ [directionDecreasing|directionIncreasing],\logicalzero\ [allowEqual],\intzero\ [table]</arguments>
     !@     <description>Return true if the table $y$-values are monotonic. Optionally, the direction of monotonicity can be specified via the {\tt direction} argument---by default either direction is allowed. By default consecutive equal values are considered non-monotonic. This behavior can be changed via the optional {\tt allowEqual} argument. If {\tt table} is specified then the {\tt table}$^{\rm th}$ table is used for the $y$-values, otherwise the first table is used.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>size</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return the size (i.e. number of $x$-values) in the table.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>x</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\intzero\ i</arguments>
     !@     <description>Return the {\tt i}$^{\rm th}$ $x$-value.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>y</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\intzero\ i,\intzero\ [table]</arguments>
     !@     <description>Return the {\tt i}$^{\rm th}$ $y$-value. If {\tt table} is specified then the {\tt table}$^{\rm th}$ table is used for the $y$-values, otherwise the first table is used.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>xs</method>
     !@     <type>\doubleone</type>
     !@     <arguments>\intzero\ i</arguments>
     !@     <description>Return an array of all $x$-values.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>ys</method>
     !@     <type>\doubleone</type>
     !@     <arguments>\intzero\ i,\intzero\ [table]</arguments>
     !@     <description>Return an array of all $y$-values. If {\tt table} is specified then the {\tt table}$^{\rm th}$ table is used for the $y$-values, otherwise the first table is used.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>xEffective</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ x</arguments>
     !@     <description>Return the effective value of $x$ to use in table interpolations.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>integrationWeights</method>
     !@     <type>\doubleone</type>
     !@     <arguments>\doublezero\ x0\argin, \doublezero\ x1\argin</arguments>
     !@     <description>Return the weights to be applied to the table to integrate (using the trapezium rule) between {\tt x0} and {\tt x1}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(Table1D_Interpolate ), deferred :: interpolate
     procedure(Table1D_Interpolate ), deferred :: interpolateGradient
     procedure                                 :: destroy            =>Table_1D_Destroy
     procedure                                 :: reverse            =>Table_1D_Reverse
     procedure                                 :: isMonotonic        =>Table1D_Is_Monotonic
     procedure                                 :: size               =>Table1D_Size
     procedure                                 :: x                  =>Table1D_X
     procedure                                 :: y                  =>Table1D_Y
     procedure                                 :: xs                 =>Table1D_Xs
     procedure                                 :: ys                 =>Table1D_Ys
     procedure                                 :: xEffective         =>Table1D_Find_Effective_X
     procedure                                 :: integrationWeights =>Table1D_Integration_Weights
  end type table1D

  interface
     double precision function Table1D_Interpolate(self,x,table)
       !% Interface to {\tt table} interpolator.
       import table1D
       implicit none
       class           (table1D), intent(inout)           :: self
       double precision         , intent(in   )           :: x
       integer                  , intent(in   ), optional :: table
     end function Table1D_Interpolate
  end interface

  type, extends(table1D) :: table1DGeneric
     !% Table type supporting generic one dimensional tables.
     type   (fgsl_interp      ) :: interpolator
     type   (fgsl_interp_accel) :: accelerator
     logical                    :: reset
   contains
     !@ <objectMethods>
     !@   <object>table1DGeneric</object>
     !@   <objectMethod>
     !@     <method>create</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ x,\intzero\ [tableCount]</arguments>
     !@     <description>Create the object with the specified {\tt x} values, and with {\tt tableCount} tables.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>populate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero|\doubleone\ y,\intzero\ [i],\intzero\ [table]</arguments>
     !@     <description>Populate the {\tt table}$^{\rm th}$ table with elements {\tt y}. If {\tt y} is a scalar, then the index, {\tt i}, of the element to set must also be specified.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: create                          =>Table_Generic_1D_Create
     procedure :: destroy                         =>Table_Generic_1D_Destroy
     procedure :: Table_Generic_1D_Populate
     procedure :: Table_Generic_1D_Populate_Single
     generic   :: populate                         => Table_Generic_1D_Populate            , &
          &                                           Table_Generic_1D_Populate_Single
     procedure :: interpolate        =>Table_Generic_1D_Interpolate
     procedure :: interpolateGradient=>Table_Generic_1D_Interpolate_Gradient
  end type table1DGeneric

  type, extends(table1D) :: table1DLinearLinear
     !% Table type supporting one dimensional table with linear spacing in $x$.
     double precision :: dxPrevious    , dyPrevious   , inverseDeltaX, xPrevious, &
          &              yPrevious
     integer          :: dTablePrevious, tablePrevious
   contains
     !@ <objectMethods>
     !@   <object>table1DLinearLinear</object>
     !@   <objectMethod>
     !@     <type>\void</type>
     !@     <method>create</method>
     !@     <arguments>\doublezero\ xMinimum,\doublezero\ xMaximum,\intzero xCount,\intzero [tableCount]</arguments>
     !@     <description>Create the object with $x$-values spanning the range {\tt xMinimum} to {\tt xMaximum} in {\tt xCount} steps, and with {\tt tableCount} tables.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>populate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero|\doubleone\ y,\intzero\ [i],\intzero\ [table]</arguments>
     !@     <description>Populate the {\tt table}$^{\rm th}$ table with elements {\tt y}. If {\tt y} is a scalar, then the index, {\tt i}, of the element to set must also be specified.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: create                          => Table_Linear_1D_Create
     procedure :: Table_Linear_1D_Populate
     procedure :: Table_Linear_1D_Populate_Single
     generic   :: populate                        => Table_Linear_1D_Populate            , Table_Linear_1D_Populate_Single
     procedure :: interpolate                     => Table_Linear_1D_Interpolate
     procedure :: interpolateGradient             => Table_Linear_1D_Interpolate_Gradient
  end type table1DLinearLinear

  type, extends(table1DLinearLinear) :: table1DLogarithmicLinear
     !% Table type supporting one dimensional table with logarithmic spacing in $x$.
     logical          :: previousSet
     double precision :: xLinearPrevious,xLogarithmicPrevious
   contains
     procedure :: create              => Table_Logarithmic_1D_Create
     procedure :: interpolate         => Table_Logarithmic_1D_Interpolate
     procedure :: interpolateGradient => Table_Logarithmic_1D_Interpolate_Gradient
     procedure :: x                   => Table_Logarithmic_1D_X
     procedure :: xs                  => Table_Logarithmic_1D_Xs
     procedure :: integrationWeights  => Table_Logarithmic_Integration_Weights
     procedure :: reverse             => Table_Logarithmic_1D_Reverse
  end type table1DLogarithmicLinear

  type, extends(table1DGeneric) :: table1DNonUniformLinearLogarithmic
     !% Table type supporting one dimensional table with non-uniform x-axis and logarithmic in $y$.
   contains
     procedure :: Table_Generic_1D_Populate        => Table_NonUniform_Linear_Logarithmic_1D_Populate
     procedure :: Table_Generic_1D_Populate_Single => Table_NonUniform_Linear_Logarithmic_1D_Populate_Single
     procedure :: interpolate                      => Table_NonUniform_Linear_Logarithmic_1D_Interpolate
     procedure :: interpolateGradient              => Table_NonUniform_Linear_Logarithmic_1D_Interpolate_Gradient
     procedure :: y                                => Table_NonUniform_Linear_Logarithmic_1D_Y
     procedure :: ys                               => Table_NonUniform_Linear_Logarithmic_1D_Ys
     procedure :: integrationWeights               => Table_NonUniform_Linear_Logarithmic_Integration_Weights
  end type table1DNonUniformLinearLogarithmic

  type, extends(table1D) :: table1DLinearCSpline
     !% Table type supporting one dimensional table with linear spacing in $x$ and cubic spline interpolation.
     double precision, allocatable, dimension(:,:) :: sv
     integer                                       :: dTablePrevious, iPrevious    , tablePrevious
     double precision                              :: deltaX        , inverseDeltaX
     double precision                              :: aPrevious     , bPrevious    , cPrevious    , dPrevious, &
          &                                           dxPrevious    , dyPrevious   , xPrevious    , yPrevious
   contains
     !@ <objectMethods>
     !@   <object>table1DLinearCSpline</object>
     !@   <objectMethod>
     !@     <method>create</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ xMinimum,\doublezero\ xMaximum,\intzero xCount,\intzero [tableCount]</arguments>
     !@     <description>Create the object with $x$-values spanning the range {\tt xMinimum} to {\tt xMaximum} in {\tt xCount} steps, and with {\tt tableCount} tables.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>populate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero|\doubleone\ y,\intzero\ [i],\intzero\ [table]</arguments>
     !@     <description>Populate the {\tt table}$^{\rm th}$ table with elements {\tt y}. If {\tt y} is a scalar, then the index, {\tt i}, of the element to set must also be specified.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: create                                 =>Table_Linear_CSpline_1D_Create
     procedure :: destroy                                =>Table_Linear_CSpline_1D_Destroy
     procedure :: Table_Linear_CSpline_1D_Populate
     procedure :: Table_Linear_CSpline_1D_Populate_Single
     generic   :: populate                                => Table_Linear_CSpline_1D_Populate            , &
          &                                                  Table_Linear_CSpline_1D_Populate_Single
     procedure :: interpolate        =>Table_Linear_CSpline_1D_Interpolate
     procedure :: interpolateGradient=>Table_Linear_CSpline_1D_Interpolate_Gradient
     procedure :: integrationWeights =>Table_Linear_CSpline_Integration_Weights
  end type table1DLinearCSpline

  type, extends(table1DLinearCSpline) :: table1DLogarithmicCSpline
     !% Table type supporting one dimensional table with logarithmic spacing in $x$ and cubic spline interpolation.
     logical          :: previousSet
     double precision :: xLinearPrevious,xLogarithmicPrevious
     double precision :: xMinimum,xMaximum
   contains
     procedure :: create             =>Table_Logarithmic_CSpline_1D_Create
     procedure :: interpolate        =>Table_Logarithmic_CSpline_1D_Interpolate
     procedure :: interpolateGradient=>Table_Logarithmic_CSpline_1D_Interpolate_Gradient
     procedure :: x                  =>Table_Logarithmic_CSpline_1D_X
     procedure :: xs                 =>Table_Logarithmic_CSpline_1D_Xs
  end type table1DLogarithmicCSpline

  abstract interface
     double precision function integrandTemplate(x)
       double precision, intent(in   ) :: x
     end function integrandTemplate
  end interface

contains

  subroutine Table_1D_Destroy(self)
    !% Destroy a 1-D table.
    use Memory_Management
    implicit none
    class(table1D), intent(inout) :: self

    if (allocated(self%xv)) call Dealloc_Array(self%xv)
    if (allocated(self%yv)) call Dealloc_Array(self%yv)
    return
  end subroutine Table_1D_Destroy

  double precision function Table1D_X(self,i)
    !% Return the {\tt i}$^{\rm th}$ $x$-value for a 1D table.
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
    implicit none
    class(table1D), intent(in   ) :: self
    double precision         , dimension(size(self%xv))  :: Table1D_Xs

    Table1D_Xs=self%xv
    return
  end function Table1D_Xs

  double precision function Table1D_Y(self,i,table)
    !% Return the {\tt i}$^{\rm th}$ $y$-value for a 1D table.
    implicit none
    class  (table1D), intent(in   )           :: self
    integer         , intent(in   )           :: i
    integer         , intent(in   ), optional :: table
    integer                                   :: ii   , tableActual

    tableActual=1
    if (present(table)) tableActual=table
    ii=i
    if (ii < 0) ii=ii+size(self%xv)+1
    Table1D_Y=self%yv(ii,tableActual)
    return
  end function Table1D_Y

  function Table1D_Ys(self)
    !% Return the $y$-values for a 1D table.
    implicit none
    class(table1D), intent(in   ) :: self
    double precision         , dimension(size(self%yv,dim=1),size(self%yv,dim=2)) :: Table1D_Ys

    Table1D_Ys=self%yv
    return
  end function Table1D_Ys

  subroutine Table_1D_Reverse(self,reversedSelf,table,precise)
    !% Reverse a 1D table (i.e. swap $x$ and $y$ components). Optionally allows specification of
    !% which $y$ table to swap with.
    use Array_Utilities
    use Galacticus_Error
    implicit none
    class  (table1D)             , intent(in   )           :: self
    class  (table1D), allocatable, intent(inout)           :: reversedSelf
    integer                      , intent(in   ), optional :: table
    logical                      , intent(in   ), optional :: precise
    integer                                                :: i           , tableActual

    if (present(precise).and.precise) call Galacticus_Error_Report('Table_1D_Reverse','table cannot be precisely reversed')
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
       ! Copy the extrapolation option from the original table.
       reversedSelf%extrapolationType=self%extrapolationType
    end select
    return
  end subroutine Table_1D_Reverse

  logical function Table1D_Is_Monotonic(self,direction,allowEqual,table)
    !% Return true if a 1D table is monotonic. Optionally allows specification of the direction,
    !% and whether or not equal elements are allowed for monotonicity.
    use Array_Utilities
    implicit none
    class  (table1D), intent(in   )           :: self
    integer         , intent(in   ), optional :: direction  , table
    logical         , intent(in   ), optional :: allowEqual
    integer                                   :: tableActual

    tableActual=1
    if (present(table)) tableActual=table
    Table1D_Is_Monotonic=Array_Is_Monotonic(self%yv(:,tableActual),direction,allowEqual)
    return
  end function Table1D_Is_Monotonic

  integer function Table1D_Size(self)
    !% Return the size of a 1D table.
    implicit none
    class(table1D), intent(in   ) :: self

    Table1D_Size=self%xCount
    return
  end function Table1D_Size

  function Table1D_Integration_Weights(self,x0,x1,integrand)
    !% Returns a set of weights for trapezoidal integration on the table between limits {\tt x0} and {\tt x1}.
    use Galacticus_Error
    implicit none
    class           (table1D          ), intent(inout)                               :: self
    double precision                   , intent(in   )                               :: x0, x1
    procedure       (integrandTemplate), intent(in   )           , pointer, optional :: integrand
    double precision                   , dimension(size(self%xv))                    :: Table1D_Integration_Weights
    double precision                                                                 :: weight, lx0, lx1
    integer                                                                          :: i

    if (x1 < x0           ) call Galacticus_Error_Report('Table1D_Integration_Weights','inverted limits'         )
    if (present(integrand)) call Galacticus_Error_Report('Table1D_Integration_Weights','integrands not supported')
    Table1D_Integration_Weights=0.0d0
    do i=2,size(self%xv)
       if (self%xv(i) <= x1 .and. self%xv(i-1) >= x0) then
          weight=self%xv(i)-self%xv(i-1)
       else if ((self%xv(i-1) < x1 .and. self%xv(i) > x1) .or. (self%xv(i) > x0 .and. self%xv(i-1) < x0)) then
          lx0=max(self%xv(i-1),x0)
          lx1=min(self%xv(i  ),x1)
          weight=lx1-lx0
       else
          weight=0.0d0
       end if
       Table1D_Integration_Weights(i-1:i)=Table1D_Integration_Weights(i-1:i)+0.5d0*weight
    end do
    return
  end function Table1D_Integration_Weights
  
  subroutine Table_Generic_1D_Create(self,x,tableCount,extrapolationType)
    !% Create a 1-D generic table.
    use Memory_Management
    implicit none
    class           (table1DGeneric)              , intent(inout)           :: self
    double precision                , dimension(:), intent(in   )           :: x
    integer                                       , intent(in   ), optional :: extrapolationType, tableCount
    integer                                                                 :: tableCountActual

    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the x-range.
    self%xCount=size(x)
    call Alloc_Array(self%xv,[size(x)                 ])
    call Alloc_Array(self%yv,[size(x),tableCountActual])
    self%xv   =x
    self%reset=.true.
    ! Set extrapolation type.
    if (present(extrapolationType)) then
       self%extrapolationType=extrapolationType
    else
       self%extrapolationType=extrapolationTypeExtrapolate
    end if
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
    class           (table1DGeneric)              , intent(inout)           :: self
    double precision                , dimension(:), intent(in   )           :: y
    integer                                       , intent(in   ), optional :: table
    integer                                                                 :: tableActual

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
    integer                                                   :: tableActual, extrapolationType

    tableActual=1
    if (present(table)) tableActual=table
    select case (self%extrapolationType)
    case (extrapolationTypeAbort)
       extrapolationType=extrapolationTypeNone
    case (extrapolationTypeExtrapolate)
       extrapolationType=extrapolationTypeLinear
    case (extrapolationTypeFix)
       extrapolationType=extrapolationTypeFixed
    end select
    Table_Generic_1D_Interpolate=Interpolate(size(self%xv),self%xv,self%yv(:,tableActual),self%interpolator,self%accelerator,self%xEffective(x),extrapolationType=extrapolationType,reset=self%reset)
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
    Table_Generic_1D_Interpolate_Gradient=Interpolate_Derivative(size(self%xv),self%xv,self%yv(:,tableActual),self%interpolator,self%accelerator,self%xEffective(x),reset=self%reset)
    return
  end function Table_Generic_1D_Interpolate_Gradient

  subroutine Table_Linear_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !% Create a 1-D linear table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class           (table1DLinearLinear), intent(inout)           :: self
    double precision                     , intent(in   )           :: xMaximum         , xMinimum
    integer                              , intent(in   )           :: xCount
    integer                              , intent(in   ), optional :: extrapolationType, tableCount
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
    self%tablePrevious=-1
    self%xPrevious    =-1.0d0
    ! Set extrapolation type.
    if (present(extrapolationType)) then
       self%extrapolationType=extrapolationType
    else
       self%extrapolationType=extrapolationTypeExtrapolate
    end if
    return
  end subroutine Table_Linear_1D_Create

  subroutine Table_Linear_1D_Populate(self,y,table)
    !% Populate a 1-D linear table.
    use Galacticus_Error
    implicit none
    class           (table1DLinearLinear)              , intent(inout)           :: self
    double precision                     , dimension(:), intent(in   )           :: y
    integer                                            , intent(in   ), optional :: table
    integer                                                                      :: tableActual

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
    integer                                                        :: i    , tableActual
    double precision                                               :: h    , xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    xEffective=self%xEffective(x)
    if (xEffective /= self%xPrevious .or. tableActual /= self%tablePrevious) then
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=max(min(int((xEffective-self%xv(1))*self%inverseDeltaX)+1,self%xCount-1),1)
       end if
       ! Compute the interpolating factor.
       h=(xEffective-self%xv(i))*self%inverseDeltaX
       ! Interpolate in the table.
       self%xPrevious    =xEffective
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
    integer                                                        :: i         , tableActual
    double precision                                               :: xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    xEffective=self%xEffective(x)
    if (xEffective /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=int((xEffective-self%xv(1))*self%inverseDeltaX)+1
       end if
       ! Interpolate in the table.
       self%dxPrevious    =xEffective
       self%dTablePrevious=tableActual
       self%    dyPrevious=(self%yv(i+1,tableActual)-self%yv(i,tableActual))*self%inverseDeltaX
    end if
    Table_Linear_1D_Interpolate_Gradient=self%dyPrevious
    return
  end function Table_Linear_1D_Interpolate_Gradient

  subroutine Table_Logarithmic_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !% Create a 1-D logarithmic table.
    implicit none
    class           (table1DLogarithmicLinear), intent(inout)           :: self
    double precision                          , intent(in   )           :: xMaximum         , xMinimum
    integer                                   , intent(in   )           :: xCount
    integer                                   , intent(in   ), optional :: extrapolationType, tableCount

    self%previousSet=.false.
    ! Call the creator for linear tables with the logarithms of the input x range.
    call self%table1DLinearLinear%create(log(xMinimum),log(xMaximum),xCount,tableCount,extrapolationType)
    return
  end subroutine Table_Logarithmic_1D_Create

  double precision function Table_Logarithmic_1D_X(self,i)
    !% Return the {\tt i}$^{\rm th}$ $x$-value for a logarithmic 1D table.
    implicit none
    class  (table1DLogarithmicLinear), intent(inout) :: self
    integer                          , intent(in   ) :: i

    Table_Logarithmic_1D_X=exp(self%table1DLinearLinear%x(i))
    return
  end function Table_Logarithmic_1D_X

  function Table_Logarithmic_1D_Xs(self)
    !% Return the $x$-values for a 1D table.
    implicit none
    class(table1DLogarithmicLinear), intent(in   ) :: self
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

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_1D_Interpolate=self%table1DLinearLinear%interpolate(self%xLogarithmicPrevious,table)
    return
  end function Table_Logarithmic_1D_Interpolate

  double precision function Table_Logarithmic_1D_Interpolate_Gradient(self,x,table)
    !% Perform linear interpolation in a logarithmic 1D table.
    implicit none
    class           (table1DLogarithmicLinear), intent(inout)           :: self
    double precision                          , intent(in   )           :: x
    integer                                   , intent(in   ), optional :: table
   
    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_1D_Interpolate_Gradient=self%table1DLinearLinear%interpolateGradient(self%xLogarithmicPrevious,table)/self%xEffective(x)
    return
  end function Table_Logarithmic_1D_Interpolate_Gradient

  function Table_Logarithmic_Integration_Weights(self,x0,x1,integrand)
    !% Returns a set of weights for trapezoidal integration on the table between limits {\tt x0} and {\tt x1}.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Integration
    use Galacticus_Error
    implicit none
    class           (table1DLogarithmicLinear ), intent(inout)                               :: self
    double precision                           , intent(in   )                               :: x0, x1
    procedure       (integrandTemplate        ), intent(in   )           , pointer, optional :: integrand
    double precision                           , dimension(size(self%xv))                    :: Table_Logarithmic_Integration_Weights
    double precision                           , parameter                                   :: logTolerance=1.0d-12
    double precision                                                                         :: gradientTerm, lx0, lx1, factor0, factor1
    integer                                                                                  :: i
    type            (fgsl_function             )                                             :: integrandFunction
    type            (fgsl_integration_workspace)                                             :: integrationWorkspace
    type            (c_ptr                     )                                             :: parameterPointer
    logical                                                                                  :: integrationReset
 
    if (x1 < x0) call Galacticus_Error_Report('Table_Logarithmic_Integration_Weights','inverted limits')
    Table_Logarithmic_Integration_Weights=0.0d0    
    do i=2,size(self%xv)
       ! Evaluate integration range for this interval of the table.
       if (self%xv(i) <= log(x1) .and. self%xv(i-1) >= log(x0)) then
          lx0=self%xv(i-1)
          lx1=self%xv(i  )
       else if ((self%xv(i-1) < log(x1) .and. self%xv(i) > log(x1)) .or. (self%xv(i) > log(x0) .and. self%xv(i-1) < log(x0))) then
          lx0=max(self%xv(i-1),log(x0))
          lx1=min(self%xv(i  ),log(x1))
       else
          cycle
       end if
       ! Proceed only for non-zero ranges. Add some tolerance to avoid attempting to evaluate for tiny ranges which arise from
       ! numerical imprecision.
       if (lx1 > lx0+logTolerance) then
          if (present(integrand)) then
             ! An integrand is given, numerically integrate the relevant terms over the integrand.
             integrationReset=.true.
             factor0=Integrate(                                       &
                  &            lx0                                  , &
                  &            lx1                                  , &
                  &            factor0Integrand                     , &
                  &            parameterPointer                     , &
                  &            integrandFunction                    , &
                  &            integrationWorkspace                 , &
                  &            toleranceRelative   =1.0d-4          , &
                  &            reset               =integrationReset  &
                  &           )
             call Integrate_Done(integrandFunction,integrationWorkspace)
             integrationReset=.true.
             factor1=Integrate(                                       &
                  &            lx0                                  , &
                  &            lx1                                  , &
                  &            factor1Integrand                     , &
                  &            parameterPointer                     , &
                  &            integrandFunction                    , &
                  &            integrationWorkspace                 , &
                  &            toleranceRelative   =1.0d-4          , &
                  &            reset               =integrationReset  &
                  &           )
             call Integrate_Done(integrandFunction,integrationWorkspace)
             Table_Logarithmic_Integration_Weights        (i-1) &
                  & =Table_Logarithmic_Integration_Weights(i-1) &
                  & +factor0                                    &
                  & -factor1
             Table_Logarithmic_Integration_Weights        (i  ) &
                  & =Table_Logarithmic_Integration_Weights(i  ) &
                  & +factor1  
          else
             ! No additional integrand, use analytic solution.
             gradientTerm=(exp(lx1)*((lx1-lx0)-1.0d0)+exp(lx0))/(lx1-lx0)
             Table_Logarithmic_Integration_Weights(i-1)=Table_Logarithmic_Integration_Weights(i-1)+max((exp(lx1)-exp(lx0))-gradientTerm,0.0d0)
             Table_Logarithmic_Integration_Weights(i  )=Table_Logarithmic_Integration_Weights(i  )+max(                   +gradientTerm,0.0d0)
          end if
       end if
    end do
    return
    
  contains
    
    function factor0Integrand(logx,parameterPointer) bind(c)
      !% Integrand used to evaluate integration weights over logarithmically spaced tables
      implicit none
      real(c_double)        :: factor0Integrand
      real(c_double), value :: logx
      type(c_ptr),    value :: parameterPointer
      real(c_double)        :: x
      
      x=exp(logx)
      factor0Integrand=x*integrand(x)
      return
    end function factor0Integrand
    
    function factor1Integrand(logx,parameterPointer) bind(c)
      !% Integrand used to evaluate integration weights over logarithmically spaced tables
      implicit none
      real(c_double)        :: factor1Integrand
      real(c_double), value :: logx
      type(c_ptr),    value :: parameterPointer
      real(c_double)        :: x
      
      x=exp(logx)
      factor1Integrand=x*integrand(x)*(logx-self%xv(i-1))/(self%xv(i)-self%xv(i-1))
      return
    end function factor1Integrand
    
  end function Table_Logarithmic_Integration_Weights

  subroutine Table_Logarithmic_1D_Reverse(self,reversedSelf,table,precise)
    !% Reverse a 1D logarithmic-linear table (i.e. swap $x$ and $y$ components). Optionally allows specification of
    !% which $y$ table to swap with.
    use Array_Utilities
    use Galacticus_Error
    implicit none
    class  (table1DLogarithmicLinear)             , intent(in   )           :: self
    class  (table1D                 ), allocatable, intent(inout)           :: reversedSelf
    integer                                       , intent(in   ), optional :: table
    logical                                       , intent(in   ), optional :: precise
    integer                                                                 :: i           , tableActual

    tableActual=1
    if (present(table)) tableActual=table
    if (.not.Array_Is_Monotonic(self%yv(:,tableActual))) call Galacticus_Error_Report('Table_Logarithmic_1D_Reverse','reversed table would not be monotonic')
    if (allocated(reversedSelf)) deallocate(reversedSelf)
    allocate(table1DNonUniformLinearLogarithmic :: reversedSelf)
    select type (reversedSelf)
    type is (table1DNonUniformLinearLogarithmic)
       call reversedSelf%create(self%yv(:,tableActual),size(self%yv,dim=2),extrapolationType=self%extrapolationType)
       reversedSelf%yv               =self%yv
       reversedSelf%yv(:,tableActual)=self%xv
        ! If the table is monotonically decreasing when reversed, we will also need to switch the order.
       if (Array_Is_Monotonic(reversedSelf%xv,direction=directionDecreasing)) then
          reversedSelf%xv(:)=Array_Reverse(reversedSelf%xv(:))
          do i=1,size(reversedSelf%yv,dim=2)
             reversedSelf%yv(:,i)=Array_Reverse(reversedSelf%yv(:,i))
          end do
       end if
    end select
    return
  end subroutine Table_Logarithmic_1D_Reverse

  subroutine Table_Linear_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !% Create a 1-D linear table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class           (table1DLinearCSpline), intent(inout)           :: self
    double precision                      , intent(in   )           :: xMaximum         , xMinimum
    integer                               , intent(in   )           :: xCount
    integer                               , intent(in   ), optional :: extrapolationType, tableCount
    integer                                                         :: tableCountActual

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
    self%tablePrevious=-1
    self%iPrevious    =-1
    ! Set extrapolation type.
    if (present(extrapolationType)) then
       self%extrapolationType=extrapolationType
    else
       self%extrapolationType=extrapolationTypeExtrapolate
    end if
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
    class           (table1DLinearCSpline)              , intent(inout)           :: self
    double precision                      , dimension(:), intent(in   )           :: y
    integer                                             , intent(in   ), optional :: table
    logical                                             , intent(in   ), optional :: computeSpline
    integer                                                                       :: tableActual
    logical                                                                       :: computeSplineActual

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
    double precision                      , allocatable  , dimension(:) :: b    , u, v
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
    integer                               , intent(in   ) :: i   , table
    double precision                      , intent(  out) :: a   , b    , c, d, dx

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
    integer                                                         :: i         , tableActual
    double precision                                                :: a         , b          , c, d, dx, &
         &                                                             xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    xEffective=self%xEffective(x)
    if (xEffective /= self%xPrevious .or. tableActual /= self%tablePrevious) then
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=int((xEffective-self%xv(1))*self%inverseDeltaX)+1
       end if
       ! Compute polynomial coefficients.
       call Table_Linear_CSpline_1D_Coefficients(self,tableActual,xEffective,i,a,b,c,d,dx)
       ! Interpolate in the table.
       self%xPrevious    =xEffective
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
    integer                                                         :: i         , tableActual
    double precision                                                :: a         , b          , c, d, dx, &
         &                                                             xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    xEffective=self%xEffective(x)
    if (xEffective /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=int((xEffective-self%xv(1))*self%inverseDeltaX)+1
       end if
       ! Compute polynomial coefficients.
       call Table_Linear_CSpline_1D_Coefficients(self,tableActual,xEffective,i,a,b,c,d,dx)
       ! Interpolate in the table.
       self%dxPrevious    =xEffective
       self%dTablePrevious=tableActual
       self%    dyPrevious=b+dx*(2.0d0*c+dx*3.0d0*d)
    end if
    Table_Linear_CSpline_1D_Interpolate_Gradient=self%dyPrevious
    return
  end function Table_Linear_CSpline_1D_Interpolate_Gradient

  subroutine Table_Logarithmic_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !% Create a 1-D logarithmic table.
    implicit none
    class           (table1DLogarithmicCSpline), intent(inout)           :: self
    double precision                           , intent(in   )           :: xMaximum         , xMinimum
    integer                                    , intent(in   )           :: xCount
    integer                                    , intent(in   ), optional :: extrapolationType, tableCount

    self%previousSet         =.false.
    self%xLinearPrevious     =-1.0d0
    self%xLogarithmicPrevious=-1.0d0
    ! Call the creator for linear tables with the logarithms of the input x range.
    call self%table1DLinearCSpline%create(log(xMinimum),log(xMaximum),xCount,tableCount,extrapolationType)
    ! Store the minimum and maximum x-values for rapid look-up.
    self%xMinimum=exp(self%xv(     1))
    self%xMaximum=exp(self%xv(xCount))
    return
  end subroutine Table_Logarithmic_CSpline_1D_Create

  double precision function Table_Logarithmic_CSpline_1D_X(self,i)
    !% Return the {\tt i}$^{\rm th}$ $x$-value for a logarithmic 1D table.
    implicit none
    class  (table1DLogarithmicCSpline), intent(inout) :: self
    integer                           , intent(in   ) :: i

    ! Check for end-points, and return stored values if possible.
    if      (i ==  1                      ) then
       Table_Logarithmic_CSpline_1D_X=self%xMinimum
    else if (i == -1 .or. i == self%xCount) then
       Table_Logarithmic_CSpline_1D_X=self%xMaximum
    else
       ! No stored value is available - simply look up the required value.
       Table_Logarithmic_CSpline_1D_X=exp(self%table1DLinearCSpline%x(i))
    end if
    return
  end function Table_Logarithmic_CSpline_1D_X

  function Table_Logarithmic_CSpline_1D_Xs(self)
    !% Return the $x$-values for a 1D table.
    implicit none
    class(table1DLogarithmicCSpline), intent(in   ) :: self
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

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_CSpline_1D_Interpolate=self%table1DLinearCSpline%interpolate(self%xLogarithmicPrevious,table)
    return
  end function Table_Logarithmic_CSpline_1D_Interpolate

  double precision function Table_Logarithmic_CSpline_1D_Interpolate_Gradient(self,x,table)
    !% Perform linear interpolation in a logarithmic 1D table.
    implicit none
    class           (table1DLogarithmicCSpline), intent(inout)           :: self
    double precision                           , intent(in   )           :: x
    integer                                    , intent(in   ), optional :: table

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_CSpline_1D_Interpolate_Gradient=self%table1DLinearCSpline%interpolateGradient(self%xLogarithmicPrevious,table)/self%xEffective(x)
    return
  end function Table_Logarithmic_CSpline_1D_Interpolate_Gradient

  function Table_Linear_CSpline_Integration_Weights(self,x0,x1,integrand)
    !% Returns a set of weights for trapezoidal integration on the table between limits {\tt x0} and {\tt x1}.
    use Galacticus_Error
    implicit none
    class           (table1DLinearCSpline), intent(inout)                               :: self
    double precision                      , intent(in   )                               :: x0, x1
    procedure       (integrandTemplate   ), intent(in   )           , pointer, optional :: integrand
    double precision                      , dimension(size(self%xv))                    :: Table_Linear_CSpline_Integration_Weights

    call Galacticus_Error_Report('Table_Linear_CSpline_Integration_Weights','integration weights not supported')
    return
  end function Table_Linear_CSpline_Integration_Weights
  
  double precision function Table1D_Find_Effective_X(self,x)
    !% Return the effective value of $x$ to use in table interpolations.
    use Galacticus_Error
    implicit none
    class           (table1D), intent(inout) :: self
    double precision         , intent(in   ) :: x

    if (x < self%x(1) .or. x > self%x(-1)) then
       select case (self%extrapolationType)
       case (extrapolationTypeExtrapolate)
          Table1D_Find_Effective_X=x
       case (extrapolationTypeFix        )
          if (x < self%x(1)) then
             Table1D_Find_Effective_X=self%x( 1)
          else
             Table1D_Find_Effective_X=self%x(-1)
          end if
       case default
          call Galacticus_Error_Report('Table_1D_Find_Effective_X','x is out of range')
       end select
    else
       Table1D_Find_Effective_X=x
    end if
    return
  end function Table1D_Find_Effective_X

  subroutine Table_NonUniform_Linear_Logarithmic_1D_Populate(self,y,table)
    !% Populate a 1-D linear-logarihtmic table.
    use Galacticus_Error
    implicit none
    class           (table1DNonUniformLinearLogarithmic)              , intent(inout)           :: self
    double precision                                    , dimension(:), intent(in   )           :: y
    integer                                                          , intent(in   ), optional :: table

    call self%table1DGeneric%populate(log(y),table)
    return
  end subroutine Table_NonUniform_Linear_Logarithmic_1D_Populate

  subroutine Table_NonUniform_Linear_Logarithmic_1D_Populate_Single(self,y,i,table)
    !% Populate a single element of a 1-D linear table.
    use Galacticus_Error
    implicit none
    class           (table1DNonUniformLinearLogarithmic), intent(inout)           :: self
    double precision                                    , intent(in   )           :: y
    integer                                             , intent(in   )           :: i
    integer                                             , intent(in   ), optional :: table

    call self%table1DGeneric%populate(log(y),i,table)
    return
  end subroutine Table_NonUniform_Linear_Logarithmic_1D_Populate_Single

  double precision function Table_NonUniform_Linear_Logarithmic_1D_Interpolate(self,x,table)
    !% Perform linear interpolation in a linear-logarithmic 1D table.
    implicit none
    class           (table1DNonUniformLinearLogarithmic), intent(inout)           :: self
    double precision                                    , intent(in   )           :: x
    integer                                             , intent(in   ), optional :: table

    Table_NonUniform_Linear_Logarithmic_1D_Interpolate=exp(self%table1DGeneric%interpolate(x,table))
    return
  end function Table_NonUniform_Linear_Logarithmic_1D_Interpolate

  double precision function Table_NonUniform_Linear_Logarithmic_1D_Interpolate_Gradient(self,x,table)
    !% Perform linear interpolation in a linear-logarithmic 1D table.
    implicit none
    class           (table1DNonUniformLinearLogarithmic), intent(inout)           :: self
    double precision                                    , intent(in   )           :: x
    integer                                             , intent(in   ), optional :: table
   
    Table_NonUniform_Linear_Logarithmic_1D_Interpolate_Gradient=self%interpolate(x,table)*self%table1DGeneric%interpolateGradient(x,table)
    return
  end function Table_NonUniform_Linear_Logarithmic_1D_Interpolate_Gradient

  function Table_NonUniform_Linear_Logarithmic_Integration_Weights(self,x0,x1,integrand)
    !% Returns a set of weights for integration on a linear-logarithmic table between limits {\tt x0} and {\tt x1}.
    use Galacticus_Error
    implicit none
    class           (table1DNonUniformLinearLogarithmic ), intent(inout)                               :: self
    double precision                                     , intent(in   )                               :: x0, x1
    procedure       (integrandTemplate                  ), intent(in   )           , pointer, optional :: integrand
    double precision                                     , dimension(size(self%xv))                    :: Table_NonUniform_Linear_Logarithmic_Integration_Weights

    call Galacticus_Error_Report('Table_NonUniform_Linear_Logarithmic_Integration_Weights','integrand is not linear in y')
    return    
  end function Table_NonUniform_Linear_Logarithmic_Integration_Weights
  
  double precision function Table_NonUniform_Linear_Logarithmic_1D_Y(self,i,table)
    !% Return the {\tt i}$^{\rm th}$ $y$-value for a 1D table.
    implicit none
    class  (table1DNonUniformLinearLogarithmic), intent(in   )           :: self
    integer                                    , intent(in   )           :: i
    integer                                    , intent(in   ), optional :: table
    integer                                                              :: ii   , tableActual

    tableActual=1
    if (present(table)) tableActual=table
    ii=i
    if (ii < 0) ii=ii+size(self%xv)+1
    Table_NonUniform_Linear_Logarithmic_1D_Y=exp(self%yv(ii,tableActual))
    return
  end function Table_NonUniform_Linear_Logarithmic_1D_Y

  function Table_NonUniform_Linear_Logarithmic_1D_Ys(self)
    !% Return the $y$-values for a 1D table.
    implicit none
    class           (table1DNonUniformLinearLogarithmic), intent(in   )                                      :: self
    double precision                                    , dimension(size(self%yv,dim=1),size(self%yv,dim=2)) :: Table_NonUniform_Linear_Logarithmic_1D_Ys

    Table_NonUniform_Linear_Logarithmic_1D_Ys=exp(self%yv)
    return
  end function Table_NonUniform_Linear_Logarithmic_1D_Ys

end module Tables
