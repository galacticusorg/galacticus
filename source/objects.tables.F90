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
Contains a module which defines a {\normalfont \ttfamily table} class with optimized interpolation operators.
!!}

module Tables
  !!{
  Defines a {\normalfont \ttfamily table} class with optimized interpolation operators.
  !!}
  use :: Numerical_Interpolation, only : interpolator
  use :: Table_Labels           , only : enumerationExtrapolationTypeType
  private
  public :: table                          , table1D                          , table1DGeneric                    , &
       &    table1DLinearLinear            , table1DLogarithmicLinear         , table1DNonUniformLinearLogarithmic, &
       &    table1DLinearCSpline           , table1DLogarithmicCSpline        , table2DLogLogLin                  , &
       &    table1DLinearMonotoneCSpline   , table1DLogarithmicMonotoneCSpline, table2DLinLinLin                  , &
       &    tablesIntegrationWeightFunction, table1DMonotoneCSpline

  !![
  <enumeration>
   <name>tableType</name>
   <description>Enumeration of table types.</description>
   <entry label="linearLinear1D"      />
   <entry label="logarithmicLinear1D" />
  </enumeration>
  !!]
  
  !![
  <stateStorable class="table">
   <table1DGeneric>
    <methodCall method="interpolatorReinitialize"/>
   </table1DGeneric>
   <table2DLinLinLin>
    <methodCall method="interpolatorReinitialize"/>
   </table2DLinLinLin>
  </stateStorable>
  !!]

  !![
  <deepCopyActions class="table">
   <table1DGeneric>
    <methodCall method="interpolatorReinitialize"/>
   </table1DGeneric>
   <table2DLinLinLin>
    <methodCall method="interpolatorReinitialize"/>
   </table2DLinLinLin>
  </deepCopyActions>
  !!]
    
  type, abstract :: table
     !!{
     Basic table type.
     !!}
   contains
     !![
     <methods>
      <method method="destroy" description="Destroy the table."/>
     </methods>
     !!]
     procedure(Table_Destroy), deferred :: destroy
  end type table

  interface
     subroutine Table_Destroy(self)
       !!{
       Interface to {\normalfont \ttfamily table} destructor.
       !!}
       import table
       implicit none
       class(table), intent(inout) :: self
     end subroutine Table_Destroy
  end interface

  type, abstract, extends(table) :: table1D
     !!{
     Basic table type.
     !!}
     integer                                                                         :: xCount
     type            (enumerationExtrapolationTypeType)             , dimension(2  ) :: extrapolationType
     double precision                                  , allocatable, dimension(:  ) :: xv
     double precision                                  , allocatable, dimension(:,:) :: yv
   contains
     !![
     <methods>
       <method description="Interpolate to {\normalfont \ttfamily x} in the {\normalfont \ttfamily table}$^\mathrm{th}$ table." method="interpolate" />
       <method description="Interpolate the gradient to {\normalfont \ttfamily x} in the {\normalfont \ttfamily table}$^\mathrm{th}$ table." method="interpolateGradient" />
       <method description="Reverse the table (i.e. swap $x$ and $y$ components) and return in {\normalfont \ttfamily reversedSelf}. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^\mathrm{th}$ table is used for the $y$-values, otherwise the first table is used. If the optional {\normalfont \ttfamily precise} argument is set to {\normalfont \ttfamily true} then the reversal must be precisely invertible---if this is not possible the method will abort." method="reverse" />
       <method description="Return true if the table $y$-values are monotonic. Optionally, the direction of monotonicity can be specified via the {\normalfont \ttfamily direction} argument---by default either direction is allowed. By default consecutive equal values are considered non-monotonic. This behavior can be changed via the optional {\normalfont \ttfamily allowEqual} argument. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^\mathrm{th}$ table is used for the $y$-values, otherwise the first table is used." method="isMonotonic" />
       <method description="Return the size (i.e. number of $x$-values) in the table." method="size" />
       <method description="Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $x$-value." method="x" />
       <method description="Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $y$-value. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^\mathrm{th}$ table is used for the $y$-values, otherwise the first table is used." method="y" />
       <method description="Return an array of all $x$-values." method="xs" />
       <method description="Return an array of all $y$-values. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^\mathrm{th}$ table is used for the $y$-values, otherwise the first table is used." method="ys" />
       <method description="Return the effective value of $x$ to use in table interpolations." method="xEffective"/>
       <method description="Return the weights to be applied to the table to integrate (using the trapezium rule) between {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}." method="integrationWeights" />
     </methods>
     !!]
     procedure(Table1D_Interpolate ), deferred :: interpolate
     procedure(Table1D_Interpolate ), deferred :: interpolateGradient
     procedure                                 :: destroy             => Table_1D_Destroy
     procedure                                 :: reverse             => Table_1D_Reverse
     procedure                                 :: isMonotonic         => Table1D_Is_Monotonic
     procedure                                 :: size                => Table1D_Size
     procedure                                 :: x                   => Table1D_X
     procedure                                 :: y                   => Table1D_Y
     procedure                                 :: xs                  => Table1D_Xs
     procedure                                 :: ys                  => Table1D_Ys
     procedure                                 :: xEffective          => Table1D_Find_Effective_X
     procedure                                 :: integrationWeights  => Table1D_Integration_Weights
     procedure                                 ::                        Table1D_Assignment
     generic                                   :: assignment(=)       => Table1D_Assignment
  end type table1D

  interface
     double precision function Table1D_Interpolate(self,x,table,status)
       !!{
       Interface to {\normalfont \ttfamily table} interpolator.
       !!}
       import table1D
       implicit none
       class           (table1D), intent(inout)           :: self
       double precision         , intent(in   )           :: x
       integer                  , intent(in   ), optional :: table
       integer                  , intent(  out), optional :: status
     end function Table1D_Interpolate
  end interface

  type, extends(table1D) :: table1DGeneric
     !!{
     Table type supporting generic one dimensional tables.
     !!}
     type   (interpolator), allocatable, dimension(:) :: interpolator_
     logical              , allocatable, dimension(:) :: interpolatorInitialized
     integer                                          :: interpolationType
   contains
     !![
     <methods>
       <method description="Create the object with the specified {\normalfont \ttfamily x} values, and with {\normalfont \ttfamily tableCount} tables." method="create" />
       <method description="Populate the {\normalfont \ttfamily table}$^\mathrm{th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified." method="populate" />
       <method description="Reinitialize the interpolator." method="interpolatorReinitialize" />
       <method description="Initialize the interpolator." method="interpolatorInitialize" />
     </methods>
     !!]
     procedure :: create                   => Table_Generic_1D_Create
     procedure :: destroy                  => Table_Generic_1D_Destroy
     procedure :: populate_                => Table_Generic_1D_Populate
     procedure :: populateSingle_          => Table_Generic_1D_Populate_Single
     generic   :: populate                 => populate_                                 , &
          &                                   populateSingle_
     procedure :: interpolate              => Table_Generic_1D_Interpolate
     procedure :: interpolateGradient      => Table_Generic_1D_Interpolate_Gradient
     procedure :: interpolatorInitialize   => Table_Generic_1D_Interpolator_Initialize
     procedure :: interpolatorReinitialize => Table_Generic_1D_Interpolator_Reinitialize
  end type table1DGeneric

  type, extends(table1D) :: table1DLinearLinear
     !!{
     Table type supporting one dimensional table with linear spacing in $x$.
     !!}
     double precision :: dxPrevious    , dyPrevious   , inverseDeltaX, xPrevious, &
          &              yPrevious
     integer          :: dTablePrevious, tablePrevious
   contains
     !![
     <methods>
       <method description="Create the object with $x$-values spanning the range {\normalfont \ttfamily xMinimum} to {\normalfont \ttfamily xMaximum} in {\normalfont \ttfamily xCount} steps, and with {\normalfont \ttfamily tableCount} tables." method="create" />
       <method description="Populate the {\normalfont \ttfamily table}$^\mathrm{th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified." method="populate" />
     </methods>
     !!]
     procedure :: create              => Table_Linear_1D_Create
     procedure ::                        Table_Linear_1D_Populate
     procedure ::                        Table_Linear_1D_Populate_Single
     generic   :: populate            => Table_Linear_1D_Populate            , Table_Linear_1D_Populate_Single
     procedure :: interpolate         => Table_Linear_1D_Interpolate
     procedure :: interpolateGradient => Table_Linear_1D_Interpolate_Gradient
  end type table1DLinearLinear

  type, extends(table1DLinearLinear) :: table1DLogarithmicLinear
     !!{
     Table type supporting one dimensional table with logarithmic spacing in $x$.
     !!}
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
     !!{
     Table type supporting one dimensional table with non-uniform x-axis and logarithmic in $y$.
     !!}
   contains
     procedure :: populate_           => Table_NonUniform_Linear_Logarithmic_1D_Populate
     procedure :: populateSingle_     => Table_NonUniform_Linear_Logarithmic_1D_Populate_Single
     procedure :: interpolate         => Table_NonUniform_Linear_Logarithmic_1D_Interpolate
     procedure :: interpolateGradient => Table_NonUniform_Linear_Logarithmic_1D_Interpolate_Gradient
     procedure :: y                   => Table_NonUniform_Linear_Logarithmic_1D_Y
     procedure :: ys                  => Table_NonUniform_Linear_Logarithmic_1D_Ys
     procedure :: integrationWeights  => Table_NonUniform_Linear_Logarithmic_Integration_Weights
  end type table1DNonUniformLinearLogarithmic

  type, extends(table1D) :: table1DLinearCSpline
     !!{
     Table type supporting one dimensional table with linear spacing in $x$ and cubic spline interpolation.
     !!}
     double precision, allocatable, dimension(:,:) :: sv            , av           , bv           , cv       , &
          &                                           dv
     integer                                       :: dTablePrevious, iPrevious    , tablePrevious
     double precision                              :: deltaX        , inverseDeltaX
     double precision                              :: aPrevious     , bPrevious    , cPrevious    , dPrevious, &
          &                                           dxPrevious    , dyPrevious   , xPrevious    , yPrevious
   contains
     !![
     <methods>
       <method description="Create the object with $x$-values spanning the range {\normalfont \ttfamily xMinimum} to {\normalfont \ttfamily xMaximum} in {\normalfont \ttfamily xCount} steps, and with {\normalfont \ttfamily tableCount} tables." method="create" />
       <method description="Populate the {\normalfont \ttfamily table}$^\mathrm{th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified." method="populate" />
     </methods>
     !!]
     procedure :: create              => Table_Linear_CSpline_1D_Create
     procedure :: destroy             => Table_Linear_CSpline_1D_Destroy
     procedure :: populateArray       => Table_Linear_CSpline_1D_Populate
     procedure :: populateScalar      => Table_Linear_CSpline_1D_Populate_Single
     procedure :: interpolate         => Table_Linear_CSpline_1D_Interpolate
     procedure :: interpolateGradient => Table_Linear_CSpline_1D_Interpolate_Gradient
     procedure :: integrationWeights  => Table_Linear_CSpline_Integration_Weights
     generic   :: populate            => populateArray, populateScalar
  end type table1DLinearCSpline

  type, extends(table1DLinearCSpline) :: table1DLogarithmicCSpline
     !!{
     Table type supporting one dimensional table with logarithmic spacing in $x$ and cubic spline interpolation.
     !!}
     logical          :: previousSet
     double precision :: xLinearPrevious, xLogarithmicPrevious
     double precision :: xMinimum       , xMaximum
   contains
     procedure :: create              => Table_Logarithmic_CSpline_1D_Create
     procedure :: interpolate         => Table_Logarithmic_CSpline_1D_Interpolate
     procedure :: interpolateGradient => Table_Logarithmic_CSpline_1D_Interpolate_Gradient
     procedure :: x                   => Table_Logarithmic_CSpline_1D_X
     procedure :: xs                  => Table_Logarithmic_CSpline_1D_Xs
  end type table1DLogarithmicCSpline

  type, extends(table1DLinearCSpline) :: table1DLinearMonotoneCSpline
     !!{
     Table type supporting one dimensional table with linear spacing in $x$ and monotonic cubic spline interpolation.
     !!}
   contains
     procedure :: create              => Table_Linear_Monotone_CSpline_1D_Create
     procedure :: destroy             => Table_Linear_Monotone_CSpline_1D_Destroy
     procedure :: populateArray       => Table_Linear_Monotone_CSpline_1D_Populate
     procedure :: populateScalar      => Table_Linear_Monotone_CSpline_1D_Populate_Single
     procedure :: interpolate         => Table_Linear_Monotone_CSpline_1D_Interpolate
     procedure :: interpolateGradient => Table_Linear_Monotone_CSpline_1D_Interpolate_Gradient
     procedure :: integrationWeights  => Table_Linear_Monotone_CSpline_Integration_Weights
  end type table1DLinearMonotoneCSpline

  type, extends(table1DLinearMonotoneCSpline) :: table1DLogarithmicMonotoneCSpline
     !!{
     Table type supporting one dimensional table with logarithmic spacing in $x$ and monotonic cubic spline interpolation.
     !!}
     logical          :: previousSet
     double precision :: xLinearPrevious, xLogarithmicPrevious
     double precision :: xMinimum       , xMaximum
   contains
     procedure :: create              => Table_Logarithmic_Monotone_CSpline_1D_Create
     procedure :: interpolate         => Table_Logarithmic_Monotone_CSpline_1D_Interpolate
     procedure :: interpolateGradient => Table_Logarithmic_Monotone_CSpline_1D_Interpolate_Gradient
     procedure :: x                   => Table_Logarithmic_Monotone_CSpline_1D_X
     procedure :: xs                  => Table_Logarithmic_Monotone_CSpline_1D_Xs
  end type table1DLogarithmicMonotoneCSpline

  type, extends(table1D) :: table1DMonotoneCSpline
     !!{
     Table type supporting one dimensional table with arbitrary spacing in $x$ and monotonic cubic spline interpolation.
     !!}
     double precision, allocatable, dimension(:,:) :: av            , bv           , cv
     integer                                       :: dTablePrevious, iPrevious    , tablePrevious
     double precision                              :: aPrevious     , bPrevious    , cPrevious    , dPrevious, &
          &                                           dxPrevious    , dyPrevious   , xPrevious    , yPrevious
   contains
     !![
     <methods>
       <method description="Create the object with {\normalfont \ttfamily xCount} points, and with {\normalfont \ttfamily tableCount} tables." method="create" />
       <method description="Populate the {\normalfont \ttfamily table}$^\mathrm{th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified." method="populate" />
     </methods>
     !!]
     procedure :: create              => Table_Monotone_CSpline_1D_Create
     procedure :: destroy             => Table_Monotone_CSpline_1D_Destroy
     procedure :: populateArray       => Table_Monotone_CSpline_1D_Populate
     procedure :: populateScalar      => Table_Monotone_CSpline_1D_Populate_Single
     procedure :: interpolate         => Table_Monotone_CSpline_1D_Interpolate
     procedure :: interpolateGradient => Table_Monotone_CSpline_1D_Interpolate_Gradient
     procedure :: integrationWeights  => Table_Monotone_CSpline_Integration_Weights
     generic   :: populate            => populateArray, populateScalar
  end type table1DMonotoneCSpline
  
  abstract interface
     double precision function tablesIntegrationWeightFunction(x)
       double precision, intent(in   ) :: x
     end function tablesIntegrationWeightFunction
  end interface

  type, extends(table) :: table2DLinLinLin
     !!{
     Table type supporting generic two dimensional tables.
     !!}
     integer                                                       :: xCount       , yCount
     double precision              , allocatable, dimension(:    ) :: xv           , yv
     double precision              , allocatable, dimension(:,:,:) :: zv
     type            (interpolator)                                :: interpolatorX, interpolatorY
   contains
     !![
     <methods>
       <method description="Create the object with the specified {\normalfont \ttfamily x} and {\normalfont \ttfamily y} values, and with {\normalfont \ttfamily tableCount} tables." method="create" />
       <method description="Populate the {\normalfont \ttfamily table}$^\mathrm{th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the indices, {\normalfont \ttfamily i}, {\normalfont \ttfamily j}, of the element to set must also be specified." method="populate" />
       <method description="Interpolate to {\normalfont \ttfamily x}, {\normalfont \ttfamily y} in the {\normalfont \ttfamily table}$^\mathrm{th}$ table." method="interpolate" />
       <method description="Return an array of all {\normalfont \ttfamily x} values." method="xs" />
       <method description="Return an array of all {\normalfont \ttfamily y} values." method="ys" />
       <method description="Return an array of all {\normalfont \ttfamily z} values." method="zs" />
       <method description="Reinitialize the interpolator." method="interpolatorReinitialize" />
     </methods>
     !!]
     procedure :: create                           => Table_2D_LinLinLin_Create
     procedure :: destroy                          => Table_2D_LinLinLin_Destroy
     procedure :: Table_2D_LinLinLin_Populate
     procedure :: Table_2D_LinLinLin_Populate_Single
     generic   :: populate                         => Table_2D_LinLinLin_Populate       , &
          &                                           Table_2D_LinLinLin_Populate_Single
     procedure :: interpolate                      => Table_2D_LinLinLin_Interpolate
     procedure :: xs                               => Table_2D_LinLinLin_Xs
     procedure :: ys                               => Table_2D_LinLinLin_Ys
     procedure :: zs                               => Table_2D_LinLinLin_Zs
     procedure :: interpolatorReinitialize         => Table_2D_LinLinLin_Interpolator_Reinitialize
  end type table2DLinLinLin

  type, extends(table) :: table2DLogLogLin
     !!{
     Two-dimensional table type with logarithmic spacing in x and y dimensions, and linear interpolation in z.
     !!}
     type            (enumerationExtrapolationTypeType)                                :: extrapolationTypeX  , extrapolationTypeY
     integer                                                                           :: xCount              , yCount              , &
          &                                                                               i                   , j                   , &
          &                                                                               tablePrevious       , dimPrevious
     double precision                                                                  :: xLinearPrevious     , yLinearPrevious     , &
          &                                                                               xLogarithmicPrevious, yLogarithmicPrevious, &
          &                                                                               hx                  , hy                  , &
          &                                                                               inverseDeltaX       , inverseDeltaY       , &
          &                                                                               zPrevious           , dzPrevious
     logical                                                                           :: xPreviousSet        , yPreviousSet
     double precision                                  , allocatable, dimension(:    ) :: xv                  , yv
     double precision                                  , allocatable, dimension(:,:,:) :: zv
   contains
     !![
     <methods>
       <method description="Compute and store interpolation factors to {\normalfont \ttfamily (x,y)}." method="interpolationFactors" />
       <method description="Interpolate to {\normalfont \ttfamily x} in the {\normalfont \ttfamily table}$^\mathrm{th}$ table." method="interpolate" />
       <method description="Interpolate the gradient to {\normalfont \ttfamily x} in the {\normalfont \ttfamily table}$^\mathrm{th}$ table." method="interpolateGradient" />
       <method description="Return the size (i.e. number of $x$ or $y$-values) in the table of the given dimension." method="size" />
       <method description="Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $x$-value." method="x" />
       <method description="Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $y$-value. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^\mathrm{th}$ table is used for the $y$-values, otherwise the first table is used." method="y" />
       <method description="Return the {\normalfont \ttfamily (i,j)}$^\mathrm{th}$ $z$-value. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^\mathrm{th}$ table is used for the $z$-values, otherwise the first table is used." method="z" />
       <method description="Return an array of all $x$-values." method="xs" />
       <method description="Return an array of all $y$-values. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^\mathrm{th}$ table is used for the $y$-values, otherwise the first table is used." method="ys" />
       <method description="Return an array of all $z$-values. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^\mathrm{th}$ table is used for the $z$-values, otherwise the first table is used." method="zs" />
       <method description="Return true if the table is initialized (this means the table is created, it may not yet have been populated)." method="isInitialized" />
       <method description="Populate the {\normalfont \ttfamily table}$^\mathrm{th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified." method="populate" />
       <method description="Create the object with $x$-values spanning the range {\normalfont \ttfamily xMinimum} to {\normalfont \ttfamily xMaximum} in {\normalfont \ttfamily xCount} steps, and with {\normalfont \ttfamily tableCount} tables." method="create" />
     </methods>
     !!]
     procedure :: create                            => Table_2DLogLogLin_Create
     procedure :: Table_2DLogLogLin_Populate
     procedure :: Table_2DLogLogLin_Populate_Single
     generic   :: populate                          => Table_2DLogLogLin_Populate             , &
          &                                            Table_2DLogLogLin_Populate_Single
     procedure :: interpolationFactors              => Table_2DLogLogLin_Interpolation_Factors
     procedure :: interpolate                       => Table_2DLogLogLin_Interpolate
     procedure :: interpolateGradient               => Table_2DLogLogLin_Interpolate_Gradient
     procedure :: destroy                           => Table_2DLogLogLin_Destroy
     procedure :: size                              => Table_2DLogLogLin_Size
     procedure :: x                                 => Table_2DLogLogLin_X
     procedure :: y                                 => Table_2DLogLogLin_Y
     procedure :: z                                 => Table_2DLogLogLin_z
     procedure :: xs                                => Table_2DLogLogLin_Xs
     procedure :: ys                                => Table_2DLogLogLin_Ys
     procedure :: zs                                => Table_2DLogLogLin_Zs
     procedure :: isInitialized                     => Table_2DLogLogLin_Is_Initialized
  end type table2DLogLogLin

contains

  subroutine Table_1D_Destroy(self)
    !!{
    Destroy a 1-D table.
    !!}
    implicit none
    class(table1D), intent(inout) :: self

    if (allocated(self%xv)) deallocate(self%xv)
    if (allocated(self%yv)) deallocate(self%yv)
    return
  end subroutine Table_1D_Destroy

  double precision function Table1D_X(self,i)
    !!{
    Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $x$-value for a 1D table.
    !!}
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
    !!{
    Return the $x$-values for a 1D table.
    !!}
    implicit none
    class(table1D), intent(in   ) :: self
    double precision         , dimension(size(self%xv))  :: Table1D_Xs

    Table1D_Xs=self%xv
    return
  end function Table1D_Xs

  double precision function Table1D_Y(self,i,table)
    !!{
    Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $y$-value for a 1D table.
    !!}
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
    !!{
    Return the $y$-values for a 1D table.
    !!}
    implicit none
    class(table1D), intent(in   ) :: self
    double precision         , dimension(size(self%yv,dim=1),size(self%yv,dim=2)) :: Table1D_Ys

    Table1D_Ys=self%yv
    return
  end function Table1D_Ys

  subroutine Table_1D_Reverse(self,reversedSelf,table,precise)
    !!{
    Reverse a 1D table (i.e. swap $x$ and $y$ components). Optionally allows specification of
    which $y$ table to swap with.
    !!}
    use :: Array_Utilities        , only : Array_Is_Monotonic, Array_Reverse, directionDecreasing
    use :: Error                  , only : Error_Report
    use :: Numerical_Interpolation, only : GSL_Interp_Linear
    implicit none
    class  (table1D)             , intent(in   )           :: self
    class  (table1D), allocatable, intent(inout)           :: reversedSelf
    integer                      , intent(in   ), optional :: table
    logical                      , intent(in   ), optional :: precise
    integer                                                :: i           , tableActual

    if (present(precise).and.precise) call Error_Report('table cannot be precisely reversed'//{introspection:location})
    tableActual=1
    if (present(table)) tableActual=table
    if (.not.Array_Is_Monotonic(self%yv(:,tableActual))) call Error_Report('reversed table would not be monotonic'//{introspection:location})
    if (allocated(reversedSelf)) then
       call reversedSelf%destroy()
       deallocate(reversedSelf)
    end if
    allocate(table1DGeneric :: reversedSelf)
    select type (reversedSelf)
    type is (table1DGeneric)
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
       ! Set linear interpolation.
       reversedSelf%interpolationType=GSL_Interp_Linear
       ! Build the interpolator.
       allocate(reversedSelf%interpolator_          (1))
       allocate(reversedSelf%interpolatorInitialized(1))
       reversedSelf%interpolator_          (1)=interpolator(reversedSelf%xv,interpolationType=GSL_Interp_Linear,extrapolationType=self%extrapolationType(1))
       reversedSelf%interpolatorInitialized(1)=.true.
    end select
    return
  end subroutine Table_1D_Reverse

  logical function Table1D_Is_Monotonic(self,direction,allowEqual,table)
    !!{
    Return true if a 1D table is monotonic. Optionally allows specification of the direction,
    and whether or not equal elements are allowed for monotonicity.
    !!}
    use :: Array_Utilities, only : Array_Is_Monotonic
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
    !!{
    Return the size of a 1D table.
    !!}
    implicit none
    class(table1D), intent(in   ) :: self

    Table1D_Size=self%xCount
    return
  end function Table1D_Size

  function Table1D_Integration_Weights(self,x0,x1,integrand)
    !!{
    Returns a set of weights for trapezoidal integration on the table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1D                        ), intent(inout)                               :: self
    double precision                                 , intent(in   )                               :: x0, x1
    procedure       (tablesIntegrationWeightFunction), intent(in   )           , pointer, optional :: integrand
    double precision                                 , dimension(size(self%xv))                    :: Table1D_Integration_Weights
    double precision                                                                               :: weight, lx0, lx1
    integer                                                                                        :: i

    if (x1 < x0           ) call Error_Report('inverted limits'//{introspection:location})
    if (present(integrand)) call Error_Report('integrands not supported'//{introspection:location})
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
  
  !![
  <workaround type="gfortran" PR="121537" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=121537">
    <description>gfortran misses a defined-assignment of a component.</description>
  </workaround>
  !!]
  subroutine Table1D_Assignment(to,from)
    !!{
    Assignment operator for the {\normalfont \ttfamily table1D} class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (table1D), intent(  out) :: to
    class  (table1D), intent(in   ) :: from
    integer                         :: i

    to%xCount           =from%xCount
    to%extrapolationType=from%extrapolationType
    if (allocated(from%xv)) then
       allocate(to%xv(size(from%xv)))
       to%xv=from%xv
    end if
    if (allocated(from%yv)) then
       allocate(to%yv(size(from%yv,dim=1),size(from%yv,dim=2)))
       to%yv=from%yv
    end if
    select type (from)
    class is (table1DGeneric)
       select type (to)
       class is (table1DGeneric)
          if (allocated(from%interpolator_)) then
             allocate(to%interpolator_(size(from%interpolator_)))
             !![
	     <workaround type="gfortran" PR="46897" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=46897">
	       <seeAlso type="gfortran" PR="57696" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=57696"/>
	       <description>Type-bound defined assignment not done because multiple part array references would occur in intermediate expressions.</description>
             !!]
             do i=1,size(from%interpolator_)
                to%interpolator_(i)=from%interpolator_(i)
             end do
             !![
	     </workaround>
             !!]
          end if
          if (allocated(from%interpolatorInitialized)) then
             allocate(to%interpolatorInitialized(size(from%interpolatorInitialized)))
             to%interpolatorInitialized=from%interpolatorInitialized
          end if
          to%xCount           =from%xCount
          to%extrapolationType=from%extrapolationType
          to%interpolationType=from%interpolationType
       class default
          call Error_Report('mismatched types in assignment'//{introspection:location})
       end select
    class is (table1DLinearLinear)
       select type (to)
       class is (table1DLinearLinear)
          to%dxPrevious    =from%dxPrevious
          to%dyPrevious    =from%dyPrevious
          to%inverseDeltaX =from%inverseDeltaX
          to%xPrevious     =from%xPrevious
          to%yPrevious     =from%yPrevious
          to%dTablePrevious=from%dTablePrevious
          to%tablePrevious =from%tablePrevious
      class default
          call Error_Report('mismatched types in assignment'//{introspection:location})
       end select
    class is (table1DLogarithmicLinear)
       select type (to)
       class is (table1DLogarithmicLinear)
          to%dxPrevious          =from%dxPrevious
          to%dyPrevious          =from%dyPrevious
          to%inverseDeltaX       =from%inverseDeltaX
          to%xPrevious           =from%xPrevious
          to%yPrevious           =from%yPrevious
          to%dTablePrevious      =from%dTablePrevious
          to%tablePrevious       =from%tablePrevious
          to%previousSet         =from%previousSet
          to%xLinearPrevious     =from%xLinearPrevious
          to%xLogarithmicPrevious=from%xLogarithmicPrevious
       class default
          call Error_Report('mismatched types in assignment'//{introspection:location})
       end select
    class is (table1DLinearCSpline)
       select type (to)
       class is (table1DLinearCSpline)
          if (allocated(to%sv)) deallocate(to%sv)
          if (allocated(to%av)) deallocate(to%av)
          if (allocated(to%bv)) deallocate(to%bv)
          if (allocated(to%cv)) deallocate(to%cv)
          if (allocated(to%dv)) deallocate(to%dv)
          if (allocated(from%sv)) allocate(to%sv,source=from%sv)
          if (allocated(from%av)) allocate(to%av,source=from%av)
          if (allocated(from%bv)) allocate(to%bv,source=from%bv)
          if (allocated(from%cv)) allocate(to%cv,source=from%cv)
          if (allocated(from%dv)) allocate(to%dv,source=from%dv)
          to%dTablePrevious=from%dTablePrevious
          to%iPrevious     =from%iPrevious
          to%tablePrevious =from%tablePrevious
          to%deltaX        =from%deltaX
          to%inverseDeltaX =from%inverseDeltaX
          to%aPrevious     =from%aPrevious
          to%bPrevious     =from%bPrevious
          to%cPrevious     =from%cPrevious
          to%dPrevious     =from%dPrevious
          to%dxPrevious    =from%dxPrevious
          to%dyPrevious    =from%dyPrevious
          to%xPrevious     =from%xPrevious
          to%yPrevious     =from%yPrevious
       class default
          call Error_Report('mismatched types in assignment'//{introspection:location})
       end select
    class is (table1DLogarithmicCSpline)
       select type (to)
       class is (table1DLogarithmicCSpline)
          if (allocated(to%sv)) deallocate(to%sv)
          if (allocated(to%av)) deallocate(to%av)
          if (allocated(to%bv)) deallocate(to%bv)
          if (allocated(to%cv)) deallocate(to%cv)
          if (allocated(to%dv)) deallocate(to%dv)
          if (allocated(from%sv)) allocate(to%sv,source=from%sv)
          if (allocated(from%av)) allocate(to%av,source=from%av)
          if (allocated(from%bv)) allocate(to%bv,source=from%bv)
          if (allocated(from%cv)) allocate(to%cv,source=from%cv)
          if (allocated(from%dv)) allocate(to%dv,source=from%dv)
          to%dTablePrevious      =from%dTablePrevious
          to%iPrevious           =from%iPrevious
          to%tablePrevious       =from%tablePrevious
          to%deltaX              =from%deltaX
          to%inverseDeltaX       =from%inverseDeltaX
          to%aPrevious           =from%aPrevious
          to%bPrevious           =from%bPrevious
          to%cPrevious           =from%cPrevious
          to%dPrevious           =from%dPrevious
          to%dxPrevious          =from%dxPrevious
          to%dyPrevious          =from%dyPrevious
          to%xPrevious           =from%xPrevious
          to%yPrevious           =from%yPrevious
          to%previousSet         =from%previousSet
          to%xLinearPrevious     =from%xLinearPrevious
          to%xLogarithmicPrevious=from%xLogarithmicPrevious
          to%xMinimum            =from%xMinimum
          to%xMaximum            =from%xMaximum
       class default
          call Error_Report('mismatched types in assignment'//{introspection:location})
       end select
    class is (table1DLogarithmicMonotoneCSpline)
       select type (to)
       class is (table1DLogarithmicMonotoneCSpline)
          if (allocated(to%sv)) deallocate(to%sv)
          if (allocated(to%av)) deallocate(to%av)
          if (allocated(to%bv)) deallocate(to%bv)
          if (allocated(to%cv)) deallocate(to%cv)
          if (allocated(to%dv)) deallocate(to%dv)
          if (allocated(from%sv)) allocate(to%sv,source=from%sv)
          if (allocated(from%av)) allocate(to%av,source=from%av)
          if (allocated(from%bv)) allocate(to%bv,source=from%bv)
          if (allocated(from%cv)) allocate(to%cv,source=from%cv)
          if (allocated(from%dv)) allocate(to%dv,source=from%dv)
          to%dTablePrevious      =from%dTablePrevious
          to%iPrevious           =from%iPrevious
          to%tablePrevious       =from%tablePrevious
          to%deltaX              =from%deltaX
          to%inverseDeltaX       =from%inverseDeltaX
          to%aPrevious           =from%aPrevious
          to%bPrevious           =from%bPrevious
          to%cPrevious           =from%cPrevious
          to%dPrevious           =from%dPrevious
          to%dxPrevious          =from%dxPrevious
          to%dyPrevious          =from%dyPrevious
          to%xPrevious           =from%xPrevious
          to%yPrevious           =from%yPrevious
          to%previousSet         =from%previousSet
          to%xLinearPrevious     =from%xLinearPrevious
          to%xLogarithmicPrevious=from%xLogarithmicPrevious
          to%xMinimum            =from%xMinimum
          to%xMaximum            =from%xMaximum
       class default
          call Error_Report('mismatched types in assignment'//{introspection:location})
       end select
    class is (table1DMonotoneCSpline)
       select type (to)
       class is (table1DMonotoneCSpline)
          if (allocated(to%av)) deallocate(to%av)
          if (allocated(to%bv)) deallocate(to%bv)
          if (allocated(to%cv)) deallocate(to%cv)
          if (allocated(from%av)) allocate(to%av,source=from%av)
          if (allocated(from%bv)) allocate(to%bv,source=from%bv)
          if (allocated(from%cv)) allocate(to%cv,source=from%cv)
          to%dTablePrevious      =from%dTablePrevious
          to%iPrevious           =from%iPrevious
          to%tablePrevious       =from%tablePrevious
          to%aPrevious           =from%aPrevious
          to%bPrevious           =from%bPrevious
          to%cPrevious           =from%cPrevious
          to%dPrevious           =from%dPrevious
          to%dxPrevious          =from%dxPrevious
          to%dyPrevious          =from%dyPrevious
          to%xPrevious           =from%xPrevious
          to%yPrevious           =from%yPrevious
       class default
          call Error_Report('mismatched types in assignment'//{introspection:location})
       end select
    class default
       call Error_Report('unsupported type in assignment'//{introspection:location})
    end select
    return
  end subroutine Table1D_Assignment
    
  subroutine Table_Generic_1D_Create(self,x,tableCount,extrapolationType,interpolationType)
    !!{
    Create a 1-D generic table.
    !!}
    use :: Error                  , only : Error_Report
    use :: Numerical_Interpolation, only : GSL_Interp_Linear
    use :: Table_Labels           , only : extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (table1DGeneric                  )              , intent(inout)           :: self
    double precision                                  , dimension(:), intent(in   )           :: x
    integer                                                         , intent(in   ), optional :: interpolationType, tableCount
    type            (enumerationExtrapolationTypeType), dimension(2), intent(in   ), optional :: extrapolationType
    !![
    <optionalArgument name="tableCount"        defaultsTo="1"                           />
    <optionalArgument name="extrapolationType" defaultsTo="extrapolationTypeExtrapolate"/>
    <optionalArgument name="interpolationType" defaultsTo="GSL_Interp_Linear"           />
    !!]
    
    ! Allocate arrays and construct the x-range.
    self%xCount=size(x)
    allocate(self%xv(size(x)            ))
    allocate(self%yv(size(x),tableCount_))
    self%xv               =x
    self%interpolationType=interpolationType_
    self%extrapolationType=extrapolationType_
    ! Validate extrapolation type.
    if (any(extrapolationType_ == extrapolationTypeZero)) call Error_Report('zero extrapolation is not supported'//{introspection:location})
    ! Allocate interpolators, and build if possible.
    if (interpolationType_ == GSL_Interp_Linear) then
       allocate(self%interpolator_          (1))
       allocate(self%interpolatorInitialized(1))
       self%interpolator_          (1)=interpolator(self%xv,extrapolationType=self%extrapolationType,interpolationType=self%interpolationType)
       self%interpolatorInitialized(1)=.true.
    else
       allocate(self%interpolator_          (tableCount_))
       allocate(self%interpolatorInitialized(tableCount_))
       self%interpolatorInitialized   =.false.
    end if
    return
  end subroutine Table_Generic_1D_Create

  subroutine Table_Generic_1D_Destroy(self)
    !!{
    Destroy a generic 1-D table.
    !!}
    implicit none
    class(table1DGeneric), intent(inout) :: self

    if (allocated(self%interpolator_          )) deallocate(self%interpolator_          )
    if (allocated(self%interpolatorInitialized)) deallocate(self%interpolatorInitialized)
    call Table_1D_Destroy(self)
    return
  end subroutine Table_Generic_1D_Destroy

  subroutine Table_Generic_1D_Populate(self,y,table)
    !!{
    Populate a 1-D generic table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DGeneric)              , intent(inout)           :: self
    double precision                , dimension(:), intent(in   )           :: y
    integer                                       , intent(in   ), optional :: table
    !![
    <optionalArgument name="table" defaultsTo="1"/>
    !!]

    ! Validate the input.
    if (.not.allocated(self%yv)       ) call Error_Report("create the table before populating it"//{introspection:location})
    if (size(self%yv,dim=1) /= size(y)) call Error_Report("provided y array is of wrong size"    //{introspection:location})
    ! Store the y values.
    self%yv(:,table_)=y
    return
  end subroutine Table_Generic_1D_Populate

  subroutine Table_Generic_1D_Populate_Single(self,y,i,table)
    !!{
    Populate a single element of a 1-D generic table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DGeneric), intent(inout)           :: self
    double precision                , intent(in   )           :: y
    integer                         , intent(in   )           :: i
    integer                         , intent(in   ), optional :: table
    !![
    <optionalArgument name="table" defaultsTo="1"/>
    !!]

    ! Validate the input.
    if (.not.allocated(self%yv)           ) call Error_Report("create the table before populating it"//{introspection:location})
    if (i < 1 .or. i > size(self%yv,dim=1)) call Error_Report("provided i value is out of bounds"    //{introspection:location})
    ! Store the y values.
    self%yv(i,table_)=y
    return
  end subroutine Table_Generic_1D_Populate_Single

  subroutine Table_Generic_1D_Interpolator_Initialize(self,table)
    !!{
    Initialize an interpolator.
    !!}
    use :: Numerical_Interpolation, only : GSL_Interp_Linear
    implicit none
    class  (table1DGeneric), intent(inout) :: self
    integer                , intent(in   ) :: table

    if (self%interpolationType == GSL_Interp_Linear) return
    if (.not.self%interpolatorInitialized(table)) then
       self%interpolator_          (table)=interpolator(self%xv,self%yv(:,table),extrapolationType=self%extrapolationType,interpolationType=self%interpolationType)
       self%interpolatorInitialized(table)=.true.
    end if
    return
  end subroutine Table_Generic_1D_Interpolator_Initialize
  
  double precision function Table_Generic_1D_Interpolate(self,x,table,status)
    !!{
    Perform generic interpolation in a generic 1D table.
    !!}
    use :: Numerical_Interpolation, only : GSL_Interp_Linear
    implicit none
    class           (table1DGeneric), intent(inout)           :: self
    double precision                , intent(in   )           :: x
    integer                         , intent(in   ), optional :: table
    integer                         , intent(  out), optional :: status
    integer                                                   :: interpolator_
    !![
    <optionalArgument name="table" defaultsTo="1"/>
    !!]

    call self%interpolatorInitialize(table_)
    if (self%interpolationType == GSL_Interp_Linear) then
       interpolator_=1
    else
       interpolator_=table_
    end if
    Table_Generic_1D_Interpolate=self%interpolator_(interpolator_)%interpolate(self%xEffective(x,status),self%yv(:,table_))
    return
  end function Table_Generic_1D_Interpolate
  
  double precision function Table_Generic_1D_Interpolate_Gradient(self,x,table,status)
    !!{
    Perform generic interpolation in a generic 1D table.
    !!}
    use :: Numerical_Interpolation, only : GSL_Interp_Linear
    implicit none
    class           (table1DGeneric), intent(inout)           :: self
    double precision                , intent(in   )           :: x
    integer                         , intent(in   ), optional :: table
    integer                         , intent(  out), optional :: status
    integer                                                   :: interpolator_
    !![
    <optionalArgument name="table" defaultsTo="1"/>
    !!]
    
    call self%interpolatorInitialize(table_)
    if (self%interpolationType == GSL_Interp_Linear) then
       interpolator_=1
    else
       interpolator_=table_
    end if
    Table_Generic_1D_Interpolate_Gradient=self%interpolator_(interpolator_)%derivative(self%xEffective(x,status),self%yv(:,table_))
    return
  end function Table_Generic_1D_Interpolate_Gradient

  subroutine Table_Generic_1D_Interpolator_Reinitialize(self)
    !!{
    Reinitialize the interpolator.
    !!}
    implicit none
    class  (table1DGeneric), intent(inout) :: self
    integer                                :: i

    if (.not.allocated(self%interpolator_)) return
    do i=1,size(self%interpolator_)
       call self%interpolator_(i)%GSLReallocate()
    end do
    return
  end subroutine Table_Generic_1D_Interpolator_Reinitialize

  subroutine Table_Linear_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !!{
    Create a 1-D linear table.
    !!}
    use :: Error            , only : Error_Report
    use :: Numerical_Ranges , only : Make_Range                  , rangeTypeLinear
    use :: Table_Labels     , only : extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (table1DLinearLinear             ), intent(inout)                         :: self
    double precision                                  , intent(in   )                         :: xMaximum         , xMinimum
    integer                                           , intent(in   )                         :: xCount
    integer                                           , intent(in   ), optional               :: tableCount
    type            (enumerationExtrapolationTypeType), intent(in   ), optional, dimension(2) :: extrapolationType
    integer                                                                                   :: tableCountActual

    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the x-range.
    self%xCount=xCount
    if (allocated(self%xv)) deallocate(self%xv)
    if (allocated(self%yv)) deallocate(self%yv)
    allocate(self%xv(xCount                 ))
    allocate(self%yv(xCount,tableCountActual))
    self%xv            =Make_Range(xMinimum,xMaximum,xCount,rangeType=rangeTypeLinear)
    self%inverseDeltaX =1.0d0/(self%xv(2)-self%xv(1))
    self%tablePrevious =-1
    self%dTablePrevious=-1
    self%xPrevious     =-1.0d0
    self%dxPrevious    =-1.0d0
    ! Set extrapolation type.
    if (present(extrapolationType)) then
       self%extrapolationType=extrapolationType
    else
       self%extrapolationType=extrapolationTypeExtrapolate
    end if
    return
  end subroutine Table_Linear_1D_Create

  subroutine Table_Linear_1D_Populate(self,y,table)
    !!{
    Populate a 1-D linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DLinearLinear)              , intent(inout)           :: self
    double precision                     , dimension(:), intent(in   )           :: y
    integer                                            , intent(in   ), optional :: table
    integer                                                                      :: tableActual

    ! Validate the input.
    if (.not.allocated(self%yv)       ) call Error_Report("create the table before populating it"//{introspection:location})
    if (size(self%yv,dim=1) /= size(y)) call Error_Report("provided y array is of wrong size"    //{introspection:location})

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
    !!{
    Populate a single element of a 1-D linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DLinearLinear), intent(inout)           :: self
    double precision                     , intent(in   )           :: y
    integer                              , intent(in   )           :: i
    integer                              , intent(in   ), optional :: table
    integer                                                        :: tableActual

    ! Validate the input.
    if (.not.allocated(self%yv)           ) call Error_Report("create the table before populating it"//{introspection:location})
    if (i < 1 .or. i > size(self%yv,dim=1)) call Error_Report("provided i value is out of bounds"    //{introspection:location})

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

  double precision function Table_Linear_1D_Interpolate(self,x,table,status)
    !!{
    Perform linear interpolation in a linear 1D table.
    !!}
    use :: Table_Labels, only : extrapolationTypeZero
    implicit none
    class           (table1DLinearLinear), intent(inout)           :: self
    double precision                     , intent(in   )           :: x
    integer                              , intent(in   ), optional :: table
    integer                              , intent(  out), optional :: status
    integer                                                        :: i    , tableActual
    double precision                                               :: h    , xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    xEffective=self%xEffective(x,status)
    if (xEffective /= self%xPrevious .or. tableActual /= self%tablePrevious) then
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          if (self%extrapolationType(1) == extrapolationTypeZero) then
             Table_Linear_1D_Interpolate=0.0d0
             return
          end if
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          if (self%extrapolationType(2) == extrapolationTypeZero) then
             Table_Linear_1D_Interpolate=0.0d0
             return
          end if
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

  double precision function Table_Linear_1D_Interpolate_Gradient(self,x,table,status)
    !!{
    Perform linear interpolation in a linear 1D table.
    !!}
    use :: Table_Labels, only : extrapolationTypeZero
    implicit none
    class           (table1DLinearLinear), intent(inout)           :: self
    double precision                     , intent(in   )           :: x
    integer                              , intent(in   ), optional :: table
    integer                              , intent(  out), optional :: status
    integer                                                        :: i         , tableActual
    double precision                                               :: xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    xEffective=self%xEffective(x,status)
    if (xEffective /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          if (self%extrapolationType(1) == extrapolationTypeZero) then
             Table_Linear_1D_Interpolate_Gradient=0.0d0
             return
          end if
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          if (self%extrapolationType(2) == extrapolationTypeZero) then
             Table_Linear_1D_Interpolate_Gradient=0.0d0
             return
          end if
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
    !!{
    Create a 1-D logarithmic table.
    !!}
    implicit none
    class           (table1DLogarithmicLinear        ), intent(inout)                         :: self
    double precision                                  , intent(in   )                         :: xMaximum         , xMinimum
    integer                                           , intent(in   )                         :: xCount
    integer                                           , intent(in   ), optional               :: tableCount
    type            (enumerationExtrapolationTypeType), intent(in   ), optional, dimension(2) :: extrapolationType

    self%previousSet         =.false.
    self%xLinearPrevious     =-1.0d0
    self%xLogarithmicPrevious=-1.0d0
    ! Call the creator for linear tables with the logarithms of the input x range.
    call self%table1DLinearLinear%create(log(xMinimum),log(xMaximum),xCount,tableCount,extrapolationType)
    return
  end subroutine Table_Logarithmic_1D_Create

  double precision function Table_Logarithmic_1D_X(self,i)
    !!{
    Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $x$-value for a logarithmic 1D table.
    !!}
    implicit none
    class  (table1DLogarithmicLinear), intent(inout) :: self
    integer                          , intent(in   ) :: i

    Table_Logarithmic_1D_X=exp(self%table1DLinearLinear%x(i))
    return
  end function Table_Logarithmic_1D_X

  function Table_Logarithmic_1D_Xs(self)
    !!{
    Return the $x$-values for a 1D table.
    !!}
    implicit none
    class           (table1DLogarithmicLinear), intent(in   )             :: self
    double precision                          , dimension(size(self%xv))  :: Table_Logarithmic_1D_Xs

    Table_Logarithmic_1D_Xs=exp(self%table1DLinearLinear%xs())
    return
  end function Table_Logarithmic_1D_Xs

  double precision function Table_Logarithmic_1D_Interpolate(self,x,table,status)
    !!{
    Perform linear interpolation in a logarithmic 1D table.
    !!}
    implicit none
    class           (table1DLogarithmicLinear), intent(inout)           :: self
    double precision                          , intent(in   )           :: x
    integer                                   , intent(in   ), optional :: table
    integer                                   , intent(  out), optional :: status

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_1D_Interpolate=self%table1DLinearLinear%interpolate(self%xLogarithmicPrevious,table,status)
    return
  end function Table_Logarithmic_1D_Interpolate

  double precision function Table_Logarithmic_1D_Interpolate_Gradient(self,x,table,status)
    !!{
    Perform linear interpolation in a logarithmic 1D table.
    !!}
    implicit none
    class           (table1DLogarithmicLinear), intent(inout)           :: self
    double precision                          , intent(in   )           :: x
    integer                                   , intent(in   ), optional :: table
    integer                                   , intent(  out), optional :: status

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_1D_Interpolate_Gradient=self%table1DLinearLinear%interpolateGradient(self%xLogarithmicPrevious,table,status)/self%xEffective(x,status)
    return
  end function Table_Logarithmic_1D_Interpolate_Gradient

  function Table_Logarithmic_Integration_Weights(self,x0,x1,integrand)
    !!{
    Returns a set of weights for trapezoidal integration on the table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
    !!}
    use :: Error     , only : Error_Report
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (table1DLogarithmicLinear       ), intent(inout)                               :: self
    double precision                                 , intent(in   )                               :: x0, x1
    procedure       (tablesIntegrationWeightFunction), intent(in   )           , pointer, optional :: integrand
    double precision                                 , dimension(size(self%xv))                    :: Table_Logarithmic_Integration_Weights
    double precision                                 , parameter                                   :: logTolerance                         =1.0d-12
    double precision                                                                               :: gradientTerm                                 , lx0        , &
         &                                                                                            lx1                                          , factor0    , &
         &                                                                                            factor1
    integer                                                                                        :: i
    type            (integrator                     )                                              :: integrator0                                  , integrator1

    if (x1 < x0) call Error_Report('inverted limits'//{introspection:location})
    integrator0=integrator(factor0Integrand,toleranceRelative=1.0d-4)
    integrator1=integrator(factor1Integrand,toleranceRelative=1.0d-4)
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
             factor0=integrator0%integrate(lx0,lx1)
             factor1=integrator1%integrate(lx0,lx1)
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

    double precision function factor0Integrand(logx)
      !!{
      Integrand used to evaluate integration weights over logarithmically spaced tables
      !!}
      implicit none
      double precision, intent(in   ) :: logx
      double precision                :: x

      x=exp(logx)
      factor0Integrand=x*integrand(x)
      return
    end function factor0Integrand

    double precision function factor1Integrand(logx)
      !!{
      Integrand used to evaluate integration weights over logarithmically spaced tables
      !!}
      implicit none
      double precision, intent(in   ) :: logx
      double precision                :: x

      x=exp(logx)
      factor1Integrand=x*integrand(x)*(logx-self%xv(i-1))/(self%xv(i)-self%xv(i-1))
      return
    end function factor1Integrand

  end function Table_Logarithmic_Integration_Weights

  subroutine Table_Logarithmic_1D_Reverse(self,reversedSelf,table,precise)
    !!{
    Reverse a 1D logarithmic-linear table (i.e. swap $x$ and $y$ components). Optionally allows specification of
    which $y$ table to swap with.
    !!}
    use :: Array_Utilities, only : Array_Is_Monotonic, Array_Reverse, directionDecreasing
    use :: Error          , only : Error_Report
    implicit none
    class  (table1DLogarithmicLinear)             , intent(in   )           :: self
    class  (table1D                 ), allocatable, intent(inout)           :: reversedSelf
    integer                                       , intent(in   ), optional :: table
    logical                                       , intent(in   ), optional :: precise
    integer                                                                 :: i           , tableActual
    !$GLC attributes unused :: precise

    tableActual=1
    if (present(table)) tableActual=table
    if (.not.Array_Is_Monotonic(self%yv(:,tableActual))) call Error_Report('reversed table would not be monotonic'//{introspection:location})
    if (allocated(reversedSelf)) then
       call reversedSelf%destroy()
       deallocate(reversedSelf)
    end if
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
    !!{
    Create a 1-D linear table.
    !!}
    use :: Error            , only : Error_Report
    use :: Numerical_Ranges , only : Make_Range                  , rangeTypeLinear
    use :: Table_Labels     , only : extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (table1DLinearCSpline            ), intent(inout)                         :: self
    double precision                                  , intent(in   )                         :: xMaximum         , xMinimum
    integer                                           , intent(in   )                         :: xCount
    integer                                           , intent(in   ), optional               :: tableCount
    type            (enumerationExtrapolationTypeType), intent(in   ), optional, dimension(2) :: extrapolationType
    integer                                                                                   :: tableCountActual

    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the x-range.
    self%xCount=xCount
    if (allocated(self%xv)) deallocate(self%xv)
    if (allocated(self%yv)) deallocate(self%yv)
    if (allocated(self%sv)) deallocate(self%sv)
    if (allocated(self%av)) deallocate(self%av)
    if (allocated(self%bv)) deallocate(self%bv)
    if (allocated(self%cv)) deallocate(self%cv)
    if (allocated(self%dv)) deallocate(self%dv)
    allocate(self%xv(xCount                   ))
    allocate(self%yv(xCount  ,tableCountActual))
    allocate(self%sv(xCount  ,tableCountActual))
    allocate(self%av(xCount-1,tableCountActual))
    allocate(self%bv(xCount-1,tableCountActual))
    allocate(self%cv(xCount-1,tableCountActual))
    allocate(self%dv(xCount-1,tableCountActual))
    self%xv            =Make_Range(xMinimum,xMaximum,xCount,rangeType=rangeTypeLinear)
    self%       deltaX =self%xv(2)-self%xv(1)
    self%inverseDeltaX =1.0d0/self%deltaX
    self%tablePrevious =-1
    self%dTablePrevious=-1
    self%iPrevious     =-1
    self%xPrevious     =-1.0d0
    self%dxPrevious    =-1.0d0
    ! Set extrapolation type.
    if (present(extrapolationType)) then
       if (any(extrapolationType == extrapolationTypeZero)) call Error_Report('zero extrapolation is not supported'//{introspection:location})
       self%extrapolationType=extrapolationType
    else
       self%extrapolationType=extrapolationTypeExtrapolate
    end if
    return
  end subroutine Table_Linear_CSpline_1D_Create

  subroutine Table_Linear_CSpline_1D_Destroy(self)
    !!{
    Destroy a linear cubic-spline 1-D table.
    !!}
    implicit none
    class(table1DLinearCSpline), intent(inout) :: self

    call Table_1D_Destroy(self)
    if (allocated(self%sv)) deallocate(self%sv)
    if (allocated(self%av)) deallocate(self%av)
    if (allocated(self%bv)) deallocate(self%bv)
    if (allocated(self%cv)) deallocate(self%cv)
    if (allocated(self%dv)) deallocate(self%dv)
    return
  end subroutine Table_Linear_CSpline_1D_Destroy

  subroutine Table_Linear_CSpline_1D_Populate(self,y,table,computeSpline)
    !!{
    Populate a 1-D linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DLinearCSpline)              , intent(inout)           :: self
    double precision                      , dimension(:), intent(in   )           :: y
    integer                                             , intent(in   ), optional :: table
    logical                                             , intent(in   ), optional :: computeSpline
    integer                                                                       :: tableActual
    logical                                                                       :: computeSplineActual

    ! Validate the input.
    if (.not.(allocated(self%yv).and.allocated(self%sv))) call Error_Report("create the table before populating it"//{introspection:location})
    if (size(self%yv,dim=1) /= size(y)                  ) call Error_Report("provided y array is of wrong size"    //{introspection:location})

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
    !!{
    Populate a single element of a 1-D linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DLinearCSpline), intent(inout)           :: self
    double precision                      , intent(in   )           :: y
    integer                               , intent(in   )           :: i
    integer                               , intent(in   ), optional :: table
    logical                               , intent(in   ), optional :: computeSpline
    integer                                                         :: tableActual
    logical                                                         :: computeSplineActual

    ! Validate the input.
    if (.not.(allocated(self%yv).and.allocated(self%sv))) call Error_Report("create the table before populating it"//{introspection:location})
    if (i < 1 .or. i > size(self%yv,dim=1)              ) call Error_Report("provided i value is out of bounds"    //{introspection:location})

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
    !!{
    Compute the interpolating spline factors for a 1-D linear spline.
    !!}
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
    ! Compute polynomial coefficients.
    do i=1,size(self%xv)-1
       self%av(i,table)=                      self%yv(i  ,table)
       self%bv(i,table)=-self%       deltaX*  self%sv(i+1,table)/6.0d0 &
            &           -self%       deltaX*  self%sv(i  ,table)/3.0d0 &
            &           +self%inverseDeltaX*( self%yv(i+1,table)       &
            &                                -self%yv(i  ,table)       &
            &                             )
       self%cv(i,table)=                      self%sv(i  ,table)/2.0d0
       self%dv(i,table)= self%inverseDeltaX*(                          &
            &                                 self%sv(i+1,table)       &
            &                                -self%sv(i  ,table)       &
            &                               )                   /6.0d0
    end do
    return
  end subroutine Table_Linear_CSpline_1D_Compute_Spline

  double precision function Table_Linear_CSpline_1D_Interpolate(self,x,table,status)
    !!{
    Perform linear interpolation in a linear 1D table.
    !!}
    implicit none
    class           (table1DLinearCSpline), intent(inout)           :: self
    double precision                      , intent(in   )           :: x
    integer                               , intent(in   ), optional :: table
    integer                               , intent(  out), optional :: status
    integer                                                         :: i         , tableActual
    double precision                                                :: xEffective, dx

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%xPrevious .or. tableActual /= self%tablePrevious) then
       xEffective=self%xEffective(x,status)
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=                                                    1
       else if (xEffective >= self%xv(self%xCount)) then
          i=                                                      self%xCount-1
       else
          i=min(int((xEffective-self%xv(1))*self%inverseDeltaX)+1,self%xCount-1)
       end if
       ! Compute offset from tabulated point.
       dx=xEffective-self%xv(i)
       ! Interpolate in the table.
       self%xPrevious    =x
       self%tablePrevious=tableActual
       self%    yPrevious=self%av(i,tableActual)+dx*(self%bv(i,tableActual)+dx*(self%cv(i,tableActual)+dx*self%dv(i,tableActual)))
    end if
    Table_Linear_CSpline_1D_Interpolate=self%yPrevious
    return
  end function Table_Linear_CSpline_1D_Interpolate

  double precision function Table_Linear_CSpline_1D_Interpolate_Gradient(self,x,table,status)
    !!{
    Perform linear interpolation in a linear 1D table.
    !!}
    implicit none
    class           (table1DLinearCSpline), intent(inout)           :: self
    double precision                      , intent(in   )           :: x
    integer                               , intent(in   ), optional :: table
    integer                               , intent(  out), optional :: status
    integer                                                         :: i         , tableActual
    double precision                                                :: xEffective, dx

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then
       xEffective=self%xEffective(x,status)
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=                                                    1
       else if (xEffective >= self%xv(self%xCount)) then
          i=                                                      self%xCount-1
       else
          i=min(int((xEffective-self%xv(1))*self%inverseDeltaX)+1,self%xCount-1)
       end if
       ! Compute offset from tabulated point.
       dx=xEffective-self%xv(i)
       ! Interpolate in the table.
       self%dxPrevious    =x
       self%dTablePrevious=tableActual
       self%    dyPrevious=self%bv(i,tableActual)+dx*(2.0d0*self%cv(i,tableActual)+dx*3.0d0*self%dv(i,tableActual))
    end if
    Table_Linear_CSpline_1D_Interpolate_Gradient=self%dyPrevious
    return
  end function Table_Linear_CSpline_1D_Interpolate_Gradient

  subroutine Table_Logarithmic_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !!{
    Create a 1-D logarithmic table.
    !!}
    implicit none
    class           (table1DLogarithmicCSpline       ), intent(inout)                         :: self
    double precision                                  , intent(in   )                         :: xMaximum         , xMinimum
    integer                                           , intent(in   )                         :: xCount
    integer                                           , intent(in   ), optional               :: tableCount
    type            (enumerationExtrapolationTypeType), intent(in   ), optional, dimension(2) :: extrapolationType

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
    !!{
    Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $x$-value for a logarithmic 1D table.
    !!}
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
    !!{
    Return the $x$-values for a 1D table.
    !!}
    implicit none
    class           (table1DLogarithmicCSpline), intent(in   )            :: self
    double precision                           , dimension(size(self%xv)) :: Table_Logarithmic_CSpline_1D_Xs

    Table_Logarithmic_CSpline_1D_Xs=exp(self%table1DLinearCSpline%xs())
    return
  end function Table_Logarithmic_CSpline_1D_Xs

  double precision function Table_Logarithmic_CSpline_1D_Interpolate(self,x,table,status)
    !!{
    Perform linear interpolation in a logarithmic 1D table.
    !!}
    implicit none
    class           (table1DLogarithmicCSpline), intent(inout)           :: self
    double precision                           , intent(in   )           :: x
    integer                                    , intent(in   ), optional :: table
    integer                                    , intent(  out), optional :: status

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_CSpline_1D_Interpolate=self%table1DLinearCSpline%interpolate(self%xLogarithmicPrevious,table,status)
    return
  end function Table_Logarithmic_CSpline_1D_Interpolate

  double precision function Table_Logarithmic_CSpline_1D_Interpolate_Gradient(self,x,table,status)
    !!{
    Perform linear interpolation in a logarithmic 1D table.
    !!}
    implicit none
    class           (table1DLogarithmicCSpline), intent(inout)           :: self
    double precision                           , intent(in   )           :: x
    integer                                    , intent(in   ), optional :: table
    integer                                    , intent(  out), optional :: status

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_CSpline_1D_Interpolate_Gradient=self%table1DLinearCSpline%interpolateGradient(self%xLogarithmicPrevious,table,status)/self%xEffective(x,status)
    return
  end function Table_Logarithmic_CSpline_1D_Interpolate_Gradient

  function Table_Linear_CSpline_Integration_Weights(self,x0,x1,integrand)
    !!{
    Returns a set of weights for trapezoidal integration on the table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DLinearCSpline           ), intent(inout)                               :: self
    double precision                                 , intent(in   )                               :: x0, x1
    procedure       (tablesIntegrationWeightFunction), intent(in   )           , pointer, optional :: integrand
    double precision                                 , dimension(size(self%xv))                    :: Table_Linear_CSpline_Integration_Weights
    !$GLC attributes unused :: self, x0, x1, integrand

    Table_Linear_CSpline_Integration_Weights=0.0d0
    call Error_Report('integration weights not supported'//{introspection:location})
    return
  end function Table_Linear_CSpline_Integration_Weights

  subroutine Table_Monotone_CSpline_1D_Create(self,x,tableCount,extrapolationType)
    !!{
    Create a 1-D monotone cubic spline table.
    !!}
    use :: Error            , only : Error_Report
    use :: Table_Labels     , only : extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (table1DMonotoneCSpline          ), intent(inout)                         :: self
    double precision                                  , intent(in   )          , dimension(:) :: x
    integer                                           , intent(in   ), optional               :: tableCount
    type            (enumerationExtrapolationTypeType), intent(in   ), optional, dimension(2) :: extrapolationType
    integer                                                                                   :: tableCountActual

    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the x-range.
    self%xCount=size(x)
    allocate(self%xv(self%xCount                   ))
    allocate(self%yv(self%xCount  ,tableCountActual))
    allocate(self%av(self%xCount  ,tableCountActual))
    allocate(self%bv(self%xCount-1,tableCountActual))
    allocate(self%cv(self%xCount-1,tableCountActual))
    self%xv            =x
    self%tablePrevious =-1
    self%dTablePrevious=-1
    self%iPrevious     =-1
    self%xPrevious     =-1.0d0
    self%dxPrevious    =-1.0d0
    ! Set extrapolation type.
    if (present(extrapolationType)) then
       if (any(extrapolationType == extrapolationTypeZero)) call Error_Report('zero extrapolation is not supported'//{introspection:location})
       self%extrapolationType=extrapolationType
    else
       self%extrapolationType=extrapolationTypeExtrapolate
    end if
    return
  end subroutine Table_Monotone_CSpline_1D_Create

  subroutine Table_Monotone_CSpline_1D_Destroy(self)
    !!{
    Destroy a monotone cubic-spline 1-D table.
    !!}
    implicit none
   class(table1DMonotoneCSpline), intent(inout) :: self

    call Table_1D_Destroy(self)
    if (allocated(self%av)) deallocate(self%av)
    if (allocated(self%bv)) deallocate(self%bv)
    if (allocated(self%cv)) deallocate(self%cv)
    return
  end subroutine Table_Monotone_CSpline_1D_Destroy

  subroutine Table_Monotone_CSpline_1D_Populate(self,y,table,computeSpline)
    !!{
    Populate a 1-D monotone cubic spline table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DMonotoneCSpline)              , intent(inout)           :: self
    double precision                        , dimension(:), intent(in   )           :: y
    integer                                               , intent(in   ), optional :: table
    logical                                               , intent(in   ), optional :: computeSpline
    integer                                                                         :: tableActual
    logical                                                                         :: computeSplineActual

    ! Validate the input.
    if (.not.(allocated(self%yv))) call Error_Report("create the table before populating it"//{introspection:location})
    if (size(self%yv,dim=1) /= size(y) ) call Error_Report("provided y array is of wrong size"    //{introspection:location})

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table

    ! Store the y values.
    self%yv(:,tableActual)=y

    ! Compute the spline interpolation for this table.
    computeSplineActual=.true.
    if (present(computeSpline)) computeSplineActual=computeSpline
    if (computeSplineActual) call Table_Monotone_CSpline_1D_Compute_Spline(self,tableActual)
    return
  end subroutine Table_Monotone_CSpline_1D_Populate

  subroutine Table_Monotone_CSpline_1D_Populate_Single(self,y,i,table,computeSpline)
    !!{
    Populate a single element of a 1-D monotone cubic spline table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DMonotoneCSpline), intent(inout)           :: self
    double precision                        , intent(in   )           :: y
    integer                                 , intent(in   )           :: i
    integer                                 , intent(in   ), optional :: table
    logical                                 , intent(in   ), optional :: computeSpline
    integer                                                           :: tableActual
    logical                                                           :: computeSplineActual

    ! Validate the input.
    if (.not.(allocated(self%yv))) call Error_Report("create the table before populating it"//{introspection:location})
    if (i < 1 .or. i > size(self%yv,dim=1)) call Error_Report("provided i value is out of bounds"    //{introspection:location})

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table

    ! Store the y values.
    self%yv(i,tableActual)=y

    ! Compute the spline interpolation for this table.
    computeSplineActual=.true.
    if (present(computeSpline)) computeSplineActual=computeSpline
    if (computeSplineActual) call Table_Monotone_CSpline_1D_Compute_Spline(self,tableActual)
    return
  end subroutine Table_Monotone_CSpline_1D_Populate_Single

  subroutine Table_Monotone_CSpline_1D_Compute_Spline(self,table)
    !!{
    Compute the interpolating spline factors for a 1-D linear spline.
    !!}
    implicit none
    type            (table1DMonotoneCSpline), intent(inout)               :: self
    integer                                 , intent(in   )               :: table
    double precision                        , allocatable  , dimension(:) :: dx   , dy    , m
    integer                                                               :: i
    double precision                                                      :: dxSum, factor

    ! Reset all previously stored values.
    self% tablePrevious=-1
    self%dTablePrevious=-1
    self%     iPrevious=-1
    ! Allocate workspace.
    allocate(dx(size(self%xv)-1))
    allocate(dy(size(self%xv)-1))
    allocate( m(size(self%xv)-1))
    ! Get consecutive differences and slopes.
    do i=1,size(self%xv)-1
       dx(i)=self%xv(i+1      )-self%xv(i      )
       dy(i)=self%yv(i+1,table)-self%yv(i,table)
       m (i)=dy(i)/dx(i)
    end do
    ! Get degree-1 coefficients.
    self%av(1,table)=m(1)
    do i=1,size(dx)-1
       if (m(i)*m(i+1) <= 0.0d0) then
          self%av(i+1,table)=0.0d0
       else
          dxSum=dx(i)+dx(i+1)
          self%av(i+1,table)=3.0d0*dxSum/((dxSum+dx(i+1))/m(i)+(dxSum+dx(i))/m(i+1))
       end if
    end do
    self%av(size(dx)+1,table)=m(size(m))
    ! Get degree-2 and degree-3 coefficients.
    do i=1,size(self%av)-1
       factor=self%av(i,table)+self%av(i+1,table)-2.0d0*m(i)
       self%bv(i,table)=(m(i)-self%av(i,table)-factor)/dx(i)
       self%cv(i,table)=factor/dx(i)**2
    end do
    ! Destroy workspace.
    deallocate(dx)
    deallocate(dy)
    deallocate( m)
    return
  end subroutine Table_Monotone_CSpline_1D_Compute_Spline

  double precision function Table_Monotone_CSpline_1D_Interpolate(self,x,table,status)
    !!{
    Perform monotonic cubic spline interpolation in a 1D table.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Arrays_Search, only : searchArray
    implicit none
    class           (table1DMonotoneCSpline), intent(inout)           :: self
    double precision                        , intent(in   )           :: x
    integer                                 , intent(in   ), optional :: table
    integer                                 , intent(  out), optional :: status
    integer                                                           :: tableActual
    integer         (c_size_t              )                          :: i
    double precision                                                  :: dx   , xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%xPrevious .or. tableActual /= self%tablePrevious) then
       xEffective=self%xEffective(x,status)
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=searchArray(self%xv,xEffective)
       end if
       ! Compute offset from tabulated point.
       dx=xEffective-self%xv(i)
       ! Interpolate in the table.
       self%xPrevious    =x
       self%tablePrevious=tableActual
       self%    yPrevious=self%yv(i,tableActual)+self%av(i,tableActual)*dx+self%bv(i,tableActual)*dx**2+self%cv(i,tableActual)*dx**3
    end if
    Table_Monotone_CSpline_1D_Interpolate=self%yPrevious
    return
  end function Table_Monotone_CSpline_1D_Interpolate

  double precision function Table_Monotone_CSpline_1D_Interpolate_Gradient(self,x,table,status)
    !!{
    Perform monotonic cubic spline interpolation in a 1D table and return the gradient.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Arrays_Search, only : searchArray
    implicit none
    class           (table1DMonotoneCSpline), intent(inout)           :: self
    double precision                        , intent(in   )           :: x
    integer                                 , intent(in   ), optional :: table
    integer                                 , intent(  out), optional :: status
    integer                                                           :: tableActual
    integer         (c_size_t              )                          :: i
    double precision                                                  :: dx   , xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then
       xEffective=self%xEffective(x,status)
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=searchArray(self%xv,xEffective)
       end if
       ! Compute offset from tabulated point.
       dx=xEffective-self%xv(i)
       ! Interpolate in the table.
       self%dxPrevious    =x
       self%dTablePrevious=tableActual
       self%    dyPrevious=self%av(i,tableActual)+2.0d0*self%bv(i,tableActual)*dx+3.0d0*self%cv(i,tableActual)*dx**2
    end if
    Table_Monotone_CSpline_1D_Interpolate_Gradient=self%dyPrevious
    return
  end function Table_Monotone_CSpline_1D_Interpolate_Gradient

  function Table_Monotone_CSpline_Integration_Weights(self,x0,x1,integrand)
    !!{
    Returns a set of weights for trapezoidal integration on the table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DMonotoneCSpline         ), intent(inout)                               :: self
    double precision                                 , intent(in   )                               :: x0, x1
    procedure       (tablesIntegrationWeightFunction), intent(in   )           , pointer, optional :: integrand
    double precision                                 , dimension(size(self%xv))                    :: Table_Monotone_CSpline_Integration_Weights
    !$GLC attributes unused :: self, x0, x1, integrand

    Table_Monotone_CSpline_Integration_Weights=0.0d0
    call Error_Report('integration weights not supported'//{introspection:location})
    return
  end function Table_Monotone_CSpline_Integration_Weights
  
  double precision function Table1D_Find_Effective_X(self,x,status)
    !!{
    Return the effective value of $x$ to use in table interpolations.
    !!}
    use :: Error       , only : Error_Report                , errorStatusOutOfRange, errorStatusSuccess
    use :: Table_Labels, only : extrapolationTypeExtrapolate, extrapolationTypeFix , extrapolationTypeZero, extrapolationTypeAbort
    implicit none
    class           (table1D), intent(inout)           :: self
    double precision         , intent(in   )           :: x
    integer                  , intent(  out), optional :: status

    if (present(status)) status=errorStatusSuccess
    if      (x < self%xv(1         )) then
       select case (self%extrapolationType(1)%ID)
       case (extrapolationTypeExtrapolate%ID,extrapolationTypeZero%ID)
          Table1D_Find_Effective_X=     x
       case (extrapolationTypeFix        %ID)
          Table1D_Find_Effective_X=self%xv(1          )
       case (extrapolationTypeAbort      %ID)
          Table1D_Find_Effective_X=     x
          if (present(status)) then
             status=errorStatusOutOfRange
          else
             call Error_Report('x is below range - extrapolation type "abort" forbids extrapolation'//{introspection:location})
          end if
       case default
          Table1D_Find_Effective_X=0.0d0
          call    Error_Report('x is below range - unknown extrapolation method'                    //{introspection:location})
       end select
    else if (x > self%x(self%xCount)) then
       select case (self%extrapolationType(2)%ID)
       case (extrapolationTypeExtrapolate%ID,extrapolationTypeZero%ID)
          Table1D_Find_Effective_X=     x
       case (extrapolationTypeFix        %ID)
          Table1D_Find_Effective_X=self%xv(self%xCount)
       case (extrapolationTypeAbort      %ID)
          Table1D_Find_Effective_X=     x
          if (present(status)) then
             status=errorStatusOutOfRange
          else
             call Error_Report('x is above range - extrapolation type "abort" forbids extrapolation'//{introspection:location})
          end if
       case default
          Table1D_Find_Effective_X=0.0d0
          call    Error_Report('x is above range - unknown extrapolation method'                    //{introspection:location})      
       end select
    else
       Table1D_Find_Effective_X=x
    end if
    return
  end function Table1D_Find_Effective_X

  subroutine Table_NonUniform_Linear_Logarithmic_1D_Populate(self,y,table)
    !!{
    Populate a 1-D linear-logarithmic table.
    !!}
    implicit none
    class           (table1DNonUniformLinearLogarithmic)              , intent(inout)           :: self
    double precision                                    , dimension(:), intent(in   )           :: y
    integer                                                           , intent(in   ), optional :: table

    call self%table1DGeneric%populate(log(y),table)
    return
  end subroutine Table_NonUniform_Linear_Logarithmic_1D_Populate

  subroutine Table_NonUniform_Linear_Logarithmic_1D_Populate_Single(self,y,i,table)
    !!{
    Populate a single element of a 1-D linear table.
    !!}
    implicit none
    class           (table1DNonUniformLinearLogarithmic), intent(inout)           :: self
    double precision                                    , intent(in   )           :: y
    integer                                             , intent(in   )           :: i
    integer                                             , intent(in   ), optional :: table

    call self%table1DGeneric%populate(log(y),i,table)
    return
  end subroutine Table_NonUniform_Linear_Logarithmic_1D_Populate_Single

  double precision function Table_NonUniform_Linear_Logarithmic_1D_Interpolate(self,x,table,status)
    !!{
    Perform linear interpolation in a linear-logarithmic 1D table.
    !!}
    implicit none
    class           (table1DNonUniformLinearLogarithmic), intent(inout)           :: self
    double precision                                    , intent(in   )           :: x
    integer                                             , intent(in   ), optional :: table
    integer                                             , intent(  out), optional :: status

    Table_NonUniform_Linear_Logarithmic_1D_Interpolate=exp(self%table1DGeneric%interpolate(x,table,status))
    return
  end function Table_NonUniform_Linear_Logarithmic_1D_Interpolate

  double precision function Table_NonUniform_Linear_Logarithmic_1D_Interpolate_Gradient(self,x,table,status)
    !!{
    Perform linear interpolation in a linear-logarithmic 1D table.
    !!}
    implicit none
    class           (table1DNonUniformLinearLogarithmic), intent(inout)           :: self
    double precision                                    , intent(in   )           :: x
    integer                                             , intent(in   ), optional :: table
    integer                                             , intent(  out), optional :: status

    Table_NonUniform_Linear_Logarithmic_1D_Interpolate_Gradient=self%interpolate(x,table,status)*self%table1DGeneric%interpolateGradient(x,table,status)
    return
  end function Table_NonUniform_Linear_Logarithmic_1D_Interpolate_Gradient

  function Table_NonUniform_Linear_Logarithmic_Integration_Weights(self,x0,x1,integrand)
    !!{
    Returns a set of weights for integration on a linear-logarithmic table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DNonUniformLinearLogarithmic ), intent(inout)                               :: self
    double precision                                     , intent(in   )                               :: x0, x1
    procedure       (tablesIntegrationWeightFunction    ), intent(in   )           , pointer, optional :: integrand
    double precision                                     , dimension(size(self%xv))                    :: Table_NonUniform_Linear_Logarithmic_Integration_Weights
    !$GLC attributes unused :: self, x0, x1, integrand

    Table_NonUniform_Linear_Logarithmic_Integration_Weights=0.0d0
    call Error_Report('integrand is not linear in y'//{introspection:location})
    return
  end function Table_NonUniform_Linear_Logarithmic_Integration_Weights

  double precision function Table_NonUniform_Linear_Logarithmic_1D_Y(self,i,table)
    !!{
    Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $y$-value for a 1D table.
    !!}
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
    !!{
    Return the $y$-values for a 1D table.
    !!}
    implicit none
    class           (table1DNonUniformLinearLogarithmic), intent(in   )                                      :: self
    double precision                                    , dimension(size(self%yv,dim=1),size(self%yv,dim=2)) :: Table_NonUniform_Linear_Logarithmic_1D_Ys

    Table_NonUniform_Linear_Logarithmic_1D_Ys=exp(self%yv)
    return
  end function Table_NonUniform_Linear_Logarithmic_1D_Ys
  
  subroutine Table_2DLogLogLin_Create(self,xMinimum,xMaximum,xCount,yMinimum,yMaximum,yCount,tableCount,extrapolationTypeX,extrapolationTypeY)
    !!{
    Create a 2-D log-log-linear table.
    !!}
    use :: Numerical_Ranges , only : Make_Range                  , rangeTypeLinear
    use :: Table_Labels     , only : extrapolationTypeExtrapolate
    implicit none
    class           (table2DLogLogLin                ), intent(inout)           :: self
    double precision                                  , intent(in   )           :: xMaximum          , xMinimum          , &
         &                                                                         yMaximum          , yMinimum
    integer                                           , intent(in   )           :: xCount            , yCount
    type            (enumerationExtrapolationTypeType), intent(in   ), optional :: extrapolationTypeX, extrapolationTypeY
    integer                                           , intent(in   ), optional :: tableCount
    integer                                                                     :: tableCountActual

    ! Initialize state.
    self%xPreviousSet        =.false.
    self%yPreviousSet        =.false.
    self%xLinearPrevious     =-1.0d0
    self%yLinearPrevious     =-1.0d0
    self%xLogarithmicPrevious=-1.0d0
    self%yLogarithmicPrevious=-1.0d0
    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the ranges.
    self%xCount=xCount
    self%yCount=yCount
    allocate(self%xv(xCount                        ))
    allocate(self%yv(yCount                 ))
    allocate(self%zv(xCount,yCount,tableCountActual))
    self%xv                  =Make_Range(log(xMinimum),log(xMaximum),xCount,rangeType=rangeTypeLinear)
    self%yv                  =Make_Range(log(yMinimum),log(yMaximum),yCount,rangeType=rangeTypeLinear)
    self%inverseDeltaX       =1.0d0/(self%xv(2)-self%xv(1))
    self%inverseDeltaY       =1.0d0/(self%yv(2)-self%yv(1))
    self%tablePrevious       =-1
    self%hx                  =-1.0d0
    self%xLinearPrevious     =-1.0d0
    self%xLogarithmicPrevious=-1.0d0
    self%hy                  =-1.0d0
    self%yLinearPrevious     =-1.0d0
    self%yLogarithmicPrevious=-1.0d0
    ! Set extrapolation type.
    if (present(extrapolationTypeX)) then
       self%extrapolationTypeX=extrapolationTypeX
    else
       self%extrapolationTypeX=extrapolationTypeExtrapolate
    end if
    if (present(extrapolationTypeY)) then
       self%extrapolationTypeY=extrapolationTypeY
    else
       self%extrapolationTypeY=extrapolationTypeExtrapolate
    end if
    return
  end subroutine Table_2DLogLogLin_Create

  double precision function Table_2DLogLogLin_X(self,i)
    !!{
    Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $x$-value for a 2D log-log table.
    !!}
    implicit none
    class  (table2DLogLogLin), intent(inout) :: self
    integer                  , intent(in   ) :: i

    Table_2DLogLogLin_X=exp(self%xv(i))
    return
  end function Table_2DLogLogLin_X

  double precision function Table_2DLogLogLin_Y(self,i)
    !!{
    Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $y$-value for a 2D log-log table.
    !!}
    implicit none
    class  (table2DLogLogLin), intent(inout) :: self
    integer                  , intent(in   ) :: i

    Table_2DLogLogLin_Y=exp(self%yv(i))
    return
  end function Table_2DLogLogLin_Y

  double precision function Table_2DLogLogLin_Z(self,i,j,table)
    !!{
    Return the {\normalfont \ttfamily (i,j)}$^\mathrm{th}$ $x$-value for a 2D log-log table.
    !!}
    implicit none
    class  (table2DLogLogLin), intent(inout)           :: self
    integer                  , intent(in   )           :: i          , j
    integer                  , intent(in   ), optional :: table
    integer                                            :: tableActual

    tableActual=1
    if (present(table)) tableActual=table
    Table_2DLogLogLin_Z=self%zv(i,j,tableActual)
    return
  end function Table_2DLogLogLin_Z

  function Table_2DLogLogLin_Xs(self)
    !!{
    Return the $x$-values for a 2D log-log table.
    !!}
    implicit none
    class(table2DLogLogLin), intent(in   )             :: self
    double precision       , dimension(size(self%xv))  :: Table_2DLogLogLin_Xs

    Table_2DLogLogLin_Xs=exp(self%xv)
    return
  end function Table_2DLogLogLin_Xs

  function Table_2DLogLogLin_Ys(self)
    !!{
    Return the $y$-values for a 2D log-log table.
    !!}
    implicit none
    class(table2DLogLogLin), intent(in   )             :: self
    double precision       , dimension(size(self%yv))  :: Table_2DLogLogLin_Ys

    Table_2DLogLogLin_Ys=exp(self%yv)
    return
  end function Table_2DLogLogLin_Ys

  function Table_2DLogLogLin_Zs(self,table)
    !!{
    Return the $y$-values for a 2D log-log table.
    !!}
    implicit none
    class           (table2DLogLogLin), intent(in   )                                    :: self
    double precision                  , dimension(size(self%xv),size(self%yv))           :: Table_2DLogLogLin_Zs
    integer                           , intent(in   )                         , optional :: table
    integer                                                                              :: tableActual

    tableActual=1
    if (present(table)) tableActual=table
    Table_2DLogLogLin_Zs=self%zv(:,:,tableActual)
    return
  end function Table_2DLogLogLin_Zs

  subroutine Table_2DLogLogLin_Populate(self,z,table)
    !!{
    Populate a 2-D log-log-linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table2DLogLogLin)                , intent(inout)           :: self
    double precision                  , dimension(:,:), intent(in   )           :: z
    integer                                           , intent(in   ), optional :: table
    integer                                                                     :: tableActual

    ! Validate the input.
    if (.not.allocated(self%zv)) call Error_Report("create the table before populating it"//{introspection:location})
    if     (                                                                                                &
         &   size(self%zv,dim=1) /= size(z,dim=1)                                                           &
         &  .or.                                                                                            &
         &   size(self%zv,dim=2) /= size(z,dim=2)                                                           &
         & ) call Error_Report("provided z array is of wrong size"//{introspection:location})
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Reset all previously stored values.
    self%tablePrevious=-1
    ! Store the z values.
    self%zv(:,:,tableActual)=z
    return
  end subroutine Table_2DLogLogLin_Populate

  subroutine Table_2DLogLogLin_Populate_Single(self,z,i,j,table)
    !!{
    Populate a single element of a 2-D log-log-linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table2DLogLogLin), intent(inout)           :: self
    double precision                  , intent(in   )           :: z
    integer                           , intent(in   )           :: i          , j
    integer                           , intent(in   ), optional :: table
    integer                                                     :: tableActual

    ! Validate the input.
    if (.not.allocated(self%zv)           ) call Error_Report("create the table before populating it"//{introspection:location})
    if (i < 1 .or. i > size(self%zv,dim=1)) call Error_Report("provided i value is out of bounds"    //{introspection:location})
    if (j < 1 .or. j > size(self%zv,dim=2)) call Error_Report("provided j value is out of bounds"    //{introspection:location})
    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Reset all previously stored values.
    self%tablePrevious=-1
    ! Store the z value.
    self%zv(i,j,tableActual)=z
    return
  end subroutine Table_2DLogLogLin_Populate_Single

  integer function Table_2DLogLogLin_Size(self,dim)
    !!{
    Return the size of a 2D log-log-linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (table2DLogLogLin), intent(in   ) :: self
    integer                  , intent(in   ) :: dim

    select case (dim)
    case (1)
       Table_2DLogLogLin_Size=self%xCount
    case (2)
       Table_2DLogLogLin_Size=self%yCount
    case default
       Table_2DLogLogLin_Size=0
       call Error_Report('1 ≤ dim ≤ 2 is required'//{introspection:location})
    end select
    return
  end function Table_2DLogLogLin_Size

  double precision function Table_2DLogLogLin_Interpolate(self,x,y,table)
    !!{
    Perform linear interpolation in a logarithmic 1D table.
    !!}
    implicit none
    class           (table2DLogLogLin), intent(inout)           :: self
    double precision                  , intent(in   )           :: x          , y
    integer                           , intent(in   ), optional :: table
    integer                                                     :: tableActual

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Test for being recalled with same values.
    if (.not.(self%xPreviousSet .and. self%yPreviousSet .and. x == self%xLinearPrevious .and. y == self%yLinearPrevious .and. tableActual == self%tablePrevious)) then
       ! Update interpolation factors.
       call self%interpolationFactors(x,y)
       ! Perform the interpolation.
       self%zPrevious=                                                                                                               &
            & +(self%zv(self%i,self%j  ,tableActual)*(1.0d0-self%hx)+self%zv(self%i+1,self%j  ,tableActual)*self%hx)*(1.0d0-self%hy) &
            & +(self%zv(self%i,self%j+1,tableActual)*(1.0d0-self%hx)+self%zv(self%i+1,self%j+1,tableActual)*self%hx)*       self%hy
    end if
    ! Return the stored value.
    Table_2DLogLogLin_Interpolate=self%zPrevious
    return
  end function Table_2DLogLogLin_Interpolate

  double precision function Table_2DLogLogLin_Interpolate_Gradient(self,x,y,dim,table)
    !!{
    Perform linear interpolation in a logarithmic 1D table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table2DLogLogLin), intent(inout)           :: self
    double precision                  , intent(in   )           :: x          , y
    integer                           , intent(in   )           :: dim
    integer                           , intent(in   ), optional :: table
    integer                                                     :: tableActual

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Test for being recalled with same values.
    if (.not.(self%xPreviousSet .and. self%yPreviousSet .and. x == self%xLinearPrevious .and. y == self%yLinearPrevious .and. tableActual == self%tablePrevious .and. dim == self%dimPrevious)) then
       ! Update interpolation factors.
       call self%interpolationFactors(x,y)
       ! Perform the interpolation.
       select case (dim)
       case (1)
          self%dzPrevious=                                                                                         &
               & +(                                                                                                &
               &   +(-self%zv(self%i,self%j  ,tableActual)+self%zv(self%i+1,self%j  ,tableActual))*(1.0d0-self%hy) &
               &   +(-self%zv(self%i,self%j+1,tableActual)+self%zv(self%i+1,self%j+1,tableActual))*       self%hy  &
               &  )                                                                                                &
               & *self%inverseDeltaX                                                                               &
               & /self%xLinearPrevious
       case (2)
          self%dzPrevious=                                                                                         &
               & +(                                                                                                &
               &   +(-self%zv(self%i  ,self%j,tableActual)+self%zv(self%i  ,self%j+1,tableActual))*(1.0d0-self%hx) &
               &   +(-self%zv(self%i+1,self%j,tableActual)+self%zv(self%i+1,self%j+1,tableActual))*       self%hx  &
               &  )                                                                                                &
               & *self%inverseDeltaY                                                                               &
               & /self%yLinearPrevious
       case default
          call Error_Report('1 ≤ dim ≤ 2 is required'//{introspection:location})
       end select
    end if
    ! Return the stored value.
    Table_2DLogLogLin_Interpolate_Gradient=self%dzPrevious
    return
  end function Table_2DLogLogLin_Interpolate_Gradient

  subroutine Table_2DLogLogLin_Interpolation_Factors(self,x,y)
    ! Update interpolation factors for x if necessary.
    implicit none
    class           (table2DLogLogLin), intent(inout)           :: self
    double precision                  , intent(in   )           :: x    , y

    ! Update interpolation factors for x if necessary.
    if (.not.self%xPreviousSet .or. x /= self%xLinearPrevious) then
       self%xPreviousSet        =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
       ! Determine the location in the table.
       if      (self%xLogarithmicPrevious <  self%xv(          1)) then
          self%i=1
       else if (self%xLogarithmicPrevious >= self%xv(self%xCount)) then
          self%i=self%xCount-1
       else
          self%i=int((self%xLogarithmicPrevious-self%xv(1))*self%inverseDeltaX)+1
       end if
       ! Compute interpolation factor.
       self%hx=(self%xLogarithmicPrevious-self%xv(self%i))*self%inverseDeltaX
    end if
    ! Update interpolation factors for y if necessary.
    if (.not.self%yPreviousSet .or. y /= self%yLinearPrevious) then
       self%yPreviousSet        =.true.
       self%yLinearPrevious     =    y
       self%yLogarithmicPrevious=log(y)
       ! Determine the location in the table.
       if      (self%yLogarithmicPrevious <  self%yv(          1)) then
          self%j=1
       else if (self%yLogarithmicPrevious >= self%yv(self%yCount)) then
          self%j=self%yCount-1
       else
          self%j=int((self%yLogarithmicPrevious-self%yv(1))*self%inverseDeltaY)+1
       end if
       ! Compute interpolation factor.
       self%hy=(self%yLogarithmicPrevious-self%yv(self%j))*self%inverseDeltaY
    end if
    return
  end subroutine Table_2DLogLogLin_Interpolation_Factors

  subroutine Table_2DLogLogLin_Destroy(self)
    !!{
    Destroy a 2D log-log-linear table.
    !!}
    implicit none
    class(table2DLogLogLin), intent(inout) :: self

    if (allocated(self%xv)) deallocate(self%xv)
    if (allocated(self%yv)) deallocate(self%yv)
    if (allocated(self%zv)) deallocate(self%zv)
    return
  end subroutine Table_2DLogLogLin_Destroy

  logical function Table_2DLogLogLin_Is_Initialized(self)
    !!{
    Return true if a 2D log-log-linear table has been created.
    !!}
    implicit none
    class(table2DLogLogLin), intent(in   ) :: self

    Table_2DLogLogLin_Is_Initialized=allocated(self%zv)
    return
  end function Table_2DLogLogLin_Is_Initialized

  subroutine Table_Linear_Monotone_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !!{
    Create a 1-D linear table.
    !!}
    use :: Error            , only : Error_Report
    use :: Numerical_Ranges , only : Make_Range                  , rangeTypeLinear
    use :: Table_Labels     , only : extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (table1DLinearMonotoneCSpline    ), intent(inout)                         :: self
    double precision                                  , intent(in   )                         :: xMaximum         , xMinimum
    integer                                           , intent(in   )                         :: xCount
    integer                                           , intent(in   ), optional               :: tableCount
    type            (enumerationextrapolationTypeType), intent(in   ), optional, dimension(2) :: extrapolationType
    integer                                                                                   :: tableCountActual

    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the x-range.
    self%xCount=xCount
    allocate(self%xv(xCount                   ))
    allocate(self%yv(xCount  ,tableCountActual))
    allocate(self%av(xCount  ,tableCountActual))
    allocate(self%bv(xCount-1,tableCountActual))
    allocate(self%cv(xCount-1,tableCountActual))
    self%xv           =Make_Range(xMinimum,xMaximum,xCount,rangeType=rangeTypeLinear)
    self%       deltaX=self%xv(2)-self%xv(1)
    self%inverseDeltaX=1.0d0/self%deltaX
    self%tablePrevious =-1
    self%dTablePrevious=-1
    self%iPrevious     =-1
    self%xPrevious     =-1.0d0
    self%dxPrevious    =-1.0d0
    ! Set extrapolation type.
    if (present(extrapolationType)) then
       if (any(extrapolationType == extrapolationTypeZero)) call Error_Report('zero extrapolation is not supported'//{introspection:location})
       self%extrapolationType=extrapolationType
    else
       self%extrapolationType=extrapolationTypeExtrapolate
    end if
    return
  end subroutine Table_Linear_Monotone_CSpline_1D_Create

  subroutine Table_Linear_Monotone_CSpline_1D_Destroy(self)
    !!{
    Destroy a linear cubic-spline 1-D table.
    !!}
    implicit none
    class(table1DLinearMonotoneCSpline), intent(inout) :: self

    call Table_1D_Destroy(self)
    if (allocated(self%av)) deallocate(self%av)
    if (allocated(self%bv)) deallocate(self%bv)
    if (allocated(self%cv)) deallocate(self%cv)
    return
  end subroutine Table_Linear_Monotone_CSpline_1D_Destroy

  subroutine Table_Linear_Monotone_CSpline_1D_Populate(self,y,table,computeSpline)
    !!{
    Populate a 1-D linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DLinearMonotoneCSpline)              , intent(inout)           :: self
    double precision                              , dimension(:), intent(in   )           :: y
    integer                                                     , intent(in   ), optional :: table
    logical                                                     , intent(in   ), optional :: computeSpline
    integer                                                                               :: tableActual
    logical                                                                               :: computeSplineActual

    ! Validate the input.
    if (.not.allocated(self%yv)       ) call Error_Report("create the table before populating it"//{introspection:location})
    if (size(self%yv,dim=1) /= size(y)) call Error_Report("provided y array is of wrong size"    //{introspection:location})

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table

    ! Store the y values.
    self%yv(:,tableActual)=y

    ! Compute the spline interpolation for this table.
    computeSplineActual=.true.
    if (present(computeSpline)) computeSplineActual=computeSpline
    if (computeSplineActual) call Table_Linear_Monotone_CSpline_1D_Compute_Spline(self,tableActual)
    return
  end subroutine Table_Linear_Monotone_CSpline_1D_Populate

  subroutine Table_Linear_Monotone_CSpline_1D_Populate_Single(self,y,i,table,computeSpline)
    !!{
    Populate a single element of a 1-D linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DLinearMonotoneCSpline), intent(inout)           :: self
    double precision                              , intent(in   )           :: y
    integer                                       , intent(in   )           :: i
    integer                                       , intent(in   ), optional :: table
    logical                                       , intent(in   ), optional :: computeSpline
    integer                                                                 :: tableActual
    logical                                                                 :: computeSplineActual

    ! Validate the input.
    if (.not.allocated(self%yv)           ) call Error_Report("create the table before populating it"//{introspection:location})
    if (i < 1 .or. i > size(self%yv,dim=1)) call Error_Report("provided i value is out of bounds"    //{introspection:location})

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table

    ! Store the y values.
    self%yv(i,tableActual)=y

    ! Compute the spline interpolation for this table.
    computeSplineActual=.true.
    if (present(computeSpline)) computeSplineActual=computeSpline
    if (computeSplineActual) call Table_Linear_Monotone_CSpline_1D_Compute_Spline(self,tableActual)
    return
  end subroutine Table_Linear_Monotone_CSpline_1D_Populate_Single

  subroutine Table_Linear_Monotone_CSpline_1D_Compute_Spline(self,table)
    !!{
    Compute the interpolating spline factors for a 1-D linear spline.
    !!}
    implicit none
    type            (table1DLinearMonotoneCSpline), intent(inout)               :: self
    integer                                       , intent(in   )               :: table
    double precision                              , allocatable  , dimension(:) :: dx   , dy    , m
    integer                                                                     :: i
    double precision                                                            :: dxSum, factor

    ! Reset all previously stored values.
    self% tablePrevious=-1
    self%dTablePrevious=-1
    self%     iPrevious=-1
    ! Allocate workspace.
    allocate(dx(size(self%xv)-1))
    allocate(dy(size(self%xv)-1))
    allocate( m(size(self%xv)-1))
    ! Get consecutive differences and slopes.
    do i=1,size(self%xv)-1
       dx(i)=self%xv(i+1      )-self%xv(i      )
       dy(i)=self%yv(i+1,table)-self%yv(i,table)
       m (i)=dy(i)/dx(i)
    end do
    ! Get degree-1 coefficients.
    self%av(1,table)=m(1)
    do i=1,size(dx)-1
       if (m(i)*m(i+1) <= 0.0d0) then
          self%av(i+1,table)=0.0d0
       else
          dxSum=dx(i)+dx(i+1)
          self%av(i+1,table)=3.0d0*dxSum/((dxSum+dx(i+1))/m(i)+(dxSum+dx(i))/m(i+1))
       end if
    end do
    self%av(size(dx)+1,table)=m(size(m))
    ! Get degree-2 and degree-3 coefficients.
    do i=1,size(self%av)-1
       factor=self%av(i,table)+self%av(i+1,table)-2.0d0*m(i)
       self%bv(i,table)=(m(i)-self%av(i,table)-factor)/dx(i)
       self%cv(i,table)=factor/dx(i)**2
    end do
    ! Destroy workspace.
    deallocate(dx)
    deallocate(dy)
    deallocate( m)
    return
  end subroutine Table_Linear_Monotone_CSpline_1D_Compute_Spline

  double precision function Table_Linear_Monotone_CSpline_1D_Interpolate(self,x,table,status)
    !!{
    Perform linear interpolation in a linear 1D table.
    !!}
    implicit none
    class           (table1DLinearMonotoneCSpline), intent(inout)           :: self
    double precision                              , intent(in   )           :: x
    integer                                       , intent(in   ), optional :: table
    integer                                       , intent(  out), optional :: status
    integer                                                                 :: i     , tableActual
    double precision                                                        :: dx    , xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%xPrevious .or. tableActual /= self%tablePrevious) then
       xEffective=self%xEffective(x,status)
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=int((xEffective-self%xv(1))*self%inverseDeltaX)+1
       end if
       ! Compute offset from tabulated point.
       dx=xEffective-self%xv(i)
       ! Interpolate in the table.
       self%xPrevious    =x
       self%tablePrevious=tableActual
       self%    yPrevious=self%yv(i,tableActual)+self%av(i,tableActual)*dx+self%bv(i,tableActual)*dx**2+self%cv(i,tableActual)*dx**3
    end if
    Table_Linear_Monotone_CSpline_1D_Interpolate=self%yPrevious
    return
  end function Table_Linear_Monotone_CSpline_1D_Interpolate

  double precision function Table_Linear_Monotone_CSpline_1D_Interpolate_Gradient(self,x,table,status)
    !!{
    Perform linear interpolation in a linear 1D table.
    !!}
    implicit none
    class           (table1DLinearMonotoneCSpline), intent(inout)           :: self
    double precision                              , intent(in   )           :: x
    integer                                       , intent(in   ), optional :: table
    integer                                       , intent(  out), optional :: status
    integer                                                                 :: i     , tableActual
    double precision                                                        :: dx    , xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then
       xEffective=self%xEffective(x,status)
       ! Determine the location in the table.
       if      (xEffective <  self%xv(          1)) then
          i=1
       else if (xEffective >= self%xv(self%xCount)) then
          i=self%xCount-1
       else
          i=int((xEffective-self%xv(1))*self%inverseDeltaX)+1
       end if
       ! Compute offset from tabulated point.
       dx=xEffective-self%xv(i)
       ! Interpolate in the table.
       self%dxPrevious    =x
       self%dTablePrevious=tableActual
       self%    dyPrevious=self%av(i,tableActual)+2.0d0*self%bv(i,tableActual)*dx+3.0d0*self%cv(i,tableActual)*dx**2
    end if
    Table_Linear_Monotone_CSpline_1D_Interpolate_Gradient=self%dyPrevious
    return
  end function Table_Linear_Monotone_CSpline_1D_Interpolate_Gradient

  function Table_Linear_Monotone_CSpline_Integration_Weights(self,x0,x1,integrand)
    !!{
    Returns a set of weights for trapezoidal integration on the table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table1DLinearMonotoneCSpline   ), intent(inout)                               :: self
    double precision                                 , intent(in   )                               :: x0, x1
    procedure       (tablesIntegrationWeightFunction), intent(in   )           , pointer, optional :: integrand
    double precision                                 , dimension(size(self%xv))                    :: Table_Linear_Monotone_CSpline_Integration_Weights
    !$GLC attributes unused :: self, x0, x1, integrand

    Table_Linear_Monotone_CSpline_Integration_Weights=0.0d0
    call Error_Report('integration weights not supported'//{introspection:location})
    return
  end function Table_Linear_Monotone_CSpline_Integration_Weights

  subroutine Table_Logarithmic_Monotone_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !!{
    Create a 1-D logarithmic table.
    !!}
    implicit none
    class           (table1DLogarithmicMonotoneCSpline), intent(inout)                         :: self
    double precision                                   , intent(in   )                         :: xMaximum         , xMinimum
    integer                                            , intent(in   )                         :: xCount
    type            (enumerationExtrapolationTypeType ), intent(in   ), optional, dimension(2) :: extrapolationType
    integer                                            , intent(in   ), optional               ::  tableCount

    self%previousSet         =.false.
    self%xLinearPrevious     =-1.0d0
    self%xLogarithmicPrevious=-1.0d0
    ! Call the creator for linear tables with the logarithms of the input x range.
    call self%table1DLinearMonotoneCSpline%create(log(xMinimum),log(xMaximum),xCount,tableCount,extrapolationType)
    ! Store the minimum and maximum x-values for rapid look-up.
    self%xMinimum=exp(self%xv(     1))
    self%xMaximum=exp(self%xv(xCount))
    return
  end subroutine Table_Logarithmic_Monotone_CSpline_1D_Create

  double precision function Table_Logarithmic_Monotone_CSpline_1D_X(self,i)
    !!{
    Return the {\normalfont \ttfamily i}$^\mathrm{th}$ $x$-value for a logarithmic 1D table.
    !!}
    implicit none
    class  (table1DLogarithmicMonotoneCSpline), intent(inout) :: self
    integer                                   , intent(in   ) :: i

    ! Check for end-points, and return stored values if possible.
    if      (i ==  1                      ) then
       Table_Logarithmic_Monotone_CSpline_1D_X=self%xMinimum
    else if (i == -1 .or. i == self%xCount) then
       Table_Logarithmic_Monotone_CSpline_1D_X=self%xMaximum
    else
       ! No stored value is available - simply look up the required value.
       Table_Logarithmic_Monotone_CSpline_1D_X=exp(self%table1DLinearMonotoneCSpline%x(i))
    end if
    return
  end function Table_Logarithmic_Monotone_CSpline_1D_X

  function Table_Logarithmic_Monotone_CSpline_1D_Xs(self)
    !!{
    Return the $x$-values for a 1D table.
    !!}
    implicit none
    class           (table1DLogarithmicMonotoneCSpline), intent(in   )            :: self
    double precision                                   , dimension(size(self%xv)) :: Table_Logarithmic_Monotone_CSpline_1D_Xs

    Table_Logarithmic_Monotone_CSpline_1D_Xs=exp(self%table1DLinearMonotoneCSpline%xs())
    return
  end function Table_Logarithmic_Monotone_CSpline_1D_Xs

  double precision function Table_Logarithmic_Monotone_CSpline_1D_Interpolate(self,x,table,status)
    !!{
    Perform linear interpolation in a logarithmic 1D table.
    !!}
    implicit none
    class           (table1DLogarithmicMonotoneCSpline), intent(inout)           :: self
    double precision                                   , intent(in   )           :: x
    integer                                            , intent(in   ), optional :: table
    integer                                            , intent(  out), optional :: status

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_Monotone_CSpline_1D_Interpolate=self%table1DLinearMonotoneCSpline%interpolate(self%xLogarithmicPrevious,table,status)
    return
  end function Table_Logarithmic_Monotone_CSpline_1D_Interpolate

  double precision function Table_Logarithmic_Monotone_CSpline_1D_Interpolate_Gradient(self,x,table,status)
    !!{
    Perform linear interpolation in a logarithmic 1D table.
    !!}
    implicit none
    class           (table1DLogarithmicMonotoneCSpline), intent(inout)           :: self
    double precision                                   , intent(in   )           :: x
    integer                                            , intent(in   ), optional :: table
    integer                                            , intent(  out), optional :: status

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_Monotone_CSpline_1D_Interpolate_Gradient=self%table1DLinearMonotoneCSpline%interpolateGradient(self%xLogarithmicPrevious,table,status)/self%xEffective(x,status)
    return
  end function Table_Logarithmic_Monotone_CSpline_1D_Interpolate_Gradient

  subroutine Table_2D_LinLinLin_Create(self,x,y,tableCount)
    !!{
    Create a 2-D generic table.
    !!}
    implicit none
    class           (table2DLinLinLin)              , intent(inout)           :: self
    double precision                  , dimension(:), intent(in   )           :: x         , y
    integer                                         , intent(in   ), optional :: tableCount
    !![
    <optionalArgument name="tableCount" defaultsTo="1"/>
    !!]

    ! Allocate arrays and construct the x-range.
    self%xCount=size(x)
    self%yCount=size(y)
    allocate(self%xv(size(x)                    ))
    allocate(self%yv(        size(y)            ))
    allocate(self%zv(size(x),size(y),tableCount_))
    self%xv=x
    self%yv=y
    ! Build the interpolators.
    self%interpolatorX=interpolator(self%xv)
    self%interpolatorY=interpolator(self%yv)
    return
  end subroutine Table_2D_LinLinLin_Create

  subroutine Table_2D_LinLinLin_Destroy(self)
    !!{
    Destroy a generic 2-D table.
    !!}
    implicit none
    class(table2DLinLinLin), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine Table_2D_LinLinLin_Destroy

  subroutine Table_2D_LinLinLin_Populate(self,z,table)
    !!{
    Populate a 2-D linear table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table2DLinLinLin)                , intent(inout)           :: self
    double precision                  , dimension(:,:), intent(in   )           :: z
    integer                                           , intent(in   ), optional :: table
    !![
    <optionalArgument name="table" defaultsTo="1"/>
    !!]

    ! Validate the input.
    if (.not.allocated(self%zv)             ) call Error_Report("create the table before populating it"//{introspection:location})
    if (size(self%zv,dim=1) /= size(z,dim=1)) call Error_Report("provided z array is of wrong size"    //{introspection:location})
    if (size(self%zv,dim=2) /= size(z,dim=2)) call Error_Report("provided z array is of wrong size"    //{introspection:location})

    ! Store the y values.
    self%zv(:,:,table_)=z
    return
  end subroutine Table_2D_LinLinLin_Populate

  subroutine Table_2D_LinLinLin_Populate_Single(self,z,i,j,table)
    !!{
    Populate a single element of a 2-D generic table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (table2DLinLinLin), intent(inout)           :: self
    double precision                  , intent(in   )           :: z
    integer                           , intent(in   )           :: i    , j
    integer                           , intent(in   ), optional :: table
    !![
    <optionalArgument name="table" defaultsTo="1"/>
    !!]

    ! Validate the input.
    if (.not.allocated(self%zv)           ) call Error_Report("create the table before populating it"//{introspection:location})
    if (i < 1 .or. i > size(self%zv,dim=1)) call Error_Report("provided i value is out of bounds"    //{introspection:location})
    if (j < 1 .or. j > size(self%zv,dim=2)) call Error_Report("provided j value is out of bounds"    //{introspection:location})

    ! Store the y values.
    self%zv(i,j,table_)=z
    return
  end subroutine Table_2D_LinLinLin_Populate_Single

  double precision function Table_2D_LinLinLin_Interpolate(self,x,y,table)
    !!{
    Perform generic interpolation in a generic 2D table.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    implicit none
    class           (table2DLinLinLin), intent(inout)            :: self
    double precision                  , intent(in   )            :: x    , y
    integer                           , intent(in   ) , optional :: table
    integer         (c_size_t        )                           :: i    , j , &
         &                                                          ii   , jj
    double precision                  , dimension(0:1)           :: hi   , hj
    !![
    <optionalArgument name="table" defaultsTo="1"/>
    !!]

    ! Compute interpolating factors.
    call self%interpolatorX%linearFactors(x,i,hi)
    call self%interpolatorY%linearFactors(y,j,hj)
    ! Perform the interpolation.
    Table_2D_LinLinLin_Interpolate=0.0d0
    do ii=0,1
       do jj=0,1
          Table_2D_LinLinLin_Interpolate=Table_2D_LinLinLin_Interpolate+self%zv(i+ii,j+jj,table_)*hi(ii)*hj(jj)
       end do
    end do
    return
  end function Table_2D_LinLinLin_Interpolate

  function Table_2D_LinLinLin_Xs(self)
    !!{
    Return the $x$-values for a 2D table.
    !!}
    implicit none
    class           (table2DLinLinLin), intent(in   )            :: self
    double precision                  , dimension(size(self%xv)) :: Table_2D_LinLinLin_Xs

    Table_2D_LinLinLin_Xs=self%xv
    return
  end function Table_2D_LinLinLin_Xs

  function Table_2D_LinLinLin_Ys(self)
    !!{
    Return the $y$-values for a 2D table.
    !!}
    implicit none
    class           (table2DLinLinLin), intent(in   )            :: self
    double precision                  , dimension(size(self%yv)) :: Table_2D_LinLinLin_Ys

    Table_2D_LinLinLin_Ys=self%yv
    return
  end function Table_2D_LinLinLin_Ys

  function Table_2D_LinLinLin_Zs(self)
    !!{
    Return the $z$-values for a 2D table.
    !!}
    implicit none
    class           (table2DLinLinLin), intent(in   )                                                          :: self
    double precision                  , dimension(size(self%zv,dim=1),size(self%zv,dim=2),size(self%zv,dim=3)) :: Table_2D_LinLinLin_Zs

    Table_2D_LinLinLin_Zs=self%zv
    return
  end function Table_2D_LinLinLin_Zs

  subroutine Table_2D_LinLinLin_Interpolator_Reinitialize(self)
    !!{
    Reinitialize the interpolator.
    !!}
    implicit none
    class(table2DLinLinLin), intent(inout) :: self

    call self%interpolatorX%GSLReallocate()
    call self%interpolatorY%GSLReallocate()
    return
  end subroutine Table_2D_LinLinLin_Interpolator_Reinitialize

end module Tables
