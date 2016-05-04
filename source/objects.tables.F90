!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which defines a {\normalfont \ttfamily table} class with optimized interpolation operators.

module Tables
  !% Defines a {\normalfont \ttfamily table} class with optimized interpolation operators.
  use FGSL
  use Table_Labels
  private
  public :: table                       , table1D                          , table1DGeneric                    , &
       &    table1DLinearLinear         , table1DLogarithmicLinear         , table1DNonUniformLinearLogarithmic, &
       &    table1DLinearCSpline        , table1DLogarithmicCSpline        , table2DLogLogLin                  , &
       &    table1DLinearMonotoneCSpline, table1DLogarithmicMonotoneCSpline

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
       !% Interface to {\normalfont \ttfamily table} destructor.
       import table
       implicit none
       class(table), intent(inout) :: self
     end subroutine Table_Destroy
  end interface

  type, abstract, extends(table) :: table1D
     !% Basic table type.
     integer                                       :: xCount
     integer                      , dimension(2  ) :: extrapolationType
     double precision, allocatable, dimension(:  ) :: xv
     double precision, allocatable, dimension(:,:) :: yv
   contains
     !@ <objectMethods>
     !@   <object>table1D</object>
     !@   <objectMethod>
     !@     <method>interpolate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ x,\intzero\ [table]</arguments>
     !@     <description>Interpolate to {\normalfont \ttfamily x} in the {\normalfont \ttfamily table}$^{\mathrm th}$ table.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>interpolateGradient</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ x,\intzero\ [table]</arguments>
     !@     <description>Interpolate the gradient to {\normalfont \ttfamily x} in the {\normalfont \ttfamily table}$^{\mathrm th}$ table.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reverse</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(table)\textgreater} reversedSelf,\intzero\ [table], \logicalzero\ [precise]</arguments>
     !@     <description>Reverse the table (i.e. swap $x$ and $y$ components) and return in {\normalfont \ttfamily reversedSelf}. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^{\mathrm th}$ table is used for the $y$-values, otherwise the first table is used. If the optional {\normalfont \ttfamily precise} argument is set to {\normalfont \ttfamily true} then the reversal must be precisely invertible---if this is not possible the method will abort.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isMonotonic</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\enum\ [directionDecreasing|directionIncreasing],\logicalzero\ [allowEqual],\intzero\ [table]</arguments>
     !@     <description>Return true if the table $y$-values are monotonic. Optionally, the direction of monotonicity can be specified via the {\normalfont \ttfamily direction} argument---by default either direction is allowed. By default consecutive equal values are considered non-monotonic. This behavior can be changed via the optional {\normalfont \ttfamily allowEqual} argument. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^{\mathrm th}$ table is used for the $y$-values, otherwise the first table is used.</description>
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
     !@     <description>Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $x$-value.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>y</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\intzero\ i,\intzero\ [table]</arguments>
     !@     <description>Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $y$-value. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^{\mathrm th}$ table is used for the $y$-values, otherwise the first table is used.</description>
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
     !@     <description>Return an array of all $y$-values. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^{\mathrm th}$ table is used for the $y$-values, otherwise the first table is used.</description>
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
     !@     <description>Return the weights to be applied to the table to integrate (using the trapezium rule) between {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.</description>
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
       !% Interface to {\normalfont \ttfamily table} interpolator.
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
     type   (fgsl_interp_type ) :: interpolationType
     logical                    :: reset            =.true.
   contains
     !@ <objectMethods>
     !@   <object>table1DGeneric</object>
     !@   <objectMethod>
     !@     <method>create</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ x,\intzero\ [tableCount]</arguments>
     !@     <description>Create the object with the specified {\normalfont \ttfamily x} values, and with {\normalfont \ttfamily tableCount} tables.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>populate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero|\doubleone\ y,\intzero\ [i],\intzero\ [table]</arguments>
     !@     <description>Populate the {\normalfont \ttfamily table}$^{\mathrm th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified.</description>
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
     !@     <description>Create the object with $x$-values spanning the range {\normalfont \ttfamily xMinimum} to {\normalfont \ttfamily xMaximum} in {\normalfont \ttfamily xCount} steps, and with {\normalfont \ttfamily tableCount} tables.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>populate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero|\doubleone\ y,\intzero\ [i],\intzero\ [table]</arguments>
     !@     <description>Populate the {\normalfont \ttfamily table}$^{\mathrm th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified.</description>
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
     double precision, allocatable, dimension(:,:) :: sv            , av           , bv           , cv       , &
          &                                           dv
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
     !@     <description>Create the object with $x$-values spanning the range {\normalfont \ttfamily xMinimum} to {\normalfont \ttfamily xMaximum} in {\normalfont \ttfamily xCount} steps, and with {\normalfont \ttfamily tableCount} tables.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>populate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero|\doubleone\ y,\intzero\ [i],\intzero\ [table]</arguments>
     !@     <description>Populate the {\normalfont \ttfamily table}$^{\mathrm th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified.</description>
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
     double precision :: xLinearPrevious, xLogarithmicPrevious
     double precision :: xMinimum       , xMaximum
   contains
     procedure :: create             =>Table_Logarithmic_CSpline_1D_Create
     procedure :: interpolate        =>Table_Logarithmic_CSpline_1D_Interpolate
     procedure :: interpolateGradient=>Table_Logarithmic_CSpline_1D_Interpolate_Gradient
     procedure :: x                  =>Table_Logarithmic_CSpline_1D_X
     procedure :: xs                 =>Table_Logarithmic_CSpline_1D_Xs
  end type table1DLogarithmicCSpline

  type, extends(table1D) :: table1DLinearMonotoneCSpline
     !% Table type supporting one dimensional table with linear spacing in $x$ and monotonic cubic spline interpolation.
     double precision, allocatable, dimension(:,:) :: c1            , c2           , c3
     integer                                       :: dTablePrevious, iPrevious    , tablePrevious
     double precision                              :: deltaX        , inverseDeltaX
     double precision                              :: dxPrevious    , dyPrevious   , xPrevious    , yPrevious
   contains
     !@ <objectMethods>
     !@   <object>table1DLinearMonotoneCSpline</object>
     !@   <objectMethod>
     !@     <method>create</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ xMinimum,\doublezero\ xMaximum,\intzero xCount,\intzero [tableCount]</arguments>
     !@     <description>Create the object with $x$-values spanning the range {\normalfont \ttfamily xMinimum} to {\normalfont \ttfamily xMaximum} in {\normalfont \ttfamily xCount} steps, and with {\normalfont \ttfamily tableCount} tables.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>populate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero|\doubleone\ y,\intzero\ [i],\intzero\ [table]</arguments>
     !@     <description>Populate the {\normalfont \ttfamily table}$^{\mathrm th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: create                                 =>Table_Linear_Monotone_CSpline_1D_Create
     procedure :: destroy                                =>Table_Linear_Monotone_CSpline_1D_Destroy
     procedure :: Table_Linear_Monotone_CSpline_1D_Populate
     procedure :: Table_Linear_Monotone_CSpline_1D_Populate_Single
     generic   :: populate                                => Table_Linear_Monotone_CSpline_1D_Populate            , &
          &                                                  Table_Linear_Monotone_CSpline_1D_Populate_Single
     procedure :: interpolate        =>Table_Linear_Monotone_CSpline_1D_Interpolate
     procedure :: interpolateGradient=>Table_Linear_Monotone_CSpline_1D_Interpolate_Gradient
     procedure :: integrationWeights =>Table_Linear_Monotone_CSpline_Integration_Weights
  end type table1DLinearMonotoneCSpline
  
  type, extends(table1DLinearMonotoneCSpline) :: table1DLogarithmicMonotoneCSpline
     !% Table type supporting one dimensional table with logarithmic spacing in $x$ and monotonic cubic spline interpolation.
     logical          :: previousSet
     double precision :: xLinearPrevious, xLogarithmicPrevious
     double precision :: xMinimum       , xMaximum
   contains
     procedure :: create             =>Table_Logarithmic_Monotone_CSpline_1D_Create
     procedure :: interpolate        =>Table_Logarithmic_Monotone_CSpline_1D_Interpolate
     procedure :: interpolateGradient=>Table_Logarithmic_Monotone_CSpline_1D_Interpolate_Gradient
     procedure :: x                  =>Table_Logarithmic_Monotone_CSpline_1D_X
     procedure :: xs                 =>Table_Logarithmic_Monotone_CSpline_1D_Xs
  end type table1DLogarithmicMonotoneCSpline

  abstract interface
     double precision function integrandTemplate(x)
       double precision, intent(in   ) :: x
     end function integrandTemplate
  end interface

  type, extends(table) :: table2DLogLogLin
     !% Two-dimensional table type with logarithmic spacing in x and y dimensions, and linear interpolation in z.
     integer                                         :: extrapolationTypeX  , extrapolationTypeY  , &
          &                                             xCount              , yCount              , &
          &                                             i                   , j                   , &
          &                                             tablePrevious       , dimPrevious
     double precision                                :: xLinearPrevious     , yLinearPrevious     , &
          &                                             xLogarithmicPrevious, yLogarithmicPrevious, &
          &                                             hx                  , hy                  , &
          &                                             inverseDeltaX       , inverseDeltaY       , &
          &                                             zPrevious           , dzPrevious
     logical                                         :: xPreviousSet        , yPreviousSet
     double precision, allocatable, dimension(:    ) :: xv                  , yv
     double precision, allocatable, dimension(:,:,:) :: zv
   contains
     !@ <objectMethods>
     !@   <object>table2DLogLogLin</object>
     !@   <objectMethod>
     !@     <method>interpolationFactors</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ x, \doublezero\ y</arguments>
     !@     <description>Compute and store interpolation factors to {\normalfont \ttfamily (x,y)}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>interpolate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ x,\intzero\ [table]</arguments>
     !@     <description>Interpolate to {\normalfont \ttfamily x} in the {\normalfont \ttfamily table}$^{\mathrm th}$ table.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>interpolateGradient</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ x,\intzero\ [table]</arguments>
     !@     <description>Interpolate the gradient to {\normalfont \ttfamily x} in the {\normalfont \ttfamily table}$^{\mathrm th}$ table.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>size</method>
     !@     <type>\intzero</type>
     !@     <arguments>\intzero\ dim</arguments>
     !@     <description>Return the size (i.e. number of $x$ or $y$-values) in the table of the given dimension.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>x</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\intzero\ i</arguments>
     !@     <description>Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $x$-value.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>y</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\intzero\ i,\intzero\ [table]</arguments>
     !@     <description>Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $y$-value. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^{\mathrm th}$ table is used for the $y$-values, otherwise the first table is used.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>z</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\intzero\ i,\intzero\ j,\intzero\ [table]</arguments>
     !@     <description>Return the {\normalfont \ttfamily (i,j)}$^{\mathrm th}$ $z$-value. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^{\mathrm th}$ table is used for the $z$-values, otherwise the first table is used.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>xs</method>
     !@     <type>\doubleone</type>
     !@     <arguments></arguments>
     !@     <description>Return an array of all $x$-values.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>ys</method>
     !@     <type>\doubleone</type>
     !@     <arguments></arguments>
     !@     <description>Return an array of all $y$-values. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^{\mathrm th}$ table is used for the $y$-values, otherwise the first table is used.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>zs</method>
     !@     <type>\doubletwo</type>
     !@     <arguments>\intzero\ [table]</arguments>
     !@     <description>Return an array of all $z$-values. If {\normalfont \ttfamily table} is specified then the {\normalfont \ttfamily table}$^{\mathrm th}$ table is used for the $z$-values, otherwise the first table is used.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isInitialized</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if the table is initialized (this means the table is created, it may not yet have been populated).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>populate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero|\doubletwo\ z,\intzero\ [i],\intzero\ [j],\intzero\ [table]</arguments>
     !@     <description>Populate the {\normalfont \ttfamily table}$^{\mathrm th}$ table with elements {\normalfont \ttfamily y}. If {\normalfont \ttfamily y} is a scalar, then the index, {\normalfont \ttfamily i}, of the element to set must also be specified.</description>
     !@   </objectMethod>  
     !@   <objectMethod>
     !@     <method>create</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ xMinimum,\doublezero\ xMaximum,\doublezero\ yMinimum,\doublezero\ yMaximum,\intzero yCount,\intzero [tableCount],\enumExtrapolationType [extrapolationTypeX],\enumExtrapolationType [extrapolationTypeY]</arguments>
     !@     <description>Create the object with $x$-values spanning the range {\normalfont \ttfamily xMinimum} to {\normalfont \ttfamily xMaximum} in {\normalfont \ttfamily xCount} steps, and with {\normalfont \ttfamily tableCount} tables.</description>
     !@   </objectMethod>
     !@ </objectMethods>
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
    !% Destroy a 1-D table.
    use Memory_Management
    implicit none
    class(table1D), intent(inout) :: self

    if (allocated(self%xv)) call Dealloc_Array(self%xv)
    if (allocated(self%yv)) call Dealloc_Array(self%yv)
    return
  end subroutine Table_1D_Destroy

  double precision function Table1D_X(self,i)
    !% Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $x$-value for a 1D table.
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
    !% Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $y$-value for a 1D table.
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
       ! Set linear interpolation.
       reversedSelf%interpolationType=FGSL_Interp_Linear
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
    !% Returns a set of weights for trapezoidal integration on the table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
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

  subroutine Table_Generic_1D_Create(self,x,tableCount,extrapolationType,interpolationType)
    !% Create a 1-D generic table.
    use Memory_Management
    use Galacticus_Error
    implicit none
    class           (table1DGeneric  )              , intent(inout)           :: self
    double precision                  , dimension(:), intent(in   )           :: x
    integer                                         , intent(in   ), optional :: tableCount
    integer                           , dimension(2), intent(in   ), optional :: extrapolationType
    type            (fgsl_interp_type)              , intent(in   ), optional :: interpolationType
    integer                                                                   :: tableCountActual

    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the x-range.
    self%xCount=size(x)
    call Alloc_Array(self%xv,[size(x)                 ])
    call Alloc_Array(self%yv,[size(x),tableCountActual])
    self%xv   =x
    self%reset=.true.
    ! Set interpoaltion type.
    if (present(interpolationType)) then
       self%interpolationType=interpolationType
    else
       self%interpolationType=FGSL_Interp_Linear
    end if
    ! Set extrapolation type.
    if (present(extrapolationType)) then
       if (any(extrapolationType == extrapolationTypeZero)) call Galacticus_Error_Report('Table_Generic_1D_Create','zero extrapolation is not supported')
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
    integer                                                   :: tableActual

    tableActual=1
    if (present(table)) tableActual=table
    Table_Generic_1D_Interpolate=Interpolate(self%xv,self%yv(:,tableActual),self%interpolator,self%accelerator,self%xEffective(x),extrapolationType=self%extrapolationType(1),interpolationType=self%interpolationType,reset=self%reset)
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
    Table_Generic_1D_Interpolate_Gradient=Interpolate_Derivative(self%xv,self%yv(:,tableActual),self%interpolator,self%accelerator,self%xEffective(x),reset=self%reset,extrapolationType=self%extrapolationType(1),interpolationType=self%interpolationType)
    return
  end function Table_Generic_1D_Interpolate_Gradient

  subroutine Table_Linear_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !% Create a 1-D linear table.
    use Memory_Management
    use Numerical_Ranges
    use Galacticus_Error
    implicit none
    class           (table1DLinearLinear), intent(inout)                         :: self
    double precision                     , intent(in   )                         :: xMaximum         , xMinimum
    integer                              , intent(in   )                         :: xCount
    integer                              , intent(in   ), optional               :: tableCount
    integer                              , intent(in   ), optional, dimension(2) :: extrapolationType
    integer                                                                      :: tableCountActual

    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the x-range.
    self%xCount=xCount
    call Alloc_Array(self%xv,[xCount                 ])
    call Alloc_Array(self%yv,[xCount,tableCountActual])
    self%xv            =Make_Range(xMinimum,xMaximum,xCount,rangeType=rangeTypeLinear)
    self%inverseDeltaX =1.0d0/(self%xv(2)-self%xv(1))
    self%tablePrevious =-1
    self%dTablePrevious=-1
    self%xPrevious     =-1.0d0
    self%dxPrevious    =-1.0d0
    ! Set extrapolation type.
    if (present(extrapolationType)) then
       if (any(extrapolationType == extrapolationTypeZero)) call Galacticus_Error_Report('Table_Linear_1D_Create','zero extrapolation is not supported')
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
    class           (table1DLogarithmicLinear), intent(inout)                         :: self
    double precision                          , intent(in   )                         :: xMaximum         , xMinimum
    integer                                   , intent(in   )                         :: xCount
    integer                                   , intent(in   ), optional               :: tableCount
    integer                                   , intent(in   ), optional, dimension(2) :: extrapolationType

    self%previousSet         =.false.
    self%xLinearPrevious     =-1.0d0
    self%xLogarithmicPrevious=-1.0d0
    ! Call the creator for linear tables with the logarithms of the input x range.
    call self%table1DLinearLinear%create(log(xMinimum),log(xMaximum),xCount,tableCount,extrapolationType)
    return
  end subroutine Table_Logarithmic_1D_Create

  double precision function Table_Logarithmic_1D_X(self,i)
    !% Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $x$-value for a logarithmic 1D table.
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
    !% Returns a set of weights for trapezoidal integration on the table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
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
    !GCC$ attributes unused :: precise
    
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
    use Galacticus_Error
    implicit none
    class           (table1DLinearCSpline), intent(inout)                         :: self
    double precision                      , intent(in   )                         :: xMaximum         , xMinimum
    integer                               , intent(in   )                         :: xCount
    integer                               , intent(in   ), optional               :: tableCount
    integer                               , intent(in   ), optional, dimension(2) :: extrapolationType
    integer                                                                       :: tableCountActual

    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the x-range.
    self%xCount=xCount
    call Alloc_Array(self%xv,[xCount                   ])
    call Alloc_Array(self%yv,[xCount  ,tableCountActual])
    call Alloc_Array(self%sv,[xCount  ,tableCountActual])
    call Alloc_Array(self%av,[xCount-1,tableCountActual])
    call Alloc_Array(self%bv,[xCount-1,tableCountActual])
    call Alloc_Array(self%cv,[xCount-1,tableCountActual])
    call Alloc_Array(self%dv,[xCount-1,tableCountActual])
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
       if (any(extrapolationType == extrapolationTypeZero)) call Galacticus_Error_Report('Table_Linear_CSpline_1D_Create','zero extrapolation is not supported')
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
    if (allocated(self%av)) call Dealloc_Array(self%av)
    if (allocated(self%bv)) call Dealloc_Array(self%bv)
    if (allocated(self%cv)) call Dealloc_Array(self%cv)
    if (allocated(self%dv)) call Dealloc_Array(self%dv)
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

  double precision function Table_Linear_CSpline_1D_Interpolate(self,x,table)
    !% Perform linear interpolation in a linear 1D table.
    implicit none
    class           (table1DLinearCSpline), intent(inout)           :: self
    double precision                      , intent(in   )           :: x
    integer                               , intent(in   ), optional :: table
    integer                                                         :: i         , tableActual
    double precision                                                :: xEffective, dx

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%xPrevious .or. tableActual /= self%tablePrevious) then
       xEffective=self%xEffective(x)
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
       self%    yPrevious=self%av(i,tableActual)+dx*(self%bv(i,tableActual)+dx*(self%cv(i,tableActual)+dx*self%dv(i,tableActual)))
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
    double precision                                                :: xEffective, dx

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then
       xEffective=self%xEffective(x)
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
       self%    dyPrevious=self%bv(i,tableActual)+dx*(2.0d0*self%cv(i,tableActual)+dx*3.0d0*self%dv(i,tableActual))
    end if
    Table_Linear_CSpline_1D_Interpolate_Gradient=self%dyPrevious
    return
  end function Table_Linear_CSpline_1D_Interpolate_Gradient

  subroutine Table_Logarithmic_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !% Create a 1-D logarithmic table.
    implicit none
    class           (table1DLogarithmicCSpline), intent(inout)                         :: self
    double precision                           , intent(in   )                         :: xMaximum         , xMinimum
    integer                                    , intent(in   )                         :: xCount
    integer                                    , intent(in   ), optional               :: tableCount
    integer                                    , intent(in   ), optional, dimension(2) :: extrapolationType

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
    !% Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $x$-value for a logarithmic 1D table.
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
    !% Returns a set of weights for trapezoidal integration on the table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
    use Galacticus_Error
    implicit none
    class           (table1DLinearCSpline), intent(inout)                               :: self
    double precision                      , intent(in   )                               :: x0, x1
    procedure       (integrandTemplate   ), intent(in   )           , pointer, optional :: integrand
    double precision                      , dimension(size(self%xv))                    :: Table_Linear_CSpline_Integration_Weights
    !GCC$ attributes unused :: self, x0, x1

    Table_Linear_CSpline_Integration_Weights=0.0d0
    call Galacticus_Error_Report('Table_Linear_CSpline_Integration_Weights','integration weights not supported')
    return
  end function Table_Linear_CSpline_Integration_Weights
  
  double precision function Table1D_Find_Effective_X(self,x)
    !% Return the effective value of $x$ to use in table interpolations.
    use Galacticus_Error
    implicit none
    class           (table1D), intent(inout) :: self
    double precision         , intent(in   ) :: x

    if      (x < self%x(+1)) then
       select case (self%extrapolationType(1))
       case (extrapolationTypeExtrapolate)
          Table1D_Find_Effective_X=x
       case (extrapolationTypeFix        )
          Table1D_Find_Effective_X=self%x(+1)
       case default
          call Galacticus_Error_Report('Table_1D_Find_Effective_X','x is below range')
       end select
    else if (x > self%x(-1)) then
       select case (self%extrapolationType(2))
       case (extrapolationTypeExtrapolate)
          Table1D_Find_Effective_X=x
       case (extrapolationTypeFix        )
          Table1D_Find_Effective_X=self%x(-1)
       case default
          call Galacticus_Error_Report('Table_1D_Find_Effective_X','x is above range')
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
    !% Returns a set of weights for integration on a linear-logarithmic table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
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
    !% Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $y$-value for a 1D table.
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
  
  subroutine Table_2DLogLogLin_Create(self,xMinimum,xMaximum,xCount,yMinimum,yMaximum,yCount,tableCount,extrapolationTypeX,extrapolationTypeY)
    !% Create a 2-D log-log-linear table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class           (table2DLogLogLin), intent(inout)           :: self
    double precision                  , intent(in   )           :: xMaximum          , xMinimum          , &
         &                                                         yMaximum          , yMinimum
    integer                           , intent(in   )           :: xCount            , yCount
    integer                           , intent(in   ), optional :: extrapolationTypeX, extrapolationTypeY, &
         &                                                         tableCount
    integer                                                     :: tableCountActual
    
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
    call Alloc_Array(self%xv,[xCount                        ])
    call Alloc_Array(self%yv,[       yCount                 ])
    call Alloc_Array(self%zv,[xCount,yCount,tableCountActual])
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
    !% Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $x$-value for a 2D log-log table.
    implicit none
    class  (table2DLogLogLin), intent(inout) :: self
    integer                  , intent(in   ) :: i

    Table_2DLogLogLin_X=exp(self%xv(i))
    return
  end function Table_2DLogLogLin_X

  double precision function Table_2DLogLogLin_Y(self,i)
    !% Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $y$-value for a 2D log-log table.
    implicit none
    class  (table2DLogLogLin), intent(inout) :: self
    integer                  , intent(in   ) :: i

    Table_2DLogLogLin_Y=exp(self%yv(i))
    return
  end function Table_2DLogLogLin_Y

  double precision function Table_2DLogLogLin_Z(self,i,j,table)
    !% Return the {\normalfont \ttfamily (i,j)}$^{\mathrm th}$ $x$-value for a 2D log-log table.
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
    !% Return the $x$-values for a 2D log-log table.
    implicit none
    class(table2DLogLogLin), intent(in   )             :: self
    double precision       , dimension(size(self%xv))  :: Table_2DLogLogLin_Xs

    Table_2DLogLogLin_Xs=exp(self%xv)
    return
  end function Table_2DLogLogLin_Xs

  function Table_2DLogLogLin_Ys(self)
    !% Return the $y$-values for a 2D log-log table.
    implicit none
    class(table2DLogLogLin), intent(in   )             :: self
    double precision       , dimension(size(self%xv))  :: Table_2DLogLogLin_Ys

    Table_2DLogLogLin_Ys=exp(self%yv)
    return
  end function Table_2DLogLogLin_Ys

  function Table_2DLogLogLin_Zs(self,table)
    !% Return the $y$-values for a 2D log-log table.
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
    !% Populate a 2-D log-log-linear table.
    use Galacticus_Error
    implicit none
    class           (table2DLogLogLin)                , intent(inout)           :: self
    double precision                  , dimension(:,:), intent(in   )           :: z
    integer                                           , intent(in   ), optional :: table
    integer                                                                     :: tableActual

    ! Validate the input.
    if (.not.allocated(self%zv)) call Galacticus_Error_Report("Table_2DLogLogLin_Populate","create the table before populating it")
    if     (                                                                                                &
         &   size(self%zv,dim=1) /= size(z,dim=1)                                                           &
         &  .or.                                                                                            &
         &   size(self%zv,dim=2) /= size(z,dim=2)                                                           &
         & ) call Galacticus_Error_Report("Table_2DLogLogLin_Populate","provided z array is of wrong size")
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
    !% Populate a single element of a 2-D log-log-linear table.
    use Galacticus_Error
    implicit none
    class           (table2DLogLogLin), intent(inout)           :: self
    double precision                  , intent(in   )           :: z
    integer                           , intent(in   )           :: i          , j
    integer                           , intent(in   ), optional :: table
    integer                                                     :: tableActual

    ! Validate the input.
    if (.not.allocated(self%zv)           ) call Galacticus_Error_Report("Table_2DLogLogLin_Populate_Single","create the table before populating it")
    if (i < 1 .or. i > size(self%zv,dim=1)) call Galacticus_Error_Report("Table_2DLogLogLin_Populate_Single","provided i value is out of bounds"    )
    if (j < 1 .or. j > size(self%zv,dim=2)) call Galacticus_Error_Report("Table_2DLogLogLin_Populate_Single","provided j value is out of bounds"    )
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
    !% Return the size of a 2D log-log-linear table.
    use Galacticus_Error
    implicit none
    class  (table2DLogLogLin), intent(in   ) :: self
    integer                  , intent(in   ) :: dim

    select case (dim)
    case (1)
       Table_2DLogLogLin_Size=self%xCount
    case (2)
       Table_2DLogLogLin_Size=self%yCount
    case default
       call Galacticus_Error_Report('Table_2DLogLogLin_Size','1 ≤ dim ≤ 2 is required')
    end select
    return
  end function Table_2DLogLogLin_Size

  double precision function Table_2DLogLogLin_Interpolate(self,x,y,table)
    !% Perform linear interpolation in a logarithmic 1D table.
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
    !% Perform linear interpolation in a logarithmic 1D table.
    use Galacticus_Error
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
          call Galacticus_Error_Report('Table_2DLogLogLin_Interpolate_Gradient','1 ≤ dim ≤ 2 is required')
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
    !% Destroy a 2D log-log-linear table.
    use Memory_Management
    implicit none
    class(table2DLogLogLin), intent(inout) :: self

    if (allocated(self%xv)) call Dealloc_Array(self%xv)
    if (allocated(self%yv)) call Dealloc_Array(self%yv)
    if (allocated(self%zv)) call Dealloc_Array(self%zv)
    return
  end subroutine Table_2DLogLogLin_Destroy

  logical function Table_2DLogLogLin_Is_Initialized(self)
    !% Return true if a 2D log-log-linear table has been created.
    implicit none
    class(table2DLogLogLin), intent(in   ) :: self

    Table_2DLogLogLin_Is_Initialized=allocated(self%zv)
    return
  end function Table_2DLogLogLin_Is_Initialized

  subroutine Table_Linear_Monotone_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !% Create a 1-D linear table.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class           (table1DLinearMonotoneCSpline), intent(inout)           :: self
    double precision                              , intent(in   )           :: xMaximum         , xMinimum
    integer                                       , intent(in   )           :: xCount
    integer                                       , intent(in   ), optional :: extrapolationType, tableCount
    integer                                                                 :: tableCountActual

    ! Determine number of tables.
    tableCountActual=1
    if (present(tableCount)) tableCountActual=tableCount
    ! Allocate arrays and construct the x-range.
    self%xCount=xCount
    call Alloc_Array(self%xv,[xCount                   ])
    call Alloc_Array(self%yv,[xCount  ,tableCountActual])
    call Alloc_Array(self%c1,[xCount+1,tableCountActual])
    call Alloc_Array(self%c2,[xCount  ,tableCountActual])
    call Alloc_Array(self%c3,[xCount  ,tableCountActual])
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
  end subroutine Table_Linear_Monotone_CSpline_1D_Create

  subroutine Table_Linear_Monotone_CSpline_1D_Destroy(self)
    !% Destroy a linear cubic-sline 1-D table.
    use Memory_Management
    implicit none
    class(table1DLinearMonotoneCSpline), intent(inout) :: self

    call Table_1D_Destroy(self)
    if (allocated(self%c1)) call Dealloc_Array(self%c1)
    if (allocated(self%c2)) call Dealloc_Array(self%c2)
    if (allocated(self%c3)) call Dealloc_Array(self%c3)
    return
  end subroutine Table_Linear_Monotone_CSpline_1D_Destroy

  subroutine Table_Linear_Monotone_CSpline_1D_Populate(self,y,table,computeSpline)
    !% Populate a 1-D linear table.
    use Galacticus_Error
    implicit none
    class           (table1DLinearMonotoneCSpline)              , intent(inout)           :: self
    double precision                              , dimension(:), intent(in   )           :: y
    integer                                                     , intent(in   ), optional :: table
    logical                                                     , intent(in   ), optional :: computeSpline
    integer                                                                               :: tableActual
    logical                                                                               :: computeSplineActual

    ! Validate the input.
    if (.not.allocated(self%yv)       ) call Galacticus_Error_Report("Table_Linear_Monotone_CSpline_1D_Populate","create the table before populating it")
    if (size(self%yv,dim=1) /= size(y)) call Galacticus_Error_Report("Table_Linear_Monotone_CSpline_1D_Populate","provided y array is of wrong size"    )

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
    !% Populate a single element of a 1-D linear table.
    use Galacticus_Error
    implicit none
    class           (table1DLinearMonotoneCSpline), intent(inout)           :: self
    double precision                              , intent(in   )           :: y
    integer                                       , intent(in   )           :: i
    integer                                       , intent(in   ), optional :: table
    logical                                       , intent(in   ), optional :: computeSpline
    integer                                                                 :: tableActual
    logical                                                                 :: computeSplineActual

    ! Validate the input.
    if (.not.allocated(self%yv)           ) call Galacticus_Error_Report("Table_Linear_Monotone_CSpline_1D_Populate_Single","create the table before populating it")
    if (i < 1 .or. i > size(self%yv,dim=1)) call Galacticus_Error_Report("Table_Linear_Monotone_CSpline_1D_Populate_Single","provided i value is out of bounds"    )

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
    !% Compute the interpolating spline factors for a 1-D linear spline.
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
    allocate(dx(size(self%xv)))
    allocate(dy(size(self%xv)))
    allocate( m(size(self%xv)))
    ! Get consecutive differences and slopes.
    do i=1,size(self%xv)-1
       dx(i)=self%xv(i+1      )-self%xv(i      )
       dy(i)=self%yv(i+1,table)-self%yv(i,table)
       m (i)=dy(i)/dx(i)
    end do
    ! Get degree-1 coefficients.
    self%c1(1,table)=m(1)
    do i=1,size(self%xv)-1
       if (m(i)*m(i+1) <= 0.0d0) then
          self%c1(i+1,table)=0.0d0
       else
          dxSum=dx(i)+dx(i+1)
          self%c1(i+1,table)=3.0d0*dxSum/((dxSum+dx(i+1))/m(i)+(dxSum+dx(i))/m(i+1))
       end if
    end do
    self%c1(size(self%xv)+1,table)=m(size(self%xv))
    ! Get degree-2 and degree-3 coefficients.
    do i=1,size(self%xv)
       factor=self%c1(i,table)+self%c1(i+1,table)-2.0d0*m(i)
       self%c2(i,table)=(m(i)-self%c1(i,table)-factor)/dx(i)
       self%c3(i,table)=factor/dx(i)**2
    end do
    ! Destroy workspace.
    deallocate(dx)
    deallocate(dy)
    deallocate( m)
    return
  end subroutine Table_Linear_Monotone_CSpline_1D_Compute_Spline

  double precision function Table_Linear_Monotone_CSpline_1D_Interpolate(self,x,table)
    !% Perform linear interpolation in a linear 1D table.
    implicit none
    class           (table1DLinearMonotoneCSpline), intent(inout)           :: self
    double precision                              , intent(in   )           :: x
    integer                                       , intent(in   ), optional :: table
    integer                                                                 :: i    , tableActual
    double precision                                                        :: dx   , xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%xPrevious .or. tableActual /= self%tablePrevious) then
       xEffective=self%xEffective(x)
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
       self%    yPrevious=self%yv(i,tableActual)+self%c1(i,tableActual)*dx+self%c2(i,tableActual)*dx**2+self%c3(i,tableActual)*dx**3
    end if
    Table_Linear_Monotone_CSpline_1D_Interpolate=self%yPrevious
    return
  end function Table_Linear_Monotone_CSpline_1D_Interpolate

  double precision function Table_Linear_Monotone_CSpline_1D_Interpolate_Gradient(self,x,table)
    !% Perform linear interpolation in a linear 1D table.
    implicit none
    class           (table1DLinearMonotoneCSpline), intent(inout)           :: self
    double precision                              , intent(in   )           :: x
    integer                                       , intent(in   ), optional :: table
    integer                                                                 :: i    , tableActual
    double precision                                                        :: dx   , xEffective

    ! Determine which table to use.
    tableActual=1
    if (present(table)) tableActual=table
    ! Check for recall with same value as previous call.
    if (x /= self%dxPrevious .or. tableActual /= self%dTablePrevious) then
       xEffective=self%xEffective(x)
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
       self%    dyPrevious=self%c1(i,tableActual)+2.0d0*self%c2(i,tableActual)*dx+3.0d0*self%c3(i,tableActual)*dx**2
    end if
    Table_Linear_Monotone_CSpline_1D_Interpolate_Gradient=self%dyPrevious
    return
  end function Table_Linear_Monotone_CSpline_1D_Interpolate_Gradient

  function Table_Linear_Monotone_CSpline_Integration_Weights(self,x0,x1,integrand)
    !% Returns a set of weights for trapezoidal integration on the table between limits {\normalfont \ttfamily x0} and {\normalfont \ttfamily x1}.
    use Galacticus_Error
    implicit none
    class           (table1DLinearMonotoneCSpline), intent(inout)                               :: self
    double precision                              , intent(in   )                               :: x0, x1
    procedure       (integrandTemplate           ), intent(in   )           , pointer, optional :: integrand
    double precision                              , dimension(size(self%xv))                    :: Table_Linear_Monotone_CSpline_Integration_Weights

    call Galacticus_Error_Report('Table_Linear_Monotone_CSpline_Integration_Weights','integration weights not supported')
    return
  end function Table_Linear_Monotone_CSpline_Integration_Weights

  subroutine Table_Logarithmic_Monotone_CSpline_1D_Create(self,xMinimum,xMaximum,xCount,tableCount,extrapolationType)
    !% Create a 1-D logarithmic table.
    implicit none
    class           (table1DLogarithmicMonotoneCSpline), intent(inout)           :: self
    double precision                                   , intent(in   )           :: xMaximum         , xMinimum
    integer                                            , intent(in   )           :: xCount
    integer                                            , intent(in   ), optional :: extrapolationType, tableCount

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
    !% Return the {\normalfont \ttfamily i}$^{\mathrm th}$ $x$-value for a logarithmic 1D table.
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
    !% Return the $x$-values for a 1D table.
    implicit none
    class           (table1DLogarithmicMonotoneCSpline), intent(in   )            :: self
    double precision                                   , dimension(size(self%xv)) :: Table_Logarithmic_Monotone_CSpline_1D_Xs

    Table_Logarithmic_Monotone_CSpline_1D_Xs=exp(self%table1DLinearMonotoneCSpline%xs())
    return
  end function Table_Logarithmic_Monotone_CSpline_1D_Xs

  double precision function Table_Logarithmic_Monotone_CSpline_1D_Interpolate(self,x,table)
    !% Perform linear interpolation in a logarithmic 1D table.
    implicit none
    class           (table1DLogarithmicMonotoneCSpline), intent(inout)           :: self
    double precision                                   , intent(in   )           :: x
    integer                                            , intent(in   ), optional :: table

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_Monotone_CSpline_1D_Interpolate=self%table1DLinearMonotoneCSpline%interpolate(self%xLogarithmicPrevious,table)
    return
  end function Table_Logarithmic_Monotone_CSpline_1D_Interpolate

  double precision function Table_Logarithmic_Monotone_CSpline_1D_Interpolate_Gradient(self,x,table)
    !% Perform linear interpolation in a logarithmic 1D table.
    implicit none
    class           (table1DLogarithmicMonotoneCSpline), intent(inout)           :: self
    double precision                                   , intent(in   )           :: x
    integer                                            , intent(in   ), optional :: table

    if (.not.self%previousSet .or. x /= self%xLinearPrevious) then
       self%previousSet         =.true.
       self%xLinearPrevious     =    x
       self%xLogarithmicPrevious=log(x)
    end if
    Table_Logarithmic_Monotone_CSpline_1D_Interpolate_Gradient=self%table1DLinearMonotoneCSpline%interpolateGradient(self%xLogarithmicPrevious,table)/self%xEffective(x)
    return
  end function Table_Logarithmic_Monotone_CSpline_1D_Interpolate_Gradient
  
end module Tables
