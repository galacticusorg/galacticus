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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which implements multilinear interpolation in tables of arbitrary dimensionality.
!!}

module Numerical_Interpolation_MultiD
  !!{RST
  Implements multilinear interpolation in tables of arbitrary dimensionality, on a grid which is separable---that is,
  one specified by an independent array of abscissae along each dimension.

  Several tabulations in Galacticus are functions of four or more variables. Interpolating in them means locating the
  bracketing pair of abscissae along each dimension and then forming the weighted sum over the :math:`2^N` corners of
  the enclosing hypercube. That pattern is mechanical but fiddly to write out by hand for each new table, and is easy
  to get subtly wrong. This class factors it out.

  Each dimension is described by an ordinary :galacticus-class:`interpolator` object, so extrapolation behaviour is
  specified per-dimension in the usual way, and a dimension which should be interpolated logarithmically is handled
  by constructing its interpolator on---and evaluating it at---the logarithm of the coordinate.

  Table values are indexed in Fortran's usual column-major order, so that the value at grid point
  :math:`(i_1,i_2,\ldots,i_N)` lies at position :math:`i_1 + n_1 (i_2-1) + n_1 n_2 (i_3-1) + \ldots` of the flattened
  array, where :math:`n_d` is the number of abscissae along dimension :math:`d`. A contiguous rank-\ :math:`N` array
  whose shape matches that returned by the ``shapeTable`` method may therefore be passed directly to the
  ``interpolate`` method.

  Where the factors along each dimension are shared between several tables, or are unchanged between successive
  interpolations, they may be found once with the ``factors`` method and then passed to the ``cornersFactors`` and
  ``interpolateFactors`` methods, so that the bracketing search along each dimension is not repeated.
  !!}
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator
  implicit none
  private
  public :: interpolatorMultiD

  ! Maximum number of dimensions supported. The limit exists only to keep the number of hypercube corners, 2^N, within
  ! a default integer; it is far above any plausible tabulation.
  integer, parameter :: countDimensionsMaximum=16

  !![
  <stateStorable class="interpolatorMultiD">
   <interpolatorMultiD>
    <methodCall method="GSLReallocate"/>
   </interpolatorMultiD>
  </stateStorable>
  <deepCopyActions class="interpolatorMultiD">
   <interpolatorMultiD>
    <methodCall method="GSLReallocate"/>
   </interpolatorMultiD>
  </deepCopyActions>
  !!]

  type :: interpolatorMultiD
     !!{RST
     A multilinear interpolator over a separable grid of arbitrary dimensionality.
     !!}
     private
     type   (interpolator), allocatable, dimension(:) :: dimensions
     integer(c_size_t    ), allocatable, dimension(:) :: countPoints
     integer(c_size_t    ), allocatable, dimension(:) :: stride
     integer                                          :: countDimensions_=0
     integer                                          :: countCorners_   =0
   contains
     !![
     <methods docformat="rst">
       <method method="countDimensions"    description="Return the number of dimensions of the interpolator."                  />
       <method method="countCorners"       description="Return the number of corners of the interpolating hypercube, i.e. 2^N."/>
       <method method="countTable"         description="Return the total number of entries in the flattened table."            />
       <method method="shapeTable"         description="Return the number of abscissae along each dimension."                  />
       <method method="factors"            description="Return the bracketing index and linear weights along each dimension."  />
       <method method="corners"            description="Return the flattened index of, and weight for, each hypercube corner." />
       <method method="cornersFactors"     description="As ``corners``, but from factors which have already been found."       />
       <method method="interpolate"        description="Interpolate in a flattened table of values at the given coordinates."  />
       <method method="interpolateFactors" description="As ``interpolate``, but from factors which have already been found."   />
       <method method="assignment(=)"      description="Assign multilinear interpolator objects."                              />
       <method method="gslReallocate"      description="Reallocate the GSL objects held along each dimension."                 />
     </methods>
     !!]
     procedure ::                           interpolatorMultiDAssign
     generic   :: assignment(=)          => interpolatorMultiDAssign
     procedure :: countDimensions    => interpolatorMultiDCountDimensions
     procedure :: countCorners       => interpolatorMultiDCountCorners
     procedure :: countTable         => interpolatorMultiDCountTable
     procedure :: shapeTable         => interpolatorMultiDShapeTable
     procedure :: factors            => interpolatorMultiDFactors
     procedure :: corners            => interpolatorMultiDCorners
     procedure :: cornersFactors     => interpolatorMultiDCornersFactors
     procedure :: interpolate        => interpolatorMultiDInterpolate
     procedure :: interpolateFactors => interpolatorMultiDInterpolateFactors
     procedure :: gslReallocate      => interpolatorMultiDGSLReallocate
  end type interpolatorMultiD

  interface interpolatorMultiD
     !!{RST
     Constructors for the ``interpolatorMultiD`` class.
     !!}
     module procedure interpolatorMultiDConstructor
  end interface interpolatorMultiD

contains

  function interpolatorMultiDConstructor(interpolators) result(self)
    !!{RST
    Constructor for ``interpolatorMultiD`` objects. ``interpolators`` provides one :galacticus-class:`interpolator`
    per dimension, in the order in which the dimensions are laid out in the table---that is, the first element
    corresponds to the most rapidly varying index of the flattened table.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (interpolatorMultiD)                              :: self
    type   (interpolator      ), intent(in   ), dimension(:) :: interpolators
    integer                                                  :: i

    if (size(interpolators) < 1                     )                                                     &
         & call Error_Report('at least one dimension is required'             //{introspection:location})
    if (size(interpolators) > countDimensionsMaximum)                                                     &
         & call Error_Report('too many dimensions'                            //{introspection:location})
    self%countDimensions_=size(interpolators)
    self%countCorners_   =2**self%countDimensions_
    allocate(self%dimensions (self%countDimensions_))
    allocate(self%countPoints(self%countDimensions_))
    allocate(self%stride     (self%countDimensions_))
    do i=1,self%countDimensions_
       self%dimensions (i)=interpolators(i)
       self%countPoints(i)=interpolators(i)%count()
       if (self%countPoints(i) < 2_c_size_t)                                                              &
            & call Error_Report('each dimension requires at least 2 abscissae'//{introspection:location})
    end do
    ! Compute strides for column-major (Fortran) ordering of the flattened table.
    self%stride(1)=1_c_size_t
    do i=2,self%countDimensions_
       self%stride(i)=self%stride(i-1)*self%countPoints(i-1)
    end do
    return
  end function interpolatorMultiDConstructor

  subroutine interpolatorMultiDGSLReallocate(self)
    !!{RST
    Reallocate the GSL objects held by the :galacticus-class:`interpolator` of each dimension.

    This is needed after an ``interpolatorMultiD`` has been copied---or restored from stored state---by the
    code-generated deep copy and state store machinery. That machinery copies this object by intrinsic assignment,
    which leaves the interpolator of each dimension sharing GSL objects with the object copied from without having
    incremented their reference count. Reallocating releases that share and gives this copy GSL objects of its own, in
    exactly the same way as for a bare :galacticus-class:`interpolator`.
    !!}
    implicit none
    class  (interpolatorMultiD), intent(inout) :: self
    integer                                    :: i

    if (.not.allocated(self%dimensions)) return
    do i=1,size(self%dimensions)
       call self%dimensions(i)%gslReallocate()
    end do
    return
  end subroutine interpolatorMultiDGSLReallocate

  subroutine interpolatorMultiDAssign(self,from)
    !!{RST
    Perform assignment of ``interpolatorMultiD`` objects.

    Each dimension is held as an :galacticus-class:`interpolator`, which is a reference-counted handle on GSL objects.
    Those handles must therefore be copied using the ``interpolator`` class' own defined assignment, so that the
    reference count of the GSL objects is incremented. A copy which did not do so would leave this object holding
    pointers to GSL objects freed when the object copied from was destroyed---as happens whenever an
    ``interpolatorMultiD`` is built from interpolators which are local to the routine building it.
    !!}
    implicit none
    class  (interpolatorMultiD), intent(inout) :: self
    class  (interpolatorMultiD), intent(in   ) :: from
    integer                                    :: i

    self%countDimensions_=from%countDimensions_
    self%countCorners_   =from%countCorners_
    if (allocated(self%dimensions )) deallocate(self%dimensions )
    if (allocated(self%countPoints)) deallocate(self%countPoints)
    if (allocated(self%stride     )) deallocate(self%stride     )
    if (allocated(from%dimensions )) then
       allocate(self%dimensions(size(from%dimensions)))
       do i=1,size(from%dimensions)
          self%dimensions(i)=from%dimensions(i)
       end do
    end if
    if (allocated(from%countPoints))   allocate(self%countPoints,source=from%countPoints)
    if (allocated(from%stride     ))   allocate(self%stride     ,source=from%stride     )
    return
  end subroutine interpolatorMultiDAssign

  integer function interpolatorMultiDCountDimensions(self) result(countDimensions)
    !!{RST
    Return the number of dimensions of an ``interpolatorMultiD``.
    !!}
    implicit none
    class(interpolatorMultiD), intent(in   ) :: self

    countDimensions=self%countDimensions_
    return
  end function interpolatorMultiDCountDimensions

  integer function interpolatorMultiDCountCorners(self) result(countCorners)
    !!{RST
    Return the number of corners, :math:`2^N`, of the hypercube in which an ``interpolatorMultiD`` interpolates.
    !!}
    implicit none
    class(interpolatorMultiD), intent(in   ) :: self

    countCorners=self%countCorners_
    return
  end function interpolatorMultiDCountCorners

  function interpolatorMultiDCountTable(self) result(countTable)
    !!{RST
    Return the total number of entries in the flattened table of values expected by an ``interpolatorMultiD``.
    !!}
    implicit none
    integer(c_size_t          )                :: countTable
    class  (interpolatorMultiD), intent(in   ) :: self

    countTable=product(self%countPoints)
    return
  end function interpolatorMultiDCountTable

  function interpolatorMultiDShapeTable(self) result(shape_)
    !!{RST
    Return the number of abscissae along each dimension of an ``interpolatorMultiD``. The product of these is the
    required size of the flattened table of values.
    !!}
    implicit none
    integer(c_size_t          ), allocatable, dimension(:) :: shape_
    class  (interpolatorMultiD), intent(in   )             :: self

    allocate(shape_(self%countDimensions_))
    shape_=self%countPoints
    return
  end function interpolatorMultiDShapeTable

  subroutine interpolatorMultiDFactors(self,x,indices,weights)
    !!{RST
    Return, for each dimension of an ``interpolatorMultiD``, the index ``indices(d)`` of the lower of the pair of
    abscissae bracketing ``x(d)``, together with the linear interpolating weights ``weights(0:1,d)`` to be applied at
    that abscissa and the next.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (interpolatorMultiD), intent(inout)                 :: self
    double precision                    , intent(in   ), dimension(:  ) :: x
    integer         (c_size_t          ), intent(  out), dimension(:  ) :: indices
    double precision                    , intent(  out), dimension(:,:) :: weights
    integer                                                             :: i

    if (size(x                ) /= self%countDimensions_)                                  &
         & call Error_Report('coordinate vector has wrong size'//{introspection:location})
    if (size(indices          ) /= self%countDimensions_)                                  &
         & call Error_Report('index array has wrong size'      //{introspection:location})
    if (size(weights,dim=2) /= self%countDimensions_ .or. size(weights,dim=1) /= 2)        &
         & call Error_Report('weight array has wrong shape'    //{introspection:location})
    do i=1,self%countDimensions_
       call self%dimensions(i)%linearFactors(x(i),indices(i),weights(:,i))
    end do
    return
  end subroutine interpolatorMultiDFactors

  subroutine interpolatorMultiDCorners(self,x,indices,weights)
    !!{RST
    Return, for each of the :math:`2^N` corners of the hypercube of an ``interpolatorMultiD`` enclosing ``x``, the
    index ``indices(c)`` of that corner in the flattened table of values, together with the multilinear weight
    ``weights(c)`` to be applied to it. The weights sum to unity, except where a dimension is extrapolating under the
    ``zero`` extrapolation type, in which case they sum to zero.

    ``indices`` and ``weights`` must each be of size ``countCorners()``.
    !!}
    implicit none
    class           (interpolatorMultiD), intent(inout)                        :: self
    double precision                    , intent(in   ), dimension(:)          :: x
    integer         (c_size_t          ), intent(  out), dimension(:)          :: indices
    double precision                    , intent(  out), dimension(:)          :: weights
    integer         (c_size_t          ), dimension(    self%countDimensions_) :: indicesDimension
    double precision                    , dimension(0:1,self%countDimensions_) :: weightsDimension

    call self%factors       (x,indicesDimension,weightsDimension                )
    call self%cornersFactors(  indicesDimension,weightsDimension,indices,weights)
    return
  end subroutine interpolatorMultiDCorners

  subroutine interpolatorMultiDCornersFactors(self,indicesDimension,weightsDimension,indices,weights)
    !!{RST
    As ``corners``, but taking the per-dimension bracketing indices and linear weights which have already been found by
    the ``factors`` method, in place of the coordinates themselves.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (interpolatorMultiD), intent(in   )                  :: self
    integer         (c_size_t          ), intent(in   ), dimension(:   ) :: indicesDimension
    double precision                    , intent(in   ), dimension(0:,:) :: weightsDimension
    integer         (c_size_t          ), intent(  out), dimension(:   ) :: indices
    double precision                    , intent(  out), dimension(:   ) :: weights
    integer                                                              :: corner          , i, &
         &                                                                  offset

    if (size(indicesDimension) /= self%countDimensions_)                                              &
         & call Error_Report('index array has wrong size'   //{introspection:location})
    if (size(weightsDimension,dim=2) /= self%countDimensions_ .or. size(weightsDimension,dim=1) /= 2) &
         & call Error_Report('weight array has wrong shape' //{introspection:location})
    if (size(indices) /= self%countCorners_ .or. size(weights) /= self%countCorners_)                 &
         & call Error_Report('corner arrays have wrong size'//{introspection:location})
    ! Enumerate the corners of the hypercube. Bit `i-1` of the corner number selects the lower (0) or upper (1)
    ! abscissa along dimension `i`.
    do corner=0,self%countCorners_-1
       indices(corner+1)=1_c_size_t
       weights(corner+1)=1.0d0
       do i=1,self%countDimensions_
          offset           =ibits(corner,i-1,1)
          indices(corner+1)=+indices(corner+1)            &
               &            +(                            &
               &              +indicesDimension(i)        &
               &              +int(offset,kind=c_size_t)  &
               &              -1_c_size_t                 &
               &             )                            &
               &            *self%stride     (i)
          weights(corner+1)=+weights         (corner+1  ) &
               &            *weightsDimension(offset  ,i)
       end do
    end do
    return
  end subroutine interpolatorMultiDCornersFactors

  double precision function interpolatorMultiDInterpolate(self,values,x) result(y)
    !!{RST
    Interpolate in the table ``values`` at the coordinates ``x``. ``values`` holds the tabulated function in Fortran's
    usual column-major ordering, so a contiguous rank-\ :math:`N` array whose shape matches that returned by the
    ``shapeTable`` method may be passed directly.
    !!}
    implicit none
    class           (interpolatorMultiD), intent(inout)                        :: self
    double precision                    , intent(in   ), dimension(*)          :: values
    double precision                    , intent(in   ), dimension(:)          :: x
    integer         (c_size_t          ), dimension(    self%countDimensions_) :: indicesDimension
    double precision                    , dimension(0:1,self%countDimensions_) :: weightsDimension

    call self%factors(x,indicesDimension,weightsDimension)
    y=self%interpolateFactors(values,indicesDimension,weightsDimension)
    return
  end function interpolatorMultiDInterpolate

  double precision function interpolatorMultiDInterpolateFactors(self,values,indicesDimension,weightsDimension) result(y)
    !!{RST
    As ``interpolate``, but taking the per-dimension bracketing indices and linear weights which have already been found
    by the ``factors`` method, in place of the coordinates themselves.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (interpolatorMultiD), intent(in   )                  :: self
    double precision                    , intent(in   ), dimension(*   ) :: values
    integer         (c_size_t          ), intent(in   ), dimension(:   ) :: indicesDimension
    double precision                    , intent(in   ), dimension(0:,:) :: weightsDimension
    integer         (c_size_t          )                                 :: index_
    double precision                                                     :: value_
    integer                                                              :: corner          , i

    if (size(indicesDimension) /= self%countDimensions_)                                              &
         & call Error_Report('index array has wrong size'  //{introspection:location})
    if (size(weightsDimension,dim=2) /= self%countDimensions_ .or. size(weightsDimension,dim=1) /= 2) &
         & call Error_Report('weight array has wrong shape'//{introspection:location})
    ! Sum over the corners of the hypercube. Bit `i-1` of the corner number selects the lower (0) or upper (1) abscissa
    ! along dimension `i`. The weights are applied to the tabulated value one dimension at a time, rather than being
    ! multiplied out into a single weight per corner and applied to the value at the end. That saves a multiply per
    ! corner, and forms each term in the same order as a corner sum written out by hand as nested loops over the
    ! dimensions---so a hand-rolled sum converted to use this class gives bit-identical results.
    y=0.0d0
    do corner=0,self%countCorners_-1
       ! Locate this corner in the flattened table.
       index_=1_c_size_t
       do i=1,self%countDimensions_
          index_=+index_                                   &
               & +(                                        &
               &   +indicesDimension(i)                    &
               &   +int(ibits(corner,i-1,1),kind=c_size_t) &
               &   -1_c_size_t                             &
               &  )                                        &
               & *self%stride       (i)
       end do
       ! Accumulate the tabulated value at this corner, weighted along each dimension in turn.
       value_=values(index_)
       do i=1,self%countDimensions_
          value_=value_*weightsDimension(ibits(corner,i-1,1),i)
       end do
       y=y+value_
    end do
    return
  end function interpolatorMultiDInterpolateFactors

end module Numerical_Interpolation_MultiD
