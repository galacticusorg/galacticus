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

!!{RST
Contains a module which implements construction of numerical ranges.
!!}

module Numerical_Ranges
  !!{RST
  Implements construction of numerical ranges.
  !!}
  implicit none
  private
  public :: Make_Range, Range_Pinned

  ! Parameters to specify type of range required.
  integer, parameter, public :: rangeTypeLinear=0, rangeTypeLogarithmic=1, rangeTypeUndefined=-1

  !![
  <enumeration docformat="rst">
   <name>gridScheme</name>
   <description>
   Used to specify the scheme by which an absolute lattice of tabulation points is laid out. A lattice with ``pointsPer``
   :math:`=N` consists of the points :math:`x_j` for integer :math:`j`, where :math:`x_j = 10^{j/N}` for ``perDecade``,
   :math:`x_j = 2^{j/N}` for ``perOctave``, and :math:`x_j = j/N` for ``perUnit``. Since the lattice is a function of
   :math:`j` and :math:`N` alone it is independent of the range over which any particular table is built, so two tables
   built on the same lattice have identical abscissae wherever their ranges overlap.
   </description>
   <entry label="undefined" description="No lattice is defined."                                                        />
   <entry label="perDecade" description="Points are spaced uniformly in :math:`\log_{10} x`, ``pointsPer`` per decade." />
   <entry label="perOctave" description="Points are spaced uniformly in :math:`\log_{2} x`, ``pointsPer`` per octave."  />
   <entry label="perUnit"   description="Points are spaced uniformly in :math:`x`, ``pointsPer`` per unit interval."    />
  </enumeration>
  !!]

  ! Tolerance used when rounding lattice coordinates to integer lattice indices. This is required because, for example,
  ! log10(1000) evaluates to 2.9999999999999996 in double precision, which would otherwise be floored to 2 instead of 3. The
  ! tolerance is applied in units of lattice steps (or of decades/octaves/unit intervals when anchoring to those), so a value of
  ! 10⁻⁶ shifts a bound by at most one part in 10⁶ of a step - far below any level at which a tabulation is meaningful.
  double precision, parameter :: toleranceLattice=1.0d-6

  type, public :: rangeLattice
     !!{RST
     A specification of a contiguous set of points on an absolute lattice. The lattice itself is fixed by ``scheme`` and
     ``pointsPer``; the particular set of points is the ``count`` points with lattice indices ``indexMinimum`` to
     ``indexMinimum``\ ``+``\ ``count``\ ``-1``.
     !!}
     type   (enumerationGridSchemeType) :: scheme      =gridSchemeUndefined
     integer                            :: pointsPer   =0
     integer                            :: indexMinimum=0
     integer                            :: count       =0
   contains
     !![
     <methods docformat="rst">
       <method description="Return true if this object specifies a usable lattice."                                            method="isDefined"        />
       <method description="Return the lattice index of the final point."                                                      method="indexMaximum"     />
       <method description="Return the value of the first point."                                                              method="minimum"          />
       <method description="Return the value of the final point."                                                              method="maximum"          />
       <method description="Return the values of all points."                                                                  method="values"           />
       <method description="Return the natural logarithms of the values of all points (logarithmic gridding schemes only)."     method="valuesLogarithmic"/>
       <method description="Return the spacing of the points (``perUnit`` gridding scheme only)."                              method="step"             />
       <method description="Return the spacing of the natural logarithms of the points (logarithmic gridding schemes only)."   method="stepLogarithmic"  />
       <method description="Return true if ``lattice`` lies on the same lattice as this object, and all of its points are points of this object." method="covers"  />
     </methods>
     !!]
     procedure :: isDefined         => rangeLatticeIsDefined
     procedure :: indexMaximum      => rangeLatticeIndexMaximum
     procedure :: minimum           => rangeLatticeMinimum
     procedure :: maximum           => rangeLatticeMaximum
     procedure :: values            => rangeLatticeValues
     procedure :: valuesLogarithmic => rangeLatticeValuesLogarithmic
     procedure :: step              => rangeLatticeStep
     procedure :: stepLogarithmic   => rangeLatticeStepLogarithmic
     procedure :: covers            => rangeLatticeCovers
  end type rangeLattice

  interface Range_Pinned
     !!{RST
     Generic interface to functions which construct a pinned range.
     !!}
     module procedure Range_Pinned_Scalar
     module procedure Range_Pinned_Array
  end interface Range_Pinned

contains

  recursive function Make_Range(rangeMinimum,rangeMaximum,rangeNumber,rangeType,rangeBinned) result (rangeValues)
    !!{RST
    Builds a numerical range between ``rangeMinimum`` and ``rangeMaximum`` using ``rangeNumber`` points and spacing as specified by ``rangeType`` (defaulting to linear spacing if no ``rangeType`` is given).
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision, intent(in   )           :: rangeMaximum                  , rangeMinimum
    integer         , intent(in   )           :: rangeNumber
    integer         , intent(in   ), optional :: rangeType
    logical         , intent(in   ), optional :: rangeBinned
    double precision                          :: rangeValues      (rangeNumber)
    integer                                   :: iRange                        , rangeTypeActual
    logical                                   :: rangeBinnedActual

    ! Find what type of range is required.
    if (present(rangeType)) then
       rangeTypeActual=rangeType
    else
       rangeTypeActual=rangeTypeLinear
    end if
    ! Determine if a binned range is required.
    if (present(rangeBinned)) then
       rangeBinnedActual=rangeBinned
    else
       rangeBinnedActual=.false.
    end if
    ! Build the range.
    select case (rangeTypeActual)
    case (rangeTypeLinear)
       ! Determine if we are returning a regular range or bin centers.
       if (rangeBinnedActual) then
          ! Build a linear binned range returning bin centers.
          forall(iRange=1:rangeNumber)
             rangeValues(iRange)=rangeMinimum+(rangeMaximum-rangeMinimum)*(dble(iRange-1)+0.5d0)/dble(rangeNumber)
          end forall
       else
          ! Check that the rangeNumber is valid.
          if (rangeNumber <= 1) call Error_Report('number of points in range must exceed 1'//{introspection:location})
          ! Build a linear range.
          forall(iRange=1:rangeNumber)
             rangeValues(iRange)=rangeMinimum+(rangeMaximum-rangeMinimum)*dble(iRange-1)/dble(rangeNumber-1)
          end forall
          ! Ensure upper limit is precise.
          rangeValues(rangeNumber)=rangeMaximum
       end if
    case (rangeTypeLogarithmic)
       ! Call ourself with logged limits and then exponentiate the result.
       rangeValues=exp(Make_Range(log(rangeMinimum),log(rangeMaximum),rangeNumber,rangeTypeLinear,rangeBinnedActual))
    case default
       call Error_Report('range type is unrecognized'//{introspection:location})
    end select
    return
  end function Make_Range

  function Range_Pinned_Scalar(valueTarget,pointsPer,scheme,marginFactor,marginOffset,latticeCurrent,rangeCurrent,limitMinimum,limitMaximum,anchorEvery) result(lattice)
    !!{RST
    Return the ``rangeLattice`` covering the single value ``valueTarget``---see ``Range_Pinned_Array`` for a description of the arguments.
    !!}
    implicit none
    type            (rangeLattice             )                                        :: lattice
    double precision                           , intent(in   )                         :: valueTarget
    integer                                    , intent(in   )                         :: pointsPer
    type            (enumerationGridSchemeType), intent(in   )                         :: scheme
    double precision                           , intent(in   ), optional               :: marginFactor  , marginOffset, &
         &                                                                                limitMinimum  , limitMaximum
    type            (rangeLattice             ), intent(in   ), optional               :: latticeCurrent
    double precision                           , intent(in   ), optional, dimension(2) :: rangeCurrent
    integer                                    , intent(in   ), optional               :: anchorEvery

    lattice=Range_Pinned_Array(                                &
         &                     [valueTarget]                 , &
         &                      pointsPer                    , &
         &                      scheme                       , &
         &                      marginFactor  =marginFactor  , &
         &                      marginOffset  =marginOffset  , &
         &                      latticeCurrent=latticeCurrent, &
         &                      rangeCurrent  =rangeCurrent  , &
         &                      limitMinimum  =limitMinimum  , &
         &                      limitMaximum  =limitMaximum  , &
         &                      anchorEvery   =anchorEvery     &
         &                    )
    return
  end function Range_Pinned_Scalar

  function Range_Pinned_Array(valueTarget,pointsPer,scheme,marginFactor,marginOffset,latticeCurrent,rangeCurrent,limitMinimum,limitMaximum,anchorEvery) result(lattice)
    !!{RST
    Return a ``rangeLattice``---a set of points on the absolute lattice specified by ``scheme`` and
    ``pointsPer``---which is guaranteed to cover every value in ``valueTarget``, with a safety margin. The bounds are found
    by, in order:

    #. applying the safety margin to the range of ``valueTarget``: the range is expanded by a factor ``marginFactor``
       (default :math:`2`) at each end for the logarithmic gridding schemes, or by an additive offset ``marginOffset``
       (default :math:`1`) for the ``perUnit`` scheme;
    #. taking the union with ``rangeCurrent``, if given (used when the current range of a table is not itself known to lie
       on the lattice);
    #. *pinning* the bounds outward to lattice points whose index is a multiple of ``anchorEvery``;
    #. clamping to the hard limits ``limitMinimum`` and/or ``limitMaximum``, if given, with the clamped edge snapped
       *inward* to the lattice so that all points remain lattice points;
    #. taking the union with ``latticeCurrent``, if given---performed in exact integer arithmetic, and applied *after* the
       clamping so that a lattice which is already in use is never shrunk (and so that
       extension of a table built on it can never fail).

    ``anchorEvery`` is the interval, measured in lattice steps, to which the bounds are pinned; it controls only how far
    beyond the requested range the table extends, and is independent of the grid density ``pointsPer`` which controls the
    accuracy of interpolation within it. It defaults to ``pointsPer``---that is, to whole decades, octaves, or unit
    intervals---which gives the coarsest set of possible ranges and hence the strongest determinism. A value of :math:`1`
    pins to the lattice points themselves. Intermediate values are useful where a whole decade of extra tabulation is
    physically or computationally extravagant: on a cosmic time axis, for example, a whole decade spans most of the history
    of the universe, whereas ``anchorEvery``\ ``=``\ ``pointsPer``\ ``/2`` pins to half decades (a factor of
    :math:`\sqrt{10}`). Determinism, exact value reuse on extension, and the mergeability of stored tables hold for *any*
    choice: only the granularity of the discrete set of possible ranges changes. A divisor of ``pointsPer`` gives the most
    natural anchor points, but is not required.

    Because the resulting bounds are lattice indices, the point count is exact, and any two ranges built on the same
    lattice have bit-identical abscissae wherever they overlap.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (rangeLattice             )                                        :: lattice
    double precision                           , intent(in   )          , dimension(:) :: valueTarget
    integer                                    , intent(in   )                         :: pointsPer
    type            (enumerationGridSchemeType), intent(in   )                         :: scheme
    double precision                           , intent(in   ), optional               :: marginFactor     , marginOffset     , &
         &                                                                                limitMinimum     , limitMaximum
    type            (rangeLattice             ), intent(in   ), optional               :: latticeCurrent
    double precision                           , intent(in   ), optional, dimension(2) :: rangeCurrent
    integer                                    , intent(in   ), optional               :: anchorEvery
    double precision                                                                   :: marginFactor_    , marginOffset_    , &
         &                                                                                coordinateMinimum, coordinateMaximum, &
         &                                                                                anchorsPerUnit
    integer                                                                            :: anchorEvery_
    integer                                                                            :: indexMinimum     , indexMaximum

    ! Validate the arguments.
    anchorEvery_=pointsPer
    if (present(anchorEvery)) anchorEvery_=anchorEvery
    if (size(valueTarget) <  1) call Error_Report('no target values were provided'                    //{introspection:location})
    if (pointsPer         <  1) call Error_Report('number of points per interval must be positive'    //{introspection:location})
    if (anchorEvery_      <  1) call Error_Report('anchor interval must be at least one lattice step' //{introspection:location})
    ! Apply the safety margin.
    coordinateMinimum=0.0d0
    coordinateMaximum=0.0d0
    select case (scheme%ID)
    case (gridSchemePerDecade%ID,gridSchemePerOctave%ID)
       if (present(marginOffset)      ) call Error_Report('`marginOffset` applies only to the `perUnit` gridding scheme - use `marginFactor`'//{introspection:location})
       marginFactor_=2.0d0
       if (present(marginFactor)      ) marginFactor_=marginFactor
       if (marginFactor_     <  1.0d0 ) call Error_Report('`marginFactor` must be at least unity'                                           //{introspection:location})
       if (any(valueTarget   <= 0.0d0)) call Error_Report('target values must be positive for logarithmic gridding schemes'                 //{introspection:location})
       coordinateMinimum=Lattice_Coordinate(minval(valueTarget)/marginFactor_,scheme)
       coordinateMaximum=Lattice_Coordinate(maxval(valueTarget)*marginFactor_,scheme)
    case (gridSchemePerUnit  %ID                       )
       if (present(marginFactor)      ) call Error_Report('`marginFactor` applies only to logarithmic gridding schemes - use `marginOffset`' //{introspection:location})
       marginOffset_=1.0d0
       if (present(marginOffset)      ) marginOffset_=marginOffset
       if (marginOffset_     <  0.0d0 ) call Error_Report('`marginOffset` must be non-negative'                                             //{introspection:location})
       coordinateMinimum=minval(valueTarget)-marginOffset_
       coordinateMaximum=maxval(valueTarget)+marginOffset_
    case default
       call Error_Report('unrecognized gridding scheme'//{introspection:location})
    end select
    ! Take the union with any current range which is not known to lie on the lattice.
    if (present(rangeCurrent)) then
       coordinateMinimum=min(coordinateMinimum,Lattice_Coordinate(minval(rangeCurrent),scheme))
       coordinateMaximum=max(coordinateMaximum,Lattice_Coordinate(maxval(rangeCurrent),scheme))
    end if
    ! Pin the bounds outward onto those lattice points whose index is a multiple of the anchor interval. The number of anchor
    ! intervals per unit of the lattice coordinate is `pointsPer/anchorEvery` - note that this evaluates exactly to unity when
    ! anchoring to whole decades/octaves/unit intervals, and exactly to `pointsPer` when anchoring to the lattice points
    ! themselves, so that those two cases involve no additional round-off.
    anchorsPerUnit=dble(pointsPer)/dble(anchorEvery_)
    indexMinimum  =Lattice_Floor  (coordinateMinimum*anchorsPerUnit)*anchorEvery_
    indexMaximum  =Lattice_Ceiling(coordinateMaximum*anchorsPerUnit)*anchorEvery_
    ! Apply any hard limits, snapping the clamped edge inward to the lattice.
    if (present(limitMinimum)) indexMinimum=max(indexMinimum,Lattice_Ceiling(Lattice_Coordinate(limitMinimum,scheme)*dble(pointsPer)))
    if (present(limitMaximum)) indexMaximum=min(indexMaximum,Lattice_Floor  (Lattice_Coordinate(limitMaximum,scheme)*dble(pointsPer)))
    ! Take the union with any current lattice. This is done in exact integer arithmetic, and after the hard limits have been
    ! applied so that a lattice already in use is never shrunk.
    if (present(latticeCurrent)) then
       if (latticeCurrent%isDefined()) then
          if (latticeCurrent%scheme    /= scheme   ) call Error_Report('current lattice uses a different gridding scheme'   //{introspection:location})
          if (latticeCurrent%pointsPer /= pointsPer) call Error_Report('current lattice uses a different density of points' //{introspection:location})
          indexMinimum=min(indexMinimum,latticeCurrent%indexMinimum  )
          indexMaximum=max(indexMaximum,latticeCurrent%indexMaximum())
       end if
    end if
    if (indexMaximum <= indexMinimum) call Error_Report('pinned range contains fewer than two points - check any imposed limits'//{introspection:location})
    ! Construct the lattice.
    lattice%scheme      =scheme
    lattice%pointsPer   =pointsPer
    lattice%indexMinimum=indexMinimum
    lattice%count       =indexMaximum-indexMinimum+1
    return
  end function Range_Pinned_Array

  logical function rangeLatticeIsDefined(self)
    !!{RST
    Return true if the object specifies a usable lattice.
    !!}
    implicit none
    class(rangeLattice), intent(in   ) :: self

    rangeLatticeIsDefined=self%scheme /= gridSchemeUndefined .and. self%pointsPer > 0 .and. self%count > 1
    return
  end function rangeLatticeIsDefined

  integer function rangeLatticeIndexMaximum(self)
    !!{RST
    Return the lattice index of the final point of the lattice range.
    !!}
    implicit none
    class(rangeLattice), intent(in   ) :: self

    rangeLatticeIndexMaximum=self%indexMinimum+self%count-1
    return
  end function rangeLatticeIndexMaximum

  double precision function rangeLatticeMinimum(self)
    !!{RST
    Return the value of the first point of the lattice range.
    !!}
    implicit none
    class(rangeLattice), intent(in   ) :: self

    rangeLatticeMinimum=Lattice_Value(dble(self%indexMinimum  )/dble(self%pointsPer),self%scheme)
    return
  end function rangeLatticeMinimum

  double precision function rangeLatticeMaximum(self)
    !!{RST
    Return the value of the final point of the lattice range.
    !!}
    implicit none
    class(rangeLattice), intent(in   ) :: self

    rangeLatticeMaximum=Lattice_Value(dble(self%indexMaximum())/dble(self%pointsPer),self%scheme)
    return
  end function rangeLatticeMaximum

  function rangeLatticeValues(self) result(values)
    !!{RST
    Return the values of all points of the lattice range. The values are computed from the integer lattice indices, and so are
    bit-identical to those of any other range on the same lattice.
    !!}
    implicit none
    class           (rangeLattice), intent(in   )         :: self
    double precision              , dimension(self%count) :: values
    integer                                                 :: i

    do i=1,self%count
       values(i)=Lattice_Value(dble(self%indexMinimum+i-1)/dble(self%pointsPer),self%scheme)
    end do
    return
  end function rangeLatticeValues

  function rangeLatticeValuesLogarithmic(self) result(values)
    !!{RST
    Return the natural logarithms of the values of all points of the lattice range. As for ``rangeLatticeValues``
    the values are computed from the integer lattice indices. Available only for the logarithmic gridding schemes.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (rangeLattice), intent(in   )         :: self
    double precision              , dimension(self%count) :: values
    double precision                                      :: logarithmBase
    integer                                               :: i

    select case (self%scheme%ID)
    case (gridSchemePerDecade%ID)
       logarithmBase=log(10.0d0)
    case (gridSchemePerOctave%ID)
       logarithmBase=log( 2.0d0)
    case default
       logarithmBase=0.0d0
       call Error_Report('logarithmic lattice values are defined only for logarithmic gridding schemes'//{introspection:location})
    end select
    do i=1,self%count
       values(i)=dble(self%indexMinimum+i-1)/dble(self%pointsPer)*logarithmBase
    end do
    return
  end function rangeLatticeValuesLogarithmic

  double precision function rangeLatticeStep(self)
    !!{RST
    Return the spacing of the points of the lattice. Available only for the ``perUnit`` gridding scheme, for which the points
    are uniformly spaced in value; for the logarithmic schemes see ``rangeLatticeStepLogarithmic``.

    Note that this is deliberately *not* evaluated as the difference of two neighbouring lattice points: that difference varies
    in its final bits with position along the lattice, so a table which used it would change its interpolation by of order one
    unit in the last place when extended. Evaluated as below the spacing is a pure function of ``pointsPer``, and so is
    invariant under extension.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(rangeLattice), intent(in   ) :: self

    rangeLatticeStep=0.0d0
    select case (self%scheme%ID)
    case (gridSchemePerUnit%ID)
       rangeLatticeStep=1.0d0/dble(self%pointsPer)
    case default
       call Error_Report('the points of a logarithmically-spaced lattice are not uniformly spaced in value'//{introspection:location})
    end select
    return
  end function rangeLatticeStep

  double precision function rangeLatticeStepLogarithmic(self)
    !!{RST
    Return the spacing of the natural logarithms of the points of the lattice. Available only for the logarithmic gridding
    schemes. As for ``rangeLatticeStep`` this is evaluated as a pure function of ``pointsPer`` so that it is invariant under
    extension of the range.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(rangeLattice), intent(in   ) :: self

    rangeLatticeStepLogarithmic=0.0d0
    select case (self%scheme%ID)
    case (gridSchemePerDecade%ID)
       rangeLatticeStepLogarithmic=log(10.0d0)/dble(self%pointsPer)
    case (gridSchemePerOctave%ID)
       rangeLatticeStepLogarithmic=log( 2.0d0)/dble(self%pointsPer)
    case default
       call Error_Report('logarithmic spacing is defined only for logarithmic gridding schemes'//{introspection:location})
    end select
    return
  end function rangeLatticeStepLogarithmic

  logical function rangeLatticeCovers(self,lattice)
    !!{RST
    Return true if ``lattice`` lies on the same lattice as ``self``, and every one of its points is also a point of ``self``.
    !!}
    implicit none
    class(rangeLattice), intent(in   ) :: self
    type (rangeLattice), intent(in   ) :: lattice

    rangeLatticeCovers=  self   %isDefined   ()                       &
         &             .and.                                          &
         &               lattice%isDefined   ()                       &
         &             .and.                                          &
         &               self   %scheme        == lattice%scheme      &
         &             .and.                                          &
         &               self   %pointsPer     == lattice%pointsPer   &
         &             .and.                                          &
         &               self   %indexMinimum  <= lattice%indexMinimum &
         &             .and.                                          &
         &               self   %indexMaximum() >= lattice%indexMaximum()
    return
  end function rangeLatticeCovers

  double precision function Lattice_Coordinate(value,scheme)
    !!{RST
    Return the (real-valued) lattice coordinate corresponding to ``value``, such that the lattice points are those at which the
    coordinate is an integer multiple of :math:`1/N`.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                           , intent(in   ) :: value
    type            (enumerationGridSchemeType), intent(in   ) :: scheme

    Lattice_Coordinate=0.0d0
    select case (scheme%ID)
    case (gridSchemePerDecade%ID)
       if (value <= 0.0d0) call Error_Report('values must be positive for logarithmic gridding schemes'//{introspection:location})
       Lattice_Coordinate=log10(value)
    case (gridSchemePerOctave%ID)
       if (value <= 0.0d0) call Error_Report('values must be positive for logarithmic gridding schemes'//{introspection:location})
       Lattice_Coordinate=log  (value)/log(2.0d0)
    case (gridSchemePerUnit  %ID)
       Lattice_Coordinate=      value
    case default
       call Error_Report('unrecognized gridding scheme'//{introspection:location})
    end select
    return
  end function Lattice_Coordinate

  double precision function Lattice_Value(coordinate,scheme)
    !!{RST
    Return the value corresponding to the lattice coordinate ``coordinate``---the inverse of ``Lattice_Coordinate``.

    This function is marked ``noinline``, which is essential to the guarantee that all lattice points are bit-identical between
    ranges built on the same lattice. If it is inlined into an array-filling loop the compiler is free to vectorize the
    exponentiation, and a vectorized ``pow()`` need not return the same result as the scalar one---so the same lattice point
    could differ in its final bits between one table and another purely because the two loops had different trip counts. Keeping
    this an opaque call forces every lattice point, everywhere, through the same scalar evaluation. The cost is negligible:
    lattices are built once, and contain at most a few thousand points.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                           , intent(in   ) :: coordinate
    type            (enumerationGridSchemeType), intent(in   ) :: scheme
!GCC$ ATTRIBUTES noinline :: Lattice_Value

    Lattice_Value=0.0d0
    select case (scheme%ID)
    case (gridSchemePerDecade%ID)
       Lattice_Value=10.0d0**coordinate
    case (gridSchemePerOctave%ID)
       Lattice_Value= 2.0d0**coordinate
    case (gridSchemePerUnit  %ID)
       Lattice_Value=        coordinate
    case default
       call Error_Report('unrecognized gridding scheme'//{introspection:location})
    end select
    return
  end function Lattice_Value

  integer function Lattice_Floor(coordinate)
    !!{RST
    Return the largest integer not exceeding ``coordinate``, allowing a small tolerance so that coordinates which should be
    exactly integral (but which are not, due to round-off in the logarithm) are not rounded down by a whole step.
    !!}
    implicit none
    double precision, intent(in   ) :: coordinate

    Lattice_Floor=floor(coordinate+toleranceLattice)
    return
  end function Lattice_Floor

  integer function Lattice_Ceiling(coordinate)
    !!{RST
    Return the smallest integer not less than ``coordinate``, allowing a small tolerance---see ``Lattice_Floor``.
    !!}
    implicit none
    double precision, intent(in   ) :: coordinate

    Lattice_Ceiling=ceiling(coordinate-toleranceLattice)
    return
  end function Lattice_Ceiling

end module Numerical_Ranges
