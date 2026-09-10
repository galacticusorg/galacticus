.. _manual-sec-growingATabulation:

Growing a Tabulation
====================

Many quantities in Galacticus are expensive enough that they are tabulated once
and interpolated thereafter. The range which must be tabulated is usually not
known in advance: it is discovered from the values actually requested as a model
runs. A tabulation therefore has to *grow*, and the question this page answers
is how to grow one without changing the values it has already returned.

The mechanism is an **absolute lattice**: a set of abscissae fixed once and for
all by the gridding scheme, entirely independent of which range any particular
table happens to span. A table built on such a lattice can be extended, and
every point it already held---and every value interpolated between those
points---is unchanged, bit for bit. The pieces are ``Range_Pinned`` and
``rangeLattice`` (in ``Numerical_Ranges``), the ``extend`` methods of the
``table1D`` and ``table2DLogLogLin`` classes, and, for tabulations expensive
enough to keep between runs, ``Table_Caches``.

Why a naively-grown tabulation is a problem
-------------------------------------------

The natural way to grow a table is to widen its bounds until they enclose the
new request, redistribute the requested number of points over the widened range,
and paste the previously computed values back in. This does not work, and the
reason is worth understanding before reaching for the machinery which avoids it.

The bounds so obtained depend on the *history* of the requests. They are grown
by repeated multiplication by some factor from an origin set by whichever value
happened to be asked for first, so two runs which end up covering the same range
of values, but reach it in a different order, tabulate at different abscissae.
Nothing about the model has changed, yet every interpolated value differs.

Worse, the abscissae pasted back in no longer belong to the table which now
holds them. Redistributing :math:`N` points over the widened range gives a new
spacing, but the preserved values retain the abscissae they were computed at, so
the table becomes a patchwork of spacings which no single interpolation factor
describes. Where such a table is indexed arithmetically---from a stored inverse
spacing, rather than by searching---the index computed drifts further from the
truth with each extension. This is not hypothetical: it is exactly the failure
which the expansion-factor tabulation in
``source/cosmology/functions/matter_lambda.F90`` used to suffer, where after
sufficiently many extensions the computed index fell beyond the end of the
table.

The absolute lattice
--------------------

A lattice is specified by a gridding scheme and a density of points,
``pointsPer`` :math:`=N`. Its points are the :math:`x_j`, for *all* integers
:math:`j`:

* ``gridSchemePerDecade`` --- :math:`x_j = 10^{j/N}`;
* ``gridSchemePerOctave`` --- :math:`x_j = 2^{j/N}`;
* ``gridSchemePerUnit`` --- :math:`x_j = j/N`.

Since :math:`x_j` is a function of :math:`j` and :math:`N` alone, the lattice
does not know, and cannot depend on, the range over which any particular table
is built. Two tables built on the same lattice therefore have identical
abscissae wherever their ranges overlap.

A ``rangeLattice`` object names a contiguous stretch of such a lattice: the
``count`` points with integer indices ``indexMinimum`` through
``indexMinimum+count-1``. It carries no abscissae of its own---they are
regenerated from the integer indices on demand by ``values()`` (or
``valuesLogarithmic()``), so they are bit-identical to those of any other range
on the same lattice. Its useful methods are:

* ``isDefined()`` --- is this a usable lattice? A default-initialized
  ``rangeLattice`` is not, and that is how "nothing is tabulated yet" is
  represented;
* ``covers(lattice)`` --- is ``lattice`` on the same lattice as this one, with
  every one of its points also a point of this one? This is the test for "is the
  existing tabulation already sufficient?";
* ``minimum()``, ``maximum()``, ``indexMaximum()``, ``values()``,
  ``valuesLogarithmic()``;
* ``step()`` (``perUnit`` only) and ``stepLogarithmic()`` (logarithmic schemes
  only) --- the spacing of the points.

``step()`` and ``stepLogarithmic()`` are evaluated as pure functions of
``pointsPer``, and deliberately *not* as the difference of two neighboring
lattice points: that difference varies in its final bits with position along the
lattice, so a table which derived its interpolation factor from it would change
every interpolated value by of order one unit in the last place whenever its
lower bound moved---defeating the very reuse the lattice exists to provide.
**Always take the spacing from the lattice.**

Pinning a range to the lattice
------------------------------

``Range_Pinned`` turns a request---one value, or an array of values, which the
tabulation must cover---into the ``rangeLattice`` to tabulate:

.. code-block:: fortran

   use :: Numerical_Ranges, only : Range_Pinned, rangeLattice, gridSchemePerDecade
   ...
   type(rangeLattice) :: lattice
   ...
   lattice=Range_Pinned(                                    &
        &                              time               , &
        &                              pointsPerDecade    , &
        &                              gridSchemePerDecade, &
        &               marginFactor  =2.0d0              , &
        &               anchorEvery   =pointsPerDecade/2  , &
        &               latticeCurrent=self%lattice         &
        &              )

The bounds are found by, in order:

#. **applying the safety margin** to the range of the requested values: a
   multiplicative ``marginFactor`` (default :math:`2`) at each end for the
   logarithmic schemes, or an additive ``marginOffset`` (default :math:`1`) for
   ``perUnit``. The margin means that a request which creeps outward by a little
   does not force a retabulation for every step it takes;
#. **taking the union with** ``rangeCurrent``, if given---a pair of values which
   the range must span, used where some bound must always be included but is not
   itself known to lie on the lattice;
#. **pinning the bounds outward** to lattice points whose index is a multiple of
   ``anchorEvery``;
#. **clamping to the hard limits** ``limitMinimum`` and ``limitMaximum``, if
   given, with the clamped edge snapped *inward* to the lattice so that all
   points remain lattice points;
#. **taking the union with** ``latticeCurrent``, if given---in exact integer
   arithmetic, and after the clamping, so that a lattice already in use is never
   shrunk and extension of a table built on it can never fail.

``anchorEvery``
~~~~~~~~~~~~~~~

``anchorEvery`` is the interval, in lattice steps, to which the bounds are
pinned. It controls only how far beyond the request the table extends, and is
independent of ``pointsPer``, which controls the accuracy of interpolation
within it. It defaults to ``pointsPer``---whole decades, octaves, or unit
intervals---which gives the coarsest set of possible ranges and hence the
strongest determinism. A value of :math:`1` pins to the lattice points
themselves.

Intermediate values are for cases where a whole decade of extra tabulation is
extravagant. On a cosmic time axis a decade spans most of the history of the
universe, so ``anchorEvery=pointsPer/2`` (half decades, a factor of
:math:`\sqrt{10}`) is the usual choice there; where each additional point costs
an ODE integration, ``anchorEvery=1`` may be the only affordable choice.
Determinism, exact reuse on extension, and the mergeability of cached tables
hold for *any* choice---only the granularity of the set of possible ranges
changes.

``rangeCurrent`` versus ``latticeCurrent``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These two look similar and are not interchangeable. The distinction is the most
common source of error in using ``Range_Pinned``, so it is worth stating
plainly:

* the **request is the target**, and only the target receives the safety margin;
* the range **already tabulated** is taken in union through ``latticeCurrent``,
  after the margin and after the hard limits.

Folding the current bounds into the target instead---taking a ``min``/``max``
against them before calling ``Range_Pinned``---applies the safety margin to a
bound to which the margin has already been applied once. The range then ratchets
outward by another factor on every retabulation, growing without limit.
``latticeCurrent`` exists precisely so that the current range can be honored
without the margin being applied to it a second time.

``rangeCurrent`` is for the different case of a bound which must always be
spanned but does not lie on the lattice: the present day, for a tabulation of
integrals evaluated to the present day, or the current range of a table which
was not itself built on a lattice. It is applied *after* the margin, so it
raises the upper bound without also lowering the lower one.

Deciding whether to retabulate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the pinned lattice, not the requested value against the current bounds:

.. code-block:: fortran

   lattice=Range_Pinned(x,pointsPer,scheme,latticeCurrent=self%table%lattice)
   if (self%table%lattice%covers(lattice)) return

Testing the request against the bounds would apply the safety margin only when
the request happened to arrive outside them, so the range finally reached would
depend on the order in which values were asked for---the very dependence pinning
exists to remove. Pinning first and comparing lattices costs one cheap integer
comparison and has no such history.

Extending the tabulation
------------------------

``extend`` grows a table onto a new lattice, creating it if it does not yet
exist, and returns an ``isComputed`` mask which is true at those points whose
values were carried over and false at those the caller must now evaluate:

.. code-block:: fortran

   logical, allocatable, dimension(:) :: isComputed
   ...
   call self%table%extend(lattice,isComputed)
   do i=1,lattice%count
      if (isComputed(i)) cycle
      call self%table%populate(expensiveFunction(self%table%x(i)),i)
   end do

Two things to note. Take the abscissae from the table (or from
``lattice%values()``), never recompute them independently, so that what is
evaluated is exactly what will later be interpolated. And skipping the computed
points is the entire purpose of the exercise---the mask is not advisory.

``extend`` and ``x`` are methods of the abstract ``table1D``, but ``populate``
is defined on the concrete table types. Where the table is held as a
``class(table1D), allocatable``---as it must be to use ``Table_Caches``---the
loop therefore needs a ``select type`` around it:

.. code-block:: fortran

   select type (table_ => self%table)
   type is (table1DLogarithmicLinear)
      do i=1,lattice%count
         if (isComputed(i)) cycle
         call table_%populate(expensiveFunction(table_%x(i)),i)
      end do
   end select

The new lattice must be commensurate with (same scheme, same ``pointsPer``) and
must *contain* the lattice on which the table is currently tabulated; anything
else is a fatal error. ``Range_Pinned`` guarantees this when it is passed the
table's current lattice as ``latticeCurrent``.

``extend`` is supported for the table types with uniformly-spaced *internal*
abscissae---the linear and logarithmic types, with linear or cubic-spline
interpolation---and is an error for the others. The two-dimensional
``table2DLogLogLin`` takes a lattice per axis, each pinned independently:

.. code-block:: fortran

   call self%table%extend(latticeX,latticeY,isComputed)

Here the preserved values form a rectangular block, so ``isComputed`` is true
over that block and false over the L-shaped remainder. Both axes are extended in
one call because the block preserved is the rectangle common to the two
tabulations.

Tabulations held in raw arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not every tabulation lives in a ``table`` object. For a plain rank-1 allocatable
array, ``Range_Lattice_Extend`` does the whole job:

.. code-block:: fortran

   use :: Numerical_Ranges, only : Range_Lattice_Extend
   ...
   call Range_Lattice_Extend(self%lattice,latticeNew,self%values,isComputed)
   self%lattice=latticeNew

``Range_Lattice_Extend_2D`` does the same for a rank-2 array, extending both
axes at once. Where the tabulation is of some other rank, or the tabulated axis
is not the first, use ``Range_Lattice_Offset``: it returns the index offset at
which the points of the current lattice appear within the new one, and leaves
the caller to move its own arrays---typically a ``Move_Alloc`` followed by a
single array-section assignment.

What is, and is not, preserved
------------------------------

The reuse guarantee is strong, and it is worth being precise about its scope.

**Preserved bit-for-bit.** The abscissae, which are regenerated from integer
lattice indices rather than carried over, and the tabulated values, which are
copied into the index window they occupy in the new table. That window is
computed from the integer lattice indices, so no searching or floating-point
comparison of abscissae is involved. For a linearly-interpolated table the
interpolation is local and the spacing comes from the lattice, so every
interpolated value on the overlap is unchanged too.

This matters to any consumer which integrates across a range of the tabulated
variable, or compares quantities evaluated either side of a retabulation. Absent
the guarantee, such a consumer sees the tabulation shift beneath it: a quantity
computed before the table grew and one computed after are no longer strictly
comparable, and differences which ought to be zero are not.

**Not preserved: cubic-spline coefficients.** A cubic spline is not local---every
coefficient depends on every tabulated value---so adding points at either end
changes the interpolant everywhere. The coefficients are therefore discarded and
rebuilt. A cubic-spline table preserves its tabulated values exactly, but *not*
the values it interpolates between them. Callers must accept this; where it is
unacceptable, use a linearly-interpolated table.

A note on why ``Lattice_Value``, which converts a lattice coordinate back to a
value, is marked ``noinline``: were it inlined into an array-filling loop the
compiler would be free to vectorize the exponentiation, and a vectorized
``pow()`` need not agree in its final bits with the scalar one. The same lattice
point could then differ between two tables merely because the loops filling them
had different trip counts. Keeping it an opaque call forces every lattice point,
everywhere, through the same scalar evaluation.

When carry-over is not valid
----------------------------

Carrying values over is correct only when each tabulated value depends on its
own abscissa and on nothing else about the range. That holds for a function
evaluated pointwise, and for an integral whose limits do not move---an integral
from each tabulated time to the present day, for instance, is unaffected when
the earliest time tabulated changes.

It does **not** hold when the values are generated by marching across the table.
The expansion-factor tabulation in ``matter_lambda.F90`` is the worked example:
the expansion factor is integrated forward from an initial condition imposed at
the earliest tabulated epoch, so every value in the table depends on where that
integration began. Extend the table downward and the initial condition moves;
every value which was carried over is then the solution of a different initial
value problem from the one the new table represents.

The idiom for such a case is to discard everything by presenting an *undefined*
lattice to the extension, which allocates afresh with ``isComputed`` false
throughout:

.. code-block:: fortran

   carryOver=self%initialized .and. self%lattice%isDefined()
   if (carryOver) carryOver=latticeNew%indexMinimum == self%lattice%indexMinimum
   if (.not.carryOver) then
      if (allocated(self%values)) deallocate(self%values)
      self%lattice=rangeLattice()
   end if
   call Range_Lattice_Extend(self%lattice,latticeNew,self%values,isComputed)

Note that the test is on the integer lattice index of the lower bound, not on
the value of the bound---exact, and free of any tolerance.

Where a tabulation is of this kind it is worth seeding its lower bound to a
fixed value rather than letting it follow the request, so that the bound moves
rarely, if ever, and the discard is rarely triggered.

Persistent caches
-----------------

Where a tabulation is expensive enough to be worth keeping between runs,
``Table_Caches`` (``source/objects/tables/caches.F90``) writes it to a file
under ``pathTypeDataDynamic``. Because the table is built on an absolute
lattice, a cached file whose range does not cover a new request need not be
discarded: the cached abscissae are identical to the new ones wherever the two
overlap, so the file can be *merged* with the table in memory, only the
genuinely missing points computed, and the union written back. Cache files
therefore grow monotonically toward the union of every range ever requested,
rather than thrashing between runs.

The file name comes from ``Table_Cache_File_Name``:

.. code-block:: fortran

   self%fileName   =Table_Cache_File_Name(                                                                                   &
        &                                 subDirectory    ='intergalacticMedium'                                           , &
        &                                 objectType      =char(self%objectType      (                                   )), &
        &                                 hashedDescriptor=char(self%hashedDescriptor(includeSourceDigest         =.true.  , &
        &                                                                             includeFileModificationTimes=.true.))  &
        &                                )
   call Directory_Make(File_Path(self%fileName))

Passing ``includeSourceDigest=.true.`` and
``includeFileModificationTimes=.true.`` makes the name change whenever either
the parameters of the object or the source code which generates the tabulation
changes---so a stale cache is never read back. Build the name with this function
rather than by hand: doing it by hand is how a discriminator comes to be
omitted, and two different tabulations come to share a file.

The full idiom is then:

.. code-block:: fortran

   lattice=Range_Pinned(x,pointsPer,scheme,latticeCurrent=self%table%lattice)
   if (self%table%lattice%covers(lattice)) return
   ! Merge in any tabulation already cached on disk, then re-evaluate the range required in the light of it.
   call Table_Cache_Restore(self%table,self%fileName,status)
   lattice=Range_Pinned(x,pointsPer,scheme,latticeCurrent=self%table%lattice)
   call self%table%extend(lattice,isComputed)
   if (any(.not.isComputed)) then
      ! Compute only the missing points. Note that no lock is held while doing so.
      select type (table_ => self%table)
      type is (table1DLogarithmicLinear)
         do i=1,lattice%count
            if (isComputed(i)) cycle
            call table_%populate(expensiveFunction(table_%x(i)),i)
         end do
      end select
      ! Merge with, and write back, the cache file.
      call Table_Cache_Store(self%table,self%fileName)
   end if

``Range_Pinned`` is called twice, and must be: restoring the cache may widen the
table's lattice, and the range finally tabulated must be the union of that with
the request.

Locking is handled internally, and the important property is what is *not*
locked. ``Table_Cache_Restore`` takes a shared lock, and ``Table_Cache_Store``
an exclusive one, but each holds it only for the duration of the file
input/output---never across the tabulation itself, which is the expensive part.
``Table_Cache_Store`` re-reads the file under its exclusive lock before writing,
so a range written by another process since our own read is merged rather than
overwritten. Note that ``Table_Cache_Store`` may extend the table, since it
returns holding the union of the cached and in-memory tabulations.

Merging in two dimensions is more restrictive than in one: the union of two
rectangles is a rectangle only if one contains the other, so unless the cached
tabulation covers that in memory (in which case it is adopted) or is covered by
it (in which case there is nothing to do), the cached tabulation is reported and
ignored.

Checklist
---------

When adding or converting a tabulation:

#. Store a ``rangeLattice`` alongside the tabulation, or use the ``lattice``
   component the ``table`` classes already provide.
#. Pin the range with ``Range_Pinned``, passing the *request* as the target and
   the current lattice as ``latticeCurrent``. Never fold the current bounds into
   the target.
#. Choose ``anchorEvery`` for the cost of over-tabulating, and ``pointsPer`` for
   the accuracy of interpolation. They are independent.
#. Return early when ``covers`` says the existing tabulation is sufficient.
#. Extend, and evaluate only the points for which ``isComputed`` is false.
#. Take abscissae, and any interpolation spacing, from the lattice or the table.
#. Check that carry-over is actually valid for the quantity being tabulated; if
   it is not, discard explicitly by resetting the lattice.
#. If the tabulation is expensive enough to keep between runs, wrap it in
   ``Table_Cache_Restore``/``Table_Cache_Store``.

The unit tests in ``source/tests/make_ranges.F90`` and
``source/tests/table_caches.F90`` exercise all of this, and are a useful source
of small, concrete examples.
