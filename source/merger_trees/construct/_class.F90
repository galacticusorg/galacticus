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
Contains a module which constructs/destructs merger trees.
!!}

module Merger_Tree_Construction
  !!{RST
  Constructs/destructs merger trees.
  !!}
  use            :: Galacticus_Nodes, only : mergerTree
  use, intrinsic :: ISO_C_Binding   , only : c_size_t
  private

  !![
  <functionClass docformat="rst">
   <name>mergerTreeConstructor</name>
   <descriptiveName>Merger Tree Constructors</descriptiveName>
   <description>
   Class providing merger tree constructors---any process by which a representation of a dark matter halo merger tree is created and made available within Galacticus for subsequent galaxy formation calculations. Implementations include stochastic Monte Carlo tree building algorithms (which generate trees on-the-fly from excursion set merger rates), reading pre-computed trees from N-body simulations stored in HDF5 or other formats, and analytic approximations. The ``construct`` method returns a :galacticus-class:`mergerTree` pointer for the given tree index, and signals completion once all trees in the set have been constructed.
   </description>
   <default>build</default>
   <method name="construct" >
    <description>
    Construct the merger tree corresponding to the given ``treeNumber``.
    </description>
    <type>type(mergerTree), pointer</type>
    <pass>yes</pass>
    <argument>integer(c_size_t), intent(in   ) :: treeNumber</argument>
    <argument>logical          , intent(  out) :: finished</argument>
   </method>
   <method name="countTrees" >
    <description>
    Return the total number of trees that will be constructed, or a negative value if this is not known. The default implementation returns a negative value; implementations that know the total tree count up-front (e.g. the ``build`` constructor) should override this to enable whole-run progress and run-time estimation.
    </description>
    <type>integer(c_size_t)</type>
    <pass>yes</pass>
    <code>
      ! By default the total number of trees is not known.
      !$GLC attributes unused :: self
      mergerTreeConstructorCountTrees=-1_c_size_t
    </code>
   </method>
   <method name="treeMasses" >
    <description>
    Return the root masses (in units of $M_\odot$) of the trees that will be constructed. The output array is left unallocated if the masses are not known up-front. The default implementation leaves the array unallocated; implementations that know the tree masses up-front (e.g. the ``build`` constructor) should override this to enable whole-run run-time estimation.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, allocatable, dimension(:), intent(  out) :: masses</argument>
    <code>
      ! By default the tree masses are not known.
      !$GLC attributes unused :: self, masses
    </code>
   </method>
   <method name="treeCountsNodes" >
    <description>
    Return the number of nodes in each of the trees that will be constructed. The output array is left unallocated if the node counts are not known up-front. The default implementation leaves the array unallocated; implementations for which the node count of each tree is known before it is evolved (e.g. the ``read`` constructor) should override this to enable node-count-based whole-run run-time estimation, which is useful where root masses are not known up-front.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer(c_size_t), allocatable, dimension(:), intent(  out) :: countsNodes</argument>
    <code>
      ! By default the tree node counts are not known.
      !$GLC attributes unused :: self, countsNodes
    </code>
   </method>
   <method name="randomSequenceNonDeterministicWarn" >
    <description>
    Display a warning if the merger tree random number generator sequence is non-deterministic.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <modules>
      <name>Display</name>
      <only>displayMessage, displayMagenta, displayReset</only>
    </modules>
    <modules>
      <name>MPI_Utilities</name>
      <only>mpiSelf</only>
    </modules>
    <modules>
      <name>OMP_Lib</name>
      <only>OMP_Get_Max_Threads</only>
    </modules>
    <argument>type(mergerTree), intent(inout) :: tree</argument>
    <code>
     logical, save :: nonDeterministicWarned=.false.

     if (.not.nonDeterministicWarned) then
        !$omp critical (treeRandomSequenceNonDeterministicWarn)
        if (.not.nonDeterministicWarned) then
           if        (                                                                                                                                                                   &amp;
                !$ &amp;   (tree%randomNumberGenerator_%openMPIndependent() .and. OMP_Get_Max_Threads() > 1)                                                                             &amp;
                !$ &amp;  .or.                                                                                                                                                           &amp;
                &amp;      (tree%randomNumberGenerator_%   mpiIndependent() .and. mpiSelf%count      () > 1)                                                                             &amp;
                &amp;    ) call displayMessage(                                                                                                                                          &amp;
                &amp;                          displayMagenta()//'WARNING:'//displayReset()//                                                                                            &amp;
                &amp;                          " per-tree random number sequences may not be deterministic - see:"                                                          //char(10)// &amp;
                &amp;                          "        https://github.com/galacticusorg/galacticus/wiki/Troubleshooting#non-deterministic-per-tree-random-number-sequences"             &amp;
                &amp;                          )
           nonDeterministicWarned=.true.
        end if
        !$omp end critical (treeRandomSequenceNonDeterministicWarn)
     end if
    </code>
   </method>
  </functionClass>
  !!]

end module Merger_Tree_Construction
