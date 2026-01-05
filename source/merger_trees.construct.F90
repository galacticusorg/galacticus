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

!!{
Contains a module which constructs/destructs merger trees.
!!}

module Merger_Tree_Construction
  !!{
  Constructs/destructs merger trees.
  !!}
  use            :: Galacticus_Nodes, only : mergerTree
  use, intrinsic :: ISO_C_Binding   , only : c_size_t
  private

  !![
  <functionClass>
   <name>mergerTreeConstructor</name>
   <descriptiveName>Merger Tree Constructors</descriptiveName>
   <description>
    Class providing merger tree constructors. Here, ``construct'' means any process of creating a representation of a merger
    tree within \glc.
   </description>
   <default>build</default>
   <method name="construct" >
    <description>Construct the merger tree corresponding to the given {\normalfont \ttfamily treeNumber}.</description>
    <type>type(mergerTree), pointer</type>
    <pass>yes</pass>
    <argument>integer(c_size_t), intent(in   ) :: treeNumber</argument>
    <argument>logical          , intent(  out) :: finished</argument>
   </method>
   <method name="randomSequenceNonDeterministicWarn" >
    <description>Display a warning if the merger tree random number generator sequence is non-deterministic.</description>
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
    <type>void</type>
    <pass>yes</pass>
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
