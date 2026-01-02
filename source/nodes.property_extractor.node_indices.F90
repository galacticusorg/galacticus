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
Implements a property extractor for basic node indices.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorNodeIndices">
   <description>
    A node property extract which extracts various indices related to the merger tree structure:
    \begin{description}
     \item [{\normalfont \ttfamily nodeIndex}] A unique\footnote{Node indices are typically unique, but there is no actual
     requirement within \protect\glc\ that this must be the case. A merger tree construction method could create nodes with
     non-unique indices.} (within a tree) integer index identifying the node;
     \item [{\normalfont \ttfamily parentIndex}] The index of this node's parent node (or $-1$ if it has no parent);
     \item [{\normalfont \ttfamily siblingIndex}] The index of this node's sibling node (or $-1$ if it has no sibling);
     \item [{\normalfont \ttfamily satelliteIndex}] The index of this node's first satellite node (or $-1$ if it has no
     satellites);
     \item [{\normalfont \ttfamily nodeIsIsolated}] Will be $0$ for a node which is a subhalo inside some other node (i.e. a
     satellite galaxy\index{satellite galaxies!identifying}) or $1$ for a node that is an isolated halo (i.e. a central
     galaxy\index{central galaxies!identifying}).
    \end{description}
    The {\normalfont \ttfamily nodeIndex} property corresponds by default to the index of the node in the original merger
    tree. This means that as a galaxy evolves through the tree and, in particular, gets promoted into a new halo the index
    associated with a galaxy will change. This is useful to identify where the galaxy resides in the original (unevolved) tree
    structure, but does not allow galaxies to be traced from one output to the next using their {\normalfont \ttfamily
    nodeIndex} value. By use of the node operator {\normalfont \ttfamily \textless nodeOperator value="indexShift"/\textgreater} this behavior
    can be changed such that the value of {\normalfont \ttfamily nodeIndex} will reflect the index of the earliest progenitor
    node along the main branch of the current node. As such, this index will remain the same for a given galaxy during its
    evolution\index{galaxies!tracing through
    outputs}\index{galaxies!indices}\index{nodes!indices}\index{indices!nodes}\index{indices!galaxies}. These two alternative
    algorithms for propagating node indices are illustrated in Figure~\ref{fig:NodePromotionIndexAlgorithms}.
    \begin{figure}
     \begin{center}
     \includegraphics[width=140mm]{Diagrams/NodePromotionIndices.pdf}
     \end{center}
     \caption{Illustration of  options for the propagation  of node indices during  node promotion events.  Two identical trees
     (top row) are evolved without (left column) and one with (right column) the node operator {\normalfont \ttfamily \textless nodeOperator value="indexShift"/\textgreater}
      The middle and lower rows indicate the resulting node indices after two stages of tree evolution.}
     \label{fig:NodePromotionIndexAlgorithms}
    \end{figure}
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerTuple) :: nodePropertyExtractorNodeIndices
     !!{
     A property extractor class for basic node indices.
     !!}
     private
   contains
     procedure :: elementCount => nodeIndicesElementCount
     procedure :: extract      => nodeIndicesExtract
     procedure :: names        => nodeIndicesNames
     procedure :: descriptions => nodeIndicesDescriptions
     procedure :: unitsInSI    => nodeIndicesUnitsInSI
  end type nodePropertyExtractorNodeIndices

  interface nodePropertyExtractorNodeIndices
     !!{
     Constructors for the \refClass{nodePropertyExtractorNodeIndices} output analysis class.
     !!}
     module procedure nodeIndicesConstructorParameters
  end interface nodePropertyExtractorNodeIndices

contains

  function nodeIndicesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorNodeIndices} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorNodeIndices)                :: self
    type(inputParameters                 ), intent(inout) :: parameters

    self=nodePropertyExtractorNodeIndices()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nodeIndicesConstructorParameters

  integer function nodeIndicesElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily nodeIndices} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorNodeIndices), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    nodeIndicesElementCount=5
    return
  end function nodeIndicesElementCount

  function nodeIndicesExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily nodeIndices} property extractor.
    !!}
    implicit none
    integer         (kind_int8                       ), dimension(:) , allocatable :: nodeIndicesExtract
    class           (nodePropertyExtractorNodeIndices), intent(inout)              :: self
    type            (treeNode                        ), intent(inout)              :: node
    double precision                                  , intent(in   )              :: time
    type            (multiCounter                    ), intent(inout), optional    :: instance
    !$GLC attributes unused :: self, time, instance

    allocate(nodeIndicesExtract(5))
    nodeIndicesExtract   (1:4)=[                             &
         &                      node               %index(), &
         &                      node%parent        %index(), &
         &                      node%sibling       %index(), &
         &                      node%firstSatellite%index()  &
         &                     ]
    select case (node%isSatellite())
    case (.true. )
       nodeIndicesExtract(5  )=0
    case (.false.)
       nodeIndicesExtract(5  )=1
    end select
    return
  end function nodeIndicesExtract

  subroutine nodeIndicesNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily nodeIndices} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorNodeIndices), intent(inout)                             :: self
    double precision                                  , intent(in   )                             :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(5))
    names(1)=var_str('nodeIndex'     )
    names(2)=var_str('parentIndex'   )
    names(3)=var_str('siblingIndex'  )
    names(4)=var_str('satelliteIndex')
    names(5)=var_str('nodeIsIsolated')
    return
  end subroutine nodeIndicesNames

  subroutine nodeIndicesDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily nodeIndices} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorNodeIndices), intent(inout)                             :: self
    double precision                                  , intent(in   )                             :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(5))
    descriptions(1)=var_str('Tree-unique ID for this node.')
    descriptions(2)=var_str('ID of parent node.'           )
    descriptions(3)=var_str('ID of sibling node.'          )
    descriptions(4)=var_str('ID of first satellite node.'  )
    descriptions(5)=var_str('Is the node isolated (0|1)?'  )
    return
  end subroutine nodeIndicesDescriptions

  function nodeIndicesUnitsInSI(self,time)
    !!{
    Return the units of the last isolated redshift property in the SI system.
    !!}
    implicit none
    double precision                                  , allocatable  , dimension(:) :: nodeIndicesUnitsInSI
    class           (nodePropertyExtractorNodeIndices), intent(inout)               :: self
    double precision                                  , intent(in   )               :: time
    !$GLC attributes unused :: self, time

    allocate(nodeIndicesUnitsInSI(5))
    nodeIndicesUnitsInSI=0.0d0
    return
  end function nodeIndicesUnitsInSI

