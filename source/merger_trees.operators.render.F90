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
Implements a merger tree operator which dumps tree data to a file suitable for 3D rendering.
!!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <mergerTreeOperator name="mergerTreeOperatorRender">
   <description>
    A merger tree operator which outputs data on the structure of a merger tree and its halos useful for rendering the tree as
    a 3-D structure to a file named {\normalfont \ttfamily
    render\_$\langle$treeIndex$\rangle$\_$\langle$outputIndex$\rangle$.hdf5} where $\langle${\normalfont \ttfamily
    treeIndex}$\rangle$ is the index of the tree and $\langle${\normalfont \ttfamily outputIndex}$\rangle$ is an incremental
    counter that tracks the number of outputs for this tree. The output is a simple HDF5 file containing the following
    datasets:
    \begin{description}
     \item [{\normalfont \ttfamily nodeIndex}] Index of the node;
     \item [{\normalfont \ttfamily parentIndex}] Index of the parent node;
     \item [{\normalfont \ttfamily childIndex}] Index of the child node;
     \item [{\normalfont \ttfamily time}] Time of the node;
     \item [{\normalfont \ttfamily expansionFactor}] Corresponding expansion factor;
     \item [{\normalfont \ttfamily radiusVirial}] Virial radius of the node;
     \item [{\normalfont \ttfamily position}] $(x,y,z)$ position of the node.
    \end{description}
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorRender
     !!{
     A merger tree operator which dumps tree data to a file suitable for 3D rendering.
     !!}
     private
     class  (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class  (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     integer(kind=kind_int8          )          :: treeIndexPrevious
     integer                                    :: outputCounter
   contains
     final     ::                        renderDestructor
     procedure :: operatePreEvolution => renderOperatePreEvolution
  end type mergerTreeOperatorRender

  interface mergerTreeOperatorRender
     !!{
     Constructors for the {\normalfont \ttfamily render} merger tree operator class.
     !!}
     module procedure renderConstructorParameters
     module procedure renderConstructorInternal
  end interface mergerTreeOperatorRender

contains

  function renderConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily render} merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(mergerTreeOperatorRender )                :: self
    type(inputParameters          ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_
    class(cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    
    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=mergerTreeOperatorRender(cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function renderConstructorParameters

  function renderConstructorInternal(cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily render} merger tree operator class.
    !!}
    implicit none
    type (mergerTreeOperatorRender)                        :: self
    class(darkMatterHaloScaleClass), intent(in   ), target :: darkMatterHaloScale_
    class(cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *cosmologyFunctions_"/>
    !!]

    self%treeIndexPrevious=-1_kind_int8
    self%outputCounter    =-1
    return
  end function renderConstructorInternal

  subroutine renderDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily render} merger tree operator class.
    !!}
    implicit none
    type(mergerTreeOperatorRender), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine renderDestructor
  
  subroutine renderOperatePreEvolution(self,tree)
    !!{
    Output the structure of {\normalfont \ttfamily tree}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic      , nodeComponentPosition, &
         &                                          treeNode
    use :: IO_HDF5                         , only : hdf5Object
    use :: Merger_Tree_Walkers             , only : mergerTreeWalkerAllNodes
    use :: Numerical_Constants_Astronomical, only : gigaYear                , megaParsec
    implicit none
    class           (mergerTreeOperatorRender), intent(inout), target         :: self
    type            (mergerTree              ), intent(inout), target         :: tree
    type            (treeNode                )               , pointer        :: node
    class           (nodeComponentBasic      )               , pointer        :: basic
    class           (nodeComponentPosition   )               , pointer        :: position
    type            (mergerTreeWalkerAllNodes)                                :: treeWalker
    integer         (kind=kind_int8          ), allocatable  , dimension(:  ) :: childIndex     , nodeIndex   , &
         &                                                                       parentIndex
    double precision                          , allocatable  , dimension(:  ) :: expansionFactor, radiusVirial, &
         &time
    double precision                          , allocatable  , dimension(:,:) :: position_
    integer                                                                   :: iNode          , nodesInTree
    character       (len=39                  )                                :: fileName
    type            (hdf5Object              )                                :: fileObject     , treeDataset

    ! Reset output incremental counter if this tree is not the same as the previous one.
    if (tree%index /= self%treeIndexPrevious) then
       self%treeIndexPrevious=tree%index
       self%outputCounter    =-1
    end if
    ! Increment the output counter.
    self%outputCounter=self%outputCounter+1
    ! Construct a file name for the output.
    write (fileName,'(a7,i16.16,a1,i10.10,a5)') "render:",tree%index,":",self%outputCounter,".hdf5"
    ! Count the number of nodes in the tree.
    nodesInTree=0
    treeWalker =mergerTreeWalkerAllNodes(tree)
    do while (treeWalker%next(node))
       nodesInTree=nodesInTree+1
    end do
    ! Allocate arrays for temporary storage.
    allocate(nodeIndex      (nodesInTree))
    allocate(parentIndex    (nodesInTree))
    allocate(childIndex     (nodesInTree))
    allocate(time           (nodesInTree))
    allocate(expansionFactor(nodesInTree))
    allocate(radiusVirial   (nodesInTree))
    allocate(position_      (3,nodesInTree))
    ! Populate arrays with data.
    iNode     =0
    treeWalker=mergerTreeWalkerAllNodes(tree)
     do while (treeWalker%next(node))
       iNode=iNode+1
       basic                    => node               %basic   ()
       position                 => node               %position()
       nodeIndex      (  iNode) =  node               %index   ()
       parentIndex    (  iNode) =  node    %parent    %index   ()
       childIndex     (  iNode) =  node    %firstChild%index   ()
       time           (  iNode) =                                                basic%time()
       expansionFactor(  iNode) =  self    %cosmologyFunctions_ %expansionFactor(basic%time())
       radiusVirial   (  iNode) =  self    %darkMatterHaloScale_%radiusVirial   (node        )
       position_      (:,iNode) =  position                     %position       (            )
    end do
    ! Open an HDF5 file.
    call fileObject%openFile(fileName,overWrite=.true.,objectsOverwritable=.true.)
    ! Write the datasets.
    call fileObject %writeDataset  (nodeIndex      ,"nodeIndex"      ,"Node index []"                                  )
    call fileObject %writeDataset  (parentIndex    ,"parentIndex"    ,"Parent index []"                                )
    call fileObject %writeDataset  (childIndex     ,"childIndex"     ,"Child index []"                                 )
    call fileObject %writeDataset  (expansionFactor,"expansionFactor","Expansion factor []"                            )
    call fileObject %writeDataset  (time           ,"time"           ,"Time [Gyr]"         ,datasetReturned=treeDataset)
    call treeDataset%writeAttribute(gigaYear       ,"unitsInSI"                                                        )
    call treeDataset%close         (                                                                                   )
    call fileObject %writeDataset  (radiusVirial   ,"radiusVirial"   ,"Virial radius [Mpc]",datasetReturned=treeDataset)
    call treeDataset%writeAttribute(megaParsec     ,"unitsInSI"                                                        )
    call treeDataset%close         (                                                                                   )
    call fileObject %writeDataset  (position_      ,"position"       ,"Position [Mpc]"     ,datasetReturned=treeDataset)
    call treeDataset%writeAttribute(megaParsec     ,"unitsInSI"                                                        )
    call treeDataset%close         (                                                                                   )
    call fileObject %close         (                                                                                   )
    ! Deallocate temporary arrays.
    deallocate(nodeIndex      )
    deallocate(parentIndex    )
    deallocate(childIndex     )
    deallocate(time           )
    deallocate(expansionFactor)
    deallocate(radiusVirial   )
    deallocate(position_      )
    return
  end subroutine renderOperatePreEvolution
