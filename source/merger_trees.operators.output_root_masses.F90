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
  Implements a merger tree operator which outputs a file of tree root masses (and weights).
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  ! Buffer size for tree data.
  !![
  <scoping>
   <module variables="outputRootMassesBufferSize"/>
  </scoping>
  !!]
  integer, parameter :: outputRootMassesBufferSize=1000

  !![
  <mergerTreeOperator name="mergerTreeOperatorOutputRootMasses">
   <description>Output a file of tree root masses (and weights).</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorOutputRootMasses
     !!{
     A merger tree operator class which outputs a file of tree root masses (and weights).
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                               :: cosmologyFunctions_         => null()
     integer                                                                          :: nodeHierarchyLevelMaximumID          , treeCount
     double precision                                                                 :: time                                 , redshift
     double precision                         , dimension(outputRootMassesBufferSize) :: mass                                 , weight
     type            (varying_string         )                                        :: fileName
     logical                                                                          :: alwaysIsolatedHalosOnly
   contains
     final     ::                        outputRootMassesDestructor
     procedure :: operatePreEvolution => outputRootMassesOperatePreEvolution
     procedure :: finalize            => outputRootMassesFinalize
  end type mergerTreeOperatorOutputRootMasses

  interface mergerTreeOperatorOutputRootMasses
     !!{
     Constructors for the tree root mass outputting merger tree operator class.
     !!}
     module procedure outputRootMassesConstructorParameters
     module procedure outputRootMassesConstructorInternal
  end interface mergerTreeOperatorOutputRootMasses

contains

  function outputRootMassesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the conditional mass function merger tree operator class which takes a parameter set as input.
    !!}
    implicit none
    type            (mergerTreeOperatorOutputRootMasses)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    double precision                                                    :: time                   , redshift
    logical                                                             :: alwaysIsolatedHalosOnly
    type            (varying_string                    )                :: fileName

    !![
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift at which to gather tree root masses.</description>
    </inputParameter>
    <inputParameter>
      <name>alwaysIsolatedHalosOnly</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>Include only always-isolated halos when gathering tree root masses?</description>
    </inputParameter>
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file to which tree masses should be written.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    ! Get time from redshift.
    time   =cosmologyFunctions_ %cosmicTime                 (          &
         &   cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                    redshift &
         &                                                   )         &
         &                                                  )
    ! Construct the instance.
    self=mergerTreeOperatorOutputRootMasses(time,alwaysIsolatedHalosOnly,fileName,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function outputRootMassesConstructorParameters

  function outputRootMassesConstructorInternal(time,alwaysIsolatedHalosOnly,fileName,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the conditional mass function merger tree operator class.
    !!}
    use :: File_Utilities  , only : File_Exists, File_Remove
    use :: HDF5_Access     , only : hdf5Access
    implicit none
    type            (mergerTreeOperatorOutputRootMasses)                        :: self
    double precision                                    , intent(in   )         :: time
    logical                                             , intent(in   )         :: alwaysIsolatedHalosOnly
    type            (varying_string                    ), intent(in   )         :: fileName
    class           (cosmologyFunctionsClass           ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="time, alwaysIsolatedHalosOnly, fileName, *cosmologyFunctions_"/>
    !!]

    ! Initialize.
    self%treeCount=0
    self%redshift =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    !![
    <addMetaProperty component="basic" name="nodeHierarchyLevelMaximum" id="self%nodeHierarchyLevelMaximumID" type="integer"/>
    !!]
    ! Remove any pre-existing file.
    !$ call hdf5Access%set()
    if (File_Exists(fileName)) call File_Remove(fileName)
    !$ call hdf5Access%unset()
    return
  end function outputRootMassesConstructorInternal

  subroutine outputRootMassesDestructor(self)
    !!{
    Destructor for  the \refClass{mergerTreeOperatorOutputRootMasses} merger tree operator class.
    !!}
    implicit none
    type(mergerTreeOperatorOutputRootMasses), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine outputRootMassesDestructor
  
  subroutine outputRootMassesOperatePreEvolution(self,tree)
    !!{
    Compute conditional mass function on {\normalfont \ttfamily tree}.
    !!}
    use :: Galacticus_Nodes    , only : mergerTree                   , nodeComponentBasic, treeNode
    use :: Merger_Tree_Walkers , only : mergerTreeWalkerIsolatedNodes
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (mergerTreeOperatorOutputRootMasses), target , intent(inout) :: self
    type            (mergerTree                        ), target , intent(inout) :: tree
    type            (treeNode                          ), pointer                :: node             , nodeChild      , &
         &                                                                          nodeSibling
    class           (nodeComponentBasic                ), pointer                :: basic            , basicChild     , &
         &                                                                          basicSibling
    type            (mergerTreeWalkerIsolatedNodes     )                         :: treeWalker
    double precision                                                             :: branchBegin      , branchEnd      , &
         &                                                                          branchMassInitial, branchMassFinal, &
         &                                                                          massRoot

    ! Iterate over nodes, searching for root nodes.
    treeWalker=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       ! Get the child node, and process if child exists.
       if (associated(node%firstChild)) then
          nodeChild => node%firstChild
       else
          nodeChild => node
       end if
       ! Check if child should be included.
       basicChild => nodeChild%basic()
       if     (                                                                                     &
            &   .not.self      %alwaysIsolatedHalosOnly                                             &
            &  .or.                                                                                 &
            &        basicChild%integerRank0MetaPropertyGet (self%nodeHierarchyLevelMaximumID) == 0 &
            & ) then
          ! Determine range of times spanned by this branch.
          basic       => node      %basic()
          branchBegin =  basicChild%time ()
          branchEnd   =  basic     %time ()
          ! Does the branch span the search time?
          if     (                                                     &
               &     branchBegin <= self%time                          &
               &  .and.                                                &
               &   (                                                   &
               &     branchEnd   >  self%time                          &
               &    .or.                                               &
               &     (                                                 &
               &       .not.associated(node%parent)                    &
               &      .and.                                            &
               &       Values_Agree(branchEnd,self%time,relTol=1.0d-6) &
               &     )                                                 &
               &   )                                                   &
               & ) then
             ! Get the masses on the branch.
             branchMassInitial=basicChild%mass()
             branchMassFinal  =basic     %mass()
             ! Remove the mass in any non-primary progenitors - we don't want to include
             ! their mass in the estimated mass growth rate of this node.
             if (associated(node%firstChild)) then
                nodeSibling => node%firstChild%sibling
                do while (associated(nodeSibling))
                   basicSibling    => nodeSibling%basic()
                   branchMassFinal =  branchMassFinal-basicSibling%mass()
                   nodeSibling     => nodeSibling%sibling
                end do
             end if
             ! Do not let the parent mass decrease along the branch.
             branchMassFinal=max(branchMassFinal,branchMassInitial)
             ! Interpolate to get the mass at the required time.
             if (branchEnd == branchBegin) then
                massRoot=branchMassFinal
             else
                massRoot=                 +branchMassInitial  &
                     &   +(branchMassFinal-branchMassInitial) &
                     &   *(self%time      -branchBegin      ) &
                     &   /(branchEnd      -branchBegin      )
             end if
             ! Store the tree data.
             self%treeCount                =self         %treeCount   +1
             self%mass     (self%treeCount)=              massRoot
             self%weight   (self%treeCount)=node%hostTree%volumeWeight
             ! Flush the buffers to file if necessary.
             if (self%treeCount == outputRootMassesBufferSize) call self%finalize()
          end if
       end if
    end do
    return
  end subroutine outputRootMassesOperatePreEvolution

  subroutine outputRootMassesFinalize(self)
    !!{
    Outputs conditional mass function.
    !!}
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5, only : hdf5Object
    implicit none
    class(mergerTreeOperatorOutputRootMasses), intent(inout) :: self
    type (hdf5Object                        ), target        :: outputFile

    ! If the buffers are empty, we have nothing to do.
    if (self%treeCount == 0) return
    ! Open the output file.
    !$ call hdf5Access%set  ()
    outputFile=hdf5Object(char(self%fileName),overWrite=.false.)
    ! Write the data.
    call outputFile%writeDataset(self%mass  (1:self%treeCount),"treeRootMass","Tree root node masses.",appendTo =.true. )
    call outputFile%writeDataset(self%weight(1:self%treeCount),"treeWeight"  ,"Tree weights."         ,appendTo =.true. )
    !$ call hdf5Access%unset()
    ! Reset the buffer counter.
    self%treeCount=0
    return
  end subroutine outputRootMassesFinalize
