!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Contains a module which implements a merger tree operator which outputs a file of tree root masses (and weights).
  
  ! Buffer size for tree data.
  integer, parameter :: outputRootMassesBufferSize=1000
  
  !# <mergerTreeOperator name="mergerTreeOperatorOutputRootMasses" defaultThreadPrivate="yes">
  !#  <description>Output a file of tree root masses (and weights).</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorOutputRootMasses
     !% A merger tree operator class which outputs a file of tree root masses (and weights).
     private
     integer                                                                 :: treeCount
     double precision                                                        :: time
     double precision                , dimension(outputRootMassesBufferSize) :: mass                   , weight
     type            (varying_string)                                        :: fileName
     logical                                                                 :: alwaysIsolatedHalosOnly
   contains
     final     ::             outputRootMassesDestructor
     procedure :: operate  => outputRootMassesOperate
     procedure :: finalize => outputRootMassesFinalize
  end type mergerTreeOperatorOutputRootMasses

  interface mergerTreeOperatorOutputRootMasses
     !% Constructors for the tree root mass outputting merger tree operator class.
     module procedure outputRootMassesConstructorParameters
     module procedure outputRootMassesConstructorInternal
  end interface mergerTreeOperatorOutputRootMasses

contains

  function outputRootMassesConstructorParameters(parameters)
    !% Constructor for the conditional mass function merger tree operator class which takes a parameter set as input.
    use Cosmology_Functions
    implicit none
    type            (mergerTreeOperatorOutputRootMasses)                :: outputRootMassesConstructorParameters
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    double precision                                                    :: time                                 , redshift
    logical                                                             :: alwaysIsolatedHalosOnly
    type            (varying_string                    )                :: fileName
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>redshift</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The redshift at which to gather tree root masses.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alwaysIsolatedHalosOnly</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>Include only always-isolated halos when gathering tree root masses?</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of the file to which tree masses should be written.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    ! Get time from redshift.
    cosmologyFunctions_ => cosmologyFunctions()
    time                =  cosmologyFunctions_ %cosmicTime                 (          &
         &                  cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                                   redshift &
         &                                                                  )         &
         &                                                                 )
    ! Construct the instance.
    outputRootMassesConstructorParameters=outputRootMassesConstructorInternal(time,alwaysIsolatedHalosOnly,fileName)
    return
  end function outputRootMassesConstructorParameters

  function outputRootMassesConstructorInternal(time,alwaysIsolatedHalosOnly,fileName)
    !% Internal constructor for the conditional mass function merger tree operator class.
    use File_Utilities
    use System_Command
    implicit none
    type            (mergerTreeOperatorOutputRootMasses)                :: outputRootMassesConstructorInternal
    double precision                                    , intent(in   ) :: time
    logical                                             , intent(in   ) :: alwaysIsolatedHalosOnly
    type            (varying_string                    ), intent(in   ) :: fileName

    ! Initialize.
    outputRootMassesConstructorInternal%time                   =time
    outputRootMassesConstructorInternal%alwaysIsolatedHalosOnly=alwaysIsolatedHalosOnly
    outputRootMassesConstructorInternal%fileName               =fileName
    outputRootMassesConstructorInternal%treeCount              =0
    ! Remove any pre-existing file.
    !$omp critical (HDF5_Access)
    if (File_Exists(fileName)) call System_Command_Do("rm -f "//fileName)
    !$omp end critical (HDF5_Access)
    return
  end function outputRootMassesConstructorInternal

  elemental subroutine outputRootMassesDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorOutputRootMasses), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine outputRootMassesDestructor
  
  subroutine outputRootMassesOperate(self,tree)
    !% Compute conditional mass function on {\normalfont \ttfamily tree}.
    use Galacticus_Nodes
    use Numerical_Comparison
    use Galacticus_Error
    implicit none
    class           (mergerTreeOperatorOutputRootMasses)         , intent(inout) :: self
    type            (mergerTree                        ), target , intent(inout) :: tree
    type            (treeNode                          ), pointer                :: node             , nodeChild      , &
         &                                                                          nodeSibling
    type            (mergerTree                        ), pointer                :: treeCurrent
    class           (nodeComponentBasic                ), pointer                :: basic            , basicChild     , &
         &                                                                          basicSibling
    class           (nodeComponentMergingStatistics    ), pointer                :: mergingStatistics
    double precision                                                             :: branchBegin      , branchEnd      , &
         &                                                                          branchMassInitial, branchMassFinal, &
         &                                                                          massRoot

    ! Iterate over trees.
    treeCurrent => tree    
    do while (associated(treeCurrent))
       ! Get root node of the tree.       
       node => treeCurrent%baseNode
       ! Walk the tree, searching for root nodes.
       do while (associated(node))
          ! Get the child node, and process if child exists.
          if (associated(node%firstChild)) then
             nodeChild => node%firstChild
          else
             nodeChild => node
          end if
          ! Check if child should be included.
          mergingStatistics => nodeChild%mergingStatistics()
          if     (                                                         &
               &   .not.self             %alwaysIsolatedHalosOnly          &
               &  .or.                                                     &
               &        mergingStatistics%nodeHierarchyLevelMaximum() == 0 &
               & ) then             
             ! Get the basic components.
             basic      => node     %basic()
             basicChild => nodeChild%basic()
             ! Determine range of times spanned by this branch.
             branchBegin=basicChild%time()
             branchEnd  =basic     %time()
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
                self%treeCount                =self       %treeCount   +1
                self%mass     (self%treeCount)=            massRoot
                self%weight   (self%treeCount)=treeCurrent%volumeWeight
                ! Flush the buffers to file if necessary.
                if (self%treeCount == outputRootMassesBufferSize) call self%finalize()
             end if
          end if
          ! Move to the next node.
          node => node%walkTree()
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine outputRootMassesOperate

  subroutine outputRootMassesFinalize(self)
    !% Outputs conditional mass function.
    use ISO_Varying_String
    use IO_HDF5
    use Galacticus_HDF5
    implicit none
    class(mergerTreeOperatorOutputRootMasses), intent(inout) :: self
    type (hdf5Object                        ), target        :: outputFile

    ! If the buffers are empty, we have nothing to do.
    if (self%treeCount == 0) return
    ! Open the output file.
    !$omp critical (HDF5_Access)
    call outputFile%openFile(char(self%fileName),overWrite=.false.)
    ! Write the data.
    call outputFile%writeDataset(self%mass  (1:self%treeCount),"treeRootMass","Tree root node masses.",appendTo=.true.)
    call outputFile%writeDataset(self%weight(1:self%treeCount),"treeWeight"  ,"Tree weights."         ,appendTo=.true.)
    ! Close the output file.
    call outputFile%close()
    !$omp end critical (HDF5_Access)
    ! Reset the buffer counter.
    self%treeCount=0
    return
  end subroutine outputRootMassesFinalize
