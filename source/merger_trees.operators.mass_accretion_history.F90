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
  Implements a merger tree operator which outputs mass accretion
  histories.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: IO_HDF5                , only : hdf5Object

  !![
  <mergerTreeOperator name="mergerTreeOperatorMassAccretionHistory">
   <description>
    A merger tree operator class which outputs mass accretion histories (i.e. the mass of the \gls{node} on the primary branch
    as a function of time). Histories are written into the \glc\ output file in a group with name given by {\normalfont
    \ttfamily [outputGroupName]}. Within that group, each merger tree has its own group named {\normalfont \ttfamily
    mergerTree\textless\ N\textgreater} where {\normalfont \ttfamily \textless\ N\textgreater} is the tree index. Within each
    such merger tree group datasets giving the node index (``{\normalfont \ttfamily nodeIndex}''), time (``{\normalfont
    \ttfamily nodeTime}''), basic mass (``{\normalfont \ttfamily nodeMass}''), expansion factor (``{\normalfont \ttfamily
    nodeExpansionFactor}'') are written. Optionally, datasets giving the spin parameter (``{\normalfont \ttfamily nodeSpin}'')
    and its vector components (``{\normalfont \ttfamily nodeSpinVector}'') are included if {\normalfont \ttfamily
    [includeSpin]} and {\normalfont \ttfamily [includeSpinVector]} respectively are set to {\normalfont \ttfamily true}.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorMassAccretionHistory
     !!{
     A merger tree operator class which outputs mass accretion histories.
     !!}
     private
     type   (hdf5Object               )          :: outputGroup
     type   (varying_string           )          :: outputGroupName
     class  (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_  => null()
     class  (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_ => null()
     logical                                     :: includeSpin                   , includeSpinVector
   contains
     final     ::                        massAccretionHistoryDestructor
     procedure :: operatePreEvolution => massAccretionHistoryOperatePreEvolution
     procedure :: finalize            => massAccretionHistoryFinalize
  end type mergerTreeOperatorMassAccretionHistory

  interface mergerTreeOperatorMassAccretionHistory
     !!{
     Constructors for the mass accretion history merger tree operator class.
     !!}
     module procedure massAccretionHistoryConstructorParameters
     module procedure massAccretionHistoryConstructorInternal
  end interface mergerTreeOperatorMassAccretionHistory

contains

  function massAccretionHistoryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the mass accretion history merger tree operator class which takes a
    parameter set as input.
    !!}
    implicit none
    type   (mergerTreeOperatorMassAccretionHistory)                :: self
    type   (inputParameters                       ), intent(inout) :: parameters
    type   (varying_string                        )                :: outputGroupName
    class  (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class  (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    logical                                                        :: includeSpin          , includeSpinVector

    !![
    <inputParameter>
      <name>outputGroupName</name>
      <source>parameters</source>
      <defaultValue>var_str('massAccretionHistories')</defaultValue>
      <description>The name of the \gls{hdf5} group to output mass accretion histories to.</description>
    </inputParameter>
    <inputParameter>
      <name>includeSpin</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, include the spin of the halo in the output.</description>
    </inputParameter>
    <inputParameter>
      <name>includeSpinVector</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, include the spin vector of the halo in the output.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=mergerTreeOperatorMassAccretionHistory(char(outputGroupName),includeSpin,includeSpinVector,cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function massAccretionHistoryConstructorParameters

  function massAccretionHistoryConstructorInternal(outputGroupName,includeSpin,includeSpinVector,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the mass accretion history merger tree operator class.
    !!}
    use :: Error           , only : Component_List      , Error_Report
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none
    type     (mergerTreeOperatorMassAccretionHistory)                        :: self
    character(len=*                                 ), intent(in   )         :: outputGroupName
    logical                                          , intent(in   )         :: includeSpin         , includeSpinVector
    class    (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class    (darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="outputGroupName, includeSpin, includeSpinVector, *cosmologyFunctions_, *darkMatterHaloScale_"/>
    !!]

    if (self%includeSpin      .and..not.defaultSpinComponent%angularMomentumIsGettable      ())                 &
         & call Error_Report                                                                                    &
         &  (                                                                                                   &
         &   'the angularMomentum property of the spin component must be gettable.'                          // &
         &   Component_List(                                                                                    &
         &                  'spin'                                                                            , &
         &                   defaultSpinComponent%angularMomentumAttributeMatch      (requireGettable=.true.)   &
         &                 )                                                                                 // &
         &   {introspection:location}                                                                           &
         &  )
    if (self%includeSpinVector.and..not.defaultSpinComponent%angularMomentumVectorIsGettable())                 &
         & call Error_Report                                                                                    &
         &  (                                                                                                   &
         &   'the angularMomentumVector property of the spin component must be gettable.'                    // &
         &   Component_List(                                                                                    &
         &                  'spin'                                                                            , &
         &                   defaultSpinComponent%angularMomentumVectorAttributeMatch(requireGettable=.true.)   &
         &                 )                                                                                 // &
         &   {introspection:location}                                                                           &
         &  )
    return
  end function massAccretionHistoryConstructorInternal

  subroutine massAccretionHistoryDestructor(self)
    !!{
    Destructor for the mass accretion history merger tree operator function class.
    !!}
    implicit none
    type(mergerTreeOperatorMassAccretionHistory), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine massAccretionHistoryDestructor

  subroutine massAccretionHistoryOperatePreEvolution(self,tree)
    !!{
    Output the mass accretion history for a merger tree.
    !!}
    use            :: Dark_Matter_Halo_Spins          , only : Dark_Matter_Halo_Angular_Momentum_Scale
    use            :: Display                         , only : displayGreen                           , displayReset
    use            :: Error                           , only : Error_Report
    use            :: Output_HDF5                     , only : outputFile
    use            :: Galacticus_Nodes                , only : mergerTree                             , nodeComponentBasic, nodeComponentSpin, treeNode
    use            :: HDF5_Access                     , only : hdf5Access
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: ISO_Varying_String              , only : varying_string
    use            :: Numerical_Constants_Astronomical, only : gigaYear                               , massSolar
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class           (mergerTreeOperatorMassAccretionHistory), intent(inout), target         :: self
    type            (mergerTree                            ), intent(inout), target         :: tree
    type            (treeNode                              )               , pointer        :: node
    integer         (kind=kind_int8                        ), allocatable  , dimension(:  ) :: nodeIndex
    double precision                                        , allocatable  , dimension(:  ) :: nodeMass             , nodeTime, &
         &                                                                                     nodeExpansionFactor  , nodeSpin
    double precision                                        , allocatable  , dimension(:,:) :: nodeSpinVector
    class           (nodeComponentBasic                    )               , pointer        :: basic
    class           (nodeComponentSpin                     )               , pointer        :: spin
    type            (mergerTree                            )               , pointer        :: treeCurrent
    integer         (c_size_t                              )                                :: accretionHistoryCount
    type            (varying_string                        )                                :: groupName
    type            (hdf5Object                            )                                :: accretionDataset     , treeGroup

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Count up number of entries expected for accretion history.
       accretionHistoryCount =  0
       node                  => treeCurrent%nodeBase
       do while (associated(node))
          accretionHistoryCount =  accretionHistoryCount           +1
          node                  => node                 %firstChild
       end do
       ! Allocate storage space.
       allocate(nodeIndex          (int(accretionHistoryCount)  ))
       allocate(nodeTime           (int(accretionHistoryCount)  ))
       allocate(nodeExpansionFactor(int(accretionHistoryCount)  ))
       allocate(nodeMass           (int(accretionHistoryCount)  ))
       if (self%includeSpin      ) allocate(nodeSpin           (int(accretionHistoryCount)  ))
       if (self%includeSpinVector) allocate(nodeSpinVector     (int(accretionHistoryCount),3))
       ! Extract accretion history.
       accretionHistoryCount =  0
       node                  => treeCurrent%nodeBase
       do while (associated(node))
          accretionHistoryCount                                               =  accretionHistoryCount+1
          basic                                                               =>                                          node%basic      (                 )
          if     (                        &
               &   self%includeSpin       &
               &  .or.                    &
               &   self%includeSpinVector &
               & )                        &
               & spin                                                         =>                                           node%spin                  (autoCreate=.true.)
          nodeIndex                                 (accretionHistoryCount  ) =                                            node %index                (                 )
          nodeTime                                  (accretionHistoryCount  ) =                                            basic%time                 (                 )
          nodeMass                                  (accretionHistoryCount  ) =                                            basic%mass                 (                 )
          nodeExpansionFactor                       (accretionHistoryCount  ) =   self%cosmologyFunctions_%expansionFactor(basic%time                 (                 ))
          if (self%includeSpin      ) nodeSpin      (accretionHistoryCount  ) =                                            spin %angularMomentum      (                 )  &
               &                                                                 /Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_)
          if (self%includeSpinVector) nodeSpinVector(accretionHistoryCount,:) =                                            spin %angularMomentumVector(                 )  &
               &                                                                 /Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_)
          node                                                                =>                                           node %firstChild
       end do
       ! Create the output group if necessary.
       !$ call hdf5Access%set()
       if (.not.self%outputGroup%isOpen()) self%outputGroup=outputFile%openGroup(char(self%outputGroupName),'Mass accretion histories of main branches in merger trees.')
       !$ call hdf5Access%unset()
       ! Output to HDF5 file.
       groupName='mergerTree'
       groupName=groupName//treeCurrent%index
       !$ call hdf5Access%set()
       if (self%outputGroup%hasGroup(char(groupName)))                                                                                                                        &
            & call Error_Report(                                                                                                                                              &
            &                   'duplicate tree index detected - mass accretion history can not be output'                                                       //char(10)// &
            &                   {introspection:location}                                                                                                                   // &
            &                   displayGreen()//'  HELP:'//displayReset()                                                                                                  // &
            &                   ' This can happen if reading merger trees which contain multiple root nodes from file. To avoid this problem, force tree indices'          // &
            &                   ' to be reset to the index of the root node by adding the following to your input parameter file:'                               //char(10)// &
            &                   '  <mergerTreeConstructor value="read">'                                                                                         //char(10)// &
            &                   '     <treeIndexToRootNodeIndex value="true"/>'                                                                                  //char(10)// &
            &                   '  </mergerTreeConstructor>'                                                                                                                  &
            &                  )
       treeGroup=self%outputGroup%openGroup(char(groupName),'Mass accretion history for main branch of merger tree.')
       call                             treeGroup       %writeDataset  (nodeIndex          ,'nodeIndex'          ,'Index of the node.'                                            )
       call                             treeGroup       %writeDataset  (nodeTime           ,'nodeTime'           ,'Time at node [Gyr].'          ,datasetReturned=accretionDataset)
       call                             accretionDataset%writeAttribute(gigaYear                                 ,'unitsInSI'                                                     )
       call                             accretionDataset%close         (                                                                                                          )
       call                             treeGroup       %writeDataset  (nodeMass           ,'nodeMass'           ,'Mass of the node [M⊙].'       ,datasetReturned=accretionDataset)
       call                             accretionDataset%writeAttribute(massSolar                                ,"unitsInSI"                                                     )
       call                             accretionDataset%close         (                                                                                                          )
       call                             treeGroup       %writeDataset  (nodeExpansionFactor,'nodeExpansionFactor','Expansion factor of the node.'                                 )
       if (self%includeSpin      ) call treeGroup       %writeDataset  (nodeSpin           ,'nodeSpin','Spin parameter of the node.'                                              )
       if (self%includeSpinVector) call treeGroup       %writeDataset  (nodeSpinVector     ,'nodeSpinVector','Spin parameter vector of the node.'                                 )
       call treeGroup       %close         (                                                                                                                                      )
       !$ call hdf5Access%unset()
       ! Deallocate storage space.
       deallocate(nodeIndex          )
       deallocate(nodeTime           )
       deallocate(nodeMass           )
       deallocate(nodeExpansionFactor)
       if (self%includeSpin      ) deallocate(nodeSpin           )
       if (self%includeSpinVector) deallocate(nodeSpinVector     )
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine massAccretionHistoryOperatePreEvolution

  subroutine massAccretionHistoryFinalize(self)
    !!{
    Close the mass accretion history group before closing the HDF5 file.
    !!}
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class(mergerTreeOperatorMassAccretionHistory), intent(inout) :: self

    !$ call hdf5Access%set()
    if (self%outputGroup%isOpen()) call self%outputGroup%close()
    !$ call hdf5Access%unset()
    return
  end subroutine massAccretionHistoryFinalize
