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
  Implements a node operator class that causes halo axis ratios to be interpolated linearly between child and parent nodes.
  !!}

  !![
  <nodeOperator name="nodeOperatorHaloAxisRatiosInterpolate">
   <description>
    A node operator class that causes halo axis ratios be interpolated linearly between child and parent nodes. For primary
    progenitor nodes then $\dot{\mathbf{a}} = (\mathbf{a}_{i+1}-\mathbf{a}_i)/(t_{i+1}-t_i)$, where $\mathbf{a}_i$ is the axis
    ratio tuple of the node in the initialized tree, $\mathbf{a}_{i+1}$ is the axis ratio tuple of its parent node, and $t_i$ and
    $t_{i+1}$ are the corresponding times.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloAxisRatiosInterpolate
     !!{
     A node operator class that causes halo axis ratios to be interpolated linearly between child and parent nodes.
     !!}
     private
     integer :: axisRatioRateID
   contains
     procedure :: nodeInitialize                      => haloAxisRatiosInterpolateNodeInitialize
     procedure :: differentialEvolutionAnalytics      => haloAxisRatiosInterpolateDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionSolveAnalytics => haloAxisRatiosInterpolateDifferentialEvolutionSolveAnalytics
     procedure :: nodePromote                         => haloAxisRatiosInterpolateNodePromote
  end type nodeOperatorHaloAxisRatiosInterpolate
  
  interface nodeOperatorHaloAxisRatiosInterpolate
     !!{
     Constructors for the \refClass{nodeOperatorHaloAxisRatiosInterpolate} node operator class.
     !!}
     module procedure haloAxisRatiosInterpolateConstructorParameters
     module procedure haloAxisRatiosInterpolateConstructorInternal
  end interface nodeOperatorHaloAxisRatiosInterpolate
  
contains
  
  function haloAxisRatiosInterpolateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorHaloAxisRatiosInterpolate} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorHaloAxisRatiosInterpolate)                :: self
    type(inputParameters                      ), intent(inout) :: parameters
     
    self=nodeOperatorHaloAxisRatiosInterpolate()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloAxisRatiosInterpolateConstructorParameters

  function haloAxisRatiosInterpolateConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorHaloAxisRatiosInterpolate} node operator class.
    !!}
    implicit none
    type(nodeOperatorHaloAxisRatiosInterpolate) :: self
     
    !![
    <addMetaProperty component="darkMatterProfile" name="axisRatioGrowth" id="self%axisRatioRateID" rank="1" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function haloAxisRatiosInterpolateConstructorInternal

  subroutine haloAxisRatiosInterpolateNodeInitialize(self,node)
    !!{
    Compute the rate of growth of halo axis ratios assuming a constant growth rate.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorHaloAxisRatiosInterpolate), intent(inout), target  :: self
    type            (treeNode                             ), intent(inout), target  :: node
    class           (nodeComponentBasic                   )               , pointer :: basic            , basicParent
    class           (nodeComponentDarkMatterProfile       )               , pointer :: darkMatterProfile, darkMatterProfileParent
    double precision                                                                :: timeInterval
    
    darkMatterProfile => node%darkMatterProfile()
    if (node%isPrimaryProgenitor()) then
       ! For primary progenitors compute and store the axis ratio growth rate.
       basic                   =>  node              %basic            ()
       basicParent             =>  node       %parent%basic            ()
       darkMatterProfileParent =>  node       %parent%darkMatterProfile()
       timeInterval            =  +basicParent       %time             () &
            &                     -basic             %time             ()
       if (timeInterval > 0.0d0) then
          call darkMatterProfile%floatRank1MetaPropertySet(                                              &
               &                                              self                   %axisRatioRateID  , &
               &                                           +(                                            &
               &                                             +darkMatterProfileParent%axisRatios     ()  &
               &                                             -darkMatterProfile      %axisRatios     ()  &
               &                                            )                                            &
               &                                           /timeInterval                                 &
               &                                          )
       else
          call darkMatterProfile%floatRank1MetaPropertySet(                                              &
               &                                              self                   %axisRatioRateID  , &
               &                                           [0.0d0,0.0d0,0.0d0]                           &
               &                                          )
       end if
    else
       ! For non-primary progenitors, assume no growth.
       call    darkMatterProfile%floatRank1MetaPropertySet(                                              &
            &                                                  self                   %axisRatioRateID , &
            &                                              [0.0d0,0.0d0,0.0d0]                           &
            &                                             )
    end if
    return
  end subroutine haloAxisRatiosInterpolateNodeInitialize
  
  subroutine haloAxisRatiosInterpolateDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorHaloAxisRatiosInterpolate), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile       ), pointer       :: darkMatterProfile
    !$GLC attributes unused :: self

    darkMatterProfile => node%darkMatterProfile()
    call darkMatterProfile%axisRatiosAnalytic()
    return
  end subroutine haloAxisRatiosInterpolateDifferentialEvolutionAnalytics

  subroutine haloAxisRatiosInterpolateDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Evolve halo angular momentum at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorHaloAxisRatiosInterpolate), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: time
    class           (nodeComponentBasic                   ), pointer       :: basicParent
    class           (nodeComponentDarkMatterProfile       ), pointer       :: darkMatterProfile, darkMatterProfileParent
    !$GLC attributes unused :: self

    if (.not.node%isPrimaryProgenitor()) return
    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    basicParent             => node%parent%basic            ()   
    call darkMatterProfile%axisRatiosSet(                                                                           &
         &                               +  darkMatterProfileParent%axisRatios               (                    ) &
         &                               +(                                                                         &
         &                                 +                        time                                            &
         &                                 -basicParent            %time                     (                    ) &
         &                                )                                                                         &
         &                               *  darkMatterProfile      %floatRank1MetaPropertyGet(self%axisRatioRateID) &
         &                              )
    return
  end subroutine haloAxisRatiosInterpolateDifferentialEvolutionSolveAnalytics

  subroutine haloAxisRatiosInterpolateNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the axis ratios of {\normalfont \ttfamily node}
    to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class(nodeOperatorHaloAxisRatiosInterpolate), intent(inout)  :: self
    type (treeNode                             ), intent(inout)  :: node
    type (treeNode                             ), pointer        :: nodeParent
    class(nodeComponentDarkMatterProfile       ), pointer        :: darkMatterProfileParent, darkMatterProfile

    nodeParent              => node      %parent
    darkMatterProfile       => node      %darkMatterProfile()
    darkMatterProfileParent => nodeParent%darkMatterProfile()
    call darkMatterProfile%axisRatiosSet            (                     darkMatterProfileParent%axisRatios               (                    ))
    call darkMatterProfile%floatRank1MetaPropertySet(self%axisRatioRateID,darkMatterProfileParent%floatRank1MetaPropertyGet(self%axisRatioRateID))
    return
  end subroutine haloAxisRatiosInterpolateNodePromote
