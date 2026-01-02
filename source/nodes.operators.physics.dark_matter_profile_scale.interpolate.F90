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
  Implements a node operator class that causes dark matter profile scale radius to be interpolated linearly between child and parent nodes.
  !!}

  !![
  <nodeOperator name="nodeOperatorDarkMatterProfileScaleInterpolate">
   <description>
    A node operator class that causes dark matter profile scale radius to be interpolated linearly between child and parent
    nodes. For primary progenitor nodes $\dot{r}_\mathrm{s} = (r_{\mathrm{s},i+1}-r_{\mathrm{s},i})/(t_{i+1}-t_i)$, where
    $r_{\mathrm{s},i}$ is the scale radius of the dark matter profile of the node in the initialized tree, $r_{\mathrm{s},i+1}$
    is the spin of its parent node, and $t_i$ and $t_{i+1}$ are the corresponding times. For non-primary progenitors the rate
    of change is set to zero, i.e. $\dot{r}_\mathrm{s}=0$.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfileScaleInterpolate
     !!{
     A node operator class that causes dark matter profile scale radius to be interpolated linearly between child and parent nodes.
     !!}
     private
     integer :: scaleGrowthRateID
   contains
     procedure :: nodeInitialize                      => dmpScaleInterpolateNodeInitialize
     procedure :: nodePromote                         => dmpScaleInterpolateNodePromote
     procedure :: differentialEvolutionAnalytics      => dmpScaleInterpolateDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionSolveAnalytics => dmpScaleInterpolateDifferentialEvolutionSolveAnalytics
  end type nodeOperatorDarkMatterProfileScaleInterpolate
  
  interface nodeOperatorDarkMatterProfileScaleInterpolate
     !!{
     Constructors for the \refClass{nodeOperatorDarkMatterProfileScaleInterpolate} node operator class.
     !!}
     module procedure dmpScaleInterpolateConstructorParameters
     module procedure dmpScaleInterpolateConstructorInternal
  end interface nodeOperatorDarkMatterProfileScaleInterpolate
  
contains
  
  function dmpScaleInterpolateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileScaleInterpolate} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorDarkMatterProfileScaleInterpolate)                :: self
    type(inputParameters                              ), intent(inout) :: parameters

    self=nodeOperatorDarkMatterProfileScaleInterpolate()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function dmpScaleInterpolateConstructorParameters

  function dmpScaleInterpolateConstructorInternal() result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileScaleInterpolate} node operator class which takes a parameter set as input.
    !!}
    implicit none
    type (nodeOperatorDarkMatterProfileScaleInterpolate) :: self

    !![
    <addMetaProperty component="darkMatterProfile" name="scaleGrowthRate" id="self%scaleGrowthRateID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function dmpScaleInterpolateConstructorInternal

  subroutine dmpScaleInterpolateNodeInitialize(self,node)
    !!{
    Compute the rate of growth of dark matter profile scale radius assuming a constant growth rate.
    !!}
    use :: Display         , only : displayBlue       , displayGreen                  , displayYellow, displayBold, &
         &                          displayReset
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileScaleInterpolate), intent(inout), target  :: self
    type            (treeNode                                     ), intent(inout), target  :: node
    class           (nodeComponentBasic                           )               , pointer :: basic            , basicParent
    class           (nodeComponentDarkMatterProfile               )               , pointer :: darkMatterProfile, darkMatterProfileParent
    double precision                                                                        :: timeInterval

    ! Set the growth rate for the scale radius.
    call self%nodeTreeInitialize(node)
    darkMatterProfile => node%darkMatterProfile()
    select type (darkMatterProfile)
    type is (nodeComponentDarkMatterProfile)
       call Error_Report(                                                                                                                                                                                            &
            &            displayBold()//'darkMatterProfile'//displayReset()//' component must be created prior to initialization scale radius interpolation'                                            //char(10)// &
            &            '   For example, by using the following nodeOperator'                                                                                                                          //char(10)// &
            &            '    <'//displayBlue()//'nodeOperator'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"darkMatterProfileScaleSet"'//displayReset()//'/>'          // &
            &            {introspection:location}                                                                                                                                                                    &
            &           )
    class default
       if (node%isPrimaryProgenitor()) then
          ! Node is the primary progenitor, so compute the scale radius growth rate.
          basic        =>  node              %basic()
          basicParent  =>  node       %parent%basic()
          timeInterval =  +basicParent       %time () &
               &          -basic             %time ()
          if (timeInterval > 0.0d0) then
             darkMatterProfileParent => node%parent%darkMatterProfile()
             call darkMatterProfile%floatRank0MetaPropertySet(                                                &
                  &                                              self                   %scaleGrowthRateID  , &
                  &                                           +(                                              &
                  &                                             +darkMatterProfileParent%scale            ()  &
                  &                                             -darkMatterProfile      %scale            ()  &
                  &                                            )                                              &
                  &                                           /                          timeInterval         &
                  &                                          )
          else
             ! Time interval is non-positive - assume zero growth rate.
             call darkMatterProfile%floatRank0MetaPropertySet(self%scaleGrowthRateID,0.0d0)
          end if
       else
          ! Node is a non-primary progenitor - assume zero growth rate.
          call    darkMatterProfile%floatRank0MetaPropertySet(self%scaleGrowthRateID,0.0d0)
       end if
    end select
    return
  end subroutine dmpScaleInterpolateNodeInitialize
  
  subroutine dmpScaleInterpolateDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileScaleInterpolate), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile               ), pointer       :: darkMatterProfile

    darkMatterProfile => node%darkMatterProfile()
    call darkMatterProfile%scaleAnalytic()
    return
  end subroutine dmpScaleInterpolateDifferentialEvolutionAnalytics

  subroutine dmpScaleInterpolateDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Evolve dark matter profile scale radius at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileScaleInterpolate), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: time
    class           (nodeComponentBasic                           ), pointer       :: basicParent
    class           (nodeComponentDarkMatterProfile               ), pointer       :: darkMatterProfile, darkMatterProfileParent
    
    if (.not.node%isPrimaryProgenitor()) return
    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    basicParent             => node%parent%basic            ()
    call darkMatterProfile%scaleSet(                                                                           &
         &                          +darkMatterProfileParent%scale                    (                      ) &
         &                          +(                                                                         &
         &                            +                      time                                              &
         &                            -basicParent          %time                     (                      ) &
         &                           )                                                                         &
         &                          *darkMatterProfile      %floatRank0MetaPropertyGet(self%scaleGrowthRateID) &
         &                         )
    return
  end subroutine dmpScaleInterpolateDifferentialEvolutionSolveAnalytics

  subroutine dmpScaleInterpolateNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the scale radius
    growth rate of {\normalfont \ttfamily node} to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileScaleInterpolate), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile               ), pointer       :: darkMatterProfile, darkMatterProfileParent
    
    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    call darkMatterProfile%floatRank0MetaPropertySet(self%scaleGrowthRateID,darkMatterProfileParent%floatRank0MetaPropertyGet(self%scaleGrowthRateID))
    call darkMatterProfile%                 scaleSet(                       darkMatterProfileParent%scale                    (                      ))
    return
  end subroutine dmpScaleInterpolateNodePromote
