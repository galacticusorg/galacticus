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
  Implements a node operator class that causes dark matter profile shape parameter to be interpolated linearly between child and parent nodes.
  !!}

  use :: Dark_Matter_Profiles_Shape, only : darkMatterProfileShapeClass
  
  !![
  <nodeOperator name="nodeOperatorDarkMatterProfileShapeInterpolate">
   <description>
    A node operator class that causes dark matter profile shape parameter to be interpolated linearly between child and parent
    nodes. For primary progenitor nodes $\dot{\alpha} = (\alpha_{i+1}-\alpha_{i})/(t_{i+1}-t_i)$, where $\alpha_{i}$ is the shape
    parameter of the dark matter profile of the node in the initialized tree, $\alpha_{i+1}$ is the spin of its parent node, and
    $t_i$ and $t_{i+1}$ are the corresponding times. For non-primary progenitors the rate of change is set to zero,
    i.e. $\dot{\alpha}=0$.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfileShapeInterpolate
     !!{
     A node operator class that causes dark matter profile shape parameter to be interpolated linearly between child and parent nodes.
     !!}
     private
     class(darkMatterProfileShapeClass), pointer :: darkMatterProfileShape_ => null()
     integer                                     :: shapeGrowthRateID
   contains
     final     ::                                        dmpShapeInterpolateConstructorDestructor
     procedure :: nodeInitialize                      => dmpShapeInterpolateNodeInitialize
     procedure :: nodePromote                         => dmpShapeInterpolateNodePromote
     procedure :: differentialEvolutionAnalytics      => dmpShapeInterpolateDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionSolveAnalytics => dmpShapeInterpolateDifferentialEvolutionSolveAnalytics
  end type nodeOperatorDarkMatterProfileShapeInterpolate
  
  interface nodeOperatorDarkMatterProfileShapeInterpolate
     !!{
     Constructors for the \refClass{nodeOperatorDarkMatterProfileShapeInterpolate} node operator class.
     !!}
     module procedure dmpShapeInterpolateConstructorParameters
     module procedure dmpShapeInterpolateConstructorInternal
  end interface nodeOperatorDarkMatterProfileShapeInterpolate
  
contains
  
  function dmpShapeInterpolateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileShapeInterpolate} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDarkMatterProfileShapeInterpolate)                :: self
    type (inputParameters                              ), intent(inout) :: parameters
    class(darkMatterProfileShapeClass                  ), pointer       :: darkMatterProfileShape_

    !![
    <objectBuilder class="darkMatterProfileShape" name="darkMatterProfileShape_" source="parameters"/>
    !!]
    self=nodeOperatorDarkMatterProfileShapeInterpolate(darkMatterProfileShape_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileShape_"/>
    !!]
    return
  end function dmpShapeInterpolateConstructorParameters

  function dmpShapeInterpolateConstructorInternal(darkMatterProfileShape_) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileShapeInterpolate} node operator class which takes a parameter set as input.
    !!}
    implicit none
    type (nodeOperatorDarkMatterProfileShapeInterpolate)                        :: self
    class(darkMatterProfileShapeClass                  ), intent(in   ), target :: darkMatterProfileShape_
    !![
    <constructorAssign variables="*darkMatterProfileShape_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="shapeGrowthRate" id="self%shapeGrowthRateID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function dmpShapeInterpolateConstructorInternal

  subroutine dmpShapeInterpolateConstructorDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorDarkMatterProfileShapeInterpolate} dark matter halo profile shape parameter class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfileShapeInterpolate), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileShape_"/>
    !!]
    return
  end subroutine dmpShapeInterpolateConstructorDestructor

  subroutine dmpShapeInterpolateNodeInitialize(self,node)
    !!{
    Compute the rate of growth of dark matter profile shape parameter assuming a constant growth rate.
    !!}
    use :: Display         , only : displayBlue       , displayGreen                  , displayYellow, displayBold, &
         &                          displayReset
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileShapeInterpolate), intent(inout), target  :: self
    type            (treeNode                                     ), intent(inout), target  :: node
    class           (nodeComponentBasic                           )               , pointer :: basic            , basicParent
    class           (nodeComponentDarkMatterProfile               )               , pointer :: darkMatterProfile, darkMatterProfileParent
    double precision                                                                        :: timeInterval

    ! Set the growth rate for the shape parameter.
    call self%nodeTreeInitialize(node)
    darkMatterProfile => node%darkMatterProfile()
    select type (darkMatterProfile)
    type is (nodeComponentDarkMatterProfile)
       call Error_Report(                                                                                                                                                                                            &
            &            displayBold()//'darkMatterProfile'//displayReset()//' component must be created prior to initialization shape parameter interpolation'                                         //char(10)// &
            &            '   For example, by using the following nodeOperator'                                                                                                                          //char(10)// &
            &            '    <'//displayBlue()//'nodeOperator'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"darkMatterProfileShapeSet"'//displayReset()//'/>'          // &
            &            {introspection:location}                                                                                                                                                                    &
            &           )
    class default
       if (node%isPrimaryProgenitor()) then
          ! Node is the primary progenitor, so compute the shape parameter growth rate.
          basic        =>  node              %basic()
          basicParent  =>  node       %parent%basic()
          timeInterval =  +basicParent       %time () &
               &          -basic             %time ()
          if (timeInterval > 0.0d0) then
             darkMatterProfileParent => node%parent%darkMatterProfile()
             call darkMatterProfile%floatRank0MetaPropertySet(                                                &
                  &                                              self                   %shapeGrowthRateID  , &
                  &                                           +(                                              &
                  &                                             +darkMatterProfileParent%shape            ()  &
                  &                                             -darkMatterProfile      %shape            ()  &
                  &                                            )                                              &
                  &                                           /                          timeInterval         &
                  &                                          )
          else
             ! Time interval is non-positive - assume zero growth rate.
             call darkMatterProfile%floatRank0MetaPropertySet(self%shapeGrowthRateID,0.0d0)
          end if
       else
          ! Node is a non-primary progenitor - assume zero growth rate.
          call    darkMatterProfile%floatRank0MetaPropertySet(self%shapeGrowthRateID,0.0d0)
       end if
    end select
    return
  end subroutine dmpShapeInterpolateNodeInitialize
  
  subroutine dmpShapeInterpolateDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileShapeInterpolate), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile               ), pointer       :: darkMatterProfile

    darkMatterProfile => node%darkMatterProfile()
    call darkMatterProfile%shapeAnalytic()
    return
  end subroutine dmpShapeInterpolateDifferentialEvolutionAnalytics

  subroutine dmpShapeInterpolateDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Evolve dark matter profile shape parameter at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileShapeInterpolate), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: time
    class           (nodeComponentBasic                           ), pointer       :: basicParent
    class           (nodeComponentDarkMatterProfile               ), pointer       :: darkMatterProfile, darkMatterProfileParent

    if (.not.node%isPrimaryProgenitor()) return
    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    basicParent             => node%parent%basic            ()
    call darkMatterProfile%shapeSet(                                                                           &
         &                          +darkMatterProfileParent%shape                    (                      ) &
         &                          +(                                                                         &
         &                            +                      time                                              &
         &                            -basicParent          %time                     (                      ) &
         &                           )                                                                         &
         &                          *darkMatterProfile      %floatRank0MetaPropertyGet(self%shapeGrowthRateID) &
         &                         )
    return
  end subroutine dmpShapeInterpolateDifferentialEvolutionSolveAnalytics

  subroutine dmpShapeInterpolateNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the shape parameter
    growth rate of {\normalfont \ttfamily node} to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileShapeInterpolate), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile               ), pointer       :: darkMatterProfile, darkMatterProfileParent
    
    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    call darkMatterProfile%floatRank0MetaPropertySet(self%shapeGrowthRateID,darkMatterProfileParent%floatRank0MetaPropertyGet(self%shapeGrowthRateID))
    call darkMatterProfile%                 shapeSet(                       darkMatterProfileParent%shape                    (                      ))
    return
  end subroutine dmpShapeInterpolateNodePromote
