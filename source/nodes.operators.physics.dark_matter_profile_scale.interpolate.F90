!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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

  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadiusClass
  
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
     class(darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     integer                                           :: scaleGrowthRateID
   contains
     final     ::                          darkMatterProfileScaleInterpolateConstructorDestructor
     procedure :: nodeInitialize        => darkMatterProfileScaleInterpolateNodeInitialize
     procedure :: nodePromote           => darkMatterProfileScaleInterpolateNodePromote
     procedure :: differentialEvolution => darkMatterProfileScaleInterpolateDifferentialEvolution
  end type nodeOperatorDarkMatterProfileScaleInterpolate
  
  interface nodeOperatorDarkMatterProfileScaleInterpolate
     !!{
     Constructors for the {\normalfont \ttfamily darkMatterProfileScaleInterpolate} node operator class.
     !!}
     module procedure darkMatterProfileScaleInterpolateConstructorParameters
     module procedure darkMatterProfileScaleInterpolateConstructorInternal
  end interface nodeOperatorDarkMatterProfileScaleInterpolate
  
contains
  
  function darkMatterProfileScaleInterpolateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily darkMatterProfileScaleInterpolate} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDarkMatterProfileScaleInterpolate)                :: self
    type (inputParameters                              ), intent(inout) :: parameters
    class(darkMatterProfileScaleRadiusClass            ), pointer       :: darkMatterProfileScaleRadius_

    !![
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    self=nodeOperatorDarkMatterProfileScaleInterpolate(darkMatterProfileScaleRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function darkMatterProfileScaleInterpolateConstructorParameters

  function darkMatterProfileScaleInterpolateConstructorInternal(darkMatterProfileScaleRadius_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily darkMatterProfileScaleInterpolate} node operator class which takes a parameter set as input.
    !!}
    implicit none
    type (nodeOperatorDarkMatterProfileScaleInterpolate)                        :: self
    class(darkMatterProfileScaleRadiusClass            ), intent(in   ), target :: darkMatterProfileScaleRadius_
    !![
    <constructorAssign variables="*darkMatterProfileScaleRadius_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="scaleGrowthRate" id="self%scaleGrowthRateID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function darkMatterProfileScaleInterpolateConstructorInternal

  subroutine darkMatterProfileScaleInterpolateConstructorDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily darkMatterProfileScaleInterpolate} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfileScaleInterpolate), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    !!]
    return
  end subroutine darkMatterProfileScaleInterpolateConstructorDestructor

  subroutine darkMatterProfileScaleInterpolateNodeInitialize(self,node)
    !!{
    Compute the rate of growth of dark matter profile scale radius assuming a constant growth rate.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileScaleInterpolate), intent(inout)          :: self
    type            (treeNode                                     ), intent(inout), target  :: node
    class           (nodeComponentBasic                           )               , pointer :: basic            , basicParent
    class           (nodeComponentDarkMatterProfile               )               , pointer :: darkMatterProfile, darkMatterProfileParent
    double precision                                                                        :: timeInterval

    ! Set the growth rate for the scale radius.
    call self%nodeTreeInitialize(node)
    darkMatterProfile => node%darkMatterProfile()
    if (node%isPrimaryProgenitor()) then
       ! Node is the primary progenitor, so compute the scale radius growth rate.
       basic        =>  node              %basic()
       basicParent  =>  node       %parent%basic()
       timeInterval =  +basicParent       %time () &
            &          -basic             %time ()
       if (timeInterval > 0.0d0) then
          darkMatterProfileParent => node%parent%darkMatterProfile()
          call darkMatterProfile%metaPropertySet(                                                &
               &                                    self                   %scaleGrowthRateID  , &
               &                                 +(                                              &
               &                                   +darkMatterProfileParent%scale            ()  &
               &                                   -darkMatterProfile      %scale            ()  &
               &                                  )                                              &
               &                                 /                          timeInterval         &
               &                                )
       else
          ! Time interval is non-positive - assume zero growth rate.
          call darkMatterProfile%metaPropertySet(self%scaleGrowthRateID,0.0d0)
       end if
    else
       ! Node is a non-primary progenitor - assume zero growth rate.
       call    darkMatterProfile%metaPropertySet(self%scaleGrowthRateID,0.0d0)
    end if
    return
  end subroutine darkMatterProfileScaleInterpolateNodeInitialize
  
  subroutine darkMatterProfileScaleInterpolateDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Evolve dark matter profile scale radius at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, propertyTypeInactive
    implicit none
    class           (nodeOperatorDarkMatterProfileScaleInterpolate), intent(inout), target  :: self
    type            (treeNode                                     ), intent(inout)          :: node
    logical                                                        , intent(inout)          :: interrupt
    procedure       (interruptTask                                ), intent(inout), pointer :: functionInterrupt
    integer                                                        , intent(in   )          :: propertyType
    class           (nodeComponentDarkMatterProfile               )               , pointer :: darkMatterProfile, darkMatterProfileParent
    double precision                                                                        :: rateScaleRadius
    !$GLC attributes unused :: interrupt, functionInterrupt
    
    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Get the dark matter profile component.
    darkMatterProfile => node             %darkMatterProfile(                      )
    rateScaleRadius   =  darkMatterProfile%metaPropertyGet  (self%scaleGrowthRateID)
    if (node%isPrimaryProgenitor()) then
       ! If necessary, limit the growth rate so that we do not exceed the scale radius of the parent halo.
       darkMatterProfileParent => node%parent%darkmatterProfile()
       if     (                                                                                              &
            &   (rateScaleRadius > 0.0d0 .and. darkMatterProfile%scale() >= darkMatterProfileParent%scale()) &
            &  .or.                                                                                          &
            &   (rateScaleRadius < 0.0d0 .and. darkMatterProfile%scale() <= darkMatterProfileParent%scale()) &
            & ) rateScaleRadius=0.0d0
    end if
    call darkMatterProfile%scaleRate(rateScaleRadius)
    return
  end subroutine darkMatterProfileScaleInterpolateDifferentialEvolution

  subroutine darkMatterProfileScaleInterpolateNodePromote(self,node)
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
    call darkMatterProfile%metaPropertySet(self%scaleGrowthRateID,darkMatterProfileParent%metaPropertyGet(self%scaleGrowthRateID))
    return
  end subroutine darkMatterProfileScaleInterpolateNodePromote
