!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  Implements a node operator class that causes halo spin be interpolated linearly between child and parent nodes.
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
   contains
     final     ::                          darkMatterProfileScaleInterpolateConstructorDestructor
     procedure :: nodeTreeInitialize    => darkMatterProfileScaleInterpolateNodeTreeInitialize
     procedure :: nodeInitialize        => darkMatterProfileScaleInterpolateNodeInitialize
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
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDarkMatterProfileScaleInterpolate)                        :: self
    class(darkMatterProfileScaleRadiusClass            ), intent(in   ), target :: darkMatterProfileScaleRadius_
    !![
    <constructorAssign variables="*darkMatterProfileScaleRadius_"/>
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

  subroutine darkMatterProfileScaleInterpolateNodeTreeInitialize(self,node)
    !!{
    Initialize dark matter profile scale radii.
    !!}
    use :: Galacticus_Nodes   , only : nodeComponentDarkMatterProfile
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    implicit none
    class(nodeOperatorDarkMatterProfileScaleInterpolate), intent(inout)          :: self
    type (treeNode                                     ), intent(inout), target  :: node
    type (treeNode                                     )               , pointer :: nodeWork
    class(nodeComponentDarkMatterProfile               )               , pointer :: darkMatterProfile, darkMatterProfileWork
    type (mergerTreeWalkerAllNodes                     )                         :: treeWalker

    ! Initialize the scale radius, if necessary.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    if (darkMatterProfile%scale() < 0.0d0) then
       ! Perform our own depth-first tree walk to set scales in all nodes of the tree. This is necessary as we require access
       ! to the parent scale to set scale growth rates, but must initialize scales in a strictly depth-first manner as some
       ! algorithms rely on knowing the progenitor structure of the tree to compute scale radii.
       treeWalker=mergerTreeWalkerAllNodes(node%hostTree,spanForest=.true.)
       do while (treeWalker%next(nodeWork))
          ! Get the scale radius - this will initialize the radius if necessary.
          darkMatterProfileWork => nodeWork%darkMatterProfile(autoCreate=.true.)
          call darkMatterProfileWork%scaleSet(self%darkMatterProfileScaleRadius_%radius(nodeWork))
      end do
    end if
    return
  end subroutine darkMatterProfileScaleInterpolateNodeTreeInitialize
  
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
          call darkMatterProfile%scaleGrowthRateSet(                                         &
               &                                    (                                        &
               &                                     +darkMatterProfileParent%scale       () &
               &                                     -darkMatterProfile      %scale       () &
               &                                    )                                        &
               &                                    /                         timeInterval   &
               &                                   )
       else
          ! Time interval is non-positive - assume zero growth rate.
          call darkMatterProfile%scaleGrowthRateSet(0.0d0)
       end if
    else
       ! Node is a non-primary progenitor - assume zero growth rate.
       call    darkMatterProfile%scaleGrowthRateSet(0.0d0)
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
    darkMatterProfile => node             %darkMatterProfile()
    rateScaleRadius   =  darkMatterProfile%scaleGrowthRate  ()
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
