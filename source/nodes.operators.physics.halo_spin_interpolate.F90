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

  !![
  <nodeOperator name="nodeOperatorHaloSpinInterpolate">
   <description>
    A node operator class that causes halo spin be interpolated linearly between child and parent nodes. For primary progenitor
    nodes, if only the scalar spin, $\lambda$, is available then $\dot{\lambda} = (\lambda_{i+1}-\lambda_i)/(t_{i+1}-t_i)$,
    where $\lambda_i$ is the spin of the node in the initialized tree, $\lambda_{i+1}$ is the spin of its parent node, and
    $t_i$ and $t_{i+1}$ are the corresponding times. If vector spin is available the same interpolation is applied to each
    individual component of spin, with the rate of change of the scalar spin computed self-consistently. For non-primary
    progenitors both scalar and vector spin are assumed to be constant, i.e. $\dot{\lambda}=0$.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloSpinInterpolate
     !!{
     A node operator class that causes halo spin be interpolated linearly between child and parent nodes.
     !!}
     private
   contains
     procedure :: nodeInitialize        => haloSpinInterpolateNodeInitialize
     procedure :: differentialEvolution => haloSpinInterpolateDifferentialEvolution
  end type nodeOperatorHaloSpinInterpolate
  
  interface nodeOperatorHaloSpinInterpolate
     !!{
     Constructors for the {\normalfont \ttfamily haloSpinInterpolate} node operator class.
     !!}
     module procedure haloSpinInterpolateConstructorParameters
  end interface nodeOperatorHaloSpinInterpolate
  
contains
  
  function haloSpinInterpolateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily haloSpinInterpolate} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorHaloSpinInterpolate)                :: self
    type(inputParameters                ), intent(inout) :: parameters
     
    self=nodeOperatorHaloSpinInterpolate()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloSpinInterpolateConstructorParameters

  subroutine haloSpinInterpolateNodeInitialize(self,node)
    !!{
    Compute the rate of growth of halo spin assuming a constant growth rate.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpin
    implicit none
    class           (nodeOperatorHaloSpinInterpolate), intent(inout)          :: self
    type            (treeNode                       ), intent(inout), target  :: node
    class           (nodeComponentBasic             )               , pointer :: basic       , basicParent
    class           (nodeComponentSpin              )               , pointer :: spin        , spinParent
    double precision                                                          :: timeInterval
    logical                                                                   :: spinIsVector
    
    if (node%isPrimaryProgenitor()) then
       ! For primary progenitors compute and store the spin growth rate.
       basic        =>  node              %basic                         ()
       basicParent  =>  node       %parent%basic                         ()
       spin         =>  node              %spin                          ()
       spinParent   =>  node       %parent%spin                          ()
       spinIsVector =   spin              %spinVectorGrowthRateIsSettable()
       timeInterval =  +basicParent       %time                          () &
            &          -basic             %time                          ()
       if (timeInterval > 0.0d0) then
          call        spin%spinGrowthRateSet      (                           &
               &                                   +(                         &
               &                                     +spinParent%spin      () &
               &                                     -spin      %spin      () &
               &                                    )                         &
               &                                   /timeInterval              &
               &                                  )
          if (spinIsVector)                                                   &
               & call spin%spinVectorGrowthRateSet(                           &
               &                                   +(                         &
               &                                     +spinParent%spinVector() &
               &                                     -spin      %spinVector() &
               &                                    )                         &
               &                                   /timeInterval              &
               &                                  )
       else
          call        spin%spinGrowthRateSet      (                           &
               &                                    +0.0d0                    &
               &                                  )
          if (spinIsVector)                                                   &
               & call spin%spinVectorGrowthRateSet(                           &
               &                                   [+0.0d0,+0.0d0,+0.0d0]     &
               &                                  )
       end if
    else
       ! For non-primary progenitors, assume no growth.
       call spin%spinGrowthRateSet(0.0d0)
    end if
    return
  end subroutine haloSpinInterpolateNodeInitialize
  
  subroutine haloSpinInterpolateDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Evolve scalar halo spin at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpin, propertyTypeInactive
    implicit none
    class    (nodeOperatorHaloSpinInterpolate), intent(inout), target  :: self
    type     (treeNode                       ), intent(inout)          :: node
    logical                                   , intent(inout)          :: interrupt
    procedure(interruptTask                  ), intent(inout), pointer :: functionInterrupt
    integer                                   , intent(in   )          :: propertyType
    class    (nodeComponentSpin              )               , pointer :: spin
    logical                                                            :: spinIsVector

    if (propertyType == propertyTypeInactive) return
    spin         => node%spin                          ()
    spinIsVector =  spin%spinVectorGrowthRateIsSettable()
    if (spinIsVector) then
       call        spin%spinVectorRate(                                  &
            &                                spin%spinVectorGrowthRate() &
            &                         )
       if (spin%spin() > 0.0d0)                                          &
            & call spin%spinRate      (                                  &
            &                          +sum(                             &
            &                               +spin%spinVector          () &
            &                               *spin%spinVectorGrowthRate() &
            &                              )                             &
            &                          /     spin%spin                () &
            &                         )
    else
       call        spin%spinRate      (                                  &
            &                                spin%spinGrowthRate      () &
            &                         )
    end if
    return
  end subroutine haloSpinInterpolateDifferentialEvolution
