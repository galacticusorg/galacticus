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
  Implements a node operator class that causes halo angular momentum be interpolated linearly between child and parent nodes.
  !!}

  !![
  <nodeOperator name="nodeOperatorHaloAngularMomentumInterpolate">
   <description>
    A node operator class that causes halo angular momentum be interpolated linearly between child and parent nodes. For primary
    progenitor nodes, if only the scalar angular momentum, $\lambda$, is available then $\dot{J} = (J_{i+1}-J_i)/(t_{i+1}-t_i)$,
    where $J_i$ is the angular momentum of the node in the initialized tree, $J_{i+1}$ is the angular momentum of its parent node,
    and $t_i$ and $t_{i+1}$ are the corresponding times. If vector angular momentum is available the same interpolation is applied
    to each individual component of angular momentum, with the rate of change of the scalar angular momentum computed
    self-consistently. For non-primary progenitors both scalar and vector angular momentum are assumed to be constant,
    i.e. $\dot{J}=0$.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloAngularMomentumInterpolate
     !!{
     A node operator class that causes halo angular momentum be interpolated linearly between child and parent nodes.
     !!}
     private
   contains
     procedure :: nodeInitialize                      => haloAngMomInterpolateNodeInitialize
     procedure :: differentialEvolutionAnalytics      => haloAngMomInterpolateDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionSolveAnalytics => haloAngMomInterpolateDifferentialEvolutionSolveAnalytics
     procedure :: nodePromote                         => haloAngMomInterpolateNodePromote
  end type nodeOperatorHaloAngularMomentumInterpolate
  
  interface nodeOperatorHaloAngularMomentumInterpolate
     !!{
     Constructors for the \refClass{nodeOperatorHaloAngularMomentumInterpolate} node operator class.
     !!}
     module procedure haloAngMomInterpolateConstructorParameters
  end interface nodeOperatorHaloAngularMomentumInterpolate
  
contains
  
  function haloAngMomInterpolateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorHaloAngularMomentumInterpolate} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorHaloAngularMomentumInterpolate)                :: self
    type(inputParameters                           ), intent(inout) :: parameters
     
    self=nodeOperatorHaloAngularMomentumInterpolate()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloAngMomInterpolateConstructorParameters

  subroutine haloAngMomInterpolateNodeInitialize(self,node)
    !!{
    Compute the rate of growth of halo angular momentum assuming a constant growth rate.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpin
    implicit none
    class           (nodeOperatorHaloAngularMomentumInterpolate), intent(inout), target  :: self
    type            (treeNode                                  ), intent(inout), target  :: node
    class           (nodeComponentBasic                        )               , pointer :: basic                  , basicParent
    class           (nodeComponentSpin                         )               , pointer :: spin                   , spinParent
    double precision                                                                     :: timeInterval
    logical                                                                              :: angularMomentumIsVector
    
    spin                    => node%spin                                     ()
    angularMomentumIsVector =  spin%angularMomentumVectorGrowthRateIsSettable()
    if (node%isPrimaryProgenitor()) then
       ! For primary progenitors compute and store the angular momentum growth rate.
       basic        =>  node              %basic()
       basicParent  =>  node       %parent%basic()
       spinParent   =>  node       %parent%spin ()
       timeInterval =  +basicParent       %time () &
            &          -basic             %time ()
       if (timeInterval > 0.0d0) then
          call        spin%angularMomentumGrowthRateSet      (                                      &
               &                                              +(                                    &
               &                                                +spinParent%angularMomentum      () &
               &                                                -spin      %angularMomentum      () &
               &                                               )                                    &
               &                                              /timeInterval                         &
               &                                             )
          if (angularMomentumIsVector)                                                              &
               & call spin%angularMomentumVectorGrowthRateSet(                                      &
               &                                              +(                                    &
               &                                                +spinParent%angularMomentumVector() &
               &                                                -spin      %angularMomentumVector() &
               &                                               )                                    &
               &                                              /timeInterval                         &
               &                                             )
       else
          call        spin%angularMomentumGrowthRateSet      (                                      &
               &                                               +0.0d0                               &
               &                                             )
          if (angularMomentumIsVector)                                                              &
               & call spin%angularMomentumVectorGrowthRateSet(                                      &
               &                                              [+0.0d0,+0.0d0,+0.0d0]                &
               &                                             )
       end if
    else
       ! For non-primary progenitors, assume no growth.
       call           spin%angularMomentumGrowthRateSet      (                                      &
            &                                                  +0.0d0                               &
            &                                                )
       if (angularMomentumIsVector)                                                                 &
            &    call spin%angularMomentumVectorGrowthRateSet(                                      &
            &                                                 [+0.0d0,+0.0d0,+0.0d0]                &
            &                                                 )
    end if
    return
  end subroutine haloAngMomInterpolateNodeInitialize
  
  subroutine haloAngMomInterpolateDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpin
    implicit none
    class  (nodeOperatorHaloAngularMomentumInterpolate), intent(inout) :: self
    type   (treeNode                                  ), intent(inout) :: node
    class  (nodeComponentSpin                         ), pointer       :: spin
    logical                                                            :: angularMomentumIsVector

    spin                    => node%spin                                     ()
    angularMomentumIsVector =  spin%angularMomentumVectorGrowthRateIsSettable()
    call        spin%angularMomentumAnalytic      ()
    if (angularMomentumIsVector)                     &
         & call spin%angularMomentumVectorAnalytic()
    return
  end subroutine haloAngMomInterpolateDifferentialEvolutionAnalytics

  subroutine haloAngMomInterpolateDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Evolve halo angular momentum at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpin
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodeOperatorHaloAngularMomentumInterpolate), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    double precision                                            , intent(in   ) :: time
    class           (nodeComponentBasic                        ), pointer       :: basicParent
    class           (nodeComponentSpin                         ), pointer       :: spin                   , spinParent
    logical                                                                     :: angularMomentumIsVector
    !$GLC attributes unused :: self

    if (.not.node%isPrimaryProgenitor()) return
    spin                    => node       %spin                                     ()
    spinParent              => node%parent%spin                                     ()
    basicParent             => node%parent%basic                                    ()
    angularMomentumIsVector =  spin       %angularMomentumVectorGrowthRateIsSettable()
    if (angularMomentumIsVector) then
       call spin%angularMomentumVectorSet(                                                               &
            &                             +                spinParent %angularMomentumVector          () &
            &                             +(                                                             &
            &                               +                          time                              &
            &                               -              basicParent%time                           () &
            &                              )                                                             &
            &                             *                spin       %angularMomentumVectorGrowthRate() &
            &                            )
       call spin%angularMomentumSet      (                                                               &
            &                             Vector_Magnitude(                                              &
            &                                              spin       %angularMomentumVector          () &
            &                                             )                                              &
            &                            )
    else
       call spin%angularMomentumSet      (                                                               &
            &                             +                spinParent %angularMomentum                () &
            &                             +(                                                             &
            &                               +                          time                              &
            &                               -              basicParent%time                           () &
            &                              )                                                             &
            &                             *                spin       %angularMomentumGrowthRate      () &
            &                            )
    end if
    return
  end subroutine haloAngMomInterpolateDifferentialEvolutionSolveAnalytics

  subroutine haloAngMomInterpolateNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the angular momentum of {\normalfont \ttfamily node}
    to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpin, treeNode
    implicit none
    class(nodeOperatorHaloAngularMomentumInterpolate), intent(inout)  :: self
    type (treeNode                                  ), intent(inout)  :: node
    type (treeNode                                  ), pointer        :: nodeParent
    class(nodeComponentSpin                         ), pointer        :: spinParent , spin
    !$GLC attributes unused :: self

    nodeParent => node      %parent
    spin       => node      %spin  ()
    spinParent => nodeParent%spin  ()
    call    spin%angularMomentumSet                (spinParent%angularMomentum                ())
    call    spin%angularMomentumGrowthRateSet      (spinParent%angularMomentumGrowthRate      ())
    if (spin%angularMomentumVectorGrowthRateIsSettable()) then
       call spin%angularMomentumVectorSet          (spinParent%angularMomentumVector          ())
       call spin%angularMomentumVectorGrowthRateSet(spinParent%angularMomentumVectorGrowthRate())
    end if
    return
  end subroutine haloAngMomInterpolateNodePromote
