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

  !!{RST
  Implements a node operator class that causes halo angular momentum be interpolated linearly between child and parent nodes.
  !!}

  !![
  <nodeOperator name="nodeOperatorHaloAngularMomentumInterpolate" docformat="rst">
   <description>
     A node operator class that causes halo angular momentum be interpolated linearly between child and parent nodes. For
     primary progenitor nodes, if only the scalar angular momentum, :math:`\lambda`, is available then :math:`\dot{J} =
     (J_{i+1}-J_i)/(t_{i+1}-t_i)`, where :math:`J_i` is the angular momentum of the node in the initialized tree,
     :math:`J_{i+1}` is the angular momentum of its parent node, and :math:`t_i` and :math:`t_{i+1}` are the corresponding
     times. If vector angular momentum is available the same interpolation is applied to each individual component of angular
     momentum, with the rate of change of the scalar angular momentum computed self-consistently. For non-primary progenitors
     both scalar and vector angular momentum are assumed to be constant, i.e. :math:`\dot{J}=0`.

     The specific angular momentum of accreted material is also computed and used to set the ``angularMomentumSpecificAccreted``
     property of the spin ``component``. If ``angularMomentumSpecificGrowthFromAccretionRate`` is true then this specific
     angular momentum is computed as :math:`j = \dot{J}/\dot{M}` where :math:`\dot{M}` is the dark matter only accretion rate
     onto the halo (which excludes resolved mergers). If ``angularMomentumSpecificGrowthFromAccretionRate`` is false then this
     specific angular momentum is computed as :math:`j = \Delta J / \Delta M` where :math:`\Delta J` is the change in angular
     momentum from child to parent halo, and :math:`\Delta M` is the change in mass from child to parent halo (and so includes
     mass growth through merging). Since :math:`\dot J` is computed from the full change in angular momentum from child to
     parent (including contributions from mergers) the latter (``angularMomentumSpecificGrowthFromAccretionRate == false``) is
     arguably more physical. The former (``angularMomentumSpecificGrowthFromAccretionRate == true``) is retained for backward
     compatibility.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloAngularMomentumInterpolate
     !!{RST
     A node operator class that causes halo angular momentum be interpolated linearly between child and parent nodes.
     !!}
     private
     logical :: angularMomentumSpecificGrowthFromAccretionRate
   contains
     procedure :: nodeInitialize                      => haloAngMomInterpolateNodeInitialize
     procedure :: differentialEvolutionAnalytics      => haloAngMomInterpolateDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionSolveAnalytics => haloAngMomInterpolateDifferentialEvolutionSolveAnalytics
     procedure :: nodePromote                         => haloAngMomInterpolateNodePromote
  end type nodeOperatorHaloAngularMomentumInterpolate
  
  interface nodeOperatorHaloAngularMomentumInterpolate
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorHaloAngularMomentumInterpolate` node operator class.
     !!}
     module procedure haloAngMomInterpolateConstructorParameters
     module procedure haloAngMomInterpolateConstructorInternal
  end interface nodeOperatorHaloAngularMomentumInterpolate
  
contains
  
  function haloAngMomInterpolateConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorHaloAngularMomentumInterpolate` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorHaloAngularMomentumInterpolate)                :: self
    type   (inputParameters                           ), intent(inout) :: parameters
    logical                                                            :: angularMomentumSpecificGrowthFromAccretionRate
          
    !![
    <inputParameter docformat="rst">
      <name>angularMomentumSpecificGrowthFromAccretionRate</name>
      <defaultValue>.true.</defaultValue>
      <description>
	Controls whether the specific angular momentum of accreted material is estimated using the dark matter only accretion
	rate (``true``; excluding growth due to mergers), or from the next change in dark matter only mass (``false``; includng
	growth due to mergers).
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorHaloAngularMomentumInterpolate(angularMomentumSpecificGrowthFromAccretionRate)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloAngMomInterpolateConstructorParameters
 
  function haloAngMomInterpolateConstructorInternal(angularMomentumSpecificGrowthFromAccretionRate) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorHaloAngularMomentumInterpolate` node operator class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorHaloAngularMomentumInterpolate)                :: self
    logical                                            , intent(in   ) :: angularMomentumSpecificGrowthFromAccretionRate
    !![
    <constructorAssign variables="angularMomentumSpecificGrowthFromAccretionRate"/>
    !!]
    
    return
  end function haloAngMomInterpolateConstructorInternal

  subroutine haloAngMomInterpolateNodeInitialize(self,node)
    !!{RST
    Compute the rate of growth of halo angular momentum assuming a constant growth rate.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpin
    implicit none
    class           (nodeOperatorHaloAngularMomentumInterpolate), intent(inout), target  :: self
    type            (treeNode                                  ), intent(inout), target  :: node
    class           (nodeComponentBasic                        )               , pointer :: basic                  , basicParent
    class           (nodeComponentSpin                         )               , pointer :: spin                   , spinParent
    double precision                                                                     :: timeInterval           , massInterval
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
          call        spin%angularMomentumGrowthRateSet      (                                          &
               &                                              +(                                        &
               &                                                +spinParent%angularMomentum          () &
               &                                                -spin      %angularMomentum          () &
               &                                               )                                        &
               &                                              /timeInterval                             &
               &                                             )
          if (angularMomentumIsVector)                                                                  &
               & call spin%angularMomentumVectorGrowthRateSet(                                          &
               &                                              +(                                        &
               &                                                +spinParent%angularMomentumVector    () &
               &                                                -spin      %angularMomentumVector    () &
               &                                               )                                        &
               &                                              /timeInterval                             &
               &                                             )
          ! Estimate the specific angular momentum of the accreted material.
          if (self%angularMomentumSpecificGrowthFromAccretionRate) then
             ! Estimate using the dark matter-only accretion rate (which excludes resolved mergers). Note that negative
             ! accretion rates are permitted here - when a halo is losing mass the gas removed from the circumgalactic
             ! medium must carry away its angular momentum. Only a precisely zero accretion rate, for which the specific
             ! angular momentum is undefined (and irrelevant, since no mass is accreted), is excluded.
             if (basic%accretionRate() /= 0.0d0) then
                ! Note that the change in the magnitude of the angular momentum is used here directly, rather than the
                ! `angularMomentumGrowthRate` property. For a vector spin component that property is a deferred function
                ! returning the instantaneous projection, d|J|/dt = (J.dJ/dt)/|J|, which varies along the branch segment as
                ! the angular momentum vector reorients. Evaluating it here would freeze its value at the start of the
                ! segment, and, since |J(t)| is convex, would systematically under-deliver angular momentum to the
                ! circumgalactic medium. Using the net change in magnitude instead integrates to exactly the change in the
                ! angular momentum of the halo across the segment, for both scalar and vector spin components.
                call  spin%angularMomentumSpecificAccretedSet(                                          &
                     &                                        +(                                        &
                     &                                          +spinParent%angularMomentum          () &
                     &                                          -spin      %angularMomentum          () &
                     &                                         )                                        &
                     &                                        /             timeInterval                &
                     &                                        /basic       %accretionRate            () &
                     &                                       )
             else
                call  spin%angularMomentumSpecificAccretedSet(                                          &
                     &                                        +0.0d0                                    &
                     &                                       )
             end if
          else
             ! Estimate using the net dark matter-only mass growth (which includes resolved mergers).
             massInterval=+basicParent%mass() &
                  &       -basic      %mass()
             if (massInterval > 0.0d0) then
                call  spin%angularMomentumSpecificAccretedSet(                                          &
                     &                                        +(                                        &
                     &                                          +spinParent%angularMomentum          () &
                     &                                          -spin      %angularMomentum          () &
                     &                                         )                                        &
                     &                                        /massInterval                             &
                     &                                       )
             else
                call  spin%angularMomentumSpecificAccretedSet(                                          &
                     &                                        +0.0d0                                    &
                     &                                       )
             end if
          end if
       else
          call        spin%angularMomentumGrowthRateSet      (                                          &
               &                                              +0.0d0                                    &
               &                                             )
          call        spin%angularMomentumSpecificAccretedSet(                                          &
               &                                              +0.0d0                                    &
               &                                             )
          if (angularMomentumIsVector)                                                                  &
               & call spin%angularMomentumVectorGrowthRateSet(                                          &
               &                                              [+0.0d0,+0.0d0,+0.0d0]                    &
               &                                             )
       end if
    else
       ! For non-primary progenitors, assume no growth.
       call           spin%angularMomentumGrowthRateSet      (                                          &
            &                                                 +0.0d0                                    &
            &                                                )
       call           spin%angularMomentumSpecificAccretedSet(                                          &
            &                                                 +0.0d0                                    &
            &                                                )
       if (angularMomentumIsVector)                                                                     &
            & call spin%angularMomentumVectorGrowthRateSet   (                                          &
            &                                                 [+0.0d0,+0.0d0,+0.0d0]                    &
            &                                                )
    end if
    return
  end subroutine haloAngMomInterpolateNodeInitialize
  
  subroutine haloAngMomInterpolateDifferentialEvolutionAnalytics(self,node)
    !!{RST
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
    !!{RST
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
    !!{RST
    Ensure that ``node`` is ready for promotion to its parent. In this case, we simply update the angular momentum of ``node`` to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpin, treeNode
    implicit none
    class(nodeOperatorHaloAngularMomentumInterpolate), intent(inout) :: self
    type (treeNode                                  ), intent(inout) :: node
    type (treeNode                                  ), pointer       :: nodeParent
    class(nodeComponentSpin                         ), pointer       :: spinParent, spin
    !$GLC attributes unused :: self

    nodeParent => node      %parent
    spin       => node      %spin  ()
    spinParent => nodeParent%spin  ()
    call    spin%angularMomentumSet                (spinParent%angularMomentum                ())
    call    spin%angularMomentumGrowthRateSet      (spinParent%angularMomentumGrowthRate      ())
    call    spin%angularMomentumSpecificAccretedSet(spinParent%angularMomentumSpecificAccreted())
    if (spin%angularMomentumVectorGrowthRateIsSettable()) then
       call spin%angularMomentumVectorSet          (spinParent%angularMomentumVector          ())
       call spin%angularMomentumVectorGrowthRateSet(spinParent%angularMomentumVectorGrowthRate())
    end if
    return
  end subroutine haloAngMomInterpolateNodePromote
