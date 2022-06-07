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
  Implements a node operator class that initializes halo angular momenta using the model of \cite[][see also
  \protect\citealt{benson_random-walk_2020}]{vitvitska_origin_2002}.
  !!}
  
  use :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profile_Scales        , only : darkMatterProfileScaleRadiusClass
  use :: Dark_Matter_Profiles_DMO          , only : darkMatterProfileDMOClass
  use :: Halo_Spin_Distributions           , only : haloSpinDistributionClass
  use :: Virial_Orbits                     , only : virialOrbitClass
  use :: Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolutionClass

  !![
  <nodeOperator name="nodeOperatorHaloAngularMomentumVitvitska2002">
   <description>
    A node operator class that initializes halo angular momenta using the model of \cite[][see also
    \protect\citealt{benson_random-walk_2020}]{vitvitska_origin_2002}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloAngularMomentumVitvitska2002
     !!{
     A node operator class that initializes halo angular momenta using the model of \cite[][see also
     \protect\citealt{benson_random-walk_2020}]{vitvitska_origin_2002}.     
     !!}
     private
     class           (haloSpinDistributionClass        ), pointer :: haloSpinDistribution_         => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_         => null()
     class           (virialOrbitClass                 ), pointer :: virialOrbit_                  => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_          => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     class           (mergerTreeMassResolutionClass    ), pointer :: mergerTreeMassResolution_     => null()
     double precision                                             :: exponentMass
   contains
     final     ::                       haloAngularMomentumVitvitska2002Destructor
     procedure :: nodeTreeInitialize => haloAngularMomentumVitvitska2002NodeTreeInitialize
  end type nodeOperatorHaloAngularMomentumVitvitska2002
  
  interface nodeOperatorHaloAngularMomentumVitvitska2002
     !!{
     Constructors for the {\normalfont \ttfamily haloAngularMomentumVitvitska2002} node operator class.
     !!}
     module procedure haloAngularMomentumVitvitska2002ConstructorParameters
     module procedure haloAngularMomentumVitvitska2002ConstructorInternal
  end interface nodeOperatorHaloAngularMomentumVitvitska2002
  
contains
  
  function haloAngularMomentumVitvitska2002ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily haloAngularMomentumVitvitska2002} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorHaloAngularMomentumVitvitska2002)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (haloSpinDistributionClass                   ), pointer       :: haloSpinDistribution_
    class           (darkMatterProfileDMOClass                   ), pointer       :: darkMatterProfileDMO_
    class           (virialOrbitClass                            ), pointer       :: virialOrbit_
    class           (darkMatterHaloScaleClass                    ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass           ), pointer       :: darkMatterProfileScaleRadius_
    class           (mergerTreeMassResolutionClass               ), pointer       :: mergerTreeMassResolution_
    double precision                                                              :: exponentMass

    !![
    <inputParameter>
      <name>exponentMass</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <description>The exponent of mass ratio appearing in the orbital angular momentum term in the Vitvitska model.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    <objectBuilder class="haloSpinDistribution"         name="haloSpinDistribution_"         source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="virialOrbit"                  name="virialOrbit_"                  source="parameters"/>
    <objectBuilder class="mergerTreeMassResolution"     name="mergerTreeMassResolution_"     source="parameters"/>
    !!]
    self=nodeOperatorHaloAngularMomentumVitvitska2002(exponentMass,darkMatterProfileScaleRadius_,haloSpinDistribution_,darkMatterProfileDMO_,darkMatterHaloScale_,virialOrbit_,mergerTreeMassResolution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloSpinDistribution_"        />
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    <objectDestructor name="virialOrbit_"                 />
    <objectDestructor name="mergerTreeMassResolution_"    />
    !!]
    return
  end function haloAngularMomentumVitvitska2002ConstructorParameters

  function haloAngularMomentumVitvitska2002ConstructorInternal(exponentMass,darkMatterProfileScaleRadius_,haloSpinDistribution_,darkMatterProfileDMO_,darkMatterHaloScale_,virialOrbit_,mergerTreeMassResolution_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily haloAngularMomentumVitvitska2002} node operator class.
    !!}
    implicit none
    type            (nodeOperatorHaloAngularMomentumVitvitska2002)                        :: self
    class           (haloSpinDistributionClass                   ), intent(in   ), target :: haloSpinDistribution_
    class           (darkMatterProfileDMOClass                   ), intent(in   ), target :: darkMatterProfileDMO_
    class           (virialOrbitClass                            ), intent(in   ), target :: virialOrbit_
    class           (darkMatterHaloScaleClass                    ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass           ), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (mergerTreeMassResolutionClass               ), intent(in   ), target :: mergerTreeMassResolution_
    double precision                                              , intent(in   )         :: exponentMass
    !![
    <constructorAssign variables="exponentMass, *darkMatterProfileScaleRadius_, *haloSpinDistribution_, *darkMatterProfileDMO_, *darkMatterHaloScale_, *virialOrbit_, *mergerTreeMassResolution_"/>
    !!]

    return
  end function haloAngularMomentumVitvitska2002ConstructorInternal

  subroutine haloAngularMomentumVitvitska2002Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily haloAngularMomentumVitvitska2002} node operator class.
    !!}
    implicit none
    type(nodeOperatorHaloAngularMomentumVitvitska2002), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    <objectDestructor name="self%haloSpinDistribution_"        />
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    <objectDestructor name="self%darkMatterHaloScale_"         />
    <objectDestructor name="self%virialOrbit_"                 />
    <objectDestructor name="self%mergerTreeMassResolution_"    />
     !!]
    return
  end subroutine haloAngularMomentumVitvitska2002Destructor

  subroutine haloAngularMomentumVitvitska2002NodeTreeInitialize(self,node)
    !!{
    Initialize the spin of {\normalfont \ttfamily node}.
    !!}
    use :: Dark_Matter_Halo_Spins  , only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes        , only : nodeComponentSpin                      , nodeComponentBasic                 , nodeComponentDarkMatterProfile, nodeComponentSatellite
    use :: Beta_Functions          , only : Beta_Function                          , Beta_Function_Incomplete_Normalized
    use :: Numerical_Constants_Math, only : Pi
    use :: Vectors                 , only : Vector_Magnitude
    use :: Merger_Tree_Walkers     , only : mergerTreeWalkerAllNodes
    implicit none
    class           (nodeOperatorHaloAngularMomentumVitvitska2002), intent(inout), target    :: self
    type            (treeNode                                    ), intent(inout), target    :: node
    type            (treeNode                                    )               , pointer   :: nodeChild                         , nodeSibling                       , &
         &                                                                                      nodeUnresolved                    , nodeWork
    class           (nodeComponentBasic                          )               , pointer   :: basicChild                        , basicSibling                      , &
         &                                                                                      basic                             , basicUnresolved
    class           (nodeComponentSpin                           )               , pointer   :: spin                              , spinSibling                       , &
         &                                                                                      spinChild
    class           (nodeComponentSatellite                      )               , pointer   :: satelliteSibling
    class           (nodeComponentDarkMatterProfile              )               , pointer   :: darkMatterProfileUnresolved
    double precision                                              , dimension(3)             :: angularMomentumUnresolved         , angularMomentumTotal
    double precision                                                             , parameter :: massFunctionSlopeLogarithmic=1.9d0
    double precision                                                                         :: angularMomentumValue              , massRatio                         , &
         &                                                                                      theta                             , phi                               , &
         &                                                                                      massUnresolved                    , radiusScaleUnresolved             , &
         &                                                                                      massResolution                    , angularMomentumSubresolutionFactor, &
         &                                                                                      a                                 , b
    type            (mergerTreeWalkerAllNodes                    )                           :: treeWalker

    spin => node%spin(autoCreate=.true.)
    if (spin%angularMomentum() == 0.0d0) then
       ! Perform our own depth-first tree walk to set scales in all nodes of the tree. This is necessary as we require access
       ! to the parent scale to set scale growth rates, but must initialize scales in a strictly depth-first manner as some
       ! algorithms rely on knowing the progenitor structure of the tree to compute scale radii.
       treeWalker=mergerTreeWalkerAllNodes(node%hostTree,spanForest=.true.)
       do while (treeWalker%next(nodeWork))
          basic => nodeWork%basic(                 )
          spin  => nodeWork%spin (autoCreate=.true.)
         ! If this node has no children, draw its spin from a distribution, and assign a direction which is isotropically
          ! distributed.
          if (.not.associated(nodeWork%firstChild)) then
             theta               =acos(2.0d0   *nodeWork%hostTree%randomNumberGenerator_%uniformSample()-1.0d0)
             phi                 =     2.0d0*Pi*nodeWork%hostTree%randomNumberGenerator_%uniformSample()
             angularMomentumValue=self%haloSpinDistribution_%sample(nodeWork)*Dark_Matter_Halo_Angular_Momentum_Scale(nodeWork,self%darkMatterProfileDMO_)
             angularMomentumTotal=angularMomentumValue*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
          else
             nodeChild      =>  nodeWork  %firstChild
             spinChild      =>  nodeChild %spin      ()
             basicChild     =>  nodeChild %basic     ()
             massUnresolved =  +basic     %mass      () &
                  &            -basicChild%mass      ()
             ! Get the resolution of the tree.
             massResolution=self%mergerTreeMassResolution_%resolution(nodeWork%hostTree)
             ! Node has multiple progenitors - iterate over them and sum their angular momenta.
             angularMomentumTotal =   0.0d0
             nodeSibling          =>  nodeChild
             do while(associated(nodeSibling%sibling))
                nodeSibling                =>  nodeSibling %sibling
                basicSibling               =>  nodeSibling %basic         (                 )
                spinSibling                =>  nodeSibling %spin          (                 )
                satelliteSibling           =>  nodeSibling %satellite     (autoCreate=.true.)
                massRatio                  =  +basicSibling%mass          (                 ) &
                     &                        /basicChild  %mass          (                 )
                massUnresolved             =  +             massUnresolved                    &
                     &                        -basicSibling%mass          (                 )
                ! Add orbital angular momentum of this sibling scaled by the reduced mass to correct to the center of mass of the
                ! sibling-child binary system.
                angularMomentumTotal=+angularMomentumTotal                &
                     &               +angularMomentumOrbital(nodeSibling) &
                     &               /(                                   &
                     &                 +1.0d0                             &
                     &                 +massRatio                         &
                     &                )**self%exponentMass
                ! Add the spin angular momentum of the sibling.
                angularMomentumTotal=+            angularMomentumTotal    &
                     &               +spinSibling%angularMomentumVector()
             end do
             ! Add in the spin angular momentum of the primary child.
             angularMomentumTotal=+           angularMomentumTotal   &
                  &               +spinChild%angularMomentumVector()
             ! Account for unresolved accretion. The assumption is that unresolved accretion has the mean specific angular momentum averaged over the distribution of virial orbits.
             nodeUnresolved                       => treeNode                                                      (                         )
             nodeUnresolved             %hostTree => node%hostTree
             basicUnresolved                      => nodeUnresolved                              %basic            (autoCreate=.true.        )
             darkMatterProfileUnresolved          => nodeUnresolved                              %darkMatterProfile(autoCreate=.true.        )
             call basicUnresolved            %massSet            (massResolution       )
             call basicUnresolved            %timeSet            (basicChild%time()    )
             call basicUnresolved            %timeLastIsolatedSet(basicChild%time()    )
             radiusScaleUnresolved                =  self          %darkMatterProfileScaleRadius_%radius           (           nodeUnresolved)
             call darkMatterProfileUnresolved%scaleSet           (radiusScaleUnresolved)
             ! Compute a correction factor to the orbital angular momentum which takes into account the mass dependence of the
             ! 1/(1+m/M)ᵅ term that is applied to the angular momentum, and the reduced mass factor that appears in the orbital
             ! angular momentum. Averaging this over a power-law mass function gives the result below. In the case that α=0 the
             ! result is identically 1 - in this case we avoid computing beta functions.
             massRatio=+basicUnresolved%mass()       &
                  &    /basicChild     %mass()
             a        =+2.0d0                        &
                  &    -massFunctionSlopeLogarithmic
             if (self%exponentMass == 0.0d0) then
                angularMomentumSubresolutionFactor=+1.0d0
             else
                b                                 =+massFunctionSlopeLogarithmic                                         &
                     &                             +self%exponentMAss                                                    &
                     &                             -2.0d0
                angularMomentumSubresolutionFactor=+Beta_Function_Incomplete_Normalized(a,b,massRatio/(1.0d0+massRatio)) &
                     &                             *Beta_Function                      (a,b                            ) &
                     &                             *           (2.0d0-massFunctionSlopeLogarithmic)                      &
                     &                             /massRatio**(2.0d0-massFunctionSlopeLogarithmic)
             end if
             ! Accumulate the angular momentum of the unresolved mass.
             angularMomentumUnresolved=+                  massUnresolved                                               &
                  &                    *self%virialOrbit_%angularMomentumVectorMean         (nodeUnresolved,nodeChild) &
                  &                    *                  angularMomentumSubresolutionFactor      
             angularMomentumTotal     =+                  angularMomentumTotal                                         &
                  &                    +                  angularMomentumUnresolved
             call nodeUnresolved%destroy()
             deallocate(nodeUnresolved)
             ! Compute the magnitude of the angular momentum.
             angularMomentumValue=Vector_Magnitude(angularMomentumTotal)
          end if
          call spin%angularMomentumSet      (angularMomentumValue)
          call spin%angularMomentumVectorSet(angularMomentumTotal)
       end do
    end if
    return
  end subroutine haloAngularMomentumVitvitska2002NodeTreeInitialize

  function angularMomentumOrbital(node)
    !!{
    Returns the orbital angular momentum vector associated with a satellite by drawing a
    random position towards the host at virial radius distance and a random velocity vector
    consistent with the orbital parameters of the satellite.
    !!}
    use :: Coordinates     , only : assignment(=)     , coordinateCartesian
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, treeNode
    use :: Kepler_Orbits   , only : keplerOrbit
    use :: Vectors         , only : Vector_Product
    implicit none
    double precision                        , dimension(3)                :: angularMomentumOrbital
    type            (treeNode              ), pointer     , intent(inout) :: node
    class           (nodeComponentSatellite), pointer                     :: satellite
    class           (nodeComponentBasic    ), pointer                     :: basic
    double precision                        , dimension(3)                :: haloVelocity            , haloPosition
    type            (keplerOrbit           )                              :: orbit
    type            (coordinateCartesian   )                              :: coordinates

    ! Get the orbital properties.
    basic        => node       %basic      (                 )
    satellite    => node       %satellite  (autoCreate=.true.)
    orbit        =  satellite  %virialOrbit(                 )
    coordinates  =  orbit      %position   (                 )
    haloPosition =  coordinates
    coordinates  =  orbit      %velocity   (                 )
    haloVelocity =  coordinates
    ! Calculate the orbital angular momentum vector.
    angularMomentumOrbital=+               basic%mass()  &
         &                 *Vector_Product(              &
         &                                 haloPosition, &
         &                                 haloVelocity  &
         &                                )
    return
  end function angularMomentumOrbital
