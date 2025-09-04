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

    In addition to the mean angular momentum vector of unresolved accretion accounted for by \cite{benson_random-walk_2020}, an
    optional stochastic contribution to the angular momentum from unresolved accretion is allowed for. This represents the fact
    that the angular momentum vector of a halo will diffuse away from zero in a random walk even if the mean angular momentum
    contributed by unresolved accretion is zero. The three components of the angular momentum vector of unresolved accretion are
    treated as independent Wiener processes with time-dependent variance. Specifically, each component of the angular momentum
    vector obeys:
    \begin{equation}
    J_\mathrm{i}(t_2) = J_\mathrm{i}(t_1) + \left[ \sigma^2 \left\{ J_\mathrm{v}^2(t_2) - J_\mathrm{v}^2(t_1) \right\} f_\mathrm{u} \right]^{1/2} N(0,1)
    \end{equation}
    where $J_\mathrm{v}(t) = M_\mathrm{v}(t) V_\mathrm{v}(t) R_\mathrm{v}(t)$ is the characteristic virial angular momentum,
    $M_\mathrm{v}(t)$, $V_\mathrm{v}(t)$, and $R_\mathrm{v}(t)$ are the virial mass, velocity, and radius respectively,
    $M_\mathrm{u}$ is the unresolved mass between times $t_1$ and $t_2$ respectively\footnote{Note that this assumes that the
    characteristic angular momentum scales in proportion to mass. In detail this is not correct, as there is also some dependence
    on the change in redshift across the timestep due to the dependence of virial densities on redshift. In practice, we ignore
    this dependence and absorb such effects into the parameter $\sigma$.}, $\sigma^2$ represents the variance in angular momentum
    per unit increase in $J_\mathrm{v}^2$, and $N(0,1)$ is a random variable distributed as a standard normal. The parameter
    $\sigma=${\normalfont \ttfamily [angularMomentumVarianceSpecific]}.

    The factor $f_\mathrm{u}$ can take on one of two forms. If {\normalfont \ttfamily
    [useOriginalSubresolutionMethod]}={\normalfont \ttfamily false} (which should be preferred) then:
    \begin{equation}
    f_\mathrm{u} = \frac{M_\mathrm{u}}{M(t_2)},
    \end{equation}
    in which case the variance is proportional to the mass in unresolved accretion (i.e. the ``time'' variable in the Wiener
    process is just mass). If {\normalfont \ttfamily [useOriginalSubresolutionMethod]}={\normalfont \ttfamily true} then the
    original form,
    \begin{equation}
    f_\mathrm{u} = \left\{ \frac{M(t_1)+M_\mathrm{r}}{M(t_1)} \right\}^2 = \left\{ \frac{M(t_2)-M_\mathrm{u}}{M(t_1)} \right\}^2
    \end{equation}
    is used. In this, flawed, approach the variance arising from unresolved mass is, incorrectly, not proportional to the
    unresolved mass accretion (and also has the incorrect scaling with mass). This approach should not be used for new
    calculations, and is retained only for backward compatibility.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloAngularMomentumVitvitska2002
     !!{
     A node operator class that initializes halo angular momenta using the model of \cite[][see also
     \protect\citealt{benson_random-walk_2020}]{vitvitska_origin_2002}.     
     !!}
     private
     class           (haloSpinDistributionClass        ), pointer :: haloSpinDistribution_          => null()
     class           (virialOrbitClass                 ), pointer :: virialOrbit_                   => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_           => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_  => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_          => null()
     class           (mergerTreeMassResolutionClass    ), pointer :: mergerTreeMassResolution_      => null()
     double precision                                             :: exponentMass                            , angularMomentumVarianceSpecific
     logical                                                      :: useOriginalSubresolutionMethod
   contains
     final     ::                       haloAngularMomentumVitvitska2002Destructor
     procedure :: nodeTreeInitialize => haloAngularMomentumVitvitska2002NodeTreeInitialize
  end type nodeOperatorHaloAngularMomentumVitvitska2002
  
  interface nodeOperatorHaloAngularMomentumVitvitska2002
     !!{
     Constructors for the \refClass{nodeOperatorHaloAngularMomentumVitvitska2002} node operator class.
     !!}
     module procedure haloAngularMomentumVitvitska2002ConstructorParameters
     module procedure haloAngularMomentumVitvitska2002ConstructorInternal
  end interface nodeOperatorHaloAngularMomentumVitvitska2002
  
contains
  
  function haloAngularMomentumVitvitska2002ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorHaloAngularMomentumVitvitska2002} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorHaloAngularMomentumVitvitska2002)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (haloSpinDistributionClass                   ), pointer       :: haloSpinDistribution_
    class           (virialOrbitClass                            ), pointer       :: virialOrbit_
    class           (darkMatterHaloScaleClass                    ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass           ), pointer       :: darkMatterProfileScaleRadius_
    class           (darkMatterProfileDMOClass                   ), pointer       :: darkMatterProfileDMO_
    class           (mergerTreeMassResolutionClass               ), pointer       :: mergerTreeMassResolution_
    double precision                                                              :: exponentMass                  , angularMomentumVarianceSpecific
    logical                                                                       :: useOriginalSubresolutionMethod

    !![
    <inputParameter>
      <name>exponentMass</name>
      <defaultValue>2.0d0</defaultValue>
      <source>parameters</source>
      <description>The exponent of mass ratio appearing in the orbital angular momentum term in the Vitvitska model.</description>
    </inputParameter>
    <inputParameter>
      <name>angularMomentumVarianceSpecific</name>
      <description>The variance in the difference in the angular momentum of a halo per unit angular momentum scale growth.</description>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>useOriginalSubresolutionMethod</name>
      <description>If true, the original (flawed) method for accounting for variance in the angular momentum of subresolution accretion will be used. Otherwise, the new method is used.</description>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
    </inputParameter>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"          source="parameters"/>
    <objectBuilder class="haloSpinDistribution"         name="haloSpinDistribution_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="virialOrbit"                  name="virialOrbit_"                  source="parameters"/>
    <objectBuilder class="mergerTreeMassResolution"     name="mergerTreeMassResolution_"     source="parameters"/>
    !!]
    self=nodeOperatorHaloAngularMomentumVitvitska2002(exponentMass,angularMomentumVarianceSpecific,useOriginalSubresolutionMethod,darkMatterProfileScaleRadius_,haloSpinDistribution_,darkMatterHaloScale_,darkMatterProfileDMO_,virialOrbit_,mergerTreeMassResolution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloSpinDistribution_"        />
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="virialOrbit_"                 />
    <objectDestructor name="mergerTreeMassResolution_"    />
    !!]
    return
  end function haloAngularMomentumVitvitska2002ConstructorParameters

  function haloAngularMomentumVitvitska2002ConstructorInternal(exponentMass,angularMomentumVarianceSpecific,useOriginalSubresolutionMethod,darkMatterProfileScaleRadius_,haloSpinDistribution_,darkMatterHaloScale_,darkMatterProfileDMO_,virialOrbit_,mergerTreeMassResolution_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorHaloAngularMomentumVitvitska2002} node operator class.
    !!}
    use :: Error           , only : Component_List      , Error_Report
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none
    type            (nodeOperatorHaloAngularMomentumVitvitska2002)                        :: self
    class           (haloSpinDistributionClass                   ), intent(in   ), target :: haloSpinDistribution_
    class           (virialOrbitClass                            ), intent(in   ), target :: virialOrbit_
    class           (darkMatterHaloScaleClass                    ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass           ), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (darkMatterProfileDMOClass                   ), intent(in   ), target :: darkMatterProfileDMO_
    class           (mergerTreeMassResolutionClass               ), intent(in   ), target :: mergerTreeMassResolution_
    double precision                                              , intent(in   )         :: exponentMass                  , angularMomentumVarianceSpecific
    logical                                                       , intent(in   )         :: useOriginalSubresolutionMethod
    !![
    <constructorAssign variables="exponentMass, angularMomentumVarianceSpecific, useOriginalSubresolutionMethod, *darkMatterProfileScaleRadius_, *haloSpinDistribution_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *virialOrbit_, *mergerTreeMassResolution_"/>
    !!]

    ! Ensure that the spin component supports vector angular momentum.
    if     (                                                                                                                               &
         &  .not.                                                                                                                          &
         &   (                                                                                                                             &
         &     defaultSpinComponent%angularMomentumVectorIsGettable()                                                                      &
         &    .and.                                                                                                                        &
         &     defaultSpinComponent%angularMomentumVectorIsSettable()                                                                      &
         &   )                                                                                                                             &
         & )                                                                                                                               &
         & call Error_Report                                                                                                               &
         &      (                                                                                                                          &
         &       'method requires a spin component that provides a gettable and settable "angularMomentumVector" property.'             // &
         &       Component_List(                                                                                                           &
         &                      'spin'                                                                                                  ,  &
         &                       defaultSpinComponent%angularMomentumVectorAttributeMatch(requireGettable=.true.,requireSettable=.true.)   &
         &                     )                                                                                                        // &
         &       {introspection:location}                                                                                                  &
         &      )
    return
  end function haloAngularMomentumVitvitska2002ConstructorInternal

  subroutine haloAngularMomentumVitvitska2002Destructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorHaloAngularMomentumVitvitska2002} node operator class.
    !!}
    implicit none
    type(nodeOperatorHaloAngularMomentumVitvitska2002), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    <objectDestructor name="self%haloSpinDistribution_"        />
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
    implicit none
    class           (nodeOperatorHaloAngularMomentumVitvitska2002), intent(inout), target    :: self
    type            (treeNode                                    ), intent(inout), target    :: node
    type            (treeNode                                    )               , pointer   :: nodeChild                         , nodeSibling                       , &
         &                                                                                      nodeUnresolved
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
         &                                                                                      a                                 , b                                 , &
         &                                                                                      angularMomentumScale              , angularMomentumScaleChild         , &
         &                                                                                      angularMomentumRootVariance       , factorMassUnresolved
    integer                                                                                  :: i
    
    basic => node%basic(                 )
    spin  => node%spin (autoCreate=.true.)
    ! If this node has no children, draw its spin from a distribution, and assign a direction which is isotropically
    ! distributed.
    if (.not.associated(node%firstChild)) then
       theta               =acos(2.0d0   *node%hostTree%randomNumberGenerator_%uniformSample()-1.0d0)
       phi                 =     2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()
       angularMomentumValue=self%haloSpinDistribution_%sample(node)*Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_,self%darkMatterProfileDMO_)
       angularMomentumTotal=angularMomentumValue*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
    else
       nodeChild      =>  node  %firstChild
       spinChild      =>  nodeChild %spin      ()
       basicChild     =>  nodeChild %basic     ()
       massUnresolved =  +basic     %mass      () &
            &            -basicChild%mass      ()
       ! Get the resolution of the tree.
       massResolution=self%mergerTreeMassResolution_%resolution(node%hostTree)
       ! Node has multiple progenitors - iterate over them and sum their angular momenta.
       angularMomentumTotal =   0.0d0
       nodeSibling          =>  nodeChild
       do while(associated(nodeSibling%sibling))
          nodeSibling                =>  nodeSibling %sibling
          basicSibling               =>  nodeSibling %basic            (                 )
          spinSibling                =>  nodeSibling %spin             (                 )
          satelliteSibling           =>  nodeSibling %satellite        (autoCreate=.true.)
          massRatio                  =  +basicSibling%mass             (                 ) &
               &                        /basicChild  %mass             (                 )
          massUnresolved             =  +             massUnresolved                       &
               &                        -basicSibling%mass             (                 ) &
               &                        *nodeSibling %subsamplingWeight(                 )
          ! Add orbital angular momentum of this sibling scaled by the reduced mass to correct to the center of mass of the
          ! sibling-child binary system.                
          angularMomentumTotal=+angularMomentumTotal                &
               &               +angularMomentumOrbital(nodeSibling) &
               &               /(                                   &
               &                 +1.0d0                             &
               &                 +massRatio                         &
               &                )**self%exponentMass                &
               &               *nodeSibling%subsamplingWeight()
          ! Add the spin angular momentum of the sibling.
          angularMomentumTotal=+            angularMomentumTotal    &
               &               +spinSibling%angularMomentumVector() &
               &               *nodeSibling%subsamplingWeight   ()
       end do
       ! Add in the spin angular momentum of the primary child.
       angularMomentumTotal=+           angularMomentumTotal   &
            &               +spinChild%angularMomentumVector()
       ! Account for unresolved accretion. The assumption is that unresolved accretion has the mean specific angular
       ! momentum averaged over the distribution of virial orbits.
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
       ! result is identically 1 - in this case we avoid computing β functions.
       massRatio=+basicUnresolved%mass()       &
            &    /basicChild     %mass()
       a        =+2.0d0                        &
            &    -massFunctionSlopeLogarithmic
       if (self%exponentMass == 0.0d0) then
          angularMomentumSubresolutionFactor=+1.0d0
       else
          b                                 =+massFunctionSlopeLogarithmic                                         &
               &                             +self%exponentMass                                                    &
               &                             -2.0d0
          angularMomentumSubresolutionFactor=+Beta_Function_Incomplete_Normalized(a,b,massRatio/(1.0d0+massRatio)) &
               &                             *Beta_Function                      (a,b                            ) &
               &                             *           (2.0d0-massFunctionSlopeLogarithmic)                      &
               &                             /massRatio**(2.0d0-massFunctionSlopeLogarithmic)
       end if
       ! Accumulate the mean angular momentum of the unresolved mass.
       angularMomentumUnresolved=+                  massUnresolved                                               &
            &                    *self%virialOrbit_%angularMomentumVectorMean         (nodeUnresolved,nodeChild) &
            &                    *                  angularMomentumSubresolutionFactor      
       angularMomentumTotal     =+                  angularMomentumTotal                                         &
            &                    +                  angularMomentumUnresolved
       call nodeUnresolved%destroy()
       deallocate(nodeUnresolved)
       ! Account for stochastic variations in the angular momentum contributed by the unresolved mass.
       angularMomentumScale       =+self      %darkMatterHaloScale_%radiusVirial  (node     ) &
            &                      *self      %darkMatterHaloScale_%velocityVirial(node     ) &
            &                      *basic                          %mass          (         )
       angularMomentumScaleChild  =+self      %darkMatterHaloScale_%radiusVirial  (nodeChild) &
            &                      *self      %darkMatterHaloScale_%velocityVirial(nodeChild) &
            &                      *basicChild                     %mass          (         )
       if (self%useOriginalSubresolutionMethod) then
          ! Use the original, flawed, calculation of the unresolved mass factor.
          factorMassUnresolved=+(                              &
               &                 +basic   %mass          ()    &
               &                 -         massUnresolved      &
               &                )                          **2 &
               &               /basicChild%mass          ()**2
       else
          ! Use the new calculation of the unresolved mass factor.
          factorMassUnresolved=+massUnresolved                 &
               &               /basic     %mass          ()
       end if
       angularMomentumRootVariance=+sqrt(                                      &
            &                            +self%angularMomentumVarianceSpecific &
            &                            *(                                    &
            &                              +max(                               &
            &                                   +angularMomentumScale     **2  &
            &                                   -angularMomentumScaleChild**2, &
            &                                   +0.0d0                         &
            &                                  )                               &
            &                              *factorMassUnresolved               &
            &                             )                                    &
            &                           )
       do i=1,3
          angularMomentumTotal(i)=+                                     angularMomentumTotal       (i) &
               &                  +                                     angularMomentumRootVariance    &
               &                  *node%hostTree%randomNumberGenerator_%standardNormalSample       ( )
       end do
       ! Compute the magnitude of the angular momentum.
       angularMomentumValue=Vector_Magnitude(angularMomentumTotal)
    end if
    call spin%angularMomentumSet      (angularMomentumValue)
    call spin%angularMomentumVectorSet(angularMomentumTotal)
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
