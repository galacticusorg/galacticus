!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
An implementation of the hot halo mass distribution class for $\beta$-profile distributions.
!!}

  use :: Hot_Halo_Mass_Distributions_Core_Radii, only : hotHaloMassDistributionCoreRadiusClass
  use :: Mass_Distributions                    , only : massDistributionBetaProfile

  !![
  <hotHaloMassDistribution name="hotHaloMassDistributionBetaProfile">
   <description>
    A hot halo mass distribution class which adopts a spherically symmetric $\beta$-profile density profile for the hot
    halo. Specifically,
    \begin{equation}
     \rho_\mathrm{hot halo}(r) \propto \left[ r^2 + r_\mathrm{core}^2 \right]^{3\beta/2},
    \end{equation}
    where the core radius, $r_\mathrm{core}$, is set using the selected cored profile core radius method (see
    \refPhysics{hotHaloMassDistributionCoreRadius}). The value of $\beta$ is specified by the {\normalfont
    \ttfamily [beta]} parameter. The profile is normalized such that the current mass in the hot gas profile is contained
    within the outer radius of the hot halo, $r_\mathrm{hot, outer}$.
   </description>
  </hotHaloMassDistribution>
  !!]
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionBetaProfile
     !!{
     A $\beta$-profile implementation of the hot halo mass distribution class.
     !!}
     private
     double precision                                                  :: beta
     type            (massDistributionBetaProfile           )          :: distribution
     class           (hotHaloMassDistributionCoreRadiusClass), pointer :: hotHaloMassDistributionCoreRadius_ => null()
   contains
     !![
     <methods>
       <method description="Initialize the $\beta$-profile density hot halo mass distribution for the given {\normalfont \ttfamily node}." method="initialize" />
     </methods>
     !!]
     final     ::                           betaProfileDestructor
     procedure :: get                    => betaProfileGet
     procedure :: initialize             => betaProfileInitialize
     procedure :: density                => betaProfileDensity
     procedure :: densityLogSlope        => betaProfileDensityLogSlope
     procedure :: enclosedMass           => betaProfileEnclosedMass
     procedure :: radialMoment           => betaProfileRadialMoment
     procedure :: densitySquaredIntegral => betaProfileDensitySquaredIntegral
     procedure :: rotationNormalization  => betaProfileRotationNormalization
  end type hotHaloMassDistributionBetaProfile

  interface hotHaloMassDistributionBetaProfile
     !!{
     Constructors for the $\beta$-profile hot halo mass distribution class.
     !!}
     module procedure betaProfileConstructorParameters
     module procedure betaProfileConstructorInternal
  end interface hotHaloMassDistributionBetaProfile

contains

  function betaProfileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the null {\normalfont \ttfamily betaProfile} hot halo mass distributionclass which builds the object from a
    parameter set.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List          , Error_Report
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    implicit none
    type            (hotHaloMassDistributionBetaProfile    )                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    logical                                                 , save          :: initialized                       =.false.
    class           (hotHaloMassDistributionCoreRadiusClass), pointer       :: hotHaloMassDistributionCoreRadius_
    double precision                                                        :: beta

    if (.not.initialized) then
       !$omp critical(betaProfileInitialized)
       if (.not.initialized) then
          ! Check that required property is gettable.
          if     (                                                                                            &
               &  .not.(                                                                                      &
               &         defaultHotHaloComponent%       massIsGettable()                                      &
               &        .and.                                                                                 &
               &         defaultHotHaloComponent%outerRadiusIsGettable()                                      &
               &       )                                                                                      &
               & ) call Error_Report                                                                          &
               & (                                                                                            &
               &  'This method requires that the "mass" property of the hot halo is gettable.'//              &
               &  Component_List(                                                                             &
               &                 'hotHalo'                                                                 ,  &
               &                  defaultHotHaloComponent%       massAttributeMatch(requireGettable=.true.)   &
               &                 .intersection.                                                               &
               &                  defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.)   &
               &                )                                                                          // &
               &  {introspection:location}                                                                    &
               & )
          initialized=.true.
       end if
       !$omp end critical(betaProfileInitialized)
    end if

    !![
    <inputParameter>
      <name>beta</name>
      <defaultValue>2.0d0/3.0d0</defaultValue>
      <description>The value of $\beta$ in $\beta$-profile hot halo mass distributions.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="hotHaloMassDistributionCoreRadius" name="hotHaloMassDistributionCoreRadius_" source="parameters"/>
    !!]
    self=hotHaloMassDistributionBetaProfile(beta,hotHaloMassDistributionCoreRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloMassDistributionCoreRadius_"/>
    !!]
    return
  end function betaProfileConstructorParameters

  function betaProfileConstructorInternal(beta,hotHaloMassDistributionCoreRadius_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily betaProfile} hot halo mass distribution class.
    !!}
    implicit none
    type            (hotHaloMassDistributionBetaProfile    )                        :: self
    double precision                                        , intent(in   )         :: beta
    class           (hotHaloMassDistributionCoreRadiusClass), intent(in   ), target :: hotHaloMassDistributionCoreRadius_
    !![
    <constructorAssign variables="beta, *hotHaloMassDistributionCoreRadius_"/>
    !!]

    return
  end function betaProfileConstructorInternal

  subroutine betaProfileDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily betaProfile} hot halo mass distribution class.
    !!}
    implicit none
    type(hotHaloMassDistributionBetaProfile), intent(inout) :: self

    !![
    <objectDestructor name="self%hotHaloMassDistributionCoreRadius_"/>
    !!]
    return
  end subroutine betaProfileDestructor

  function betaProfileGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the $\beta$-profile hot halo mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo, treeNode
    use :: Galactic_Structure_Options, only : componentTypeHotHalo, massTypeGaseous, weightByMass
    implicit none
    class           (massDistributionClass             ), pointer                 :: massDistribution_
    class           (hotHaloMassDistributionBetaProfile), intent(inout)           :: self
    type            (treeNode                          ), intent(inout)           :: node
    type            (enumerationWeightByType           ), intent(in   ), optional :: weightBy
    integer                                             , intent(in   ), optional :: weightIndex
    class           (nodeComponentHotHalo              ), pointer                 :: hotHalo
    double precision                                                              :: radiusScale      , radiusOuter, &
         &                                                                           mass
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Get properties of the hot halo.
    radiusScale =  self   %hotHaloMassDistributionCoreRadius_%radius     (node)
    hotHalo     => node                                      %hotHalo    (    )
    radiusOuter =  hotHalo                                   %outerRadius(    )
    mass        =  hotHalo                                   %mass       (    )
    ! If outer radius is non-positive return a null profile.
    if (radiusOuter <= 0.0d0 .or. mass <= 0.0d0) return
    ! Create the beta-profile distribution.
    allocate(massDistributionBetaProfile :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionBetaProfile)
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionBetaProfile(                                                 &amp;
             &amp;                     beta                 =self%beta                , &amp;
             &amp;                     coreRadius           =     radiusScale         , &amp;
             &amp;                     mass                 =     mass                , &amp;
             &amp;                     outerRadius          =     radiusOuter         , &amp;
             &amp;                     truncateAtOuterRadius=     .true.              , &amp;
             &amp;                     componentType        =     componentTypeHotHalo, &amp;
             &amp;                     massType             =     massTypeGaseous       &amp;
             &amp;                    )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function betaProfileGet

  subroutine betaProfileInitialize(self,node)
    !!{
    Initialize the $\beta$-profile hot halo density profile for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloMassDistributionBetaProfile    ), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    class           (nodeComponentHotHalo                  ), pointer       :: hotHalo
    double precision                                                        :: radiusScale, radiusOuter, &
         &                                                                     mass

    radiusScale =  self%hotHaloMassDistributionCoreRadius_%radius     (node)
    hotHalo     =>      node                              %hotHalo    (    )
    radiusOuter =       hotHalo                           %outerRadius(    )
    mass        =       hotHalo                           %mass       (    )
    if (radiusOuter <= 0.0d0) then
       ! If outer radius is non-positive, set mass to zero and outer radius to an arbitrary value.
       mass=0.0d0
       radiusOuter=1.0d0
    end if
    self%distribution=massDistributionBetaProfile(                           &
         &                                        beta         =self%beta  , &
         &                                        coreRadius   =radiusScale, &
         &                                        mass         =mass       , &
         &                                        outerRadius  =radiusOuter  &
         &                                       )
    return
  end subroutine betaProfileInitialize

  double precision function betaProfileDensity(self,node,radius)
    !!{
    Return the density in a single-betaProfile hot halo mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class           (hotHaloMassDistributionBetaProfile), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                    , intent(in   ) :: radius
    type            (coordinateSpherical               )                :: position

    call self%initialize(node)
    position          =[radius,0.0d0,0.0d0]
    betaProfileDensity=self%distribution%density(position)
    return
  end function betaProfileDensity

  double precision function betaProfileDensityLogSlope(self,node,radius)
    !!{
    Return the logarithmic slope of the density of the hot halo at the given {\normalfont \ttfamily radius}.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class           (hotHaloMassDistributionBetaProfile), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                    , intent(in   ) :: radius
    type            (coordinateSpherical               )                :: position

    call self%initialize(node)
    position                  =[radius,0.0d0,0.0d0]
    betaProfileDensityLogSlope=self%distribution%densityGradientRadial(position,logarithmic=.true.)
    return
  end function betaProfileDensityLogSlope

  double precision function betaProfileEnclosedMass(self,node,radius)
    !!{
    Return the mass enclosed in the hot halo at the given {\normalfont \ttfamily radius}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloMassDistributionBetaProfile), intent(inout)          :: self
    type            (treeNode                          ), intent(inout), target  :: node
    double precision                                    , intent(in   )          :: radius
    class           (nodeComponentHotHalo              )               , pointer :: hotHalo

    hotHalo => node%hotHalo()
    if (radius > hotHalo%outerRadius()) then
       betaProfileEnclosedMass=hotHalo%mass()
    else
       call self%initialize(node)
       betaProfileEnclosedMass=self%distribution%massEnclosedBySphere(radius)
    end if
    return
  end function betaProfileEnclosedMass

  double precision function betaProfileRadialMoment(self,node,moment,radius)
    !!{
    Return the radial moment of the density profile of the hot halo to the given {\normalfont \ttfamily radius}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloMassDistributionBetaProfile), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                    , intent(in   ) :: moment , radius
    class           (nodeComponentHotHalo              ), pointer       :: hotHalo

    call self%initialize(node)
    hotHalo                => node%hotHalo()
    betaProfileRadialMoment=                                              &
         & self%                                                          &
         &  distribution%                                                 &
         &   densityRadialMoment(                                         &
         &                       moment                                 , &
         &                       radiusMinimum=0.0d0                    , &
         &                       radiusMaximum=min(                       &
         &                                         radius               , &
         &                                         hotHalo%outerRadius()  &
         &                                        )                       &
         &                      )
    return
  end function betaProfileRadialMoment

  double precision function betaProfileDensitySquaredIntegral(self,node,radius)
    !!{
    Return the integral of the square of the density profile of the hot halo to the given {\normalfont \ttfamily radius}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloMassDistributionBetaProfile), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                    , intent(in   ) :: radius
    class           (nodeComponentHotHalo              ), pointer       :: hotHalo

    call self%initialize(node)
    hotHalo                => node%hotHalo()
    betaProfileDensitySquaredIntegral=                                      &
         & self%                                                            &
         &  distribution%                                                   &
         &   densitySquareIntegral(                                         &
         &                         radiusMinimum=0.0d0                    , &
         &                         radiusMaximum=min(                       &
         &                                           radius               , &
         &                                           hotHalo%outerRadius()  &
         &                                          )                       &
         &                         )
    return
  end function betaProfileDensitySquaredIntegral

  double precision function betaProfileRotationNormalization(self,node)
    !!{
    Returns the relation between specific angular momentum and rotation velocity (assuming a
    rotation velocity that is constant in radius) for {\normalfont \ttfamily node}. Specifically, the
    normalization, $A$, returned is such that $V_\mathrm{rot} = A J/M$.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class(hotHaloMassDistributionBetaProfile), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node
    class(nodeComponentHotHalo              ), pointer       :: hotHalo

    hotHalo                         => node%hotHalo()
    betaProfileRotationNormalization=                           &
         &  self%radialMoment(node,2.0d0,hotHalo%outerRadius()) &
         & /self%radialMoment(node,3.0d0,hotHalo%outerRadius())
    return
  end function betaProfileRotationNormalization
