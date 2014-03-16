!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% An implementation of the hot halo mass distribution class for $\beta$-profile distributions.

  double precision :: betaProfileBeta
  logical          :: betaProfileInitialized=.false.

  !# <hotHaloMassDistribution name="hotHaloMassDistributionBetaProfile">
  !#  <description>Provides a $\beta$-profile implementation of the hot halo mass distribution class.</description>
  !# </hotHaloMassDistribution>
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionBetaProfile
     !% A $\beta$-profile implementation of the hot halo mass distribution class.
     private
     double precision                              :: beta
     type            (massDistributionBetaProfile) :: distribution
   contains
     !@ <objectMethods>
     !@   <object>hotHaloMassDistributionBetaProfile</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <type>void</type>
     !@     <arguments>\textcolor{red}{\textless *type(treeNode)\textgreater} node\argin</arguments>
     !@     <description>Initialize the $\beta$-profile density hot halo mass distribution for the given {\tt node}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                          betaProfileDestructor
     procedure :: initialize            => betaProfileInitialize
     procedure :: density               => betaProfileDensity
     procedure :: densityLogSlope       => betaProfileDensityLogSlope
     procedure :: enclosedMass          => betaProfileEnclosedMass
     procedure :: radialMoment          => betaProfileRadialMoment
     procedure :: rotationNormalization => betaProfileRotationNormalization
  end type hotHaloMassDistributionBetaProfile

  interface hotHaloMassDistributionBetaProfile
     !% Constructors for the $\beta$-profile hot halo mass distribution class.
     module procedure betaProfileConstructor
     module procedure betaProfileDefaultConstructor
  end interface hotHaloMassDistributionBetaProfile

contains

  function betaProfileDefaultConstructor()
    !% Default constructor for the betaProfile hot halo mass distribution class.
    use Galacticus_Error
    use Array_Utilities
    use Input_Parameters
    implicit none
    type(hotHaloMassDistributionBetaProfile) :: betaProfileDefaultConstructor

    if (.not.betaProfileInitialized) then
       !$omp critical(betaProfileInitialized)
       if (.not.betaProfileInitialized) then
          ! Check that required propert is gettable.
          if     (                                                                                                      &
               &  .not.(                                                                                                &
               &         defaultHotHaloComponent%       massIsGettable()                                                &
               &        .and.                                                                                           &
               &         defaultHotHaloComponent%outerRadiusIsGettable()                                                &
               &       )                                                                                                &
               & ) call Galacticus_Error_Report                                                                         &
               & (                                                                                                      &
               &  'betaProfileDefaultConstructor'                                                                 , &
               &  'This method requires that the "mass" property of the hot halo is gettable.'//                        &
               &  Galacticus_Component_List(                                                                            &
               &                            'hotHalo'                                                                 , &
               &                             defaultHotHaloComponent%       massAttributeMatch(requireGettable=.true.)  &
               &                            .intersection.                                                              &
               &                             defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.)  &
               &                           )                                                                            &
               & )
          !@ <inputParameter>
          !@   <name>hotHaloMassDistributionBeta</name>
          !@   <defaultValue>$2/3$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The value of $\beta$ in $\beta$-profile hot halo mass distributions.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloMassDistributionBeta',betaProfileBeta,defaultValue=2.0d0/3.0d0)
          ! Record that implementation is now initialized.
          betaProfileInitialized=.true.
       end if
       !$omp end critical(betaProfileInitialized)
    end if
    betaProfileDefaultConstructor=betaProfileConstructor(betaProfileBeta)
    return
  end function betaProfileDefaultConstructor
  
  function betaProfileConstructor(beta)
    !% Default constructor for the {\tt betaProfile} hot halo mass distribution class.
    use Input_Parameters
    implicit none
    type            (hotHaloMassDistributionBetaProfile)                :: betaProfileConstructor
    double precision                                    , intent(in   ) :: beta
    
    betaProfileConstructor%beta=beta
    return
  end function betaProfileConstructor
  
  elemental subroutine betaProfileDestructor(self)
    !% Destructor for the {\tt betaProfile} hot halo mass distribution class.
    implicit none
    type(hotHaloMassDistributionBetaProfile), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine betaProfileDestructor

  subroutine betaProfileInitialize(self,node)
    !% Initialize the $\beta$-profile hot halo density profile for the given {\tt node}.
    use Hot_Halo_Mass_Distributions_Core_Radii
    implicit none
    class           (hotHaloMassDistributionBetaProfile    ), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo
    class           (hotHaloMassDistributionCoreRadiusClass)               , pointer :: hotHaloMassDistributionCoreRadius_
    double precision                                                                 :: radiusScale, radiusOuter, &
         &                                                                              mass

    hotHaloMassDistributionCoreRadius_ => hotHaloMassDistributionCoreRadius             (    )
    radiusScale                        =  hotHaloMassDistributionCoreRadius_%radius     (node)
    hotHalo                            => node                              %hotHalo    (    )
    radiusOuter                        =  hotHalo                           %outerRadius(    )
    mass                               =  hotHalo                           %mass       (    )
    if (radiusOuter <= 0.0d0) then
       ! If outer radius is non-positive, set mass to zero and outer radius to an arbitrary value.
       mass=0.0d0
       radiusOuter=1.0d0
    end if
    call self%                                  &
         & distribution%                        &
         &  initialize(                         &
         &             beta       =self%beta  , &
         &             coreRadius =radiusScale, &
         &             mass       =mass       , &
         &             outerRadius=radiusOuter  &
         &            )
    return
  end subroutine betaProfileInitialize

  double precision function betaProfileDensity(self,node,radius)
    !% Return the density in a single-betaProfile hot halo mass distribution.
    use Coordinates
    implicit none
    class           (hotHaloMassDistributionBetaProfile), intent(inout)          :: self
    type            (treeNode                          ), intent(inout), pointer :: node
    double precision                                    , intent(in   )          :: radius
    type            (coordinateSpherical               )                         :: position

    call self%initialize(node)
    position              =[radius,0.0d0,0.0d0]
    betaProfileDensity=self%distribution%density(position)
    return
  end function betaProfileDensity

  double precision function betaProfileDensityLogSlope(self,node,radius)
    !% Return the logarithmic slope of the density of the hot halo at the given {\tt radius}.
    use Coordinates
    implicit none
    class           (hotHaloMassDistributionBetaProfile), intent(inout)          :: self
    type            (treeNode                          ), intent(inout), pointer :: node
    double precision                                    , intent(in   )          :: radius
    type            (coordinateSpherical               )                         :: position

    call self%initialize(node)
    position                      =[radius,0.0d0,0.0d0]
    betaProfileDensityLogSlope=self%distribution%densityGradientRadial(position,logarithmic=.true.)
    return
  end function betaProfileDensityLogSlope
  
  double precision function betaProfileEnclosedMass(self,node,radius)
    !% Return the mass enclosed in the hot halo at the given {\tt radius}.
    implicit none
    class           (hotHaloMassDistributionBetaProfile), intent(inout)          :: self
    type            (treeNode                          ), intent(inout), pointer :: node
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
    !% Return the radial moment of the density profile of the hot halo to the given {\tt radius}.
    implicit none
    class           (hotHaloMassDistributionBetaProfile), intent(inout)          :: self
    type            (treeNode                          ), intent(inout), pointer :: node
    double precision                                    , intent(in   )          :: moment , radius
    class           (nodeComponentHotHalo              )               , pointer :: hotHalo

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

  double precision function betaProfileRotationNormalization(self,node)
    !% Returns the relation between specific angular momentum and rotation velocity (assuming a
    !% rotation velocity that is constant in radius) for {\tt node}. Specifically, the
    !% normalization, $A$, returned is such that $V_{\rm rot} = A J/M$.
    implicit none
    class(hotHaloMassDistributionBetaProfile), intent(inout)          :: self
    type (treeNode                          ), intent(inout), pointer :: node
    class(nodeComponentHotHalo              )               , pointer :: hotHalo

    hotHalo                         => node%hotHalo()
    betaProfileRotationNormalization=                           &
         &  self%radialMoment(node,2.0d0,hotHalo%outerRadius()) &
         & /self%radialMoment(node,3.0d0,hotHalo%outerRadius())
    return
  end function betaProfileRotationNormalization
