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

!% An implementation of the hot halo mass distribution class for cored isothermal distributions.

  !# <hotHaloMassDistribution name="hotHaloMassDistributionCoredIsothermal">
  !#  <description>Provides a cored isothermal implementation of the hot halo mass distribution class.</description>
  !# </hotHaloMassDistribution>
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionCoredIsothermal
     !% A coredIsothermal implementation of the hot halo mass distribution class.
     private
     type(massDistributionBetaProfile) :: distribution
   contains
     !@ <objectMethods>
     !@   <object>hotHaloMassDistributionCoredIsothermal</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <type>void</type>
     !@     <arguments>\textcolor{red}{\textless *type(treeNode)\textgreater} node\argin</arguments>
     !@     <description>Initialize the cored isothermal density hot halo mass distribution for the given {\tt node}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                          coredIsothermalDestructor
     procedure :: initialize            => coredIsothermalInitialize
     procedure :: density               => coredIsothermalDensity
     procedure :: densityLogSlope       => coredIsothermalDensityLogSlope
     procedure :: enclosedMass          => coredIsothermalEnclosedMass
     procedure :: radialMoment          => coredIsothermalRadialMoment
     procedure :: rotationNormalization => coredIsothermalRotationNormalization
  end type hotHaloMassDistributionCoredIsothermal

  interface hotHaloMassDistributionCoredIsothermal
     !% Constructors for the cored isothermal hot halo mass distribution class.
     module procedure coredIsothermalDefaultConstructor
  end interface hotHaloMassDistributionCoredIsothermal

  logical :: coredIsothermalInitialized=.false.

contains

  function coredIsothermalDefaultConstructor()
    !% Default constructor for the coredIsothermal hot halo mass distribution class.
    use Galacticus_Error
    use Array_Utilities
    implicit none
    type(hotHaloMassDistributionCoredIsothermal) :: coredIsothermalDefaultConstructor

    if (.not.coredIsothermalInitialized) then
       !$omp critical(coredIsothermalInitialized)
       if (.not.coredIsothermalInitialized) then
          ! Check that required propert is gettable.
          if     (                                                                                                      &
               &  .not.(                                                                                                &
               &         defaultHotHaloComponent%       massIsGettable()                                                &
               &        .and.                                                                                           &
               &         defaultHotHaloComponent%outerRadiusIsGettable()                                                &
               &       )                                                                                                &
               & ) call Galacticus_Error_Report                                                                         &
               & (                                                                                                      &
               &  'coredIsothermalDefaultConstructor'                                                                 , &
               &  'This method requires that the "mass" property of the hot halo is gettable.'//                        &
               &  Galacticus_Component_List(                                                                            &
               &                            'hotHalo'                                                                 , &
               &                             defaultHotHaloComponent%       massAttributeMatch(requireGettable=.true.)  &
               &                            .intersection.                                                              &
               &                             defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.)  &
               &                           )                                                                            &
               & )
          ! Record that implementation is now initialized.
          coredIsothermalInitialized=.true.
       end if
       !$omp end critical(coredIsothermalInitialized)
    end if
    return
  end function coredIsothermalDefaultConstructor

  elemental subroutine coredIsothermalDestructor(self)
    !% Destructor for the single-coredIsothermal hot halo mass distribution class.
    implicit none
    type(hotHaloMassDistributionCoredIsothermal), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine coredIsothermalDestructor

  subroutine coredIsothermalInitialize(self,node)
    !% Initialize the cored isothermal hot halo density profile for the given {\tt node}.
    use Hot_Halo_Mass_Distributions_Core_Radii
    implicit none
    class           (hotHaloMassDistributionCoredIsothermal), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo
    class           (hotHaloMassDistributionCoreRadiusClass)               , pointer :: hotHaloMassDistributionCoreRadius_
    double precision                                        , parameter              :: betaIsothermal=2.0d0/3.0d0
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
    call self%                                     &
         & distribution%                           &
         &  initialize(                            &
         &             beta       =betaIsothermal, &
         &             coreRadius =radiusScale   , &
         &             mass       =mass          , &
         &             outerRadius=radiusOuter     &
         &            )
    return
  end subroutine coredIsothermalInitialize

  double precision function coredIsothermalDensity(self,node,radius)
    !% Return the density in a single-coredIsothermal hot halo mass distribution.
    use Coordinates
    implicit none
    class           (hotHaloMassDistributionCoredIsothermal), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    double precision                                        , intent(in   )          :: radius
    type            (coordinateSpherical                   )                         :: position

    call self%initialize(node)
    position              =[radius,0.0d0,0.0d0]
    coredIsothermalDensity=self%distribution%density(position)
    return
  end function coredIsothermalDensity

  double precision function coredIsothermalDensityLogSlope(self,node,radius)
    !% Return the logarithmic slope of the density of the hot halo at the given {\tt radius}.
    use Coordinates
    implicit none
    class           (hotHaloMassDistributionCoredIsothermal), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    double precision                                        , intent(in   )          :: radius
    type            (coordinateSpherical                   )                         :: position

    call self%initialize(node)
    position                      =[radius,0.0d0,0.0d0]
    coredIsothermalDensityLogSlope=self%distribution%densityGradientRadial(position,logarithmic=.true.)
    return
  end function coredIsothermalDensityLogSlope
  
  double precision function coredIsothermalEnclosedMass(self,node,radius)
    !% Return the mass enclosed in the hot halo at the given {\tt radius}.
    implicit none
    class           (hotHaloMassDistributionCoredIsothermal), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    double precision                                        , intent(in   )          :: radius
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo

    hotHalo => node%hotHalo()
    if (radius > hotHalo%outerRadius()) then
       coredIsothermalEnclosedMass=hotHalo%mass()
    else
       call self%initialize(node)
       coredIsothermalEnclosedMass=self%distribution%massEnclosedBySphere(radius)
    end if
    return
  end function coredIsothermalEnclosedMass
  
  double precision function coredIsothermalRadialMoment(self,node,moment,radius)
    !% Return the radial moment of the density profile of the hot halo to the given {\tt radius}.
    implicit none
    class           (hotHaloMassDistributionCoredIsothermal), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    double precision                                        , intent(in   )          :: moment , radius
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo

    call self%initialize(node)
    hotHalo                    => node%hotHalo()
    coredIsothermalRadialMoment=                                          &
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
  end function coredIsothermalRadialMoment

  double precision function coredIsothermalRotationNormalization(self,node)
    !% Returns the relation between specific angular momentum and rotation velocity (assuming a
    !% rotation velocity that is constant in radius) for {\tt node}. Specifically, the
    !% normalization, $A$, returned is such that $V_{\rm rot} = A J/M$.
    implicit none
    class(hotHaloMassDistributionCoredIsothermal), intent(inout)          :: self
    type (treeNode                              ), intent(inout), pointer :: node
    class(nodeComponentHotHalo                  )               , pointer :: hotHalo

    hotHalo                             => node%hotHalo()
    coredIsothermalRotationNormalization=                       &
         &  self%radialMoment(node,2.0d0,hotHalo%outerRadius()) &
         & /self%radialMoment(node,3.0d0,hotHalo%outerRadius())
    return
  end function coredIsothermalRotationNormalization
