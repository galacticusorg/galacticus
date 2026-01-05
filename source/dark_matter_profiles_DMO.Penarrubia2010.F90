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
  An implementation of \cite{penarrubia_impact_2010} dark matter halo profiles.
  !!}

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOPenarrubia2010">
   <description>
    A dark matter profile DMO class which builds \refClass{massDistributionZhao1996} to implement the
    \cite{penarrubia_impact_2010} density profile.
   </description>
   <deepCopy>
    <functionClass variables="massDistributionStripped, massDistributionUnstripped"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="massDistributionStripped, massDistributionUnstripped"/>
   </stateStorable>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOZhao1996) :: darkMatterProfileDMOPenarrubia2010
     !!{
     A dark matter halo profile class implementing \cite{penarrubia_impact_2010} dark matter halos.
     !!}
     private
     double precision                                   :: betaStripped                                     , muRadius                                       , &
          &                                                etaRadius                                        , muVelocity                                     , &
          &                                                etaVelocity                                      , ratioRadiusMaximumRadiusScaleStripped          , &
          &                                                ratioRadiusMaximumRadiusScaleUnstripped
     type           (massDistributionZhao1996), pointer :: massDistributionStripped                => null(), massDistributionUnstripped            => null()
   contains
     final     ::                  penarrubia2010Destructor
     procedure :: exponents     => penarrubia2010Exponents
     procedure :: scaleRadius   => penarrubia2010ScaleRadius
     procedure :: normalization => penarrubia2010Normalization
  end type darkMatterProfileDMOPenarrubia2010

  interface darkMatterProfileDMOPenarrubia2010
     !!{
     Constructors for the \refClass{darkMatterProfileDMOPenarrubia2010} dark matter halo profile class.
     !!}
     module procedure penarrubia2010ConstructorParameters
     module procedure penarrubia2010ConstructorInternal
  end interface darkMatterProfileDMOPenarrubia2010

  ! Mass fraction at which the transition from "unstripped" to "stripped" profile occurs.
  double precision, parameter :: fractionMassTransition=0.9d0
  
contains

  function penarrubia2010ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily penarrubia2010} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOPenarrubia2010)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_
    double precision                                                    :: alpha               , beta        , &
         &                                                                 gamma               , betaStripped, &
         &                                                                 muRadius            , etaRadius   , &
         &                                                                 muVelocity          , etaVelocity
    
    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <description>The parameter $\alpha$ of the \cite{zhao_analytical_1996} dark matter density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <description>The parameter $\beta$ of the \cite{zhao_analytical_1996} dark matter density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>betaStripped</name>
      <source>parameters</source>
      <defaultValue>5.0d0</defaultValue>
      <description>The parameter $\beta_\mathrm{stripped}$ of the \cite{penarrubia_impact_2010} dark matter density profile. This is the $\beta$ exponent of the \cite{zhao_analytical_1996} dark matter density profile for cases where significant stripping of the profile has occurred.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <source>parameters</source>
      <description>The parameter $\gamma$ of the \cite{zhao_analytical_1996} dark matter density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>muRadius</name>
      <source>parameters</source>
      <description>The parameter $\mu$ of the \cite{penarrubia_impact_2010} tidal track for $r_\mathrm{max}$.</description>
    </inputParameter>
    <inputParameter>
      <name>etaRadius</name>
      <source>parameters</source>
      <description>The parameter $\eta$ of the \cite{penarrubia_impact_2010} tidal track for $r_\mathrm{max}$.</description>
    </inputParameter>
    <inputParameter>
      <name>muVelocity</name>
      <source>parameters</source>
      <description>The parameter $\mu$ of the \cite{penarrubia_impact_2010} tidal track for $V_\mathrm{max}$.</description>
    </inputParameter>
    <inputParameter>
      <name>etaVelocity</name>
      <source>parameters</source>
      <description>The parameter $\eta$ of the \cite{penarrubia_impact_2010} tidal track for $V_\mathrm{max}$.</description>
    </inputParameter>
    <inputParameter>
      <name>etaVelocity</name>
      <source>parameters</source>
      <description>The parameter $\eta$ of the \cite{penarrubia_impact_2010} tidal track for $V_\mathrm{max}$.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOPenarrubia2010(alpha,beta,gamma,betaStripped,muRadius,etaRadius,muVelocity,etaVelocity,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function penarrubia2010ConstructorParameters

  function penarrubia2010ConstructorInternal(alpha,beta,gamma,betaStripped,muRadius,etaRadius,muVelocity,etaVelocity,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOPenarrubia2010} dark matter halo profile class.
    !!}
    use :: Mass_Distributions, only : massDistributionZhao1996
    implicit none
    type            (darkMatterProfileDMOPenarrubia2010)                        :: self
    class           (darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                    , intent(in   )         :: alpha               , beta        , &
         &                                                                         gamma               , betaStripped, &
         &                                                                         muRadius            , etaRadius   , &
         &                                                                         muVelocity          , etaVelocity
    !![
    <constructorAssign variables="alpha, beta, gamma, betaStripped, muRadius, etaRadius, muVelocity, etaVelocity, *darkMatterHaloScale_"/>
    !!]

    ! Compute the mapping between scale radius and radius of peak velocity in the scale-free stripped and unstripped profiles. We
    ! use a halo of arbitrary mass here, but do set a somewhat physical scale radius to ensure that numerical solvers have a
    ! reasonable profile to work with.
    allocate(self%massDistributionStripped  )
    allocate(self%massDistributionUnstripped)
    !![
    <referenceConstruct isResult="yes" owner="self" object="massDistributionStripped"   constructor="massDistributionZhao1996(alpha,betaStripped,gamma,scaleLength=1.0d0,mass=1.0d0,radiusOuter=1.0d0)"/>
    <referenceConstruct isResult="yes" owner="self" object="massDistributionUnstripped" constructor="massDistributionZhao1996(alpha,beta        ,gamma,scaleLength=1.0d0,mass=1.0d0,radiusOuter=1.0d0)"/>
    !!]
    self%ratioRadiusMaximumRadiusScaleStripped  =+self%massDistributionStripped  %radiusRotationCurveMaximum()
    self%ratioRadiusMaximumRadiusScaleUnstripped=+self%massDistributionUnstripped%radiusRotationCurveMaximum()
    return
  end function penarrubia2010ConstructorInternal

  subroutine penarrubia2010Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOPenarrubia2010} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOPenarrubia2010), intent(inout) :: self
    
    !![
    <objectDestructor name="self%massDistributionStripped"  />
    <objectDestructor name="self%massDistributionUnstripped"/>
    !!]
    return
  end subroutine penarrubia2010Destructor
  
  subroutine penarrubia2010Exponents(self,node,alpha,beta,gamma)
    !!{
    Compute the exponents of the {\normalfont \ttfamily penarrubia2010} dark matter halo profile.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentBasic
    implicit none
    class           (darkMatterProfileDMOPenarrubia2010), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                    , intent(  out) :: alpha    , beta, &
         &                                                                 gamma
    class           (nodeComponentBasic                ), pointer       :: basic
    class           (nodeComponentSatellite            ), pointer       :: satellite

    basic     => node%basic    ()
    satellite => node%satellite()
    alpha     =  self%alpha
    gamma     =  self%gamma
    if (satellite%boundMass() < fractionMassTransition*basic%mass()) then
       ! Use exponent for the stripped case.
       beta   =  self%betaStripped
    else
       ! Use exponent for the unstripped case.
       beta   =  self%beta
    end if
    return
  end subroutine penarrubia2010Exponents
  
  double precision function penarrubia2010ScaleRadius(self,node) result(radiusScale)
    !!{
    Compute the scale radius of the {\normalfont \ttfamily penarrubia2010} dark matter halo profile.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (darkMatterProfileDMOPenarrubia2010), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    class           (nodeComponentBasic                ), pointer       :: basic
    class           (nodeComponentSatellite            ), pointer       :: satellite
    class           (nodeComponentDarkMatterProfile    ), pointer       :: darkMatterProfile
    double precision                                                    :: fractionMassBound            , fractionRadiusMaximum, &
         &                                                                 ratioRadiusMaximumRadiusScale
    
    basic                 =>  node     %basic            ()
    satellite             =>  node     %satellite        ()
    darkMatterProfile     =>  node     %darkMatterProfile()
    fractionMassBound     =  +satellite%boundMass        () &
         &                   /basic    %mass             ()
    fractionRadiusMaximum =  +2.0d0              **self%muRadius  &
         &                   *  fractionMassBound**self%etaRadius &
         &                   /(                                   &
         &                     +1.0d0                             &
         &                     +fractionMassBound                 &
         &                    )                  **self%muRadius
    if (fractionMassBound >= fractionMassTransition) then
       ratioRadiusMaximumRadiusScale=self%ratioRadiusMaximumRadiusScaleUnstripped
    else
       ratioRadiusMaximumRadiusScale=self%ratioRadiusMaximumRadiusScaleStripped
    end if
    radiusScale=+                  fractionRadiusMaximum                     &
         &      *self             %ratioRadiusMaximumRadiusScaleUnstripped   &
         &      /                  ratioRadiusMaximumRadiusScale             &
         &      *darkMatterProfile%scale                                  ()
    return
  end function penarrubia2010ScaleRadius

  double precision function penarrubia2010Normalization(self,node)
    !!{
    Compute the mass normalization of the {\normalfont \ttfamily penarrubia2010} dark matter halo profile.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class(darkMatterProfileDMOPenarrubia2010       ), intent(inout) :: self
    type (treeNode                                 ), intent(inout) :: node
    class(nodeComponentBasic                       ), pointer       :: basic
    class           (nodeComponentSatellite        ), pointer       :: satellite
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    double precision                                                :: fractionMassBound     , fractionVelocityMaximum, &
         &                                                             radiusScale           , radiusVirial           , &
         &                                                             velocityRotationFactor
    
    basic                   =>  node     %basic            ()
    satellite               =>  node     %satellite        ()
    darkMatterProfile       =>  node     %darkMatterProfile()
    fractionMassBound       =  +satellite%boundMass        () &
         &                     /basic    %mass             ()
    fractionVelocityMaximum =  +2.0d0              **self%muVelocity  &
         &                     *  fractionMassBound**self%etaVelocity &
         &                     /(                                     &
         &                       +1.0d0                               &
         &                       +fractionMassBound                   &
         &                      )                  **self%muVelocity
    radiusScale             =darkMatterProfile                     %scale       (    )
    radiusVirial            =self             %darkMatterHaloScale_%radiusVirial(node)
    if (fractionMassBound >= fractionMassTransition) then
       velocityRotationFactor=+self%massDistributionUnstripped%        rotationCurve       (radius=radiusVirial/self%scaleRadius(node)) &
            &                 /self%massDistributionUnstripped%velocityRotationCurveMaximum(                                          )
    else
       velocityRotationFactor=+self%massDistributionStripped  %        rotationCurve       (radius=radiusVirial/self%scaleRadius(node)) &
            &                 /self%massDistributionStripped  %velocityRotationCurveMaximum(                                          )
    end if
    penarrubia2010Normalization=+basic%mass()                                                                                    &
         &                      *fractionVelocityMaximum**2                                                                      &
         &                      *velocityRotationFactor **2                                                                      &
         &                      /(                                                                                               &
         &                        +self%massDistributionUnstripped%rotationCurve               (radius=radiusVirial/radiusScale) &
         &                        /self%massDistributionUnstripped%velocityRotationCurveMaximum(                               ) &
         &                       )**2
    return
  end function penarrubia2010Normalization
