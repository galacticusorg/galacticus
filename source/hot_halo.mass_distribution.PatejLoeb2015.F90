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
An implementation of the hot halo mass distribution class which uses the model of \cite{patej_simple_2015}.
!!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <hotHaloMassDistribution name="hotHaloMassDistributionPatejLoeb2015">
   <description>Provides an implementation of the hot halo mass distribution class which uses the model of \cite{patej_simple_2015}.</description>
  </hotHaloMassDistribution>
  !!]
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionPatejLoeb2015
     !!{
     An implementation of the hot halo mass distribution class which uses the model of \cite{patej_simple_2015}.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     double precision                                     :: gamma                          , radiusShock
   contains
     final     ::                          patejLoeb2015Destructor
     procedure :: density               => patejLoeb2015Density
     procedure :: densityLogSlope       => patejLoeb2015DensityLogSlope
     procedure :: enclosedMass          => patejLoeb2015EnclosedMass
     procedure :: radialMoment          => patejLoeb2015RadialMoment
     procedure :: rotationNormalization => patejLoeb2015RotationNormalization
  end type hotHaloMassDistributionPatejLoeb2015

  interface hotHaloMassDistributionPatejLoeb2015
     !!{
     Constructors for the {\normalfont \ttfamily patejLoeb2015} hot halo mass distribution class.
     !!}
     module procedure patejLoeb2015ConstructorParameters
     module procedure patejLoeb2015ConstructorInternal
  end interface hotHaloMassDistributionPatejLoeb2015

contains

  function patejLoeb2015ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily patejLoeb2015} hot halo mass distribution class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hotHaloMassDistributionPatejLoeb2015)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass           ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), pointer       :: darkMatterHaloScale_
    double precision                                                      :: gamma               , radiusShock

    !![
    <inputParameter>
      <name>gamma</name>
      <defaultValue>1.15d0</defaultValue>
      <description>The parameter $\Gamma$ in the \cite{patej_simple_2015} hot halo gas mass distribution model.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusShock</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The shock radius, $s$, (in units of the halo virial radius) in the \cite{patej_simple_2015} hot halo gas mass distribution model.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=hotHaloMassDistributionPatejLoeb2015(gamma,radiusShock,darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function patejLoeb2015ConstructorParameters

  function patejLoeb2015ConstructorInternal(gamma,radiusShock,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily patejLoeb2015} hot halo mass distribution class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List                   , Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent, defaultHotHaloComponent
    implicit none
    type            (hotHaloMassDistributionPatejLoeb2015)                        :: self
    double precision                                      , intent(in   )         :: gamma                         , radiusShock
    class           (darkMatterProfileDMOClass           ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), intent(in   ), target :: darkMatterHaloScale_
    logical                                               , save                  :: initialized           =.false.
    !![
    <constructorAssign variables="gamma, radiusShock, *darkMatterProfileDMO_, *darkMatterHaloScale_"/>
    !!]

    ! Check that required properties are gettable.
    if (.not.initialized) then
       !$omp critical(patejLoeb2015Initialized)
       if (.not.initialized) then
          if     (                                                                                                      &
               &  .not.(                                                                                                &
               &         defaultHotHaloComponent          %       massIsGettable()                                      &
               &        .and.                                                                                           &
               &         defaultHotHaloComponent          %outerRadiusIsGettable()                                      &
               &       )                                                                                                &
               & ) call Error_Report                                                                                    &
               & (                                                                                                      &
               &  'This method requires that the "mass" and "outerRadius" properties of the hot halo are gettable.'//   &
               &  Component_List(                                                                                       &
               &                 'hotHalo'                                                                           ,  &
               &                  defaultHotHaloComponent          %       massAttributeMatch(requireGettable=.true.)   &
               &                 .intersection.                                                                         &
               &                  defaultHotHaloComponent          %outerRadiusAttributeMatch(requireGettable=.true.)   &
               &                )                                                                                    // &
               &  {introspection:location}                                                                              &
               & )
          if     (                                                                                                      &
               &  .not.(                                                                                                &
               &         defaultDarkMatterProfileComponent%      scaleIsGettable()                                      &
               &       )                                                                                                &
               & ) call Error_Report                                                                                    &
               & (                                                                                                      &
               &  'This method requires that the "scale" property of the dark matter profile is gettable.'//            &
               &  Component_List(                                                                                       &
               &                 'darkMatterProfile'                                                                 ,  &
               &                  defaultDarkMatterProfileComponent%      scaleAttributeMatch(requireGettable=.true.)   &
               &                )                                                                                    // &
               &  {introspection:location}                                                                              &
               & )
          initialized=.true.
       end if
       !$omp end critical(patejLoeb2015Initialized)
    end if
    return
  end function patejLoeb2015ConstructorInternal

  subroutine patejLoeb2015Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily patejLoeb2015} hot halo mass distribution class.
    !!}
    implicit none
    type(hotHaloMassDistributionPatejLoeb2015), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine patejLoeb2015Destructor

  double precision function patejLoeb2015Density(self,node,radius)
    !!{
    Return the density in a {\normalfont \ttfamily patejLoeb2015} hot halo mass distribution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloMassDistributionPatejLoeb2015), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    class           (nodeComponentHotHalo                ), pointer       :: hotHalo
    class           (nodeComponentDarkMatterProfile      ), pointer       :: darkMatterHaloProfile
    double precision                                                      :: radiusShock          , densityNormalization, &
         &                                                                   radiusOuter          , radiusDarkMatter

    ! Get the hot halo and dark matter profile components.
    hotHalo               => node%hotHalo          ()
    darkMatterHaloProfile => node%darkMatterProfile()
    ! Find the shock and outer radii.
    radiusShock         =+self                        %radiusShock                         &
         &               *self   %darkMatterHaloScale_%radiusVirial(node                 )
    radiusOuter         =+hotHalo                     %outerRadius (                     )
    ! Find the density normalization.
    radiusDarkMatter    =+radiusShock                                                         &
         &               *(radiusOuter/radiusShock)**self%gamma
    densityNormalization=+hotHalo                        %mass        (                     ) &
         &               /self   %darkMatterProfileDMO_  %enclosedMass(node,radiusDarkMatter)
    ! Compute the density.
    patejLoeb2015Density=+self%gamma                                                         &
         &               *densityNormalization                                               &
         &               *(                                                                  &
         &                 +radius                                                           &
         &                 /radiusShock                                                      &
         &                )**(3.0d0*self%gamma-3.0d0)                                        &
         &               *self%darkMatterProfileDMO_%density(                                &
         &                                                   node                          , &
         &                                                   +radiusShock                    &
         &                                                   *(                              &
         &                                                     +radius                       &
         &                                                     /radiusShock                  &
         &                                                   )**self%gamma                   &
         &                                                  )
    return
  end function patejLoeb2015Density

  double precision function patejLoeb2015DensityLogSlope(self,node,radius)
    !!{
    Return the density in a {\normalfont \ttfamily patejLoeb2015} hot halo mass distribution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (hotHaloMassDistributionPatejLoeb2015), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    class           (nodeComponentDarkMatterProfile      ), pointer       :: darkMatterHaloProfile
    double precision                                                      :: radiusShock

    ! Get the dark matter profile component.
    darkMatterHaloProfile => node%darkMatterProfile()
    ! Find the shock radius.
    radiusShock         =+self                     %radiusShock                    &
         &               *self%darkMatterHaloScale_%radiusVirial(node            )
    ! Compute the log slope of density.
    patejLoeb2015DensityLogSlope=+3.0d0                                                       &
         &                       *(self%gamma-1.0d0)                                          &
         &                       + self%gamma                                                 &
         &                       * self%darkMatterProfileDMO_%densityLogSlope(                &
         &                                                                    node          , &
         &                                                                    +radiusShock    &
         &                                                                    *(              &
         &                                                                      +radius       &
         &                                                                      /radiusShock  &
         &                                                                     )**self%gamma  &
         &                                                                   )
    return
  end function patejLoeb2015DensityLogSlope

  double precision function patejLoeb2015EnclosedMass(self,node,radius)
    !!{
    Return the enclosed mass in a {\normalfont \ttfamily patejLoeb2015} hot halo mass distribution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloMassDistributionPatejLoeb2015), intent(inout)          :: self
    type            (treeNode                            ), intent(inout), target  :: node
    double precision                                      , intent(in   )          :: radius
    class           (nodeComponentHotHalo                )               , pointer :: hotHalo
    class           (nodeComponentDarkMatterProfile      )               , pointer :: darkMatterHaloProfile
    double precision                                                               :: radiusShock          , densityNormalization, &
         &                                                                            radiusOuter          , radiusScale         , &
         &                                                                            radiusDarkMatter

    ! Get the hot halo and dark matter profile components.
    hotHalo               => node%hotHalo          ()
    darkMatterHaloProfile => node%darkMatterProfile()
    ! Find the shock, outer, and scale radii.
    radiusShock              =+self                     %radiusShock                         &
         &                    *self%darkMatterHaloScale_%radiusVirial(node                 )
    radiusOuter              =     hotHalo              %outerRadius (                     )
    radiusScale              =     darkMatterHaloProfile%scale       (                     )
    ! Find the density normalization.
    radiusDarkMatter         =+radiusShock                                                   &
         &                    *(radiusOuter/radiusShock)**self%gamma
    densityNormalization     =+     hotHalo                %mass        (                     ) &
         &                    /self%darkMatterProfileDMO_  %enclosedMass(node,radiusDarkMatter)
    ! Compute the corresponding radius in the dark matter halo.
    radiusDarkMatter         =+radiusShock                                                   &
         &                    *(radius     /radiusShock)**self%gamma
    ! Compute the enclosed mass (eqn. 4 of Patej & Loeb 2015).
    patejLoeb2015EnclosedMass=+densityNormalization                                      &
         &                    *self%darkMatterProfileDMO_%enclosedMass(                  &
         &                                                             node            , &
         &                                                             radiusDarkMatter  &
         &                                                           )
    return
  end function patejLoeb2015EnclosedMass

  double precision function patejLoeb2015RadialMoment(self,node,moment,radius)
    !!{
    Compute a radial moment in a {\normalfont \ttfamily patejLoeb2015} hot halo mass distribution.
    For this profile we have:
    \begin{equation}
    \rho_\mathrm{g}(r) = f \Gamma (r/s)^{3 \Gamma - 3} \rho_\mathrm{DM}(s[r/s]^\Gamma).
    \end{equation}
    Defining $R=s[r/s]^\Gamma$, such that $r/s = (R/s)^{1/\Gamma}$, and $\mathrm{d}r = \Gamma^{-1} (R/s)^{1/\Gamma-1} \mathrm{d}R$, then
    \begin{equation}
    \int r^m \rho_\mathrm{g}(r) \mathrm{d}r = f s^{(m-2)(\Gamma-1)/\Gamma} \int R^{(2\Gamma-2+m)/\Gamma} \rho_\mathrm{DM}(R) \mathrm{d}R,
    \end{equation}
    or
    \begin{equation}
    \mathcal{R}_\mathrm{g}(r;m) = f s^{(m-2)(\Gamma-1)/\Gamma} \mathcal{R}_\mathrm{DM}(r;(2\Gamma-2+m)/\Gamma),
    \end{equation}
    where $\mathcal{R}(r;m)$ is the $m^\mathrm{th}$ radial moment of the density profile.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloMassDistributionPatejLoeb2015), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: moment               , radius
    class           (nodeComponentHotHalo                ), pointer       :: hotHalo
    class           (nodeComponentDarkMatterProfile      ), pointer       :: darkMatterHaloProfile
    double precision                                                      :: radiusShock          , densityNormalization, &
         &                                                                   radiusOuter          , radiusScale         , &
         &                                                                   radiusDarkMatter

   ! Get the hot halo and dark matter profile components.
    hotHalo               => node%hotHalo          ()
    darkMatterHaloProfile => node%darkMatterProfile()
    ! Find the shock, outer, and scale radii.
    radiusShock         =+self                     %radiusShock                         &
         &               *self%darkMatterHaloScale_%radiusVirial(node                 )
    radiusOuter         =     hotHalo              %outerRadius (                     )
    radiusScale         =     darkMatterHaloProfile%scale       (                     )
    ! Find the density normalization.
    radiusDarkMatter    =+radiusShock                                                   &
         &               *(radiusOuter/radiusShock)**self%gamma
    densityNormalization=+     hotHalo                %mass        (                     ) &
         &               /self%darkMatterProfileDMO_  %enclosedMass(node,radiusDarkMatter)
    ! Compute the corresponding radius in the dark matter halo.
    radiusDarkMatter=+(                                   &
         &             +radiusShock                       &
         &             /radiusScale                       &
         &            )**(1.0d0-self%gamma)               &
         &            *(                                  &
         &             +min(radius,hotHalo%outerRadius()) &
         &             /radiusScale                       &
         &            )**       self%gamma
    ! Compute the radial moment.
    patejLoeb2015RadialMoment=+densityNormalization                                               &
         &                    *radiusShock**((self%gamma-1.0d0)*(moment-2.0d0)/self%gamma)        &
         &                    *self%darkMatterProfileDMO_%radialMoment(node,moment,radiusDarkMatter)
    return
  end function patejLoeb2015RadialMoment

  double precision function patejLoeb2015RotationNormalization(self,node)
    !!{
    Returns the relation between specific angular momentum and rotation velocity (assuming a
    rotation velocity that is constant in radius) for {\normalfont \ttfamily node}. Specifically, the
    normalization, $A$, returned is such that $V_\mathrm{rot} = A J/M$.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class(hotHaloMassDistributionPatejLoeb2015), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    class(nodeComponentHotHalo                ), pointer       :: hotHalo

    hotHalo                            =>  node%hotHalo     (                                )
    patejLoeb2015RotationNormalization =  +self%radialMoment(node,2.0d0,hotHalo%outerRadius()) &
         &                                /self%radialMoment(node,3.0d0,hotHalo%outerRadius())
    return
  end function patejLoeb2015RotationNormalization
