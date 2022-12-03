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
An implementation of the hot halo mass distribution class which uses the model of \cite{ricotti_feedback_2000}.
!!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <hotHaloMassDistribution name="hotHaloMassDistributionRicotti2000">
   <description>
    A hot halo mass distribution class which adopts a spherically symmetric $\beta$-profile density profile for the hot halo
    with parameters selected using the fitting function of \cite{ricotti_feedback_2000} who found these by solving for
    hydrostatic equilibrium in an \gls{nfw} density profile. Specifically, $\beta = 0.9 b$ where
    \begin{equation}
     b = {2 c \over 9 \Gamma} \left[ \log(1+c) - {c\over 1+c}\right]^{-1},
    \end{equation}
    $c$ is the concentration parameter of the dark matter halo, and $\Gamma$ is the ratio of virial and gas temperatures, which
    is assumed to be unity. The core radius is $r_\mathrm{core} = 0.22 r_\mathrm{s}$ where $r_\mathrm{s}$ is the scale radius
    of the dark matter profile. The profile is normalized such that the current mass in the hot gas profile is contained within
    the outer radius of the hot halo, $r_\mathrm{hot, outer}$.
   </description>
  </hotHaloMassDistribution>
  !!]
  type, extends(hotHaloMassDistributionBetaProfile) :: hotHaloMassDistributionRicotti2000
     !!{
     An implementation of the hot halo mass distribution class which uses the model of \cite{ricotti_feedback_2000}.
     !!}
     private
     class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     class(darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
   contains
     final     ::               ricotti2000Destructor
     procedure :: initialize => ricotti2000Initialize
  end type hotHaloMassDistributionRicotti2000

  interface hotHaloMassDistributionRicotti2000
     !!{
     Constructors for the {\normalfont \ttfamily ricotti2000} hot halo mass distribution class.
     !!}
     module procedure ricotti2000ConstructorParameters
     module procedure ricotti2000ConstructorInternal
  end interface hotHaloMassDistributionRicotti2000

contains

  function ricotti2000ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily ricotti2000} hot halo mass distribution class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloMassDistributionRicotti2000)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_
    class(darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=hotHaloMassDistributionRicotti2000(darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function ricotti2000ConstructorParameters

  function ricotti2000ConstructorInternal(darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily ricotti2000} hot halo mass distribution class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List                   , Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent, defaultHotHaloComponent
    implicit none
    type   (hotHaloMassDistributionRicotti2000)                        :: self
    class  (darkMatterProfileDMOClass         ), intent(in   ), target :: darkMatterProfileDMO_
    class  (darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    logical                                    , save                  :: initialized          =.false.
    !![
    <constructorAssign variables="*darkMatterProfileDMO_, *darkMatterHaloScale_"/>
    !!]

    if (.not.initialized) then
       !$omp critical(ricotti2000Initialized)
       if (.not.initialized) then
          ! Check that required properties are gettable.
          if     (                                                                                                      &
               &  .not.(                                                                                                &
               &         defaultHotHaloComponent          %       massIsGettable()                                      &
               &        .and.                                                                                           &
               &         defaultHotHaloComponent          %outerRadiusIsGettable()                                      &
               &       )                                                                                                &
               & ) call Error_Report                                                                                    &
               & (                                                                                                      &
               &  'This method requires that the "mass" property of the hot halo is gettable.'//                        &
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
          ! Record that implementation is now initialized.
          initialized=.true.
       end if
       !$omp end critical(ricotti2000Initialized)
    end if
    ! Set the value of β. This is arbitrary as it will be computed as needed, but avoids compiler complaints that this
    ! constructor is not initialized
    self%beta=-1.0d0
    return
  end function ricotti2000ConstructorInternal

  subroutine ricotti2000Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily ricotti2000} hot halo mass distribution class.
    !!}
    implicit none
    type(hotHaloMassDistributionRicotti2000), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine ricotti2000Destructor

  subroutine ricotti2000Initialize(self,node)
    !!{
    Initialize the {\normalfont \ttfamily ricotti2000} hot halo density profile for the given {\normalfont \ttfamily
    node}. Parameterizations of $\beta$ and core radius are taken from section 2.1 of \cite{ricotti_feedback_2000}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloMassDistributionRicotti2000    ), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    class           (nodeComponentHotHalo                  ), pointer       :: hotHalo
    class           (nodeComponentDarkMatterProfile        ), pointer       :: darkMatterProfile
    double precision                                        , parameter     :: virialToGasTemperatureRatio=1.0d0
    double precision                                                        :: mass                             , radiusOuter  , &
         &                                                                     radiusScale                      , radiusVirial , &
         &                                                                     radiusCore                       , concentration, &
         &                                                                     b                                , beta

    ! Compute parameters of the profile.
    hotHalo           =>  node                        %hotHalo          (    )
    darkMatterProfile =>  node                        %darkMatterProfile(    )
    radiusOuter       =   hotHalo                     %outerRadius      (    )
    mass              =   hotHalo                     %mass             (    )
    radiusScale       =           darkMatterProfile   %scale            (    )
    radiusVirial      =   self   %darkMatterHaloScale_%radiusVirial     (node)
    concentration     =  +radiusVirial              &
         &               /radiusScale
    b                 =  +(                             &
         &                 +2.0d0                       &
         &                 *concentration               &
         &                 /9.0d0                       &
         &                 /virialToGasTemperatureRatio &
         &                )                             &
         &               /(                             &
         &                 +log(1.0d0+concentration)    &
         &                 -          concentration     &
         &                 /(                           &
         &                   +1.0d0                     &
         &                   +concentration             &
         &                  )                           &
         &                )
    beta              =  +0.90d0      &
         &               *b
    radiusCore        =  +0.22d0      &
         &               *radiusScale
    ! Construct the mass distribution.
    if (radiusOuter <= 0.0d0) then
       ! If outer radius is non-positive, set mass to zero and outer radius to an arbitrary value.
       mass=0.0d0
       radiusOuter=1.0d0
    end if
     self%distribution=massDistributionBetaProfile(                         &
         &                                         beta       =beta       , &
         &                                         coreRadius =radiusCore , &
         &                                         mass       =mass       , &
         &                                         outerRadius=radiusOuter  &
         &                                        )
    return
  end subroutine ricotti2000Initialize
