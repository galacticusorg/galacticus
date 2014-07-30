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

!% An implementation of the hot halo mass distribution class which uses the model of \cite{ricotti_feedback_2000}.

  !# <hotHaloMassDistribution name="hotHaloMassDistributionRicotti2000">
  !#  <description>Provides an implementation of the hot halo mass distribution class which uses the model of \cite{ricotti_feedback_2000}.</description>
  !# </hotHaloMassDistribution>
  type, extends(hotHaloMassDistributionBetaProfile) :: hotHaloMassDistributionRicotti2000
     !% An implementation of the hot halo mass distribution class which uses the model of \cite{ricotti_feedback_2000}.
     private
   contains
     procedure :: initialize => ricotti2000Initialize
  end type hotHaloMassDistributionRicotti2000

  interface hotHaloMassDistributionRicotti2000
     !% Constructors for the {\tt ricotti2000} hot halo mass distribution class.
     module procedure ricotti2000DefaultConstructor
  end interface hotHaloMassDistributionRicotti2000

  logical :: ricotti2000Initialized=.false.

contains

  function ricotti2000DefaultConstructor()
    !% Default constructor for the ricotti2000 hot halo mass distribution class.
    use Galacticus_Error
    use Array_Utilities
    implicit none
    type(hotHaloMassDistributionRicotti2000) :: ricotti2000DefaultConstructor

    if (.not.ricotti2000Initialized) then
       !$omp critical(ricotti2000Initialized)
       if (.not.ricotti2000Initialized) then
          ! Check that required properties are gettable.
          if     (                                                                                                                &
               &  .not.(                                                                                                          &
               &         defaultHotHaloComponent          %       massIsGettable()                                                &
               &        .and.                                                                                                     &
               &         defaultHotHaloComponent          %outerRadiusIsGettable()                                                &
               &       )                                                                                                          &
               & ) call Galacticus_Error_Report                                                                                   &
               & (                                                                                                                &
               &  'ricotti2000DefaultConstructor'                                                                               , &
               &  'This method requires that the "mass" property of the hot halo is gettable.'//                                  &
               &  Galacticus_Component_List(                                                                                      &
               &                            'hotHalo'                                                                           , &
               &                             defaultHotHaloComponent          %       massAttributeMatch(requireGettable=.true.)  &
               &                            .intersection.                                                                        &
               &                             defaultHotHaloComponent          %outerRadiusAttributeMatch(requireGettable=.true.)  &
               &                           )                                                                                      &
               & )
          if     (                                                                                                                &
               &  .not.(                                                                                                          &
               &         defaultDarkMatterProfileComponent%      scaleIsGettable()                                                &
               &       )                                                                                                          &
               & ) call Galacticus_Error_Report                                                                                   &
               & (                                                                                                                &
               &  'ricotti2000DefaultConstructor'                                                                               , &
               &  'This method requires that the "scale" property of the dark matter profile is gettable.'//                      &
               &  Galacticus_Component_List(                                                                                      &
               &                            'darkMatterProfile'                                                                 , &
               &                             defaultDarkMatterProfileComponent%      scaleAttributeMatch(requireGettable=.true.)  &
               &                           )                                                                                      &
               & )
          ! Record that implementation is now initialized.
          ricotti2000Initialized=.true.
       end if
       !$omp end critical(ricotti2000Initialized)
    end if
    return
  end function ricotti2000DefaultConstructor

  subroutine ricotti2000Initialize(self,node)
    !% Initialize the {\tt ricotti2000} hot halo density profile for the given {\tt node}. Parameterizations of $\beta$ and core
    !% radius are taken from section 2.1 of \cite{ricotti_feedback_2000}.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (hotHaloMassDistributionRicotti2000    ), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo
    class           (nodeComponentDarkMatterProfile        )               , pointer :: darkMatterProfile
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    double precision                                        , parameter              :: virialToGasTemperatureRatio=1.0d0
    double precision                                                                 :: mass                             , radiusOuter  , &
         &                                                                              radiusScale                      , radiusVirial , &
         &                                                                              radiusCore                       , concentration, &
         &                                                                              b                                , beta
    
    ! Compute parameters of the profile.
    darkMatterHaloScale_ => darkMatterHaloScale                   (    )
    hotHalo              => node                %hotHalo          (    )
    darkMatterProfile    => node                %darkMatterProfile(    )
    radiusOuter          =  hotHalo             %outerRadius      (    )
    mass                 =  hotHalo             %mass             (    )
    radiusScale          =  darkMatterProfile   %scale            (    )
    radiusVirial         =  darkMatterHaloScale_%virialRadius     (node)
    concentration        =  radiusVirial/radiusScale
    b                    =   (                             &
         &                     2.0d0                       &
         &                    *concentration               &
         &                    /9.0d0                       &
         &                    /virialToGasTemperatureRatio &
         &                   )                             &
         &                  /(                             &
         &                     log(1.0d0+concentration)    &
         &                    -concentration               &
         &                    /(                           &
         &                       1.0d0                     &
         &                      +concentration             &
         &                     )                           &
         &                   )
    beta                 =  0.90d0*b
    radiusCore           =  0.22d0*radiusScale
    ! Construct the mass distribution.
    if (radiusOuter <= 0.0d0) then
       ! If outer radius is non-positive, set mass to zero and outer radius to an arbitrary value.
       mass=0.0d0
       radiusOuter=1.0d0
    end if
    call self%                                  &
         & distribution%                        &
         &  initialize(                         &
         &             beta       =beta       , &
         &             coreRadius =radiusCore , &
         &             mass       =mass       , &
         &             outerRadius=radiusOuter  &
         &            )
    return
  end subroutine ricotti2000Initialize
