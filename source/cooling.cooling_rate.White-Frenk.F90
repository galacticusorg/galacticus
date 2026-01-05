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
  Implementation of a cooling rate class for the \cite{white_galaxy_1991} cooling rate calculation.
  !!}

  use :: Cooling_Infall_Radii   , only : coolingInfallRadiusClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <coolingRate name="coolingRateWhiteFrenk1991">
   <description>
    A cooling rate class that uses the expression given by \cite{white_galaxy_1991}, namely
    \begin{equation}
    \dot{M}_\mathrm{cool} = \left\{ \begin{array}{ll} 4 \pi r_\mathrm{infall}^2 \rho(r_\mathrm{infall})
    \dot{r}_\mathrm{infall}&amp; \hbox{ if } r_\mathrm{infall} &lt; r_\mathrm{hot, outer} \\ M_\mathrm{hot}/\tau_\mathrm{halo,
    dynamical} &amp; \hbox{ if } r_\mathrm{infall} \ge r_\mathrm{hot, outer}, \end{array} \right. ,
    \end{equation}
    where $r_\mathrm{infall}$ is the infall radius (see \refPhysics{coolingInfallRadius}) in the hot halo and $\rho(r)$ is the
    density profile of the hot halo.
   </description>
  </coolingRate>
  !!]
  type, extends(coolingRateClass) :: coolingRateWhiteFrenk1991
     !!{
     Implementation of cooling rate class for the \cite{white_galaxy_1991} cooling rate calculation.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     class           (coolingInfallRadiusClass), pointer :: coolingInfallRadius_ => null()
     double precision                                    :: velocityCutOff
   contains
     final     ::         whiteFrenk1991Destructor
     procedure :: rate => whiteFrenk1991Rate
  end type coolingRateWhiteFrenk1991

  interface coolingRateWhiteFrenk1991
     !!{
     Constructors for the \cite{white_galaxy_1991} cooling rate class.
     !!}
     module procedure whiteFrenk1991ConstructorParameters
     module procedure whiteFrenk1991ConstructorInternal
  end interface coolingRateWhiteFrenk1991

contains

  function whiteFrenk1991ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{white_galaxy_1991} cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coolingRateWhiteFrenk1991)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass ), pointer       :: darkMatterHaloScale_
    class           (coolingInfallRadiusClass ), pointer       :: coolingInfallRadius_
    double precision                                           :: velocityCutOff

    !![
    <inputParameter>
      <name>velocityCutOff</name>
      <source>parameters</source>
      <defaultValue>1.0d4</defaultValue>
      <description>The halo virial velocity (in km/s) above which cooling rates are forced to zero in the \cite{white_galaxy_1991} cooling rate model.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="coolingInfallRadius" name="coolingInfallRadius_" source="parameters"/>
    !!]
    self=coolingRateWhiteFrenk1991(velocityCutOff,darkMatterHaloScale_,coolingInfallRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="coolingInfallRadius_"/>
    !!]
    return
  end function whiteFrenk1991ConstructorParameters

  function whiteFrenk1991ConstructorInternal(velocityCutOff,darkMatterHaloScale_,coolingInfallRadius_) result(self)
    !!{
    Internal constructor for the \cite{white_galaxy_1991} cooling rate class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List          , Error_Report
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none
    type            (coolingRateWhiteFrenk1991)                        :: self
    double precision                           , intent(in   )         :: velocityCutOff
    class           (darkMatterHaloScaleClass ), intent(in   ), target :: darkMatterHaloScale_
    class           (coolingInfallRadiusClass ), intent(in   ), target :: coolingInfallRadius_
    !![
    <constructorAssign variables="velocityCutOff, *darkMatterHaloScale_, *coolingInfallRadius_"/>
    !!]

    ! Check that the properties we need are gettable.
    if     (                                                                                                   &
         &  .not.(                                                                                             &
         &         defaultHotHaloComponent%       massIsGettable()                                             &
         &        .and.                                                                                        &
         &         defaultHotHaloComponent%outerRadiusIsGettable()                                             &
         &       )                                                                                             &
         & ) call Error_Report                                                                                 &
         &        (                                                                                            &
         &         'mass and outerRadius properties of hot halo component must be gettable.'//                 &
         &         Component_List(                                                                             &
         &                        'hotHalo'                                                                 ,  &
         &                         defaultHotHaloComponent%       massAttributeMatch(requireGettable=.true.)   &
         &                        .intersection.                                                               &
         &                         defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.)   &
         &                       )                                                                          // &
         &         {introspection:location}                                                                    &
         &        )
    return
  end function whiteFrenk1991ConstructorInternal

  subroutine whiteFrenk1991Destructor(self)
    !!{
    Destructor for the \cite{white_galaxy_1991} cooling rate class.
    !!}
    implicit none
    type(coolingRateWhiteFrenk1991), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%coolingInfallRadius_"/>
    !!]
    return
  end subroutine whiteFrenk1991Destructor

  double precision function whiteFrenk1991Rate(self,node)
    !!{
    Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for the \cite{white_galaxy_1991} cooling rate
    model.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo , treeNode
    use :: Numerical_Constants_Math  , only : Pi
    use :: Galactic_Structure_Options, only : componentTypeHotHalo , massTypeGaseous
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Coordinates               , only : coordinateSpherical  , assignment(=)
    implicit none
    class           (coolingRateWhiteFrenk1991), intent(inout) :: self
    type            (treeNode                 ), intent(inout) :: node
    class           (nodeComponentHotHalo     ), pointer       :: hotHalo
    class           (massDistributionClass    ), pointer       :: massDistribution_
    type            (coordinateSpherical      )                :: coordinates
    double precision                                           :: densityCooling        , radiusInfall  , &
         &                                                        radiusOuter           , velocityVirial, &
         &                                                        radiusInfallGrowthRate

    ! Get the virial velocity.
    velocityVirial=self%darkMatterHaloScale_%velocityVirial(node)
    ! Return zero cooling rate if virial velocity exceeds critical value.
    if (velocityVirial > self%velocityCutOff) then
       whiteFrenk1991Rate=0.0d0
       return
    end if
    ! Get the outer radius of the hot halo.
    hotHalo     => node   %hotHalo    ()
    radiusOuter =  hotHalo%outerRadius()
    ! Get the infall radius.
    radiusInfall=self%coolingInfallRadius_%radius(node)
    if (radiusInfall >= radiusOuter) then
       ! Infall radius exceeds the outer radius. Limit infall to the dynamical timescale.
       whiteFrenk1991Rate=+hotHalo                     %mass              (    ) &
            &             /self   %darkMatterHaloScale_%timescaleDynamical(node)
    else
       ! Find the density at the cooling radius.  
       coordinates              =  [radiusInfall,0.0d0,0.0d0]
       massDistribution_        => node                                  %massDistribution  (componentTypeHotHalo,massTypeGaseous)
       densityCooling           =  massDistribution_                     %density           (coordinates                         )
       !![
       <objectDestructor name="massDistribution_"/>
       !!]       
       ! Find infall radius growth rate.
       radiusInfallGrowthRate   =  self             %coolingInfallRadius_%radiusIncreaseRate(node                                )
       ! Compute the infall rate.
       whiteFrenk1991Rate       =  +4.0d0                     &
            &                      *Pi                        &
            &                      *radiusInfall          **2 &
            &                      *densityCooling            &
            &                      *radiusInfallGrowthRate
    end if
    return
  end function whiteFrenk1991Rate
