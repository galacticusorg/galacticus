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
  Implementation of a cooling rate class for the \cite{cole_hierarchical_2000} cooling rate calculation.
  !!}

  use :: Cooling_Infall_Radii, only : coolingInfallRadiusClass

  !![
  <coolingRate name="coolingRateCole2000">
   <description>
    A cooling rate class that uses the algorithm of \cite{cole_hierarchical_2000}. The cooling rate is given by
    \begin{equation}
    \dot{M}_\mathrm{cool} = \left\{ \begin{array}{ll} 4 \pi r_\mathrm{infall}^2 \rho(r_\mathrm{infall}) \dot{r}_\mathrm{infall}
    &amp; \hbox{ if } r_\mathrm{infall} &lt; r_\mathrm{hot, outer} \\ 0 &amp; \hbox{ if } r_\mathrm{infall} \ge r_\mathrm{hot,
    outer}, \end{array} \right.
    \end{equation}
    where $\rho(r)$ is the density profile of the hot halo, and $r_\mathrm{infall}$ is the infall radius (see
    \refPhysics{coolingInfallRadius}).
   </description>
  </coolingRate>
  !!]
  type, extends(coolingRateClass) :: coolingRateCole2000
     !!{
     Implementation of cooling rate class for the \cite{cole_hierarchical_2000} cooling rate calculation.
     !!}
     private
     class(coolingInfallRadiusClass), pointer :: coolingInfallRadius_ => null()
   contains
     final     ::         cole2000Destructor
     procedure :: rate => cole2000Rate
  end type coolingRateCole2000

  interface coolingRateCole2000
     !!{
     Constructors for the \cite{cole_hierarchical_2000} cooling rate class.
     !!}
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface coolingRateCole2000

contains

  function cole2000ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{cole_hierarchical_2000} cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingRateCole2000     )                :: self
    type (inputParameters         ), intent(inout) :: parameters
    class(coolingInfallRadiusClass), pointer       :: coolingInfallRadius_

    !![
    <objectBuilder class="coolingInfallRadius" name="coolingInfallRadius_" source="parameters"/>
    !!]
    self=coolingRateCole2000(coolingInfallRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingInfallRadius_"/>
    !!]
    return
  end function cole2000ConstructorParameters

  function cole2000ConstructorInternal(coolingInfallRadius_) result(self)
    !!{
    Internal constructor for the \cite{cole_hierarchical_2000} cooling rate class.
    !!}
    implicit none
    type (coolingRateCole2000     )                        :: self
    class(coolingInfallRadiusClass), intent(in   ), target :: coolingInfallRadius_
    !![
    <constructorAssign variables="*coolingInfallRadius_"/>
    !!]

    return
  end function cole2000ConstructorInternal

  subroutine cole2000Destructor(self)
    !!{
    Destructor for the \cite{cole_hierarchical_2000} cooling rate class.
    !!}
    implicit none
    type(coolingRateCole2000), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingInfallRadius_"/>
    !!]
    return
  end subroutine cole2000Destructor

  double precision function cole2000Rate(self,node)
    !!{
    Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for the \cite{white_galaxy_1991} cooling rate
    model.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic   , nodeComponentHotHalo, treeNode
    use :: Numerical_Constants_Math  , only : Pi
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Coordinates               , only : coordinateSpherical  , assignment(=)
    use :: Galactic_Structure_Options, only : componentTypeHotHalo , massTypeGaseous
    implicit none
    class           (coolingRateCole2000  ), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    class           (nodeComponentBasic   ), pointer       :: basicFormation
    class           (nodeComponentHotHalo ), pointer       :: hotHaloFormation
    class           (massDistributionClass), pointer       :: massDistribution_
    type            (coordinateSpherical  )                :: coordinates
    double precision                                       :: densityInfall            , radiusInfall, &
         &                                                    radiusInfallGrowthRate   , radiusOuter
    !$GLC attributes unused :: self

    ! Get formation node components.
    basicFormation   => node%formationNode%basic  ()
    hotHaloFormation => node%formationNode%hotHalo()
    ! Check for empty halos.
    if (hotHaloFormation%mass() <= 0.0d0) then
       cole2000Rate=0.0d0
       return
    end if
    ! Get the outer radius of the hot halo.
    radiusOuter =hotHaloFormation%outerRadius()
    ! Get the infall radius.
    radiusInfall=self%coolingInfallRadius_%radius(node%formationNode)
    if (radiusInfall >= radiusOuter) then
       ! Infall radius exceeds the outer radius - zero infall rate.
       cole2000Rate=0.0d0
    else
       ! Find the density at the infall radius.
       coordinates              =  [radiusInfall,0.0d0,0.0d0]
       massDistribution_        => node             %formationNode       %massDistribution  (componentTypeHotHalo,massTypeGaseous)
       densityInfall            =  massDistribution_                     %density           (coordinates                         )
       !![
       <objectDestructor name="massDistribution_"/>
       !!]          
       ! Find infall radius growth rate.
       radiusInfallGrowthRate   =  self             %coolingInfallRadius_%radiusIncreaseRate(node%formationNode                  )
       ! Compute the infall rate.
       cole2000Rate             = +4.0d0                     &
            &                     *Pi                        &
            &                     *radiusInfall          **2 &
            &                     *densityInfall             &
            &                     *radiusInfallGrowthRate
    end if
    return
  end function cole2000Rate
