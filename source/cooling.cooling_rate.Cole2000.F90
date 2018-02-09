!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implementation of a cooling rate class for the \cite{cole_hierarchical_2000} cooling rate calculation.

  use Cooling_Infall_Radii
  use Hot_Halo_Mass_Distributions

  !# <coolingRate name="coolingRateCole2000" defaultThreadPrivate="yes">
  !#  <description>Computes the mass cooling rate in a hot gas halo utilizing the \cite{cole_hierarchical_2000} method. This is based on the
  !# properties of the halo at formation time, and gives a zero cooling rate when the cooling radius exceeds the virial radius.</description>
  !# </coolingRate>
  type, extends(coolingRateClass) :: coolingRateCole2000
     !% Implementation of cooling rate class for the \cite{cole_hierarchical_2000} cooling rate calculation.
     private
     class(coolingInfallRadiusClass    ), pointer :: coolingInfallRadius_
     class(hotHaloMassDistributionClass), pointer :: hotHaloMassDistribution_
   contains
     final     ::         cole2000Destructor
     procedure :: rate => cole2000Rate
  end type coolingRateCole2000

  interface coolingRateCole2000
     !% Constructors for the \cite{cole_hierarchical_2000} cooling rate class.
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface coolingRateCole2000

contains

  function cole2000ConstructorParameters(parameters) result(self)
    !% Constructor for the \cite{cole_hierarchical_2000} cooling rate class which builds the object from a parameter set.
    use Input_Parameters
    implicit none
    type (coolingRateCole2000         )                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(coolingInfallRadiusClass    ), pointer       :: coolingInfallRadius_
    class(hotHaloMassDistributionClass), pointer       :: hotHaloMassDistribution_
    
    !# <objectBuilder class="coolingInfallRadius"     name="coolingInfallRadius_"     source="parameters"/>
    !# <objectBuilder class="hotHaloMassDistribution" name="hotHaloMassDistribution_" source="parameters"/>
    self=coolingRateCole2000(coolingInfallRadius_,hotHaloMassDistribution_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function cole2000ConstructorParameters

  function cole2000ConstructorInternal(coolingInfallRadius_,hotHaloMassDistribution_) result(self)
    !% Internal constructor for the \cite{cole_hierarchical_2000} cooling rate class.
    implicit none
    type (coolingRateCole2000         )                        :: self
    class(coolingInfallRadiusClass    ), intent(in   ), target :: coolingInfallRadius_
    class(hotHaloMassDistributionClass), intent(in   ), target :: hotHaloMassDistribution_
    !# <constructorAssign variables="*coolingInfallRadius_, *hotHaloMassDistribution_"/>

    return
  end function cole2000ConstructorInternal

  subroutine cole2000Destructor(self)
    !% Destructor for the \cite{cole_hierarchical_2000} cooling rate class.
    implicit none
    type(coolingRateCole2000), intent(inout) :: self

    !# <objectDestructor name="self%coolingInfallRadius_"    />
    !# <objectDestructor name="self%hotHaloMassDistribution_"/>
    return
  end subroutine cole2000Destructor

  double precision function cole2000Rate(self,node)
    !% Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for the \cite{white_galaxy_1991} cooling rate
    !% model.
    use Numerical_Constants_Math
    implicit none
    class           (coolingRateCole2000 ), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node    
    class           (nodeComponentBasic  ), pointer       :: basicFormation
    class           (nodeComponentHotHalo), pointer       :: hotHaloFormation
    double precision                                      :: densityInfall            , radiusInfall, &
         &                                                   radiusInfallGrowthRate   , radiusOuter
    !GCC$ attributes unused :: self

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
       densityInfall            =  self%hotHaloMassDistribution_%density           (node%formationNode,radiusInfall)
       ! Find infall radius growth rate.
       radiusInfallGrowthRate   =  self%coolingInfallRadius_    %radiusIncreaseRate(node%formationNode             )
       ! Compute the infall rate.
       cole2000Rate             = +4.0d0                     &
            &                     *Pi                        &
            &                     *radiusInfall          **2 &
            &                     *densityInfall             &
            &                     *radiusInfallGrowthRate
    end if
    return
  end function cole2000Rate
