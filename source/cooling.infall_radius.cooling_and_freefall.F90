!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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
  
  !% Implementation of an infall radius calculation in which the infall radius is the smaller of the cooling and freefall radii.
  
  use Cooling_Radii
  use Freefall_Radii

  !# <coolingInfallRadius name="coolingInfallRadiusCoolingFreefall" defaultThreadPrivate="yes">
  !#  <description>An infall radius calculation in which the infall radius is the smaller of the cooling and freefall radii.</description>
  !# </coolingInfallRadius>
  type, extends(coolingInfallRadiusClass) :: coolingInfallRadiusCoolingFreefall
     !% Implementation of an infall radius calculation in which the infall radius is the smaller of the cooling and freefall radii.
     private
     class(coolingRadiusClass), pointer :: coolingRadius_
   contains
     final     ::                       coolingFreefallDestructor
     procedure :: radius             => coolingFreefallRadius
     procedure :: radiusIncreaseRate => coolingFreefallRadiusIncreaseRate
  end type coolingInfallRadiusCoolingFreefall

  interface coolingInfallRadiusCoolingFreefall
     !% Constructors for the cooling radius infall radii class.
     module procedure coolingFreefallConstructorParameters
     module procedure coolingFreefallConstructorInternal
  end interface coolingInfallRadiusCoolingFreefall

contains

  function coolingFreefallConstructorParameters(parameters) result(self)
    !% Constructor for the cooling radius infall radii class which builds the object from a parameter set.
    use Input_Parameters2
    implicit none
    type (coolingInfallRadiusCoolingFreefall)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(coolingRadiusClass                ), pointer       :: coolingRadius_

    !# <objectBuilder class="coolingRadius" name="coolingRadius_" source="parameters"/>
    self=coolingInfallRadiusCoolingFreefall(coolingRadius_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function coolingFreefallConstructorParameters

  function coolingFreefallConstructorInternal(coolingRadius_) result(self)
    !% Internal constructor for the cooling radius infall radii class.
    use Galacticus_Error
    implicit none
    type (coolingInfallRadiusCoolingFreefall)                        :: self
    class(coolingRadiusClass                ), intent(in   ), target :: coolingRadius_
    !# <constructorAssign variables="*coolingRadius_"/>

    return
  end function coolingFreefallConstructorInternal

  subroutine coolingFreefallDestructor(self)
    !% Destructor for the cooling radius infall radii class.
    implicit none
    type(coolingInfallRadiusCoolingFreefall), intent(inout) :: self

    !# <objectDestructor name="self%coolingRadius_" />
    return
  end subroutine coolingFreefallDestructor

  double precision function coolingFreefallRadius(self,node)
    !% Return the infall radius in the ``cooling radius'' model in Mpc/Gyr.
    implicit none
    class           (coolingInfallRadiusCoolingFreefall), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                                    :: radiusCooling, radiusFreefall

 
    radiusCooling =self%coolingRadius_%radius(node)
    radiusFreefall=Freefall_Radius           (node)
    if (radiusCooling < radiusFreefall) then
       coolingFreefallRadius=radiusCooling
    else
       coolingFreefallRadius=radiusFreefall
    end if
    return
  end function coolingFreefallRadius

  double precision function coolingFreefallRadiusIncreaseRate(self,node)
    !% Return the growth rate of the infall radius in the ``cooling radius'' model in Mpc/Gyr.
    implicit none
    class           (coolingInfallRadiusCoolingFreefall), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                                    :: radiusCooling, radiusFreefall

    radiusCooling =self%coolingRadius_%radius(node)
    radiusFreefall=Freefall_Radius           (node)
    if (radiusCooling < radiusFreefall) then
       coolingFreefallRadiusIncreaseRate=self%coolingRadius_%radiusGrowthRate(node)
    else
       coolingFreefallRadiusIncreaseRate=Freefall_Radius_Growth_Rate         (node)
    end if
    return
  end function coolingFreefallRadiusIncreaseRate
