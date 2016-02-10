!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a \cite{press_formation_1974} dark matter halo mass function class.
  use Cosmological_Mass_Variance

  !# <haloMassFunction name="haloMassFunctionPressSchechter">
  !#  <description>The halo mass function is computed from the function given by \cite{press_formation_1974}.</description>
  !# </haloMassFunction>
  type, extends(haloMassFunctionClass) :: haloMassFunctionPressSchechter
     !% A halo mass function class using the model of \cite{press_formation_1974}.
     private
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_
    contains
     final     ::                 pressSchechterDestructor
     procedure :: differential => pressSchechterDifferential
  end type haloMassFunctionPressSchechter

  interface haloMassFunctionPressSchechter
     !% Constructors for the {\normalfont \ttfamily pressSchechter} halo mass function class.
     module procedure pressSchechterConstructorParameters
     module procedure pressSchechterConstructorInternal
  end interface haloMassFunctionPressSchechter

contains

  function pressSchechterConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily pressSchechter} halo mass function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(haloMassFunctionPressSchechter)                :: pressSchechterConstructorParameters
    type(inputParameters               ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <objectBuilder class="cosmologyParameters"      name="pressSchechterConstructorParameters%cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="pressSchechterConstructorParameters%cosmologicalMassVariance_" source="parameters"/>
   return
  end function pressSchechterConstructorParameters

  function pressSchechterConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_)
    !% Internal constructor for the {\normalfont \ttfamily pressSchechter} halo mass function class.
    implicit none
    type (haloMassFunctionPressSchechter)                        :: pressSchechterConstructorInternal
    class(cosmologyParametersClass      ), target, intent(in   ) :: cosmologyParameters_    
    class(cosmologicalMassVarianceClass ), target, intent(in   ) :: cosmologicalMassVariance_

    pressSchechterConstructorInternal%cosmologyParameters_      => cosmologyParameters_
    pressSchechterConstructorInternal%cosmologicalMassVariance_ => cosmologicalMassVariance_
    return
  end function pressSchechterConstructorInternal

  subroutine pressSchechterDestructor(self)
    !% Destructor for the {\normalfont \ttfamily pressSchechter} halo mass function class.
    implicit none
    type(haloMassFunctionPressSchechter), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"      />
    !# <objectDestructor name="self%cosmologicalMassVariance_" />
    return
  end subroutine pressSchechterDestructor

  double precision function pressSchechterDifferential(self,time,mass)
    !% Return the differential halo mass function at the given time and mass.
    use Excursion_Sets_First_Crossings
    implicit none
    class           (haloMassFunctionPressSchechter), intent(inout) :: self
    double precision                                , intent(in   ) :: time , mass    
    double precision                                                :: alpha, variance

    alpha                     =abs(self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass))
    variance                  =    self%cosmologicalMassVariance_%rootVariance                   (mass) **2
    pressSchechterDifferential=+2.0d0&
         &                     *self%cosmologyParameters_%OmegaMatter    (             ) &
         &                     *self%cosmologyParameters_%densityCritical(             ) &
         &                     /mass**2                                                  &
         &                     *alpha                                                    &
         &                     *variance                                                 &
         &                     *Excursion_Sets_First_Crossing_Probability(variance,time)
    return
  end function pressSchechterDifferential
