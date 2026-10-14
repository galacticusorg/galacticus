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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implementation of a gaseous dynamical friction class which uses the collisionless formula of
  :cite:t:`chandrasekhar_dynamical_1943`.
  !!}

  !![
  <dynamicalFrictionGaseous name="dynamicalFrictionGaseousChandrasekhar1943" docformat="rst">
   <description>
   A gaseous dynamical friction class which treats the gas as a collisionless medium, using the
   :cite:t:`chandrasekhar_dynamical_1943` formula with the sound speed in place of the velocity dispersion:

   .. math::

      \mathcal{I} = \ln\Lambda \left[ \mathrm{erf}(x) - \frac{2x}{\sqrt{\pi}} \exp(-x^2) \right],

   where :math:`x=\mathcal{M}/\sqrt{2}`. This is not a physically-appropriate model for a gaseous medium, but is provided for
   comparison with the collisionless case.
   </description>
  </dynamicalFrictionGaseous>
  !!]
  type, extends(dynamicalFrictionGaseousClass) :: dynamicalFrictionGaseousChandrasekhar1943
     !!{RST
     Implementation of a gaseous dynamical friction class which uses the collisionless formula of
     :cite:t:`chandrasekhar_dynamical_1943`.
     !!}
     private
   contains
     procedure :: forceFactor => chandrasekhar1943ForceFactor
  end type dynamicalFrictionGaseousChandrasekhar1943

  interface dynamicalFrictionGaseousChandrasekhar1943
     !!{RST
     Constructors for the :galacticus-class:`dynamicalFrictionGaseousChandrasekhar1943` gaseous dynamical friction class.
     !!}
     module procedure chandrasekhar1943ConstructorParameters
  end interface dynamicalFrictionGaseousChandrasekhar1943

contains

  function chandrasekhar1943ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dynamicalFrictionGaseousChandrasekhar1943` gaseous dynamical friction class which
    builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(dynamicalFrictionGaseousChandrasekhar1943)                :: self
    type(inputParameters                          ), intent(inout) :: parameters

    self=dynamicalFrictionGaseousChandrasekhar1943()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function chandrasekhar1943ConstructorParameters

  double precision function chandrasekhar1943ForceFactor(self,machNumber,logarithmCoulomb) result(forceFactor)
    !!{RST
    Return the dimensionless drag force for the :cite:t:`chandrasekhar_dynamical_1943` formula.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (dynamicalFrictionGaseousChandrasekhar1943), intent(inout) :: self
    double precision                                           , intent(in   ) :: machNumber       , logarithmCoulomb
    double precision                                           , parameter     :: xMaximum  =10.0d0
    double precision                                                           :: x
    !$GLC attributes unused :: self

    ! For large x the bracketed term is unity to within machine precision.
    x=machNumber/sqrt(2.0d0)
    if (x <= xMaximum) then
       forceFactor=+logarithmCoulomb      &
            &      *(                     &
            &        +erf ( x   )         &
            &        -2.0d0               &
            &        *      x             &
            &        *exp (-x**2)         &
            &        /sqrt( Pi  )         &
            &       )
    else
       forceFactor=+logarithmCoulomb
    end if
    return
  end function chandrasekhar1943ForceFactor
