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
  Implements a \cite{tinker_towardhalo_2008} dark matter halo mass function class.
  !!}
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Linear_Growth             , only : linearGrowthClass

  !![
  <haloMassFunction name="haloMassFunctionTinker2008Form" abstract="yes">
   <description>The halo mass function is computed from the function given by \cite{tinker_towardhalo_2008}.</description>
  </haloMassFunction>
  !!]
  type, abstract, extends(haloMassFunctionClass) :: haloMassFunctionTinker2008Form
     !!{
     A halo mass function class using the fitting function of \cite{tinker_towardhalo_2008}.
     !!}
     private
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (linearGrowthClass            ), pointer :: linearGrowth_             => null()
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     double precision                                         :: time                               , mass, &
          &                                                      massFunction
   contains
     !![
     <methods>
       <method description="Return the parameter $a$ in the \cite{tinker_towardhalo_2008} halo mass function fit." method="a"             />
       <method description="Return the parameter $b$ in the \cite{tinker_towardhalo_2008} halo mass function fit." method="b"             />
       <method description="Return the parameter $c$ in the \cite{tinker_towardhalo_2008} halo mass function fit." method="c"             />
       <method description="Return the parameter $A$ in the \cite{tinker_towardhalo_2008} halo mass function fit." method="normalization" />
     </methods>
     !!]
     procedure                                    :: differential  => tinker2008FormDifferential
     procedure(tinker2008FormParameter), deferred :: normalization
     procedure(tinker2008FormParameter), deferred :: a
     procedure(tinker2008FormParameter), deferred :: b
     procedure(tinker2008FormParameter), deferred :: c
  end type haloMassFunctionTinker2008Form

  abstract interface
     !!{
     Interface to parameter functions.
     !!}
     double precision function tinker2008FormParameter(self,time,mass)
       import haloMassFunctionTinker2008Form
       class           (haloMassFunctionTinker2008Form), intent(inout) :: self
       double precision                                , intent(in   ) :: time, mass
     end function tinker2008FormParameter
  end interface

contains

  double precision function tinker2008FormDifferential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008Form), intent(inout), target    :: self
    double precision                                , intent(in   )            :: time , mass
    type            (treeNode                      ), intent(inout) , optional :: node
    double precision                                                           :: sigma, alpha
    !$GLC attributes unused :: node

    ! Update fitting function parameters if the time differs from that on the previous call.
    if (time /= self%time .or. mass /= self%mass) then
       ! Compute and store the mass function.
       sigma            =+    self%cosmologicalMassVariance_%rootVariance                   (mass,time)
       alpha            =+abs(self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass,time))
       self%massFunction=+    self%cosmologyParameters_     %OmegaMatter                    (         )  &
            &            *    self%cosmologyParameters_     %densityCritical                (         )  &
            &            /mass**2                                                                        &
            &            *alpha                                                                          &
            &            *self%normalization(time,mass)                                                  &
            &            *exp(                                                                           &
            &                 -self%c(time,mass)                                                         &
            &                 /sigma**2                                                                  &
            &                )                                                                           &
            &            *(                                                                              &
            &              +1.0d0                                                                        &
            &              +(                                                                            &
            &                +self%b(time,mass)                                                          &
            &                /sigma                                                                      &
            &               )**self%a(time,mass)                                                         &
            &             )
       ! Store the time and mass.
       self%time=time
       self%mass=mass
    end if
    ! Return the stored mass function.
    tinker2008FormDifferential=self%massFunction
    return
  end function tinker2008FormDifferential
