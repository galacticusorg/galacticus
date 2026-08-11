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

!!{RST
Contains a module which computes the reference mass and radius, :math:`M_\mathrm{200m}` and :math:`R_\mathrm{200m}`, used by
the splashback radius classes.

This duplicates a small amount of the functionality of
:galacticus-function:`Dark_Matter_Profile_Mass_Definition`. That is deliberate: the
:galacticus-class:`darkMatterProfileDMOAccretionFlowDiemerKravtsov2014` and
:galacticus-class:`darkMatterProfileDMOAccretionFlowShi2016` classes make use of splashback radii, so the splashback radius
classes must not themselves depend on the ``darkMatterProfileDMO`` class---which
:galacticus-function:`Dark_Matter_Profile_Mass_Definition` does---as that would introduce a circular module dependency.
!!}

module Dark_Matter_Halo_Splashback_Radii_Reference
  !!{RST
  Computes the reference mass and radius used by the splashback radius classes.
  !!}
  !+ Implemented by Andrew Benson with assistance from Claude.
  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Galacticus_Nodes       , only : treeNode
  use :: Mass_Distributions     , only : massDistributionClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass
  implicit none
  private
  public :: splashbackReferenceScales

  ! The density contrast, relative to the mean density of the universe, defining the reference mass and radius. All of the
  ! splashback models implemented in Galacticus are formulated relative to R₂₀₀ₘ and M₂₀₀ₘ.
  double precision, parameter, public :: splashbackDensityContrastReference=2.0d2

contains

  subroutine splashbackReferenceScales(node,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,mass,radius,massDistribution_)
    !!{RST
    Return the mass, :math:`M_\mathrm{200m}`, and radius, :math:`R_\mathrm{200m}`, enclosing a mean density of 200 times the
    mean density of the universe for the halo in ``node``.

    If ``massDistribution_`` is supplied it is used to find the radius enclosing that density---this allows a caller which
    has already built a mass distribution for this node to obtain results consistent with it, and avoids the infinite
    recursion which would otherwise occur if the ``darkMatterProfileDMO`` class in use itself makes use of splashback radii.
    Otherwise, the dark matter only mass distribution is obtained from ``node``.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkMatterOnly, massTypeDark
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Math_Exponentiation       , only : cubeRoot
    use :: Numerical_Comparison      , only : Values_Agree
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    type            (treeNode                  ), intent(inout)           :: node
    class           (cosmologyParametersClass  ), intent(inout)           :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), intent(inout)           :: cosmologyFunctions_
    class           (virialDensityContrastClass), intent(inout)           :: virialDensityContrast_
    double precision                            , intent(  out)           :: mass                  , radius
    class           (massDistributionClass     ), intent(inout), optional :: massDistribution_
    class           (massDistributionClass     ), pointer                 :: massDistributionNode_
    class           (nodeComponentBasic        ), pointer                 :: basic
    double precision                                                      :: density               , time

    basic   =>  node                 %basic                             (    )
    time    =   basic                %time                              (    )
    density =  +                      splashbackDensityContrastReference          &
         &     *cosmologyParameters_%omegaMatter                        (    )    &
         &     *cosmologyParameters_%densityCritical                    (    )    &
         &     /cosmologyFunctions_ %expansionFactor                    (time)**3
    ! If the requested density contrast matches the virial density contrast of this halo then the mass is simply the mass of
    ! the halo, and no mass distribution is needed.
    if     (                                                                                                    &
         &   Values_Agree(                                                                                      &
         &                       virialDensityContrast_%densityContrast(basic%mass(),basic%timeLastIsolated()), &
         &                       splashbackDensityContrastReference                                           , &
         &                relTol=1.0d-4                                                                         &
         &               )                                                                                      &
         &  .and.                                                                                               &
         &   time == basic%timeLastIsolated()                                                                   &
         & ) then
       mass  =basic%mass()
       radius=cubeRoot(               &
            &          +3.0d0         &
            &          /4.0d0         &
            &          /Pi            &
            &          *mass          &
            &          /density       &
            &         )
    else
       if (present(massDistribution_)) then
          radius               =  massDistribution_     %radiusEnclosingDensity(density                                 )
       else
          massDistributionNode_ => node                 %massDistribution      (componentTypeDarkMatterOnly,massTypeDark)
          radius                =  massDistributionNode_%radiusEnclosingDensity(density                                 )
          !![
          <objectDestructor name="massDistributionNode_"/>
          !!]
       end if
       mass  =+4.0d0     &
            & *Pi        &
            & *density   &
            & *radius**3 &
            & /3.0d0
    end if
    return
  end subroutine splashbackReferenceScales

end module Dark_Matter_Halo_Splashback_Radii_Reference
