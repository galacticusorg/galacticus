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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module of utilities needed by ``virialOrbit`` classes.
!!}

module Virial_Orbit_Utilities
  !!{RST
  Provides utilities needed by ``virialOrbit`` classes.

  Each takes the virial ``densityContrast`` defining the orbit class as a value, rather than the object which
  evaluates it, since for some classes (``lossCone``, for example) that object *is* the ``virialDensityContrast_``
  which these also require, and passing one object as two arguments which may be redefined is not permitted.
  !!}
  private
  public :: Virial_Orbit_Energy_Mean             , Virial_Orbit_Angular_Momentum_Magnitude_Mean, &
       &    Virial_Orbit_Angular_Momentum_Reduced, Virial_Orbit_Host_Mass_Radius               , &
       &    Virial_Orbit_Density_Contrast

contains

  double precision function Virial_Orbit_Density_Contrast(host,virialDensityContrastDefinition_) result(densityContrast)
    !!{RST
    Return the virial density contrast defining an orbit class, evaluated for ``host`` at the time it was last isolated.
    !!}
    use :: Galacticus_Nodes       , only : nodeComponentBasic, treeNode
    use :: Virial_Density_Contrast, only : virialDensityContrastClass
    implicit none
    type     (treeNode                  ), intent(inout) :: host
    class    (virialDensityContrastClass), intent(inout) :: virialDensityContrastDefinition_
    class    (nodeComponentBasic        ), pointer       :: basicHost

    basicHost       =>  host%basic()
    densityContrast =   virialDensityContrastDefinition_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated())
    return
  end function Virial_Orbit_Density_Contrast

  subroutine Virial_Orbit_Host_Mass_Radius(host,densityContrast,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_,massHost,radiusHost)
    !!{RST
    Return the mass and radius of ``host`` in the definition specified by the given ``densityContrast``.
    !!}
    use :: Cosmology_Functions                 , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters                , only : cosmologyParametersClass
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Dark_Matter_Profiles_DMO            , only : darkMatterProfileDMOClass
    use :: Galacticus_Nodes                    , only : treeNode
    use :: Virial_Density_Contrast             , only : virialDensityContrastClass
    implicit none
    type            (treeNode                  ), intent(inout) :: host
    double precision                            , intent(in   ) :: densityContrast
    class           (cosmologyParametersClass  ), intent(inout) :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), intent(inout) :: cosmologyFunctions_
    class           (virialDensityContrastClass), intent(inout) :: virialDensityContrast_
    class           (darkMatterProfileDMOClass ), intent(inout) :: darkMatterProfileDMO_
    double precision                            , intent(  out) :: massHost               , radiusHost
    double precision                                            :: velocityHost

    massHost=Dark_Matter_Profile_Mass_Definition(                                                &
         &                                                              host                  , &
         &                                                              densityContrast       , &
         &                                                              radiusHost            , &
         &                                                              velocityHost          , &
         &                                       cosmologyParameters_  =cosmologyParameters_  , &
         &                                       cosmologyFunctions_   =cosmologyFunctions_   , &
         &                                       virialDensityContrast_=virialDensityContrast_, &
         &                                       darkMatterProfileDMO_ =darkMatterProfileDMO_   &
         &                                      )
    return
  end subroutine Virial_Orbit_Host_Mass_Radius

  double precision function Virial_Orbit_Energy_Mean(node,host,velocityTotalRootMeanSquared,densityContrast,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_) result(energyMean)
    !!{RST
    Return the mean energy of orbits of ``node`` about ``host``, given the root mean squared total velocity of those orbits,
    :math:`v_\mathrm{rms}`:

    .. math::
     \langle E \rangle = \frac{v_\mathrm{rms}^2}{2 \left( 1 + M_\mathrm{node}/M_\mathrm{host} \right)} - \frac{\mathrm{G} M_\mathrm{host}}{r_\mathrm{host}},

    where the factor involving the masses accounts for the reduced mass, and the mass and radius of the host are those in the
    definition specified by the given ``densityContrast``.
    !!}
    use :: Cosmology_Functions             , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters            , only : cosmologyParametersClass
    use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOClass
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Virial_Density_Contrast         , only : virialDensityContrastClass
    implicit none
    type            (treeNode                  ), intent(inout) :: node                        , host
    double precision                            , intent(in   ) :: velocityTotalRootMeanSquared, densityContrast
    class           (cosmologyParametersClass  ), intent(inout) :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), intent(inout) :: cosmologyFunctions_
    class           (virialDensityContrastClass), intent(inout) :: virialDensityContrast_
    class           (darkMatterProfileDMOClass ), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic        ), pointer       :: basic                       , basicHost
    double precision                                            :: massHost                    , radiusHost

    basic      =>  node%basic()
    basicHost  =>  host%basic()
    call Virial_Orbit_Host_Mass_Radius(host,densityContrast,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_,massHost,radiusHost)
    energyMean =  +0.5d0                            &
         &        *velocityTotalRootMeanSquared**2  &
         &        /(                                & ! Account for reduced mass.
         &          +1.0d0                          &
         &          +basic    %mass()               &
         &          /basicHost%mass()               &
         &         )                                &
         &        -gravitationalConstant_internal   &
         &        *massHost                         &
         &        /radiusHost
    return
  end function Virial_Orbit_Energy_Mean

  double precision function Virial_Orbit_Angular_Momentum_Reduced(node,host,velocityTangentialMagnitudeMean,radiusHost) result(angularMomentumMagnitudeMean)
    !!{RST
    Return the mean magnitude of the specific angular momentum of orbits of ``node`` about ``host``, given the mean magnitude
    of the tangential velocity of those orbits, :math:`v_\mathrm{t}`, and the radius of the host, :math:`r_\mathrm{host}`:

    .. math::
     \langle |j| \rangle = \frac{v_\mathrm{t} r_\mathrm{host}}{1 + M_\mathrm{node}/M_\mathrm{host}},

    where the denominator accounts for the reduced mass. Classes which must test the mass of the host before evaluating the
    tangential velocity call this together with ``Virial_Orbit_Host_Mass_Radius``; those needing no such test call
    ``Virial_Orbit_Angular_Momentum_Magnitude_Mean`` instead.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    type            (treeNode          ), intent(inout) :: node                           , host
    double precision                    , intent(in   ) :: velocityTangentialMagnitudeMean, radiusHost
    class           (nodeComponentBasic), pointer       :: basic                          , basicHost

    basic                        =>  node%basic()
    basicHost                    =>  host%basic()
    angularMomentumMagnitudeMean =  +velocityTangentialMagnitudeMean &
         &                          *radiusHost                      &
         &                          /(                               & ! Account for reduced mass.
         &                            +1.0d0                         &
         &                            +basic    %mass()              &
         &                            /basicHost%mass()              &
         &                           )
    return
  end function Virial_Orbit_Angular_Momentum_Reduced

  double precision function Virial_Orbit_Angular_Momentum_Magnitude_Mean(node,host,velocityTangentialMagnitudeMean,densityContrast,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_) result(angularMomentumMagnitudeMean)
    !!{RST
    Return the mean magnitude of the specific angular momentum of orbits of ``node`` about ``host``, using the radius of the
    host in the definition specified by the given ``densityContrast``.
    !!}
    use :: Cosmology_Functions     , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters    , only : cosmologyParametersClass
    use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
    use :: Galacticus_Nodes        , only : treeNode
    use :: Virial_Density_Contrast , only : virialDensityContrastClass
    implicit none
    type            (treeNode                  ), intent(inout) :: node                           , host
    double precision                            , intent(in   ) :: velocityTangentialMagnitudeMean, densityContrast
    class           (cosmologyParametersClass  ), intent(inout) :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), intent(inout) :: cosmologyFunctions_
    class           (virialDensityContrastClass), intent(inout) :: virialDensityContrast_
    class           (darkMatterProfileDMOClass ), intent(inout) :: darkMatterProfileDMO_
    double precision                                            :: massHost                       , radiusHost

    call Virial_Orbit_Host_Mass_Radius(host,densityContrast,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_,massHost,radiusHost)
    angularMomentumMagnitudeMean=Virial_Orbit_Angular_Momentum_Reduced(node,host,velocityTangentialMagnitudeMean,radiusHost)
    return
  end function Virial_Orbit_Angular_Momentum_Magnitude_Mean

end module Virial_Orbit_Utilities
