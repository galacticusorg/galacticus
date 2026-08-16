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
Contains a module which provides a class implementing splashback radii of dark matter halos.
!!}

module Dark_Matter_Halo_Splashback_Radii
  !!{RST
  Provides a class implementing splashback radii of dark matter halos.
  !!}
  !+ Contributions to this file made by: Andrew Benson, Claude.
  use :: Galacticus_Nodes  , only : treeNode
  use :: Mass_Distributions, only : massDistributionClass
  private

  !![
  <functionClass docformat="rst">
   <name>darkMatterHaloSplashbackRadius</name>
   <descriptiveName>Dark Matter Halo Splashback Radii</descriptiveName>
   <description>
   Class providing the splashback radius, :math:`R_\mathrm{sp}`, of dark matter halos---the radius at which accreted material
   reaches the apocenter of its first orbit after collapse. Physically this marks the true outer edge of the collapsed halo,
   and manifests observationally as a sharp steepening of the density profile. It is not a spherical overdensity radius, and
   so does not correspond to any fixed density contrast: it depends primarily on the halo's mass accretion rate,
   :math:`\Gamma`, with rapidly accreting halos having smaller :math:`R_\mathrm{sp}/R_\mathrm{200m}`.

   Implementations of this class differ both in the model used and, for those models calibrated against the distribution of
   particle apocenters, in which statistic of that distribution is being predicted (see the ``definition`` parameter of the
   :galacticus-class:`darkMatterHaloSplashbackRadiusDiemer2017` and
   :galacticus-class:`darkMatterHaloSplashbackRadiusDiemer2020` classes).
   </description>
   <default>diemer2020</default>
   <method name="radius" >
    <description>
    Returns the splashback radius, :math:`R_\mathrm{sp}` (in Mpc), of the halo in ``node``.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type (treeNode             ), intent(inout)           :: node</argument>
    <argument>class(massDistributionClass), intent(inout), optional :: massDistribution_</argument>
   </method>
   <method name="radiusRatio" >
    <description>
    Returns the ratio :math:`R_\mathrm{sp}/R_\mathrm{200m}` for the halo in ``node``. Splashback models are formulated in terms
    of this ratio, so this method allows a caller which has already computed :math:`R_\mathrm{200m}` for its own density
    profile to obtain a splashback radius consistent with that profile, instead of with the profile used internally by this
    class.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type (treeNode             ), intent(inout)           :: node</argument>
    <argument>class(massDistributionClass), intent(inout), optional :: massDistribution_</argument>
   </method>
   <method name="mass" >
    <description>
    Returns the splashback mass, :math:`M_\mathrm{sp}` (in :math:`\mathrm{M}_\odot`), of the halo in ``node``---that is, the
    mass enclosed by the splashback radius. Note that this is predicted directly by the fitting functions of those models
    which provide it, and is in general *not* equal to the mass enclosed within :math:`R_\mathrm{sp}` in any given density
    profile model. An error is reported by models which do not predict this quantity.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Error</modules>
    <argument>type (treeNode             ), intent(inout)           :: node</argument>
    <argument>class(massDistributionClass), intent(inout), optional :: massDistribution_</argument>
    <code>
     !$GLC attributes unused :: self, node, massDistribution_
     darkMatterHaloSplashbackRadiusMass=0.0d0
     call Error_Report('splashback mass is not provided by the "'//char(self%objectType())//'" splashback radius class'//{introspection:location})
    </code>
   </method>
   <method name="massRatio" >
    <description>
    Returns the ratio :math:`M_\mathrm{sp}/M_\mathrm{200m}` for the halo in ``node``. An error is reported by models which do
    not predict this quantity.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Error</modules>
    <argument>type (treeNode             ), intent(inout)           :: node</argument>
    <argument>class(massDistributionClass), intent(inout), optional :: massDistribution_</argument>
    <code>
     !$GLC attributes unused :: self, node, massDistribution_
     darkMatterHaloSplashbackRadiusMassRatio=0.0d0
     call Error_Report('splashback mass is not provided by the "'//char(self%objectType())//'" splashback radius class'//{introspection:location})
    </code>
   </method>
   <method name="radiusScatter" >
    <description>
    Returns the 68% scatter in :math:`\log_{10}(R_\mathrm{sp}/R_\mathrm{200m})` for halos matching those in ``node``. An
    error is reported by models which do not predict this scatter.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Error</modules>
    <argument>type (treeNode             ), intent(inout)           :: node</argument>
    <argument>class(massDistributionClass), intent(inout), optional :: massDistribution_</argument>
    <code>
     !$GLC attributes unused :: self, node, massDistribution_
     darkMatterHaloSplashbackRadiusRadiusScatter=0.0d0
     call Error_Report('scatter in splashback radius is not provided by the "'//char(self%objectType())//'" splashback radius class'//{introspection:location})
    </code>
   </method>
   <method name="massScatter" >
    <description>
    Returns the 68% scatter in :math:`\log_{10}(M_\mathrm{sp}/M_\mathrm{200m})` for halos matching those in ``node``. An
    error is reported by models which do not predict this scatter.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Error</modules>
    <argument>type (treeNode             ), intent(inout)           :: node</argument>
    <argument>class(massDistributionClass), intent(inout), optional :: massDistribution_</argument>
    <code>
     !$GLC attributes unused :: self, node, massDistribution_
     darkMatterHaloSplashbackRadiusMassScatter=0.0d0
     call Error_Report('scatter in splashback mass is not provided by the "'//char(self%objectType())//'" splashback radius class'//{introspection:location})
    </code>
   </method>
  </functionClass>
  !!]

end module Dark_Matter_Halo_Splashback_Radii
