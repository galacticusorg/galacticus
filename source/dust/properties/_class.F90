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
Contains a module which provides a class implementing the properties of interstellar dust.
!!}

module Dust_Properties
  !!{RST
  Provides a class implementing the properties of interstellar dust---how much of it a galaxy has, and how strongly it
  extinguishes light.

  Dust enters a model in two places: it attenuates the light of stars and gas, and it re-emits the energy it absorbs in
  the infrared. Both depend on the same quantities---the mass of dust, set by a dust-to-metals ratio, and its opacity per
  unit mass---and were the two calculations to take those from different places they would disagree silently, the
  infrared emission of a galaxy corresponding to a different mass of dust from the one which reddened its starlight.
  Dust attenuation and emission models therefore obtain these quantities from an object of this class rather than
  carrying their own.

  The opacity here is a property of the dust grains. It is distinct from the extinction curve of a ``dustAttenuation``
  model, which may be an *effective* attenuation law folding in the geometry of the dust and stars---the
  :math:`\lambda^{-0.7}` law of :cite:t:`charlot_simple_2000`, for example---and so need not share the wavelength
  dependence of the grains themselves.
  !!}
  use :: Galactic_Structure_Options      , only : enumerationComponentTypeType
  use :: Galacticus_Nodes                , only : treeNode
  use :: Numerical_Constants_Astronomical, only : hydrogenByMassSolar         , massSolar, opticalDepthToMagnitudes, parsec
  use :: Numerical_Constants_Atomic      , only : atomicMassUnit
  use :: Numerical_Constants_Prefixes    , only : hecto                       , kilo
  private
  ! Made public so that it survives into the object file, and so that dust attenuation classes, which need the same
  ! quantities, can use it too.
  public :: componentGasProperties

  ! Metallicity of the local interstellar medium, by mass.
  double precision, parameter, public :: metallicityISMLocal                         =2.0d-2
  ! Dust-to-metals ratio of the Milky Way (Popping et al. 2017; MNRAS; 471; 3152).
  double precision, parameter, public :: dustToMetalsRatioMilkyWay                   =0.44d0
  ! A_V/E(B-V), and N_H/E(B-V) in atoms/cm²/mag (Savage & Mathis 1979; ARA&A; 17; 73).
  double precision, parameter         :: AVToEBV                                     =3.10d+00
  double precision, parameter         :: NHToEBV                                     =5.80d+21
  ! V-band optical depth per unit surface density of metals, in cm² g⁻¹, for Milky Way dust: the relation between column
  ! density and reddening of Savage & Mathis (1979), for gas of the metallicity of the local interstellar medium.
  double precision, parameter, public :: depthOpticalVPerSurfaceDensityMetalsMilkyWay=+AVToEBV                  &
       &                                                                              /NHToEBV                  &
       &                                                                              *hydrogenByMassSolar      &
       &                                                                              /atomicMassUnit           &
       &                                                                              /kilo                     &
       &                                                                              /metallicityISMLocal      &
       &                                                                              /opticalDepthToMagnitudes
  ! Surface density of gas, in M☉ pc⁻², for which gas of the metallicity of the local interstellar medium, containing
  ! Milky Way dust, has unit V-band optical depth.
  double precision, parameter, public :: densitySurfaceGasDepthOpticalVUnitMilkyWay  =+1.0d0                                        &
       &                                                                              /depthOpticalVPerSurfaceDensityMetalsMilkyWay &
       &                                                                              /metallicityISMLocal                          &
       &                                                                              *(                                            &
       &                                                                                +parsec                                     &
       &                                                                                *hecto                                      &
       &                                                                               )**2                                         &
       &                                                                              /massSolar                                    &
       &                                                                              /kilo

  !![
  <functionClass docformat="rst">
   <name>dustProperties</name>
   <descriptiveName>Dust Properties</descriptiveName>
   <description>
   Class providing the properties of interstellar dust: its abundance relative to metals, its mass, and its opacity.
   </description>
   <default>simple</default>
   <method name="dustToMetalsRatio" >
    <description>
    Return the fraction of the mass of metals in the gas of the given component of ``node`` which is in dust.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout), target :: node         </argument>
    <argument>type(enumerationComponentTypeType), intent(in   )         :: componentType</argument>
   </method>
   <method name="opacityExtinctionV" >
    <description>
    Return the :math:`V`-band extinction opacity per unit mass of dust, in cm² g⁻¹: the optical depth, due to
    absorption and scattering together, of a column containing one gram of dust per square centimeter.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="massDust" >
    <description>
    Return the mass of dust, in :math:`M_\odot`, in the gas of the given component of ``node``. The default is the
    product of the dust-to-metals ratio and the mass of metals in the gas, and is zero for a component with no gas or
    no size.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode                    ), intent(inout), target :: node         </argument>
    <argument>type(enumerationComponentTypeType), intent(in   )         :: componentType</argument>
    <code>
     double precision :: massGas, radius, metallicity
     call componentGasProperties(node,componentType,massGas,radius,metallicity)
     dustPropertiesMassDust=+self%dustToMetalsRatio(node,componentType)*metallicity*massGas
    </code>
   </method>
  </functionClass>
  !!]

contains

  subroutine componentGasProperties(node,componentType,massGas,radius,metallicity)
    !!{RST
    Return the gas mass (:math:`M_\odot`), scale radius (Mpc), and gas-phase metallicity (linear, by mass) of the
    given component of the given ``node``.

    Dust properties and dust attenuation classes which scale the dust content with the gas content of a component all
    need these three quantities, so the extraction lives here rather than in any one implementation.

    A component with no gas or no size returns zero for all three, which callers should treat as containing no dust.
    !!}
    use :: Abundances_Structure      , only : abundances       , metallicityTypeLinearByMass
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeNuclearStarCluster, componentTypeSpheroid
    use :: Galacticus_Nodes          , only : nodeComponentDisk, nodeComponentNSC               , nodeComponentSpheroid
    implicit none
    type            (treeNode                    ), intent(inout), target  :: node
    type            (enumerationComponentTypeType), intent(in   )          :: componentType
    double precision                              , intent(  out)          :: massGas           , radius, &
         &                                                                    metallicity
    class           (nodeComponentDisk           )               , pointer :: disk
    class           (nodeComponentSpheroid       )               , pointer :: spheroid
    class           (nodeComponentNSC            )               , pointer :: nuclearStarCluster
    type            (abundances                  )                         :: abundancesGas

    select case (componentType%ID)
    case (componentTypeDisk              %ID)
       disk               => node              %disk         ()
       massGas            =  disk              %massGas      ()
       radius             =  disk              %radius       ()
       abundancesGas      =  disk              %abundancesGas()
    case (componentTypeSpheroid          %ID)
       spheroid           => node              %spheroid     ()
       massGas            =  spheroid          %massGas      ()
       radius             =  spheroid          %radius       ()
       abundancesGas      =  spheroid          %abundancesGas()
    case (componentTypeNuclearStarCluster%ID)
       nuclearStarCluster => node              %NSC          ()
       massGas            =  nuclearStarCluster%massGas      ()
       radius             =  nuclearStarCluster%radius       ()
       abundancesGas      =  nuclearStarCluster%abundancesGas()
    case default
       massGas            =  0.0d0
       radius             =  0.0d0
       metallicity        =  0.0d0
       call Error_Report('component can not host dust'//{introspection:location})
       return
    end select
    if (massGas <= 0.0d0 .or. radius <= 0.0d0) then
       massGas    =0.0d0
       radius     =0.0d0
       metallicity=0.0d0
       return
    end if
    call abundancesGas%massToMassFraction(massGas)
    metallicity=abundancesGas%metallicity(metallicityTypeLinearByMass)
    return
  end subroutine componentGasProperties

end module Dust_Properties
