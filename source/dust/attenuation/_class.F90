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
Contains a module which provides a class implementing attenuation of galactic emission by dust.
!!}

module Dust_Attenuations
  !!{RST
  Provides a class implementing attenuation of galactic emission by dust.

  A ``dustAttenuation`` object answers one question: of the light emitted by a given parcel of a given galaxy, what
  fraction escapes? It is handed a ``node`` and a list of ``emissionDescriptor`` parcels---each one
  wavelength, one component, one source type, one age range---and returns a transmission factor in :math:`[0,1]` for
  each.

  Two design points are worth stating explicitly.

  **Transmission, not magnitudes.** The method returns a multiplicative factor rather than an attenuation in
  magnitudes or an optical depth. This is unambiguous, composes by multiplication (so a sequence of attenuators is
  simply a product), and remains well defined for models---such as radiative transfer through a specified
  geometry---whose result is not proportional to any single optical depth.

  **Vectorized over parcels.** ``transmission`` receives every parcel at once so that an implementation can hoist its
  per-node work---the :math:`V`-band optical depth of each component, the inclination, interpolation factors in a
  tabulated atlas---out of the loop over wavelength, without caching state on the object and risking staleness.

  How finely a producer must split its luminosity is negotiated: ``request`` reports the age bin boundaries an
  attenuator distinguishes and which other axes it depends upon, so that an age-independent screen applied to a
  broad-band luminosity is handed one parcel per component, rather than one per wavelength, age and metallicity.

  Note that any constant which must be shared between implementations of this class has to be declared here rather
  than in an implementation file, because each implementation is generated into its own submodule and so cannot see
  module-level declarations made by its siblings.
  !!}
  use :: Dust_Attenuation_Descriptors, only : decompositionRequest, emissionDescriptor
  use :: Galactic_Structure_Options  , only : componentTypeAll    , enumerationComponentTypeType
  use :: Galacticus_Nodes            , only : treeNode
  private
  ! Made public so that it survives into the object file: it is called only from the submodules into which the
  ! implementations of this class are generated, never from this module itself, and a private procedure with no
  ! caller in its own module can be discarded before those submodules are linked against it.
  public :: componentGasProperties

  ! Metallicity of the local interstellar medium, by mass. Dust-to-gas ratios are scaled relative to this value, on
  ! the assumption that the dust-to-metals ratio is universal. Declared here, rather than in an implementation file,
  ! because implementations are generated into sibling submodules which can see this module but not each other.
  double precision, parameter, public :: metallicityISMLocal=2.0d-2

  !![
  <functionClass docformat="rst">
   <name>dustAttenuation</name>
   <descriptiveName>Dust Attenuation</descriptiveName>
   <description>
   Class computing the fraction of the emission from a galaxy which is transmitted through its dust.
   </description>
   <default>zero</default>
   <method name="transmission" >
    <description>
    Return the fraction of the emission transmitted through dust, in the range :math:`[0,1]`, for each of the given
    parcels of emission from the given ``node``. A value of unity indicates no attenuation. The result has one element
    per element of ``descriptors``, in the same order.

    ``inclination`` overrides the angle at which an orientation-dependent attenuator is evaluated, in radians. It is
    how :galacticus-class:`dustAttenuationInclinationAveraged` drives its quadrature: the wrapped attenuator is
    asked for the transmission at each angle in turn, without any shared state being mutated, so the arrangement is
    safe under threading. An implementation which does not depend on orientation ignores it; one which does falls
    back to its own ``galacticInclination`` object when it is absent.
    </description>
    <type>double precision, dimension(size(descriptors))</type>
    <pass>yes</pass>
    <argument>type(treeNode          ), intent(inout), target       :: node       </argument>
    <argument>type(emissionDescriptor), intent(in   ), dimension(:) :: descriptors</argument>
    <argument>double precision        , intent(in   ), optional     :: inclination</argument>
   </method>
   <method name="request" >
    <description>
    Return the resolution which a decomposition must have for this attenuator to be applied correctly: the boundaries
    of the age bins it distinguishes, and whether it depends upon component, metallicity, or radius. A producer of
    luminosities uses this to split its output no more finely than necessary.

    The default requests no age resolution and no metallicity or radius resolution, but does request splitting by
    component---which is what any attenuator depending on the properties of an individual component needs.
    </description>
    <type>type(decompositionRequest)</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     dustAttenuationRequest%resolveComponents =.true.
     dustAttenuationRequest%resolveMetallicity=.false.
     dustAttenuationRequest%resolveRadius     =.false.
    </code>
   </method>
   <method name="supportsComponent" >
    <description>
    Return true if this attenuator can be applied to emission from the given component.

    The default accepts any individual component but rejects ``componentTypeAll``: different components are attenuated
    differently, so a luminosity summed over components can not meaningfully be attenuated---it must be decomposed,
    attenuated, and only then summed. This is checked when the attenuating property extractor is constructed, so that
    a misconfiguration is reported immediately rather than part way through a run.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(enumerationComponentTypeType), intent(in   ) :: componentType</argument>
    <code>
     !$GLC attributes unused :: self
     dustAttenuationSupportsComponent=(componentType /= componentTypeAll)
    </code>
   </method>
  </functionClass>
  !!]

contains

  subroutine componentGasProperties(node,componentType,massGas,radius,metallicity)
    !!{RST
    Return the gas mass (:math:`M_\odot`), scale radius (Mpc), and gas-phase metallicity (linear, by mass) of the
    given component of the given ``node``.

    Attenuators which scale the dust content with the gas content of a component all need these three quantities, so
    the extraction lives here rather than in any one implementation: each implementation of this class is generated
    into its own submodule, and while a submodule can see its parent module it can not see its siblings.

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

end module Dust_Attenuations
