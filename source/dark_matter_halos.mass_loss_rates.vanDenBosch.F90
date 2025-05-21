!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implementation of mass loss rates from dark matter halos using the prescription of \cite{van_den_bosch_mass_2005}.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass

  !![
  <darkMatterHaloMassLossRate name="darkMatterHaloMassLossRateVanDenBosch">
   <description>
    A dark matter halo mass loss rate class which uses the algorithm of \cite{van_den_bosch_mass_2005} to compute the rate of
    mass loss. Specifically:
    \begin{equation}
    \dot{M}_\mathrm{node,bound} = -{M_\mathrm{node,bound}\over \tau} \left({M_\mathrm{node,bound} / M_\mathrm{node,parent}}\right)^\zeta,
    \end{equation}
    where $M_\mathrm{node,parent}$ is the mass of the parent \gls{node} in which the halo lives and
    \begin{equation}
    \tau = \tau_0 \left({\Delta_\mathrm{vir}(t) \over \Delta(t_0)}\right)^{-1/2} a^{3/2},
    \end{equation}
    where $\Delta_\mathrm{vir}(t)$ is the virial overdensity of halos at time $t$ and $a$ is the expansion factor. The fitting
    parameters, $\tau_0$ and $\zeta$ have values of 0.13~Gyr and 0.36 respectively as determined by
    \cite{van_den_bosch_mass_2005}. Note that \cite{van_den_bosch_mass_2005} write this expression in a slightly different form
    since their $\Delta_\mathrm{vir}$ is defined relative to the critical density rather than the mean density as it is in
    \glc. In both cases, the timescale $\tau$ simply scales as $\langle \rho_\mathrm{vir} \rangle ^{-1/2}$ where $\langle
    \rho_\mathrm{vir} \rangle$ is the mean virial overdensity of halos.
   </description>
  </darkMatterHaloMassLossRate>
  !!]
  type, extends(darkMatterHaloMassLossRateClass) :: darkMatterHaloMassLossRateVanDenBosch
     !!{
     Implementation of a dark matter halo mass loss rate class which uses the prescription of \cite{van_den_bosch_mass_2005}.
     !!}
     private
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_ => null()
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
     double precision                                      :: timescaleNormalization
     double precision                                      :: zeta
   contains
     final     ::         vanDenBoschDestructor
     procedure :: rate => vanDenBoschRate
  end type darkMatterHaloMassLossRateVanDenBosch

  interface darkMatterHaloMassLossRateVanDenBosch
     !!{
     Constructors for the vanDenBosch dark matter halo mass loss rate class.
     !!}
     module procedure vanDenBoschConstructorParameters
     module procedure vanDenBoschConstructorInternal
  end interface darkMatterHaloMassLossRateVanDenBosch

contains

  function vanDenBoschConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloMassLossRateVanDenBosch} dark matter halo mass loss rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterHaloMassLossRateVanDenBosch)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class           (virialDensityContrastClass           ), pointer       :: virialDensityContrast_
    double precision                                                       :: timescaleNormalization, zeta

    !![
    <inputParameter>
      <name>timescaleNormalization</name>
      <source>parameters</source>
      <defaultValue>0.13d0</defaultValue>
      <description>The mass loss timescale normalization (in Gyr) for the \cite{van_den_bosch_mass_2005} dark matter halo mass loss rate algorithm.</description>
    </inputParameter>
    <inputParameter>
      <name>zeta</name>
      <source>parameters</source>
      <defaultValue>0.36d0</defaultValue>
      <description>The mass loss scaling with halo mass for the \cite{van_den_bosch_mass_2005} dark matter halo mass loss rate algorithm.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=darkMatterHaloMassLossRateVanDenBosch(timescaleNormalization,zeta,cosmologyFunctions_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function vanDenBoschConstructorParameters

  function vanDenBoschConstructorInternal(timescaleNormalization,zeta,cosmologyFunctions_,virialDensityContrast_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloMassLossRateVanDenBosch} dark matter halo mass loss rate class.
    !!}
    implicit none
    type            (darkMatterHaloMassLossRateVanDenBosch)                        :: self
    double precision                                       , intent(in   )         :: timescaleNormalization, zeta
    class           (cosmologyFunctionsClass              ), intent(in   ), target :: cosmologyFunctions_
    class           (virialDensityContrastClass           ), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="timescaleNormalization, zeta, *cosmologyFunctions_, *virialDensityContrast_"/>
    !!]

    return
  end function vanDenBoschConstructorInternal

  subroutine vanDenBoschDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloMassLossRateVanDenBosch} dark matter halo mass loss rate class.
    !!}
    implicit none
    type(darkMatterHaloMassLossRateVanDenBosch), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"    />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine vanDenBoschDestructor

  double precision function vanDenBoschRate(self,node)
    !!{
    Returns the mass loss rate from the dark matter halo of the given \gls{node} in units of $M_\odot$/Gyr.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, treeNode
    implicit none
    class           (darkMatterHaloMassLossRateVanDenBosch), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    class           (nodeComponentBasic                   ), pointer       :: basicParent           , basic
    class           (nodeComponentSatellite               ), pointer       :: satellite
    double precision                                                       :: timescaleMassLoss     , massSatelliteBound, &
         &                                                                    ratioMassSatelliteHost, timeSatellite

    satellite          => node     %satellite()
    massSatelliteBound =  satellite%boundMass()
    if (massSatelliteBound > 0.0d0) then
       basic                  =>  node        %basic()
       basicParent            =>  node %parent%basic()
       timeSatellite          =   basic       %time ()
       timescaleMassLoss      =  +     self%timescaleNormalization                                                                                          &
            &                    *sqrt(self%virialDensityContrast_%densityContrast(basic%mass(),self%cosmologyFunctions_%cosmicTime(1.0d0        )))        &
            &                    *     self%cosmologyFunctions_   %expansionFactor(                                                 timeSatellite ) **1.5d0 &
            &                    /sqrt(self%virialDensityContrast_%densityContrast(basic%mass(),                                    timeSatellite ))
       ratioMassSatelliteHost =  +massSatelliteBound        &
            &                    /basicParent       %mass()
       vanDenBoschRate        =  -massSatelliteBound                &
            &                    *ratioMassSatelliteHost**self%zeta &
            &                    /timescaleMassLoss
    else
       vanDenBoschRate        =  +0.0d0
    end if
    return
  end function vanDenBoschRate
