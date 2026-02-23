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
  Implements a dark matter profile scale radius class using the energy-based model of \cite{johnson_random_2021}.
  !!}

  use :: Cosmology_Functions               , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field        , only : cosmologicalMassVarianceClass , criticalOverdensityClass
  use :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO          , only : darkMatterProfileDMOClass
  use :: Galacticus_Nodes                  , only : nodeComponentDarkMatterProfile
  use :: Virial_Orbits                     , only : virialOrbitClass
  use :: Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolutionClass
  use :: Kind_Numbers                      , only : kind_int8

  type :: branch
     !!{
     Type used to store radius offsets along a branch.
     !!}
     private
     integer         (kind_int8)                            :: uniqueID
     double precision           , allocatable, dimension(:) :: radiusOffsetLogarithmic
     type            (branch   ), pointer                   :: next                    => null()
  end type branch

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusJohnson2021">
    <description>
      A dark matter profile scale radius class that computes scale radii based on the energy conservation approach of
      \cite{johnson_random_2021}, with some modifications and improvements.

      Specifically, for ``well-resolved'' halos (see below for discussion of this point), the concentration is found by computing
      the total energy of the halo as the sum over the energies (internal and orbital) of its progenitor halos. The halo is
      assumed to have virialized at this energy, and the scale radius is then solved for such that the energy of the assumed
      density profile (e.g. \glc{nfw}) is equal to the computed energy of the halo.

      In detail, the energy, $E_\mathrm{int}$, of a node is assumed to be given by:
      \begin{equation}
       E_\mathrm{int} = E_{\mathrm{int}, 0} + \sum_{i=1}^N \left( E_{\mathrm{int}, i} +  E_{\mathrm{orb}, i,0} \right) (1 + \mu)^{-\alpha} (1+b \nu^\beta) w_i 10^{\sigma \mathcal{N}(0,1)},
      \end{equation}
      where $E_{\mathrm{int}, 0}$ is the internal energy of the primary progenitor halo, $E_{\mathrm{int}, i}$ is the internal
      energy of the $i^\mathrm{th}$ non-primary progenitor halo, $E_{\mathrm{orb}, i,0}$ is the orbital energy of the
      $i^\mathrm{th}$ non-primary progenitor halo about the primary progenitor halo, $\mu = M_i/M_0$ is the mass ratio of the
      $i^\mathrm{th}$ non-primary progenitor and the primary progenitor ratio, $\nu$ is the peak height parameter for the primary
      progenitor halo, $w_i$ is the subsampling weight of the $i^\mathrm{th}$ non-primary progenitor, $\mathcal{N}(0,1)$ is a
      standard normal deviate, $\alpha=${\normalfont \ttfamily [massExponent]}, $\beta=${\normalfont \ttfamily
      [peakHeightExponent]}, $b=${\normalfont \ttfamily [energyBoost]}, and $\sigma=${\normalfont \ttfamily [scatterExcess]}.

      To account for the contribution to the energy from unresolved accretion, we proceed as follows. First, the unresolved mass
      is determined by subtracting the mass of all progenitors from the halo mass:
      \begin{equation}
      M_{\mathrm unres} = M - \sum_{i=0}^N M_i.
      \end{equation}
      This mass will be accreted in halos spanning a range of masses below the unresolved mass scale, $M_\mathrm{unres}$. We assume the
      mass function of these unresolved halos follows a power-law, $n(M) \propto M^\delta$, with $\delta = -1.8$. We further assume
      that the mean orbit energy of an unresolved halo is simply proportional to $M$ (i.e. the distribution of orbital parameters
      is independent of mass), and that the internal energy of an unresolved halo scales as $E_\mathrm{int} \propto
      M^{5/3+\epsilon}$ where the $5/3$ exponent is the expected scaling for halos assuming a self-similar structure, and
      $\epsilon = -0.02$ accounts for the non-self-similarity (e.g. that concentration is a function of mass) and was estimated
      for typical CDM halos.

      The mean energy (orbital or internal) per unit mass of unresolved halos is then found by averaging these scalings over the
      mass function. Dividing these by the energy of an unresolved halo of mass $M_\mathrm{unres}$ then gives a correction factor
      that can be applied to the energy computed for halos of $M_\mathrm{res}$ to account for the spectrum of unresolved halo
      masses. That is, for an energy that scales as $x^p$ where $x = M/M_0$, the correction factor is:
      \begin{equation}
      c_p = \frac{1}{x_\mathrm{unres}^{p-1} (1+x_\mathrm{unres})^{-\alpha}} \frac{\int_0^{x_\mathrm{unres}} n(x) x^p (1+x)^{-\alpha} \mathrm{d}x}{\int_0^{x_\mathrm{unres}} n(x) x \mathrm{d}x}.
      \end{equation}
      This evaluates to:
      \begin{equation}
      c_p = \frac{2+\delta}{1+\delta+p} (1+x_\mathrm{unres})^\alpha\, _2\mathrm{F}_1({1+\delta+p,\alpha},{2+\delta+p},-x_\mathrm{unres}).
      \end{equation}

      We further compute correction factors to account for the fact that this unresolved mass is accreted over some interval of
      time, during which mean halo properties (e.g. virial density) will evolve. For internal energies we have $E_\mathrm{int}
      \propto M^2/r \propto M^{5/3}/\rho \propto M^{5/3}/a(t)$ where $a(t)$ is the expansion factor. If the current halo exists at
      time $t$, and the primary progenitor at time $t_0$, we therefore compute a correction factor to the midpoint of this
      interval,
      \begin{equation}
      f_\mathrm{int} = \frac{1}{2} \left( 1 + \frac{a(t_0)}{a(t)} \right),
      \end{equation}
      where the second term in the parentheses is the ratio of the internal energy of a halo of fixed mass at $t$ and $t_0$. For
      orbital energy we expect a scaling $E_\mathrm{orb} \propto M M_\mathrm{host}/r_\mathrm{host} \propto M
      M_\mathrm{host}^{2/3}/\rho \propto M M_\mathrm{host}^{2/3}/a(t)$ where $M_\mathrm{host}$ and $r_\mathrm{host}$ are the
      virial mass and radius of the host (current) halo, respectively. Computing a correction factor to the midpoint of the time
      interval we find,
      \begin{equation}
      f_\mathrm{orb} = \frac{1}{2} \left( 1 + \frac{a(t_0)}{a(t)} \left[\frac{M}{M_0}\right]^{2/3} \right).
      \end{equation}

      We next estimate the mean and root-variance, $\bar{E}_\mathrm{unres}$ and $\sigma_\mathrm{unres}$, respectively, of the
      energy of a halo of mass $M_\mathrm{unres}$ via a Monte Carlo approach. We generate $N_\mathrm{MC}=${\normalfont \ttfamily
      [countSampleEnergyUnresolved]} such halos, each with scale radii set using the fall-back \refClass{darkMatterHaloScaleClass}
      object with an added scatter of $\sigma^\prime=${\normalfont \ttfamily [scatter]} dex, and a randomly selected orbit. For
      each such halo, the energy is computed as
      \begin{equation}
      E_\mathrm{unres} = u (E_\mathrm{orb} c_\mathrm{orb} f_\mathrm{orb} + E_\mathrm{int} c_\mathrm{int} f_\mathrm{int}) (1+b \nu^\beta),
      \end{equation}
      where $u = ${\normalfont \ttfamily [unresolvedEnergy]}.

      To estimate the deviation from the mean unresolved energy we again consider the contributions from a spectrum of unresolved
      halo masses. For simplicity we here ignore the internal energies (which are typically small relative to the orbital energy
      assuming that the unresolved halo masses are small compared to that of the primary progenitor. The mean energy of unresolved
      halos is then
      \begin{equation}
      \epsilon = \bar{E}_\mathrm{unres} \int_0^1 x^{1+a} \mathrm{d}x = \frac{\bar{E}_\mathrm{unres}}{2+a},
      \end{equation}
      while the variance in this energy is
      \begin{equation}
      \sigma^2_\epsilon = \sigma^2_\mathrm{unres} \int_0^1 x^{2+a} \mathrm{d}x = \frac{\sigma^2_\mathrm{unres}}{3+a}.
      \end{equation}
      The fractional variance is then
      \begin{equation}
      \frac{\sigma^2_\epsilon}{\epsilon^2} =\frac{2+a}{3+a} \left(\frac{\sigma_\mathrm{unres}}{\bar{E}_\mathrm{unres}}\right)^2.
      \end{equation}

      We therefore add to the energy of the halo an unresolved contribution of
      \begin{equation}
      \bar{E}_\mathrm{unres} \exp\left( \left[ \frac{2+a}{3+a} \left(\frac{\sigma_\mathrm{unres}}{\bar{E}_\mathrm{unres}}\right)^2 + (\sigma_e \log_\mathrm{e}10)^2 \right]^{1/2} \mathcal{N}(0,1) \right).
      \end{equation}
      where $\sigma_\mathrm{e}=${\normalfont \ttfamily [scatterExcess]} accounts for scatter missed by this model.

      The scale radius which corresponds to this energy is then solved for.

      For halos with mass less than $f_\mathrm{res} M_\mathrm{res}$, where $f_\mathrm{res}=${\normalfont \ttfamily
      [factorMassResolution]} and $M_\mathrm{res}$ is the mass resolution of the merger tree, and for any halo which has no
      progenitors (a leaf node), the scale radius is instead computed using an alternative method\footnote{For leaf nodes, there
      are no progenitors for which to apply the above energy calculation, and for halos sufficiently close to the mass resolution
      the energy calculation may not be reliable due to the poorly-resolved formation history of the node.}. In these cases, the
      entire extent of the branch for which this criterion applies is first determined. Each halo in this sub-branch is first
      assigned a scale radius from a fall-back \refClass{darkMatterHaloScaleClass} object which should be configured to return a
      \emph{scatter-free} scale radius for halos of given mass and redshift\footnote{As scatter will be added directly by the
      present class.} Then, a correlated set of random, log-normal deviates are applied to the scale radii of these nodes. That
      is, the scale radius of the $i^\mathrm{th}$ node in such a sub-branch will be $r_\mathrm{s} = \bar{r}_{\mathrm{s}, i} 10^{x_i}$
      where $x_i$ is a normally-distributed random variate with mean zero and dispersion $\sigma^\prime=${\normalfont \ttfamily
      [scatter]}. The deviates $x_i$ are assumed to be correlated with correlation matrix:
      \begin{equation}
       C_{i,j} = \exp\left( -\gamma \left| \log_{10} \frac{M_i}{M_j} \right|^\mu \right),
      \end{equation}

      where\footnote{These values were found by fitting to results from this class.} $\gamma = ${\normalfont \ttfamily
      [correlationRateDecay]}, $\mu =${\normalfont \ttfamily [correlationExponent]} and $M_i$ is the mass of the $i^\mathrm{th}$
      halo in the sub-branch. This results in a scale radius along the sub-branch with the correct mean and scatter, but
      correlated over mass increment scales in a way that matches the predictions of this algorithm.
    </description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusJohnson2021
     !!{
     A dark matter profile scale radius class that assigns dark matter profile scale radii using the energy-based model of
     \cite{johnson_random_2021}.
     !!}
     private
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_           => null()
     class           (criticalOverdensityClass         ), pointer :: criticalOverdensity_          => null()
     class           (cosmologicalMassVarianceClass    ), pointer :: cosmologicalMassVariance_     => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_          => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_         => null()
     class           (virialOrbitClass                 ), pointer :: virialOrbit_                  => null()
     class           (mergerTreeMassResolutionClass    ), pointer :: mergerTreeMassResolution_     => null()
     type            (branch                           ), pointer :: branches                      => null()
     double precision                                             :: massExponent                           , energyBoost              , &
          &                                                          unresolvedEnergy                       , factorMassResolution     , &
          &                                                          scatter                                , scatterExcess            , &
          &                                                          peakHeightExponent                     , correlationRateDecay     , &
          &                                                          correlationExponent
     integer                                                      :: countSampleEnergyUnresolved
     logical                                                      :: mainBranchOnly                         , applySubsamplingWeights  , &
          &                                                          acceptUnboundOrbits                    , includeUnresolvedVariance
   contains
     final     ::           darkMatterProfileScaleJohnson2021Destructor
     procedure :: radius => darkMatterProfileScaleJohnson2021Radius
  end type darkMatterProfileScaleRadiusJohnson2021
  
  interface darkMatterProfileScaleRadiusJohnson2021
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusJohnson2021} node operator class.
     !!}
     module procedure darkMatterProfileScaleJohnson2021ConstructorParameters
     module procedure darkMatterProfileScaleJohnson2021ConstructorInternal
  end interface darkMatterProfileScaleRadiusJohnson2021

  ! Sub-module-scope variables used in root finding.
  double precision                                                   :: energyTotal
  class           (darkMatterProfileScaleRadiusJohnson2021), pointer :: self_
  class           (nodeComponentDarkMatterProfile         ), pointer :: darkMatterProfile_
  type            (treeNode                               ), pointer :: node_
  !$omp threadprivate(energyTotal,self_,node_,darkMatterProfile_)

contains
  
  function darkMatterProfileScaleJohnson2021ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusJohnson2021} dark matter profile scale radius class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileScaleRadiusJohnson2021)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass               ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass          ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterProfileScaleRadiusClass      ), pointer       :: darkMatterProfileScaleRadius_
    class           (darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass              ), pointer       :: darkMatterProfileDMO_
    class           (virialOrbitClass                       ), pointer       :: virialOrbit_
    class           (mergerTreeMassResolutionClass          ), pointer       :: mergerTreeMassResolution_
    double precision                                                         :: massExponent                 , energyBoost           , &
         &                                                                      unresolvedEnergy             , factorMassResolution  , &
         &                                                                      scatter                      , scatterExcess         , &
         &                                                                      peakHeightExponent           , correlationRateDecay  , &
         &                                                                      correlationExponent
    integer                                                                  :: countSampleEnergyUnresolved
    logical                                                                  :: mainBranchOnly               , applySubsamplingWeights, &
         &                                                                      acceptUnboundOrbits          , includeUnresolvedVariance
    
    !![
    <inputParameter>
      <name>energyBoost</name>
      <defaultValue>0.797d0</defaultValue>
      <defaultSource>\citep{johnson_random_2021}</defaultSource>
      <source>parameters</source>
      <description>A boost to the energy.</description>
    </inputParameter>
    <inputParameter>
      <name>massExponent</name>
      <defaultValue>2.168d0</defaultValue>
      <defaultSource>\citep{johnson_random_2021}</defaultSource>
      <source>parameters</source>
      <description>The exponent of mass ratio appearing in the orbital energy term.</description>
    </inputParameter>
    <inputParameter>
      <name>peakHeightExponent</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The exponent of peak height in the orbital energy term.</description>
    </inputParameter>
    <inputParameter>
      <name>unresolvedEnergy</name>
      <defaultValue>0.550d0</defaultValue>
      <defaultSource>\citep{johnson_random_2021}</defaultSource>
      <source>parameters</source>
      <description>Factor multiplying the estimate of the internal energy of unresolved accretion.</description>
    </inputParameter>
    <inputParameter>
      <name>factorMassResolution</name>
      <defaultValue>1.0d2</defaultValue>
      <source>parameters</source>
      <description>The \cite{johnson_random_2021} model is applied only for halos with mass greater than $f M_\mathrm{res}$ where $f=${\normalfont \ttfamily [factorMassResolution]}. Below this mass the fall-back method is used, with correlated scatter along the branch.</description>
    </inputParameter>
    <inputParameter>
      <name>scatter</name>
      <defaultValue>0.16d0</defaultValue>
      <source>parameters</source>
      <description>The scatter in scale radius (in dex) to be applied to halos that are too low mass for the \cite{johnson_random_2021} model to be applied.</description>
    </inputParameter>
    <inputParameter>
      <name>scatterExcess</name>
      <defaultValue>0.116d0</defaultValue>
      <source>parameters</source>
      <description>The additional scatter radius (in dex) to be applied to the energies of merging halos to account for the fact that the \cite{johnson_random_2021} underpredicts the scatter in concentrations.</description>
    </inputParameter>
    <inputParameter>
      <name>correlationRateDecay</name>
      <defaultValue>3.8361d0</defaultValue>
      <source>parameters</source>
      <description>The decay rate for the exponential correlation model for scale radii.</description>
    </inputParameter>
    <inputParameter>
      <name>correlationExponent</name>
      <defaultValue>1.6198d0</defaultValue>
      <source>parameters</source>
      <description>The exponent for the exponential correlation model for scale radii.</description>
    </inputParameter>
    <inputParameter>
      <name>countSampleEnergyUnresolved</name>
      <defaultValue>100</defaultValue>
      <source>parameters</source>
      <description>The number of samples to use in Monte Carlo estimates of the mean and scatter in the energy of unresolved halos.</description>
    </inputParameter>
    <inputParameter>
      <name>mainBranchOnly</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, the \cite{johnson_random_2021} algorithm is applied to the main branch only, with other branches using the fall-back method.</description>
    </inputParameter>
    <inputParameter>
      <name>applySubsamplingWeights</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>If true, account for halo subsampling weights.</description>
    </inputParameter>
    <inputParameter>
      <name>acceptUnboundOrbits</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, allow unbound orbits when sampling orbits for unresolved halos.</description>
    </inputParameter>
    <inputParameter>
      <name>includeUnresolvedVariance</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, account for the variance in the energy of unresolved halos.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="criticalOverdensity"          name="criticalOverdensity_"          source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"     name="cosmologicalMassVariance_"     source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="virialOrbit"                  name="virialOrbit_"                  source="parameters"/>
    <objectBuilder class="mergerTreeMassResolution"     name="mergerTreeMassResolution_"     source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    !!]
    self=darkMatterProfileScaleRadiusJohnson2021(massExponent,peakHeightExponent,energyBoost,unresolvedEnergy,factorMassResolution,scatter,scatterExcess,correlationRateDecay,correlationExponent,countSampleEnergyUnresolved,mainBranchOnly,applySubsamplingWeights,acceptUnboundOrbits,includeUnresolvedVariance,cosmologyFunctions_,darkMatterProfileScaleRadius_,darkMatterHaloScale_,darkMatterProfileDMO_,virialOrbit_,mergerTreeMassResolution_,criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="criticalOverdensity_"         />
    <objectDestructor name="cosmologicalMassVariance_"    />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="virialOrbit_"                 />
    <objectDestructor name="mergerTreeMassResolution_"    />
    <objectDestructor name="darkMatterProfileDMO_"        />
    !!]
    return
  end function darkMatterProfileScaleJohnson2021ConstructorParameters

  function darkMatterProfileScaleJohnson2021ConstructorInternal(massExponent,peakHeightExponent,energyBoost,unresolvedEnergy,factorMassResolution,scatter,scatterExcess,correlationRateDecay,correlationExponent,countSampleEnergyUnresolved,mainBranchOnly,applySubsamplingWeights,acceptUnboundOrbits,includeUnresolvedVariance,cosmologyFunctions_,darkMatterProfileScaleRadius_,darkMatterHaloScale_,darkMatterProfileDMO_,virialOrbit_,mergerTreeMassResolution_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileScaleRadiusJohnson2021} dark matter profile scale radius class.
    !!}
    implicit none
    type            (darkMatterProfileScaleRadiusJohnson2021)                        :: self
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass               ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass          ), intent(in   ), target :: cosmologicalMassVariance_
    class           (darkMatterProfileScaleRadiusClass      ), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    class           (virialOrbitClass                       ), intent(in   ), target :: virialOrbit_
    class           (mergerTreeMassResolutionClass          ), intent(in   ), target :: mergerTreeMassResolution_
    class           (darkMatterProfileDMOClass              ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                         , intent(in   )         :: massExponent                 , energyBoost             , &
         &                                                                              unresolvedEnergy             , factorMassResolution    , &
         &                                                                              scatter                      , scatterExcess           , &
         &                                                                              peakHeightExponent           , correlationRateDecay    , &
         &                                                                              correlationExponent
    integer                                                  , intent(in   )         :: countSampleEnergyUnresolved
    logical                                                  , intent(in   )         :: mainBranchOnly               , applySubsamplingWeights  , &
         &                                                                              acceptUnboundOrbits          , includeUnresolvedVariance
    !![
    <constructorAssign variables="massExponent, peakHeightExponent, energyBoost, unresolvedEnergy, factorMassResolution, scatter, scatterExcess, correlationRateDecay, correlationExponent, countSampleEnergyUnresolved, mainBranchOnly, applySubsamplingWeights, acceptUnboundOrbits, includeUnresolvedVariance, *cosmologyFunctions_, *darkMatterProfileScaleRadius_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *virialOrbit_, *mergerTreeMassResolution_, *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function darkMatterProfileScaleJohnson2021ConstructorInternal

  subroutine darkMatterProfileScaleJohnson2021Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileScaleRadiusJohnson2021} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(darkMatterProfileScaleRadiusJohnson2021), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"          />
    <objectDestructor name="self%criticalOverdensity_"         />
    <objectDestructor name="self%cosmologicalMassVariance_"    />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    <objectDestructor name="self%darkMatterHaloScale_"         />
    <objectDestructor name="self%virialOrbit_"                 />
    <objectDestructor name="self%mergerTreeMassResolution_"    />
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    !!]
    return
  end subroutine darkMatterProfileScaleJohnson2021Destructor

  double precision function darkMatterProfileScaleJohnson2021Radius(self,node) result(radiusScale)
    !!{
    Initialize dark matter profile scale radii.
    !!}
    use :: Calculations_Resets             , only : Calculations_Reset
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , nodeComponentDarkMatterProfile, nodeComponentSatellite
    use :: Root_Finder                     , only : rootFinder                    , rangeExpandMultiplicative     , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    use :: Kepler_Orbits                   , only : keplerOrbit
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Hypergeometric_Functions        , only : Hypergeometric_2F1
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Linear_Algebra                  , only : matrix                        , assignment(=)
    implicit none
    class           (darkMatterProfileScaleRadiusJohnson2021), intent(inout) , target      :: self
    type            (treeNode                               ), intent(inout) , target      :: node
    type            (treeNode                               )                , pointer     :: nodeChild                                       , nodeSibling                     , &
         &                                                                                    nodeUnresolved                                  , nodeSiblingLast                 , &
         &                                                                                    nodeBranchTip
    class           (nodeComponentBasic                     )                , pointer     :: basicChild                                      , basicSibling                    , &
         &                                                                                    basic                                           , basicUnresolved
    class           (nodeComponentDarkMatterProfile         )                , pointer     :: darkMatterProfile                               , darkMatterProfileSibling        , &
         &                                                                                    darkMatterProfileChild                          , darkMatterProfileUnresolved
    class           (nodeComponentSatellite                 )                , pointer     :: satelliteSibling                                , satelliteUnresolved
    class           (massDistributionClass                  )                , pointer     :: massDistribution_
    type            (branch                                 )                , pointer     :: branch_                                         , branchLast
    double precision                                         , dimension(:  ), allocatable :: massLogarithmic                                 , deviates
    double precision                                         , dimension(:,:), allocatable :: covarianceNodes                                 , cholesky
    logical                                                  , dimension(:  ), allocatable :: independent
    double precision                                                         , parameter   :: massFunctionSlopeLogarithmic            =-1.80d0
    double precision                                                         , parameter   :: energyInternalFormFactorSlopeLogarithmic=-0.02d0
    double precision                                                         , parameter   :: scatterFractionalMaximum                =+5.00d0
    logical                                                                                :: haveBranch
    integer         (c_size_t                               )                              :: countNodes                                     , indexAlongBranch                 , &
         &                                                                                    countNodesUniqueMasses
    integer         (kind_int8                              )                              :: indexBranch                                    , i                                , &
         &                                                                                    j                                              , ii                               , &
         &                                                                                    jj
    double precision                                                                       :: energyOrbital                                  , massRatio                        , &
         &                                                                                    radiusScaleChild                               , radiusScaleOriginal              , &
         &                                                                                    massUnresolved                                 , radiusVirial                     , &
         &                                                                                    radiusScaleUnresolved                          , massResolution                   , &
         &                                                                                    energyOrbitalSubresolutionFactor               , energyInternalSubresolutionFactor, &
         &                                                                                    energyMean                                     , energyScatter                    , &
         &                                                                                    factorGrowthOrbital                            , factorGrowthInternal             , &
         &                                                                                    energySample                                   , peakHeight                       , &
         &                                                                                    energyVariance                                 , weightSubsampling
    type            (rootFinder                             )                              :: finder
    type            (keplerOrbit                            )                              :: orbit
    type            (matrix                                 )                              :: cholesky_
    
    ! Get the resolution of the tree.
    massResolution=self%mergerTreeMassResolution_%resolution(node%hostTree)
    ! Build a root finder to find scale radii.
    finder=rootFinder(                                                             &
         &            rootFunction                 =radiusScaleRoot              , &
         &            toleranceAbsolute            =1.0d-6                       , &
         &            toleranceRelative            =1.0d-3                       , &
         &            rangeExpandDownward          =0.5d+0                       , &
         &            rangeExpandUpward            =2.0d+0                       , &
         &            rangeExpandType              =rangeExpandMultiplicative    , &
         &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
         &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative  &
         &           )
    ! Create the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    ! If this node has no children, is below the mass at which we apply our energy model, or is not on the main branch and our
    ! method is being applied to the main branch only, set its dark matter profile scale radius using the fallback method plus
    ! correlated scatter along the branch.
    basic => node%basic()
    if     (                                                         &
         &   .not.associated(node%firstChild)                        &
         &  .or.                                                     &
         &   basic%mass() < self%factorMassResolution*massResolution &
         &  .or.                                                     &
         &   (                                                       &
         &          self%mainBranchOnly                              &
         &    .and.                                                  &
         &     .not.node%isOnMainBranch()                            &
         &   )                                                       &
         & ) then
       ! Find the unique ID identifying this branch (which we take to be the unique ID of the node at the branch tip), and also
       ! the position of this node along the branch.
       indexAlongBranch =  1_c_size_t
       nodeChild        => node
       do while (associated(nodeChild%firstChild))
          indexAlongBranch =  indexAlongBranch           +1_c_size_t
          nodeChild        => nodeChild       %firstChild
       end do
       nodeBranchTip => nodeChild
       indexBranch   =  nodeBranchTip%uniqueID()
       ! Walk up the branch to find its extent.
       countNodes =  indexAlongBranch-1_c_size_t
       nodeChild  => node
       basicChild => nodeChild%basic()
       do while (                                                                &
            &     associated(nodeChild)                                          &
            &    .and.                                                           &
            &     (                                                              &
            &       basicChild%mass() < self%factorMassResolution*massResolution &
            &      .or.                                                          &
            &       (self%mainBranchOnly .and. .not.node%isOnMainBranch())       &
            &     )                                                              &
            &   )
          countNodes=countNodes+1_c_size_t
          if (nodeChild%isPrimaryProgenitor()) then
             nodeChild  => nodeChild%parent
             basicChild => nodeChild%basic ()
          else
             nodeChild => null()
          end if
       end do
       if (countNodes <= 1_c_size_t) then
          ! With only a single node we can do simple sampling.
          radiusScale=self%darkMatterProfileScaleRadius_%radius(node)*10.0d0**(self%scatter*node%hostTree%randomNumberGenerator_%standardNormalSample())
       else
          ! Check if we have already computed radius offsets for this branch.
          branch_ => self%branches
          do while (associated(branch_) .and. branch_%uniqueID /= indexBranch)
             branch_ => branch_%next
          end do
          haveBranch=associated(branch_) .and. branch_%uniqueID == indexBranch
          if (.not.haveBranch) then
             ! We have not - compute offsets for the entire branch now. Begin by construct the covariance matrix.
             allocate(massLogarithmic(countNodes))
             allocate(independent    (countNodes))
             countNodesUniqueMasses =  0_c_size_t
             nodeChild              => nodeBranchTip
             do i=1,countNodes
                basicChild         =>       nodeChild %basic ()
                massLogarithmic(i) =  log10(basicChild%mass  ())
                nodeChild          =>       nodeChild %parent
                independent    (i) =                   i  == 1                    &
                     &                .or.                                        &
                     &                 massLogarithmic(i) >  massLogarithmic(i-1)
                if (independent(i)) countNodesUniqueMasses=+countNodesUniqueMasses &
                     &                                     +1_c_size_t
             end do
             allocate(deviates       (countNodesUniqueMasses                       ))
             allocate(covarianceNodes(countNodesUniqueMasses,countNodesUniqueMasses))
             allocate(cholesky       (countNodesUniqueMasses,countNodesUniqueMasses))
             ii=0_c_size_t
             do i=1,countNodes
                if (.not.independent(i)) cycle
                ii          =ii+1_c_size_t
                jj          =  +0_c_size_t
                deviates(ii)=node%hostTree%randomNumberGenerator_%standardNormalSample()
                do j=1,countNodes
                   if (.not.independent(j)) cycle
                   jj=jj+1_c_size_t
                   covarianceNodes(ii,jj)=+self%scatter**2                     &
                        &                 *exp(                                &
                        &                      -self%correlationRateDecay      &
                        &                      *abs(                           &
                        &                           +massLogarithmic(i)        &
                        &                           -massLogarithmic(j)        &
                        &                          )**self%correlationExponent &
                        &                     )
                end do
             end do
             ! Create a new branch to store this set of offsets.
             allocate(branch_)
             branch_%uniqueID =  indexBranch
             branch_%next     => null()
             allocate(branch_%radiusOffsetLogarithmic(countNodes))
             ! Generate random deviates for scale radii.
             cholesky_                      =covarianceNodes
             call cholesky_%choleskyDecomposition()
             cholesky                       =cholesky_
             branch_%radiusOffsetLogarithmic=0.0d0
             ii                             =0_c_size_t
             do i=1,countNodes
                if (independent(i)) ii=ii+1_c_size_t
                jj=0_c_size_t
                do j=1,i
                   if (.not.independent(j)) cycle
                   jj=jj+1_c_size_t
                   branch_%radiusOffsetLogarithmic(i)=+branch_ %radiusOffsetLogarithmic(i    ) &
                        &                             +cholesky                        (ii,jj) &
                        &                             *deviates                        (   jj)
                end do
             end do
             ! Store the branch.
             if (associated(self%branches)) then
                branchLast => self%branches
                do while (associated(branchLast%next))
                   branchLast => branchLast%next
                end do
                branchLast%next => branch_
             else
                self%branches => branch_
             end if
          end if
          ! Set the radius for this node.
          radiusScale=self%darkMatterProfileScaleRadius_%radius(node)*10.0d0**branch_%radiusOffsetLogarithmic(indexAlongBranch)
          ! Remove the branch if we have reached the end of it.
          if (indexAlongBranch == countNodes) then
             if (associated(self%branches,branch_)) then
                self%branches => branch_%next
             else
                branchLast => self%branches
                do while (.not.associated(branchLast%next,branch_))
                   branchLast => branchLast%next
                end do
                branchLast%next => branch_%next
             end if
             deallocate(branch_)
          end if
       end if
    else
       ! Use the energy model to compute the scale radius for this node.
       basic                   =>  node                                            %basic            (                                                            )
       nodeChild               =>  node                                            %firstChild
       darkMatterProfileChild  =>  nodeChild                                       %darkMatterProfile(                                                            )
       basicChild              =>  nodeChild                                       %basic            (                                                            )
       radiusScaleChild        =   darkMatterProfileChild                          %scale            (                                                            )
       massUnresolved          =  +basic                                           %mass             (                                                            ) &
            &                     -basicChild                                      %mass             (                                                            )
       peakHeight              =  +self                  %criticalOverdensity_     %value            (time=basicChild%time(),mass=basicChild%mass(),node=nodeChild) &
            &                     /self                  %cosmologicalMassVariance_%rootVariance     (time=basicChild%time(),mass=basicChild%mass()               )
       ! Iterate over progenitors and sum their energies.
       nodeSibling       =>                                                       nodeChild
       massDistribution_ => self             %darkMatterProfileDMO_ %get         (nodeSibling                   )
       radiusVirial      =  self             %darkMatterHaloScale_  %radiusVirial(nodeSibling                   )
       energyTotal       =  massDistribution_%energy                             (radiusVirial,massDistribution_)
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
       do while (associated(nodeSibling%sibling))
          nodeSibling              =>  nodeSibling       %sibling
          basicSibling             =>  nodeSibling       %basic                             (                              )
          darkMatterProfileSibling =>  nodeSibling       %darkMatterProfile                 (                              )
          satelliteSibling         =>  nodeSibling       %satellite                         (autoCreate=.true.             )
          massDistribution_        =>  self              %darkMatterProfileDMO_%get         (nodeSibling                   )
          radiusVirial             =   self              %darkMatterHaloScale_ %radiusVirial(nodeSibling                   )
          orbit                    =   satelliteSibling  %virialOrbit                       (                              )
          if (self%applySubsamplingWeights) then
             weightSubsampling     =   nodeSibling       %subsamplingWeight                 (                              )
          else
             weightSubsampling     =  +1.0d0
          end if
          massRatio                =  +basicSibling      %mass                              (                              ) &
               &                      /basicChild        %mass                              (                              )
          energyOrbital            =  +orbit             %energy                            (                              ) &
               &                      *basicSibling      %mass                              (                              ) &
               &                      /(                                                                                     &
               &                        +1.0d0                                                                               &
               &                        +massRatio                                                                           &
               &                       )**self%massExponent
          massUnresolved           =  +massUnresolved                                                                        &
               &                      -basicSibling      %mass                              (                              ) &
               &                      *weightSubsampling
          ! Add orbital and energy of this sibling.
          energyTotal              = +energyTotal                                                                            &
               &                     +(                                                                                      &
               &                       +energyOrbital                                                                        &
               &                       +massDistribution_%energy                            (radiusVirial,massDistribution_) &
               &                       /(                                                                                    &
               &                         +1.0d0                                                                              &
               &                         +massRatio                                                                          &
               &                        )**self%massExponent                                                                 &
               &                      )                                                                                      &
               &                     *(                                                                                      &
               &                       +1.0d0                                                                                &
               &                       +self%energyBoost                                                                     &
               &                      )                                                                                      &
               &                     *peakHeight**self%peakHeightExponent                                                    &
               &                     *weightSubsampling                                                                      &
               &                     *10.0d0**(                                                                              &
               &                               +self         %scatterExcess                                                  &
               &                               *node%hostTree%randomNumberGenerator_%standardNormalSample()                  &
               &                              )
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]
       end do
       ! Account for unresolved accretion. We assume that unresolved halos are accreted with the mean orbital energy of
       ! the virial orbital parameter distribution, plus an internal energy corresponding to that of a halo with mass
       ! equal to the total unresolved mass scaled by some correction factor (to account for the fact that the unresolved
       ! accretion will not in fact be in a single halo).
       if (massUnresolved > 0.0d0) then
          ! Construct a halo at the resolution limit.
          nodeUnresolved          => treeNode         ()
          nodeUnresolved%hostTree => node    %hostTree
          nodeUnresolved%parent   => node
          ! Find the last sibling in the children of our node. This allows us to attach our unresolved node into the tree
          ! structure temporarily.
          nodeSiblingLast => node%firstChild
          do while (associated(nodeSiblingLast%sibling))
             nodeSiblingLast => nodeSiblingLast%sibling
          end do
          nodeSiblingLast%sibling => nodeUnresolved
          ! Set fixed properties of our unresolved node.
          basicUnresolved             => nodeUnresolved%basic            (autoCreate=.true.)
          darkMatterProfileUnresolved => nodeUnresolved%darkMatterProfile(autoCreate=.true.)
          call basicUnresolved%massSet            (min(massResolution   ,massUnresolved))
          call basicUnresolved%timeSet            (    basicChild%time()                )
          call basicUnresolved%timeLastIsolatedSet(    basicChild%time()                )
          ! Compute a correction factor to the orbital energy which takes into account the mass dependence of the 1/(1+m/M)ᵅ term
          ! that is applied to the orbital energy. Averaging this over a power-law mass function gives the result below.
          massRatio                       =+basicUnresolved%mass() &
               &                           /basicChild     %mass()
          energyOrbitalSubresolutionFactor=+(1.0d0+massRatio)**self%massExponent                                       &
               &                           *Hypergeometric_2F1(                                                        &
               &                                               [2.0d0+massFunctionSlopeLogarithmic,self%massExponent], &
               &                                               [3.0d0+massFunctionSlopeLogarithmic                  ], &
               &                                               -massRatio                                              &
               &                                              )
          ! Compute a correction factor to the internal energy which takes into account the mass dependence of the 1/(1+m/M)ᵅ
          ! term that is applied to the orbital energy. Averaging this over a power-law mass function gives the result
          ! below.
          energyInternalSubresolutionFactor=+(1.0d0+massRatio)**self%massExponent                                                                                            &
               &                            *(2.0d0+massFunctionSlopeLogarithmic                                                     )                                       &
               &                            /(1.0d0+massFunctionSlopeLogarithmic+5.0d0/3.0d0+energyInternalFormFactorSlopeLogarithmic)                                       &
               &                            *Hypergeometric_2F1(                                                                                                             &
               &                                                [1.0d0+massFunctionSlopeLogarithmic+5.0d0/3.0d0+energyInternalFormFactorSlopeLogarithmic,self%massExponent], &
               &                                                [2.0d0+massFunctionSlopeLogarithmic+5.0d0/3.0d0+energyInternalFormFactorSlopeLogarithmic                  ], &
               &                                                -massRatio                                                                                                   &
               &                                               )
          ! Compute factors that account for the finite interval of time over which unresolved growth occurs.
          factorGrowthInternal=+0.5d0                                                                                                                                &
               &               *(                                                                                                                                    &
               &                 +1.0d0                                                                                                                              &
               &                 +1.0d0                                                                                                                              &
               &                 /(self%cosmologyFunctions_%expansionFactor(basic%time())/self%cosmologyFunctions_%expansionFactor(basicChild%time()))               &
               &               )
          factorGrowthOrbital =+0.5d0                                                                                                                                &
               &               *(                                                                                                                                    &
               &                 +1.0d0                                                                                                                              &
               &                 +(                                         basic%mass() /                                         basicChild%mass() )**(2.0d0/3.d0) &
               &                 /(self%cosmologyFunctions_%expansionFactor(basic%time())/self%cosmologyFunctions_%expansionFactor(basicChild%time()))               &
               &               )
          ! Get the virial radius of our unresolved node.
          radiusVirial=self%darkMatterHaloScale_%radiusVirial(nodeUnresolved)
          ! Monte Carlo sample orbits and internal structures for our unresolved halo. Use these to estimate the mean and scatter
          ! in the energy contributed by this halo.
          energyMean   =0.0d0
          energyScatter=0.0d0
          do j=1,self%countSampleEnergyUnresolved
             call Calculations_Reset(nodeUnresolved)
             ! Set internal structure of the unresolved halo.
             radiusScaleUnresolved=+          self         %darkMatterProfileScaleRadius_%radius              (nodeUnresolved) &
                  &                *10.0d0**(                                                                                  &
                  &                          +self                                       %scatter                              &
                  &                          *node%hostTree%randomNumberGenerator_       %standardNormalSample(              ) &
                  &                         )
             call darkMatterProfileUnresolved%scaleSet(radiusScaleUnresolved)
             massDistribution_ => self%darkMatterProfileDMO_%get(nodeUnresolved)
             ! Get an orbit for the unresolved halo.
             call nodeUnresolved%satelliteDestroy()
             satelliteUnresolved => nodeUnresolved     %satellite  (autoCreate=.true.)
             orbit               =  satelliteUnresolved%virialOrbit(                 )
             if (.not.orbit%isDefined()) orbit=self%virialOrbit_%orbit(nodeUnresolved,node%firstChild,self%acceptUnboundOrbits)
             ! Determine the orbital and internal energies.
             energySample=+massUnresolved                                                                                                                           &
                  &       *self%unresolvedEnergy                                                                                                                    &
                  &       *(                                                                                                                                        &
                  &         +orbit            %energy(                              )                       *energyOrbitalSubresolutionFactor *factorGrowthOrbital  &
                  &         +massDistribution_%energy(radiusVirial,massDistribution_)/basicUnresolved%mass()*energyInternalSubresolutionFactor*factorGrowthInternal &
                  &        )                                                                                                                                        &
                  &       /(                                                                                                                                        &
                  &         +1.0d0                                                                                                                                  &
                  &         +massRatio                                                                                                                              &
                  &        )**self%massExponent                                                                                                                     &
                  &       *(                                                                                                                                        &
                  &         +1.0d0                                                                                                                                  &
                  &         +self%energyBoost                                                                                                                       &
                  &        )                                                                                                                                        &
                  &       *peakHeight**self%peakHeightExponent
             !![
	     <objectDestructor name="massDistribution_"/>
	     !!]
             ! Accumulate to the mean and scatter in energy.
             energyMean   =+energyMean   +energySample
             energyScatter=+energyScatter+energySample**2
          end do
          ! Destroy the unresolved node and decouple it from the merger tree.
          call nodeUnresolved%destroy()
          deallocate(nodeUnresolved)
          nodeSiblingLast%sibling => null()
          ! Evaluate the mean and scatter in energy.
          energyMean    =+energyMean   /dble(self%countSampleEnergyUnresolved)
          energyVariance=+energyScatter/dble(self%countSampleEnergyUnresolved)-energyMean**2
          if (energyVariance > 0.0d0 .and. self%includeUnresolvedVariance) then
             energyScatter=sqrt(energyVariance)
          else
             energyScatter=0.0d0
          end if
          ! Accumulate the energy of unresolved halos to the total, including random scatter.
          energyTotal=+energyTotal                                                                 &
               &      +energyMean                                                                  &
               &      *exp(                                                                        &
               &           +node%hostTree%randomNumberGenerator_%standardNormalSample()            &
               &           *sqrt(                                                                  &
               &                 +(min(energyScatter/abs(energyMean),scatterFractionalMaximum**2)) &
               &                 *(2.0d0+massFunctionSlopeLogarithmic)                             &
               &                 /(3.0d0+massFunctionSlopeLogarithmic)                             &
               &                 +(self%scatterExcess*log(10.0d0))**2                              &
               &                )                                                                  &
               &          )
       end if
       ! Check for positive energy.
       if (energyTotal >= 0.0d0) then
          ! Energy is positive - the model has failed - simply assume no change in the scale radius.
          radiusScale         =  radiusScaleChild
       else
          ! Convert energy back to scale radius.
          self_               => self
          node_               => node
          darkMatterProfile_  => darkMatterProfile
          radiusScaleOriginal =  darkMatterProfile%scale(                          )
          radiusScale         =  finder           %find (rootGuess=radiusScaleChild)
       end if
       call darkMatterProfile%scaleSet(radiusScaleOriginal)
       call Calculations_Reset(node)
    end if
    return
  end function darkMatterProfileScaleJohnson2021Radius

  double precision function radiusScaleRoot(radiusScale)
    !!{
    Function used in root-finding to compute the scale radius of a dark matter profile as a given energy.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    use :: Mass_Distributions , only : massDistributionClass
    implicit none
    double precision                       , intent(in   ) :: radiusScale
    class           (massDistributionClass), pointer       :: massDistribution_
    
    call darkMatterProfile_%scaleSet(radiusScale)
    call Calculations_Reset(node_)
    massDistribution_ =>  self_%darkMatterProfileDMO_%get        (                                        node_                   )
    radiusScaleRoot   =   +                           energyTotal                                                                   &
         &                -     massDistribution_    %energy     (self_%darkMatterHaloScale_%radiusVirial(node_),massDistribution_)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function radiusScaleRoot
