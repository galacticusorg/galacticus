!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  An implementation of a dark matter density profile which includes the accretion flow surrounding the halo.
  !!}

  use :: Cosmology_Functions         , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters        , only : cosmologyParametersClass
  use :: Cosmological_Density_Field  , only : cosmologicalMassVarianceClass      , criticalOverdensityClass
  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOAccretionFlow">
    <description>
      An implementation of a dark matter density profile which includes the accretion flow surrounding the halo. The density
      profile is modelled as
      \begin{equation}
      \rho(r) = f_\mathrm{trans}(r) \rho_\mathrm{halo}(r) + \rho_\mathrm{accretion}(r),
      \end{equation}
      where $\rho_\mathrm{halo}(r)$ is the halo density profile (provided by a \refClass{darkMatterProfileDMOClass} object) and
      $\rho_\mathrm{accretion}(r)$ is the accretion flow density profile (provided by a \refClass{accretionFlowsClass} object),
      and $f_\mathrm{trans}(r)$ is the transition function as defined by equation~(7) of \cite{diemer_dependence_2014}.

      Note that some \refClass{accretionFlowsClass} objects make use of an \refClass{darkMatterProfileClass} object. As such, to
      avoid an infinite recursive loop, it is recommended to use an explicit, separate \refClass{darkMatterProfileClass} for the
      relevant \refClass{accretionFlowsClass} object when using this implementation.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOAccretionFlow
     !!{
     A dark matter halo profile class which implements a dark matter density profile which includes the accretion flow surrounding the halo.
     !!}
     private
     class(cosmologyParametersClass     ), pointer :: cosmologyParameters_      => null()
     class(cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class(*                            ), pointer :: accretionFlows_           => null()
     class(darkMatterProfileDMOClass    ), pointer :: darkMatterProfileDMO_     => null()
   contains
     final     ::                                      accretionFlowDestructor
     procedure :: density                           => accretionFlowDensity
     procedure :: densityLogSlope                   => accretionFlowDensityLogSlope
     procedure :: radiusEnclosingDensity            => accretionFlowRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => accretionFlowRadiusEnclosingMass
     procedure :: radialMoment                      => accretionFlowRadialMoment
     procedure :: enclosedMass                      => accretionFlowEnclosedMass
     procedure :: potential                         => accretionFlowPotential
     procedure :: circularVelocity                  => accretionFlowCircularVelocity
     procedure :: radiusCircularVelocityMaximum     => accretionFlowRadiusCircularVelocityMaximum
     procedure :: circularVelocityMaximum           => accretionFlowCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => accretionFlowRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => accretionFlowRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => accretionFlowRotationNormalization
     procedure :: energy                            => accretionFlowEnergy
     procedure :: kSpace                            => accretionFlowKSpace
     procedure :: freefallRadius                    => accretionFlowFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => accretionFlowFreefallRadiusIncreaseRate
     procedure :: deepCopy                          => accretionFlowDeepCopy
     procedure :: deepCopyReset                     => accretionFlowDeepCopyReset
     procedure :: deepCopyFinalize                  => accretionFlowDeepCopyFinalize
  end type darkMatterProfileDMOAccretionFlow

  interface darkMatterProfileDMOAccretionFlow
     !!{
     Constructors for the {\normalfont \ttfamily accretionFlow} dark matter halo profile class.
     !!}
     module procedure accretionFlowConstructorParameters
     module procedure accretionFlowConstructorInternal
  end interface darkMatterProfileDMOAccretionFlow

contains

  function accretionFlowConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily accretionFlow} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter          , inputParameters
    use :: Functions_Global, only : accretionFlowsConstruct_, accretionFlowsDestruct_
    implicit none
    type            (darkMatterProfileDMOAccretionFlow)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (*                                ), pointer       :: accretionFlows_
    class           (darkMatterProfileDMOClass        ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
    class           (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass         ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass    ), pointer       :: cosmologicalMassVariance_
    double precision                                                   :: toleranceRelativePotential

    !![
    <inputParameter>
      <name>toleranceRelativePotential</name>
      <defaultValue>1.0d-6</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the potential in dark-matter-only density profiles.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"     name="darkMatterProfileDMO_"     source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    !!]
    call accretionFlowsConstruct_(parameters,accretionFlows_)
    self=darkMatterProfileDMOAccretionFlow(toleranceRelativePotential,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,accretionFlows_,darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters" extraAllowedNames="accretionFlows"/>
    <objectDestructor name="darkMatterHaloScale_"     />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="darkMatterProfileDMO_"    />
    <objectDestructor name="darkMatterHaloScale_"     />
    !!]
    if (associated(accretionFlows_)) call accretionFlowsDestruct_(accretionFlows_)
    return
  end function accretionFlowConstructorParameters

  function accretionFlowConstructorInternal(toleranceRelativePotential,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,accretionFlows_,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily accretionFlow} dark matter profile class.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (darkMatterProfileDMOAccretionFlow)                        :: self
    class           (*                                ), intent(in   ), target :: accretionFlows_
    class           (darkMatterProfileDMOClass        ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyParametersClass         ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass         ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass    ), intent(in   ), target :: cosmologicalMassVariance_
    double precision                                   , intent(in   )         :: toleranceRelativePotential
    !![
    <constructorAssign variables="toleranceRelativePotential, *cosmologyParameters_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *accretionFlows_, *darkMatterProfileDMO_, *darkMatterHaloScale_"/>
    !!]

    return
  end function accretionFlowConstructorInternal

  subroutine accretionFlowDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily accretionFlow} dark matter halo profile class.
    !!}
    use :: Functions_Global, only : accretionFlowsDestruct_
    implicit none
    type(darkMatterProfileDMOAccretionFlow), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"    />
    <objectDestructor name="self%darkMatterHaloScale_"     />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    if (associated(self%accretionFlows_)) call accretionFlowsDestruct_(self%accretionFlows_)
    return
  end subroutine accretionFlowDestructor

  double precision function accretionFlowDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Functions_Global, only : accretionFlowsDensity_
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , intent(in   ) :: radius
    class           (nodeComponentBasic               ), pointer       :: basic
    double precision                                                   :: radius200Mean   , density200Mean    , &
         &                                                                radiusTransition, fractionTransition, &
         &                                                                peakHeight

    if (node%isSatellite()) then
       ! No accretion flow for satellites.
       accretionFlowDensity =  +self%darkMatterProfileDMO_%density(node,radius) 
    else
       ! Include accretion flow for non-satellites
       basic                =>  node%basic()
       peakHeight           =  +self%criticalOverdensity_     %value                 (time=basic%time(),mass=basic%mass()  ) &
            &                  /self%cosmologicalMassVariance_%rootVariance          (time=basic%time(),mass=basic%mass()  )
       density200Mean       =  +200.0d0                                                                                      &
            &                  *self%cosmologyFunctions_      %matterDensityEpochal  (time=basic%time()                    )
       radius200Mean        =  +self%darkMatterProfileDMO_    %radiusEnclosingDensity(     node        ,     density200Mean)
       radiusTransition     =  +(             &
            &                    +1.90d0      &
            &                    -0.18d0      &
            &                    *peakHeight  &
            &                   )             &
            &                  *radius200Mean
       fractionTransition   =  +1.0d0                &
            &                  /(                    &
            &                    +1.0d0              &
            &                    +(                  &
            &                      +radius           &
            &                      /radiusTransition &
            &                     )**4               &
            &                   )**2
       accretionFlowDensity =  +self%darkMatterProfileDMO_%density           (                     node,radius) &
            &                  *                           fractionTransition                                   &
            &                  +accretionFlowsDensity_                       (self%accretionFlows_,node,radius)
    end if
    return
  end function accretionFlowDensity

  double precision function accretionFlowDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily
    node} at the given {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , intent(in   ) :: radius
    
    if (node%isSatellite()) then
       accretionFlowDensityLogSlope=self%darkMatterProfileDMO_%densityLogSlope         (node,radius)
    else
       accretionFlowDensityLogSlope=self                      %densityLogSlopeNumerical(node,radius)
    end if
    return
  end function accretionFlowDensityLogSlope

  double precision function accretionFlowRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout), target :: self
    type            (treeNode                         ), intent(inout), target :: node
    double precision                                   , intent(in   )         :: density

    if (node%isSatellite()) then
       accretionFlowRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity         (node,density)
    else
       accretionFlowRadiusEnclosingDensity=self                      %radiusEnclosingDensityNumerical(node,density)
    end if
    return
  end function accretionFlowRadiusEnclosingDensity

  double precision function accretionFlowRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout), target :: self
    type            (treeNode                         ), intent(inout), target :: node
    double precision                                   , intent(in   )         :: mass

    if (node%isSatellite()) then
       accretionFlowRadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass         (node,mass)
    else
       accretionFlowRadiusEnclosingMass=self                      %radiusEnclosingMassNumerical(node,mass)
    end if
    return
  end function accretionFlowRadiusEnclosingMass

  double precision function accretionFlowRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the radial moment of the density profile.
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout)           :: self
    type            (treeNode                         ), intent(inout)           :: node
    double precision                                   , intent(in   )           :: moment
    double precision                                   , intent(in   ), optional :: radiusMinimum, radiusMaximum
  
    if (node%isSatellite()) then
       accretionFlowRadialMoment=self%darkMatterProfileDMO_%radialMoment         (node,moment,radiusMinimum,radiusMaximum)
    else
       accretionFlowRadialMoment=self                      %radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    end if
    return
  end function accretionFlowRadialMoment

  double precision function accretionFlowEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , intent(in   ) :: radius

    if (node%isSatellite()) then
       accretionFlowEnclosedMass=self%darkMatterProfileDMO_%enclosedMass         (node,radius)
    else
       accretionFlowEnclosedMass=self                      %enclosedMassNumerical(node,radius)
    end if
    return
  end function accretionFlowEnclosedMass

  double precision function accretionFlowPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout)           :: self
    type            (treeNode                         ), intent(inout), target   :: node
    double precision                                   , intent(in   )           :: radius
    integer                                            , intent(  out), optional :: status

    if (node%isSatellite()) then
       accretionFlowPotential=self%darkMatterProfileDMO_%potential         (node,radius,status)
    else
       accretionFlowPotential=self                      %potentialNumerical(node,radius,status)
    end if
    return
  end function accretionFlowPotential

  double precision function accretionFlowCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , intent(in   ) :: radius

    if (radius > 0.0d0) then
       accretionFlowCircularVelocity=sqrt(                                 &
            &                             +gravitationalConstantGalacticus &
            &                             *self%enclosedMass(node,radius)  &
            &                             /                       radius   &
            &                            )
    else
       accretionFlowCircularVelocity=0.0d0
    end if
    return
  end function accretionFlowCircularVelocity

  double precision function accretionFlowRadiusCircularVelocityMaximum(self,node)
    !!{
    Returns the radius (in Mpc) at which the maximum circular velocity is acheived in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node

    if (node%isSatellite()) then
       accretionFlowRadiusCircularVelocityMaximum=self%darkMatterProfileDMO_%radiusCircularVelocityMaximum         (node)
    else
       accretionFlowRadiusCircularVelocityMaximum=self                      %radiusCircularVelocityMaximumNumerical(node)
    end if
    return
  end function accretionFlowRadiusCircularVelocityMaximum

  double precision function accretionFlowCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node

    if (node%isSatellite()) then
       accretionFlowCircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum         (node)
    else
       accretionFlowCircularVelocityMaximum=self                      %circularVelocityMaximumNumerical(node)
    end if
    return
  end function accretionFlowCircularVelocityMaximum

  double precision function accretionFlowRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , intent(in   ) :: radius
 
    if (node%isSatellite()) then
       accretionFlowRadialVelocityDispersion=self%darkMatterProfileDMO_%radialVelocityDispersion         (node,radius)
    else
       accretionFlowRadialVelocityDispersion=self                      %radialVelocityDispersionNumerical(node,radius)
    end if
    return
  end function accretionFlowRadialVelocityDispersion

  double precision function accretionFlowRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , intent(in   ) :: specificAngularMomentum

    if (node%isSatellite()) then
       accretionFlowRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum         (node,specificAngularMomentum)
    else
       accretionFlowRadiusFromSpecificAngularMomentum=self                      %radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    end if
    return
  end function accretionFlowRadiusFromSpecificAngularMomentum

  double precision function accretionFlowRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                                   :: radiusVirial

    radiusVirial                      =+self%darkMatterHaloScale_%radiusVirial(node                                                           )
    accretionFlowRotationNormalization=+self                     %radialMoment(node,moment=2.0d0,radiusMinimum=0.0d0,radiusMaximum=radiusVirial) &
         &                             /self                     %radialMoment(node,moment=3.0d0,radiusMinimum=0.0d0,radiusMaximum=radiusVirial)
    return
  end function accretionFlowRotationNormalization

  double precision function accretionFlowEnergy(self,node)
    !!{
    Return the energy of a accretionFlow halo density profile.
    !!}
    implicit none
    class(darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node

    if (node%isSatellite()) then
       accretionFlowEnergy=self%darkMatterProfileDMO_%energy         (node)
    else
       accretionFlowEnergy=self                      %energyNumerical(node)
    end if
    return
  end function accretionFlowEnergy

  double precision function accretionFlowKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the accretionFlow density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$), using the expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout)         :: self
    type            (treeNode                         ), intent(inout), target :: node
    double precision                                   , intent(in   )         :: waveNumber

    if (node%isSatellite()) then
       accretionFlowKSpace=self%darkMatterProfileDMO_%kSpace         (node,waveNumber)
    else
       accretionFlowKSpace=self                      %kSpaceNumerical(node,waveNumber)
    end if
    return
  end function accretionFlowKSpace

  double precision function accretionFlowFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the accretionFlow density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout), target :: self
    type            (treeNode                         ), intent(inout), target :: node
    double precision                                   , intent(in   )         :: time

    if (node%isSatellite()) then
       accretionFlowFreefallRadius=self%darkMatterProfileDMO_%freefallRadius         (node,time)
    else
       accretionFlowFreefallRadius=self                      %freefallRadiusNumerical(node,time)
    end if
    return
  end function accretionFlowFreefallRadius

  double precision function accretionFlowFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the accretionFlow density profile at the specified {\normalfont
    \ttfamily time} (given in Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOAccretionFlow), intent(inout), target :: self
    type            (treeNode                         ), intent(inout), target :: node
    double precision                                   , intent(in   )         :: time

    if (node%isSatellite()) then
       accretionFlowFreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate         (node,time)
    else
       accretionFlowFreefallRadiusIncreaseRate=self                      %freefallRadiusIncreaseRateNumerical(node,time)
    end if
    return
  end function accretionFlowFreefallRadiusIncreaseRate

  subroutine accretionFlowDeepCopyReset(self)
    !!{
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : accretionFlowsDeepCopyReset_
    implicit none
    class(darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    
    self%copiedSelf => null()
    if (associated(self%cosmologyParameters_     )) call self%cosmologyParameters_     %deepCopyReset()
    if (associated(self%cosmologyFunctions_      )) call self%cosmologyFunctions_      %deepCopyReset()
    if (associated(self%criticalOverdensity_     )) call self%criticalOverdensity_     %deepCopyReset()
    if (associated(self%cosmologicalMassVariance_)) call self%cosmologicalMassVariance_%deepCopyReset()
    if (associated(self%darkMatterHaloScale_     )) call self%darkMatterHaloScale_     %deepCopyReset() 
    if (associated(self%darkMatterProfileDMO_    )) call self%darkMatterProfileDMO_    %deepCopyReset()
    if (associated(self%accretionFlows_          )) call accretionFlowsDeepCopyReset_(self%accretionFlows_)
    return
  end subroutine accretionFlowDeepCopyReset
  
  subroutine accretionFlowDeepCopyFinalize(self)
    !!{
    Finalize a deep reset of the object.
    !!}
    use :: Functions_Global, only : accretionFlowsDeepCopyFinalize_
    implicit none
    class(darkMatterProfileDMOAccretionFlow), intent(inout) :: self

    if (associated(self%cosmologyParameters_     )) call self%cosmologyParameters_     %deepCopyFinalize()
    if (associated(self%cosmologyFunctions_      )) call self%cosmologyFunctions_      %deepCopyFinalize()
    if (associated(self%criticalOverdensity_     )) call self%criticalOverdensity_     %deepCopyFinalize()
    if (associated(self%cosmologicalMassVariance_)) call self%cosmologicalMassVariance_%deepCopyFinalize()
    if (associated(self%darkMatterHaloScale_     )) call self%darkMatterHaloScale_     %deepCopyFinalize()
    if (associated(self%darkMatterProfileDMO_    )) call self%darkMatterProfileDMO_    %deepCopyFinalize()
    if (associated(self%accretionFlows_          )) call accretionFlowsDeepCopyFinalize_(self%accretionFlows_)
    return
  end subroutine accretionFlowDeepCopyFinalize
  
  subroutine accretionFlowDeepCopy(self,destination)
    !!{
    Perform a deep copy of the object.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Functions_Global, only : accretionFlowsDeepCopy_
    implicit none
    class(darkMatterProfileDMOAccretionFlow), intent(inout) :: self
    class(darkMatterProfileDMOClass        ), intent(inout) :: destination

    call self%darkMatterProfileDMOClass%deepCopy(destination)
    select type (destination)
    type is (darkMatterProfileDMOAccretionFlow)
       nullify(destination%cosmologyParameters_)
       if (associated(self%cosmologyParameters_)) then
          if (associated(self%cosmologyParameters_%copiedSelf)) then
             select type(s => self%cosmologyParameters_%copiedSelf)
                class is (cosmologyParametersClass)
                destination%cosmologyParameters_ => s
                class default
                call Galacticus_Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%cosmologyParameters_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%cosmologyParameters_,mold=self%cosmologyParameters_)
             call self%cosmologyParameters_%deepCopy(destination%cosmologyParameters_)
             self%cosmologyParameters_%copiedSelf => destination%cosmologyParameters_
             call destination%cosmologyParameters_%autoHook()
          end if
       end if
       nullify(destination%cosmologyFunctions_)
       if (associated(self%cosmologyFunctions_)) then
          if (associated(self%cosmologyFunctions_%copiedSelf)) then
             select type(s => self%cosmologyFunctions_%copiedSelf)
                class is (cosmologyFunctionsClass)
                destination%cosmologyFunctions_ => s
                class default
                call Galacticus_Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%cosmologyFunctions_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%cosmologyFunctions_,mold=self%cosmologyFunctions_)
             call self%cosmologyFunctions_%deepCopy(destination%cosmologyFunctions_)
             self%cosmologyFunctions_%copiedSelf => destination%cosmologyFunctions_
             call destination%cosmologyFunctions_%autoHook()
          end if
       end if
       nullify(destination%criticalOverdensity_)
       if (associated(self%criticalOverdensity_)) then
          if (associated(self%criticalOverdensity_%copiedSelf)) then
             select type(s => self%criticalOverdensity_%copiedSelf)
                class is (criticalOverdensityClass)
                destination%criticalOverdensity_ => s
                class default
                call Galacticus_Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%criticalOverdensity_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%criticalOverdensity_,mold=self%criticalOverdensity_)
             call self%criticalOverdensity_%deepCopy(destination%criticalOverdensity_)
             self%criticalOverdensity_%copiedSelf => destination%criticalOverdensity_
             call destination%criticalOverdensity_%autoHook()
          end if
       end if
       nullify(destination%cosmologicalMassVariance_)
       if (associated(self%cosmologicalMassVariance_)) then
          if (associated(self%cosmologicalMassVariance_%copiedSelf)) then
             select type(s => self%cosmologicalMassVariance_%copiedSelf)
                class is (cosmologicalMassVarianceClass)
                destination%cosmologicalMassVariance_ => s
                class default
                call Galacticus_Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%cosmologicalMassVariance_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%cosmologicalMassVariance_,mold=self%cosmologicalMassVariance_)
             call self%cosmologicalMassVariance_%deepCopy(destination%cosmologicalMassVariance_)
             self%cosmologicalMassVariance_%copiedSelf => destination%cosmologicalMassVariance_
             call destination%cosmologicalMassVariance_%autoHook()
          end if
       end if
       nullify(destination%darkMatterHaloScale_)
       if (associated(self%darkMatterHaloScale_)) then
          if (associated(self%darkMatterHaloScale_%copiedSelf)) then
             select type(s => self%darkMatterHaloScale_%copiedSelf)
                class is (darkMatterHaloScaleClass)
                destination%darkMatterHaloScale_ => s
                class default
                call Galacticus_Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%darkMatterHaloScale_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%darkMatterHaloScale_,mold=self%darkMatterHaloScale_)
             call self%darkMatterHaloScale_%deepCopy(destination%darkMatterHaloScale_)
             self%darkMatterHaloScale_%copiedSelf => destination%darkMatterHaloScale_
             call destination%darkMatterHaloScale_%autoHook()
          end if
       end if
       nullify(destination%darkMatterProfileDMO_)
       if (associated(self%darkMatterProfileDMO_)) then
          if (associated(self%darkMatterProfileDMO_%copiedSelf)) then
             select type(s => self%darkMatterProfileDMO_%copiedSelf)
                class is (darkMatterProfileDMOClass)
                destination%darkMatterProfileDMO_ => s
                class default
                call Galacticus_Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%darkMatterProfileDMO_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%darkMatterProfileDMO_,mold=self%darkMatterProfileDMO_)
             call self%darkMatterProfileDMO_%deepCopy(destination%darkMatterProfileDMO_)
             self%darkMatterProfileDMO_%copiedSelf => destination%darkMatterProfileDMO_
             call destination%darkMatterProfileDMO_%autoHook()
          end if
       end if
       nullify(destination%accretionFlows_)
       if (associated(self%accretionFlows_)) then
          allocate(destination%accretionFlows_,mold=self%accretionFlows_)
          call accretionFlowsDeepCopy_(self%accretionFlows_,destination%accretionFlows_)
       end if
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine accretionFlowDeepCopy
