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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !!{
  Contains a class which implements the tidal heating rate model of \cite{gnedin_tidal_1999}.
  !!}

  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <satelliteTidalHeatingRate name="satelliteTidalHeatingRateGnedin1999">
   <description>
    A satellite tidal heating rate class which uses the formalism of \cite{gnedin_tidal_1999} to compute the heating rate:
    \begin{equation}
    \dot{Q}_\mathrm{tidal}=\frac{1}{3}\epsilon\left[1+\left(\frac{T_\mathrm{shock}}{T_\mathrm{orb}}\right)^2\right]^{-\gamma}
    g_{ij} G^{ij}
    \end{equation}
    where $T_\mathrm{orb}$ and $T_\mathrm{shock}$ are the orbital period and shock duration, respectively, of the satellite,
    $\epsilon=${\normalfont \ttfamily [epsilon]} and $\gamma=${\normalfont \ttfamily [gamma]} are model parameters, $g_{ij}$ is
    the tidal tensor, and $G_{ij}$ is the integral with respect to time of $g_{ij}$ along the orbit of the satellite.  Upon
    tidal heating, a mass element at radius $r_\mathrm{i}$ expands to radius $r_\mathrm{f}$, according to the equation
    \begin{equation}
    \frac{1}{r_\mathrm{f}}=\frac{1}{r_\mathrm{i}}-\frac{2r_\mathrm{i}^3Q_\mathrm{tidal}}{\mathrm{G}M_\mathrm{sat}(&lt;r_\mathrm{i})}.
    \end{equation}
   </description>
  </satelliteTidalHeatingRate>
  !!]
  type, extends(satelliteTidalHeatingRateClass) :: satelliteTidalHeatingRateGnedin1999
     !!{
     A satellite tidal heating rate class which implements the tidal heating rate model of \cite{gnedin_tidal_1999}.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: epsilon                       , gamma
   contains
     final     ::                gnedin1999Destructor
     procedure :: heatingRate => gnedin1999HeatingRate
  end type satelliteTidalHeatingRateGnedin1999

  interface satelliteTidalHeatingRateGnedin1999
     !!{
     Constructors for the gnedin1999 satellite tidal heating rate class.
     !!}
     module procedure gnedin1999ConstructorParameters
     module procedure gnedin1999ConstructorInternal
  end interface satelliteTidalHeatingRateGnedin1999

contains

  function gnedin1999ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily gnedin1999} satellite tidal heating rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteTidalHeatingRateGnedin1999)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyParametersClass           ), pointer       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    double precision                                                     :: epsilon             , gamma

    !![
    <inputParameter>
      <name>epsilon</name>
      <defaultValue>3.0d0</defaultValue>
      <description>Parameter, $\epsilon$, controlling the tidal heating rate of satellites in the {\normalfont \ttfamily Gnedin1999} method.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <defaultValue>2.5d0</defaultValue>
      <description>Parameter, $\gamma$, controlling the tidal heating rate of satellites in the {\normalfont \ttfamily Gnedin1999} method.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=satelliteTidalHeatingRateGnedin1999(epsilon,gamma,cosmologyParameters_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function gnedin1999ConstructorParameters

  function gnedin1999ConstructorInternal(epsilon,gamma,cosmologyParameters_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily gnedin1999} satellite tidal heating rate class.
    !!}
    implicit none
    type            (satelliteTidalHeatingRateGnedin1999)                        :: self
    class           (cosmologyParametersClass           ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                     , intent(in)            :: epsilon             , gamma
    !![
    <constructorAssign variables="epsilon, gamma, *cosmologyParameters_, *darkMatterHaloScale_"/>
    !!]

    return
  end function gnedin1999ConstructorInternal

  subroutine gnedin1999Destructor(self)
    !!{
    Default constructor for the {\normalfont \ttfamily gnedin1999} satellite tidal heating rate class.
    !!}
    implicit none
    type(satelliteTidalHeatingRateGnedin1999), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine gnedin1999Destructor

  double precision function gnedin1999HeatingRate(self,node)
    !!{
    Return the tidal heating rate for satellite halos assuming the model of \cite{gnedin_tidal_1999}.
    !!}
    use :: Galactic_Structure_Options      , only : componentTypeAll     , coordinateSystemCartesian, massTypeDark
    use :: Galacticus_Nodes                , only : nodeComponentBasic   , nodeComponentSatellite   , treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear             , megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Tensors                         , only : assignment(=)        , max                      , operator(*) , tensorRank2Dimension3Symmetric
    use :: Vectors                         , only : Vector_Magnitude
    use :: Coordinates                     , only : coordinateCartesian  , assignment(=)
    implicit none
    class           (satelliteTidalHeatingRateGnedin1999), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    class           (nodeComponentSatellite             ), pointer       :: satellite
    type            (treeNode                           ), pointer       :: nodeHost
    class           (nodeComponentBasic                 ), pointer       :: basic
    class           (massDistributionClass              ), pointer       :: massDistribution_        , massDistributionHost_
    double precision                                     , dimension(3)  :: velocity
    double precision                                     , parameter     :: radiusHalfMassSatelliteTiny=1.0d-30, fractionMassTiny         =1.0d-6
    type            (coordinateCartesian                )                :: position
    double precision                                                     :: massSatellite            , velocityCircularSatellite, &
         &                                                                  radius                   , speed                    , &
         &                                                                  timescaleShock           , heatingRateNormalized    , &
         &                                                                  orbitalFrequencySatellite, radiusHalfMassSatellite  , &
         &                                                                  massHalfSatellite        , fractionDarkMatter
    logical                                                              :: useFrequencyOrbital
    type            (tensorRank2Dimension3Symmetric     )                :: tidalTensor              , tidalTensorPathIntegrated
    
    ! Construct required properties of satellite and host.
    nodeHost                  => node     %mergesWith               (        )
    satellite                 => node     %satellite                (        )
    massSatellite             =  satellite%boundMass                (        )
    position                  =  satellite%position                 (        )
    velocity                  =  satellite%velocity                 (        )
    tidalTensorPathIntegrated =  satellite%tidalTensorPathIntegrated(        )
    radius                    =  position %rSpherical               (        )
    speed                     =  Vector_Magnitude                   (velocity)
    ! Find the universal dark matter fraction.
    fractionDarkMatter        =  +(                                         &
         &                         +self%cosmologyParameters_%OmegaMatter() &
         &                         -self%cosmologyParameters_%OmegaBaryon() &
         &                        )                                         &
         &                       /  self%cosmologyParameters_%OmegaMatter()
    ! Find the gravitational tidal tensor.
    massDistributionHost_    => nodeHost%massDistribution()
    tidalTensor              =  massDistributionHost_%tidalTensor(position)
    !![
    <objectDestructor name="massDistributionHost_"/>
    !!]
    ! Find the orbital frequency at the half mass radius of the satellite.
    massDistribution_        => node%massDistribution(componentTypeAll,massTypeDark)
    basic                    => node%basic()
    massHalfSatellite        =  +0.50d0                                                                                                       &
         &                      *min(                                                                                                         &
         &                           +                  fractionDarkMatter                                                                    &
         &                           *                  massSatellite                                                                       , &
         &                           +massDistribution_%massEnclosedBySphere(radius=self%darkMatterHaloScale_%radiusVirial           (node))  &
         &                      )
    useFrequencyOrbital      = massHalfSatellite > fractionMassTiny*basic%mass()
    if (useFrequencyOrbital) then
       radiusHalfMassSatellite  =     massDistribution_%radiusEnclosingMass (mass  =                                massHalfSatellite      )
       velocityCircularSatellite=     massDistribution_%rotationCurve       (radius=                          radiusHalfMassSatellite      )
       ! Compute the orbital frequency.
       orbitalFrequencySatellite =  +velocityCircularSatellite &
            &                       /radiusHalfMassSatellite   &
            &                       *gigaYear                  &
            &                       *kilo                      &
            &                       /megaParsec
    else
       ! No well-defined half mass radius exists in the satellite. Use the virial orbital frequency instead.
       orbitalFrequencySatellite =  +self%darkMatterHaloScale_%velocityVirial(node) &
            &                       /self%darkMatterHaloScale_%radiusVirial  (node) &
            &                       *gigaYear                                       &
            &                       *kilo                                           &
            &                       /megaParsec
    end if
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Catch non-positive speeds.
    if (speed > 0.0d0) then
       ! Find the shock timescale (i.e. crossing time in the radial direction).
       timescaleShock=+megaParsec &
            &         /kilo       &
            &         /gigaYear   &
            &         *radius     &
            &         /speed
       ! Compute the heating rate.
       heatingRateNormalized=+self%epsilon                                          &
            &                /(                                                     &
            &                  +1.0d0                                               &
            &                  +(                                                   &
            &                    +timescaleShock                                    &
            &                    *orbitalFrequencySatellite                         &
            &                   )**2                                                &
            &                 )**self%gamma                                         &
            &                /3.0d0                                                 &
            &                *tidalTensor%doubleContract(tidalTensorPathIntegrated) &
            &                *(kilo*gigaYear/megaParsec)**2
       ! Limit the heating rate to be non-negative.
       gnedin1999HeatingRate=max(heatingRateNormalized,0.0d0)
    else
       ! Speed is non-positive - assume zero heating rate.
       gnedin1999HeatingRate=0.0d0
    end if
    return
  end function gnedin1999HeatingRate


