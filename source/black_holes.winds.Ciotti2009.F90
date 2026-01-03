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
  Implements a black hole winds class based on the model of \cite{ciotti_feedbackcentral_2009}.
  !!}

  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  use :: Accretion_Disks           , only : accretionDisksClass

  !![
  <blackHoleWind name="blackHoleWindCiotti2009">
   <description>
     A black hole winds class based (loosely) on the model of \cite{ciotti_feedbackcentral_2009}. The wind power is given by:
     \begin{equation}
     L_\mathrm{w} = \epsilon_\mathrm{w} \epsilon_\mathrm{r} f_\mathrm{w} \dot{M}_\bullet \mathrm{c}^2,
     \end{equation}
     where $\dot{M}_\bullet$ is the black hole accretion rate, $\epsilon_\mathrm{w}=${\normalfont \ttfamily [efficiencyWind]} is
     an overall efficiency parameter, $\epsilon_\mathrm{r}$ is the radiative efficiency of the accretion flow \emph{if}
     {\normalfont \ttfamily [efficiencyWindScalesWithEfficiencyRadiative]=true}, and $1$ otherwise, and $f_\mathrm{w}$ represents
     the fraction of the wind power that is coupled to the surrounding galaxy.

     The model for $f_\mathrm{w}$ is inspired by \cite{ciotti_feedbackcentral_2009} who state:
     \begin{quotation}     
     If the pressure corresponding to the momentum ﬂow within the jet or wind is much greater than the pressure in the ambient
     gas, very little mass, momentum and kinetic energy is taken from it and deposited in that ambient gas. But when the [ratio of
     ISM to wind pressure] approaches unity, the ``working surface'' has been reached and the jet or wind discharges its content.
     \end{quotation}
     
     The energy density (pressure) in the wind at some radius $r$ in the galaxy is simply $\epsilon_\mathrm{w} \epsilon_\mathrm{r}
     \dot{M} \mathrm{c}^2 / 4 \pi r^2 v_\mathrm{w}$ (i.e. the energy input into a shell over time $\delta t$, $\epsilon_\mathrm{w}
     \epsilon_\mathrm{r} \dot{M} \mathrm{c}^2 \delta t$, divided by the volume of the shell occupied by the wind in time $\delta
     t$, $4 \pi r^2 v_\mathrm{w} \delta t$) where $v_\mathrm{w}$ is the wind velocity (assumed fixed at $10^4$~km/s). The
     corresponding ISM pressure is just $(3/2) \mathrm{k}T \rho(r)/m_\mathrm{H}$ where $T$ is the ISM temperature (assumed fixed
     at $10^4$~K) and $\rho$ is the ISM density. We approximate that the spheroid ISM density as $3 M/4/\pi/r^3$, where $M$ is the
     total gas mass in the spheroid, such that we find a ratio of ISM to wind pressures of:     
     \begin{equation}
      \frac{P_\mathrm{ISM}}{P_\mathrm{w}} = (3/2) \mathrm{k}T \rho(r)/m_\mathrm{H} \left/ \frac{\epsilon_\mathrm{w} \epsilon_\mathrm{r} \dot{M} \mathrm{c}^2}{4 \pi r^2 v_\mathrm{w}} \right. .
     \end{equation}
     We then smoothly interpolate $f_\mathrm{w}$ across the transition as:
     \begin{equation}
     f_\mathrm{w} = \left\{ \begin{array}{ll} 0 &amp; \hbox{if } x \le 0 \\ 3 x^2 - 2 x^3 &amp; \hbox{if } 0 &lt; x &lt; 1 \\ 1 &amp; \hbox{if } x \ge 1, \end{array} \right.
     \end{equation}
     where $x=P_\mathrm{ISM}/P_\mathrm{w}-1/2$.
    </description>
  </blackHoleWind>
  !!]
  type, extends(blackHoleWindClass) :: blackHoleWindCiotti2009
     !!{
     A black hole winds class based on the model of \cite{ciotti_feedbackcentral_2009}.
     !!}
     private
     class           (blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_                     => null()
     class           (accretionDisksClass        ), pointer :: accretionDisks_                             => null()
     double precision                                       :: efficiencyWind
     logical                                                :: efficiencyWindScalesWithEfficiencyRadiative
   contains
     final     ::          ciotti2009Destructor
     procedure :: power => ciotti2009Power
  end type blackHoleWindCiotti2009
  
  interface blackHoleWindCiotti2009
     !!{
     Constructors for the \refClass{blackHoleWindCiotti2009} black hole winds class.
     !!}
     module procedure ciotti2009ConstructorParameters
     module procedure ciotti2009ConstructorInternal
  end interface blackHoleWindCiotti2009

contains

  function ciotti2009ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleWindCiotti2009} black hole winds class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleWindCiotti2009    )                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    class           (blackHoleAccretionRateClass), pointer       :: blackHoleAccretionRate_
    class           (accretionDisksClass        ), pointer       :: accretionDisks_
    double precision                                             :: efficiencyWind
    logical                                                      :: efficiencyWindScalesWithEfficiencyRadiative

    !![
    <inputParameter>
      <name>efficiencyWind</name>
      <defaultValue>2.4d-3</defaultValue>
      <description>The efficiency of the black hole accretion-driven wind: $L_\mathrm{wind} = \epsilon_\mathrm{wind} \dot{M}_\bullet \clight^2$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>efficiencyWindScalesWithEfficiencyRadiative</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether the black hole wind efficiency should scale with the radiative efficiency of the accretion disk.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    !!]
    self=blackHoleWindCiotti2009(efficiencyWind,efficiencyWindScalesWithEfficiencyRadiative,blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function ciotti2009ConstructorParameters

  function ciotti2009ConstructorInternal(efficiencyWind,efficiencyWindScalesWithEfficiencyRadiative,blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Internal constructor for the \refClass{blackHoleWindCiotti2009} node operator class.
    !!}
    implicit none
    type            (blackHoleWindCiotti2009    )                        :: self
    class           (blackHoleAccretionRateClass), target, intent(in   ) :: blackHoleAccretionRate_
    class           (accretionDisksClass        ), target, intent(in   ) :: accretionDisks_
    double precision                                     , intent(in   ) :: efficiencyWind
    logical                                              , intent(in   ) :: efficiencyWindScalesWithEfficiencyRadiative
    !![
    <constructorAssign variables="efficiencyWind, efficiencyWindScalesWithEfficiencyRadiative, *blackHoleAccretionRate_, *accretionDisks_"/>
    !!]

    return
  end function ciotti2009ConstructorInternal

  subroutine ciotti2009Destructor(self)
    !!{
    Destructor for the ciotti2009 black hole winds class.
    !!}
    implicit none
    type(blackHoleWindCiotti2009), intent(inout) :: self
    
    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    <objectDestructor name="self%accretionDisks_"        />
    !!]
    return
  end subroutine ciotti2009Destructor
  
  double precision function ciotti2009Power(self,blackHole) result(power)
    !!{
    Compute the power of a black hole-driven wind that couples to the surrounding galaxy.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentSpheroid
    use :: Numerical_Constants_Astronomical, only : gigaYear             , megaParsec, massSolar
    use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant   , speedLight
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (blackHoleWindCiotti2009), intent(inout) :: self
    class           (nodeComponentBlackHole ), intent(inout) :: blackHole
    double precision                         , parameter     :: velocityWind         =1.0d4
    double precision                         , parameter     :: temperatureISM       =1.0d4
    class           (nodeComponentSpheroid  ), pointer       :: spheroid
    double precision                                         :: rateAccretionSpheroid                , rateAccretionHotHalo, &
         &                                                      rateAccretionNuclearStarCluster      , rateAccretion       , &
         &                                                      pressureWind                         , pressureISM         , &
         &                                                      fractionWindCoupled                  , massGasSpheroid     , &
         &                                                      radiusSpheroid                       , x                   , &
         &                                                      efficiencyWind
    
    ! Find the accretion rate and wind production efficiency.
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
    rateAccretion =+     rateAccretionSpheroid           &
         &         +     rateAccretionHotHalo            &
         &         +     rateAccretionNuclearStarCluster
    efficiencyWind=+self%efficiencyWind
    ! Apply any scaling of the wind production efficiency with the radiative efficiency of the accretion disk.
    if (self%efficiencyWindScalesWithEfficiencyRadiative)                                    &
         & efficiencyWind=+                     efficiencyWind                               &
         &                *self%accretionDisks_%efficiencyRadiative(blackHole,rateAccretion)
    ! Return if there is no wind.
    if     (                         &
         &   rateAccretion  <= 0.0d0 &
         &  .or.                     &
         &   efficiencyWind <= 0.0d0 &
         & ) then
       power=0.0d0
       return
    end if
    ! Compute the fraction of the wind that couples of the ISM of the spheroid component.
    spheroid            => blackHole%hostNode%spheroid()
    massGasSpheroid     =  spheroid          %massGas ()
    fractionWindCoupled =  0.0d0
    if (massGasSpheroid > 0.0d0) then
       radiusSpheroid=spheroid%radius()
       if (radiusSpheroid > 0.0d0) then
          ! Evaluate wind and ISM pressures.
          pressureWind=+efficiencyWind                       &
               &       *rateAccretion    *massSolar/gigaYear &
               &       *speedLight**2                        &
               &       /4.0d0                                &
               &       /Pi                                   &
               &       /velocityWind     /kilo               &
               &       /radiusSpheroid**2/megaParsec**2
          pressureISM =+3.0d0                                &
               &       /4.0d0                                &
               &       /Pi                                   &
               &       *massGasSpheroid  *massSolar          &
               &       /radiusSpheroid**3/megaParsec**3      &
               &       /massHydrogenAtom                     &
               &       *3.0d0                                &
               &       /2.0d0                                &
               &       *boltzmannsConstant                   &
               &       *temperatureISM
          ! Construct an interpolating factor such that the energy input from the wind drops to zero below half of the ISM
          ! pressure.
          x      =+pressureISM  &
               &  /pressureWind &
               &  -0.50d0
          if (x <= 0.0d0) then
             ! No energy input below half of critical density.
             fractionWindCoupled=+0.0d0
          else if (x >= 1.0d0) then
             ! Full energy input above 1.5 times critical density.
             fractionWindCoupled=+1.0d0
          else
             ! Smooth polynomial interpolating function between these limits.
             fractionWindCoupled=+3.0d0*x**2 &
                  &              -2.0d0*x**3
          end if
       end if
    end if
    efficiencyWind=+efficiencyWind      &
         &         *fractionWindCoupled
    ! Evaluate the wind power.
    power  =+efficiencyWind    &
         &  *rateAccretion     &
         &  *speedLight    **2 &
         &  /kilo          **2
    return
  end function ciotti2009Power
