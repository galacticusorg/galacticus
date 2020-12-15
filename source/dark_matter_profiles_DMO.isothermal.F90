!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% An implementation of isothermal dark matter halo profiles.

  !# <darkMatterProfileDMO name="darkMatterProfileDMOIsothermal">
  !#  <description>
  !#   A dark matter profile DMO class in which the density profile is given by:
  !#   \begin{equation}
  !#    \rho_\mathrm{dark matter}(r) \propto r^{-2},
  !#   \end{equation}
  !#   normalized such that the total mass of the \gls{node} is enclosed with the virial radius.
  !#  </description>
  !# </darkMatterProfileDMO>
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOIsothermal
     !% A dark matter halo profile class implementing isothermal dark matter halos.
     private
   contains
     final     ::                                      isothermalDestructor
     procedure :: density                           => isothermalDensity
     procedure :: densityLogSlope                   => isothermalDensityLogSlope
     procedure :: radialMoment                      => isothermalRadialMoment
     procedure :: enclosedMass                      => isothermalEnclosedMass
     procedure :: radiusEnclosingDensity            => isothermalRadiusEnclosingDensity
     procedure :: potential                         => isothermalPotential
     procedure :: circularVelocity                  => isothermalCircularVelocity
     procedure :: circularVelocityMaximum           => isothermalCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => isothermalRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => isothermalRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => isothermalRotationNormalization
     procedure :: energy                            => isothermalEnergy
     procedure :: energyGrowthRate                  => isothermalEnergyGrowthRate
     procedure :: kSpace                            => isothermalKSpace
     procedure :: freefallRadius                    => isothermalFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => isothermalFreefallRadiusIncreaseRate
  end type darkMatterProfileDMOIsothermal

  interface darkMatterProfileDMOIsothermal
     !% Constructors for the {\normalfont \ttfamily isothermal} dark matter halo profile class.
     module procedure isothermalConstructorParameters
     module procedure isothermalConstructorInternal
  end interface darkMatterProfileDMOIsothermal

contains

  function isothermalConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily isothermal} dark matter halo profile class.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDMOIsothermal)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_

    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=darkMatterProfileDMOIsothermal(darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function isothermalConstructorParameters

  function isothermalConstructorInternal(darkMatterHaloScale_) result(self)
    !% Generic constructor for the {\normalfont \ttfamily isothermal} dark matter halo profile class.
    implicit none
    type (darkMatterProfileDMOIsothermal)                        :: self
    class(darkMatterHaloScaleClass      ), intent(in   ), target :: darkMatterHaloScale_
    !# <constructorAssign variables="*darkMatterHaloScale_"/>

    return
  end function isothermalConstructorInternal

  subroutine isothermalDestructor(self)
    !% Destructor for the {\normalfont \ttfamily isothermal} dark matter halo profile class.
    implicit none
    type(darkMatterProfileDMOIsothermal), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_" />
    return
  end subroutine isothermalDestructor

  double precision function isothermalDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given
    !% in units of Mpc).
    use :: Galacticus_Nodes        , only : nodeComponentBasic, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentBasic            ), pointer       :: basic

    basic             => node %basic()
    isothermalDensity =  basic%mass ()/4.0d0/Pi/self%darkMatterHaloScale_%virialRadius(node)/radius**2
    return
  end function isothermalDensity

  double precision function isothermalDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, radius

    isothermalDensityLogSlope=-2.0d0
    return
  end function isothermalDensityLogSlope

  double precision function isothermalRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given
    !% in units of Mpc).
    use :: Galacticus_Nodes        , only : nodeComponentBasic, treeNode
    use :: Numerical_Comparison    , only : Values_Agree
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout)           :: self
    type            (treeNode                      ), intent(inout)           :: node
    double precision                                , intent(in   )           :: moment
    double precision                                , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    class           (nodeComponentBasic            )               , pointer  :: basic
    double precision                                                          :: radiusMinimumActual, radiusMaximumActual

    radiusMinimumActual=0.0d0
    radiusMaximumActual=self%darkMatterHaloScale_%virialRadius(node)
    if (present(radiusMinimum)) radiusMinimumActual=radiusMinimum
    if (present(radiusMaximum)) radiusMaximumActual=radiusMaximum
    basic   => node%basic         ()
    if (Values_Agree(moment,1.0d0,absTol=1.0d-6)) then
       isothermalRadialMoment=+basic%mass()                                 &
            &                 /4.0d0                                        &
            &                 /Pi                                           &
            &                 /self%darkMatterHaloScale_%virialRadius(node) &
            &                 *log(                                         &
            &                      +radiusMaximumActual                     &
            &                      /radiusMinimumActual                     &
            &                     )
    else
       isothermalRadialMoment=+basic%mass()                                 &
            &                 /4.0d0                                        &
            &                 /Pi                                           &
            &                 /self%darkMatterHaloScale_%virialRadius(node) &
            &                 /                       (moment-1.0d0)        &
            &                 *(                                            &
            &                   +radiusMaximumActual**(moment-1.0d0)        &
            &                   -radiusMinimumActual**(moment-1.0d0)        &
            &                  )
    end if
    return
  end function isothermalRadialMoment

  double precision function isothermalEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentBasic            ), pointer       :: basic

    basic   => node%basic     ()
    isothermalEnclosedMass=basic%mass()*(radius/self%darkMatterHaloScale_%virialRadius(node))
    return
  end function isothermalEnclosedMass

  double precision function isothermalPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use :: Galactic_Structure_Options, only : structureErrorCodeInfinite, structureErrorCodeSuccess
    use :: Galacticus_Error          , only : Galacticus_Error_Report
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout)           :: self
    type            (treeNode                      ), intent(inout), target   :: node
    double precision                                , intent(in   )           :: radius
    integer                                         , intent(  out), optional :: status
    double precision                                , parameter               :: radiusFractionalMinimum=1.0d-30
    double precision                                                       :: radiusFractional

    radiusFractional      =  radius/self%darkMatterHaloScale_%virialRadius(node)
    if (radiusFractional <= 0.0d0) then
       isothermalPotential=0.0d0
       if (present(status)) status=structureErrorCodeInfinite
    else
       isothermalPotential=log(radiusFractional)*self%darkMatterHaloScale_%virialVelocity(node)**2
       if (present(status)) status=structureErrorCodeSuccess
    end if
    return
  end function isothermalPotential

  double precision function isothermalCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc). For an isothermal halo this is independent of radius and therefore equal to the virial velocity.
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    !$GLC attributes unused :: radius

    isothermalCircularVelocity=self%darkMatterHaloScale_%virialVelocity(node)
    return
  end function isothermalCircularVelocity

  double precision function isothermalCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}. For an isothermal halo circular
    !% velocity is independent of radius.
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node

    isothermalCircularVelocityMaximum=self%circularVelocity(node,0.0d0)
    return
  end function isothermalCircularVelocityMaximum

  double precision function isothermalRadialVelocityDispersion(self,node,radius)
    !% Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius}
    !% (given in units of Mpc). For an isothermal halo this is independent of radius and equal to the virial velocity divided by $\sqrt(2)$.
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    !$GLC attributes unused :: radius

    isothermalRadialVelocityDispersion=self%darkMatterHaloScale_%virialVelocity(node)/sqrt(2.0d0)
    return
  end function isothermalRadialVelocityDispersion

  double precision function isothermalRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc). For an isothermal halo, the circular velocity is constant (and therefore equal to the virial
    !% velocity). Therefore, $r = j/V_\mathrm{virial}$ where $j$(={\normalfont \ttfamily specificAngularMomentum}) is the specific angular momentum and
    !% $r$ the required radius.
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: specificAngularMomentum

    isothermalRadiusFromSpecificAngularMomentum=specificAngularMomentum/self%darkMatterHaloScale_%virialVelocity(node)
    return
  end function isothermalRadiusFromSpecificAngularMomentum

  double precision function isothermalRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    implicit none
    class(darkMatterProfileDMOIsothermal), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node

    isothermalRotationNormalization=2.0d0/self%darkMatterHaloScale_%virialRadius(node)
    return
  end function isothermalRotationNormalization

  double precision function isothermalEnergy(self,node)
    !% Return the energy of an isothermal halo density profile.
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(darkMatterProfileDMOIsothermal), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentBasic            ), pointer       :: basic

    basic            =>       node%basic ()
    isothermalEnergy =  -0.5d0*basic%mass()*self%darkMatterHaloScale_%virialVelocity(node)**2
    return
  end function isothermalEnergy

  double precision function isothermalEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of an isothermal halo density profile.
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(darkMatterProfileDMOIsothermal), intent(inout)          :: self
    type (treeNode                      ), intent(inout), target  :: node
    class(nodeComponentBasic            )               , pointer :: basic

    basic   => node%basic     ()
    isothermalEnergyGrowthRate=self%energy(node)                                                                                                &
         &                     *(                                                                                                               &
         &                       +basic%accretionRate()/basic%mass()                                                                            &
         &                       +2.0d0*self%darkMatterHaloScale_%virialVelocityGrowthRate(node)/self%darkMatterHaloScale_%virialVelocity(node) &
         &                      )
    return
  end function isothermalEnergyGrowthRate

  double precision function isothermalKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the isothermal density profile at the specified {\normalfont \ttfamily waveNumber} (given in Mpc$^{-1}$), using the
    !% expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; table~1).
    use :: Exponential_Integrals, only : Sine_Integral
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), target  :: node
    double precision                                , intent(in   )          :: waveNumber
    double precision                                                         :: radiusScale, waveNumberScaleFree

    ! Get the scale radius (for which we use the virial radius).
    radiusScale          =  self%darkMatterHaloScale_%virialRadius(node)

    ! Get the dimensionless wavenumber.
    waveNumberScaleFree=waveNumber*radiusScale

    ! Compute the Fourier transformed profile.
    isothermalKSpace=Sine_Integral(waveNumberScaleFree)/waveNumberScaleFree

    return
  end function isothermalKSpace

  double precision function isothermalFreefallRadius(self,node,time)
    !% Returns the freefall radius in the isothermal density profile at the specified {\normalfont \ttfamily time} (given in Gyr). For an isothermal
    !% potential, the freefall radius, $r_\mathrm{ff}(t)$, is:
    !% \begin{equation}
    !% r_\mathrm{ff}(t) = \sqrt{{2 \over \pi}} V_\mathrm{virial} t.
    !% \end{equation}
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: time

    isothermalFreefallRadius=sqrt(2.0d0/Pi)*self%darkMatterHaloScale_%virialVelocity(node)*time&
         &/Mpc_per_km_per_s_To_Gyr
    return
  end function isothermalFreefallRadius

  double precision function isothermalFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the isothermal density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr). For an isothermal potential, the rate of increase of the freefall radius, $\dot{r}_\mathrm{ff}(t)$, is:
    !% \begin{equation}
    !% \dot{r}_\mathrm{ff}(t) = \sqrt{{2 \over \pi}} V_\mathrm{virial}.
    !% \end{equation}
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: time
    !$GLC attributes unused :: time

    isothermalFreefallRadiusIncreaseRate=sqrt(2.0d0/Pi)*self%darkMatterHaloScale_%virialVelocity(node) &
         & /Mpc_per_km_per_s_To_Gyr
    return
  end function isothermalFreefallRadiusIncreaseRate

  double precision function isothermalRadiusEnclosingDensity(self,node,density)
    !% Null implementation of function to compute the radius enclosing a given density for isothermal dark matter halo profiles.
    use :: Galacticus_Nodes        , only : nodeComponentBasic, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOIsothermal), intent(inout), target :: self
    type            (treeNode                      ), intent(inout), target :: node
    double precision                                , intent(in   )         :: density
    class           (nodeComponentBasic            ), pointer               :: basic

    basic                            =>             node                %basic       (    )
    isothermalRadiusEnclosingDensity =  +sqrt(                                              &
         &                                    +     basic               %mass        (    ) &
         &                                    /4.0d0                                        &
         &                                    *3.0d0                                        &
         &                                    /Pi                                           &
         &                                    /self%darkMatterHaloScale_%virialRadius(node) &
         &                                    /density                                      &
         &                                   )
    return
  end function isothermalRadiusEnclosingDensity

