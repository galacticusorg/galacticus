!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !# <darkMatterProfile name="darkMatterProfileIsothermal">
  !#  <description>Isothermal dark matter halo profiles</description>
  !# </darkMatterProfile>

  use Dark_Matter_Halo_Scales

  type, extends(darkMatterProfileClass) :: darkMatterProfileIsothermal
     !% A dark matter halo profile class implementing isothermal dark matter halos.
     private
     class(darkMatterHaloScaleClass), pointer :: scale
   contains
     final                                             isothermalDestructor
     procedure :: calculationReset                  => isothermalCalculationReset
     procedure :: density                           => isothermalDensity
     procedure :: enclosedMass                      => isothermalEnclosedMass
     procedure :: potential                         => isothermalPotential
     procedure :: circularVelocity                  => isothermalCircularVelocity
     procedure :: circularVelocityMaximum           => isothermalCircularVelocityMaximum
     procedure :: radiusFromSpecificAngularMomentum => isothermalRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => isothermalRotationNormalization
     procedure :: energy                            => isothermalEnergy
     procedure :: energyGrowthRate                  => isothermalEnergyGrowthRate
     procedure :: kSpace                            => isothermalKSpace
     procedure :: freefallRadius                    => isothermalFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => isothermalFreefallRadiusIncreaseRate
  end type darkMatterProfileIsothermal

  interface darkMatterProfileIsothermal
     !% Constructors for the {\tt isothermal} dark matter halo profile class.
     module procedure isothermalDefaultConstructor
     module procedure isothermalConstructor
  end interface darkMatterProfileIsothermal

contains

  function isothermalDefaultConstructor()
    !% Default constructor for the {\tt isothermal} dark matter halo profile class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileIsothermal), target :: isothermalDefaultConstructor

    isothermalDefaultConstructor%scale => darkMatterHaloScale()
    return
  end function isothermalDefaultConstructor

  function isothermalConstructor(scale)
    !% Generic constructor for the {\tt isothermal} dark matter halo profile class.
    use Input_Parameters
    implicit none
    type (darkMatterProfileIsothermal), target :: isothermalConstructor
    class(darkMatterHaloScaleClass   ), target :: scale

    isothermalConstructor%scale => scale
    return
  end function isothermalConstructor
  
  subroutine isothermalDestructor(self)
    !% Destructor for the {\tt isothermal} dark matter halo profile class.
    implicit none
    type(darkMatterProfileIsothermal), intent(inout) :: self

    if (self%scale%isFinalizable()) deallocate(self%scale)
    return
  end subroutine isothermalDestructor
  
  subroutine isothermalCalculationReset(self,thisNode)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileIsothermal), intent(inout)          :: self
    type (treeNode                   ), intent(inout), pointer :: thisNode

    call self%scale%calculationReset(thisNode)
    return
  end subroutine isothermalCalculationReset

  double precision function isothermalDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\tt node} at the given {\tt radius} (given
    !% in units of Mpc).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileIsothermal), intent(inout)          :: self
    type            (treeNode                   ), intent(inout), pointer :: node
    double precision                             , intent(in   )          :: radius
    class           (nodeComponentBasic         )               , pointer :: thisBasicComponent

    thisBasicComponent   => node%basic         ()
    isothermalDensity=thisBasicComponent%mass()/4.0d0/Pi/self%scale%virialRadius(node)/radius**2
    return
  end function isothermalDensity

  double precision function isothermalEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt node} at the given {\tt radius} (given in
    !% units of Mpc).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileIsothermal), intent(inout)          :: self
    type            (treeNode                   ), intent(inout), pointer :: node
    double precision                             , intent(in   )          :: radius
    class           (nodeComponentBasic         )               , pointer :: thisBasicComponent

    thisBasicComponent   => node%basic     ()
    isothermalEnclosedMass=thisBasicComponent%mass()*(radius/self%scale%virialRadius(node))
    return
  end function isothermalEnclosedMass

  double precision function isothermalPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\tt node} at the given {\tt radius} (given in
    !% units of Mpc).
    use Galacticus_Nodes
    use Dark_Matter_Profiles_Error_Codes
    use Dark_Matter_Halo_Scales
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileIsothermal), intent(inout)           :: self
    type            (treeNode                   ), intent(inout), pointer  :: node
    double precision                             , intent(in   )           :: radius
    integer                                      , intent(  out), optional :: status
    double precision                             , parameter               :: radiusFractionalMinimum=1.0d-30
    double precision                                                       :: radiusFractional

    radiusFractional      =  radius/self%scale%virialRadius(node)
    if (radiusFractional <= 0.0d0) then
       isothermalPotential=0.0d0
       if (present(status)) status=darkMatterProfileErrorInfinite
    else
       isothermalPotential=(-1.0d0+log(radiusFractional))&
            &*self%scale%virialVelocity(node)**2
       if (present(status)) status=darkMatterProfileSuccess
    end if
    return
  end function isothermalPotential

  double precision function isothermalCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\tt node} at the given {\tt radius} (given in
    !% units of Mpc). For an isothermal halo this is independent of radius and therefore equal to the virial velocity.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileIsothermal), intent(inout)          :: self
    type            (treeNode                   ), intent(inout), pointer :: node
    double precision                             , intent(in   )          :: radius

    isothermalCircularVelocity=self%scale%virialVelocity(node)
    return
  end function isothermalCircularVelocity

  double precision function isothermalCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\tt node}. For an isothermal halo circular
    !% velocity is independent of radius.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileIsothermal), intent(inout)          :: self
    type            (treeNode                   ), intent(inout), pointer :: node

    isothermalCircularVelocityMaximum=self%circularVelocity(node,0.0d0)
    return
  end function isothermalCircularVelocityMaximum

  double precision function isothermalRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\tt node} at which a circular orbit has the given {\tt specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc). For an isothermal halo, the circular velocity is constant (and therefore equal to the virial
    !% velocity). Therefore, $r = j/V_{\rm virial}$ where $j$(={\tt specificAngularMomentum}) is the specific angular momentum and
    !% $r$ the required radius.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileIsothermal), intent(inout)          :: self
    type            (treeNode                   ), intent(inout), pointer :: node
    double precision                             , intent(in   )          :: specificAngularMomentum

    isothermalRadiusFromSpecificAngularMomentum=specificAngularMomentum/self%scale%virialVelocity(node)
    return
  end function isothermalRadiusFromSpecificAngularMomentum

  double precision function isothermalRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileIsothermal), intent(inout)          :: self
    type (treeNode                   ), intent(inout), pointer :: node

    isothermalRotationNormalization=2.0d0/self%scale%virialRadius(node)
    return
  end function isothermalRotationNormalization

  double precision function isothermalEnergy(self,node)
    !% Return the energy of an isothermal halo density profile.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileIsothermal), intent(inout)          :: self
    type (treeNode                   ), intent(inout), pointer :: node
    class(nodeComponentBasic         )               , pointer :: thisBasicComponent

    thisBasicComponent    => node%basic     ()
    isothermalEnergy=-0.5d0*thisBasicComponent%mass()*self%scale%virialVelocity(node)**2
    return
  end function isothermalEnergy

  double precision function isothermalEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of an isothermal halo density profile.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileIsothermal), intent(inout)          :: self
    type (treeNode                   ), intent(inout), pointer :: node
    class(nodeComponentBasic         )               , pointer :: thisBasicComponent

    thisBasicComponent   => node%basic     ()
    isothermalEnergyGrowthRate=self%energy(node)&
         &*(thisBasicComponent%accretionRate()/thisBasicComponent%mass()+2.0d0&
         &*self%scale%virialVelocityGrowthRate(node)/self%scale%virialVelocity(node))
    return
  end function isothermalEnergyGrowthRate

  double precision function isothermalKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the isothermal density profile at the specified {\tt waveNumber} (given in Mpc$^{-1}$), using the
    !% expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; table~1).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Exponential_Integrals
    implicit none
    class           (darkMatterProfileIsothermal), intent(inout)          :: self
    type            (treeNode                   ), intent(inout), pointer :: node
    double precision                             , intent(in   )          :: waveNumber
    double precision                                                      :: radiusScale, waveNumberScaleFree

    ! Get the scale radius (for which we use the virial radius).
    radiusScale          =  self%scale%virialRadius(node)

    ! Get the dimensionless wavenumber.
    waveNumberScaleFree=waveNumber*radiusScale

    ! Compute the Fourier transformed profile.
    isothermalKSpace=Sine_Integral(waveNumberScaleFree)/waveNumberScaleFree

    return
  end function isothermalKSpace

  double precision function isothermalFreefallRadius(self,node,time)
    !% Returns the freefall radius in the isothermal density profile at the specified {\tt time} (given in Gyr). For an isothermal
    !% potential, the freefall radius, $r_{\rm ff}(t)$, is:
    !% \begin{equation}
    !% r_{\rm ff}(t) = \sqrt{{2 \over \pi}} V_{\rm virial} t.
    !% \end{equation}
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    class           (darkMatterProfileIsothermal), intent(inout)          :: self
    type            (treeNode                   ), intent(inout), pointer :: node
    double precision                             , intent(in   )          :: time

    isothermalFreefallRadius=sqrt(2.0d0/Pi)*self%scale%virialVelocity(node)*time&
         &/Mpc_per_km_per_s_To_Gyr
    return
  end function isothermalFreefallRadius

  double precision function isothermalFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the isothermal density profile at the specified {\tt time} (given in
    !% Gyr). For an isothermal potential, the rate of increase of the freefall radius, $\dot{r}_{\rm ff}(t)$, is:
    !% \begin{equation}
    !% \dot{r}_{\rm ff}(t) = \sqrt{{2 \over \pi}} V_{\rm virial}.
    !% \end{equation}
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    class           (darkMatterProfileIsothermal), intent(inout)          :: self
    type            (treeNode                   ), intent(inout), pointer :: node
    double precision                             , intent(in   )          :: time

    isothermalFreefallRadiusIncreaseRate=sqrt(2.0d0/Pi)*self%scale%virialVelocity(node)&
         &/Mpc_per_km_per_s_To_Gyr
    return
  end function isothermalFreefallRadiusIncreaseRate
